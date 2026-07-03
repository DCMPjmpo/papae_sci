"""
extract_site.py 【5蛋白全量最终定稿｜版本锁定不再迭代】
精细化修复：1.差值二次强转float16  2.非法位点拦截防CLS误取
已完成全部生产级加固：
- WT numpy.copy 隔离内存视图
- 断点续跑计数对齐修复
- 浮点类型全程严格锁定float16
- 非法突变位置前置拦截
- 2000间隔flush降低IO压力
- 跑完自动清理断点避免下次误续跑
- 显式1-based位点注释杜绝索引混淆
- 全程信任mutation_position字段，不解析突变字符串
适配P6五蛋白SCI跨蛋白泛化研究，上下游路径匹配config.py
输入：config.MERGE_EXPANDED(all_proteins_expanded.csv)、config.WT_SEQUENCES(wt_sequences.csv)
输出：
    all_proteins_delta_site.dat (memmap float16 (N,33,1280))
    all_proteins_site_metadata.csv
"""
import pandas as pd
import numpy as np
import torch
import os
import sys
import random
import json
from tqdm import tqdm
from collections import defaultdict

# 导入全局配置
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config

# ===================== 全局运行配置（正式全量固定配置）=====================
DEBUG_ONLY_LAST_LAYER = False               # 正式实验固定False，提取全部33层
ENABLE_FREQUENT_CUDA_CACHE_CLEAR = False    # 关闭频繁empty_cache避免性能损耗
SAMPLE_ASSERT_CHECK_NUM = 100               # 抽样100条自动校验氨基酸匹配
SAVE_CHECKPOINT = True                       # 开启断点续跑应对长时间任务中断
CHECKPOINT_INTERVAL = 500                   # 每500个唯一序列保存断点
FLUSH_INTERVAL = 2000                        # 优化：每2000序列刷新memmap，降低频繁磁盘IO
CLEAR_CACHE_SEQ_INTERVAL = 1000
CHECKPOINT_PATH = os.path.join(config.OUTPUT_DIR, "checkpoint.json")
# ======================================================================

print("=" * 80)
print("extract_site.py | 五蛋白Delta位点提取 最终定型版（类型锁死+非法位点拦截）")
if DEBUG_ONLY_LAST_LAYER:
    print("运行模式：调试模式（仅ESM第33层，快速验证Pipeline正确性）")
else:
    print("运行模式：正式全量模式（提取ESM 1~33全部层，适配P1层间SCI跨蛋白泛化研究）")
print("=" * 80)

# ========== 加载展开后数据集 ==========
df = pd.read_csv(config.MERGE_EXPANDED)
total_samples = len(df)
print(f"展开后总样本条数: {total_samples:,}")

unique_seq_df = df[["protein", "mutated_sequence"]].drop_duplicates()
unique_seq_cnt = len(unique_seq_df)
print(f"唯一(蛋白+突变序列)总数: {unique_seq_cnt:,}")

# 统计去重前向推理节省比例
forward_saving_ratio = 1 - (unique_seq_cnt / total_samples)
print(f"序列去重前向推理算力节省比例: {forward_saving_ratio:.2%}")

# 加载WT序列映射表
wt_df = pd.read_csv(config.WT_SEQUENCES)
wt_map = dict(zip(wt_df['protein'], wt_df['wt_sequence']))
wt_protein_list = list(wt_map.keys())
print(f"\n参与计算蛋白列表: {wt_protein_list}")

# ========== ESM2模型初始化 ==========
print(f"\n加载ESM2模型: {config.ESM_MODEL_NAME}...")
import esm
model, alphabet = esm.pretrained.load_model_and_alphabet(config.ESM_MODEL_NAME)
model = model.to(config.DEVICE)
model.eval()
batch_converter = alphabet.get_batch_converter()

# 动态设置提取层数
if DEBUG_ONLY_LAST_LAYER:
    extract_layers = [config.ESM_N_LAYERS]
    effective_n_layers = 1
else:
    extract_layers = list(range(1, config.ESM_N_LAYERS + 1))
    effective_n_layers = config.ESM_N_LAYERS

# ========== 工具函数（显式位点对齐注释，杜绝长期索引疑惑） ==========
def extract_full_embedding(sequence, label=""):
    """提取整条序列分层Embedding，返回CPU tensor (n_layer, seq_len, 1280)"""
    data = [(label, sequence)]
    _, _, tokens = batch_converter(data)
    tokens = tokens.to(config.DEVICE)
    with torch.no_grad():
        results = model(tokens, repr_layers=extract_layers, return_contacts=False)
    layer_reps = []
    for lid in extract_layers:
        layer_reps.append(results["representations"][lid][0].cpu())
    return torch.stack(layer_reps)

def get_site_from_full(full_emb, position):
    """
    截取指定突变位点表征
    position: 传入CSV内mutation_position，为1-based位点
    ESM Token规则：
        token 0 → CLS占位符，无氨基酸对应
        token 1 → 序列第1个氨基酸
        token N → 序列第N个氨基酸
    因此mutation_position(1-based)可直接作为索引取值，无需±1偏移
    return: float16 tensor shape (n_layer, 1280)
    """
    return full_emb[:, position, :].half()

# ========== Step1 WT推理：新增非法位点校验 + 独立numpy缓存防内存共享 ==========
print(f"\n{'='*70}")
print("Step 1: WT前向推理 → 位点合法性校验 + 独立Numpy缓存，杜绝CLS误取&内存共享隐患")
print(f"{'='*70}")
wt_site_cache = dict()  # key=(protein, pos), value=np.float16 (n_layer,1280)

for prot_name, wt_seq in wt_map.items():
    seq_len_wt = len(wt_seq)
    print(f"\n  {prot_name} | WT序列长度: {seq_len_wt}")
    prot_mask = df["protein"] == prot_name
    needed_pos = df.loc[prot_mask, "mutation_position"].unique()
    if len(needed_pos) == 0:
        print(f"    该蛋白无突变位点，跳过推理")
        continue
    wt_full_emb = extract_full_embedding(wt_seq, f"{prot_name}_WT")
    print(f"    WT完整Embedding shape: {wt_full_emb.shape}")
    for pos in needed_pos:
        # 新增修复：拦截小于1的非法位点，防止误取CLS token
        if pos < 1:
            raise ValueError(f"【数据致命错误】蛋白{prot_name}检测到非法突变位点 position={pos}，小于1会读取CLS占位，终止运行！")
        site_tensor = get_site_from_full(wt_full_emb, pos)
        # .copy() 生成独立内存数组，杜绝numpy视图共享底层buffer隐性错乱
        wt_site_cache[(prot_name, pos)] = site_tensor.numpy().copy()
    del wt_full_emb
    print(f"    缓存位点数量: {len(needed_pos)}，释放整条WT Embedding ✅")
print(f"\n✅ WT 独立Numpy缓存完成，总缓存位点数量: {len(wt_site_cache)}")

# ========== Step2 初始化Memmap输出 & 断点续跑判断 ==========
print(f"\n{'='*70}")
print("Step 2: 初始化Delta Memmap磁盘映射文件")
print(f"{'='*70}")
output_shape = (total_samples, effective_n_layers, config.ESM_DIM)
memmap_path = config.DELTA_SITE.replace(".npy", ".dat")
resume_from_seq_idx = 0

# 断点续跑逻辑
if SAVE_CHECKPOINT and os.path.exists(CHECKPOINT_PATH):
    with open(CHECKPOINT_PATH, "r", encoding="utf-8") as f:
        ckpt_data = json.load(f)
    resume_from_seq_idx = ckpt_data.get("last_seq_idx", 0)
    print(f"检测到断点文件，将从第 {resume_from_seq_idx} 个唯一序列接续运行")
    delta_mem = np.memmap(memmap_path, dtype=np.float16, mode="r+", shape=output_shape)
else:
    delta_mem = np.memmap(memmap_path, dtype=np.float16, mode="w+", shape=output_shape)

est_gb = np.prod(output_shape) * 2 / (1024 ** 3)
print(f"Memmap路径: {memmap_path}")
print(f"输出张量Shape: {output_shape}")
print(f"预估磁盘占用: {est_gb:.2f} GB")

# ========== Step3 构建唯一序列索引（itertuples提速，替换低效iterrows） ==========
print(f"\n{'='*70}")
print("Step 3: 构建唯一序列-样本索引映射")
print(f"{'='*70}")
seq_to_indices = defaultdict(list)
for idx, row in enumerate(df.itertuples(index=False)):
    key = (row.protein, row.mutated_sequence)
    pos_int = int(row.mutation_position)
    seq_to_indices[key].append((idx, pos_int))
seq_items = list(seq_to_indices.items())
print(f"索引构建完成，唯一序列总数: {len(seq_items)}")

# 全局位点越界前置预检（校验mutation_position不超过序列总长）
print("\n全局突变位点越界预检:")
for prot_name, wt_seq in wt_map.items():
    wt_len = len(wt_seq)
    pos_arr = df[df["protein"] == prot_name]["mutation_position"].unique()
    if len(pos_arr) == 0:
        continue
    max_p = pos_arr.max()
    if max_p > wt_len:
        raise ValueError(f"数据异常：{prot_name} 最大突变位置{max_p} > WT序列长度{wt_len}")
    print(f"  {prot_name}: 最大位点{max_p}, WT总长{wt_len} | 合法 ✅")

# ========== Step4 主推理循环：断点计数修复 + 浮点类型全程锁死 + 周期Flush ==========
print(f"\n{'='*70}")
print("Step 4: Mut序列推理 → Delta计算 → 实时写入Memmap（断点续跑启用）")
print(f"{'='*70}")
processed_sample_cnt = 0
# 断点修复：续跑初始化对齐序号，保证断点间隔判断逻辑永远正确
seq_process_cnt = resume_from_seq_idx

for seq_idx, ((prot_name, mut_seq), idx_pos_list) in enumerate(tqdm(seq_items, desc="序列推理进度", total=len(seq_items))):
    # 跳过已处理完成的序列
    if seq_idx < resume_from_seq_idx:
        for sample_idx, pos in idx_pos_list:
            processed_sample_cnt += 1
        continue
    
    mut_full_emb = extract_full_embedding(mut_seq, f"{prot_name}_MUT")
    for sample_idx, pos in idx_pos_list:
        mut_tensor = get_site_from_full(mut_full_emb, pos)
        wt_np = wt_site_cache[(prot_name, pos)]
        # 本次修复：差值计算后二次强转float16，彻底锁死类型，规避numpy版本自动升精度
        delta_np = (mut_tensor.numpy().astype(np.float16) - wt_np).astype(np.float16)
        delta_mem[sample_idx] = delta_np
        processed_sample_cnt += 1
    del mut_full_emb

    seq_process_cnt += 1

    # 周期性刷新磁盘，降低高频IO开销
    if seq_process_cnt % FLUSH_INTERVAL == 0:
        delta_mem.flush()
    # 可选显存清理开关
    if ENABLE_FREQUENT_CUDA_CACHE_CLEAR and torch.cuda.is_available():
        if seq_process_cnt % CLEAR_CACHE_SEQ_INTERVAL == 0:
            torch.cuda.empty_cache()
    # 周期性落地断点文件
    if SAVE_CHECKPOINT and seq_process_cnt % CHECKPOINT_INTERVAL == 0:
        ckpt_data = {"last_seq_idx": seq_idx + 1}
        with open(CHECKPOINT_PATH, "w", encoding="utf-8") as f:
            json.dump(ckpt_data, f, indent=2)

# 全部循环结束强制落盘缓冲区
delta_mem.flush()

# ========== 运行完成自动删除断点文件，防止下次启动误续跑 ==========
if SAVE_CHECKPOINT:
    if os.path.exists(CHECKPOINT_PATH):
        os.remove(CHECKPOINT_PATH)
        print(f"\n✅ 任务全部完成，自动清理断点文件: {CHECKPOINT_PATH}，避免下次运行直接跳过全部数据")

print(f"\n✅ 处理完成，总处理样本: {processed_sample_cnt:,}")

# ========== Step5 导出Metadata表格 ==========
print(f"\n{'='*70}")
print("Step 5: 导出位点Metadata CSV")
print(f"{'='*70}")
metadata_df = df[[
    "protein", "mutant", "DMS_score", "DMS_score_bin",
    "mutation_position", "mutation_index", "n_mutations", "dataset"
]].copy()
if "parent_mutant" in df.columns:
    metadata_df["parent_mutant"] = df["parent_mutant"]
else:
    metadata_df["parent_mutant"] = df["mutant"]
metadata_df = metadata_df[[
    "protein", "mutant", "parent_mutant", "mutation_position",
    "mutation_index", "n_mutations", "DMS_score", "DMS_score_bin", "dataset"
]]
metadata_df.to_csv(config.SITE_METADATA, index=False)
print(f"Metadata保存路径: {config.SITE_METADATA}")
print(f"表格总行数: {len(metadata_df):,}")

# ========== Step6 全自动正确性断言校验（校验上游expand拆分逻辑） ==========
print(f"\n{'='*70}")
print(f"Step 6: 自动校验（抽样{SAMPLE_ASSERT_CHECK_NUM}条氨基酸位点匹配断言，校验上游数据正确性）")
print(f"{'='*70}")
delta_check = np.memmap(memmap_path, dtype=np.float16, mode="r", shape=output_shape)
assert delta_check.shape[0] == len(metadata_df), f"样本数量不匹配 memmap:{delta_check.shape[0]} vs metadata:{len(metadata_df)}"
assert delta_check.shape[1] == effective_n_layers, "提取层数维度不匹配"
assert delta_check.shape[2] == config.ESM_DIM, "Embedding特征维度不匹配"

random.seed(42)
total_meta = len(metadata_df)
check_indices = random.sample(range(total_meta), min(SAMPLE_ASSERT_CHECK_NUM, total_meta))
error_flag = False

for sid in check_indices:
    row_meta = metadata_df.iloc[sid]
    prot = row_meta["protein"]
    pos_1 = int(row_meta["mutation_position"])
    mut_anno = row_meta["mutant"]
    wt_seq = wt_map[prot]
    wt_aa_exp = wt_seq[pos_1 - 1]
    wt_aa_anno = mut_anno[0]
    mut_aa_anno = mut_anno[-1]
    original_df_row = df.iloc[sid]
    curr_mut_seq = original_df_row["mutated_sequence"]
    mut_aa_seq = curr_mut_seq[pos_1 - 1]
    try:
        assert wt_aa_exp == wt_aa_anno, f"【上游数据错误】样本{sid} {mut_anno}: WT注释{wt_aa_anno} 实际序列{wt_aa_exp} 位置{pos_1}"
        assert mut_aa_seq == mut_aa_anno, f"【上游数据错误】样本{sid} {mut_anno}: Mut注释{mut_aa_anno} 实际序列{mut_aa_seq} 位置{pos_1}"
    except AssertionError as e:
        print(e)
        error_flag = True

if error_flag:
    raise RuntimeError("校验失败：突变位置/氨基酸解析错位，问题来源于expand_multi_mutations.py上游脚本，终止运行！")
else:
    print(f"✅ {len(check_indices)}条抽样断言全部通过，上游位点拆分逻辑无误，本脚本取值逻辑安全")

rand_idx = random.randint(0, total_meta - 1)
vec_mean = delta_check[rand_idx].mean()
vec_std = delta_check[rand_idx].std()
print(f"随机样本{rand_idx} Delta向量均值: {vec_mean:.6f}, 标准差: {vec_std:.6f}")

del delta_mem
del delta_check

# ========== 收尾汇总 ==========
print(f"\n{'='*80}")
print("✅ extract_site.py 完全定型，可提交服务器跑全量；下一步转向P1 SCI计算+跨蛋白泛化实验")
print(f"输出文件汇总：")
print(f"  Delta Memmap文件: {memmap_path} | 大小: {os.path.getsize(memmap_path)/1024**3:.2f} GB")
print(f"  Metadata表格: {config.SITE_METADATA}")
print(f"\nP1固定读取写法（规避后缀路径坑）：")
print(f'delta_dat_path = config.DELTA_SITE.replace(".npy", ".dat")')
print(f'import numpy as np')
print(f'delta = np.memmap(delta_dat_path, dtype=np.float16, mode="r", shape={output_shape})')
print("=" * 80)