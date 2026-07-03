"""
P1_build_sci.py 【五蛋白95142样本 最终封版生产脚本｜永久定型不再迭代】
本次仅做可观测性增强，算法内核完全不变，符合科研严谨性建议：
1. 新增统计：存在零方差层样本数量、相关系数NaN/Inf异常样本数量
2. 出现异常样本实时打印告警，便于溯源，不会静默掩盖数据细节
3. 保留原有nan_to_num兜底：无定义相关系数直接置0，拒绝人为加随机噪声
4. 不做clip、不做标准化、不修改原始Delta数值，保留极端突变信号
已固化全部生产级修复与优化：
✅ 修复delta memmap一维解包崩溃，由metadata获取样本量显式指定三维shape
✅ SCI矩阵动态r+/w+模式，真正支持矩阵断点续跑，不会覆盖已有运算结果
✅ 三类SCI得分独立memmap实现真断点，重启不会清零历史计算值
✅ np.corrcoef批量相关矩阵，替代双重pearsonr循环提速数十倍
✅ np.nan_to_num清洗全零层NaN/Inf异常值，避免得分污染
✅ np.partition局部极值筛选，替代全局全排序节约运算开销
✅ SCI矩阵强制float32，杜绝半精度相关系数精度损失
✅ 运算结束自动导出.npy文件，P2~P5旧代码无需修改直接兼容读取
✅ 循环外预计算上三角索引，减少重复冗余运算
✅ 断点写入前全局flush所有memmap，杜绝断点超前数据导致样本丢失
✅ 周期断点+运行收尾自动清理断点文件，规避下次误续跑
"""
import numpy as np
import pandas as pd
import os
import sys
import json
from tqdm import tqdm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config

# ===================== 全局运行配置（生产固定参数，禁止修改）=====================
SAVE_FULL_SCI_MATRIX = True        # 强制存储完整SCI矩阵，支撑P2重布线、P3蛋白对比、P6跨蛋白泛化
SAVE_CHECKPOINT = True              # 开启断点续跑应对长时间任务中断
CHECKPOINT_INTERVAL = 1000         # 每1000个样本落地一次断点
CHECKPOINT_PATH = os.path.join(config.OUTPUT_DIR, "p1_checkpoint.json")
# ===================================================================================

print("=" * 80)
print("P1_build_sci.py | 层间SCI计算 最终封版定稿（新增异常统计监控，算法无改动）")
print("模式：强制生成完整SCI 33×33矩阵 + Memmap断点得分 + 末尾兼容npy导出")
print("=" * 80)

# 路径定义（严格匹配extract_site后缀替换规则，杜绝路径找不到bug）
delta_dat_path = config.DELTA_SITE.replace(".npy", ".dat")
sci_mat_out_path = config.SCI_MATRICES.replace(".npy", ".dat")
# SCI得分memmap路径（断点持久化用）
sci_mean_path = config.SCI_SCORES_MEAN.replace(".npy", ".dat")
sci_top20_path = config.SCI_SCORES_TOP20.replace(".npy", ".dat")
sci_top50_path = config.SCI_SCORES_TOP50.replace(".npy", ".dat")

# ========== 读取metadata获取样本总数，显式指定delta memmap三维shape，修复一维解包崩溃 ==========
meta_csv_path = config.SITE_METADATA
metadata = pd.read_csv(meta_csv_path)
n_sample = len(metadata)
n_layer = config.ESM_N_LAYERS
dim_emb = config.ESM_DIM
print(f"总样本量: {n_sample:,}, ESM层数: {n_layer}, 嵌入维度: {dim_emb}")

delta_mm = np.memmap(
    delta_dat_path,
    dtype=np.float16,
    mode="r",
    shape=(n_sample, n_layer, dim_emb)
)

# ========== SCI矩阵动态判断打开模式，实现矩阵完整断点续跑 ==========
if os.path.exists(sci_mat_out_path):
    sci_mat_mode = "r+"
    print(f"检测到已有SCI矩阵文件，以续写模式打开，保留已计算数据")
else:
    sci_mat_mode = "w+"
    print(f"新建SCI矩阵memmap文件")

sci_mat_mm = np.memmap(
    sci_mat_out_path,
    dtype=np.float32,
    mode=sci_mat_mode,
    shape=(n_sample, n_layer, n_layer)
)
print(f"SCI完整矩阵输出路径: {sci_mat_out_path}")

# -------------------------- SCI得分Memmap断点初始化 --------------------------
def init_score_memmap(file_path, total_n):
    """断点感知创建得分memmap：存在则续写，不存在则新建初始化"""
    if os.path.exists(file_path):
        return np.memmap(file_path, dtype=np.float32, mode="r+", shape=total_n)
    else:
        return np.memmap(file_path, dtype=np.float32, mode="w+", shape=total_n)

sci_mean_mm = init_score_memmap(sci_mean_path, n_sample)
sci_top20_mm = init_score_memmap(sci_top20_path, n_sample)
sci_top50_mm = init_score_memmap(sci_top50_path, n_sample)

# 断点续跑起始下标
start_idx = 0
if SAVE_CHECKPOINT and os.path.exists(CHECKPOINT_PATH):
    with open(CHECKPOINT_PATH, "r", encoding="utf-8") as f:
        ckpt_info = json.load(f)
    start_idx = ckpt_info.get("processed_sample_idx", 0)
    print(f"\n检测P1断点文件，从第 {start_idx:,} 号样本接续运算")

# 循环外预计算上三角掩码，避免重复生成
triu_mask = np.triu_indices(n_layer, k=1)

# ===================== 新增统计计数器 =====================
zero_var_sample_cnt = 0    # 存在至少一层方差近似为0的样本总数
nan_corr_sample_cnt = 0    # 相关矩阵出现NaN/Inf的样本总数

print(f"\n开始逐样本计算层间SCI相关系数...")
for idx in tqdm(range(start_idx, n_sample), desc="SCI计算进度"):
# end_limit = min(start_idx + 1000, n_sample)
# for idx in tqdm(range(start_idx,end_limit),desc="SCI计算进度【前1000样本测试】"):
    delta_single = delta_mm[idx].astype(np.float32)  # shape (33, 1280)

    # 统计：判断是否存在零方差层
    layer_var = np.var(delta_single, axis=1)
    if np.any(layer_var < 1e-15):
        zero_var_sample_cnt += 1

    # 高速批量生成33×33相关矩阵
    corr_mat = np.corrcoef(delta_single)

    # 异常检测：相关矩阵存在NaN/Inf则告警并计数
    if not np.all(np.isfinite(corr_mat)):
        nan_corr_sample_cnt += 1
        print(f"\n[WARN] 样本 {idx} 相关矩阵检测到 NaN / Inf 非法相关系数")

    # 容错清洗NaN/正负无穷（单层方差为0、全零Delta产生非法相关值）
    # 策略：无定义相关系数直接置0，不添加随机噪声，论文解释性更严谨
    corr_mat = np.nan_to_num(corr_mat, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)

    sci_mat_mm[idx] = corr_mat

    triu_vals = corr_mat[triu_mask]

    # ===================== SCI取值逻辑说明 =====================
    # 方案1【当前默认，与前期PTEN/TEM1论文基线完全一致】：原始相关系数，区分正负相关
    sci_mean_mm[idx] = np.mean(triu_vals)
    top20_vals = np.partition(triu_vals, -20)[-20:]
    sci_top20_mm[idx] = np.mean(top20_vals)
    top50_vals = np.partition(triu_vals, -50)[-50:]
    sci_top50_mm[idx] = np.mean(top50_vals)

    # 方案2【如需表征层间相关强度/结构冲突幅度，取消下方注释启用绝对值模式】
    # abs_vals = np.abs(triu_vals)
    # sci_mean_mm[idx] = np.mean(abs_vals)
    # top20_vals = np.partition(abs_vals, -20)[-20:]
    # sci_top20_mm[idx] = np.mean(top20_vals)
    # top50_vals = np.partition(abs_vals, -50)[-50:]
    # sci_top50_mm[idx] = np.mean(top50_vals)
    # ==========================================================

    # ===================== 核心时序保护：先全部落盘，再写入断点，杜绝断点超前数据 =====================
    if SAVE_CHECKPOINT and ((idx + 1) % CHECKPOINT_INTERVAL == 0):
        # 强制刷新所有memmap缓冲区写入磁盘
        sci_mat_mm.flush()
        sci_mean_mm.flush()
        sci_top20_mm.flush()
        sci_top50_mm.flush()
        # 确认数据落地后，再更新断点下标
        ckpt_data = {"processed_sample_idx": idx + 1}
        with open(CHECKPOINT_PATH, "w", encoding="utf-8") as f:
            json.dump(ckpt_data, f, indent=2)

# 全部运算结束，最终全局落盘
sci_mat_mm.flush()
sci_mean_mm.flush()
sci_top20_mm.flush()
sci_top50_mm.flush()

# 导出传统npy文件，完全兼容原有P2~P5旧代码，无需修改下游逻辑
np.save(config.SCI_SCORES_MEAN, np.array(sci_mean_mm))
np.save(config.SCI_SCORES_TOP20, np.array(sci_top20_mm))
np.save(config.SCI_SCORES_TOP50, np.array(sci_top50_mm))
print(f"\n✅ 已导出传统.npy得分文件，原有P2/P3/P4/P5脚本可直接np.load读取，无需改动代码")

# 运算完成自动清理断点，防止下次启动误跳过全部样本
if SAVE_CHECKPOINT and os.path.exists(CHECKPOINT_PATH):
    os.remove(CHECKPOINT_PATH)
    print(f"✅ 全量计算完成，自动清理断点文件: {CHECKPOINT_PATH}")

# 释放内存映射句柄
del delta_mm, sci_mat_mm, sci_mean_mm, sci_top20_mm, sci_top50_mm

# ===================== 收尾输出统计信息（可直接写入论文描述） =====================
print(f"\n========== 全局异常统计汇总 ==========")
print(f"总样本数量：{n_sample:,}")
print(f"存在至少一层零方差的样本数：{zero_var_sample_cnt:,}，占比：{zero_var_sample_cnt / n_sample:.4%}")
print(f"相关矩阵出现NaN/Inf异常的样本数：{nan_corr_sample_cnt:,}，占比：{nan_corr_sample_cnt / n_sample:.4%}")
print(f"======================================")

print(f"\n✅ P1 SCI层间相关计算全部完成（95142样本五蛋白正式封版，不再迭代修改）")
print("输出文件清单：")
print(f"  SCI 33×33完整矩阵Memmap: {sci_mat_out_path}")
print(f"  SCI均值得分Memmap: {sci_mean_path}")
print(f"  SCI Top20得分Memmap: {sci_top20_path}")
print(f"  SCI Top50得分Memmap: {sci_top50_path}")
print(f"  兼容旧脚本 SCI均值.npy: {config.SCI_SCORES_MEAN}")
print(f"  兼容旧脚本 SCI Top20.npy: {config.SCI_SCORES_TOP20}")
print(f"  兼容旧脚本 SCI Top50.npy: {config.SCI_SCORES_TOP50}")
print("\n下游读取示例：")
print('import numpy as np')
print(f'n = {n_sample}')
print(f'# Memmap低内存读取大矩阵')
print(f'sci_mat = np.memmap("{sci_mat_out_path}", dtype=np.float32, mode="r", shape=(n,{n_layer},{n_layer}))')
print(f'# 原有旧脚本直接使用np.load，无需修改')
print(f'sci_mean = np.load("{config.SCI_SCORES_MEAN}")')
print("=" * 80)