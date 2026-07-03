import pandas as pd
import os
from collections import Counter

# ========== 配置路径 ==========
DATA_DIR = r"/mnt/sda/gws_1020251255/data/raw/proteingym/DMS_ProteinGym_substitutions"
OUTPUT_DIR = r"/mnt/sda/gws_1020251255/data/processed"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ========== 文件列表 ==========
files = [
    # BRCA1
    ("BRCA1", "BRCA1_Findlay2018", "BRCA1_HUMAN_Findlay_2018.csv"),
    
    # PTEN
    ("PTEN", "PTEN_Mighell2018", "PTEN_HUMAN_Mighell_2018.csv"),
    ("PTEN", "PTEN_Matreyek2021", "PTEN_HUMAN_Matreyek_2021.csv"),
    
    # TEM-1 (BLAT_ECOLX) — 4个数据集
    ("TEM1", "TEM1_Deng2012", "BLAT_ECOLX_Deng_2012.csv"),
    ("TEM1", "TEM1_Firnberg2014", "BLAT_ECOLX_Firnberg_2014.csv"),
    ("TEM1", "TEM1_Jacquier2013", "BLAT_ECOLX_Jacquier_2013.csv"),
    ("TEM1", "TEM1_Stiffler2015", "BLAT_ECOLX_Stiffler_2015.csv"),
]

# ========== 读取并合并 ==========
all_data = []

for protein, dataset, filename in files:
    filepath = os.path.join(DATA_DIR, filename)
    print(f"正在读取: {filename}")
    
    df = pd.read_csv(filepath)
    df["protein"] = protein
    df["dataset"] = dataset
    
    # 只保留需要的列
    df = df[["protein", "dataset", "mutant", "mutated_sequence", "DMS_score", "DMS_score_bin"]]
    
    all_data.append(df)
    print(f"  → {len(df)} 条记录")

# 合并
merged = pd.concat(all_data, ignore_index=True)
print(f"\n{'='*50}")
print(f"合并后总计: {len(merged)} 条")

# ========== 检查重复（按 protein + mutant）==========
print(f"\n{'='*50}")
print("检查重复情况:")
print(f"{'='*50}")

# 统计各蛋白的重复情况
for protein in merged["protein"].unique():
    sub = merged[merged["protein"] == protein]
    dup_mask = sub.duplicated(subset=["mutant"], keep=False)
    n_dup = dup_mask.sum()
    n_total = len(sub)
    print(f"{protein}: {n_total}条, 重复记录{n_dup}条 ({n_dup/n_total*100:.1f}%)")

# ========== 去重：每个 (protein, mutant) 只保留一条 ==========
print(f"\n{'='*50}")
print("开始去重 (按 protein + mutant 聚合):")
print(f"{'='*50}")

def aggregate_group(group):
    """对一个 (protein, mutant) 分组进行聚合"""
    # DMS_score: 取中位数
    score_median = group["DMS_score"].median()
    
    # DMS_score_bin: 多数投票
    labels = group["DMS_score_bin"].tolist()
    label_counts = Counter(labels)
    majority_label = label_counts.most_common(1)[0][0]  # 出现最多的标签
    
    # mutated_sequence: 取第一个（理论上应该一样）
    seq = group["mutated_sequence"].iloc[0]
    
    # dataset: 记录所有来源，用逗号分隔
    datasets = ",".join(sorted(group["dataset"].unique()))
    
    return pd.Series({
        "mutated_sequence": seq,
        "DMS_score": score_median,
        "DMS_score_bin": majority_label,
        "dataset": datasets,
        "n_sources": len(group),  # 来源数量
    })

# 按 protein + mutant 分组聚合
deduped = merged.groupby(["protein", "mutant"]).apply(aggregate_group).reset_index()

# 添加聚合标记
deduped["is_aggregated"] = deduped["n_sources"] > 1

print(f"去重后: {len(deduped)} 条")
print(f"聚合来源数分布:\n{deduped['n_sources'].value_counts().sort_index()}")

# 检查标签冲突处理情况
print(f"\n{'='*50}")
print("标签冲突处理检查:")
print(f"{'='*50}")
conflict_check = merged.groupby(["protein", "mutant"])["DMS_score_bin"].nunique()
n_conflict = (conflict_check > 1).sum()
print(f"原始标签冲突的突变数: {n_conflict}")
print(f"聚合后全部解决 (取多数投票)")

# ========== 保存去重后的全量数据 ==========
# 保留需要的列
final_cols = ["protein", "dataset", "mutant", "mutated_sequence", 
              "DMS_score", "DMS_score_bin", "n_sources", "is_aggregated"]
deduped = deduped[final_cols]

all_path = os.path.join(OUTPUT_DIR, "all_mutations.csv")
deduped.to_csv(all_path, index=False)
print(f"\n已保存去重全量数据: {all_path}")
print(f"  → {len(deduped)} 条记录")

# ========== 统计 ==========
print(f"\n{'='*50}")
print("各蛋白数据统计:")
print(f"{'='*50}")
print(deduped.groupby(["protein", "dataset"]).size())

print(f"\n{'='*50}")
print("各蛋白标签分布:")
print(f"{'='*50}")
print(deduped.groupby(["protein", "DMS_score_bin"]).size())

# ========== 分层采样 ==========
print(f"\n{'='*50}")
print("分层采样 (每蛋白每标签500条):")
print(f"{'='*50}")

sampled_list = []

for protein in deduped["protein"].unique():
    protein_df = deduped[deduped["protein"] == protein]
    
    for label in [0, 1]:
        label_df = protein_df[protein_df["DMS_score_bin"] == label]
        n_available = len(label_df)
        n_sample = min(500, n_available)
        
        if n_available < 500:
            print(f"  ⚠️ {protein} {label}: 只有{n_available}条，全部取用")
        else:
            print(f"  {protein} {'pathogenic(0)' if label==0 else 'neutral(1)'}: 可用{n_available}条, 采样{n_sample}条")
        
        sampled = label_df.sample(n=n_sample, random_state=42)
        sampled_list.append(sampled)

pilot = pd.concat(sampled_list, ignore_index=True)

# 打乱顺序
pilot = pilot.sample(frac=1, random_state=42).reset_index(drop=True)

print(f"\n采样总计: {len(pilot)} 条")

# 保存
pilot_path = os.path.join(OUTPUT_DIR, "pilot_3000.csv")
pilot.to_csv(pilot_path, index=False)
print(f"\n已保存采样数据: {pilot_path}")

print(f"\n{'='*50}")
print("采样后分布:")
print(f"{'='*50}")
print(pilot.groupby(["protein", "DMS_score_bin"]).size())