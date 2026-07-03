import pandas as pd
import os

# ========== 路径配置 ==========
DATA_DIR = r"/mnt/sda/gws_1020251255/data/raw/proteingym/DMS_ProteinGym_substitutions"

# ========== 读取四个 TEM-1 文件 ==========
files = [
    ("TEM1_Deng2012", "BLAT_ECOLX_Deng_2012.csv"),
    ("TEM1_Firnberg2014", "BLAT_ECOLX_Firnberg_2014.csv"),
    ("TEM1_Jacquier2013", "BLAT_ECOLX_Jacquier_2013.csv"),
    ("TEM1_Stiffler2015", "BLAT_ECOLX_Stiffler_2015.csv"),
]

all_tem1 = []
for dataset, filename in files:
    filepath = os.path.join(DATA_DIR, filename)
    df = pd.read_csv(filepath)
    df["dataset"] = dataset
    df["protein"] = "TEM1"
    df = df[["protein", "dataset", "mutant", "mutated_sequence", "DMS_score", "DMS_score_bin"]]
    all_tem1.append(df)
    print(f"{dataset}: {len(df)} 条")

merged_tem1 = pd.concat(all_tem1, ignore_index=True)
print(f"\n{'='*50}")
print(f"TEM-1 总计: {len(merged_tem1)} 条")

# ========== 检查重复突变 ==========
print(f"\n{'='*50}")
print("检查重复突变 (按 protein + mutant):")
print(f"{'='*50}")

# 只看 TEM1 内部重复
dup_mask = merged_tem1.duplicated(subset=["protein", "mutant"], keep=False)
dup_df = merged_tem1[dup_mask]

if len(dup_df) == 0:
    print("✅ 无重复突变")
else:
    n_dup_records = len(dup_df)
    n_dup_mutants = dup_df["mutant"].nunique()
    total = len(merged_tem1)
    dup_rate = n_dup_records / total * 100
    
    print(f"重复记录数: {n_dup_records}")
    print(f"重复突变种类数: {n_dup_mutants}")
    print(f"重复率: {dup_rate:.2f}%")
    
    print(f"\n{'='*50}")
    print("重复突变示例 (前20条):")
    print(f"{'='*50}")
    print(dup_df.sort_values("mutant").head(20)[["dataset", "mutant", "DMS_score", "DMS_score_bin"]].to_string())
    
    # 检查同一突变在不同数据集中的标签是否一致
    print(f"\n{'='*50}")
    print("检查同一突变在不同数据集的标签一致性:")
    print(f"{'='*50}")
    
    conflict = []
    for mutant in dup_df["mutant"].unique():
        sub = dup_df[dup_df["mutant"] == mutant]
        labels = sub["DMS_score_bin"].unique()
        if len(labels) > 1:
            conflict.append((mutant, labels.tolist()))
    
    if conflict:
        print(f"⚠️ 标签冲突的突变数: {len(conflict)}")
        print("示例:")
        for m, labels in conflict[:5]:
            print(f"  {m}: 标签 {labels}")
    else:
        print("✅ 所有重复突变的标签一致")

# ========== 建议处理方案 ==========
print(f"\n{'='*50}")
print("处理建议:")
print(f"{'='*50}")

if len(dup_df) == 0:
    print("无需处理，直接合并即可。")
elif dup_rate < 5:
    print(f"重复率 {dup_rate:.2f}% < 5%，建议：直接合并，保留所有记录。")
    print("理由：重复样本在训练/分析中会被自然平均掉，不影响结论。")
else:
    print(f"重复率 {dup_rate:.2f}% >= 5%，建议：去重，每个突变只保留一条。")
    print("去重方法：按 (protein, mutant) 分组，取 DMS_score 中位数。")