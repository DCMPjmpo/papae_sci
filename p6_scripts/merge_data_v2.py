"""
P6 多蛋白数据合并 (最终版)
"""
import pandas as pd
import numpy as np
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config

print("=" * 70)
print("P6 多蛋白数据合并 (最终版)")
print("=" * 70)

all_data = []
skipped = []

for protein, dataset, filename in config.PROTEIN_DATASETS:
    filepath = os.path.join(config.RAW_DATA_DIR, filename)
    if not os.path.exists(filepath):
        print(f"  [SKIP] 不存在: {filename}")
        continue
    
    df = pd.read_csv(filepath)
    
    # 检查必备列
    required = ["mutant", "mutated_sequence", "DMS_score", "DMS_score_bin"]
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"{filename}: 缺失 {missing}, 可用: {list(df.columns)}")
    
    df["protein"] = protein
    df["dataset"] = dataset
    
    # 提前过滤超长
    max_len = df["mutated_sequence"].str.len().max()
    if max_len > config.ESM_MAX_SEQ_LEN:
        print(f"  [SKIP] {protein} 超长: {max_len} > {config.ESM_MAX_SEQ_LEN}")
        skipped.append((protein, dataset, max_len))
        continue
    
    print(f"  {protein:6s} | {dataset:20s} | {len(df):6d} 条 | len={max_len}")
    all_data.append(df)

if not all_data:
    raise ValueError("没有读取到任何数据!")

merged = pd.concat(all_data, ignore_index=True)
print(f"\n合并: {len(merged)} 条, {merged['protein'].nunique()} 个蛋白")
if skipped:
    print(f"跳过: {skipped}")

# 添加辅助字段
merged["protein_length"] = merged["mutated_sequence"].str.len()
merged["mutation_count"] = merged["mutant"].str.count(r"\:") + 1

# 重复检查
print(f"\n{'='*50}")
print("重复检查:")
for p in sorted(merged["protein"].unique()):
    sub = merged[merged["protein"] == p]
    dup = sub.duplicated(subset=["protein", "mutant"], keep=False).sum()
    print(f"  {p:6s}: {len(sub):6d} 条, 重复 {dup:5d} ({dup/len(sub)*100:.1f}%)")

# 去重聚合（含序列冲突检查 + 平票处理）
print(f"\n{'='*50}")
print("去重聚合:")

def aggregate_group(group):
    # 序列冲突检查：同一个 mutant 必须对应同一个 sequence
    seqs = group["mutated_sequence"].unique()
    if len(seqs) > 1:
        raise ValueError(
            f"序列冲突: {group.name}\n"
            f"  发现 {len(seqs)} 种序列:\n" +
            "\n".join([f"    {i+1}: {s[:60]}..." for i, s in enumerate(seqs[:3])])
        )
    
    # Bin平票处理
    bins = group["DMS_score_bin"]
    votes = Counter(bins)
    
    if len(votes) > 1 and votes[0] == votes[1]:
        median = group["DMS_score"].median()
        # 用全局中位数作为阈值重新分bin
        global_median = group["DMS_score"].median()
        bin_label = 1 if median >= global_median else 0
        print(f"  [平票] {group.name}: 0={votes[0]}, 1={votes[1]} → 用DMS_score中位数={median:.4f} → bin={bin_label}")
    else:
        bin_label = votes.most_common(1)[0][0]
    
    return pd.Series({
        "mutated_sequence": group["mutated_sequence"].iloc[0],
        "DMS_score": group["DMS_score"].median(),
        "DMS_score_bin": bin_label,
        "dataset": ",".join(sorted(group["dataset"].unique())),
        "n_sources": len(group),
        "protein_length": group["protein_length"].iloc[0],
        "mutation_count": group["mutation_count"].iloc[0],
    })

try:
    deduped = merged.groupby(["protein", "mutant"], as_index=False).apply(aggregate_group)
except ValueError as e:
    print(f"\n❌ 聚合失败: {e}")
    sys.exit(1)

deduped["is_aggregated"] = deduped["n_sources"] > 1

print(f"去重: {len(merged)} → {len(deduped)} 条")

# 统计
print(f"\n{'='*50}")
print("各蛋白统计:")
for p in sorted(deduped["protein"].unique()):
    sub = deduped[deduped["protein"] == p]
    n0 = (sub["DMS_score_bin"] == 0).sum()
    n1 = (sub["DMS_score_bin"] == 1).sum()
    n_multi = (sub["mutation_count"] > 1).sum()
    print(f"  {p:6s}: {len(sub):6d} 条 | 0:{n0:4d} | 1:{n1:4d} | 多点:{n_multi:4d}")

# 多点统计
print(f"\n{'='*50}")
print("多点突变统计:")
mut_dist = deduped["mutation_count"].value_counts().sort_index()
print(mut_dist.to_string())

# 保存
deduped.to_csv(config.MERGE_CSV, index=False)
print(f"\n✅ 已保存: {config.MERGE_CSV}")
print(f"   {len(deduped)} 条, {deduped['protein'].nunique()} 个蛋白")

# 最终校验
assert deduped["protein"].nunique() > 1, "至少需要2个蛋白"
assert deduped["mutant"].notna().all(), "存在空mutant"
assert deduped["mutated_sequence"].notna().all(), "存在空sequence"
print("✅ 最终校验通过")