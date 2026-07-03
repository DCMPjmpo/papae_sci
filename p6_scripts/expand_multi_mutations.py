"""
P6 多点突变展开 (最终版)
输入: all_proteins_with_wt.csv
输出: all_proteins_expanded.csv
新增: mutation_position 列 (提取突变位点位置)
"""
import pandas as pd
import re
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config

print("=" * 70)
print("P6 多点突变展开 (最终版)")
print("=" * 70)

df = pd.read_csv(config.MERGE_WITH_WT)
print(f"输入: {len(df)} 条")

MUTATION_PATTERN = re.compile(r'([A-Z])(\d+)([A-Z\\*])')
records = []

for _, row in df.iterrows():
    mutant = row["mutant"]
    parts = mutant.split(config.MULTI_SEPARATOR)

    if len(parts) == 1:
        # 单点
        rec = row.to_dict()
        rec["parent_mutant"] = mutant
        rec["mutation_index"] = 1
        rec["n_mutations"] = 1
        records.append(rec)
    else:
        # 多点展开
        for idx, part in enumerate(parts, 1):
            part = part.strip()
            m = MUTATION_PATTERN.match(part)
            if not m:
                print(f"  [WARN] 跳过: {mutant} 中的 {part}")
                continue

            rec = row.to_dict()
            rec["mutant"] = part
            rec["parent_mutant"] = mutant
            rec["mutation_index"] = idx
            rec["n_mutations"] = len(parts)
            records.append(rec)

expanded = pd.DataFrame(records)
print(f"展开: {len(df)} → {len(expanded)} 条")

# ===== 新增: mutation_position =====
print("\n提取 mutation_position...")
# 从 mutant 列提取数字位置，例如 "A123V" → 123
pos_match = expanded["mutant"].str.extract(r"[A-Z](\d+)")
expanded["mutation_position"] = pos_match.astype(int)

print(f"mutation_position 范围: {expanded['mutation_position'].min()} - {expanded['mutation_position'].max()}")

# 统计
print(f"\n{'='*50}")
print("展开统计:")
for p in sorted(expanded["protein"].unique()):
    sub = expanded[expanded["protein"] == p]
    orig = df[df["protein"] == p]
    n_multi = (orig["mutation_count"] > 1).sum() if "mutation_count" in orig.columns else 0
    print(f"  {p:6s}: {len(orig):5d} → {len(sub):5d} (多点:{n_multi})")

# 保存
expanded.to_csv(config.MERGE_EXPANDED, index=False)
print(f"\n✅ 保存: {config.MERGE_EXPANDED}")

# 校验
assert expanded["mutant"].notna().all(), "mutant 存在空值"
assert expanded["wt_sequence"].notna().all(), "wt_sequence 存在空值"
assert expanded["mutated_sequence"].notna().all(), "mutated_sequence 存在空值"
assert expanded["mutation_position"].notna().all(), "mutation_position 存在空值"
assert (expanded["mutation_position"] > 0).all(), "mutation_position 必须 > 0"
print("✅ 校验通过")