"""
P6 多蛋白 WT 恢复 (最终版)
策略: 多候选投票, 警告模式(不严格失败), 逐位点全量验证
"""
import pandas as pd
import re
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config

print("=" * 70)
print("P6 多蛋白 WT 恢复 (最终版)")
print("=" * 70)

df = pd.read_csv(config.MERGE_CSV)
print(f"输入: {len(df)} 条, {df['protein'].nunique()} 个蛋白")

MUTATION_PATTERN = re.compile(r'([A-Z])(\d+)([A-Z\*])')

def recover_wt(mutant_str, mutated_seq, strict=False):
    """
    恢复WT。strict=False时，序列不匹配只警告不失败。
    """
    parts = mutant_str.split(config.MULTI_SEPARATOR)
    wt_seq = list(mutated_seq)
    warnings = 0
    
    for part in parts:
        part = part.strip()
        if not part:
            continue
        
        m = MUTATION_PATTERN.match(part)
        if not m:
            return None, 0
        
        wt_aa, pos_str, mut_aa = m.groups()
        pos = int(pos_str) - 1
        
        if pos >= len(wt_seq):
            return None, 0
        
        # 警告模式：不匹配时记录但不失败
        if wt_seq[pos] != mut_aa:
            warnings += 1
            if strict:
                return None, warnings
        
        wt_seq[pos] = wt_aa
    
    return ''.join(wt_seq), warnings


wt_sequences = {}
protein_stats = {}

for protein in sorted(df["protein"].unique()):
    protein_df = df[df["protein"] == protein]
    print(f"\n{'='*60}")
    print(f"【{protein}】{len(protein_df)} 条")
    
    candidates = []
    warning_count = 0
    fail_count = 0
    
    for _, row in protein_df.iterrows():
        wt_seq, warns = recover_wt(row["mutant"], row["mutated_sequence"], strict=False)
        warning_count += warns
        
        if wt_seq is not None:
            candidates.append(wt_seq)
        else:
            fail_count += 1
    
    print(f"  成功: {len(candidates)} | 失败: {fail_count} | 警告: {warning_count}")
    
    if len(candidates) == 0:
        print(f"  ❌ 无法恢复!")
        protein_stats[protein] = {"recovered": False, "wt_len": 0}
        continue
    
    # 投票选WT
    wt_counter = Counter(candidates)
    wt_seq, votes = wt_counter.most_common(1)[0]
    
    if len(wt_counter) > 1:
        print(f"  ⚠️ {len(wt_counter)} 种WT，取最多({votes}票)")
    else:
        print(f"  ✅ WT一致，验证{votes}次")
    
    wt_sequences[protein] = wt_seq
    protein_stats[protein] = {
        "recovered": True,
        "wt_len": len(wt_seq),
        "warnings": warning_count,
        "n_unique_wt": len(wt_counter),
    }

# 附加WT列
df["wt_sequence"] = df["protein"].map(wt_sequences)

missing = df["wt_sequence"].isna().sum()
if missing > 0:
    print(f"\n❌ {missing} 条缺失WT")
else:
    print(f"\n✅ 全部恢复成功")

# 全量验证
print(f"\n{'='*60}")
print("全量验证:")
total = pass_count = fail_count = 0
fail_examples = []

for _, row in df.iterrows():
    if pd.isna(row["wt_sequence"]):
        continue
    
    parts = row["mutant"].split(config.MULTI_SEPARATOR)
    ok = True
    
    for part in parts:
        part = part.strip()
        m = MUTATION_PATTERN.match(part)
        if not m:
            ok = False
            break
        
        wt_aa, pos_str, mut_aa = m.groups()
        pos = int(pos_str) - 1
        wt_seq = row["wt_sequence"]
        mut_seq = row["mutated_sequence"]
        
        if pos >= len(wt_seq) or pos >= len(mut_seq):
            ok = False
            break
        
        if wt_seq[pos] != wt_aa or mut_seq[pos] != mut_aa:
            ok = False
            if len(fail_examples) < 5:
                fail_examples.append({
                    "protein": row["protein"],
                    "mutant": row["mutant"],
                    "part": part,
                    "pos": pos+1,
                    "wt_exp": wt_aa, "wt_act": wt_seq[pos],
                    "mut_exp": mut_aa, "mut_act": mut_seq[pos],
                })
            break
    
    total += 1
    if ok:
        pass_count += 1
    else:
        fail_count += 1

print(f"  检查: {total} | ✅: {pass_count} ({pass_count/total*100:.1f}%) | ❌: {fail_count}")

if fail_examples:
    print("\n  失败样例:")
    for ex in fail_examples:
        print(f"    {ex['protein']} {ex['mutant']} 位点{ex['part']}")
        print(f"      WT[{ex['pos']}]: {ex['wt_exp']}→{ex['wt_act']}")
        print(f"      MUT[{ex['pos']}]: {ex['mut_exp']}→{ex['mut_act']}")

# 保存
df.to_csv(config.MERGE_WITH_WT, index=False)
print(f"\n✅ 保存: {config.MERGE_WITH_WT}")

# WT序列表
wt_df = pd.DataFrame([
    {"protein": p, "wt_sequence": s, "length": len(s)}
    for p, s in wt_sequences.items()
])
wt_df.to_csv(config.WT_SEQUENCES, index=False)
print(f"✅ 保存WT: {config.WT_SEQUENCES}")
print(wt_df.to_string(index=False))

# 校验
assert len(wt_sequences) > 0, "WT为空"
assert df["wt_sequence"].isna().sum() == 0, "存在缺失WT"
print("✅ 校验通过")