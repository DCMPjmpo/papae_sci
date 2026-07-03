import pandas as pd
import re
import os

# ========== 路径 ==========
PILOT_PATH = r"/mnt/sda/gws_1020251255/data/processed/pilot_3000.csv"
OUTPUT_DIR = r"/mnt/sda/gws_1020251255/data/processed"

df = pd.read_csv(PILOT_PATH)
print(f"数据量: {len(df)} 条")

# ========== 恢复 WT 序列 ==========
def recover_wt(mutant, mutated_seq):
    """
    mutant 格式: "A123V" (原始氨基酸 + 位置 + 突变氨基酸)
    """
    match = re.match(r'([A-Z])(\d+)([A-Z])', mutant)
    if not match:
        return None
    
    wt_aa, pos, mut_aa = match.groups()
    pos = int(pos) - 1  # 0-based 索引
    
    if pos >= len(mutated_seq):
        return None
    
    # 验证 mutated_sequence 在 pos 位置是否是 mut_aa
    if mutated_seq[pos] != mut_aa:
        # 有些数据格式不一致，尝试不验证直接替换
        pass
    
    # 恢复 WT
    wt_seq = list(mutated_seq)
    wt_seq[pos] = wt_aa
    return ''.join(wt_seq)

# 对每个蛋白，找一个突变恢复 WT（同一蛋白所有突变共享 WT）
wt_sequences = {}

for protein in df["protein"].unique():
    protein_df = df[df["protein"] == protein]
    
    # 尝试恢复 WT，直到成功
    for _, row in protein_df.iterrows():
        wt_seq = recover_wt(row["mutant"], row["mutated_sequence"])
        if wt_seq is not None:
            wt_sequences[protein] = wt_seq
            break
    
    print(f"{protein}: WT 序列长度 = {len(wt_sequences.get(protein, ''))}")

# ========== 为每条记录添加 WT 序列 ==========
df["wt_sequence"] = df["protein"].map(wt_sequences)

# 验证：随机抽查几个突变
print(f"\n{'='*50}")
print("WT 恢复验证 (抽查):")
print(f"{'='*50}")
for _, row in df.sample(3, random_state=42).iterrows():
    print(f"\n蛋白: {row['protein']}")
    print(f"突变: {row['mutant']}")
    print(f"WT  : {row['wt_sequence'][:50]}...")
    print(f"MUT : {row['mutated_sequence'][:50]}...")

# 检查是否有缺失
missing = df["wt_sequence"].isna().sum()
print(f"\n缺失 WT 序列: {missing} 条")

# ========== 保存 ==========
output_path = os.path.join(OUTPUT_DIR, "pilot_3000_with_wt.csv")
df.to_csv(output_path, index=False)
print(f"\n已保存: {output_path}")

# 同时保存 WT 序列字典，方便后面批量提取
wt_df = pd.DataFrame([
    {"protein": p, "wt_sequence": s, "seq_len": len(s)}
    for p, s in wt_sequences.items()
])
wt_path = os.path.join(OUTPUT_DIR, "wt_sequences.csv")
wt_df.to_csv(wt_path, index=False)
print(f"已保存 WT 序列表: {wt_path}")
print(wt_df)

# ========== 深度验证：检查突变位置是否正确 ==========
print(f"\n{'='*50}")
print("深度验证：检查具体突变位置")
print(f"{'='*50}")

for _, row in df.sample(5, random_state=42).iterrows():
    mutant = row["mutant"]
    wt_seq = row["wt_sequence"]
    mut_seq = row["mutated_sequence"]
    
    match = re.match(r'([A-Z])(\d+)([A-Z])', mutant)
    if match:
        wt_aa, pos, mut_aa = match.groups()
        pos = int(pos) - 1
        
        print(f"\n突变: {mutant}")
        print(f"  WT 位置 {pos}: {wt_seq[pos]} (应为 {wt_aa})")
        print(f"  MUT位置 {pos}: {mut_seq[pos]} (应为 {mut_aa})")
        
        # 检查是否只有这一个位置不同
        diff_count = sum(1 for a, b in zip(wt_seq, mut_seq) if a != b)
        print(f"  总差异位置数: {diff_count}")
        if diff_count == 1:
            print("  ✅ 单点突变正确")
        else:
            print(f"  ⚠️ 差异位置数异常（应为1）")