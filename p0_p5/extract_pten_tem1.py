import torch
import os

# ========== 关键：强制模型从项目目录加载 ==========
PROJECT_DIR = r"/mnt/sda/gws_1020251255/model/hub"
os.makedirs(PROJECT_DIR, exist_ok=True)
os.environ["TORCH_HOME"] = PROJECT_DIR
torch.hub.set_dir(PROJECT_DIR)
# ================================================

import esm
import pandas as pd
import numpy as np
from tqdm import tqdm

# ========== 配置 ==========
DATA_PATH = r"/mnt/sda/gws_1020251255/data/processed/pilot_3000_with_wt.csv"
OUTPUT_DIR = r"/mnt/sda/gws_1020251255/data/processed"
os.makedirs(OUTPUT_DIR, exist_ok=True)

BATCH_SIZE = 4
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"使用设备: {DEVICE}")

# ========== 加载模型 ==========
print("加载 ESM2-650M 模型...")
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
model = model.to(DEVICE)
model.eval()
batch_converter = alphabet.get_batch_converter()

# ========== 读取数据，只保留 PTEN 和 TEM1 ==========
df = pd.read_csv(DATA_PATH)
df = df[df["protein"].isin(["PTEN", "TEM1"])].reset_index(drop=True)
print(f"\n筛选后数据量: {len(df)} 条")
print(df.groupby("protein").size())

# 检查长度是否安全
print(f"\n序列长度检查:")
for protein in df["protein"].unique():
    max_len = df[df["protein"] == protein]["mutated_sequence"].str.len().max()
    print(f"  {protein}: 最大长度 {max_len} (ESM2限制: 1022)")
    if max_len > 1022:
        print(f"  ⚠️ 警告: {protein} 有截断风险!")

# ========== 提取函数 ==========
def extract_layers(sequence, model, alphabet, device):
    data = [("protein", sequence)]
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_tokens = batch_tokens.to(device)
    
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=list(range(1, 34)), return_contacts=False)
    
    seq_len = len(sequence)
    layers = []
    for layer_idx in range(1, 34):
        hidden = results["representations"][layer_idx]
        hidden = hidden[0, 1:1+seq_len, :]  # 去掉 <cls> 和 <eos>
        layers.append(hidden.cpu().numpy())
    
    return np.stack(layers, axis=0)

# ========== 提取 Mut embedding ==========
print(f"\n{'='*50}")
print("提取 Mut 33层 embedding...")
print(f"{'='*50}")

all_mut_emb = []
all_meta = []

n = len(df)
for i in tqdm(range(0, n, BATCH_SIZE), desc="Mut提取"):
    batch_df = df.iloc[i:i+BATCH_SIZE]
    for _, row in batch_df.iterrows():
        seq = row["mutated_sequence"]
        protein = row["protein"]
        mutant = row["mutant"]
        
        emb = extract_layers(seq, model, alphabet, DEVICE)  # (33, seq_len, 1280)
        emb_mean = emb.mean(axis=1)  # (33, 1280)
        
        all_mut_emb.append(emb_mean)
        all_meta.append({
            "protein": protein,
            "mutant": mutant,
            "DMS_score_bin": row["DMS_score_bin"],
            "DMS_score": row["DMS_score"],
        })

mut_embeddings = np.stack(all_mut_emb, axis=0)  # (n, 33, 1280)
print(f"\nMut embedding 形状: {mut_embeddings.shape}")

# 保存
np.save(os.path.join(OUTPUT_DIR, "pten_tem1_mut_embeddings.npy"), mut_embeddings)
pd.DataFrame(all_meta).to_csv(os.path.join(OUTPUT_DIR, "pten_tem1_metadata.csv"), index=False)

# ========== 提取 WT embedding（只有3条唯一序列）==========
print(f"\n{'='*50}")
print("提取 WT 33层 embedding...")
print(f"{'='*50}")

wt_df = df[["protein", "wt_sequence"]].drop_duplicates()
print(wt_df)

wt_embeddings = {}
for _, row in wt_df.iterrows():
    protein = row["protein"]
    seq = row["wt_sequence"]
    print(f"\n  {protein}: 长度 {len(seq)}")
    emb = extract_layers(seq, model, alphabet, DEVICE)
    emb_mean = emb.mean(axis=1)
    wt_embeddings[protein] = emb_mean
    print(f"  → 形状: {emb_mean.shape}")

# 保存 WT
np.save(os.path.join(OUTPUT_DIR, "pten_tem1_wt_embeddings_dict.npy"), wt_embeddings)

# ========== 构建 WT 对齐矩阵（每条样本对应其蛋白的 WT）==========
print(f"\n{'='*50}")
print("构建 Delta = Mut - WT...")
print(f"{'='*50}")

wt_matrix = np.zeros_like(mut_embeddings)  # (n, 33, 1280)
for i, meta in enumerate(all_meta):
    protein = meta["protein"]
    wt_matrix[i] = wt_embeddings[protein]

delta = mut_embeddings - wt_matrix  # (n, 33, 1280)
print(f"Delta 形状: {delta.shape}")

# 保存 Delta
np.save(os.path.join(OUTPUT_DIR, "pten_tem1_delta_embeddings.npy"), delta)

print(f"\n{'='*50}")
print("完成！文件清单:")
print(f"{'='*50}")
print(f"  - pten_tem1_mut_embeddings.npy")
print(f"  - pten_tem1_metadata.csv")
print(f"  - pten_tem1_wt_embeddings_dict.npy")
print(f"  - pten_tem1_delta_embeddings.npy")
print(f"\n下一步: 计算 SCI (33x33 similarity matrix)")