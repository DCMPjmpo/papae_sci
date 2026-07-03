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
DATA_PATH = r"/mnt/sda/gws_1020251255/data/processed/pilot_3000.csv"
OUTPUT_DIR = r"/mnt/sda/gws_1020251255/data/processed"
os.makedirs(OUTPUT_DIR, exist_ok=True)

BATCH_SIZE = 4  # CPU 用 4，有 GPU 可以调大
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

print(f"使用设备: {DEVICE}")
print(f"模型加载目录: {torch.hub.get_dir()}")

# ========== 加载模型（从本地）==========
print("\n加载 ESM2-650M 模型...")

model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
model = model.to(DEVICE)
model.eval()

batch_converter = alphabet.get_batch_converter()

# ========== 读取数据 ==========
df = pd.read_csv(DATA_PATH)
print(f"\n数据量: {len(df)} 条")
print(df.groupby("protein").size())

# ========== 提取函数 ==========
def extract_layers(sequence, model, alphabet, device):
    """
    提取 ESM2 全部33层的 hidden states
    返回: (33, seq_len, 1280) 的 numpy 数组
    """
    data = [("protein", sequence)]
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_tokens = batch_tokens.to(device)
    
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=list(range(1, 34)), return_contacts=False)
    
    seq_len = len(sequence)
    layers = []
    for layer_idx in range(1, 34):
        hidden = results["representations"][layer_idx]
        hidden = hidden[0, 1:1+seq_len, :]
        layers.append(hidden.cpu().numpy())
    
    return np.stack(layers, axis=0)

# ========== 批量提取 ==========
print(f"\n开始提取 33 层 embedding...")
print(f"Batch size: {BATCH_SIZE}")

all_embeddings = []
all_metadata = []

n = len(df)
for i in tqdm(range(0, n, BATCH_SIZE), desc="提取进度"):
    batch_df = df.iloc[i:i+BATCH_SIZE]
    
    for _, row in batch_df.iterrows():
        seq = row["mutated_sequence"]
        protein = row["protein"]
        mutant = row["mutant"]
        label = row["DMS_score_bin"]
        
        try:
            emb = extract_layers(seq, model, alphabet, DEVICE)
            emb_mean = emb.mean(axis=1)  # (33, 1280)
            
            all_embeddings.append(emb_mean)
            all_metadata.append({
                "protein": protein,
                "mutant": mutant,
                "DMS_score_bin": label,
                "seq_len": len(seq),
                "DMS_score": row["DMS_score"],
            })
            
        except Exception as e:
            print(f"\n❌ 错误: {protein} {mutant} - {e}")
            continue

# ========== 保存 ==========
print(f"\n提取完成: {len(all_embeddings)} 条")

embeddings = np.stack(all_embeddings, axis=0)
print(f"Embedding 形状: {embeddings.shape}")

emb_path = os.path.join(OUTPUT_DIR, "esm2_33layer_embeddings.npy")
np.save(emb_path, embeddings)
print(f"已保存: {emb_path}")

meta_df = pd.DataFrame(all_metadata)
meta_path = os.path.join(OUTPUT_DIR, "esm2_metadata.csv")
meta_df.to_csv(meta_path, index=False)
print(f"已保存: {meta_path}")

print(f"\n{'='*50}")
print("完成！下一步: 计算 SCI (Semantic Conflict Index)")
print(f"{'='*50}")