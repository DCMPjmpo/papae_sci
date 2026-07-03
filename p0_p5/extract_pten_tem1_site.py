import torch
import os
import re

# ========== 强制模型从项目目录加载 ==========
PROJECT_DIR = r"/mnt/sda/gws_1020251255/model/hub"
os.makedirs(PROJECT_DIR, exist_ok=True)
os.environ["TORCH_HOME"] = PROJECT_DIR
torch.hub.set_dir(PROJECT_DIR)

import esm
import pandas as pd
import numpy as np
from tqdm import tqdm
from scipy import stats

# ========== 配置 ==========
DATA_PATH = r"/mnt/sda/gws_1020251255/data/processed/pilot_3000_with_wt.csv"
OUTPUT_DIR = r"/mnt/sda/gws_1020251255/data/processed"
os.makedirs(OUTPUT_DIR, exist_ok=True)

DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"使用设备: {DEVICE}")

# ========== 加载模型 ==========
print("加载 ESM2-650M 模型...")
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
model = model.to(DEVICE)
model.eval()
batch_converter = alphabet.get_batch_converter()

# ========== 读取数据 ==========
df = pd.read_csv(DATA_PATH)
df = df[df["protein"].isin(["PTEN", "TEM1"])].reset_index(drop=True)
print(f"数据量: {len(df)} 条")

# ========== 提取 WT 所有位点（每个蛋白只跑一次）==========
def extract_wt_all_positions(wt_seq, used_positions, model, alphabet, device):
    """
    提取 WT 序列**用到过的位点**的 33 层表示
    返回: {position: (33, 1280)}
    """
    data = [("protein", wt_seq)]
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_tokens = batch_tokens.to(device)
    
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=list(range(1, 34)), return_contacts=False)
    
    seq_len = len(wt_seq)
    all_embeddings = {}
    
    for pos in used_positions:
        if pos >= seq_len:
            raise ValueError(f"位点 {pos} 超出序列长度 {seq_len}")
        
        layers = []
        for layer_idx in range(1, 34):
            hidden = results["representations"][layer_idx]
            hidden = hidden[0, 1:1+seq_len, :]
            site_vec = hidden[pos, :]
            layers.append(site_vec.cpu().numpy())
        all_embeddings[pos] = np.stack(layers, axis=0)  # (33, 1280)
    
    return all_embeddings

# ========== 提取 Mut 突变位点 ==========
def extract_mutation_site(sequence, mutant, model, alphabet, device):
    """
    提取突变位点的 33 层表示
    """
    data = [("protein", sequence)]
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_tokens = batch_tokens.to(device)
    
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=list(range(1, 34)), return_contacts=False)
    
    # 解析突变位置
    match = re.match(r'([A-Z])(\d+)([A-Z])', mutant)
    if not match:
        raise ValueError(f"无法解析突变: {mutant}")
    
    wt_aa, pos, mut_aa = match.groups()
    pos = int(pos) - 1  # 0-based
    
    seq_len = len(sequence)
    if pos >= seq_len:
        raise ValueError(f"突变位置 {pos} 超出序列长度 {seq_len}")
    
    # 严格校验：mutated_sequence 在 pos 位置应该是 mut_aa
    actual_aa = sequence[pos]
    if actual_aa != mut_aa:
        raise ValueError(
            f"突变 {mutant} 校验失败: "
            f"期望 mutated_sequence 在位置 {pos+1} 为 '{mut_aa}', "
            f"实际为 '{actual_aa}'"
        )
    
    # 提取突变位点的 33 层表示: (33, 1280)
    layers = []
    for layer_idx in range(1, 34):
        hidden = results["representations"][layer_idx]
        hidden = hidden[0, 1:1+seq_len, :]
        site_vec = hidden[pos, :]
        layers.append(site_vec.cpu().numpy())
    
    return np.stack(layers, axis=0)

# ========== Step 1: 统计每个蛋白用到的突变位点 ==========
print(f"\n{'='*50}")
print("Step 1: 统计每个蛋白用到的突变位点")
print(f"{'='*50}")

used_positions_by_protein = {}
for protein in ["PTEN", "TEM1"]:
    sub = df[df["protein"] == protein]
    positions = set()
    for mutant in sub["mutant"]:
        match = re.match(r'([A-Z])(\d+)([A-Z])', mutant)
        pos = int(match.group(2)) - 1
        positions.add(pos)
    used_positions_by_protein[protein] = sorted(positions)
    print(f"  {protein}: {len(positions)} 个唯一突变位点")

# ========== Step 2: 提取 WT 用到过的位点（每个蛋白只跑一次）==========
print(f"\n{'='*50}")
print("Step 2: 提取 WT 用到过的位点")
print(f"{'='*50}")

wt_df = df[["protein", "wt_sequence"]].drop_duplicates()
wt_all_embeddings = {}

for _, row in wt_df.iterrows():
    protein = row["protein"]
    seq = row["wt_sequence"]
    used_pos = used_positions_by_protein[protein]
    
    print(f"\n  {protein}: 序列长度 {len(seq)}, 提取 {len(used_pos)} 个位点")
    wt_all_embeddings[protein] = extract_wt_all_positions(
        seq, used_pos, model, alphabet, DEVICE
    )
    print(f"  → 已缓存")

# ========== Step 3: 提取 Mut 突变位点 ==========
print(f"\n{'='*50}")
print("Step 3: 提取 Mut 突变位点")
print(f"{'='*50}")

mut_site_emb = []
all_meta = []

for _, row in tqdm(df.iterrows(), total=len(df), desc="Mut-site提取"):
    seq = row["mutated_sequence"]
    protein = row["protein"]
    mutant = row["mutant"]
    
    emb = extract_mutation_site(seq, mutant, model, alphabet, DEVICE)
    mut_site_emb.append(emb)
    
    match = re.match(r'([A-Z])(\d+)([A-Z])', mutant)
    all_meta.append({
        "protein": protein,
        "mutant": mutant,
        "DMS_score_bin": row["DMS_score_bin"],
        "DMS_score": row["DMS_score"],
        "mutation_position": int(match.group(2)) if match else -1,
    })

mut_site_embeddings = np.stack(mut_site_emb, axis=0)  # (2000, 33, 1280)
print(f"\nMut-site embedding: {mut_site_embeddings.shape}")

# ========== Step 4: 构建 WT-site 对齐矩阵 ==========
print(f"\n{'='*50}")
print("Step 4: 构建 WT-site 对齐矩阵")
print(f"{'='*50}")

wt_site_matrix = np.zeros_like(mut_site_embeddings)
for i, meta in enumerate(all_meta):
    protein = meta["protein"]
    pos = meta["mutation_position"] - 1  # 0-based
    wt_site_matrix[i] = wt_all_embeddings[protein][pos]

# ========== Step 5: 计算 Delta ==========
print(f"\n{'='*50}")
print("Step 5: 计算 Delta = Mut-site - WT-site")
print(f"{'='*50}")

delta_site = mut_site_embeddings - wt_site_matrix
print(f"Delta-site: {delta_site.shape}")

np.save(os.path.join(OUTPUT_DIR, "pten_tem1_delta_site.npy"), delta_site)

# ========== Step 6: 计算 SCI（保存完整矩阵 + 多种 score）==========
print(f"\n{'='*50}")
print("Step 6: 计算 SCI（完整矩阵 + mean/top20/top50）")
print(f"{'='*50}")

def compute_sci_full(delta):
    """
    计算完整 SCI 矩阵和多种 SCI score
    delta: (33, 1280)
    返回: sci_matrix (33,33), sci_mean, sci_top20, sci_top50
    """
    n_layers, dim = delta.shape
    sci_matrix = np.zeros((n_layers, n_layers))
    
    for j in range(n_layers):
        for k in range(j+1, n_layers):
            d_j = delta[j]
            d_k = delta[k]
            
            norm_j = np.linalg.norm(d_j)
            norm_k = np.linalg.norm(d_k)
            
            # 跳过近零向量
            if norm_j < 1e-6 or norm_k < 1e-6:
                continue
            
            sim = np.dot(d_j, d_k) / (norm_j * norm_k)
            sci = 1 - sim
            sci_matrix[j, k] = sci
            sci_matrix[k, j] = sci
    
    # 上三角
    upper_tri = sci_matrix[np.triu_indices(n_layers, k=1)]
    upper_tri = upper_tri[upper_tri > 0]  # 去掉跳过产生的0
    
    # 多种 score
    sci_mean = upper_tri.mean()
    sci_top20 = np.sort(upper_tri)[-20:].mean() if len(upper_tri) >= 20 else sci_mean
    sci_top50 = np.sort(upper_tri)[-50:].mean() if len(upper_tri) >= 50 else sci_mean
    
    return sci_matrix, sci_mean, sci_top20, sci_top50

# 批量计算
n_samples = len(delta_site)
sci_matrices = []
sci_scores_mean = []
sci_scores_top20 = []
sci_scores_top50 = []

for i in tqdm(range(n_samples), desc="计算SCI"):
    d = delta_site[i]
    mat, sm, s20, s50 = compute_sci_full(d)
    sci_matrices.append(mat)
    sci_scores_mean.append(sm)
    sci_scores_top20.append(s20)
    sci_scores_top50.append(s50)

sci_matrices = np.array(sci_matrices)  # (2000, 33, 33)
sci_scores_mean = np.array(sci_scores_mean)
sci_scores_top20 = np.array(sci_scores_top20)
sci_scores_top50 = np.array(sci_scores_top50)

print(f"\nSCI matrices: {sci_matrices.shape}")
print(f"SCI mean:     {sci_scores_mean.shape}")
print(f"SCI top20:    {sci_scores_top20.shape}")
print(f"SCI top50:    {sci_scores_top50.shape}")

# 保存
np.save(os.path.join(OUTPUT_DIR, "pten_tem1_sci_site_matrices.npy"), sci_matrices)
np.save(os.path.join(OUTPUT_DIR, "pten_tem1_sci_site_scores_mean.npy"), sci_scores_mean)
np.save(os.path.join(OUTPUT_DIR, "pten_tem1_sci_site_scores_top20.npy"), sci_scores_top20)
np.save(os.path.join(OUTPUT_DIR, "pten_tem1_sci_site_scores_top50.npy"), sci_scores_top50)
pd.DataFrame(all_meta).to_csv(os.path.join(OUTPUT_DIR, "pten_tem1_site_metadata.csv"), index=False)

# ========== Step 7: 统计检验 ==========
print(f"\n{'='*50}")
print("Step 7: 统计检验")
print(f"{'='*50}")

labels = np.array([m["DMS_score_bin"] for m in all_meta])

for score_name, scores in [
    ("SCI-mean", sci_scores_mean),
    ("SCI-top20", sci_scores_top20),
    ("SCI-top50", sci_scores_top50),
]:
    pathogenic = scores[labels == 0]
    neutral = scores[labels == 1]
    
    u, p = stats.mannwhitneyu(pathogenic, neutral, alternative='greater')
    
    # Cliff's Delta
    nx, ny = len(pathogenic), len(neutral)
    u_cd, _ = stats.mannwhitneyu(pathogenic, neutral, alternative="two-sided")
    cliff_d = (2 * u_cd) / (nx * ny) - 1
    
    # AUC
    from sklearn.metrics import roc_auc_score
    auc = roc_auc_score(1 - labels, scores)
    
    print(f"\n{score_name}:")
    print(f"  Pathogenic: {pathogenic.mean():.6f} ± {pathogenic.std():.6f}")
    print(f"  Neutral:    {neutral.mean():.6f} ± {neutral.std():.6f}")
    print(f"  Δ:          {pathogenic.mean() - neutral.mean():.6f}")
    print(f"  p-value:    {p:.6f}")
    print(f"  Cliff's d:  {cliff_d:.4f}")
    print(f"  ROC-AUC:    {auc:.4f}")

print(f"\n{'='*50}")
print("完成！文件清单:")
print(f"{'='*50}")
print(f"  - pten_tem1_delta_site.npy")
print(f"  - pten_tem1_sci_site_matrices.npy      (33×33 matrices)")
print(f"  - pten_tem1_sci_site_scores_mean.npy")
print(f"  - pten_tem1_sci_site_scores_top20.npy")
print(f"  - pten_tem1_sci_site_scores_top50.npy")
print(f"  - pten_tem1_site_metadata.csv")