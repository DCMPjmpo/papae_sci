import numpy as np
import pandas as pd
from scipy import stats
import os

# ========== 路径 ==========
OUTPUT_DIR = r"/mnt/sda/gws_1020251255/data/processed"

# ========== 读取数据 ==========
delta = np.load(os.path.join(OUTPUT_DIR, "pten_tem1_delta_embeddings.npy"))  # (2000, 33, 1280)
meta = pd.read_csv(os.path.join(OUTPUT_DIR, "pten_tem1_metadata.csv"))
labels = meta["DMS_score_bin"].values

print(f"Delta 形状: {delta.shape}")
print(f"Labels 分布: {np.bincount(labels)}")

# ========== 计算 SCI-score ==========
def compute_sci(delta):
    n_samples, n_layers, dim = delta.shape
    sci_scores = []
    
    for i in range(n_samples):
        d = delta[i]
        pairs = []
        for j in range(n_layers):
            for k in range(j+1, n_layers):
                sim = np.dot(d[j], d[k]) / (np.linalg.norm(d[j]) * np.linalg.norm(d[k]) + 1e-8)
                pairs.append(1 - sim)
        sci_scores.append(np.mean(pairs))
    
    return np.array(sci_scores)

print("\n计算 SCI-score...")
sci = compute_sci(delta)

# ========== 统计检验 ==========
pathogenic = sci[labels == 0]
neutral = sci[labels == 1]

u_stat, p_value = stats.mannwhitneyu(pathogenic, neutral, alternative='greater')

print(f"\n{'='*50}")
print("统计检验: Pathogenic vs Neutral")
print(f"{'='*50}")
print(f"Pathogenic SCI: {pathogenic.mean():.6f} ± {pathogenic.std():.6f}")
print(f"Neutral SCI:    {neutral.mean():.6f} ± {neutral.std():.6f}")
print(f"Mann-Whitney U: p = {p_value:.6f}")
print(f"{'*** 显著!' if p_value < 0.001 else '** 显著!' if p_value < 0.01 else '* 显著!' if p_value < 0.05 else '不显著'}")

# 保存
np.save(os.path.join(OUTPUT_DIR, "pten_tem1_sci_scores.npy"), sci)
print(f"\n已保存: pten_tem1_sci_scores.npy")