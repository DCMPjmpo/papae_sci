import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
import os

# ========== 路径 ==========
OUTPUT_DIR = r"/mnt/sda/gws_1020251255/data/processed"
PLOT_DIR = r"/mnt/sda/gws_1020251255/plots"
os.makedirs(PLOT_DIR, exist_ok=True)

# ========== 读取数据 ==========
delta = np.load(os.path.join(OUTPUT_DIR, "pten_tem1_delta_embeddings.npy"))
meta = pd.read_csv(os.path.join(OUTPUT_DIR, "pten_tem1_metadata.csv"))
labels = meta["DMS_score_bin"].values

path_delta = delta[labels == 0]
neut_delta = delta[labels == 1]

# ========== 矩阵运算计算 SCI 矩阵 ==========
def compute_sci_matrices_batch(deltas):
    norms = np.linalg.norm(deltas, axis=2, keepdims=True) + 1e-8
    normalized = deltas / norms
    sim = np.matmul(normalized, normalized.transpose(0, 2, 1))
    sci = 1 - sim
    return sci

print("计算 Pathogenic SCI 矩阵...")
path_sci_matrices = compute_sci_matrices_batch(path_delta)
path_mean = path_sci_matrices.mean(axis=0)

print("计算 Neutral SCI 矩阵...")
neut_sci_matrices = compute_sci_matrices_batch(neut_delta)
neut_mean = neut_sci_matrices.mean(axis=0)

diff = path_mean - neut_mean

# ========== 计算原始 p-value + FDR 校正 ==========
print("计算 P-value 矩阵...")
n_layers = 33

# 收集上三角的 p-value
upper_ps = []
upper_indices = []

for i in range(n_layers):
    for j in range(i+1, n_layers):
        path_values = path_sci_matrices[:, i, j]
        neut_values = neut_sci_matrices[:, i, j]
        _, p = stats.mannwhitneyu(path_values, neut_values, alternative='greater')
        upper_ps.append(p)
        upper_indices.append((i, j))

# FDR 校正 (Benjamini-Hochberg)
_, p_adj, _, _ = multipletests(upper_ps, method="fdr_bh")

# 重构校正后的 p-value 矩阵
pvalue_matrix = np.ones((n_layers, n_layers))
pvalue_adj_matrix = np.ones((n_layers, n_layers))

for (i, j), p_raw, p_fdr in zip(upper_indices, upper_ps, p_adj):
    pvalue_matrix[i, j] = p_raw
    pvalue_matrix[j, i] = p_raw
    pvalue_adj_matrix[i, j] = p_fdr
    pvalue_adj_matrix[j, i] = p_fdr

# ========== 画图 ==========
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Pathogenic
sns.heatmap(path_mean, ax=axes[0,0], cmap='RdYlBu_r', vmin=0, vmax=1, 
            square=True, cbar_kws={'label': 'SCI'})
axes[0,0].set_title("Pathogenic: Mean SCI Matrix")

# Neutral
sns.heatmap(neut_mean, ax=axes[0,1], cmap='RdYlBu_r', vmin=0, vmax=1,
            square=True, cbar_kws={'label': 'SCI'})
axes[0,1].set_title("Neutral: Mean SCI Matrix")

# Difference
max_diff = np.abs(diff).max()
sns.heatmap(diff, ax=axes[1,0], cmap='RdBu_r', vmin=-max_diff, vmax=max_diff,
            square=True, cbar_kws={'label': 'ΔSCI'})
axes[1,0].set_title("Difference: Pathogenic - Neutral")

# FDR 校正后的显著性
log_fdr = -np.log10(pvalue_adj_matrix + 1e-300)
sns.heatmap(log_fdr, ax=axes[1,1], cmap='Reds', vmin=0, vmax=6,
            square=True, cbar_kws={'label': '-log10(FDR)'})
axes[1,1].set_title("Significance (FDR-corrected): -log10(q)")

# 标注显著区域
for idx, (i, j) in enumerate(upper_indices):
    if p_adj[idx] < 0.001:
        axes[1,1].text(j+0.5, i+0.5, '***', ha='center', va='center', 
                      fontsize=4, color='white')
    elif p_adj[idx] < 0.01:
        axes[1,1].text(j+0.5, i+0.5, '**', ha='center', va='center', 
                      fontsize=4, color='white')
    elif p_adj[idx] < 0.05:
        axes[1,1].text(j+0.5, i+0.5, '*', ha='center', va='center', 
                      fontsize=4, color='black')

for ax in axes.flat:
    ax.set_xlabel("Layer")
    ax.set_ylabel("Layer")

plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, "sci_heatmap.png"), dpi=300, bbox_inches="tight")
plt.show()

print(f"\n已保存: {os.path.join(PLOT_DIR, 'sci_heatmap.png')}")

# ========== 差异最大的显著层对 ==========
print(f"\n{'='*50}")
print("差异最大的显著层对 (Top 20, FDR-corrected):")
print(f"{'='*50}")

# 按校正后的 p-value 排序
sorted_idx = np.argsort(p_adj)
for rank in range(min(20, len(sorted_idx))):
    idx = sorted_idx[rank]
    i, j = upper_indices[idx]
    q = p_adj[idx]
    d = diff[i, j]
    print(f"  Rank {rank+1}: Layer {i+1} vs Layer {j+1}: Δ = {d:.6f}, q = {q:.6f} {'***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else ''}")

# 统计显著层对数量
sig_raw = sum(1 for p in upper_ps if p < 0.05)
sig_fdr = sum(1 for q in p_adj if q < 0.05)
total_pairs = len(upper_ps)

print(f"\n显著层对统计:")
print(f"  Raw p < 0.05:  {sig_raw} / {total_pairs} ({sig_raw/total_pairs*100:.1f}%)")
print(f"  FDR q < 0.05:  {sig_fdr} / {total_pairs} ({sig_fdr/total_pairs*100:.1f}%)")