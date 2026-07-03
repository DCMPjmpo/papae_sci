import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

matrices = np.load('data/processed/pten_tem1_sci_site_matrices.npy')
metadata = pd.read_csv('data/processed/pten_tem1_site_metadata.csv')

labels = metadata["DMS_score_bin"].values
path_mask = labels == 0
neut_mask = labels == 1

path_mean = matrices[path_mask].mean(axis=0)
neut_mean = matrices[neut_mask].mean(axis=0)
diff_matrix = path_mean - neut_mean

triu_idx = np.triu_indices(33, k=1)
diff_values = diff_matrix[triu_idx]
positive_idx = np.where(diff_values > 0)[0]
sorted_idx = positive_idx[np.argsort(diff_values[positive_idx])[::-1]][:20]

top20_pairs = []
print("="*60)
print("Top 20 正差异层对 (P > N)")
print("="*60)
print(f"{'Rank':<6}{'Layer i':<10}{'Layer j':<10}{'ΔSCI':<12}")
print("-"*60)

for rank, idx in enumerate(sorted_idx, 1):
    i, j = triu_idx[0][idx], triu_idx[1][idx]
    d = diff_values[idx]
    top20_pairs.append((i, j, d))
    print(f"{rank:<6}{i+1:<10}{j+1:<10}{d:<12.4f}")

layer_counter = Counter()
early, middle, late = 0, 0, 0
for i, j, d in top20_pairs:
    for layer in [i+1, j+1]:
        layer_counter[layer] += 1
        if layer <= 11: early += 1
        elif layer <= 22: middle += 1
        else: late += 1

print(f"\n区间分布: Early(1-11)={early}, Middle(12-22)={middle}, Late(23-33)={late}")
if late > early + middle:
    print("🔥 高层主导！低层↔高层冲突机制成立")
elif early > late:
    print("⚠️ 低层更多，机制需调整")
else:
    print("✅ 分布均衡")

np.save('top20_positive_layer_pairs.npy', np.array(top20_pairs, dtype=object))

fig, axes = plt.subplots(2, 2, figsize=(14, 12))
for ax, data, title in [(axes[0,0], path_mean, 'Pathogenic'),
                        (axes[0,1], neut_mean, 'Neutral'),
                        (axes[1,0], diff_matrix, 'Difference (P-N)')]:
    sns.heatmap(data, cmap='RdBu_r', center=0, square=True,
                xticklabels=range(1,34), yticklabels=range(1,34), ax=ax)
    ax.set_title(title, fontweight='bold')

vmax = np.abs(diff_matrix).max()
sns.heatmap(diff_matrix, cmap='RdBu_r', center=0, vmin=-vmax, vmax=vmax, square=True,
            xticklabels=range(1,34), yticklabels=range(1,34), ax=axes[1,1])
for rank, (i, j, d) in enumerate(top20_pairs, 1):
    axes[1,1].add_patch(plt.Rectangle((j, i), 1, 1, fill=False, edgecolor='lime', lw=2))
    axes[1,1].add_patch(plt.Rectangle((i, j), 1, 1, fill=False, edgecolor='lime', lw=2))
    if rank <= 5:
        axes[1,1].text(j+0.5, i+0.5, str(rank), ha='center', va='center',
                       fontsize=8, fontweight='bold', bbox=dict(boxstyle='round', facecolor='white'))
axes[1,1].set_title('Top 20 Positive Pairs', fontweight='bold')
plt.tight_layout()
plt.savefig('data/layer/layer_pair_mining.png', dpi=300, bbox_inches='tight')
plt.show()
print("\n✅ Saved: layer_pair_mining.png")