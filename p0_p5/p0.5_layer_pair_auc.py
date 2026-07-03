import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns

# 路径
matrices = np.load('data/p0_5processed/pten_tem1_sci_site_matrices.npy')
metadata = pd.read_csv('data//p0_5processed/pten_tem1_site_metadata.csv')
labels = metadata["DMS_score_bin"].values

y = (labels == 0).astype(int)
path_mask = labels == 0
neut_mask = labels == 1

print("="*60)
print("528 层对 AUC 扫描 + Permutation + Bootstrap")
print("="*60)

results = []
triu_idx = np.triu_indices(33, k=1)

# 真实 AUC 扫描
for idx in range(len(triu_idx[0])):
    i, j = triu_idx[0][idx], triu_idx[1][idx]
    sci_ij = matrices[:, i, j]
    try:
        auc = roc_auc_score(y, sci_ij)
    except:
        auc = 0.5
    diff_pn = matrices[path_mask, i, j].mean() - matrices[neut_mask, i, j].mean()
    results.append({'Layer_i': i+1, 'Layer_j': j+1, 'AUC': auc, 'Diff_P_N': diff_pn})

results_df = pd.DataFrame(results)
results_df['AUC_abs'] = np.abs(results_df['AUC'] - 0.5)
results_df = results_df.sort_values(['AUC_abs', 'Diff_P_N'], ascending=[False, False])

# ① Permutation Test（修正 p-value）
print("\n" + "="*60)
print("Permutation Test (n=200, 先确认趋势)")
print("="*60)

n_perm = 200
perm_best_aucs = []

for p in range(n_perm):
    y_perm = np.random.permutation(y)
    best_auc_perm = 0.5
    
    for idx in range(len(triu_idx[0])):
        i, j = triu_idx[0][idx], triu_idx[1][idx]
        sci_ij = matrices[:, i, j]
        try:
            auc_p = roc_auc_score(y_perm, sci_ij)
        except:
            auc_p = 0.5
        auc_p_flip = max(auc_p, 1 - auc_p)
        best_auc_perm = max(best_auc_perm, auc_p_flip)
    
    perm_best_aucs.append(best_auc_perm)
    
    if (p + 1) % 50 == 0:
        print(f"  {p+1}/{n_perm} done")

perm_best_aucs = np.array(perm_best_aucs)

# 真实最佳
best_real = results_df.iloc[0]
best_auc_raw = best_real['AUC']
best_auc_flip = max(best_auc_raw, 1 - best_auc_raw)

# 修正 p-value: (count + 1) / (n + 1)
p_value = (np.sum(perm_best_aucs >= best_auc_flip) + 1) / (n_perm + 1)

print(f"\n真实最佳 FlipAUC: {best_auc_flip:.4f}")
print(f"Permutation 均值: {perm_best_aucs.mean():.4f} ± {perm_best_aucs.std():.4f}")
print(f"Permutation 95% 分位: {np.percentile(perm_best_aucs, 95):.4f}")
print(f"p-value: {p_value:.4f}")

if p_value < 0.001:
    print("🔥🔥🔥 极显著！真实机制信号")
elif p_value < 0.01:
    print("🔥 显著")
elif p_value < 0.05:
    print("✅ 边缘显著")
else:
    print("⚠️ 不显著")

# ② Bootstrap（修正：Top10 而不是 Top1）
print("\n" + "="*60)
print("Bootstrap 稳定性 (n=50, Top10 层对)")
print("="*60)

n_boot = 50
top10_boot_pairs = []

for b in range(n_boot):
    idx_boot = np.random.choice(len(y), size=len(y), replace=True)
    y_boot = y[idx_boot]
    matrices_boot = matrices[idx_boot]
    
    boot_results = []
    for idx in range(len(triu_idx[0])):
        i, j = triu_idx[0][idx], triu_idx[1][idx]
        sci_b = matrices_boot[:, i, j]
        try:
            auc_b = roc_auc_score(y_boot, sci_b)
        except:
            auc_b = 0.5
        auc_b_flip = max(auc_b, 1 - auc_b)
        boot_results.append({
            'Layer_i': i+1, 'Layer_j': j+1,
            'AUC': auc_b, 'FlipAUC': auc_b_flip
        })
    
    boot_df = pd.DataFrame(boot_results)
    boot_df['AUC_abs'] = np.abs(boot_df['AUC'] - 0.5)
    boot_df = boot_df.sort_values('AUC_abs', ascending=False)
    
    # 保存 Top10
    for _, row in boot_df.head(10).iterrows():
        pair = (int(row['Layer_i']), int(row['Layer_j']))
        top10_boot_pairs.append(pair)
    
    if (b + 1) % 10 == 0:
        print(f"  {b+1}/{n_boot} done")

pair_counter = Counter(top10_boot_pairs)
print(f"\nBootstrap Top10 稳定性:")
for pair, count in pair_counter.most_common(10):
    pct = count / n_boot * 100
    bar = "█" * int(pct / 5)
    print(f"  L{pair[0]}↔L{pair[1]}: {count}/{n_boot} ({pct:.0f}%) {bar}")

# ③ 输出 Top10
print("\n" + "="*60)
print("Top 10 层对")
print("="*60)
print(f"{'Rank':<6}{'Layer i':<10}{'Layer j':<10}{'AUC':<10}{'FlipAUC':<10}{'Δ(P-N)':<10}")
print("-"*60)
for rank, (_, row) in enumerate(results_df.head(10).iterrows(), 1):
    flip = max(row['AUC'], 1-row['AUC'])
    print(f"{rank:<6}{int(row['Layer_i']):<10}{int(row['Layer_j']):<10}"
          f"{row['AUC']:<10.4f}{flip:<10.4f}{row['Diff_P_N']:<10.4f}")

# 保存
results_df.head(20).to_csv('data/layer/layer_pair_top20.csv', index=False)
results_df.head(50).to_csv('data/layer/layer_pair_top50.csv', index=False)
results_df.to_csv('data/layer/layer_pair_auc_ranking.csv', index=False)

# ④ 反向强层对
print("\n" + "="*60)
print("反向强层对")
print("="*60)
reverse = results_df[results_df['AUC'] < 0.5].sort_values('AUC_abs', ascending=False).head(5)
print(f"{'Layer i':<10}{'Layer j':<10}{'AUC':<10}{'FlipAUC':<10}")
print("-"*60)
for _, row in reverse.iterrows():
    print(f"{int(row['Layer_i']):<10}{int(row['Layer_j']):<10}{row['AUC']:<10.4f}{1-row['AUC']:<10.4f}")

# ⑤ Top50 层频率 + 33×33 frequency map
top50 = results_df.head(50)
layer_counter = Counter()
freq_matrix = np.zeros((33, 33))

for _, row in top50.iterrows():
    li, lj = int(row['Layer_i'])-1, int(row['Layer_j'])-1
    layer_counter[int(row['Layer_i'])] += 1
    layer_counter[int(row['Layer_j'])] += 1
    freq_matrix[li, lj] += 1
    freq_matrix[lj, li] += 1

print("\n" + "="*60)
print("Top50 层频率")
print("="*60)
for layer, count in layer_counter.most_common(10):
    print(f"  Layer {layer:<2}: {count} 次")

# 33×33 frequency heatmap
fig, ax = plt.subplots(figsize=(10, 8))
sns.heatmap(freq_matrix, cmap='YlOrRd', square=True, annot=False,
            xticklabels=range(1,34), yticklabels=range(1,34), ax=ax)
ax.set_title('Top50 Layer Pair Frequency Map', fontsize=14, fontweight='bold')
ax.set_xlabel('Layer'); ax.set_ylabel('Layer')
plt.tight_layout()
plt.savefig('data/layer/layer_pair_frequency_map.png', dpi=300, bbox_inches='tight')
plt.show()

# ⑥ 最佳层对分布图（自动翻转）
best_i = int(best_real['Layer_i']) - 1
best_j = int(best_real['Layer_j']) - 1
path_vals = matrices[path_mask, best_i, best_j]
neut_vals = matrices[neut_mask, best_i, best_j]

if best_auc_raw < 0.5:
    path_vals, neut_vals = neut_vals, path_vals

fig, ax = plt.subplots(figsize=(8, 6))
bp = ax.boxplot([neut_vals, path_vals], labels=['Neutral', 'Pathogenic'],
                patch_artist=True, medianprops=dict(color='red', linewidth=2))
bp['boxes'][0].set_facecolor('lightblue')
bp['boxes'][1].set_facecolor('lightcoral')
ax.set_ylabel(f'SCI (L{best_i+1} ↔ L{best_j+1})', fontsize=12)
ax.set_title(f'Best Pair FlipAUC={best_auc_flip:.3f} (p={p_value:.4f})', 
             fontsize=14, fontweight='bold')
ax.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig('data/layer/best_layer_pair_distribution.png', dpi=300, bbox_inches='tight')
plt.show()

# 最终判断
print(f"\n🔥 最佳层对: L{int(best_real['Layer_i'])} ↔ L{int(best_real['Layer_j'])}")
print(f"   FlipAUC: {best_auc_flip:.4f} | p-value: {p_value:.4f}")

if best_auc_flip > 0.75 and p_value < 0.001:
    print("🚀🚀🚀 情况 A: 强机制信号，直接写论文")
elif best_auc_flip > 0.70:
    print("✅ 情况 B: 有信号，建议 Top10 ensemble")
else:
    print("⚠️ 情况 C: 单层对不够，需转向 Graph SCI")

print("\n✅ Done. Upload: layer_pair_mining.png, layer_pair_frequency_map.png, best_layer_pair_distribution.png")