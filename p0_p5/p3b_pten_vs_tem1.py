import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ===================== 数据加载 =====================
matrices = np.load('data/processed/pten_tem1_sci_site_matrices.npy')
metadata = pd.read_csv('data/processed/pten_tem1_site_metadata.csv')

print(f"[Info] 总样本数: {matrices.shape[0]}")
print(f"[Info] Metadata列: {metadata.columns.tolist()}")

# ===================== 手动指定蛋白标识列（可复现）====================
protein_col = "protein"  # <-- 手动指定，确保可复现
# 如果你的列名不同，改成对应的，例如：protein_col = "gene"

if protein_col not in metadata.columns:
    raise ValueError(f"列 '{protein_col}' 不存在！可用列: {metadata.columns.tolist()}")

unique_proteins = metadata[protein_col].dropna().astype(str).unique()
print(f"[Info] 检测到蛋白: {unique_proteins}")

# 分离两个蛋白（自动识别 PTEN/TEM1）
ptene_name = [u for u in unique_proteins if 'PTEN' in u.upper()][0] if any('PTEN' in u.upper() for u in unique_proteins) else unique_proteins[0]
tem1_name = [u for u in unique_proteins if 'TEM' in u.upper()][0] if any('TEM' in u.upper() for u in unique_proteins) else unique_proteins[1]

print(f"[Info] PTEN标识: '{ptene_name}', TEM-1标识: '{tem1_name}'")

pten_mask = metadata[protein_col].astype(str) == ptene_name
tem1_mask = metadata[protein_col].astype(str) == tem1_name

pten_matrices = matrices[pten_mask]
pten_meta = metadata[pten_mask].copy().reset_index(drop=True)
tem1_matrices = matrices[tem1_mask]
tem1_meta = metadata[tem1_mask].copy().reset_index(drop=True)

print(f"[Info] PTEN样本: {pten_matrices.shape[0]}, TEM-1样本: {tem1_matrices.shape[0]}")

# ===================== 模块定义 =====================
EARLY = list(range(0, 11))    # Layer 1-11
MIDDLE = list(range(11, 22))  # Layer 12-22
LATE = list(range(22, 33))    # Layer 23-33

MODULES = {'Early': EARLY, 'Middle': MIDDLE, 'Late': LATE}

MODULE_PAIRS = [
    ('Early', 'Early'), ('Early', 'Middle'), ('Early', 'Late'),
    ('Middle', 'Middle'), ('Middle', 'Late'), ('Late', 'Late')
]

# ===================== 核心函数 =====================
def compute_module_rewiring(mats, meta, label_col='DMS_score_bin',
                            path_label=0, neut_label=1):
    """
    计算单个数据集的模块重组 (Full Mean only)
    返回: results_dict, diff_matrix, ranking_list
    """
    labels = meta[label_col].values
    path_mask = labels == path_label
    neut_mask = labels == neut_label

    n_path = path_mask.sum()
    n_neut = neut_mask.sum()

    if n_path == 0 or n_neut == 0:
        print(f"[Warning] Pathogenic={n_path}, Neutral={n_neut}, 跳过")
        return None, None, []

    path_mean = mats[path_mask].mean(axis=0)
    neut_mean = mats[neut_mask].mean(axis=0)
    diff_matrix = path_mean - neut_mean

    results = {}

    for m1_name, m2_name in MODULE_PAIRS:
        m1_idx = MODULES[m1_name]
        m2_idx = MODULES[m2_name]

        pair_diffs = []
        for i in m1_idx:
            for j in m2_idx:
                if m1_name == m2_name:
                    if i < j:
                        pair_diffs.append(diff_matrix[i, j])
                else:
                    if i < j:
                        pair_diffs.append(diff_matrix[i, j])

        delta = np.mean(pair_diffs) if pair_diffs else 0.0

        key = f"{m1_name}-{m2_name}"
        results[key] = {
            'delta': delta,
            'n_pairs': len(pair_diffs),
        }

    # 排序：哪个模块对最大
    ranking = sorted(results.items(), key=lambda x: x[1]['delta'], reverse=True)

    return results, diff_matrix, ranking


def bootstrap_module(mats, meta, label_col='DMS_score_bin',
                     path_label=0, neut_label=1,
                     n_bootstrap=1000, target_pair='Middle-Late'):
    """
    Bootstrap 验证 Middle-Late 是否显著大于0
    返回: mean, ci_low, ci_high, boot_vals
    """
    labels = meta[label_col].values
    path_idx = np.where(labels == path_label)[0]
    neut_idx = np.where(labels == neut_label)[0]

    if len(path_idx) == 0 or len(neut_idx) == 0:
        return None, None, None, None

    boot_vals = []
    for b in range(n_bootstrap):
        p_sample = np.random.choice(path_idx, size=len(path_idx), replace=True)
        n_sample = np.random.choice(neut_idx, size=len(neut_idx), replace=True)
        sample_idx = np.concatenate([p_sample, n_sample])

        sample_mats = mats[sample_idx]
        sample_meta = pd.DataFrame({label_col: labels[sample_idx]})

        res, _, _ = compute_module_rewiring(
            sample_mats, sample_meta,
            label_col=label_col, path_label=path_label, neut_label=neut_label
        )

        if res and target_pair in res:
            boot_vals.append(res[target_pair]['delta'])

    boot_vals = np.array(boot_vals)
    mu = boot_vals.mean()
    ci_low, ci_high = np.percentile(boot_vals, [2.5, 97.5])
    return mu, ci_low, ci_high, boot_vals


# ===================== 主分析 =====================
print("\n" + "="*70)
print("P3: Independent Protein Validation (PTEN vs TEM-1)")
print("="*70)

# --- 分别计算两个蛋白 ---
print("\n--- PTEN ---")
pten_full, pten_diff, pten_rank = compute_module_rewiring(pten_matrices, pten_meta)

print("\n【PTEN 模块排名 (Full Mean)】")
for rank, (name, info) in enumerate(pten_rank, 1):
    marker = " 🔥" if name == 'Middle-Late' else ""
    print(f"  #{rank}: {name:<15} ΔSCI = {info['delta']:.4f}{marker}")

print("\n--- TEM-1 ---")
tem1_full, tem1_diff, tem1_rank = compute_module_rewiring(tem1_matrices, tem1_meta)

print("\n【TEM-1 模块排名 (Full Mean)】")
for rank, (name, info) in enumerate(tem1_rank, 1):
    marker = " 🔥" if name == 'Middle-Late' else ""
    print(f"  #{rank}: {name:<15} ΔSCI = {info['delta']:.4f}{marker}")

# --- 跨蛋白对比表 ---
print("\n" + "="*70)
print("【跨蛋白对比: 所有模块 ΔSCI (Full Mean)】")
print("="*70)
print(f"{'Module Pair':<15} {'PTEN':<12} {'TEM-1':<12} {'一致?':<8}")
print("-"*70)

for mp in MODULE_PAIRS:
    key = f"{mp[0]}-{mp[1]}"
    p_val = pten_full[key]['delta'] if pten_full else 0
    t_val = tem1_full[key]['delta'] if tem1_full else 0
    consistent = "✅" if (p_val > 0 and t_val > 0) or (p_val < 0 and t_val < 0) else "❌"
    print(f"{key:<15} {p_val:<12.4f} {t_val:<12.4f} {consistent:<8}")

# --- Bootstrap: Middle-Late ---
print("\n" + "="*70)
print("【Bootstrap: Middle-Late ΔSCI > 0 (n=1000)】")
print("="*70)

pten_mu, pten_l, pten_h, pten_boot = bootstrap_module(
    pten_matrices, pten_meta, n_bootstrap=1000
)
tem1_mu, tem1_l, tem1_h, tem1_boot = bootstrap_module(
    tem1_matrices, tem1_meta, n_bootstrap=1000
)

print(f"\nPTEN:")
print(f"  Middle-Late ΔSCI : {pten_mu:.4f}  95% CI [{pten_l:.4f}, {pten_h:.4f}]")
print(f"  {'✅ CI excludes 0' if pten_l > 0 else '❌ CI crosses 0'}")

print(f"\nTEM-1:")
print(f"  Middle-Late ΔSCI : {tem1_mu:.4f}  95% CI [{tem1_l:.4f}, {tem1_h:.4f}]")
print(f"  {'✅ CI excludes 0' if tem1_l > 0 else '❌ CI crosses 0'}")

# --- 最终结论 ---
print("\n" + "="*70)
print("【P3 最终结论】")
print("="*70)

pten_rank_ml = [i for i, (k, _) in enumerate(pten_rank, 1) if k == 'Middle-Late'][0] if pten_rank else -1
tem1_rank_ml = [i for i, (k, _) in enumerate(tem1_rank, 1) if k == 'Middle-Late'][0] if tem1_rank else -1

print(f"\nPTEN  Middle-Late 排名: #{pten_rank_ml}")
print(f"TEM-1 Middle-Late 排名: #{tem1_rank_ml}")

if pten_rank_ml == 1 and tem1_rank_ml == 1 and pten_l > 0 and tem1_l > 0:
    print("\n🔥🔥🔥 独立蛋白验证成功！")
    print("   两个蛋白的 Middle-Late 均为 #1 排名")
    print("   且 Bootstrap 95% CI 均 excludes 0")
    print("   这是论文级跨蛋白泛化证据！")
    conclusion = "STRONG"
elif pten_rank_ml <= 2 and tem1_rank_ml <= 2 and pten_l > 0 and tem1_l > 0:
    print("\n✅ 方向一致，Middle-Late 排名靠前")
    print("   但排名非 #1，建议作为探索性结果")
    conclusion = "MODERATE"
else:
    print("\n⚠️ 跨蛋白复现不完全，需检查数据或考虑蛋白特异性")
    conclusion = "WEAK"

# ===================== 画图: Figure 5 =====================
fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.3)

# --- A: PTEN Module Matrix ---
ax_a = fig.add_subplot(gs[0, 0])
pten_matrix = np.zeros((3, 3))
name_map = {'Early-Early': (0,0), 'Early-Middle': (0,1), 'Early-Late': (0,2),
            'Middle-Middle': (1,1), 'Middle-Late': (1,2), 'Late-Late': (2,2)}
for k, (i, j) in name_map.items():
    pten_matrix[i, j] = pten_full[k]['delta']
    if i != j:
        pten_matrix[j, i] = pten_full[k]['delta']

mask = np.triu(np.ones((3,3), dtype=bool), k=1)
sns.heatmap(pten_matrix, mask=mask, annot=True, fmt='.4f', cmap='RdBu_r',
            center=0, vmin=-0.01, vmax=0.06, square=True, ax=ax_a,
            xticklabels=['Early','Middle','Late'], yticklabels=['Early','Middle','Late'],
            cbar_kws={'shrink': 0.8})
ax_a.set_title('A. PTEN Module Rewiring', fontweight='bold', fontsize=12)

# --- B: TEM-1 Module Matrix ---
ax_b = fig.add_subplot(gs[0, 1])
tem1_matrix = np.zeros((3, 3))
for k, (i, j) in name_map.items():
    tem1_matrix[i, j] = tem1_full[k]['delta']
    if i != j:
        tem1_matrix[j, i] = tem1_full[k]['delta']
sns.heatmap(tem1_matrix, mask=mask, annot=True, fmt='.4f', cmap='RdBu_r',
            center=0, vmin=-0.01, vmax=0.06, square=True, ax=ax_b,
            xticklabels=['Early','Middle','Late'], yticklabels=['Early','Middle','Late'],
            cbar_kws={'shrink': 0.8})
ax_b.set_title('B. TEM-1 Module Rewiring', fontweight='bold', fontsize=12)

# --- C: Ranking Comparison ---
ax_c = fig.add_subplot(gs[0, 2])
rank_data = []
for rank, (name, info) in enumerate(pten_rank, 1):
    rank_data.append({'Protein': 'PTEN', 'Module': name, 'Rank': rank, 'Delta': info['delta']})
for rank, (name, info) in enumerate(tem1_rank, 1):
    rank_data.append({'Protein': 'TEM-1', 'Module': name, 'Rank': rank, 'Delta': info['delta']})

rank_df = pd.DataFrame(rank_data)
rank_pivot = rank_df.pivot(index='Module', columns='Protein', values='Rank')
rank_pivot = rank_pivot.reindex([r[0] for r in pten_rank])

sns.heatmap(rank_pivot, annot=True, fmt='.0f', cmap='RdYlGn_r',
            vmin=1, vmax=6, square=True, ax=ax_c, 
            cbar_kws={'shrink': 0.8, 'label': 'Rank (#1=best)'})
ax_c.set_title('C. Module Ranking', fontweight='bold', fontsize=12)
ax_c.set_xlabel('')

# --- D: Bootstrap PTEN ---
ax_d = fig.add_subplot(gs[1, 0])
ax_d.hist(pten_boot, bins=50, color='steelblue', edgecolor='black', alpha=0.7)
ax_d.axvline(x=0, color='red', linestyle='--', linewidth=2)
ax_d.axvline(x=pten_mu, color='blue', linestyle='-', linewidth=2, label=f'Mean={pten_mu:.4f}')
ax_d.axvline(x=pten_l, color='blue', linestyle=':', alpha=0.7)
ax_d.axvline(x=pten_h, color='blue', linestyle=':', alpha=0.7)
ax_d.set_xlabel('Middle-Late ΔSCI')
ax_d.set_ylabel('Frequency')
ax_d.set_title(f'D. PTEN Bootstrap\nCI=[{pten_l:.4f}, {pten_h:.4f}]', fontweight='bold', fontsize=11)
ax_d.legend(fontsize=8)

# --- E: Bootstrap TEM-1 ---
ax_e = fig.add_subplot(gs[1, 1])
ax_e.hist(tem1_boot, bins=50, color='coral', edgecolor='black', alpha=0.7)
ax_e.axvline(x=0, color='red', linestyle='--', linewidth=2)
ax_e.axvline(x=tem1_mu, color='darkred', linestyle='-', linewidth=2, label=f'Mean={tem1_mu:.4f}')
ax_e.axvline(x=tem1_l, color='darkred', linestyle=':', alpha=0.7)
ax_e.axvline(x=tem1_h, color='darkred', linestyle=':', alpha=0.7)
ax_e.set_xlabel('Middle-Late ΔSCI')
ax_e.set_ylabel('Frequency')
ax_e.set_title(f'E. TEM-1 Bootstrap\nCI=[{tem1_l:.4f}, {tem1_h:.4f}]', fontweight='bold', fontsize=11)
ax_e.legend(fontsize=8)

# --- F: Cross-Protein Barplot ---
ax_f = fig.add_subplot(gs[1, 2])
modules = [r[0] for r in pten_rank]
pten_vals = [pten_full[m]['delta'] for m in modules]
tem1_vals = [tem1_full[m]['delta'] for m in modules]

x = np.arange(len(modules))
width = 0.35
ax_f.bar(x - width/2, pten_vals, width, label='PTEN', color='steelblue', edgecolor='black')
ax_f.bar(x + width/2, tem1_vals, width, label='TEM-1', color='coral', edgecolor='black')

# 标记 Middle-Late
for i, m in enumerate(modules):
    if m == 'Middle-Late':
        ax_f.annotate('★', xy=(x[i], max(pten_vals[i], tem1_vals[i])), 
                     ha='center', fontsize=20, color='red', fontweight='bold')

ax_f.set_ylabel('ΔSCI')
ax_f.set_title('F. Cross-Protein Comparison', fontweight='bold', fontsize=12)
ax_f.set_xticks(x)
ax_f.set_xticklabels(modules, rotation=15, ha='right')
ax_f.legend()
ax_f.axhline(y=0, color='black', linewidth=0.5)
ax_f.grid(True, alpha=0.3, axis='y')

plt.suptitle('Figure 5: Independent Protein Validation of Middle-Late Module Rewiring',
             fontsize=14, fontweight='bold', y=0.98)

plt.savefig('figure5_independent_protein_validation.png', dpi=300, bbox_inches='tight')
plt.show()

print("\n✅ Saved: figure5_independent_protein_validation.png")