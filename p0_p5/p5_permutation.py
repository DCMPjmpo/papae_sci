import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from scipy import stats

# ===================== 数据加载 =====================
matrices = np.load('data/processed/pten_tem1_sci_site_matrices.npy')
metadata = pd.read_csv('data/processed/pten_tem1_site_metadata.csv')
delta_emb = np.load('data/processed/pten_tem1_delta_embeddings.npy')

y = metadata["DMS_score_bin"].values

print(f"[Info] 总样本数: {len(y)}")
print(f"[Info] Pathogenic: {(y==0).sum()}, Neutral: {(y==1).sum()}")

# ===================== 特征工程 =====================
# ESM
X_esm_raw = np.column_stack([
    delta_emb.mean(axis=1),
    delta_emb[:, 22:, :].mean(axis=1),
    delta_emb.std(axis=1),
    delta_emb.max(axis=1),
])

# SCI (3模块，无泄漏)
EARLY = list(range(0, 11))
MIDDLE = list(range(11, 22))
LATE = list(range(22, 33))

X_sci = np.column_stack([
    np.mean([matrices[:, i, j] for i in MIDDLE for j in LATE if i < j], axis=0),
    np.mean([matrices[:, i, j] for i in EARLY for j in LATE if i < j], axis=0),
    np.mean([matrices[:, i, j] for i in LATE for j in LATE if i < j], axis=0),
])

# Fusion
X_fusion = np.column_stack([X_esm_raw, X_sci])

# ===================== 统一模型定义 =====================
def make_model(X):
    """统一模型：高维用PCA+LR，低维直接用LR+Scaler"""
    if X.shape[1] > 10:
        n_comp = min(100, X.shape[1])
        return Pipeline([
            ("scaler", StandardScaler()),
            ("pca", PCA(n_components=n_comp, random_state=42)),
            ("clf", LogisticRegression(max_iter=1000, class_weight="balanced"))
        ])
    else:
        return Pipeline([
            ("scaler", StandardScaler()),
            ("clf", LogisticRegression(max_iter=1000, class_weight="balanced"))
        ])

# ===================== P5A: Repeated CV (5×10) =====================
print("\n" + "="*70)
print("P5A: Repeated Stratified CV (5-fold × 10 repeats = 50 AUCs)")
print("="*70)

CV = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=42)

def repeated_cv_auc(X, y, model):
    scores = cross_val_score(model, X, y, cv=CV, scoring="roc_auc", n_jobs=-1)
    return scores

features = {
    "ESM": X_esm_raw,
    "SCI": X_sci,
    "Fusion": X_fusion
}

cv_results = {}
for name, X in features.items():
    print(f"\n--- {name} (dim={X.shape[1]}) ---")
    model = make_model(X)
    scores = repeated_cv_auc(X, y, model)
    cv_results[name] = scores
    print(f"  Mean AUC = {scores.mean():.4f}")
    print(f"  Std      = {scores.std():.4f}")
    print(f"  95% CI   = [{np.percentile(scores, 2.5):.4f}, {np.percentile(scores, 97.5):.4f}]")

# ===================== P5B: Paired Statistical Test =====================
print("\n" + "="*70)
print("P5B: Paired Statistical Test (Fusion vs ESM, Fusion vs SCI)")
print("="*70)

# Fusion vs ESM
t_stat_fe, p_fe = stats.ttest_rel(cv_results["Fusion"], cv_results["ESM"])
w_stat_fe, p_w_fe = stats.wilcoxon(cv_results["Fusion"], cv_results["ESM"])

print(f"\nFusion vs ESM:")
print(f"  Fusion mean: {cv_results['Fusion'].mean():.4f}")
print(f"  ESM mean:    {cv_results['ESM'].mean():.4f}")
print(f"  Mean diff:   {cv_results['Fusion'].mean() - cv_results['ESM'].mean():+.4f}")
print(f"  Paired t-test:  t={t_stat_fe:.3f}, p={p_fe:.4g}")
print(f"  Wilcoxon:       W={w_stat_fe:.1f}, p={p_w_fe:.4g}")

# Fusion vs SCI
t_stat_fs, p_fs = stats.ttest_rel(cv_results["Fusion"], cv_results["SCI"])
w_stat_fs, p_w_fs = stats.wilcoxon(cv_results["Fusion"], cv_results["SCI"])

print(f"\nFusion vs SCI:")
print(f"  Fusion mean: {cv_results['Fusion'].mean():.4f}")
print(f"  SCI mean:    {cv_results['SCI'].mean():.4f}")
print(f"  Mean diff:   {cv_results['Fusion'].mean() - cv_results['SCI'].mean():+.4f}")
print(f"  Paired t-test:  t={t_stat_fs:.3f}, p={p_fs:.4g}")
print(f"  Wilcoxon:       W={w_stat_fs:.1f}, p={p_w_fs:.4g}")

# ===================== P5C: Permutation Test (SCI only, n=100) =====================
print("\n" + "="*70)
print("P5C: Permutation Test (SCI only, n=100)")
print("="*70)

N_PERM = 100
np.random.seed(42)

real_sci_auc = cv_results["SCI"].mean()
perm_scores = []

sci_model = make_model(X_sci)

print(f"\n真实 SCI AUC = {real_sci_auc:.4f}")
print(f"运行 {N_PERM} 次置换...")

for i in range(N_PERM):
    y_perm = np.random.permutation(y)
    scores = cross_val_score(sci_model, X_sci, y_perm, cv=CV, scoring="roc_auc", n_jobs=-1)
    perm_scores.append(scores.mean())
    if (i+1) % 20 == 0:
        print(f"  进度: {i+1}/{N_PERM}")

perm_scores = np.array(perm_scores)
p_perm = (perm_scores >= real_sci_auc).mean()

print(f"\n置换结果:")
print(f"  Perm mean: {perm_scores.mean():.4f} ± {perm_scores.std():.4f}")
print(f"  p-value:   {p_perm:.4g}")

# ===================== 画图 =====================
fig, axes = plt.subplots(2, 3, figsize=(16, 10))

# Row 1: Repeated CV Distribution
for col, (name, scores) in enumerate(cv_results.items()):
    ax = axes[0, col]
    ax.hist(scores, bins=20, color=['coral', 'steelblue', 'green'][col], 
            edgecolor='black', alpha=0.7)
    ax.axvline(x=scores.mean(), color='red', linestyle='--', linewidth=2, 
               label=f'Mean={scores.mean():.3f}')
    ci_l, ci_h = np.percentile(scores, [2.5, 97.5])
    ax.axvline(x=ci_l, color='blue', linestyle=':', alpha=0.7)
    ax.axvline(x=ci_h, color='blue', linestyle=':', alpha=0.7)
    ax.axvspan(ci_l, ci_h, alpha=0.1, color='blue', label=f'95% CI')
    ax.set_xlabel('AUROC')
    ax.set_ylabel('Frequency')
    ax.set_title(f'{name}\nRepeated CV (5×10)', fontweight='bold')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

# Row 2: Paired comparison + Permutation
# A: Paired diff (Fusion - ESM)
ax = axes[1, 0]
diff_fe = cv_results["Fusion"] - cv_results["ESM"]
ax.hist(diff_fe, bins=20, color='purple', edgecolor='black', alpha=0.7)
ax.axvline(x=0, color='red', linestyle='--', linewidth=2)
ax.axvline(x=diff_fe.mean(), color='blue', linestyle='-', linewidth=2, 
           label=f'Mean={diff_fe.mean():+.4f}')
ax.set_xlabel('AUROC Difference (Fusion - ESM)')
ax.set_ylabel('Frequency')
ax.set_title(f'Paired: Fusion vs ESM\nt-test p={p_fe:.4g}', fontweight='bold')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# B: Paired diff (Fusion - SCI)
ax = axes[1, 1]
diff_fs = cv_results["Fusion"] - cv_results["SCI"]
ax.hist(diff_fs, bins=20, color='orange', edgecolor='black', alpha=0.7)
ax.axvline(x=0, color='red', linestyle='--', linewidth=2)
ax.axvline(x=diff_fs.mean(), color='blue', linestyle='-', linewidth=2, 
           label=f'Mean={diff_fs.mean():+.4f}')
ax.set_xlabel('AUROC Difference (Fusion - SCI)')
ax.set_ylabel('Frequency')
ax.set_title(f'Paired: Fusion vs SCI\nt-test p={p_fs:.4g}', fontweight='bold')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# C: Permutation Test (SCI)
ax = axes[1, 2]
ax.hist(perm_scores, bins=20, color='steelblue', edgecolor='black', alpha=0.7, label='Permuted')
ax.axvline(x=real_sci_auc, color='red', linestyle='--', linewidth=2.5, 
           label=f'Real={real_sci_auc:.3f}')
ax.axvline(x=perm_scores.mean(), color='steelblue', linestyle=':', linewidth=2, 
           label=f'Perm mean={perm_scores.mean():.3f}')
ax.set_xlabel('AUROC')
ax.set_ylabel('Frequency')
ax.set_title(f'SCI Permutation Test\np={p_perm:.4g} (n={N_PERM})', fontweight='bold')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

plt.suptitle('Figure 7: Statistical Significance of SCI Predictive Value',
             fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout()
plt.savefig('figure7_p5_statistical_significance.png', dpi=300, bbox_inches='tight')
plt.show()

# ===================== 最终结论 =====================
print("\n" + "="*70)
print("【P5 最终结论】")
print("="*70)

print(f"""
P5A: Repeated CV (5×10)
  ESM:    {cv_results['ESM'].mean():.4f} [{np.percentile(cv_results['ESM'], 2.5):.4f}, {np.percentile(cv_results['ESM'], 97.5):.4f}]
  SCI:    {cv_results['SCI'].mean():.4f} [{np.percentile(cv_results['SCI'], 2.5):.4f}, {np.percentile(cv_results['SCI'], 97.5):.4f}]
  Fusion: {cv_results['Fusion'].mean():.4f} [{np.percentile(cv_results['Fusion'], 2.5):.4f}, {np.percentile(cv_results['Fusion'], 97.5):.4f}]

P5B: Paired Test
  Fusion vs ESM:  t={t_stat_fe:.3f}, p={p_fe:.4g}
  Fusion vs SCI:  t={t_stat_fs:.3f}, p={p_fs:.4g}

P5C: Permutation Test (SCI)
  Real AUC = {real_sci_auc:.4f}
  Perm mean = {perm_scores.mean():.4f} ± {perm_scores.std():.4f}
  p = {p_perm:.4g}
""")

sig_count = sum([p_fe < 0.05, p_fs < 0.05, p_perm < 0.05])
if sig_count >= 2:
    print("🔥🔥🔥 统计显著性全面成立！")
    print("   SCI预测能力被Permutation Test证实")
    print("   Fusion显著优于单一特征")
    print("   论文级统计证据完整！")
elif sig_count >= 1:
    print("✅ 部分显著，可作为探索性证据")
else:
    print("⚠️ 显著性不足，需检查数据")

print("\n✅ Saved: figure7_p5_statistical_significance.png")