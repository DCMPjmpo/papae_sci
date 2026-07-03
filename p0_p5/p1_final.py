import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.inspection import permutation_importance

try:
    import xgboost as xgb
    XGB_AVAILABLE = True
except ImportError:
    XGB_AVAILABLE = False
    print("[Warning] xgboost not installed. Run: pip install xgboost")

# ===================== 数据接入区（接入你的P0数据）====================
matrices = np.load('data/processed/pten_tem1_sci_site_matrices.npy')
metadata = pd.read_csv('data/processed/pten_tem1_site_metadata.csv')
y = metadata["DMS_score_bin"].values

print(f"[Info] 样本数: {matrices.shape[0]}, SCI矩阵维度: {matrices.shape[1]}x{matrices.shape[2]}")
# =====================================================

# --------------------- 1. 层对定义（按Bootstrap频率排序）---------------------
top_pairs = [
    (14, 32),  # 96%
    (15, 32),  # 86%
    (5, 32),   # 62%
    (7, 32),   # 62%
    (29, 32),  # 62%
    (10, 32),
    (16, 32),
    (17, 32),
    (10, 31),
    (6, 32),
    (29, 33),
    (14, 33),
    (15, 33),
    (30, 32),
    (31, 32)
]

def build_features(matrices, pairs):
    """pairs 为 1-indexed"""
    return np.column_stack([matrices[:, i-1, j-1] for i, j in pairs])

X_ensemble = build_features(matrices, top_pairs)
X_baseline = matrices[:, 13, 31].reshape(-1, 1)  # L14↔L32

pair_names = [f"L{i}↔L{j}" for i, j in top_pairs]
print(f"[Info] Ensemble维度: {X_ensemble.shape[1]}")

# --------------------- 2. 统一模型定义 ---------------------
def make_lr():
    return Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(max_iter=1000, class_weight="balanced"))
    ])

def make_rf():
    return RandomForestClassifier(
        n_estimators=300, max_depth=5,
        class_weight="balanced", random_state=42
    )

def make_xgb():
    return xgb.XGBClassifier(
        n_estimators=500, max_depth=2, learning_rate=0.03,
        subsample=0.8, colsample_bytree=0.8,
        reg_lambda=3, reg_alpha=1,
        eval_metric="logloss", random_state=42
    )

MODELS = {
    "LR": make_lr,
    "RF": make_rf,
}
if XGB_AVAILABLE:
    MODELS["XGB"] = make_xgb

# --------------------- 3. 评估函数（Repeated CV）---------------------
CV = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=42)

def cv_auc(X, y, model, name="Model"):
    scores = cross_val_score(model, X, y, cv=CV, scoring="roc_auc", n_jobs=-1)
    print(f"{name:38s}: AUC = {scores.mean():.4f} ± {scores.std():.4f}")
    return scores.mean(), scores.std()

# --------------------- 4. 基线：Single Pair L14↔L32 ---------------------
print("\n" + "="*60)
print("=== Baseline: Single Pair L14↔L32 (Repeated CV 5×10) ===")
print("="*60)

baseline_results = {}
for name, fn in MODELS.items():
    mu, std = cv_auc(X_baseline, y, fn(), name=f"{name}  (L14↔L32)")
    baseline_results[name] = (mu, std)

# --------------------- 5. Ensemble：Top15 ---------------------
print("\n" + "="*60)
print("=== P1: Top-15 Layer-Pair Ensemble (Repeated CV 5×10) ===")
print("="*60)

ensemble_results = {}
for name, fn in MODELS.items():
    mu, std = cv_auc(X_ensemble, y, fn(), name=f"{name}  (Top15)")
    ensemble_results[name] = (mu, std)

# --------------------- 6. 公平对比表 ---------------------
print("\n" + "="*60)
print("=== Fair Comparison: Baseline vs Ensemble ===")
print("="*60)
print(f"{'Model':<6} {'Baseline':<16} {'Ensemble':<16} {'Delta':<10}")
print("-"*60)
for name in MODELS.keys():
    b_mu, b_std = baseline_results[name]
    e_mu, e_std = ensemble_results[name]
    print(f"{name:<6} {b_mu:.4f}±{b_std:.4f}     {e_mu:.4f}±{e_std:.4f}     {e_mu-b_mu:+.4f}")
print("="*60)

# --------------------- 7. Permutation Importance ---------------------
if XGB_AVAILABLE:
    print("\n" + "="*60)
    print("=== Permutation Importance (XGBoost, All-Data Fit) ===")
    print("="*60)
    xgb_model = make_xgb()
    xgb_model.fit(X_ensemble, y)

    perm = permutation_importance(
        xgb_model, X_ensemble, y,
        scoring="roc_auc", n_repeats=20, random_state=42, n_jobs=-1
    )

    perm_df = pd.DataFrame({
        "pair": pair_names,
        "importance_mean": perm.importances_mean,
        "importance_std": perm.importances_std
    }).sort_values("importance_mean", ascending=False)
    print(perm_df.to_string(index=False))

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Permutation Importance plot
    ax = axes[0]
    y_pos = np.arange(len(perm_df))
    ax.barh(y_pos, perm_df["importance_mean"][::-1], 
            xerr=perm_df["importance_std"][::-1], color="steelblue", capsize=3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(perm_df["pair"][::-1])
    ax.set_xlabel("Permutation Importance (AUC)")
    ax.set_title("Permutation Importance (XGBoost)")
    ax.axvline(x=0, color="black", linewidth=0.5)

    # XGB native importance (exploratory)
    ax = axes[1]
    native_imp = xgb_model.feature_importances_
    native_df = pd.DataFrame({"pair": pair_names, "importance": native_imp})
    native_df = native_df.sort_values("importance", ascending=False)
    ax.barh(range(len(native_df)), native_df["importance"][::-1], color="coral")
    ax.set_yticks(range(len(native_df)))
    ax.set_yticklabels(native_df["pair"][::-1])
    ax.set_xlabel("Native Feature Importance")
    ax.set_title("Native XGB Importance (Exploratory)")

    plt.tight_layout()
    plt.savefig("p1_permutation_importance.png", dpi=150)
    plt.show()

# --------------------- 8. Feature Correlation Heatmap ---------------------
print("\n" + "="*60)
print("=== Feature Correlation Heatmap ===")
print("="*60)

X_df = pd.DataFrame(X_ensemble, columns=pair_names)
corr = X_df.corr()

plt.figure(figsize=(12, 10))
mask = np.triu(np.ones_like(corr, dtype=bool), k=1)
sns.heatmap(corr, mask=mask, annot=True, fmt=".2f", cmap="RdBu_r",
            center=0, vmin=-1, vmax=1, square=True, linewidths=0.5,
            cbar_kws={"shrink": 0.8})
plt.title("Top-15 Layer-Pair Feature Correlation (Lower Triangle)")
plt.tight_layout()
plt.savefig("p1_feature_correlation.png", dpi=150)
plt.show()

# 打印高相关对
high_corr = []
for i in range(len(corr.columns)):
    for j in range(i+1, len(corr.columns)):
        c = corr.iloc[i, j]
        if abs(c) > 0.8:
            high_corr.append((corr.columns[i], corr.columns[j], c))
high_corr.sort(key=lambda x: abs(x[2]), reverse=True)
print(f"\n高相关对 (|r| > 0.8): {len(high_corr)} 对")
for p1, p2, c in high_corr[:10]:
    print(f"  {p1} ↔ {p2}: r = {c:.3f}")

# --------------------- 9. Progressive K（三种模型分别跑）---------------------
print("\n" + "="*60)
print("=== Progressive K Experiment (Repeated CV 5×10) ===")
print("="*60)

k_list = [1, 3, 5, 7, 10, 15]
prog_results = {name: [] for name in MODELS.keys()}

for k in k_list:
    if k > len(top_pairs):
        break
    X_k = build_features(matrices, top_pairs[:k])
    print(f"\n--- K = {k} ---")
    for name, fn in MODELS.items():
        mu, std = cv_auc(X_k, y, fn(), name=f"{name}  (Top{k})")
        prog_results[name].append((k, mu, std))

# --------------------- 10. K-AUC 曲线（三种模型）---------------------
fig, ax = plt.subplots(figsize=(9, 6))
colors = {"LR": "blue", "RF": "purple", "XGB": "crimson"}
markers = {"LR": "o", "RF": "s", "XGB": "D"}

for name in MODELS.keys():
    ks, mus, stds = zip(*prog_results[name])
    ax.errorbar(ks, mus, yerr=stds, fmt="-" + markers[name], 
                color=colors[name], linewidth=2, markersize=7,
                capsize=4, label=f"{name} Ensemble", alpha=0.9)
    # 画基线
    b_mu, b_std = baseline_results[name]
    ax.axhline(y=b_mu, color=colors[name], linestyle="--", 
               linewidth=1.5, alpha=0.5, label=f"{name} Baseline = {b_mu:.3f}")

ax.axhline(y=0.78, color="orange", linestyle=":", linewidth=2, label="Target AUC = 0.78")
ax.set_xlabel("Number of Top Layer Pairs (K)", fontsize=12)
ax.set_ylabel("Repeated CV AUC (5×10)", fontsize=12)
ax.set_title("P1: Ensemble Performance vs. K (All Models)", fontsize=14)
ax.legend(loc="lower right", fontsize=8, ncol=2)
ax.grid(True, alpha=0.3)
ax.set_xticks(k_list)
plt.tight_layout()
plt.savefig("p1_k_vs_auc_all_models.png", dpi=150)
plt.show()

# --------------------- 11. 结论与下一步 ---------------------
print("\n" + "="*60)
best_model = max(ensemble_results.keys(), key=lambda k: ensemble_results[k][0])
best_auc, best_std = ensemble_results[best_model]
best_base, _ = baseline_results[best_model]
print(f"Best Model         : {best_model}")
print(f"Best Ensemble AUC  : {best_auc:.4f} ± {best_std:.4f}")
print(f"Best Baseline AUC  : {best_base:.4f}")
print(f"Delta              : {best_auc - best_base:+.4f}")
print("="*60)

# 饱和判断
if len(prog_results[best_model]) >= 3:
    last = prog_results[best_model][-1][1]
    prev = prog_results[best_model][-2][1]
    print(f"\nK={prog_results[best_model][-2][0]}→K={prog_results[best_model][-1][0]} 增益: {last-prev:+.4f}")

print("\n" + "="*60)
if best_auc >= 0.76:
    print("✅ 结果优秀，建议直接进入 P1.5：SCI + Mutation Features")
elif best_auc >= 0.74:
    print("⚠️  有提升但边际递减，建议：1) 特征选择去冗余；2) 进入P1.5")
else:
    print("❌ Ensemble提升有限，建议转做 Layer32 Single-Layer 多统计量")
print("="*60)