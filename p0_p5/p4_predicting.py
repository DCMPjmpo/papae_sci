import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.decomposition import PCA
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_auc_score

try:
    import xgboost as xgb
    XGB_AVAILABLE = True
except ImportError:
    XGB_AVAILABLE = False
    print("[Warning] xgboost not installed. Run: pip install xgboost")

# ===================== 数据加载 =====================
matrices = np.load('data/processed/pten_tem1_sci_site_matrices.npy')
metadata = pd.read_csv('data/processed/pten_tem1_site_metadata.csv')

# 真实 ESM-2 delta embedding
delta_emb = np.load('data/processed/pten_tem1_delta_embeddings.npy')

y = metadata["DMS_score_bin"].values
protein = metadata["protein"].astype(str).values

print(f"[Info] 总样本数: {len(y)}")
print(f"[Info] Delta embedding shape: {delta_emb.shape}")
print(f"[Info] Pathogenic: {(y==0).sum()}, Neutral: {(y==1).sum()}")
print(f"[Info] Protein分布: {dict(zip(*np.unique(protein, return_counts=True)))}")

# ===================== 特征工程 =====================

# --- 1. Position Baseline ---
positions = metadata["mutation_position"].values
X_position = positions.reshape(-1, 1) / positions.max()
print(f"[Info] Position 特征维度: {X_position.shape[1]}")

# --- 2. ESM Delta Embedding (降维) ---
# 全局平均 + Late层平均 + 统计量
X_esm_global = delta_emb.mean(axis=1)
X_esm_late = delta_emb[:, 22:, :].mean(axis=1)
X_esm_raw = np.column_stack([
    X_esm_global,
    X_esm_late,
    delta_emb.std(axis=1),
    delta_emb.max(axis=1),
])
print(f"[Info] ESM raw 维度: {X_esm_raw.shape[1]}")

# --- 3. SCI 特征 (仅模块均值，无泄漏) ---
# 模块定义来自理论，非数据驱动
EARLY = list(range(0, 11))
MIDDLE = list(range(11, 22))
LATE = list(range(22, 33))

def compute_module_sci(mats):
    ml_pairs = [mats[:, i, j] for i in MIDDLE for j in LATE if i < j]
    el_pairs = [mats[:, i, j] for i in EARLY for j in LATE if i < j]
    ll_pairs = [mats[:, i, j] for i in LATE for j in LATE if i < j]

    return np.column_stack([
        np.mean(ml_pairs, axis=0),
        np.mean(el_pairs, axis=0),
        np.mean(ll_pairs, axis=0),
    ])

X_sci = compute_module_sci(matrices)
print(f"[Info] SCI 特征维度: {X_sci.shape[1]} (Middle-Late, Early-Late, Late-Late)")

# --- 4. Fusion ---
X_fusion_raw = np.column_stack([X_esm_raw, X_sci])
print(f"[Info] Fusion raw 维度: {X_fusion_raw.shape[1]}")

# ===================== 模型定义 (PCA降维) =====================
def make_lr_pca(n_comp=100):
    return Pipeline([
        ("scaler", StandardScaler()),
        ("pca", PCA(n_components=n_comp, random_state=42)),
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

# ESM/Fusion 用 PCA+LR，避免高维灾难
# SCI/Position 维度低，直接用原始模型
MODELS_LOW = {"LR": LogisticRegression(max_iter=1000, class_weight="balanced"), 
              "RF": make_rf()}
if XGB_AVAILABLE:
    MODELS_LOW["XGB"] = make_xgb()

MODELS_HIGH = {"LR_PCA100": make_lr_pca(100), 
               "RF": make_rf()}
if XGB_AVAILABLE:
    MODELS_HIGH["XGB"] = make_xgb()

# ===================== 评估函数 =====================
CV_5FOLD = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

def cv_auc(X, y, model, cv, name="Model"):
    # 自动调整 PCA n_components，避免低维数据报错
    model_copy = model
    if hasattr(model, 'named_steps') and 'pca' in model.named_steps:
        from sklearn.decomposition import PCA
        n_comp = min(100, X.shape[1])
        model_copy = Pipeline([
            ("scaler", StandardScaler()),
            ("pca", PCA(n_components=n_comp, random_state=42)),
            ("clf", model.named_steps['clf'])
        ])
    scores = cross_val_score(model_copy, X, y, cv=cv, scoring="roc_auc", n_jobs=-1)
    print(f"  {name:35s}: AUC = {scores.mean():.4f} ± {scores.std():.4f}")
    return scores.mean(), scores.std()

def lopo_auc(X, y, protein, model, train_prot, test_prot, name="Model"):
    train_mask = protein == train_prot
    test_mask = protein == test_prot
    if train_mask.sum() == 0 or test_mask.sum() == 0:
        return None
    # 自动调整 PCA n_components
    model_copy = model
    if hasattr(model, 'named_steps') and 'pca' in model.named_steps:
        from sklearn.decomposition import PCA
        n_comp = min(100, X.shape[1])
        model_copy = Pipeline([
            ("scaler", StandardScaler()),
            ("pca", PCA(n_components=n_comp, random_state=42)),
            ("clf", model.named_steps['clf'])
        ])
    model_copy.fit(X[train_mask], y[train_mask])
    y_pred = model_copy.predict_proba(X[test_mask])[:, 1]
    auc = roc_auc_score(y[test_mask], y_pred)
    print(f"  {name:35s}: AUC = {auc:.4f}  (train={train_prot}, test={test_prot})")
    return auc

# ===================== P4A: Position vs SCI =====================
print("\n" + "="*70)
print("P4A: Position Baseline vs SCI (Low-dim, no PCA needed)")
print("="*70)

p4a_results = {}
for model_name, model in MODELS_LOW.items():
    print(f"\n--- {model_name} ---")
    pos_mu, pos_std = cv_auc(X_position, y, model, CV_5FOLD, f"{model_name} (Position)")
    sci_mu, sci_std = cv_auc(X_sci, y, model, CV_5FOLD, f"{model_name} (SCI)")
    p4a_results[model_name] = {"Position": (pos_mu, pos_std), "SCI": (sci_mu, sci_std)}

# ===================== P4B: ESM vs SCI vs Fusion =====================
print("\n" + "="*70)
print("P4B: ESM vs SCI vs Fusion (ESM/Fusion use PCA+LR)")
print("="*70)

p4b_results = {}
for model_name, model in MODELS_HIGH.items():
    print(f"\n--- {model_name} ---")
    esm_mu, esm_std = cv_auc(X_esm_raw, y, model, CV_5FOLD, f"{model_name} (ESM)")
    sci_mu, sci_std = cv_auc(X_sci, y, model, CV_5FOLD, f"{model_name} (SCI)")
    fusion_mu, fusion_std = cv_auc(X_fusion_raw, y, model, CV_5FOLD, f"{model_name} (Fusion)")
    p4b_results[model_name] = {
        "ESM": (esm_mu, esm_std),
        "SCI": (sci_mu, sci_std),
        "Fusion": (fusion_mu, fusion_std)
    }

# ===================== P4C: Leave-One-Protein-Out =====================
print("\n" + "="*70)
print("P4C: Leave-One-Protein-Out (Pilot, n=2 proteins)")
print("="*70)

unique_prots = np.unique(protein)
p4c_results = {p: {} for p in unique_prots}

# 用 LR_PCA100 做LOPO (ESM维度高需要PCA)
lopo_model = make_lr_pca(100)

for test_prot in unique_prots:
    train_prot = [p for p in unique_prots if p != test_prot][0]
    print(f"\n--- Test on {test_prot}, Train on {train_prot} ---")
    for feat_name, X_feat in [("ESM", X_esm_raw), ("SCI", X_sci), ("Fusion", X_fusion_raw)]:
        auc = lopo_auc(X_feat, y, protein, lopo_model, train_prot, test_prot, feat_name)
        if auc is not None:
            p4c_results[test_prot][feat_name] = auc

# ===================== P4D: SCI Module Ablation =====================
print("\n" + "="*70)
print("P4D: SCI Module Ablation (3 features only)")
print("="*70)

X_ml = X_sci[:, 0:1]   # Middle-Late
X_el = X_sci[:, 1:2]   # Early-Late
X_ll = X_sci[:, 2:3]   # Late-Late
X_ml_el = X_sci[:, 0:2]

abl_features = {
    "ML_only": X_ml,
    "EL_only": X_el,
    "LL_only": X_ll,
    "ML+EL": X_ml_el,
    "All_3": X_sci,
}

p4d_results = {}
for feat_name, X_feat in abl_features.items():
    print(f"\n--- {feat_name} (dim={X_feat.shape[1]}) ---")
    for model_name, model in MODELS_LOW.items():
        mu, std = cv_auc(X_feat, y, model, CV_5FOLD, f"{model_name}")
        p4d_results[f"{feat_name}_{model_name}"] = (mu, std)

# ===================== 结果汇总 =====================
print("\n" + "="*70)
print("【P4 完整结果汇总】")
print("="*70)

# P4A
print("\n--- P4A: Position vs SCI ---")
print(f"{'Model':<10} {'Position':<12} {'SCI':<12} {'SCI>Pos?':<10}")
print("-"*50)
for m in MODELS_LOW.keys():
    pos = p4a_results[m]["Position"][0]
    sci = p4a_results[m]["SCI"][0]
    print(f"{m:<10} {pos:<12.4f} {sci:<12.4f} {'✅' if sci > pos else '❌'}")

# P4B
print("\n--- P4B: ESM vs SCI vs Fusion ---")
print(f"{'Model':<12} {'ESM':<12} {'SCI':<12} {'Fusion':<12} {'Fus>ESM?':<10} {'Fus>SCI?':<10}")
print("-"*70)
for m in MODELS_HIGH.keys():
    esm = p4b_results[m]["ESM"][0]
    sci = p4b_results[m]["SCI"][0]
    fus = p4b_results[m]["Fusion"][0]
    print(f"{m:<12} {esm:<12.4f} {sci:<12.4f} {fus:<12.4f} {'✅' if fus > esm else '❌':<10} {'✅' if fus > sci else '❌':<10}")

# P4C
print("\n--- P4C: Leave-One-Protein-Out ---")
for test_prot in unique_prots:
    print(f"\nTest on {test_prot}:")
    for k, v in sorted(p4c_results[test_prot].items(), key=lambda x: x[1], reverse=True):
        print(f"  {k:<20s}: {v:.4f}")

# P4D
print("\n--- P4D: Module Ablation ---")
for feat_name in abl_features.keys():
    best_mu = 0
    best_model = ""
    for m in MODELS_LOW.keys():
        key = f"{feat_name}_{m}"
        if key in p4d_results:
            mu, _ = p4d_results[key]
            if mu > best_mu:
                best_mu = mu
                best_model = m
    print(f"  {feat_name:<15s}: {best_mu:.4f}  (best model: {best_model})")

# ===================== 画图 =====================
fig = plt.figure(figsize=(18, 12))
gs = fig.add_gridspec(3, 3, hspace=0.4, wspace=0.35)

# A: P4A Position vs SCI
ax = fig.add_subplot(gs[0, 0])
models = list(MODELS_LOW.keys())
x = np.arange(len(models))
width = 0.35
pos_vals = [p4a_results[m]["Position"][0] for m in models]
sci_vals = [p4a_results[m]["SCI"][0] for m in models]
ax.bar(x - width/2, pos_vals, width, label='Position', color='gray', edgecolor='black')
ax.bar(x + width/2, sci_vals, width, label='SCI', color='steelblue', edgecolor='black')
ax.set_ylabel('AUROC')
ax.set_title('A. Position vs SCI', fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(models)
ax.legend()
ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.3)
ax.grid(True, alpha=0.3, axis='y')

# B: P4B ESM vs SCI vs Fusion
ax = fig.add_subplot(gs[0, 1])
width = 0.25
for i, feat in enumerate(["ESM", "SCI", "Fusion"]):
    vals = [p4b_results[m][feat][0] for m in MODELS_HIGH.keys()]
    errs = [p4b_results[m][feat][1] for m in MODELS_HIGH.keys()]
    ax.bar(x + i*width, vals, width, yerr=errs, label=feat, capsize=3)
ax.set_ylabel('AUROC')
ax.set_title('B. ESM vs SCI vs Fusion', fontweight='bold')
ax.set_xticks(x + width)
ax.set_xticklabels(list(MODELS_HIGH.keys()))
ax.legend()
ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.3)
ax.grid(True, alpha=0.3, axis='y')

# C: P4C LOPO
ax = fig.add_subplot(gs[0, 2])
lopo_data = []
for test_prot in unique_prots:
    for feat, auc in p4c_results[test_prot].items():
        lopo_data.append({'Test': test_prot, 'Feature': feat, 'AUC': auc})
lopo_df = pd.DataFrame(lopo_data)
lopo_pivot = lopo_df.pivot(index='Test', columns='Feature', values='AUC')
lopo_pivot.plot(kind='bar', ax=ax, color=['coral', 'steelblue', 'green'])
ax.set_ylabel('AUROC')
ax.set_title('C. LOPO (Pilot, n=2)', fontweight='bold')
ax.legend(title='Feature')
ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.3)
plt.setp(ax.xaxis.get_majorticklabels(), rotation=0)

# D: P4D Ablation
ax = fig.add_subplot(gs[1, :])
abl_names = list(abl_features.keys())
abl_vals = [max(p4d_results.get(f"{n}_{m}", (0,0))[0] for m in MODELS_LOW.keys()) for n in abl_names]
colors = plt.cm.RdYlGn(np.linspace(0.2, 0.8, len(abl_names)))
ax.bar(abl_names, abl_vals, color=colors, edgecolor='black')
ax.set_ylabel('AUROC (Best Model)')
ax.set_title('D. SCI Module Ablation', fontweight='bold')
ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.3)
plt.setp(ax.xaxis.get_majorticklabels(), rotation=15, ha='right')
ax.grid(True, alpha=0.3, axis='y')

# E: 增益矩阵
ax = fig.add_subplot(gs[2, 0])
gain_data = []
for m in MODELS_HIGH.keys():
    esm = p4b_results[m]["ESM"][0]
    sci = p4b_results[m]["SCI"][0]
    fus = p4b_results[m]["Fusion"][0]
    gain_data.append({
        'Model': m,
        'SCI vs ESM': sci - esm,
        'Fusion vs ESM': fus - esm,
        'Fusion vs SCI': fus - sci
    })
gain_df = pd.DataFrame(gain_data).set_index('Model')
sns.heatmap(gain_df, annot=True, fmt='.4f', cmap='RdYlGn', center=0,
            vmin=-0.05, vmax=0.1, square=True, ax=ax,
            cbar_kws={'shrink': 0.8})
ax.set_title('E. Performance Gain', fontweight='bold')

# F: 关键数字总结
ax = fig.add_subplot(gs[2, 1:])
ax.axis('off')

best_fus = max(p4b_results[m]["Fusion"][0] for m in MODELS_HIGH.keys())
best_esm = max(p4b_results[m]["ESM"][0] for m in MODELS_HIGH.keys())
best_sci = max(p4b_results[m]["SCI"][0] for m in MODELS_HIGH.keys())
best_pos = max(p4a_results[m]["Position"][0] for m in MODELS_LOW.keys())

summary_text = f"""
P4 关键结果 (Best across models)
{'='*45}
Position Baseline:   {best_pos:.4f}
ESM (PCA100):        {best_esm:.4f}
SCI (3 modules):     {best_sci:.4f}
ESM + SCI Fusion:    {best_fus:.4f}

Fusion Gain vs ESM:   {best_fus - best_esm:+.4f}
Fusion Gain vs SCI:   {best_fus - best_sci:+.4f}
SCI Gain vs Position: {best_sci - best_pos:+.4f}

{'='*45}
"""

if best_fus > best_esm + 0.03 and best_sci > best_pos + 0.05:
    summary_text += "🔥🔥🔥 核心结论成立!\n"
    summary_text += "   SCI > Position (排除位置偏置)\n"
    summary_text += "   Fusion > ESM (3维模块提升高维ESM)\n"
elif best_fus > best_esm:
    summary_text += "✅ Fusion 优于单一特征\n"
else:
    summary_text += "⚠️ 增益有限\n"

ax.text(0.1, 0.5, summary_text, transform=ax.transAxes, fontsize=11,
        verticalalignment='center', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.suptitle('Figure 6: SCI Predictive Value (No-Leakage Version)',
             fontsize=14, fontweight='bold', y=0.98)

plt.savefig('data/layer/figure6_p4_no_leakage.png', dpi=300, bbox_inches='tight')
plt.show()

print("\n✅ Saved: figure6_p4_no_leakage.png")