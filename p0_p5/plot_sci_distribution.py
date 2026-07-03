import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

# ========== 路径 ==========
OUTPUT_DIR = r"/mnt/sda/gws_1020251255/data/processed"
PLOT_DIR = r"/mnt/sda/gws_1020251255/plots"
os.makedirs(PLOT_DIR, exist_ok=True)

# ========== 读取数据 ==========
sci = np.load(os.path.join(OUTPUT_DIR, "pten_tem1_sci_scores.npy"))
meta = pd.read_csv(os.path.join(OUTPUT_DIR, "pten_tem1_metadata.csv"))
labels = meta["DMS_score_bin"].values
proteins = meta["protein"].values

df_plot = pd.DataFrame({
    "SCI": sci,
    "Label": ["Pathogenic" if l == 0 else "Neutral" for l in labels],
    "Protein": proteins
})

# ========== 辅助函数：修正 Cliff's Delta ==========
def cliffs_delta(x, y):
    """
    标准 Cliff's Delta 实现
    |d| < 0.147: 微小
    |d| < 0.33: 小
    |d| < 0.474: 中等
    |d| >= 0.474: 大
    """
    nx = len(x)
    ny = len(y)
    u, _ = stats.mannwhitneyu(x, y, alternative="two-sided")
    return (2 * u) / (nx * ny) - 1

def effect_size_label(d):
    abs_d = abs(d)
    if abs_d < 0.147:
        return "微小"
    elif abs_d < 0.33:
        return "小"
    elif abs_d < 0.474:
        return "中等"
    else:
        return "大"

# ========== 1. 总体分布 ==========
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

sns.boxplot(data=df_plot, x="Label", y="SCI", ax=axes[0], 
            palette=["#e74c3c", "#2ecc71"], showfliers=False)
axes[0].set_title("SCI Distribution: Boxplot")
axes[0].set_ylabel("SCI Score")

sns.violinplot(data=df_plot, x="Label", y="SCI", ax=axes[1], 
               palette=["#e74c3c", "#2ecc71"], inner="box")
axes[1].set_title("SCI Distribution: Violin Plot")
axes[1].set_ylabel("SCI Score")

plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, "sci_distribution.png"), dpi=300, bbox_inches="tight")
plt.show()

# ========== 2. 总体统计 + 效应量 ==========
path_all = df_plot[df_plot["Label"] == "Pathogenic"]["SCI"].values
neut_all = df_plot[df_plot["Label"] == "Neutral"]["SCI"].values

u, p = stats.mannwhitneyu(path_all, neut_all, alternative='greater')
cd = cliffs_delta(path_all, neut_all)

print(f"\n{'='*50}")
print("总体统计:")
print(f"{'='*50}")
print(f"Pathogenic (n={len(path_all)}): {path_all.mean():.6f} ± {path_all.std():.6f}")
print(f"Neutral    (n={len(neut_all)}): {neut_all.mean():.6f} ± {neut_all.std():.6f}")
print(f"Δ:          {path_all.mean() - neut_all.mean():.6f}")
print(f"Mann-Whitney U: p = {p:.6f}")
print(f"Cliff's Delta:  {cd:.4f} ({effect_size_label(cd)})")

# ========== 3. 按蛋白拆分 ==========
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for idx, protein in enumerate(["PTEN", "TEM1"]):
    sub = df_plot[df_plot["Protein"] == protein]
    path = sub[sub["Label"] == "Pathogenic"]["SCI"].values
    neut = sub[sub["Label"] == "Neutral"]["SCI"].values
    u, p = stats.mannwhitneyu(path, neut, alternative='greater')
    cd = cliffs_delta(path, neut)
    
    sns.violinplot(data=sub, x="Label", y="SCI", ax=axes[idx], 
                   palette=["#e74c3c", "#2ecc71"], inner="box")
    axes[idx].set_title(f"{protein} (n={len(sub)}): SCI Distribution")
    axes[idx].set_ylabel("SCI Score")
    axes[idx].text(0.5, 0.95, 
                   f"Δ = {path.mean() - neut.mean():.4f}\np = {p:.6f}\nCliff's d = {cd:.3f} ({effect_size_label(cd)})", 
                   transform=axes[idx].transAxes, ha='center', va='top',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, "sci_per_protein.png"), dpi=300, bbox_inches="tight")
plt.show()

# ========== 4. 按蛋白统计摘要 ==========
print(f"\n{'='*50}")
print("按蛋白统计摘要:")
print(f"{'='*50}")
for protein in ["PTEN", "TEM1"]:
    sub = df_plot[df_plot["Protein"] == protein]
    path = sub[sub["Label"] == "Pathogenic"]["SCI"].values
    neut = sub[sub["Label"] == "Neutral"]["SCI"].values
    u, p = stats.mannwhitneyu(path, neut, alternative='greater')
    cd = cliffs_delta(path, neut)
    
    print(f"\n{protein}:")
    print(f"  Pathogenic (n={len(path)}): {path.mean():.4f} ± {path.std():.4f}")
    print(f"  Neutral    (n={len(neut)}): {neut.mean():.4f} ± {neut.std():.4f}")
    print(f"  Δ:          {path.mean() - neut.mean():.4f}")
    print(f"  p-value:    {p:.6f}")
    print(f"  Cliff's d:  {cd:.4f} ({effect_size_label(cd)})")
    print(f"  {'*** 显著!' if p < 0.001 else '** 显著!' if p < 0.01 else '* 显著!' if p < 0.05 else '不显著'}")