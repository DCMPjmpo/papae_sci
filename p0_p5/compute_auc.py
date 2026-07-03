import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_curve, auc
import matplotlib.pyplot as plt
import os

# ========== 路径 ==========
OUTPUT_DIR = r"/mnt/sda/gws_1020251255/data/processed"
PLOT_DIR = r"/mnt/sda/gws_1020251255/plots"
os.makedirs(PLOT_DIR, exist_ok=True)

# ========== 读取数据 ==========
sci = np.load(os.path.join(OUTPUT_DIR, "pten_tem1_sci_scores.npy"))
meta = pd.read_csv(os.path.join(OUTPUT_DIR, "pten_tem1_metadata.csv"))
labels = meta["DMS_score_bin"].values  # 0 = pathogenic, 1 = neutral
proteins = meta["protein"].values

# 预测分数：SCI 越高越可能是 pathogenic
# 所以 labels 需要翻转：pathogenic = 1
y_true = 1 - labels

# ========== 总体 AUC ==========
auc_roc = roc_auc_score(y_true, sci)
precision, recall, _ = precision_recall_curve(y_true, sci)
auc_pr = auc(recall, precision)  # 修复：x=recall, y=precision

print(f"{'='*50}")
print("总体预测性能:")
print(f"{'='*50}")
print(f"ROC-AUC: {auc_roc:.4f}")
print(f"PR-AUC:  {auc_pr:.4f}")

# ========== 按蛋白 AUC ==========
print(f"\n{'='*50}")
print("按蛋白预测性能:")
print(f"{'='*50}")

for protein in ["PTEN", "TEM1"]:
    mask = proteins == protein
    sub_y_true = y_true[mask]
    sub_sci = sci[mask]
    
    auc_roc_p = roc_auc_score(sub_y_true, sub_sci)
    precision_p, recall_p, _ = precision_recall_curve(sub_y_true, sub_sci)
    auc_pr_p = auc(recall_p, precision_p)
    
    print(f"\n{protein}:")
    print(f"  ROC-AUC: {auc_roc_p:.4f}")
    print(f"  PR-AUC:  {auc_pr_p:.4f}")

# ========== Negative Control：随机标签 ==========
print(f"\n{'='*50}")
print("Negative Control: 随机标签")
print(f"{'='*50}")

np.random.seed(42)
n_runs = 100
random_aucs = []

for _ in range(n_runs):
    shuffled = np.random.permutation(y_true)
    random_aucs.append(roc_auc_score(shuffled, sci))

random_mean = np.mean(random_aucs)
random_std = np.std(random_aucs)

print(f"随机标签 AUC: {random_mean:.4f} ± {random_std:.4f}")
print(f"真实 AUC:     {auc_roc:.4f}")
print(f"差异:         {auc_roc - random_mean:.4f}")
print(f"{'✅ 信号真实!' if auc_roc > random_mean + 2*random_std else '⚠️ 需要检查'}")

# ========== ROC 曲线 ==========
fig, ax = plt.subplots(figsize=(8, 6))

# 总体
fpr, tpr, _ = roc_curve(y_true, sci)
ax.plot(fpr, tpr, label=f"Overall (AUC={auc_roc:.3f})", linewidth=2, color='black')

# 按蛋白
for protein, color in zip(["PTEN", "TEM1"], ["#e74c3c", "#3498db"]):
    mask = proteins == protein
    sub_y_true = y_true[mask]
    sub_sci = sci[mask]
    fpr_p, tpr_p, _ = roc_curve(sub_y_true, sub_sci)
    auc_p = roc_auc_score(sub_y_true, sub_sci)
    ax.plot(fpr_p, tpr_p, label=f"{protein} (AUC={auc_p:.3f})", linewidth=2, color=color)

# 随机标签参考线
ax.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Random (AUC=0.5)')

ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
ax.set_title("ROC Curve: SCI for Pathogenic Prediction")
ax.legend(loc='lower right')
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, "roc_curve.png"), dpi=300, bbox_inches="tight")
plt.show()

print(f"\n已保存: {os.path.join(PLOT_DIR, 'roc_curve.png')}")