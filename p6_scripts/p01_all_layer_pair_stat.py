import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import os

# 全局配置
plt.rcParams['font.size'] = 9
# 修复：删除多余注释行、引号修正，Windows原始字符串标准写法
INPUT_CSV = r"D:\文件\工作室\website\data\processed\p0_output\All\All_layer_pair_metrics.csv"
OUT_DIR = r"D:\文件\工作室\website\data\processed\p0_output\P01_full_statistics"
os.makedirs(OUT_DIR, exist_ok=True)

# 读取全部528条层对数据
df = pd.read_csv(INPUT_CSV)
df["abs_spear"] = df["Spearman_r"].abs()
category_list = ["EE", "EM", "EL", "MM", "ML", "LL"]

# 1. 六类层对汇总统计表
stat_res = []
for cat in category_list:
    sub_df = df[df["Pair_Category"] == cat]
    row = {
        "Category": cat,
        "Sample_Count": len(sub_df),
        # 修复变量名 sub → sub_df
        "AUROC_abs_Mean": round(sub_df["AUROC_abs"].mean(), 4),
        "AUROC_abs_Median": round(sub_df["AUROC_abs"].median(), 4),
        "AUROC_abs_Std": round(sub_df["AUROC_abs"].std(), 4),
        "Composite_Mean": round(sub_df["CompositeScore"].mean(), 4),
        "Composite_Median": round(sub_df["CompositeScore"].median(), 4),
        "Abs_Spearman_Mean": round(sub_df["abs_spear"].mean(), 4)
    }
    stat_res.append(row)
stat_df = pd.DataFrame(stat_res)
stat_csv_path = os.path.join(OUT_DIR, "all_6category_stats.csv")
stat_df.to_csv(stat_csv_path, index=False)
print("==== 全部528层对六类指标汇总 ====")
print(stat_df)

# 2. CompositeScore 箱线图
plt.figure(figsize=(10, 6), dpi=150)
composite_data = [df[df["Pair_Category"] == c]["CompositeScore"].values for c in category_list]
plt.boxplot(composite_data, labels=category_list)
plt.ylabel("Composite Score")
plt.xlabel("Layer Pair Category")
plt.title("Composite Score Distribution (All 528 Layer Pairs)")
plt.grid(axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "boxplot_composite_all.png"))
plt.close()

# 3. Corrected AUROC 箱线图
plt.figure(figsize=(10, 6), dpi=150)
auroc_data = [df[df["Pair_Category"] == c]["AUROC_abs"].values for c in category_list]
plt.boxplot(auroc_data, labels=category_list)
plt.ylabel("Corrected AUROC")
plt.xlabel("Layer Pair Category")
plt.title("AUROC Distribution (All 528 Layer Pairs)")
plt.grid(axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "boxplot_auroc_all.png"))
plt.close()

# 4. Mann-Whitney U 检验：EL组分别对比其余五类
el_group = df[df["Pair_Category"] == "EL"]["CompositeScore"].values
test_log = []
test_log.append("Mann-Whitney U Test: EL vs other layer pair categories\n")
for other_cat in ["EE", "EM", "MM", "ML", "LL"]:
    other_group = df[df["Pair_Category"] == other_cat]["CompositeScore"].values
    stat_val, p_val = mannwhitneyu(el_group, other_group, alternative="greater")
    test_log.append(f"EL vs {other_cat}: U = {stat_val:.2f}, p-value = {p_val:.4e}")

# 保存检验文本
txt_path = os.path.join(OUT_DIR, "stat_test_result.txt")
with open(txt_path, "w", encoding="utf-8") as f:
    f.write("\n".join(test_log))
print("\n".join(test_log))
print(f"\n所有图表、统计文件已输出至 {OUT_DIR}")