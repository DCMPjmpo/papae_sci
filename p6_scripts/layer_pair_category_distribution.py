import pandas as pd
import matplotlib.pyplot as plt
import os

# 配置
# [旧·服务器路径]
# OUTPUT_DIR = "//mnt/sda/gws_1020251255/data/processed/p0_output"
# [新·本地路径]
OUTPUT_DIR = "D:/文件/工作室/website/data/p0_output"
PROTEIN_LIST = ["All", "CBS", "GAL4", "PABP", "PTEN", "TEM1"]
CATEGORIES = ["EE", "EM", "EL", "MM", "ML", "LL"]
CATEGORY_NAMES = {
    "EE": "Early-Early",
    "EM": "Early-Middle",
    "EL": "Early-Late",
    "MM": "Middle-Middle",
    "ML": "Middle-Late",
    "LL": "Late-Late"
}

# 1. 读取所有蛋白的层对占比数据
all_data = []
for prot in PROTEIN_LIST:
    csv_path = os.path.join(OUTPUT_DIR, prot, f"{prot}_pair_category_ratio_Composite.csv")
    if not os.path.exists(csv_path):
        print(f"跳过{prot}，文件不存在")
        continue
    df = pd.read_csv(csv_path)
    prot_data = {"Protein": prot}
    for _, row in df.iterrows():
        prot_data[row["Pair_Type"]] = row["Ratio"]
    all_data.append(prot_data)

# 2. 转为DataFrame
df_plot = pd.DataFrame(all_data).set_index("Protein")
df_plot = df_plot[CATEGORIES]  # 固定顺序

# 3. 绘制堆叠柱状图（论文标准格式）
plt.figure(figsize=(12, 7), dpi=150)
ax = df_plot.plot(kind="bar", stacked=True, colormap="RdBu_r", ax=plt.gca())

# 美化
plt.xlabel("Protein Subset", fontsize=12)
plt.ylabel("Proportion in Composite Top50 Layer Pairs", fontsize=12)
plt.title("Distribution of Layer-Pair Categories in High-Contribution SCI Pairs", fontsize=14)
plt.xticks(rotation=0, fontsize=10)
plt.legend([CATEGORY_NAMES[c] for c in CATEGORIES], title="Layer-Pair Category", bbox_to_anchor=(1.01, 1), loc="upper left")
plt.tight_layout()

# 保存
plot_path = os.path.join(OUTPUT_DIR, "layer_pair_category_distribution.png")
plt.savefig(plot_path, bbox_inches="tight")
plt.close()

# 4. 输出全局统计结果
print("===== Global Layer-Pair Category Distribution (All Proteins) =====")
df_all = df_plot.loc["All"]
for cat in CATEGORIES:
    print(f"{CATEGORY_NAMES[cat]}: {df_all[cat]:.2%}")

print(f"\n图表已保存至: {plot_path}")