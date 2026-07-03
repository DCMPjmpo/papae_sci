import pandas as pd
import matplotlib.pyplot as plt
import os

# 配置
# [旧·服务器路径]
# OUTPUT_DIR = "/mnt/sda/gws_1020251255/data/processed/p0_output"
# [新·本地路径]
OUTPUT_DIR = "D:/文件/工作室/website/data/p0_output"
TOP_N = 10

# 1. 读取跨蛋白层对数据
csv_path = os.path.join(OUTPUT_DIR, "cross_protein_universal_pairs.csv")
df = pd.read_csv(csv_path)

# 2. 筛选Top10高频层对
df_top = df.head(TOP_N).copy()
# 格式化层对名称（1-based）
df_top["Layer_Pair"] = df_top.apply(lambda x: f"{x['layer_i_1']}-{x['layer_j_1']}", axis=1)

# 3. 绘制横向柱状图（论文标准格式）
plt.figure(figsize=(10, 6), dpi=150)
plt.barh(df_top["Layer_Pair"], df_top["Cross_Protein_Count"], color="#4472c4")

# 美化
plt.xlabel("Number of Proteins the Pair Appears in (Top50 Composite)", fontsize=12)
plt.ylabel("Layer Pair (1-based ESM Index)", fontsize=12)
plt.title(f"Top{TOP_N} Recurrent Layer Pairs Across 5 Proteins", fontsize=14)
plt.gca().invert_yaxis()  # 最高频在最上面
plt.grid(axis="x", linestyle="--", alpha=0.7)
plt.tight_layout()

# 保存
plot_path = os.path.join(OUTPUT_DIR, "top_recurrent_cross_protein_pairs.png")
plt.savefig(plot_path, bbox_inches="tight")
plt.close()

# 输出结果
print(f"===== Top{TOP_N} Recurrent Cross-Protein Layer Pairs =====")
print(df_top[["Layer_Pair", "Cross_Protein_Count"]])
print(f"\n图表已保存至: {plot_path}")