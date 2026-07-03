import pandas as pd
import matplotlib.pyplot as plt
import os

plt.rcParams['font.size'] = 9
# ===================== 全局路径统一配置（适配你本地文件夹结构）=====================
BASE_DATA = r"D:\文件\工作室\website\data\processed"
PABP_CSV = os.path.join(BASE_DATA, "p0_output", "PABP", "PABP_layer_pair_metrics.csv")
COMPARE_PROTEINS = ["CBS", "GAL4", "PTEN", "TEM1"]
OUT_DIR = os.path.join(BASE_DATA, "p0_output", "P02_PABP_analysis")
os.makedirs(OUT_DIR, exist_ok=True)
category_list = ["EE", "EM", "EL", "MM", "ML", "LL"]

# 读取PABP全量层对数据
df_pabp = pd.read_csv(PABP_CSV)
df_pabp["abs_spear"] = df_pabp["Spearman_r"].abs()

# 1. PABP自身六类指标统计表
pabp_stat = []
for cat in category_list:
    sub = df_pabp[df_pabp["Pair_Category"] == cat]
    pabp_stat.append({
        "Category": cat,
        "Count": len(sub),
        "AUROC_Mean": round(sub["AUROC_abs"].mean(), 4),
        "Composite_Mean": round(sub["CompositeScore"].mean(), 4)
    })
df_pabp_stat = pd.DataFrame(pabp_stat)
df_pabp_stat.to_csv(os.path.join(OUT_DIR, "PABP_6category_stats.csv"), index=False)

# 2. PABP六类综合得分箱线图
plt.figure(figsize=(10, 6), dpi=150)
pabp_com_data = [df_pabp[df_pabp["Pair_Category"] == c]["CompositeScore"].values for c in category_list]
plt.boxplot(pabp_com_data, labels=category_list)
plt.ylabel("Composite Score")
plt.xlabel("Layer Pair Category")
plt.title("PABP All Layer Pairs Composite Score Distribution")
plt.grid(axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "boxplot_composite_PABP.png"))
plt.close()

# 3. PABP与其余四类蛋白均值横向对比（修复你截图里残缺代码）
compare_result = []
# 先存入PABP数据
base_row = {"Protein": "PABP"}
for cat in category_list:
    base_row[cat] = df_pabp[df_pabp["Pair_Category"] == cat]["CompositeScore"].mean()
compare_result.append(base_row)

# 循环读取其他蛋白
for prot in COMPARE_PROTEINS:
    prot_csv_path = os.path.join(BASE_DATA, "p0_output", prot, f"{prot}_layer_pair_metrics.csv")
    df_prot = pd.read_csv(prot_csv_path)
    temp_row = {"Protein": prot}
    for cat in category_list:
        temp_row[cat] = df_prot[df_prot["Pair_Category"] == cat]["CompositeScore"].mean()
    compare_result.append(temp_row)

df_compare = pd.DataFrame(compare_result)
df_compare.to_csv(os.path.join(OUT_DIR, "PABP_vs_other_proteins.csv"), index=False)

# 4. 提取PABP Top50内EE配对清单（修正路径拼接方式，杜绝斜杠错误）
top50_pabp_path = os.path.join(BASE_DATA, "p0_output", "PABP", "top50_composite_PABP.csv")
top50_pabp = pd.read_csv(top50_pabp_path)
top_ee_pairs = top50_pabp[top50_pabp["Pair_Category"] == "EE"]
top_ee_pairs.to_csv(os.path.join(OUT_DIR, "PABP_top50_EE_pairs.csv"), index=False)

print(f"PABP Top50中EE配对数量：{len(top_ee_pairs)}")
print("==== PABP与其余蛋白综合得分对比 ====")
print(df_compare.round(4))