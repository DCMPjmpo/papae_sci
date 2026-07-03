import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import roc_auc_score, average_precision_score

# ===================== 路径配置（适配本地文件） =====================
# [旧·服务器路径]
# SCI_1D_PATH = "/mnt/sda/gws_1020251255/data/processed/all_proteins_sci_site_scores_mean.npy"
# META_PATH = "/mnt/sda/gws_1020251255/data/processed/all_proteins_site_metadata.csv"
# SAVE_FIG = "/mnt/sda/gws_1020251255/data/processed/P1.5_sci_signal_validation.png"
# SAVE_TABLE = "/mnt/sda/gws_1020251255/data/processed/P1.5_signal_summary.csv"
# [新·本地路径]
SCI_1D_PATH = "D:/文件/工作室/website/data/processed/all_proteins_sci_site_scores_mean.npy"
META_PATH = "D:/文件/工作室/website/data/processed/all_proteins_site_metadata.csv"
SAVE_FIG = "D:/文件/工作室/website/data/processed/P1.5_sci_signal_validation.png"
SAVE_TABLE = "D:/文件/工作室/website/data/processed/P1.5_signal_summary.csv"

# ===================== 加载数据 =====================
# 读取P1输出：仅上三角528个层对均值，严格对齐SCI定义
sci_mean = np.load(SCI_1D_PATH)
meta = pd.read_csv(META_PATH)
meta["SCI_mean"] = sci_mean

# 【优化1：分蛋白标准化DMS_score，消除跨蛋白尺度差异】
meta["DMS_score_z"] = meta.groupby("protein")["DMS_score"].transform(
    lambda x: (x - x.mean()) / x.std() if x.std() != 0 else 0
)

# 字段校验
required_cols = ["protein", "DMS_score", "DMS_score_z", "DMS_score_bin"]
for col in required_cols:
    if col not in meta.columns:
        raise ValueError(f"Metadata缺失关键字段：{col}")

protein_list = sorted(meta["protein"].unique())
print(f"检测到5个目标蛋白：{protein_list}")
print(f"全局总样本量：{len(meta)}")
protein_sample_count = {}
for p in protein_list:
    cnt = len(meta[meta["protein"] == p])
    protein_sample_count[p] = cnt
    print(f"  {p}: {cnt} 个突变位点")

# 自动获取所有蛋白最小样本量，作为均衡采样上限
min_protein_n = min(protein_sample_count.values())
print(f"\n均衡采样基准：取所有蛋白最小样本量 = {min_protein_n}")

# ===================== 统计工具函数（已修复所有变量Bug） =====================
def calc_signal_stats(df_sub, dataset_label):
    x = df_sub["SCI_mean"].values
    y_raw = df_sub["DMS_score"].values
    y_z = df_sub["DMS_score_z"].values
    y_bin = df_sub["DMS_score_bin"].values

    # 原始DMS相关（仅单蛋白有参考意义）
    pear_r_raw, pear_p_raw = pearsonr(x, y_raw)
    spear_r_raw, spear_p_raw = spearmanr(x, y_raw)
    # Z标准化DMS相关（全局/均衡数据集优先使用）
    pear_r_z, pear_p_z = pearsonr(x, y_z)
    spear_r_z, spear_p_z = spearmanr(x, y_z)

    # 双向AUROC
    auc_forward = roc_auc_score(y_bin, x)
    auc_reverse = roc_auc_score(y_bin, -x)
    if auc_reverse > auc_forward:
        auc_correct = auc_reverse
        x_opt = -x
        trend_desc = "SCI越大 → 致病/功能损伤越强"
    else:
        auc_correct = auc_forward
        x_opt = x
        trend_desc = "SCI越大 → 中性突变越强"

    # 【优化3：分开保存原始/校正AUPRC，和AUC对齐】
    auprc_forward = average_precision_score(y_bin, x)
    auprc_reverse = average_precision_score(y_bin, -x)
    auprc_correct = max(auprc_forward, auprc_reverse)

    return {
        "dataset": dataset_label,
        "n_sample": len(df_sub),
        # 原始DMS相关
        "pearson_r_raw": round(pear_r_raw, 4),
        "pearson_p_raw": pear_p_raw,
        "spearman_r_raw": round(spear_r_raw, 4),
        "spearman_p_raw": spear_p_raw,
        # Z标准化DMS相关（全局推荐）
        "pearson_r_z": round(pear_r_z, 4),
        "pearson_p_z": pear_p_z,
        "spearman_r_z": round(spear_r_z, 4),
        "spearman_p_z": spear_p_z,
        # 分类指标
        "AUROC_raw": round(auc_forward, 4),
        "AUROC_corrected": round(auc_correct, 4),
        "AUPRC_raw": round(auprc_forward, 4),
        "AUPRC_corrected": round(auprc_correct, 4),
        "SCI_trend": trend_desc
    }

# ===================== 均衡采样函数【优化2：分层抽样，维持bin比例】 =====================
def build_balanced_dataset(df, sample_limit, seed=42):
    np.random.seed(seed)
    subset_list = []
    for prot in protein_list:
        prot_df = df[df["protein"] == prot].copy()
        target_size = min(len(prot_df), sample_limit)
        # 分层抽样，保证0/1比例不变
        frac = target_size / len(prot_df)
        sample_df = prot_df.groupby("DMS_score_bin", group_keys=False).apply(
            lambda g: g.sample(frac=frac, random_state=seed)
        )
        subset_list.append(sample_df)
    balanced_df = pd.concat(subset_list, axis=0, ignore_index=True)
    return balanced_df

# 构造均衡分层数据集
balanced_meta = build_balanced_dataset(meta, sample_limit=min_protein_n)
print(f"均衡分层数据集构造完成，每个蛋白最多{min_protein_n}条，总样本：{len(balanced_meta)}")

# ===================== 批量计算所有分组指标 =====================
result_rows = []
# 1. 原始全量数据集
stats_all = calc_signal_stats(meta, "All_Proteins_Full")
result_rows.append(stats_all)

# 2. 均衡分层数据集（核心参考，Z分数为主要相关指标）
stats_balanced = calc_signal_stats(balanced_meta, "Balanced_Stratified")
result_rows.append(stats_balanced)

# 3. 分蛋白独立统计（最可信结果）
for prot_name in protein_list:
    single_df = meta[meta["protein"] == prot_name].copy()
    prot_stats = calc_signal_stats(single_df, prot_name)
    result_rows.append(prot_stats)

# 导出汇总表格
summary_df = pd.DataFrame(result_rows)
summary_df.to_csv(SAVE_TABLE, index=False, encoding="utf-8-sig")
print("\n========== SCI单值均值 与功能信号统计汇总表 ==========")
print("【关键解读提示】")
print("1. 跨蛋白混合数据集(All/Balanced)优先查看 DMS_score_z 标准化相关系数；")
print("2. 单蛋白分组可直接使用原始DMS_score相关，打分体系统一；")
print("3. SCI_mean为528层对全局平均，易稀释局部层对特异性信号，结果仅作前置筛查。")
print(summary_df.to_string(index=False))

# ===================== 绘图模块 =====================
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
fig = plt.figure(figsize=(16, 12))
grid = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.32)

# A 全局原始SCI分箱分布
ax1 = fig.add_subplot(grid[0, 0])
sns.boxplot(data=meta, x="DMS_score_bin", y="SCI_mean", palette="RdBu", ax=ax1)
ax1.set_xticklabels(["DMS_bin=0 Pathogenic", "DMS_bin=1 Neutral"])
ax1.set_title("A. Global SCI Distribution by Mutation Category", weight="bold")
ax1.set_xlabel("Mutation Bin")
ax1.set_ylabel("SCI (Upper Triangle Average)")
ax1.grid(axis="y", alpha=0.3)

# B 均衡数据集 SCI vs 标准化DMS（核心散点）
ax2 = fig.add_subplot(grid[0, 1])
sns.scatterplot(data=balanced_meta, x="SCI_mean", y="DMS_score_z", s=12, alpha=0.5, ax=ax2)
fit_coeff = np.polyfit(balanced_meta["SCI_mean"], balanced_meta["DMS_score_z"], 1)
fit_line = np.poly1d(fit_coeff)
ax2.plot(balanced_meta["SCI_mean"], fit_line(balanced_meta["SCI_mean"]), c="red", lw=2)
ax2.set_title(f"B. Balanced Stratified: SCI vs Z-DMS\nSpearman r_z = {stats_balanced['spearman_r_z']}", weight="bold")
ax2.set_xlabel("SCI Mean")
ax2.set_ylabel("Z-Normalized DMS Score")
ax2.grid(alpha=0.3)

# C 各组校正AUROC对比
ax3 = fig.add_subplot(grid[1, 0])
sns.barplot(data=summary_df, x="dataset", y="AUROC_corrected", ax=ax3, palette="viridis")
ax3.axhline(y=0.5, linestyle="--", c="red", label="Random Baseline AUC=0.5")
ax3.set_title("C. Corrected AUROC (All Groups)", weight="bold")
ax3.set_xlabel("Dataset Group")
ax3.set_ylabel("AUROC")
ax3.tick_params(axis='x', rotation=45)
ax3.legend()
ax3.grid(axis="y", alpha=0.3)

# D 各组标准化Spearman相关系数
ax4 = fig.add_subplot(grid[1, 1])
sns.barplot(data=summary_df, x="dataset", y="spearman_r_z", ax=ax4, palette="coolwarm")
ax4.axhline(y=0, linestyle="--", c="black")
ax4.set_title("D. Spearman Correlation (SCI ~ Z-DMS)", weight="bold")
ax4.set_xlabel("Dataset Group")
ax4.set_ylabel("Spearman r (Z-score DMS)")
ax4.tick_params(axis='x', rotation=45)
ax4.grid(axis="y", alpha=0.3)

plt.suptitle("P1.5 SCI Single Value Signal Check (Upper-Triangle Average, Z-Norm DMS)", fontsize=14, weight="bold", y=0.98)
plt.tight_layout()
plt.savefig(SAVE_FIG, dpi=300, bbox_inches="tight")
plt.show()

print(f"\n✅ 可视化图表已保存：{SAVE_FIG}")
print(f"✅ 统计汇总表格已保存：{SAVE_TABLE}")

# ===================== 分层结论输出 =====================
print("\n" + "="*80)
print("【关键前置说明】")
print("1. 当前仅检验528层对压缩后的一维SCI均值，大量局部层对特异性信号被平均稀释；")
print("2. 跨蛋白DMS实验打分尺度不一致，混合数据集仅采信Z标准化DMS相关；")
print("3. 无论本脚本结果强弱，均不判定SCI矩阵整体无效，跑完直接进入p0_layer_pair_mining。")
print("="*80)

p_threshold = 0.05
auc_min_threshold = 0.65
corr_min_threshold = 0.1

# 均衡分层数据集（核心判断依据，使用Z标准化指标）
bal_significant = stats_balanced["spearman_p_z"] < p_threshold
bal_auc_ok = stats_balanced["AUROC_corrected"] >= auc_min_threshold
bal_corr_ok = abs(stats_balanced["spearman_r_z"]) >= corr_min_threshold

# 全量数据集仅作对照
full_significant = stats_all["spearman_p_z"] < p_threshold
full_auc_ok = stats_all["AUROC_corrected"] >= auc_min_threshold
full_corr_ok = abs(stats_all["spearman_r_z"]) >= corr_min_threshold

print(f"【均衡分层数据集（消除样本偏差，Z-DMS为基准）】")
print(f"  Spearman r_z = {stats_balanced['spearman_r_z']}, p-value = {stats_balanced['spearman_p_z']:.2e}")
print(f"  校正AUROC = {stats_balanced['AUROC_corrected']}, 趋势：{stats_balanced['SCI_trend']}")
print(f"【原始全量数据集（PABP样本占比极高，仅对照）】")
print(f"  Spearman r_z = {stats_all['spearman_r_z']}, p-value = {stats_all['spearman_p_z']:.2e}")
print(f"  校正AUROC = {stats_all['AUROC_corrected']}, 趋势：{stats_all['SCI_trend']}")
print("-"*80)

# 三档解读，全部引导下一步层对挖掘
if bal_significant and bal_auc_ok and bal_corr_ok:
    print("🔥 情况A：一维平均SCI存在明显跨蛋白功能关联信号（标准化后）")
    print("解读：全局层间扰动均值即可区分突变损伤，后续层对筛选预期效果更佳；")
    print("行动：直接运行 p0_layer_pair_mining 挖掘高区分层对。")
elif bal_significant:
    print("✅ 情况B：一维平均SCI存在弱但统计显著信号")
    print("解读：全局均值信息量有限，局部特定层对会携带更强生物学信号；")
    print("行动：运行p0筛选核心层组合，使用多维度SCI矩阵分析。")
else:
    print("ℹ️ 情况C：一维平均SCI无显著全局关联信号")
    print("解读：528组层对平均稀释特异性扰动，不能否定33×33原始SCI矩阵价值；")
    print("行动：必须进入p0单独筛选、验证高贡献局部层对。")

print("\n统一执行路线：本脚本跑完，直接运行 p0_layer_pair_mining.py")
print("="*80)