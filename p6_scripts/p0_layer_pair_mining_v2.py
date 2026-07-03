import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
import os

# ====================== 【配置区】 ======================
# [旧·服务器路径]
# MATRIX_PATH = "/mnt/sda/gws_1020251255/data/processed/all_proteins_sci_site_matrices.dat"
# META_CSV_PATH = "/mnt/sda/gws_1020251255/data/processed/all_proteins_site_metadata.csv"
# OUTPUT_DIR = "/mnt/sda/gws_1020251255/data/processed/p0_output"
# [新·本地路径]
MATRIX_PATH = "D:/文件/工作室/website/data/processed/all_proteins_sci_site_matrices.dat"
META_CSV_PATH = "D:/文件/工作室/website/data/processed/all_proteins_site_metadata.csv"
OUTPUT_DIR = "D:/文件/工作室/website/data/p0_output"
N_LAYERS = 33
TOP_N = 50
# 修正：下标0~32分层边界
EARLY_END = 10
MIDDLE_END = 21
PROTEIN_LIST = ["All", "CBS", "GAL4", "PABP", "PTEN", "TEM1"]
WEIGHT_AUC = 0.5
WEIGHT_SPEARMAN_ABS = 0.5
# =======================================================

def cohen_d(x, y):
    nx = len(x)
    ny = len(y)
    mean_x = np.mean(x)
    mean_y = np.mean(y)
    var_x = np.var(x, ddof=1)
    var_y = np.var(y, ddof=1)
    pooled_std = np.sqrt(((nx-1)*var_x + (ny-1)*var_y) / (nx + ny - 2))
    if pooled_std == 0:
        return 0.0
    return (mean_x - mean_y) / pooled_std

def get_group(l):
    """判断单层归属 E/M/L，输入0~32下标"""
    if l <= EARLY_END:
        return "E"
    elif l <= MIDDLE_END:
        return "M"
    else:
        return "L"

def get_pair_category(i, j):
    """生成六类配对标签 EE/EM/EL/MM/ML/LL"""
    g1 = get_group(i)
    g2 = get_group(j)
    sorted_pair = sorted([g1, g2])
    return "".join(sorted_pair)

def compute_all_metrics(scores, y_bin, dms_z_vals):
    """
    指标计算：Spearman 传入 Z归一化DMS
    复合得分做值域对齐，消除AUC与相关系数量级差异
    """
    auc_forward = roc_auc_score(y_bin, scores)
    auroc_abs = max(auc_forward, 1 - auc_forward)

    if auc_forward >= 0.5:
        direction = "Neutral ↑ when SCI ↑"
    else:
        direction = "Pathogenic ↑ when SCI ↑"

    # 修正：使用Z-score DMS做相关性计算
    spr_r, spr_p = spearmanr(scores, dms_z_vals)

    idx_patho = (y_bin == 0)
    idx_neutral = (y_bin == 1)
    s_path = scores[idx_patho]
    s_neut = scores[idx_neutral]
    mean_path = np.mean(s_path)
    mean_neut = np.mean(s_neut)
    mean_diff = mean_path - mean_neut
    cd = cohen_d(s_path, s_neut)

    # AUC归一化到 [0, 1]
    auc_norm = (auroc_abs - 0.5) / 0.5
    composite_score = WEIGHT_AUC * auc_norm + WEIGHT_SPEARMAN_ABS * abs(spr_r)

    return {
        "AUROC_raw": auc_forward,
        "AUROC_abs": auroc_abs,
        "AUROC_normalized": auc_norm,
        "Direction": direction,
        "Spearman_r": spr_r,
        "Spearman_p": spr_p,
        "Cohen_d": cd,
        "pathogenic_mean": mean_path,
        "neutral_mean": mean_neut,
        "path_minus_neut_diff": mean_diff,
        "CompositeScore": composite_score
    }

def calc_freq_and_category(df_top, subset_name, sort_tag, out_subdir):
    """
    复用函数：对任意Top50表格统计
    1. 层出现频次
    2. 六类层对占比
    输出独立文件，三套排序分开统计
    """
    tag = f"_{sort_tag}"
    # 层频次统计
    all_layers = list(df_top["layer_i"]) + list(df_top["layer_j"])
    freq_counter = Counter(all_layers)
    freq_df = pd.DataFrame(list(freq_counter.items()), columns=["Layer_idx_1based", "Occurrence_Freq"])
    freq_df = freq_df.sort_values("Occurrence_Freq", ascending=False).reset_index(drop=True)
    freq_csv = os.path.join(out_subdir, f"{subset_name}_layer_frequency{tag}.csv")
    freq_df.to_csv(freq_csv, index=False)

    # 频次柱状图
    plt.figure(figsize=(10,5), dpi=150)
    plt.bar(freq_df["Layer_idx_1based"].astype(str), freq_df["Occurrence_Freq"], color="#4472c4")
    plt.xlabel("Layer Index (1~33, ESM convention)")
    plt.ylabel("Occurrence Count in Top50 Pairs")
    plt.title(f"{subset_name}: Layer Frequency | Top50 {sort_tag}")
    plt.xticks(rotation=90, fontsize=7)
    plt.tight_layout()
    fig_path = os.path.join(out_subdir, f"{subset_name}_layer_frequency{tag}.png")
    plt.savefig(fig_path, bbox_inches="tight")
    plt.close()

    # 早中晚期总占比
    total_freq = sum(freq_counter.values())
    early_total = sum(v for k,v in freq_counter.items() if (k-1) <= EARLY_END)
    mid_total = sum(v for k,v in freq_counter.items() if EARLY_END < (k-1) <= MIDDLE_END)
    late_total = sum(v for k,v in freq_counter.items() if (k-1) > MIDDLE_END)
    ratio_early = early_total / total_freq if total_freq>0 else 0
    ratio_mid = mid_total / total_freq if total_freq>0 else 0
    ratio_late = late_total / total_freq if total_freq>0 else 0

    # 六类配对占比
    cat_counter = Counter(df_top["Pair_Category"])
    cat_total = sum(cat_counter.values())
    cat_ratio = {}
    for cat in ["EE","EM","EL","MM","ML","LL"]:
        cnt_cat = cat_counter.get(cat, 0)
        cat_ratio[cat] = cnt_cat / cat_total if cat_total>0 else 0
    cat_df = pd.DataFrame(list(cat_ratio.items()), columns=["Pair_Type", "Ratio"])
    cat_df["Count"] = cat_df["Pair_Type"].map(cat_counter)
    cat_csv = os.path.join(out_subdir, f"{subset_name}_pair_category_ratio{tag}.csv")
    cat_df.to_csv(cat_csv, index=False)

    return (ratio_early, ratio_mid, ratio_late), cat_df

def process_single_subset(sci_mat, df_meta, subset_name, out_subdir):
    # 样本筛选，列名使用 protein
    if subset_name == "All":
        mask = np.ones(len(df_meta), dtype=bool)
    else:
        mask = (df_meta["protein"] == subset_name)
    df_sub = df_meta.loc[mask].reset_index(drop=True)
    mat_sub = sci_mat[mask, :, :]
    n_sample = len(df_sub)
    if n_sample < 10:
        print(f"警告：{subset_name}样本过少({n_sample})，跳过")
        return None, None

    y_bin = df_sub["DMS_score_bin"].values
    dms_z = df_sub["DMS_score_z"].values
    total_pairs = int(N_LAYERS * (N_LAYERS - 1) / 2)
    pair_results = []
    cnt = 0

    print(f"\n==== Processing subset: {subset_name}, samples={n_sample}, total pairs={total_pairs} ====")
    for i in range(N_LAYERS):
        for j in range(i + 1, N_LAYERS):
            cnt += 1
            pair_score = mat_sub[:, i, j]
            metrics = compute_all_metrics(pair_score, y_bin, dms_z)
            pair_type = get_pair_category(i, j)
            # 修正：输出层编号+1，转为论文常用1~33编号
            row = {
                "layer_i_0based": i,
                "layer_j_0based": j,
                "layer_i": i + 1,
                "layer_j": j + 1,
                "Pair_Category": pair_type
            }
            row.update(metrics)
            pair_results.append(row)
            if cnt % 50 == 0:
                print(f"{subset_name}: processed {cnt}/{total_pairs}")

    df_all_pairs = pd.DataFrame(pair_results)
    csv_all = os.path.join(out_subdir, f"{subset_name}_layer_pair_metrics.csv")
    df_all_pairs.to_csv(csv_all, index=False)

    # 生成三套独立Top50
    df_top50_auc = df_all_pairs.sort_values("AUROC_abs", ascending=False).head(TOP_N).reset_index(drop=True)
    df_top50_spear = df_all_pairs.sort_values("Spearman_r", key=abs, ascending=False).head(TOP_N).reset_index(drop=True)
    df_top50_comp = df_all_pairs.sort_values("CompositeScore", ascending=False).head(TOP_N).reset_index(drop=True)

    df_top50_auc.to_csv(os.path.join(out_subdir, f"top50_auroc_{subset_name}.csv"), index=False)
    df_top50_spear.to_csv(os.path.join(out_subdir, f"top50_spearman_{subset_name}.csv"), index=False)
    df_top50_comp.to_csv(os.path.join(out_subdir, f"top50_composite_{subset_name}.csv"), index=False)

    # 三套分别统计频次+配对占比
    seg_ratio_auc, cat_auc = calc_freq_and_category(df_top50_auc, subset_name, "AUC", out_subdir)
    seg_ratio_spear, cat_spear = calc_freq_and_category(df_top50_spear, subset_name, "Spearman", out_subdir)
    seg_ratio_comp, cat_comp = calc_freq_and_category(df_top50_comp, subset_name, "Composite", out_subdir)

    print(f"\n[{subset_name}] Segment Ratio (AUC Top50): Early={seg_ratio_auc[0]:.2%}, Mid={seg_ratio_auc[1]:.2%}, Late={seg_ratio_auc[2]:.2%}")
    print(f"[{subset_name}] Segment Ratio (Spearman Top50): Early={seg_ratio_spear[0]:.2%}, Mid={seg_ratio_spear[1]:.2%}, Late={seg_ratio_spear[2]:.2%}")
    print(f"[{subset_name}] Segment Ratio (Composite Top50): Early={seg_ratio_comp[0]:.2%}, Mid={seg_ratio_comp[1]:.2%}, Late={seg_ratio_comp[2]:.2%}")

    # AUROC_abs 对称热力图（图内显示1~33层）
    auroc_mat = np.full((N_LAYERS, N_LAYERS), np.nan)
    for row in df_all_pairs.itertuples():
        i0, j0, auc = row.layer_i_0based, row.layer_j_0based, row.AUROC_abs
        auroc_mat[i0, j0] = auc
        auroc_mat[j0, i0] = auc
    plt.figure(figsize=(12,10), dpi=150)
    im = plt.imshow(auroc_mat, cmap="RdBu_r", vmin=0.5, vmax=0.8)
    plt.colorbar(im, label="Absolute Corrected AUROC")
    plt.xticks(range(N_LAYERS), [x+1 for x in range(N_LAYERS)], fontsize=7)
    plt.yticks(range(N_LAYERS), [x+1 for x in range(N_LAYERS)], fontsize=7)
    plt.xlabel("Layer j (1~33)")
    plt.ylabel("Layer i (1~33)")
    plt.title(f"{subset_name} AUROC_abs Heatmap of All Layer-Pair SCI Scores", fontsize=11)
    plt.tight_layout()
    heatmap_path = os.path.join(out_subdir, f"{subset_name}_auroc_heatmap.png")
    plt.savefig(heatmap_path, bbox_inches="tight")
    plt.close()

    return df_all_pairs, df_top50_comp, seg_ratio_comp

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    # 加载memmap矩阵
    print("Loading SCI memmap matrix (95142, 33, 33) ...")
    sci_mat = np.memmap(
        MATRIX_PATH,
        dtype=np.float32,
        mode="r",
        shape=(95142, N_LAYERS, N_LAYERS)
    )
    print(f"Matrix shape loaded: {sci_mat.shape}")

    # 读取元数据
    print("Loading metadata csv...")
    df_meta = pd.read_csv(META_CSV_PATH)
    print("Metadata columns:", df_meta.columns.tolist())
    required_raw_cols = ["protein", "DMS_score_bin", "DMS_score"]
    for col in required_raw_cols:
        if col not in df_meta.columns:
            raise KeyError(f"缺失必要列：{col}")

    # 按需生成 按蛋白分组Z归一化DMS，后续Spearman专用
    if "DMS_score_z" not in df_meta.columns:
        print("DMS_score_z not found, generate per-protein Z-score normalization...")
        df_meta["DMS_score_z"] = df_meta.groupby("protein")["DMS_score"].transform(
            lambda x: (x - x.mean()) / x.std()
        )

    assert sci_mat.shape[0] == len(df_meta), "样本行数不匹配：矩阵行数 ≠ 元数据行数"

    protein_top_pair_set = dict()
    summary_ratio_list = []

    for prot in PROTEIN_LIST:
        sub_dir = os.path.join(OUTPUT_DIR, prot)
        os.makedirs(sub_dir, exist_ok=True)
        res = process_single_subset(sci_mat, df_meta, prot, sub_dir)
        if res is not None:
            df_all, df_t50_comp, seg_ratios = res
            summary_ratio_list.append({
                "Protein_Subset": prot,
                "Ratio_Early_Composite": seg_ratios[0],
                "Ratio_Middle_Composite": seg_ratios[1],
                "Ratio_Late_Composite": seg_ratios[2]
            })
            pair_set = set(zip(df_t50_comp["layer_i_0based"], df_t50_comp["layer_j_0based"]))
            protein_top_pair_set[prot] = pair_set

    # 汇总全局分层占比总表
    df_summary_ratio = pd.DataFrame(summary_ratio_list)
    df_summary_ratio.to_csv(os.path.join(OUTPUT_DIR, "layer_segment_ratio_summary.csv"), index=False)
    print("\n=== All subsets calculation finished ===")
    print(df_summary_ratio.round(4))

    # 跨蛋白普适层对统计（论文亮点）
    print("\n==== Cross-protein universal layer-pair analysis ====")
    single_prot_names = [p for p in PROTEIN_LIST if p != "All"]
    pair_occurrence_count = defaultdict(int)
    for p_name in single_prot_names:
        pairs = protein_top_pair_set[p_name]
        for pr in pairs:
            pair_occurrence_count[pr] += 1

    cross_pair_df = pd.DataFrame(list(pair_occurrence_count.items()), columns=["Layer_Pair_0based", "Cross_Protein_Count"])
    cross_pair_df[["layer_i_0", "layer_j_0"]] = pd.DataFrame(cross_pair_df["Layer_Pair_0based"].tolist(), index=cross_pair_df.index)
    cross_pair_df["layer_i_1"] = cross_pair_df["layer_i_0"] + 1
    cross_pair_df["layer_j_1"] = cross_pair_df["layer_j_0"] + 1
    cross_pair_df = cross_pair_df.sort_values("Cross_Protein_Count", ascending=False).reset_index(drop=True)
    cross_csv_path = os.path.join(OUTPUT_DIR, "cross_protein_universal_pairs.csv")
    cross_pair_df.to_csv(cross_csv_path, index=False)
    print(f"Cross-protein pair statistics saved to {cross_csv_path}")
    print("Top recurring cross-protein layer pairs:")
    print(cross_pair_df.head(10))

if __name__ == "__main__":
    main()