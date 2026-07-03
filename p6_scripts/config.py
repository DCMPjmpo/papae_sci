"""
P6 多蛋白 SCI 扩展 — 全局配置
所有脚本统一从此读取，不再硬编码路径/蛋白名
"""

import os

# ========== 项目根目录 ==========
# config.py 在 p6_scripts/ 下，ROOT_DIR 要上跳一级到 website/
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ========== 数据目录 ==========
RAW_DATA_DIR = os.path.join(ROOT_DIR, "data", "p6_processed")
OUTPUT_DIR = os.path.join(ROOT_DIR, "data", "processed")
PLOT_DIR = os.path.join(ROOT_DIR, "p6_plots")
MODEL_DIR = os.path.join(ROOT_DIR, "model", "hub")

# 确保目录存在（只创建输出和图表目录）
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)
# MODEL_DIR 已有，不创建


# ============================================================
# 蛋白定义
# ============================================================
# 格式: (统一蛋白名, 数据集名, 文件名)
# 一个蛋白可由多个数据集组成（如 TEM1 有4个来源）
PROTEIN_DATASETS = [
    # PTEN
    ("PTEN", "PTEN_Mighell2018",    "PTEN_HUMAN_Mighell_2018.csv"),
    ("PTEN", "PTEN_Matreyek2021",   "PTEN_HUMAN_Matreyek_2021.csv"),

    # TEM-1 (BLAT_ECOLX)
    ("TEM1", "TEM1_Deng2012",       "BLAT_ECOLX_Deng_2012.csv"),
    ("TEM1", "TEM1_Firnberg2014",   "BLAT_ECOLX_Firnberg_2014.csv"),
    ("TEM1", "TEM1_Jacquier2013",   "BLAT_ECOLX_Jacquier_2013.csv"),
    ("TEM1", "TEM1_Stiffler2015",   "BLAT_ECOLX_Stiffler_2015.csv"),

    # CBS
    ("CBS",  "CBS_Sun2020",         "CBS_HUMAN_Sun_2020.csv"),

    # GAL4
    ("GAL4", "GAL4_Kitzman2015",    "GAL4_YEAST_Kitzman_2015.csv"),

    # PABP (含双点突变)
    ("PABP", "PABP_Melamed2013",    "PABP_YEAST_Melamed_2013.csv"),
]

# 蛋白列表（去重）
PROTEINS = sorted(set(p for p, _, _ in PROTEIN_DATASETS))

# 排除列表（如 BRCA1 — 序列过长 >1022aa）
EXCLUDED_PROTEINS = []

# ============================================================
# 多点突变处理
# ============================================================
EXPAND_MULTI = True              # 是否将多点突变展开为多条单点记录
MULTI_SEPARATOR = ":"            # 多点突变分隔符

# ============================================================
# ESM2 模型配置
# ============================================================
ESM_MODEL_NAME = "esm2_t33_650M_UR50D"
ESM_N_LAYERS = 33
ESM_DIM = 1280
ESM_MAX_SEQ_LEN = 1022           # ESM2-650M 最大序列长度
DEVICE = "cuda"                  # 或 "cpu"
BATCH_SIZE = 4

# ============================================================
# 输出文件前缀（统一命名）
# ============================================================
PREFIX = "all_proteins"

# 中间产物路径
MERGE_CSV       = os.path.join(OUTPUT_DIR, f"{PREFIX}.csv")
MERGE_WITH_WT   = os.path.join(OUTPUT_DIR, f"{PREFIX}_with_wt.csv")
WT_SEQUENCES    = os.path.join(OUTPUT_DIR, "wt_sequences.csv")

MERGE_EXPANDED = os.path.join(OUTPUT_DIR, f"{PREFIX}_expanded.csv")

# Phase 2 输出 (global delta)
MUT_EMBEDDINGS  = os.path.join(OUTPUT_DIR, f"{PREFIX}_mut_embeddings.npy")
WT_EMBEDDINGS   = os.path.join(OUTPUT_DIR, f"{PREFIX}_wt_embeddings_dict.npy")
DELTA_EMBEDDINGS = os.path.join(OUTPUT_DIR, f"{PREFIX}_delta_embeddings.npy")
GLOBAL_METADATA  = os.path.join(OUTPUT_DIR, f"{PREFIX}_metadata.csv")

# Phase 2 输出 (site-level)
DELTA_SITE      = os.path.join(OUTPUT_DIR, f"{PREFIX}_delta_site.npy")
SCI_MATRICES    = os.path.join(OUTPUT_DIR, f"{PREFIX}_sci_site_matrices.npy")
SCI_SCORES_MEAN = os.path.join(OUTPUT_DIR, f"{PREFIX}_sci_site_scores_mean.npy")
SCI_SCORES_TOP20 = os.path.join(OUTPUT_DIR, f"{PREFIX}_sci_site_scores_top20.npy")
SCI_SCORES_TOP50 = os.path.join(OUTPUT_DIR, f"{PREFIX}_sci_site_scores_top50.npy")
SITE_METADATA   = os.path.join(OUTPUT_DIR, f"{PREFIX}_site_metadata.csv")

# Phase 3 输出
SCI_SCORES      = os.path.join(OUTPUT_DIR, f"{PREFIX}_sci_scores.npy")

# ============================================================
# 分层采样配置（可选：控制数据量）
# ============================================================
MAX_SAMPLES_PER_PROTEIN_LABEL = None  # None = 不采样，用全部数据；设 500 = 每蛋白每标签最多500条

# ============================================================
# 打印配置摘要
# ============================================================
def print_config():
    print("=" * 60)
    print("P6 多蛋白 SCI 扩展 — 全局配置")
    print("=" * 60)
    print(f"  蛋白数: {len(PROTEINS)}")
    for p in PROTEINS:
        datasets = [d for prot, d, _ in PROTEIN_DATASETS if prot == p]
        print(f"    {p}: {len(datasets)} 个数据集 {datasets}")
    print(f"  多点展开: {EXPAND_MULTI} (分隔符='{MULTI_SEPARATOR}')")
    print(f"  ESM2 模型: {ESM_MODEL_NAME} ({ESM_N_LAYERS}层, {ESM_DIM}维)")
    print(f"  设备: {DEVICE}")
    print(f"  输出前缀: {PREFIX}")
    print(f"  输出目录: {OUTPUT_DIR}")
    print(f"  最大序列长度: {ESM_MAX_SEQ_LEN}")
    print("=" * 60)


if __name__ == "__main__":
    print_config()
