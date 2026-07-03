"""
p5_permutation.py  (v3.1 — fig_D hotfix)
=========================================

修复了 v3 中 fig_D (QQ significance plot) 因画布尺寸爆炸导致的
ValueError: Image size ... too large.

变更点（仅 fig_D_qq_significance）：
  1. 过滤 NaN/Inf z-score，避免 matplotlib 坐标非有限值错误
  2. figsize 从 (6.6, 6.2) 增大到 (10, 10)
  3. 标签采用左右交替 + 颜色区分，减少重叠
  4. 移除 fig.tight_layout()；依赖 rcParams 中已有的 savefig.bbox='tight'
  5. 绘图边界 lim 放大 1.15 倍，给右侧标签留足空间

其余所有代码与 v3 完全一致。
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy import stats as sps


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
ROOT = Path("D:/文件/工作室/website")
DATA_DIR = ROOT / "data" / "processed"
OUT_DIR = ROOT / "data" / "p0_output" / "p5_permutation"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PROTEINS = ["CBS", "GAL4", "PABP", "PTEN", "TEM1"]
N_LAYERS = 33
N_MUT_EXPECTED = 95142
N_PAIRS_EXPECTED = 528
TOP_K = 50

PERM_STAGES = [1000, 5000, 10000]
PERM_STOP_P = 0.001
N_PERM_PROT = 1000

N_BOOT_T1 = 1000
N_BOOT_T23 = 1000
N_BOOT_VALIDATION = 100

RNG_SEED = 42

PRIMARY_STATS = ["T1", "T2", "T3", "T3_std"]
ALL_STATS = ["T1", "T2", "T3", "T3_std", "T_max"]
N_BONFERRONI_TESTS = len(PROTEINS) * len(PRIMARY_STATS)

EARLY = set(range(1, 9))
MIDDLE = set(range(9, 25))
LATE = set(range(25, 34))


def layer_band(L: int) -> str:
    if L in EARLY:
        return "E"
    if L in MIDDLE:
        return "M"
    return "L"


def pair_category(i: int, j: int) -> str:
    return "".join(sorted([layer_band(i), layer_band(j)]))


PAIR_LAYERS: List[Tuple[int, int]] = [
    (i, j) for i in range(1, N_LAYERS + 1) for j in range(i + 1, N_LAYERS + 1)
]
PAIR_CATEGORIES: List[str] = [pair_category(i, j) for (i, j) in PAIR_LAYERS]
PAIR_CATEGORIES_ARR = np.array(PAIR_CATEGORIES)
N_PAIRS = len(PAIR_LAYERS)
EL_MASK = PAIR_CATEGORIES_ARR == "EL"
N_EL = int(EL_MASK.sum())
N_NONEL = N_PAIRS - N_EL
assert N_PAIRS == N_PAIRS_EXPECTED

STAT_LABELS = {
    "T1":     r"$T_1$ : SCI mean diff",
    "T2":     fr"$T_2$ : # EL pairs in Top{TOP_K}",
    "T3":     r"$T_3$ : EL dominance |AUROC-0.5|",
    "T3_std": r"$T_3^{std}$ : standardised EL dominance",
    "T_max":  r"$T_{max}$ : max |AUROC-0.5| over 528 pairs",
}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "axes.linewidth": 0.9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8,
    "legend.frameon": False,
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_inputs() -> Tuple[np.memmap, np.ndarray, pd.DataFrame]:
    sci_full_path = DATA_DIR / "all_proteins_sci_site_matrices.dat"
    sci_top50_path = DATA_DIR / "all_proteins_sci_site_scores_top50.npy"
    meta_path = DATA_DIR / "all_proteins_site_metadata.csv"

    for p in (sci_full_path, sci_top50_path, meta_path):
        if not p.exists():
            raise FileNotFoundError(p)

    sci_full = np.memmap(sci_full_path, dtype=np.float32, mode="r")
    sci_full = sci_full.reshape(N_MUT_EXPECTED, N_LAYERS, N_LAYERS)
    sci_top50 = np.load(sci_top50_path, mmap_mode="r")
    meta = pd.read_csv(meta_path)

    if "protein" not in meta.columns or "DMS_score_bin" not in meta.columns:
        raise ValueError("metadata missing \'protein\' or \'DMS_score_bin\'")
    if len(meta) != N_MUT_EXPECTED or len(sci_top50) != N_MUT_EXPECTED:
        raise ValueError("row-count mismatch")

    return sci_full, np.asarray(sci_top50, dtype=np.float32), meta


def materialise_pair_matrix(sci_full: np.memmap) -> np.ndarray:
    pair_idx_i = np.array([i - 1 for (i, j) in PAIR_LAYERS])
    pair_idx_j = np.array([j - 1 for (i, j) in PAIR_LAYERS])
    n = sci_full.shape[0]
    out = np.empty((n, N_PAIRS), dtype=np.float32)
    chunk = 4096
    for start in range(0, n, chunk):
        stop = min(start + chunk, n)
        block = np.asarray(sci_full[start:stop])
        out[start:stop] = block[:, pair_idx_i, pair_idx_j]
    return out


# ---------------------------------------------------------------------------
# AUROC via rank trick
# ---------------------------------------------------------------------------
def precompute_ranks(pair_matrix: np.ndarray) -> np.ndarray:
    n, k = pair_matrix.shape
    ranks = np.empty_like(pair_matrix, dtype=np.float32)
    for col in range(k):
        ranks[:, col] = sps.rankdata(pair_matrix[:, col],
                                     method="average").astype(np.float32)
    return ranks


def aurocs_from_ranks(ranks: np.ndarray, pos_mask: np.ndarray) -> np.ndarray:
    n_pos = int(pos_mask.sum())
    n_neg = ranks.shape[0] - n_pos
    if n_pos == 0 or n_neg == 0:
        return np.full(ranks.shape[1], np.nan)
    sum_pos_ranks = ranks[pos_mask].sum(axis=0, dtype=np.float64)
    aurocs = (sum_pos_ranks - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
    return aurocs


# ---------------------------------------------------------------------------
# Test statistics
# ---------------------------------------------------------------------------
def stat_T1(sci_top50: np.ndarray, labels: np.ndarray) -> float:
    pos = labels == 1
    if pos.sum() == 0 or (~pos).sum() == 0:
        return np.nan
    return float(sci_top50[pos].mean() - sci_top50[~pos].mean())


def stat_T2(aurocs: np.ndarray) -> float:
    abs_auc = np.abs(aurocs - 0.5)
    top_idx = np.argpartition(-abs_auc, TOP_K - 1)[:TOP_K]
    return float(EL_MASK[top_idx].sum())


def stat_T3(aurocs: np.ndarray) -> float:
    abs_auc = np.abs(aurocs - 0.5)
    return float(abs_auc[EL_MASK].mean() - abs_auc[~EL_MASK].mean())


def stat_T3_std(aurocs: np.ndarray) -> float:
    abs_auc = np.abs(aurocs - 0.5)
    mu_non = abs_auc[~EL_MASK].mean()
    sd_non = abs_auc[~EL_MASK].std(ddof=1)
    sd_non = max(sd_non, 1e-8)
    if np.isnan(sd_non):
        return np.nan
    return float((abs_auc[EL_MASK].mean() - mu_non) / sd_non)


def stat_T_max(aurocs: np.ndarray) -> float:
    return float(np.nanmax(np.abs(aurocs - 0.5)))


def compute_all_stats(sci_top50: np.ndarray, labels: np.ndarray,
                      ranks: np.ndarray) -> Dict[str, float]:
    aurocs = aurocs_from_ranks(ranks, labels == 1)
    return {
        "T1":     stat_T1(sci_top50, labels),
        "T2":     stat_T2(aurocs),
        "T3":     stat_T3(aurocs),
        "T3_std": stat_T3_std(aurocs),
        "T_max":  stat_T_max(aurocs),
    }


# ---------------------------------------------------------------------------
# Permutation infrastructure
# ---------------------------------------------------------------------------
def within_protein_shuffle(labels: np.ndarray, protein_groups: List[np.ndarray],
                           rng: np.random.Generator) -> np.ndarray:
    shuffled = labels.copy()
    for idx in protein_groups:
        shuffled[idx] = rng.permutation(labels[idx])
    return shuffled


def empirical_p(t_null: np.ndarray, t_obs: float, two_sided: bool) -> float:
    t_null = np.asarray(t_null, dtype=np.float64)
    if two_sided:
        mu = float(np.nanmean(t_null))
        extreme = np.abs(t_null - mu) >= np.abs(t_obs - mu)
    else:
        extreme = t_null >= t_obs
    return (1.0 + int(np.nansum(extreme))) / (len(t_null) + 1.0)


def summarise(t_null: np.ndarray, t_obs: float, two_sided: bool,
              n_perm_effective: int) -> Dict[str, float]:
    t_null = np.asarray(t_null, dtype=np.float64)
    mu = float(np.nanmean(t_null))
    sd = float(np.nanstd(t_null, ddof=1))
    z = (t_obs - mu) / sd if sd > 0 else np.nan
    cohens_d = z
    p = empirical_p(t_null, t_obs, two_sided)
    ci95 = np.nanpercentile(t_null, [2.5, 97.5]).tolist()
    ci99 = np.nanpercentile(t_null, [0.5, 99.5]).tolist()
    return {
        "T_obs": float(t_obs),
        "null_mean": mu,
        "null_std": sd,
        "z_score": z,
        "cohens_d": cohens_d,
        "p_empirical": p,
        "effect_raw": float(t_obs - mu),
        "null_ci95_low": ci95[0],
        "null_ci95_high": ci95[1],
        "null_ci99_low": ci99[0],
        "null_ci99_high": ci99[1],
        "n_perm_effective": n_perm_effective,
        "two_sided": two_sided,
    }


# ---------------------------------------------------------------------------
# Adaptive permutation runner
# ---------------------------------------------------------------------------
@dataclass
class RunCfg:
    label: str
    two_sided: Dict[str, bool]
    stages: List[int]
    stop_p: float
    n_perm_fixed: int = None


def _fmt_eta(elapsed: float, done: int, total: int) -> str:
    if done == 0:
        return "?"
    rate = done / elapsed
    remain = (total - done) / rate
    return f"{remain/60:.1f} min"


def run_permutation_block(
    ranks: np.ndarray,
    sci_top50: np.ndarray,
    labels: np.ndarray,
    protein_groups: List[np.ndarray] | None,
    obs: Dict[str, float],
    cfg: RunCfg,
    rng: np.random.Generator,
) -> Tuple[Dict[str, np.ndarray], int]:
    stages = ([cfg.n_perm_fixed] if cfg.n_perm_fixed is not None else cfg.stages)
    n_target_final = stages[-1]

    null_acc = {s: [] for s in ALL_STATS}
    t0 = time.time()
    n_so_far = 0
    effective_n = stages[-1]

    for stage_idx, target_n in enumerate(stages):
        n_to_run = target_n - n_so_far
        for k in range(n_to_run):
            if protein_groups is not None:
                y_perm = within_protein_shuffle(labels, protein_groups, rng)
            else:
                y_perm = rng.permutation(labels)
            s = compute_all_stats(sci_top50, y_perm, ranks)
            for stat in ALL_STATS:
                null_acc[stat].append(s[stat])
            n_so_far += 1

            if n_so_far % 100 == 0 or n_so_far == target_n:
                p_running = {
                    stat: empirical_p(np.asarray(null_acc[stat]),
                                       obs[stat], cfg.two_sided[stat])
                    for stat in PRIMARY_STATS
                }
                p_max = max(p_running.values())
                elapsed = time.time() - t0
                eta = _fmt_eta(elapsed, n_so_far, n_target_final)
                print(f"      [{cfg.label}] perm={n_so_far:>5d}/{n_target_final}  "
                      f"p_max~{p_max:.4f}  elapsed={elapsed/60:.1f}min  ETA={eta}")

        if cfg.n_perm_fixed is None:
            p_running = {
                stat: empirical_p(np.asarray(null_acc[stat]),
                                   obs[stat], cfg.two_sided[stat])
                for stat in PRIMARY_STATS
            }
            p_max = max(p_running.values())
            print(f"      [{cfg.label}] stage {stage_idx+1}/{len(stages)} done "
                  f"(n={n_so_far}): p_max = {p_max:.4f} "
                  f"({'STOP' if p_max < cfg.stop_p else 'continue'})")
            if p_max < cfg.stop_p:
                effective_n = n_so_far
                break
    else:
        effective_n = n_so_far

    null_dist = {s: np.asarray(null_acc[s], dtype=np.float64) for s in ALL_STATS}
    return null_dist, effective_n


# ---------------------------------------------------------------------------
# Bootstrap
# ---------------------------------------------------------------------------
def bootstrap_T1(sci_top50: np.ndarray, labels: np.ndarray,
                 n_boot: int, rng: np.random.Generator) -> np.ndarray:
    N = len(labels)
    out = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, N, N)
        out[b] = stat_T1(sci_top50[idx], labels[idx])
    return out


def bootstrap_pair_stats(ranks: np.ndarray, labels: np.ndarray,
                         n_boot: int, rng: np.random.Generator,
                         label: str = "boot") -> Dict[str, np.ndarray]:
    N, K = ranks.shape
    out = {s: np.empty(n_boot) for s in ["T2", "T3", "T3_std", "T_max"]}
    t0 = time.time()
    for b in range(n_boot):
        idx = rng.integers(0, N, N)
        y = labels[idx]
        pos_mask = y == 1
        n_pos = int(pos_mask.sum())
        n_neg = N - n_pos
        if n_pos == 0 or n_neg == 0:
            for s in out:
                out[s][b] = np.nan
            continue
        pos_orig_idx = idx[pos_mask]
        sum_pos_ranks = ranks[pos_orig_idx].sum(axis=0, dtype=np.float64)
        aurocs = (sum_pos_ranks - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
        out["T2"][b] = stat_T2(aurocs)
        out["T3"][b] = stat_T3(aurocs)
        out["T3_std"][b] = stat_T3_std(aurocs)
        out["T_max"][b] = stat_T_max(aurocs)
        if (b + 1) % 100 == 0 or (b + 1) == n_boot:
            elapsed = time.time() - t0
            eta = _fmt_eta(elapsed, b + 1, n_boot)
            print(f"      [{label}] boot={b+1:>4d}/{n_boot}  "
                  f"elapsed={elapsed/60:.1f}min  ETA={eta}")
    return out


def bootstrap_validation_exact_vs_fast(
    pair_matrix: np.ndarray,
    ranks: np.ndarray,
    labels: np.ndarray,
    n_valid: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    N, K = pair_matrix.shape
    records = []
    t0 = time.time()

    for b in range(n_valid):
        idx = rng.integers(0, N, N)
        y = labels[idx]
        pos_mask = y == 1
        n_pos = int(pos_mask.sum())
        n_neg = N - n_pos
        if n_pos == 0 or n_neg == 0:
            continue

        sum_pos_exact = np.empty(K, dtype=np.float64)
        for col in range(K):
            rcol = sps.rankdata(pair_matrix[idx, col], method="average")
            sum_pos_exact[col] = rcol[pos_mask].sum()
        auroc_exact = (sum_pos_exact - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)

        sum_pos_fast = ranks[idx[pos_mask]].sum(axis=0, dtype=np.float64)
        auroc_fast = (sum_pos_fast - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)

        rec = {"boot_id": b + 1}
        for stat_name, stat_fn in [("T2", stat_T2), ("T3", stat_T3), ("T3_std", stat_T3_std)]:
            v_exact = stat_fn(auroc_exact)
            v_fast = stat_fn(auroc_fast)
            rec[f"{stat_name}_exact"] = v_exact
            rec[f"{stat_name}_fast"] = v_fast
            rec[f"{stat_name}_diff"] = v_fast - v_exact
        records.append(rec)

        if (b + 1) % 10 == 0 or (b + 1) == n_valid:
            elapsed = time.time() - t0
            eta = _fmt_eta(elapsed, b + 1, n_valid)
            print(f"      [validation] {b+1:>3d}/{n_valid}  "
                  f"elapsed={elapsed/60:.1f}min  ETA={eta}")

    df = pd.DataFrame(records)
    summary = {"boot_id": "SUMMARY"}
    for stat_name in ["T2", "T3", "T3_std"]:
        diffs = df[f"{stat_name}_diff"].astype(float)
        summary[f"{stat_name}_exact"] = np.nan
        summary[f"{stat_name}_fast"] = np.nan
        summary[f"{stat_name}_diff"] = (
            f"mean|diff|={diffs.abs().mean():.4g}; "
            f"max|diff|={diffs.abs().max():.4g}; "
            f"rmse={np.sqrt((diffs**2).mean()):.4g}"
        )
    df = pd.concat([df, pd.DataFrame([summary])], ignore_index=True)
    return df


def _sanity_check_rank_reuse_bootstrap(
    pair_matrix: np.ndarray, ranks: np.ndarray, labels: np.ndarray,
    n_check: int = 10, rng: np.random.Generator | None = None,
) -> Dict[str, float]:
    if rng is None:
        rng = np.random.default_rng(0)
    N, K = pair_matrix.shape
    devs_max = []
    devs_mean = []
    for _ in range(n_check):
        idx = rng.integers(0, N, N)
        y = labels[idx]
        pos_mask = y == 1
        n_pos = int(pos_mask.sum())
        n_neg = N - n_pos
        if n_pos == 0 or n_neg == 0:
            continue
        sum_pos_legacy = np.empty(K, dtype=np.float64)
        for col in range(K):
            ranks_col = sps.rankdata(pair_matrix[idx, col], method="average")
            sum_pos_legacy[col] = ranks_col[pos_mask].sum()
        auroc_legacy = (sum_pos_legacy - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
        sum_pos_fast = ranks[idx[pos_mask]].sum(axis=0, dtype=np.float64)
        auroc_fast = (sum_pos_fast - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
        d = np.abs(auroc_legacy - auroc_fast)
        devs_max.append(float(d.max()))
        devs_mean.append(float(d.mean()))
    return {
        "n_check":        n_check,
        "max_abs_dAUROC": float(np.mean(devs_max)) if devs_max else np.nan,
        "max_max_dAUROC": float(np.max(devs_max)) if devs_max else np.nan,
        "mean_abs_dAUROC": float(np.mean(devs_mean)) if devs_mean else np.nan,
    }


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
PROT_COLORS = {
    "CBS":  "#1f77b4",
    "GAL4": "#2ca02c",
    "PABP": "#D7263D",
    "PTEN": "#9467bd",
    "TEM1": "#ff7f0e",
}


def fig_A_global_null(global_summary, null_dist, out: Path) -> None:
    keys = PRIMARY_STATS
    fig, axes = plt.subplots(1, len(keys), figsize=(4.0 * len(keys), 3.6))
    for ax, key in zip(axes, keys):
        t_null = null_dist[key]
        t_obs = global_summary[key]["T_obs"]
        z = global_summary[key]["z_score"]
        p = global_summary[key]["p_empirical"]
        ax.hist(t_null, bins=60, color="#9AA1A8", edgecolor="white", linewidth=0.3)
        ax.axvline(t_obs, color="#D7263D", linewidth=1.8, label="observed")
        ax.set_xlabel(STAT_LABELS[key])
        ax.set_ylabel("count")
        ax.set_title(f"{key}  |  z={z:.2f}, p={p:.1e}")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(loc="upper right")
    fig.suptitle(f"Global permutation null distributions "
                 f"(within-protein shuffle, N_eff = {global_summary['T1']['n_perm_effective']:,})",
                 y=1.02)
    fig.tight_layout()
    fig.savefig(out.with_suffix(".png"))
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def fig_B_per_protein_null(per_prot, null_per_prot, out: Path) -> None:
    keys = PRIMARY_STATS
    fig, axes = plt.subplots(len(PROTEINS), len(keys),
                             figsize=(4.0 * len(keys), 2.4 * len(PROTEINS)))
    for r, prot in enumerate(PROTEINS):
        for c, key in enumerate(keys):
            ax = axes[r, c]
            t_null = null_per_prot[prot][key]
            t_obs = per_prot[prot][key]["T_obs"]
            z = per_prot[prot][key]["z_score"]
            p = per_prot[prot][key]["p_empirical"]
            ax.hist(t_null, bins=40, color="#CCCCCC",
                    edgecolor="white", linewidth=0.3)
            ax.axvline(t_obs, color=PROT_COLORS[prot], linewidth=1.7)
            if r == 0:
                ax.set_title(STAT_LABELS[key], fontsize=10)
            if c == 0:
                ax.set_ylabel(f"{prot}\ncount", color=PROT_COLORS[prot])
            ax.text(0.97, 0.92, f"z={z:.2f}\np={p:.1e}",
                    transform=ax.transAxes, ha="right", va="top",
                    fontsize=8,
                    bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
    fig.suptitle(f"Per-protein permutation null distributions "
                 f"(within-protein shuffle, N = {N_PERM_PROT:,} each)", y=1.005)
    fig.tight_layout()
    fig.savefig(out.with_suffix(".png"))
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def fig_C_observed_vs_null_bar(global_summary, per_prot, out: Path) -> None:
    keys = PRIMARY_STATS
    fig, ax = plt.subplots(figsize=(10, 4.4))
    groups = ["Global"] + PROTEINS
    n_g = len(groups)
    n_s = len(keys)
    bar_w = 0.8 / n_s
    xs = np.arange(n_g)
    palette = ["#1B4965", "#D7263D", "#F4A261", "#7B2D7E"]
    for k, key in enumerate(keys):
        zs = [global_summary[key]["z_score"]] + \
             [per_prot[p][key]["z_score"] for p in PROTEINS]
        offset = (k - (n_s - 1) / 2) * bar_w
        ax.bar(xs + offset, zs, width=bar_w, label=key,
               edgecolor="black", linewidth=0.4, color=palette[k])
    ax.set_xticks(xs)
    ax.set_xticklabels(groups)
    ax.axhline(0, color="black", linewidth=0.6)
    ax.axhline(1.96, color="grey", linewidth=0.6, linestyle=":")
    ax.axhline(-1.96, color="grey", linewidth=0.6, linestyle=":")
    z_bonf = sps.norm.ppf(1 - 0.05 / (2 * N_BONFERRONI_TESTS))
    ax.axhline(z_bonf, color="red", linewidth=0.6, linestyle="--",
               label=f"Bonferroni |z|≈{z_bonf:.2f}")
    ax.axhline(-z_bonf, color="red", linewidth=0.6, linestyle="--")
    ax.set_ylabel("z-score (T_obs − null mean) / null SD")
    ax.set_title("Observed test-statistic z-scores vs permutation null")
    ax.legend(loc="upper left", ncol=n_s + 1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out.with_suffix(".png"))
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


# ---------------------------------------------------------------------------
# FIG D  —  HOTFIX for v3 canvas explosion
# ---------------------------------------------------------------------------
def fig_D_qq_significance(global_summary, per_prot, out: Path) -> None:
    """
    QQ-style half-normal plot of observed |z| vs expected under H0.

    HOTFIX (v3.1):
      - Filter out NaN / Inf z-scores (can arise when null_std == 0).
      - Increase figsize from (6.6, 6.2) -> (10, 10) to accommodate 24 labels.
      - Remove fig.tight_layout() which interacts badly with savefig(bbox='tight')
        and many text artists, causing canvas height to explode to 73k px.
      - Add horizontal stagger to text labels and extra x-margin (lim * 1.2).
    """
    labels, zs, kinds = [], [], []
    for key in PRIMARY_STATS:
        labels.append(f"Global·{key}")
        zs.append(global_summary[key]["z_score"])
        kinds.append("Global")
        for prot in PROTEINS:
            labels.append(f"{prot}·{key}")
            zs.append(per_prot[prot][key]["z_score"])
            kinds.append(prot)

    zs = np.asarray(zs, dtype=np.float64)
    abs_zs = np.abs(zs)

    # ---- 1. Filter non-finite values -------------------------------------
    finite_mask = np.isfinite(abs_zs)
    if not finite_mask.all():
        n_drop = int((~finite_mask).sum())
        print(f"      [fig_D] Warning: dropping {n_drop} non-finite z-score(s)")
    abs_zs = abs_zs[finite_mask]
    labels = [labels[i] for i in np.where(finite_mask)[0]]
    kinds  = [kinds[i]  for i in np.where(finite_mask)[0]]

    # ---- 2. Sort by |z| --------------------------------------------------
    order = np.argsort(abs_zs)
    sorted_abs_zs = abs_zs[order]
    sorted_labels = [labels[i] for i in order]
    sorted_kinds  = [kinds[i]  for i in order]
    n_tests = len(abs_zs)
    p_emp = (np.arange(1, n_tests + 1) - 0.5) / n_tests
    expected = sps.norm.ppf((p_emp + 1) / 2.0)

    # ---- 3. Plot ---------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 10))   # enlarged from (6.6, 6.2)
    palette = {"Global": "#000000", **PROT_COLORS}

    for x, y, lab, kind in zip(expected, sorted_abs_zs, sorted_labels, sorted_kinds):
        ax.scatter(x, y, s=90, c=palette[kind], edgecolor="white",
                   linewidth=0.6, zorder=3)

    # Staggered text labels: alternate above/below to reduce overlap
    for i, (x, y, lab, kind) in enumerate(
        zip(expected, sorted_abs_zs, sorted_labels, sorted_kinds)
    ):
        x_text = x + 0.06
        # alternate vertical alignment; nudge y slightly for close points
        va = "bottom" if i % 2 == 0 else "top"
        y_text = y + (0.04 if i % 2 == 0 else -0.04)
        ax.text(
            x_text, y_text, lab,
            fontsize=8,
            va=va,
            ha="left",
            color=palette[kind],
            fontweight="bold" if kind == "Global" else "normal",
            zorder=3,
        )

    lim = max(expected.max(), sorted_abs_zs.max()) * 1.20   # extra room for labels
    ax.plot([0, lim], [0, lim], color="grey", linewidth=0.8, linestyle="--",
            label="H0 (y=x)", zorder=1)
    ax.axhline(1.96, color="black", linewidth=0.6, linestyle=":",
               label="|z|=1.96 (α=0.05)", zorder=2)
    z_bonf = sps.norm.ppf(1 - 0.05 / (2 * N_BONFERRONI_TESTS))
    ax.axhline(z_bonf, color="red", linewidth=0.6, linestyle=":",
               label=f"|z|={z_bonf:.2f} (Bonferroni)", zorder=2)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("Expected |z| under H₀ (half-normal quantile)")
    ax.set_ylabel("Observed |z|")
    ax.set_title(f"QQ-style significance plot ({n_tests} tests)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    handles = [Patch(facecolor=c, label=name) for name, c in palette.items()]
    ax.legend(handles=handles + ax.get_legend_handles_labels()[0],
              loc="lower right", fontsize=8, frameon=False)

    # ---- 4. Save (NO tight_layout here) ---------------------------------
    # rcParams already has savefig.bbox='tight', which is sufficient.
    # fig.tight_layout() was removed because it fights with text labels
    # and causes canvas size explosion (>2^16 px).
    fig.savefig(out.with_suffix(".png"))
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


# ---------------------------------------------------------------------------
# Helpers for CSV assembly
# ---------------------------------------------------------------------------
def bonferroni_adjust(p: float, n_tests: int) -> float:
    return min(1.0, float(p) * n_tests)


def bh_fdr(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=np.float64)
    out = np.full_like(p, np.nan)
    finite_mask = np.isfinite(p)
    if not finite_mask.any():
        return out
    p_finite = p[finite_mask]
    n = len(p_finite)
    order = np.argsort(p_finite)
    ranks = np.arange(1, n + 1, dtype=np.float64)
    sorted_q = p_finite[order] * n / ranks
    for i in range(n - 2, -1, -1):
        sorted_q[i] = min(sorted_q[i], sorted_q[i + 1])
    sorted_q = np.minimum(sorted_q, 1.0)
    q_finite = np.empty(n)
    q_finite[order] = sorted_q
    out[finite_mask] = q_finite
    return out


def fmt_ci(lo: float, hi: float, digits: int = 3) -> str:
    return f"[{lo:.{digits}g}, {hi:.{digits}g}]"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    print("=" * 70)
    print("p5_permutation.py  v3.1  —  permutation + bootstrap validation of F1-F5")
    print("=" * 70)

    rng = np.random.default_rng(RNG_SEED)
    boot_rng = np.random.default_rng(RNG_SEED + 1)

    # 1. Load
    print("\n[1/9] Loading SCI matrices and metadata ...")
    sci_full, sci_top50, meta = load_inputs()
    labels = meta["DMS_score_bin"].astype(int).to_numpy()
    proteins_arr = meta["protein"].astype(str).to_numpy()
    print(f"      n_mutations = {len(labels):,}")
    print(f"      n_pos = {int(labels.sum()):,}  n_neg = {int((1-labels).sum()):,}")

    protein_groups = []
    print("      per-protein balance:")
    for prot in PROTEINS:
        idx = np.where(proteins_arr == prot)[0]
        protein_groups.append(idx)
        pos = int(labels[idx].sum())
        print(f"        {prot:5s}: n={len(idx):6d}  pos={pos:5d}  neg={len(idx)-pos:5d}")

    # 2. Materialise pair matrix and ranks (float32)
    print("\n[2/9] Materialising (n, 528) pair matrix in float32 ...")
    pair_matrix = materialise_pair_matrix(sci_full)
    del sci_full
    print(f"      pair_matrix: shape={pair_matrix.shape}, "
          f"dtype={pair_matrix.dtype}, mem={pair_matrix.nbytes/1e6:.0f} MB")

    print("\n[3/9] Precomputing per-pair ranks ...")
    t0 = time.time()
    ranks = precompute_ranks(pair_matrix)
    print(f"      ranks: shape={ranks.shape}, dtype={ranks.dtype}, "
          f"mem={ranks.nbytes/1e6:.0f} MB ({time.time()-t0:.1f}s)")

    # 3. Observed values
    print("\n[4/9] Observed test statistics ...")
    obs_global = compute_all_stats(sci_top50, labels, ranks)
    for s in ALL_STATS:
        print(f"      {s:6s} obs = {obs_global[s]:+.4f}")

    two_sided_flags = {"T1": True, "T2": False, "T3": False, "T3_std": False,
                       "T_max": False}

    # 4. Global adaptive permutation
    print("\n[5/9] Global adaptive permutation ...")
    cfg = RunCfg(label="global", two_sided=two_sided_flags,
                 stages=PERM_STAGES, stop_p=PERM_STOP_P)
    global_null, eff_n_global = run_permutation_block(
        ranks, sci_top50, labels, protein_groups, obs_global, cfg, rng)

    global_summary = {s: summarise(global_null[s], obs_global[s],
                                    two_sided_flags[s], eff_n_global)
                      for s in ALL_STATS}

    # 5. Per-protein permutation (fixed N)
    print(f"\n[6/9] Per-protein permutation (N = {N_PERM_PROT:,}) ...")
    per_prot: Dict[str, Dict[str, Dict[str, float]]] = {p: {} for p in PROTEINS}
    null_per_prot: Dict[str, Dict[str, np.ndarray]] = {p: {} for p in PROTEINS}

    for prot, idx in zip(PROTEINS, protein_groups):
        print(f"\n   --- {prot} (n = {len(idx):,}) ---")
        sub_ranks = ranks[idx]
        sub_top50 = sci_top50[idx]
        sub_labels = labels[idx]
        obs_p = compute_all_stats(sub_top50, sub_labels, sub_ranks)
        cfg_p = RunCfg(label=prot, two_sided=two_sided_flags,
                       stages=[N_PERM_PROT], stop_p=0,
                       n_perm_fixed=N_PERM_PROT)
        null_p, eff_p = run_permutation_block(
            sub_ranks, sub_top50, sub_labels, None, obs_p, cfg_p, rng)
        per_prot[prot] = {s: summarise(null_p[s], obs_p[s],
                                       two_sided_flags[s], eff_p)
                          for s in ALL_STATS}
        null_per_prot[prot] = null_p

    # 6. Bonferroni + BH-FDR
    for prot in PROTEINS:
        for s in PRIMARY_STATS:
            per_prot[prot][s]["p_bonferroni"] = bonferroni_adjust(
                per_prot[prot][s]["p_empirical"], N_BONFERRONI_TESTS)

    pp_keys = [(prot, s) for prot in PROTEINS for s in PRIMARY_STATS]
    pp_ps = np.array([per_prot[prot][s]["p_empirical"] for (prot, s) in pp_keys])
    pp_qs = bh_fdr(pp_ps)
    for (prot, s), q in zip(pp_keys, pp_qs):
        per_prot[prot][s]["p_fdr"] = float(q)
    for prot in PROTEINS:
        per_prot[prot]["T_max"]["p_fdr"] = np.nan
        per_prot[prot]["T_max"]["p_bonferroni"] = np.nan

    g_ps = np.array([global_summary[s]["p_empirical"] for s in PRIMARY_STATS])
    g_qs = bh_fdr(g_ps)
    for s, q in zip(PRIMARY_STATS, g_qs):
        global_summary[s]["p_fdr"] = float(q)
    global_summary["T_max"]["p_fdr"] = np.nan

    # 7. Bootstrap
    print(f"\n[7/9] Bootstrap (T1: N={N_BOOT_T1}, T2/T3/T3_std/T_max: N={N_BOOT_T23}) ...")
    print("   T1 bootstrap ...")
    t1_boot = bootstrap_T1(sci_top50, labels, N_BOOT_T1, boot_rng)
    print("   T2/T3/T3_std/T_max bootstrap (rank-reuse fast path) ...")
    pair_boot = bootstrap_pair_stats(ranks, labels, N_BOOT_T23,
                                     boot_rng, label="boot")

    boot_rows = []
    for stat, vec in [("T1", t1_boot),
                      ("T2", pair_boot["T2"]),
                      ("T3", pair_boot["T3"]),
                      ("T3_std", pair_boot["T3_std"]),
                      ("T_max", pair_boot["T_max"])]:
        ci95 = np.nanpercentile(vec, [2.5, 97.5])
        ci99 = np.nanpercentile(vec, [0.5, 99.5])
        boot_rows.append({
            "Statistic": stat,
            "T_obs": obs_global[stat],
            "boot_mean": float(np.nanmean(vec)),
            "boot_std": float(np.nanstd(vec, ddof=1)),
            "bootstrap_ci_low_95": float(ci95[0]),
            "bootstrap_ci_high_95": float(ci95[1]),
            "bootstrap_ci_low_99": float(ci99[0]),
            "bootstrap_ci_high_99": float(ci99[1]),
            "n_boot": len(vec),
        })
    df_boot = pd.DataFrame(boot_rows)
    df_boot.to_csv(OUT_DIR / "bootstrap_summary.csv", index=False,
                   float_format="%.6g")

    # 8. Bootstrap validation: Exact vs Rank-Reuse
    print(f"\n[8/9] Bootstrap validation (Exact vs Rank-Reuse, N={N_BOOT_VALIDATION}) ...")
    df_valid = bootstrap_validation_exact_vs_fast(
        pair_matrix, ranks, labels, N_BOOT_VALIDATION, boot_rng
    )
    df_valid.to_csv(OUT_DIR / "bootstrap_validation.csv", index=False,
                    float_format="%.6g")
    print("      Validation summary:")
    for stat in ["T2", "T3", "T3_std"]:
        diffs = pd.to_numeric(df_valid[df_valid["boot_id"] != "SUMMARY"][f"{stat}_diff"],
                              errors="coerce")
        print(f"        {stat}: mean |diff| = {diffs.abs().mean():.4g}, "
              f"max |diff| = {diffs.abs().max():.4g}, "
              f"RMSE = {np.sqrt((diffs**2).mean()):.4g}")

    # 9. Write CSVs and figures
    print("\n[9/9] Writing CSVs and figures ...")

    # Global summary
    g_rows = []
    for s in ALL_STATS:
        rec = global_summary[s]
        g_rows.append({
            "Statistic": s,
            "Description": STAT_LABELS[s].replace("$", ""),
            "T_obs": rec["T_obs"],
            "null_mean": rec["null_mean"],
            "null_std": rec["null_std"],
            "effect_raw": rec["effect_raw"],
            "z_score": rec["z_score"],
            "cohens_d": rec["cohens_d"],
            "p_empirical": rec["p_empirical"],
            "p_fdr": rec.get("p_fdr", np.nan),
            "two_sided": rec["two_sided"],
            "null_ci95_low": rec["null_ci95_low"],
            "null_ci95_high": rec["null_ci95_high"],
            "null_ci99_low": rec["null_ci99_low"],
            "null_ci99_high": rec["null_ci99_high"],
            "n_perm_effective": rec["n_perm_effective"],
        })
    pd.DataFrame(g_rows).to_csv(OUT_DIR / "permutation_global_summary.csv",
                                 index=False, float_format="%.6g")

    # Per-protein summary
    p_rows = []
    for prot in PROTEINS:
        for s in ALL_STATS:
            rec = per_prot[prot][s]
            row = {
                "Protein": prot,
                "Statistic": s,
                "T_obs": rec["T_obs"],
                "null_mean": rec["null_mean"],
                "null_std": rec["null_std"],
                "effect_raw": rec["effect_raw"],
                "z_score": rec["z_score"],
                "cohens_d": rec["cohens_d"],
                "p_empirical": rec["p_empirical"],
                "p_bonferroni": rec.get("p_bonferroni", np.nan),
                "p_fdr": rec.get("p_fdr", np.nan),
                "two_sided": rec["two_sided"],
                "null_ci95_low": rec["null_ci95_low"],
                "null_ci95_high": rec["null_ci95_high"],
                "null_ci99_low": rec["null_ci99_low"],
                "null_ci99_high": rec["null_ci99_high"],
                "n_perm_effective": rec["n_perm_effective"],
            }
            p_rows.append(row)
    pd.DataFrame(p_rows).to_csv(OUT_DIR / "permutation_per_protein_summary.csv",
                                 index=False, float_format="%.6g")

    # Permutation effect sizes
    effect_rows = []
    for s in ALL_STATS:
        rec = global_summary[s]
        effect_rows.append({
            "Scope": "Global",
            "Statistic": s,
            "Cohen_d": rec["cohens_d"],
            "Z": rec["z_score"],
            "EmpiricalP": rec["p_empirical"],
        })
    for prot in PROTEINS:
        for s in ALL_STATS:
            rec = per_prot[prot][s]
            effect_rows.append({
                "Scope": prot,
                "Statistic": s,
                "Cohen_d": rec["cohens_d"],
                "Z": rec["z_score"],
                "EmpiricalP": rec["p_empirical"],
            })
    pd.DataFrame(effect_rows).to_csv(OUT_DIR / "permutation_effect_sizes.csv",
                                      index=False, float_format="%.6g")

    # Publication summary
    def _fmt_q(v):
        return "—" if not np.isfinite(v) else f"{v:.2e}"

    pub_rows = []
    for s in ALL_STATS:
        rec = global_summary[s]
        boot_row = next(b for b in boot_rows if b["Statistic"] == s)
        pub_rows.append({
            "Scope":         "Global",
            "Statistic":     s,
            "Observed":      f"{rec['T_obs']:.4g}",
            "NullMean":      f"{rec['null_mean']:.4g}",
            "NullSD":        f"{rec['null_std']:.4g}",
            "EffectRaw":     f"{rec['effect_raw']:.4g}",
            "Z":             f"{rec['z_score']:.2f}",
            "CohensD":       f"{rec['cohens_d']:.2f}",
            "EmpiricalP":    f"{rec['p_empirical']:.2e}",
            "Bonferroni":    "—",
            "FDR":           _fmt_q(rec.get("p_fdr", np.nan)),
            "CI95":          fmt_ci(rec['null_ci95_low'], rec['null_ci95_high']),
            "CI99":          fmt_ci(rec['null_ci99_low'], rec['null_ci99_high']),
            "BootCI95":      fmt_ci(boot_row['bootstrap_ci_low_95'],
                                    boot_row['bootstrap_ci_high_95']),
            "N_perm":        rec["n_perm_effective"],
        })
    for prot in PROTEINS:
        for s in ALL_STATS:
            rec = per_prot[prot][s]
            pub_rows.append({
                "Scope":         prot,
                "Statistic":     s,
                "Observed":      f"{rec['T_obs']:.4g}",
                "NullMean":      f"{rec['null_mean']:.4g}",
                "NullSD":        f"{rec['null_std']:.4g}",
                "EffectRaw":     f"{rec['effect_raw']:.4g}",
                "Z":             f"{rec['z_score']:.2f}",
                "CohensD":       f"{rec['cohens_d']:.2f}",
                "EmpiricalP":    f"{rec['p_empirical']:.2e}",
                "Bonferroni":    _fmt_q(rec.get("p_bonferroni", np.nan)),
                "FDR":           _fmt_q(rec.get("p_fdr", np.nan)),
                "CI95":          fmt_ci(rec['null_ci95_low'], rec['null_ci95_high']),
                "CI99":          fmt_ci(rec['null_ci99_low'], rec['null_ci99_high']),
                "BootCI95":      "—",
                "N_perm":        rec["n_perm_effective"],
            })
    pd.DataFrame(pub_rows).to_csv(OUT_DIR / "publication_summary.csv", index=False)

    # Persist raw nulls
    np.savez_compressed(
        OUT_DIR / "permutation_null_distributions.npz",
        **{f"global_{s}": global_null[s] for s in ALL_STATS},
        **{f"{prot}_{s}": null_per_prot[prot][s]
           for prot in PROTEINS for s in ALL_STATS},
        t1_boot=t1_boot,
        **{f"boot_{s}": pair_boot[s] for s in ["T2", "T3", "T3_std", "T_max"]},
    )

    # Figures
    print("   plotting figures ...")
    fig_A_global_null(global_summary, global_null,
                      OUT_DIR / "fig_A_global_null")
    fig_B_per_protein_null(per_prot, null_per_prot,
                           OUT_DIR / "fig_B_per_protein_null")
    fig_C_observed_vs_null_bar(global_summary, per_prot,
                               OUT_DIR / "fig_C_observed_vs_null_bar")
    fig_D_qq_significance(global_summary, per_prot,
                          OUT_DIR / "fig_D_qq_significance")

    # Console summary
    print("\n" + "=" * 70)
    print("Summary  (also written to publication_summary.csv)")
    print("=" * 70)
    print(pd.DataFrame(pub_rows).to_string(index=False))
    print(f"\nAll outputs written to: {OUT_DIR}")


if __name__ == "__main__":
    main()