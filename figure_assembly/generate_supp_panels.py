#!/usr/bin/env python3
"""
Generate Missing Supplementary Figure Panels
==============================================

Creates the panel images that do not yet exist on disk:

- data/p0_5plots/sci_per_protein.png   (Supp Fig S1 — per-protein histograms)
- data/p0_5plots/sci_distribution.png   (Supp Fig S1 — pooled histogram)
- data/p0_5plots/sci_heatmap.png        (Supp Fig S3 — representative 33×33 matrix)

Usage:
    cd d:/文件/工作室/website
    python figure_assembly/generate_supp_panels.py

Requirements:
    pip install matplotlib numpy pandas scipy
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# =========================================================================
# Paths
# =========================================================================
ROOT = Path("d:/文件/工作室/website")
DATA_DIR = ROOT / "data/processed"
OUT_DIR = ROOT / "data/p0_5plots"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SCORES_PATH = DATA_DIR / "all_proteins_sci_site_scores_top50.npy"
META_PATH = DATA_DIR / "all_proteins_site_metadata.csv"
MATRICES_PATH = DATA_DIR / "all_proteins_sci_site_matrices.dat"

# =========================================================================
# Global style — Nature Communications
# =========================================================================
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
})

PROTEINS = ["CBS", "GAL4", "PABP", "PTEN", "TEM1"]
PROTEIN_COLORS = {
    "CBS": "#1f77b4", "GAL4": "#2ca02c", "PABP": "#D7263D",
    "PTEN": "#9467bd", "TEM1": "#ff7f0e",
}
FUNC_COLOR = "#2166AC"      # blue — functional
NONFUNC_COLOR = "#B2182B"   # red — non-functional


# =========================================================================
# Helper
# =========================================================================
def load_data():
    """Load SCI top50 scores and metadata, aligned by row index."""
    scores = np.load(SCORES_PATH).astype(np.float64)
    meta = pd.read_csv(META_PATH)
    assert len(scores) == len(meta), f"Length mismatch: {len(scores)} vs {len(meta)}"
    return scores, meta


# =========================================================================
# Figure S1 — Per-protein SCI distribution histograms
# =========================================================================
def generate_sci_per_protein(scores, meta):
    """
    Supplementary Figure S1.

    Multi-panel figure (2×3 grid): 5 per-protein panels + 1 pooled panel.
    Each panel shows overlaid histograms of Top50-mean SCI score stratified
    by DMS_score_bin (blue = functional, red = non-functional).
    """
    n_cols, n_rows = 3, 2
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(10, 6.5))

    # Flatten axes for iteration
    axes_flat = axes.flatten()

    # Determine global bin edges for consistency
    all_func = scores[meta["DMS_score_bin"] == 1]
    all_nonfunc = scores[meta["DMS_score_bin"] == 0]
    combined = np.concatenate([all_func, all_nonfunc])
    bins = np.linspace(combined.min(), combined.max(), 50)

    # Plot each protein + pooled
    for idx, prot in enumerate(PROTEINS + ["Pooled"]):
        ax = axes_flat[idx]

        if prot == "Pooled":
            mask_func = meta["DMS_score_bin"] == 1
            mask_nonfunc = meta["DMS_score_bin"] == 0
        else:
            prot_mask = meta["protein"] == prot
            mask_func = prot_mask & (meta["DMS_score_bin"] == 1)
            mask_nonfunc = prot_mask & (meta["DMS_score_bin"] == 0)

        s_func = scores[mask_func.values] if hasattr(mask_func, 'values') else scores[mask_func]
        s_nonfunc = scores[mask_nonfunc.values] if hasattr(mask_nonfunc, 'values') else scores[mask_nonfunc]

        # Histograms — semi-transparent
        ax.hist(s_func, bins=bins, alpha=0.6, color=FUNC_COLOR,
                edgecolor="white", linewidth=0.2, label=f"Functional (n={len(s_func):,})")
        ax.hist(s_nonfunc, bins=bins, alpha=0.6, color=NONFUNC_COLOR,
                edgecolor="white", linewidth=0.2, label=f"Non-functional (n={len(s_nonfunc):,})")

        # Vertical dashed lines for means
        mean_func = np.mean(s_func)
        mean_nonfunc = np.mean(s_nonfunc)
        ax.axvline(mean_func, color=FUNC_COLOR, linestyle="--", linewidth=0.8)
        ax.axvline(mean_nonfunc, color=NONFUNC_COLOR, linestyle="--", linewidth=0.8)

        # Styling
        title = prot if prot != "Pooled" else "Pooled (All proteins)"
        ax.set_title(title, fontweight="bold", pad=6)
        ax.set_xlabel("Top50-mean SCI score")
        ax.set_ylabel("Count")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if idx == 0:
            ax.legend(loc="upper right", frameon=True, fontsize=6,
                      facecolor="white", edgecolor="#CCCCCC", framealpha=0.9)

    # Remove any leftover subplot
    for idx in range(len(PROTEINS) + 1, n_rows * n_cols):
        axes_flat[idx].set_visible(False)

    fig.suptitle("Supplementary Figure S1 - Per-protein SCI distributions",
                 fontsize=11, fontweight="bold", y=1.01)

    fig.tight_layout()
    out_path = OUT_DIR / "sci_per_protein.png"
    fig.savefig(out_path, dpi=600)
    print(f"  [OK] {out_path.name}")
    plt.close(fig)


# =========================================================================
# Figure S1 (pooled) — Pooled SCI distribution
# =========================================================================
def generate_sci_distribution(scores, meta):
    """
    Pooled SCI distribution — a single-panel version for cross-reference.
    Same structure as each panel in sci_per_protein but larger and with
    more detailed annotations.
    """
    fig, ax = plt.subplots(figsize=(7, 5))

    mask_func = meta["DMS_score_bin"] == 1
    mask_nonfunc = meta["DMS_score_bin"] == 0

    s_func = scores[mask_func.values]
    s_nonfunc = scores[mask_nonfunc.values]

    bins = np.linspace(scores.min(), scores.max(), 60)

    ax.hist(s_func, bins=bins, alpha=0.6, color=FUNC_COLOR,
            edgecolor="white", linewidth=0.3, label=f"Functional (n={len(s_func):,})")
    ax.hist(s_nonfunc, bins=bins, alpha=0.6, color=NONFUNC_COLOR,
            edgecolor="white", linewidth=0.3, label=f"Non-functional (n={len(s_nonfunc):,})")

    mean_func = np.mean(s_func)
    mean_nonfunc = np.mean(s_nonfunc)
    ax.axvline(mean_func, color=FUNC_COLOR, linestyle="--", linewidth=1.0,
               label=f"Mean functional = {mean_func:.4f}")
    ax.axvline(mean_nonfunc, color=NONFUNC_COLOR, linestyle="--", linewidth=1.0,
               label=f"Mean non-functional = {mean_nonfunc:.4f}")

    ax.set_xlabel("Top50-mean SCI score")
    ax.set_ylabel("Count")
    ax.set_title("Pooled SCI distribution (n = 95,142 mutations)", fontweight="bold")
    ax.legend(loc="upper right", frameon=True, facecolor="white",
              edgecolor="#CCCCCC", framealpha=0.9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Annotate T1
    t1 = mean_func - mean_nonfunc
    ax.text(0.97, 0.95, f"T1 = {t1:.5f}\n(z = 57.0, p_FDR = 10-3)",
            transform=ax.transAxes, ha="right", va="top", fontsize=8,
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"))

    fig.tight_layout()
    out_path = OUT_DIR / "sci_distribution.png"
    fig.savefig(out_path, dpi=600)
    print(f"  [OK] {out_path.name}")
    plt.close(fig)


# =========================================================================
# Figure S3 — Representative SCI 33×33 matrix heatmap
# =========================================================================
def generate_sci_heatmap(scores, meta):
    """
    Supplementary Figure S3 (left panel).

    Representative per-mutation 33×33 SCI matrix heatmaps.
    Shows one functional mutation and one non-functional mutation
    side by side, with a colorbar.
    """
    # Select representative examples: median SCI functional and non-functional
    # to avoid outliers
    mask_func = meta["DMS_score_bin"] == 1
    mask_nonfunc = meta["DMS_score_bin"] == 0

    func_indices = np.where(mask_func.values)[0]
    nonfunc_indices = np.where(mask_nonfunc.values)[0]

    # Pick mutations near the 50th percentile of SCI scores for each group
    s_func = scores[func_indices]
    s_nonfunc = scores[nonfunc_indices]

    # Compute distance from median
    med_func = np.median(s_func)
    med_nonfunc = np.median(s_nonfunc)

    idx_func = func_indices[np.argmin(np.abs(s_func - med_func))]
    idx_nonfunc = nonfunc_indices[np.argmin(np.abs(s_nonfunc - med_nonfunc))]

    # Load matrices via memmap
    shape = (95142, 33, 33)
    dtype = np.dtype('float32')
    mm = np.memmap(str(MATRICES_PATH), dtype=dtype, mode='r', shape=shape)

    mat_func = mm[idx_func].copy()
    mat_nonfunc = mm[idx_nonfunc].copy()
    del mm

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

    vmin = 0.0
    vmax = 1.0

    for ax, mat, label, idx in zip(
        axes,
        [mat_func, mat_nonfunc],
        ["Functional mutation", "Non-functional mutation"],
        [idx_func, idx_nonfunc]
    ):
        im = ax.imshow(mat, aspect="equal", cmap="RdBu_r",
                       vmin=vmin, vmax=vmax, interpolation="nearest")

        ax.set_xticks(np.arange(0, 33, 4))
        ax.set_xticklabels(np.arange(1, 34, 4))
        ax.set_yticks(np.arange(0, 33, 4))
        ax.set_yticklabels(np.arange(1, 34, 4))
        ax.set_xlabel("Layer j")
        ax.set_ylabel("Layer i")
        ax.set_title(f"{label}\n(row index {idx})", fontsize=9, fontweight="bold")

        # Add band separators
        for x in (7.5, 23.5):
            ax.axvline(x, color="black", linewidth=0.4, linestyle=":", alpha=0.5)
            ax.axhline(x, color="black", linewidth=0.4, linestyle=":", alpha=0.5)

    # Shared colorbar
    fig.subplots_adjust(right=0.88, wspace=0.25)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.015, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("SCI value (Pearson r)", rotation=270, labelpad=14)
    cbar.outline.set_linewidth(0.6)

    fig.suptitle("Supplementary Figure S3 - Representative per-mutation SCI matrices",
                 fontsize=11, fontweight="bold", y=0.98)

    fig.savefig(OUT_DIR / "sci_heatmap.png", dpi=600)
    print(f"  [OK] sci_heatmap.png")
    plt.close(fig)


# =========================================================================
# Main
# =========================================================================
def main():
    print("=" * 60)
    print("Generating missing supplementary figure panels")
    print("=" * 60)

    print("\nLoading SCI scores and metadata...")
    scores, meta = load_data()
    print(f"  {len(scores)} mutations loaded")

    print("\n--- Supp Fig S1: Per-protein SCI distributions ---")
    generate_sci_per_protein(scores, meta)
    generate_sci_distribution(scores, meta)

    print("\n--- Supp Fig S3: SCI 33x33 matrix heatmap ---")
    generate_sci_heatmap(scores, meta)

    print("\n" + "=" * 60)
    print(f"All panels generated in: {OUT_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
