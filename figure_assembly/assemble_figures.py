#!/usr/bin/env python3
"""
Figure Assembly — ESM-2 Layer-Pair Manuscript
===============================================

Generates publication-quality composite figures (Fig1–Fig5) for
Nature Communications submission from existing panel PNGs.

Usage:
    cd d:/文件/工作室/website
    python figure_assembly/assemble_figures.py

Output:
    figures/Fig{1,2,3,4,5}.{png,tiff,pdf}   at 300 DPI

Requirements:
    pip install matplotlib pillow numpy

Layout specification: see figure_assembly/LAYOUT_SPECIFICATION.md
"""

import os, sys
from pathlib import Path
from PIL import Image
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox

# =========================================================================
# Paths
# =========================================================================
ROOT = Path("d:/文件/工作室/website")
OUT_DIR = ROOT / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Source panel directories
PROC = ROOT / "data/processed"
P5_PERM = ROOT / "data/p0_output/p5_permutation"
P_CLUST = ROOT / "data/p0_output/protein_clustering"
P_PROC_OUT = ROOT / "data/processed/p0_output"

# =========================================================================
# Global matplotlibrc — Nature Communications style
# =========================================================================
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "standard",          # no tight bbox — we control the canvas
    "savefig.pad_inches": 0,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
})

# Panel-label styling
LABEL_STYLE = dict(fontsize=8, fontweight="bold", fontfamily="sans-serif")


# =========================================================================
# Helper functions
# =========================================================================

def crop_panel(image: Image.Image, left: float, top: float,
               right: float, bottom: float) -> Image.Image:
    """
    Crop a sub-panel from a multi-panel figure using fractional coordinates
    (0–1 relative to image dimensions).
    """
    w, h = image.size
    return image.crop((left * w, top * h, right * w, bottom * h))


def load_image(path: Path, target_dpi: int = 300) -> Image.Image:
    """Load a PNG, convert to RGB, return PIL Image."""
    img = Image.open(path).convert("RGB")
    return img


def img_to_array(img: Image.Image) -> np.ndarray:
    """PIL Image → numpy array for matplotlib imshow."""
    return np.asarray(img)


def place_panel(fig, ax, image: Image.Image):
    """Display a PIL image on the given Axes, filling it exactly."""
    ax.imshow(np.asarray(image), aspect="equal", interpolation="bilinear")
    ax.axis("off")


def label_panel(ax, label: str, x: float = 0.03, y: float = 0.95):
    """
    Add a panel label ('A', 'B', …) in the top-left corner of *ax*,
    in axes-fraction coordinates (0–1).
    """
    ax.text(x, y, label, transform=ax.transAxes, **LABEL_STYLE,
            va="top", ha="left")


def save_figure(fig: plt.Figure, name: str):
    """Save figure as PNG, TIFF (LZW), and PDF."""
    png_path = OUT_DIR / f"{name}.png"
    tiff_path = OUT_DIR / f"{name}.tiff"
    pdf_path = OUT_DIR / f"{name}.pdf"

    fig.savefig(png_path, dpi=300)
    print(f"  [OK] {png_path}")

    fig.savefig(tiff_path, dpi=300, format="tiff",
                pil_kwargs={"compression": "tiff_lzw"})
    print(f"  [OK] {tiff_path}")

    fig.savefig(pdf_path, dpi=300, format="pdf")
    print(f"  [OK] {pdf_path}")


# =========================================================================
# Figure 1 — SCI carries mutation-effect signal (single column)
# =========================================================================
def build_figure_1():
    print("\n=== Building Figure 1 ===")
    fig_w = 3.346       # 85 mm single column
    fig_h = 5.315       # ≈ 135 mm
    height_ratios = [0.60, 0.40]
    gap = 0.04          # inches between panels

    fig = plt.figure(figsize=(fig_w, fig_h))
    total_gap = gap * (len(height_ratios) - 1)
    avail_h = fig_h - total_gap
    panel_heights = [avail_h * r / sum(height_ratios) for r in height_ratios]

    # --- Panel A: SCI distribution (from P1.5, crop top-left) ---
    src_a = load_image(PROC / "P1.5_sci_signal_validation.png", 300)
    # 2×2 gridspec, top-left quadrant; fraction estimated from the
    # code: figsize=(16,12), gridspec(2,2,hspace=0.32,wspace=0.32)
    # Crop generously: (left, top, right, bottom)
    panel_a_img = crop_panel(src_a, 0.03, 0.03, 0.47, 0.49)
    ax_a = fig.add_axes([0, (panel_heights[1] + gap) / fig_h, 1.0, panel_heights[0] / fig_h])
    place_panel(fig, ax_a, panel_a_img)
    label_panel(ax_a, "A", x=0.02, y=0.96)

    # --- Panel B: T1 null (from fig_A_global_null, T1 leftmost 25 %) ---
    src_b = load_image(P5_PERM / "fig_A_global_null.png", 300)
    # 4 equal subplots in a 1×4 row; T1 is the leftmost ~25 %
    panel_b_img = crop_panel(src_b, 0.00, 0.00, 0.25, 1.00)
    ax_b = fig.add_axes([0, 0, 1.0, panel_heights[1] / fig_h])
    place_panel(fig, ax_b, panel_b_img)
    label_panel(ax_b, "B", x=0.02, y=0.96)

    save_figure(fig, "Fig1")
    plt.close(fig)


# =========================================================================
# Figure 2 — Top-ranked layer pairs and the 528-pair landscape
# =========================================================================
def build_figure_2():
    print("\n=== Building Figure 2 ===")
    fig_w = 7.008       # 178 mm double column
    fig_h = 6.0         # ≈ 152 mm
    gap = 0.05

    fig = plt.figure(figsize=(fig_w, fig_h))

    # Layout: top panel A, bottom row B + C
    top_h = 0.38 * fig_h
    bot_h = fig_h - top_h - gap

    # --- Panel A: stacked bar (wide, full-width) ---
    src_a = load_image(P_PROC_OUT / "layer_pair_category_distribution.png", 300)
    ax_a = fig.add_axes([0, (bot_h + gap) / fig_h, 1.0, top_h / fig_h])
    place_panel(fig, ax_a, src_a)
    label_panel(ax_a, "A")

    # --- Panel B: T2 null (second 25 % of fig_A_global_null) ---
    src_bc = load_image(P5_PERM / "fig_A_global_null.png", 300)
    panel_b_img = crop_panel(src_bc, 0.25, 0.00, 0.50, 1.00)
    ax_b = fig.add_axes([0, 0, 0.49, bot_h / fig_h])
    place_panel(fig, ax_b, panel_b_img)
    label_panel(ax_b, "B")

    # --- Panel C: T3 null (third 25 % of fig_A_global_null) ---
    panel_c_img = crop_panel(src_bc, 0.50, 0.00, 0.75, 1.00)
    ax_c = fig.add_axes([0.51, 0, 0.49, bot_h / fig_h])
    place_panel(fig, ax_c, panel_c_img)
    label_panel(ax_c, "C")

    save_figure(fig, "Fig2")
    plt.close(fig)


# =========================================================================
# Figure 3 — PABP forms an independent cluster (C4)
# =========================================================================
def build_figure_3():
    print("\n=== Building Figure 3 ===")
    fig_w = 7.008       # double column
    fig_h = 5.5         # ≈ 140 mm
    gap = 0.05

    fig = plt.figure(figsize=(fig_w, fig_h))

    # Layout: A left ~57 % full height, B + C right ~43 % stacked
    left_w = 0.57 * fig_w
    right_w = fig_w - left_w - gap
    top_h_right = 0.48 * fig_h
    bot_h_right = fig_h - top_h_right - gap

    # --- Panel A: heatmap ---
    src_a = load_image(P_CLUST / "heatmap_layer_pair_categories.png", 300)
    ax_a = fig.add_axes([0, 0, left_w / fig_w, 1.0])
    place_panel(fig, ax_a, src_a)
    label_panel(ax_a, "A")

    # --- Panel B: dendrogram ---
    src_b = load_image(P_CLUST / "dendrogram_protein_clustering.png", 300)
    ax_b = fig.add_axes([(left_w + gap) / fig_w, (bot_h_right + gap) / fig_h, right_w / fig_w, top_h_right / fig_h])
    place_panel(fig, ax_b, src_b)
    label_panel(ax_b, "B")

    # --- Panel C: PCA ---
    src_c = load_image(P_CLUST / "pca_protein_layer_pair.png", 300)
    ax_c = fig.add_axes([(left_w + gap) / fig_w, 0, right_w / fig_w, bot_h_right / fig_h])
    place_panel(fig, ax_c, src_c)
    label_panel(ax_c, "C")

    save_figure(fig, "Fig3")
    plt.close(fig)


# =========================================================================
# Figure 4 — Cross-protein recurrent layer pairs (C5)
# =========================================================================
def build_figure_4():
    print("\n=== Building Figure 4 ===")
    fig_w = 7.008       # double column
    fig_h = 5.0         # ≈ 127 mm
    gap = 0.05

    fig = plt.figure(figsize=(fig_w, fig_h))
    half_h = (fig_h - gap) / 2

    # --- Panel A: recurrence barplot ---
    src_a = load_image(P_CLUST / "layer_recurrence_barplot.png", 300)
    ax_a = fig.add_axes([0, (half_h + gap) / fig_h, 1.0, half_h / fig_h])
    place_panel(fig, ax_a, src_a)
    label_panel(ax_a, "A")

    # --- Panel B: frequency histogram ---
    src_b = load_image(P_CLUST / "layer_frequency_histogram.png", 300)
    ax_b = fig.add_axes([0, 0, 1.0, half_h / fig_h])
    place_panel(fig, ax_b, src_b)
    label_panel(ax_b, "B")

    save_figure(fig, "Fig4")
    plt.close(fig)


# =========================================================================
# Figure 5 — Statistical robustness (C6 + C7)
# =========================================================================
def build_figure_5():
    print("\n=== Building Figure 5 ===")
    fig_w = 7.008       # double column
    fig_h = 8.27        # ≈ 210 mm
    gap = 0.05

    fig = plt.figure(figsize=(fig_w, fig_h))

    # Layout: A full-width top ~58 %, B + C side-by-side bottom ~42 %
    top_h = 0.58 * fig_h
    bot_h = fig_h - top_h - gap
    half_w = (fig_w - gap) / 2

    # --- Panel A: per-protein null grid (5×4) ---
    src_a = load_image(P5_PERM / "fig_B_per_protein_null.png", 300)
    ax_a = fig.add_axes([0, (bot_h + gap) / fig_h, 1.0, top_h / fig_h])
    place_panel(fig, ax_a, src_a)
    label_panel(ax_a, "A")

    # --- Panel B: z-score bar chart ---
    src_b = load_image(P5_PERM / "fig_C_observed_vs_null_bar.png", 300)
    ax_b = fig.add_axes([0, 0, half_w / fig_w, bot_h / fig_h])
    place_panel(fig, ax_b, src_b)
    label_panel(ax_b, "B")

    # --- Panel C: QQ significance plot ---
    src_c = load_image(P5_PERM / "fig_D_qq_significance.png", 300)
    ax_c = fig.add_axes([(half_w + gap) / fig_w, 0, half_w / fig_w, bot_h / fig_h])
    place_panel(fig, ax_c, src_c)
    label_panel(ax_c, "C")

    save_figure(fig, "Fig5")
    plt.close(fig)


# =========================================================================
# Main
# =========================================================================
def main():
    print("=" * 60)
    print("ESM-2 Layer-Pair Manuscript — Figure Assembly")
    print(f"Output directory: {OUT_DIR}")
    print("=" * 60)

    build_figure_1()
    build_figure_2()
    build_figure_3()
    build_figure_4()
    build_figure_5()

    print("\n" + "=" * 60)
    print("All figures assembled successfully.")
    print(f"Output files in: {OUT_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
