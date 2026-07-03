#!/usr/bin/env python3
"""
Assemble Supplementary Figures
===============================

Creates the 5 supplementary figure composite images:

- Supp_Fig_S1.png — Per-protein SCI distributions (copy of sci_per_protein.png)
- Supp_Fig_S2.png — 528-pair AUROC landscape + Composite ranking
- Supp_Fig_S3.png — SCI matrix heatmap + best-pair distribution
- Supp_Fig_S4.png — Per-protein permutation grid (copy of fig_B_per_protein_null.png)
- Supp_Fig_S5.png — QQ significance summary (copy of fig_D_qq_significance.png)

Usage:
    cd d:/文件/工作室/website
    python figure_assembly/assemble_supp_figures.py

Output:
    figures/supplementary/Supp_Fig_S{1,2,3,4,5}.png
"""

import sys
from pathlib import Path
from PIL import Image
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# =========================================================================
# Paths
# =========================================================================
ROOT = Path("d:/文件/工作室/website")
OUT_DIR = ROOT / "figures" / "supplementary"
OUT_DIR.mkdir(parents=True, exist_ok=True)

P5_PERM = ROOT / "data/p0_output/p5_permutation"
P_CLUST = ROOT / "data/p0_output/protein_clustering"
P_PLOTS = ROOT / "data/p0_5plots"
P_LAYER = ROOT / "data/p0_5layer"

# =========================================================================
# Style
# =========================================================================
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 8,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
})

LABEL_STYLE = dict(fontsize=10, fontweight="bold", fontfamily="sans-serif")


def place_panel(fig, ax, image: Image.Image):
    ax.imshow(np.asarray(image.convert("RGB")), aspect="equal", interpolation="bilinear")
    ax.axis("off")


def label_panel(ax, label: str, x: float = 0.02, y: float = 0.96):
    ax.text(x, y, label, transform=ax.transAxes, **LABEL_STYLE,
            va="top", ha="left")


def save_figure(fig: plt.Figure, name: str):
    out_path = OUT_DIR / f"{name}.png"
    fig.savefig(out_path, dpi=300)
    print(f"  [OK] {out_path.name}")
    plt.close(fig)


# =========================================================================
# S1 — Per-protein SCI distributions (single existing composite)
# =========================================================================
def build_supp_s1():
    """Copy sci_per_protein.png (already a complete 2x3 grid)."""
    print("\n=== Supp Fig S1 ===")
    src = P_PLOTS / "sci_per_protein.png"
    out = OUT_DIR / "Supp_Fig_S1.png"
    Image.open(src).convert("RGB").save(out)
    print(f"  [OK] {out.name}")


# =========================================================================
# S2 — 528-pair AUROC landscape + Composite ranking
# =========================================================================
def build_supp_s2():
    """
    Two panels stacked vertically:
      Top:    layer_pair_frequency_map.png  (AUROC heatmap)
      Bottom: layer_pair_mining.png        (ranking visualisation)
    """
    print("\n=== Supp Fig S2 ===")
    fig_w, fig_h = 7.2, 8.0
    gap = 0.04
    half_h = (fig_h - gap) / 2

    fig = plt.figure(figsize=(fig_w, fig_h))

    # Panel A: frequency map
    img_a = Image.open(P_LAYER / "layer_pair_frequency_map.png").convert("RGB")
    ax_a = fig.add_axes([0, half_h + gap, fig_w, half_h])
    place_panel(fig, ax_a, img_a)
    label_panel(ax_a, "A")

    # Panel B: layer_pair_mining
    img_b = Image.open(P_LAYER / "layer_pair_mining.png").convert("RGB")
    ax_b = fig.add_axes([0, 0, fig_w, half_h])
    place_panel(fig, ax_b, img_b)
    label_panel(ax_b, "B")

    save_figure(fig, "Supp_Fig_S2")


# =========================================================================
# S3 — SCI matrix heatmap + best-pair distribution
# =========================================================================
def build_supp_s3():
    """
    Two panels side-by-side:
      Left:  sci_heatmap.png  (newly generated 33x33 matrix)
      Right: best_layer_pair_distribution.png (existing)
    """
    print("\n=== Supp Fig S3 ===")
    fig_w, fig_h = 10.0, 5.2
    gap = 0.06
    half_w = (fig_w - gap) / 2

    fig = plt.figure(figsize=(fig_w, fig_h))

    # Panel A: SCI heatmap
    img_a = Image.open(P_PLOTS / "sci_heatmap.png").convert("RGB")
    ax_a = fig.add_axes([0, 0, half_w, fig_h])
    place_panel(fig, ax_a, img_a)
    label_panel(ax_a, "A")

    # Panel B: best layer-pair distribution
    img_b = Image.open(P_LAYER / "best_layer_pair_distribution.png").convert("RGB")
    ax_b = fig.add_axes([half_w + gap, 0, half_w, fig_h])
    place_panel(fig, ax_b, img_b)
    label_panel(ax_b, "B")

    save_figure(fig, "Supp_Fig_S3")


# =========================================================================
# S4 — Per-protein permutation grid (copy of fig_B_per_protein_null.png)
# =========================================================================
def build_supp_s4():
    """Copy fig_B_per_protein_null.png (identical to Fig 5a, larger format)."""
    print("\n=== Supp Fig S4 ===")
    src = P5_PERM / "fig_B_per_protein_null.png"
    out = OUT_DIR / "Supp_Fig_S4.png"
    Image.open(src).convert("RGB").save(out)
    print(f"  [OK] {out.name}")


# =========================================================================
# S5 — QQ significance summary (copy of fig_D_qq_significance.png)
# =========================================================================
def build_supp_s5():
    """Copy fig_D_qq_significance.png (identical to Fig 5c)."""
    print("\n=== Supp Fig S5 ===")
    src = P5_PERM / "fig_D_qq_significance.png"
    out = OUT_DIR / "Supp_Fig_S5.png"
    Image.open(src).convert("RGB").save(out)
    print(f"  [OK] {out.name}")


# =========================================================================
# Main
# =========================================================================
def main():
    print("=" * 60)
    print("Assembling Supplementary Figures")
    print("=" * 60)

    build_supp_s1()
    build_supp_s2()
    build_supp_s3()
    build_supp_s4()
    build_supp_s5()

    print("\n" + "=" * 60)
    print(f"All supplementary figures in: {OUT_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
