#!/usr/bin/env python3
"""
Calibration Utility — Figure Panel Crop Regions
=================================================

Run this script to inspect the crop regions used in assemble_figures.py.
It displays the extracted panel images so you can visually verify that
each crop captures the correct sub-panel.

Usage:
    cd d:/文件/工作室/website
    python figure_assembly/calibrate_crops.py

Adjust the crop_fractions in the PANEL_CROPS dict below until
each cropped image shows exactly the intended panel content.
Then update the corresponding values in assemble_figures.py.
"""

import sys
from pathlib import Path
from PIL import Image

ROOT = Path("d:/文件/工作室/website")
P5_PERM = ROOT / "data/p0_output/p5_permutation"
PROC = ROOT / "data/processed"
P_PROC_OUT = ROOT / "data/processed/p0_output"

OUT = ROOT / "figure_assembly/calibration_crops"
OUT.mkdir(parents=True, exist_ok=True)

# =========================================================================
# CROP REGIONS — adjust these fractions until the output looks correct
# =========================================================================
# Each crop is (left, top, right, bottom) as fractions of the source image.
PANEL_CROPS = {
    # Figure 1
    "Fig1_PanelA": {
        "source": PROC / "P1.5_sci_signal_validation.png",
        "crop":  (0.03, 0.03, 0.47, 0.49),
        "desc":  "Fig1a — SCI distribution (top-left of 2×2 grid)",
    },
    "Fig1_PanelB": {
        "source": P5_PERM / "fig_A_global_null.png",
        "crop":  (0.00, 0.00, 0.25, 1.00),
        "desc":  "Fig1b — T1 null (leftmost of 4 horizontal panels)",
    },

    # Figure 2 — panels B and C are cropped from the same fig_A_global_null
    "Fig2_PanelA": {
        "source": P_PROC_OUT / "layer_pair_category_distribution.png",
        "crop":  (0.0, 0.0, 1.0, 1.0),   # full image
        "desc":  "Fig2a — Category distribution (full image)",
    },
    "Fig2_PanelB": {
        "source": P5_PERM / "fig_A_global_null.png",
        "crop":  (0.25, 0.00, 0.50, 1.00),
        "desc":  "Fig2b — T2 null (second of 4 horizontal panels)",
    },
    "Fig2_PanelC": {
        "source": P5_PERM / "fig_A_global_null.png",
        "crop":  (0.50, 0.00, 0.75, 1.00),
        "desc":  "Fig2c — T3 null (third of 4 horizontal panels)",
    },

    # Figure 3 — standalone panels, full images
    "Fig3_PanelA": {
        "source": ROOT / "data/p0_output/protein_clustering/heatmap_layer_pair_categories.png",
        "crop":  (0.0, 0.0, 1.0, 1.0),
        "desc":  "Fig3a — Heatmap (full image)",
    },
    "Fig3_PanelB": {
        "source": ROOT / "data/p0_output/protein_clustering/dendrogram_protein_clustering.png",
        "crop":  (0.0, 0.0, 1.0, 1.0),
        "desc":  "Fig3b — Dendrogram (full image)",
    },
    "Fig3_PanelC": {
        "source": ROOT / "data/p0_output/protein_clustering/pca_protein_layer_pair.png",
        "crop":  (0.0, 0.0, 1.0, 1.0),
        "desc":  "Fig3c — PCA (full image)",
    },

    # Figure 4 — standalone panels
    "Fig4_PanelA": {
        "source": ROOT / "data/p0_output/protein_clustering/layer_recurrence_barplot.png",
        "crop":  (0.0, 0.0, 1.0, 1.0),
        "desc":  "Fig4a — Recurrence barplot (full image)",
    },
    "Fig4_PanelB": {
        "source": ROOT / "data/p0_output/protein_clustering/layer_frequency_histogram.png",
        "crop":  (0.0, 0.0, 1.0, 1.0),
        "desc":  "Fig4b — Frequency histogram (full image)",
    },

    # Figure 5 — standalone panels
    "Fig5_PanelA": {
        "source": P5_PERM / "fig_B_per_protein_null.png",
        "crop":  (0.0, 0.0, 1.0, 1.0),
        "desc":  "Fig5a — Per-protein null grid (full image)",
    },
    "Fig5_PanelB": {
        "source": P5_PERM / "fig_C_observed_vs_null_bar.png",
        "crop":  (0.0, 0.0, 1.0, 1.0),
        "desc":  "Fig5b — z-score bar chart (full image)",
    },
    "Fig5_PanelC": {
        "source": P5_PERM / "fig_D_qq_significance.png",
        "crop":  (0.0, 0.0, 1.0, 1.0),
        "desc":  "Fig5c — QQ significance plot (full image)",
    },
}


def main():
    print("=" * 60)
    print("Crop Region Calibration")
    print("=" * 60)

    for key, info in PANEL_CROPS.items():
        src_path = info["source"]
        left, top, right, bottom = info["crop"]
        desc = info["desc"]

        print(f"\n{key}: {desc}")
        print(f"  Source: {src_path}")
        print(f"  Crop:   left={left:.3f}, top={top:.3f}, right={right:.3f}, bottom={bottom:.3f}")

        if not src_path.exists():
            print(f"  ⚠ MISSING — skipping")
            continue

        img = Image.open(src_path).convert("RGB")
        w, h = img.size
        print(f"  Size:   {w}×{h} px @ {img.info.get('dpi', '?')} DPI")

        # Crop
        box = (left * w, top * h, right * w, bottom * h)
        cropped = img.crop(box)
        out_path = OUT / f"{key}.png"
        cropped.save(out_path)
        print(f"  ✓ Saved: {out_path} ({cropped.size[0]}×{cropped.size[1]} px)")

    print(f"\n{'=' * 60}")
    print(f"All calibration crops saved to: {OUT}")
    print(f"Inspect each .png to verify the crop region is correct.")
    print(f"If a crop is wrong, adjust the fractions in PANEL_CROPS")
    print(f"and re-run this script, then update assemble_figures.py.")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
