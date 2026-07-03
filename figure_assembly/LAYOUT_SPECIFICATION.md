# Publication Figure Layout Specification
## ESM-2 Layer-Pair Manuscript — Nature Communications Format

---

## General Specifications

| Parameter | Value |
|-----------|-------|
| Journal width (single column) | 85 mm = 3.346 in |
| Journal width (double column) | 178 mm = 7.008 in |
| Maximum figure height | 230 mm = 9.055 in |
| Output DPI | 300 (halftone/ composite) |
| Source DPI | Panels at 600 DPI, downsampled to 300 DPI |
| Font family | Arial (Helvetica fallback) |
| Panel label font | 8 pt bold Arial |
| Axis label font | 6.5–7 pt Arial |
| Tick label font | 6 pt Arial |
| Annotation font | 6 pt Arial |
| Line weight | 0.5–0.75 pt |
| Panel label offset | Top-left corner of each panel, dx = −0.05 in, dy = +0.05 in from panel edge |
| Colour space | RGB |
| File formats | TIFF (LZW compression) + PNG (for preview) |
| Naming | `Fig1.tiff`, `Fig2.tiff`, …, `Fig5.tiff` |

---

## Figure 1 — SCI carries mutation-effect signal (C1)

| Property | Value |
|----------|-------|
| Width | Single column: 85 mm (3.346 in) |
| Height | ≈ 135 mm (5.315 in) |
| Aspect | Portrait |
| Panels | 2 (A over B, vertical stack) |

```
┌──────────────────┐
│  A               │  ← Overlaid histograms of Top50-mean SCI
│  (60 % height)   │    (cropped from P1.5_sci_signal_validation.png)
│                  │
├──────────────────┤
│  B               │  ← T1 permutation null histogram
│  (40 % height)   │    (cropped from fig_A_global_null.png, T1 panel)
│                  │
└──────────────────┘
    3.346 in (85 mm)
```

### Panel A — SCI distribution
- **Source:** `data/processed/P1.5_sci_signal_validation.png` (3961×3599 px, 300 DPI)
- **Crop region:** Top-left panel of the 2×2 gridspec
- **Display width:** 3.346 in (1004 px at 300 DPI)
- **Display height:** ≈ 2.9 in (870 px)

### Panel B — T1 permutation null
- **Source:** `data/p0_output/p5_permutation/fig_A_global_null.png` (9540×2211 px, 600 DPI)
- **Crop region:** T1 subplot, leftmost ~25 % of image
- **Display width:** 3.346 in (1004 px at 300 DPI)
- **Display height:** ≈ 1.9 in (570 px)

---

## Figure 2 — Top-ranked layer pairs and the 528-pair landscape (C2 + C3)

| Property | Value |
|----------|-------|
| Width | Double column: 178 mm (7.008 in) |
| Height | ≈ 152 mm (6.0 in) |
| Aspect | Landscape |
| Panels | 3 (A full-width top, B + C side-by-side below) |

```
┌──────────────────────────────────────┐
│  A                                   │  ← Stacked bar chart, category
│  (40 % height)                       │    composition of Top50 layer pairs
│                                      │    (layer_pair_category_distribution.png)
├──────────────────┬───────────────────┤
│  B               │  C                │  B: T2 (AUROC-Top50 EL count) null
│  (60 % height)   │  (60 % height)    │  C: T3 (EL category dominance) null
│                  │                   │  (both cropped from fig_A_global_null.png)
└──────────────────┴───────────────────┘
    ───── 50 % ───── ──── 50 % ────
    7.008 in (178 mm)
```

### Panel A — Category distribution
- **Source:** `data/processed/p0_output/layer_pair_category_distribution.png` (1785×1035 px, 150 DPI)
- **Display width:** 7.008 in (2102 px at 300 DPI)
- **Display height:** ≈ 2.2 in (660 px)
- **Note:** Source is 150 DPI; upscale to 300 DPI for output

### Panel B — T2 null (AUROC-Top50 EL count)
- **Source:** `data/p0_output/p5_permutation/fig_A_global_null.png` (9540×2211 px, 600 DPI)
- **Crop:** T2 subplot (second quartile)
- **Display width:** 3.35 in (1005 px at 300 DPI)
- **Display height:** ≈ 3.2 in (960 px)

### Panel C — T3 null (EL dominance)
- **Source:** `data/p0_output/p5_permutation/fig_A_global_null.png`
- **Crop:** T3 subplot (third quartile)
- **Display width:** 3.35 in (1005 px at 300 DPI)
- **Display height:** ≈ 3.2 in (960 px)

---

## Figure 3 — PABP forms an independent cluster (C4)

| Property | Value |
|----------|-------|
| Width | Double column: 178 mm (7.008 in) |
| Height | ≈ 140 mm (5.5 in) |
| Aspect | Landscape |
| Panels | 3 (A left column full height, B + C stacked right column) |

```
┌───────────────────┬────────────────────┐
│                   │  B                 │  ← Hierarchical clustering
│  A                │  (dendrogram)      │    dendrogram
│  (heatmap)        │                    │    (dendrogram_protein_clustering.png)
│  ~57 % width      ├────────────────────┤
│  full height      │  C                 │  ← PCA scatter
│                   │  (PCA)             │    (pca_protein_layer_pair.png)
│                   │                    │
└───────────────────┴────────────────────┘
    ── 57 % ──────── ─── 43 % ─────────
    7.008 in (178 mm)
```

### Panel A — Heatmap
- **Source:** `data/p0_output/protein_clustering/heatmap_layer_pair_categories.png` (3820×2768 px, 600 DPI)
- **Display width:** ≈ 4.0 in (1200 px at 300 DPI)
- **Display height:** ≈ 5.0 in (1500 px)

### Panel B — Dendrogram
- **Source:** `data/p0_output/protein_clustering/dendrogram_protein_clustering.png` (3410×2480 px, 600 DPI)
- **Display width:** ≈ 3.0 in (900 px at 300 DPI)
- **Display height:** ≈ 2.4 in (720 px)

### Panel C — PCA
- **Source:** `data/p0_output/protein_clustering/pca_protein_layer_pair.png` (3456×2852 px, 600 DPI)
- **Display width:** ≈ 3.0 in (900 px at 300 DPI)
- **Display height:** ≈ 2.4 in (720 px)

---

## Figure 4 — Recurrent layer pairs at layers 30–32 (C5)

| Property | Value |
|----------|-------|
| Width | Double column: 178 mm (7.008 in) |
| Height | ≈ 127 mm (5.0 in) |
| Aspect | Landscape |
| Panels | 2 (A over B, vertical stack) |

```
┌──────────────────────────────────────┐
│  A                                   │  ← Recurrence barplot (Top20 set)
│  (50 % height)                       │    (layer_recurrence_barplot.png)
│                                      │
├──────────────────────────────────────┤
│  B                                   │  ← Layer-frequency histogram
│  (50 % height)                       │    (Top50 union)
│                                      │    (layer_frequency_histogram.png)
└──────────────────────────────────────┘
    7.008 in (178 mm)
```

### Panel A — Recurrence barplot
- **Source:** `data/p0_output/protein_clustering/layer_recurrence_barplot.png` (4539×2851 px, 600 DPI)
- **Display width:** 7.008 in (2102 px at 300 DPI)
- **Display height:** ≈ 2.3 in (690 px)

### Panel B — Frequency histogram
- **Source:** `data/p0_output/protein_clustering/layer_frequency_histogram.png` (4696×2739 px, 600 DPI)
- **Display width:** 7.008 in (2102 px at 300 DPI)
- **Display height:** ≈ 2.3 in (690 px)

---

## Figure 5 — Statistical robustness (C6 + C7)

| Property | Value |
|----------|-------|
| Width | Double column: 178 mm (7.008 in) |
| Height | ≈ 210 mm (8.27 in) |
| Aspect | Tall landscape |
| Panels | 3 (A full-width top ~58 %, B + C side-by-side bottom ~42 %) |

```
┌──────────────────────────────────────┐
│                                      │
│  A                                   │  ← 5×4 per-protein null grid
│  (58 % height)                       │    (fig_B_per_protein_null.png)
│                                      │
├──────────────────┬───────────────────┤
│  B               │  C                │  B: z-score bar chart
│  (42 % height)   │  (42 % height)    │    (fig_C_observed_vs_null_bar.png)
│                  │                   │  C: QQ significance plot
│                  │                   │    (fig_D_qq_significance.png)
└──────────────────┴───────────────────┘
    ──── 50 % ──── ───── 50 % ──────
    7.008 in (178 mm)
```

### Panel A — Per-protein null grid
- **Source:** `data/p0_output/p5_permutation/fig_B_per_protein_null.png` (9542×7262 px, 600 DPI)
- **Display width:** 7.008 in (2102 px at 300 DPI)
- **Display height:** ≈ 4.6 in (1380 px)
- **Note:** Source is large; scale down to ~30 % linear for 300 DPI output

### Panel B — z-score bar chart
- **Source:** `data/p0_output/p5_permutation/fig_C_observed_vs_null_bar.png` (5940×2576 px, 600 DPI)
- **Display width:** 3.35 in (1005 px at 300 DPI)
- **Display height:** ≈ 3.1 in (930 px)

### Panel C — QQ significance plot
- **Source:** `data/p0_output/p5_permutation/fig_D_qq_significance.png` (5041×5109 px, 600 DPI)
- **Display width:** 3.35 in (1005 px at 300 DPI)
- **Display height:** ≈ 3.1 in (930 px)

---

## Panel Label Placement

All panel labels follow the same convention:

```
┌────────────────────────┐
│ A                      │  ← Label: 8 pt bold, horizontal baseline
│   (panel content)      │     Position: (x0 − 2 pt, y1 + 2 pt)
│                        │     where (x0, y0) is panel bottom-left,
│                        │     (x1, y1) is panel top-right in figure coords
└────────────────────────┘
```

- **Font:** Arial, 8 pt, bold, black
- **Position:** Top-left corner of each panel's border, offset slightly inside
- **No box/background:** Label sits directly over the panel content
- **Letter style:** Single uppercase letter (A, B, C)
- **No period** after letter

---

## Output Specifications

| Format | Options | Use |
|--------|---------|-----|
| TIFF | LZW compression, RGB, 300 DPI | Journal submission |
| PNG | sRGB, 300 DPI | Preview / review |
| PDF | Vector labels + embedded raster panels | Editor revision |

---

## File Naming

```
figures/
├── Fig1.tiff          # Publication-ready TIFF
├── Fig1.png           # Preview PNG
├── Fig1.pdf           # Vector-editable PDF
├── Fig2.tiff
├── Fig2.png
├── Fig2.pdf
├── Fig3.tiff
├── Fig3.png
├── Fig3.pdf
├── Fig4.tiff
├── Fig4.png
├── Fig4.pdf
├── Fig5.tiff
├── Fig5.png
├── Fig5.pdf
└── assemble_figures.py  # Master assembly script
```
