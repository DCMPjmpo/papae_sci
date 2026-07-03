# Integration Report — Manuscript First Draft

## 1. Scanned Manuscript Files

| File | Path | Size (chars) | Sections |
|---|---|---:|---|
| final_draft.md | `paper/final_draft.md` | 185,874 | Title, Abstract, Introduction, Results (2.1–2.8), Discussion, Methods (4.1–4.12), Figure Legends |
| figure_map.md | `paper/figure_map.md` | 12,179 | Figure mapping for C1–C7 claims |
| supplementary.md | `paper/supplementary.md` | 6,478 | Tables S1–S2, Figures S1–S5 |
| manuscript_first_draft.md | `paper/manuscript_first_draft.md` | 190,342 | **Generated first draft** |

### Manuscript Structure Summary

```
Title
├── Abstract
├── Introduction (1)
├── Results (2.1–2.8)
│   ├── 2.1 SCI discriminates functional from non-functional mutations (C1)
│   ├── 2.2 Top-ranked layer pairs are dominated by Early–Late interactions (C2)
│   ├── 2.3 Early–Late pairs dominate the full 528-pair landscape (C3)
│   ├── 2.4 PABP forms an independent cluster on the layer-pair composition axis (C4)
│   ├── 2.5 Recurrent cross-protein layer pairs concentrate at layers 30–32 (C5)
│   ├── 2.6 Per-protein reproducibility (C6)
│   ├── 2.7 Statistical robustness (C7)
│   └── 2.8 SCI vs. ESM-2 zero-shot LLR baseline
├── Discussion (3)
├── Methods (4.1–4.12)
│   ├── 4.1 Dataset and preprocessing
│   ├── 4.2 ESM-2 feature extraction
│   ├── 4.3 SCI computation
│   ├── 4.4 Per-mutation SCI summaries
│   ├── 4.5 Discrimination metrics
│   ├── 4.6 Layer-pair scoring
│   ├── 4.7 Permutation test
│   ├── 4.8 Multiple-testing correction
│   ├── 4.9 Non-parametric bootstrap
│   ├── 4.10 Statistical reporting
│   ├── 4.11 Code and data availability
│   └── 4.12 ESM-2 zero-shot LLR baseline
├── Figure Legends (5)
│   ├── Figure 1 (C1)
│   ├── Figure 2 (C2 + C3)
│   ├── Figure 3 (C4)
│   ├── Figure 4 (C5)
│   ├── Figure 5 (C6 + C7)
│   └── Table 1 (§2.8)
└── Supplementary Information (6)
    ├── Table S1 — Dataset characteristics
    ├── Table S2 — Full statistical results
    └── Figures S1–S5
```

## 2. Scanned Figure Files

### Figure Inventory

| Figure | Panel | Source File | Legend Mapping | Current Status |
|---|---|---|---|---|
| Figure 1 | Panel A | `figure_assembly/calibration_crops/Fig1_PanelA.png` | C1: SCI distribution stratified by DMS label | ✅ Present |
| Figure 1 | Panel B | `figure_assembly/calibration_crops/Fig1_PanelB.png` | C1: Permutation null of T₁ | ✅ Present |
| Figure 1 | Composite | — | — | ❌ Missing (Fig1.png) |
| Figure 2 | Panel A | `figure_assembly/calibration_crops/Fig2_PanelA.png` | C2: Layer-pair category composition | ✅ Present |
| Figure 2 | Panel B | `figure_assembly/calibration_crops/Fig2_PanelB.png` | C2: T₂ permutation null | ✅ Present |
| Figure 2 | Panel C | `figure_assembly/calibration_crops/Fig2_PanelC.png` | C3: T₃ permutation null | ✅ Present |
| Figure 2 | Composite | — | — | ❌ Missing (Fig2.png) |
| Figure 3 | Panel A | `figure_assembly/calibration_crops/Fig3_PanelA.png` | C4: Composition heatmap | ✅ Present |
| Figure 3 | Panel B | `figure_assembly/calibration_crops/Fig3_PanelB.png` | C4: Hierarchical clustering | ✅ Present |
| Figure 3 | Panel C | `figure_assembly/calibration_crops/Fig3_PanelC.png` | C4: PCA | ✅ Present |
| Figure 3 | Composite | — | — | ❌ Missing (Fig3.png) |
| Figure 4 | Panel A | `figure_assembly/calibration_crops/Fig4_PanelA.png` | C5: Recurrent layer-end counts | ✅ Present |
| Figure 4 | Panel B | `figure_assembly/calibration_crops/Fig4_PanelB.png` | C5: Layer-frequency histogram | ✅ Present |
| Figure 4 | Composite | — | — | ❌ Missing (Fig4.png) |
| Figure 5 | Panel A | `figure_assembly/calibration_crops/Fig5_PanelA.png` | C6: Per-protein permutation grid | ✅ Present |
| Figure 5 | Panel B | `figure_assembly/calibration_crops/Fig5_PanelB.png` | C7: Z-score bar chart | ✅ Present |
| Figure 5 | Panel C | `figure_assembly/calibration_crops/Fig5_PanelC.png` | C7: QQ significance plot | ✅ Present |
| Figure 5 | Composite | — | — | ❌ Missing (Fig5.png) |

### Additional Figure Sources (Raw PNGs)

- `figures/raw/` — Contains individual analysis figures including:
  - `sci_distribution.png` (Fig1A source)
  - `permutation_null_T1.png` (Fig1B source)
  - `layer_pair_category_distribution.png` (Fig2A source)
  - `permutation_null_T2.png` (Fig2B source)
  - `permutation_null_T3.png` (Fig2C source)
  - `composition_heatmap.png` (Fig3A source)
  - `hierarchical_clustering.png` (Fig3B source)
  - `pca_plot.png` (Fig3C source)
  - `recurrent_layer_counts.png` (Fig4A source)
  - `layer_frequency_histogram.png` (Fig4B source)
  - `per_protein_permutation_grid.png` (Fig5A source)
  - `z_score_bar_chart.png` (Fig5B source)
  - `qq_significance_plot.png` (Fig5C source)

## 3. Figure Insertion Status

| Figure | Target Section | Insertion Marker | Panel Crop References | Status |
|---|---|---|---|---|
| Figure 1 | Results 2.1 | `[INSERT FIGURE 1 HERE]` | Fig1_PanelA.png, Fig1_PanelB.png | ✅ Marked |
| Figure 2 | Results 2.2–2.3 | `[INSERT FIGURE 2 HERE]` | Fig2_PanelA.png, Fig2_PanelB.png, Fig2_PanelC.png | ✅ Marked |
| Figure 3 | Results 2.4 | `[INSERT FIGURE 3 HERE]` | Fig3_PanelA.png, Fig3_PanelB.png, Fig3_PanelC.png | ✅ Marked |
| Figure 4 | Results 2.5 | `[INSERT FIGURE 4 HERE]` | Fig4_PanelA.png, Fig4_PanelB.png | ✅ Marked |
| Figure 5 | Results 2.6–2.7 | `[INSERT FIGURE 5 HERE]` | Fig5_PanelA.png, Fig5_PanelB.png, Fig5_PanelC.png | ✅ Marked |

**Note:** Composite figures (Fig1.png–Fig5.png) are not present. Individual panel crops are referenced in the manuscript with relative paths to `figure_assembly/calibration_crops/`.

## 4. Missing Content

### Critical Missing Items

| Item | Location | Description | Action |
|---|---|---|---|
| Fig1.png | `figure_assembly/` | Composite Figure 1 (Panels A+B) | Documented — waiting for assembly pipeline |
| Fig2.png | `figure_assembly/` | Composite Figure 2 (Panels A+B+C) | Documented — waiting for assembly pipeline |
| Fig3.png | `figure_assembly/` | Composite Figure 3 (Panels A+B+C) | Documented — waiting for assembly pipeline |
| Fig4.png | `figure_assembly/` | Composite Figure 4 (Panels A+B) | Documented — waiting for assembly pipeline |
| Fig5.png | `figure_assembly/` | Composite Figure 5 (Panels A+B+C) | Documented — waiting for assembly pipeline |

### Potential Gaps

| Item | Description | Status |
|---|---|---|
| References | No dedicated References section | ❌ Not present |
| Author list | No author affiliations or correspondence | ❌ Not present |
| Keywords | No manuscript keywords | ❌ Not present |

## 5. Format Issues

### Identified Inconsistencies

| Issue | Location | Details | Severity |
|---|---|---|---|
| PABP distance claim | `figure_map.md` line 115 | States PABP distance "exceeds" within-other-four distance (0.389 < 0.397) | Medium |
| ΔAUROC value | `figure_map.md` Table 1 caption | Lists ΔAUROC = 0.008; `final_draft.md` uses 0.007 | Low |
| Caption alignment | `figure_map.md` | Some figure captions contain incomplete sentences (e.g., Fig3 caption cut off) | Low |

### Formatting Recommendations

| Category | Issue | Recommendation |
|---|---|---|
| Section numbering | Results uses 2.x, Discussion uses 3, Methods uses 4.x | Standard academic format — ✅ OK |
| Equation numbering | Subscript notation (T₁, T₂, T₃) | Consistent throughout — ✅ OK |
| Figure references | Text references match Figure Legend numbers | Verified — ✅ OK |
| Supplementary tables | Tables S1–S2 present and formatted | ✅ OK |

## 6. Submission Readiness Check

### Checklist

| Check | Status | Notes |
|---|---|---|
| Title | ✅ Present | "Structural Context Inconsistency (SCI) captures mutation-effect signal in the hidden states of ESM-2" |
| Abstract | ✅ Present | Full abstract with background, methods, results, conclusions |
| Introduction | ✅ Present | Context-setting, gap identification, paper structure |
| Results | ✅ Present | 8 subsections (2.1–2.8) covering all 7 claims (C1–C7) |
| Discussion | ✅ Present | Interpretation, limitations, future directions |
| Methods | ✅ Present | 12 subsections covering all computational details |
| Figure Legends | ✅ Present | Complete legends for all 5 figures + Table 1 |
| Supplementary | ✅ Present | Tables S1–S2, Figures S1–S5 |
| **References** | ❌ Missing | Requires addition before submission |
| **Author list** | ❌ Missing | Requires addition before submission |
| **Composite Figures** | ❌ Missing | Requires assembly from panel crops |

### Overall Assessment

**Manuscript Framework:** ✅ Complete  
**Scientific Content:** ✅ Complete (all claims C1–C7 covered)  
**Statistical Reporting:** ✅ Complete (permutation, bootstrap, multiple testing)  
**Figure Infrastructure:** ⚠️ Partial (panels present, composites missing)  
**Submission Ready:** ❌ No (missing References, Authors, Composite Figures)

### Next Steps for Submission

1. Run figure assembly pipeline to generate Fig1.png–Fig5.png
2. Add References section with proper citation formatting
3. Add Author list with affiliations and correspondence information
4. Verify all figure captions are complete and accurate
5. Resolve identified text inconsistencies (PABP distance claim, ΔAUROC value)

---

*Generated: 2025-07-01*  
*Source files: final_draft.md, figure_map.md, supplementary.md*  
*Mode: READ-ONLY — no modifications made to source content*