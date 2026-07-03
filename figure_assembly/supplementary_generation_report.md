# Supplementary Figure Generation Report

**Date:** 2026-07-02
**Scope:** Resolve confirmed manuscript artefacts from audit (9 action items)

---

## 1. Files Created

| # | File | Size | Source data | Purpose |
|---|------|------|-------------|---------|
| 1 | `data/p0_5plots/sci_per_protein.png` | 439 KB | `all_proteins_sci_site_scores_top50.npy` + `all_proteins_site_metadata.csv` | Supp Fig S1 — per-protein SCI histograms (2×3 grid, 5 proteins + pooled) |
| 2 | `data/p0_5plots/sci_distribution.png` | 222 KB | Same as above | Supp Fig S1 — pooled SCI distribution (single panel, cross-reference) |
| 3 | `data/p0_5plots/sci_heatmap.png` | 371 KB | `all_proteins_sci_site_matrices.dat` (memmap) | Supp Fig S3 — representative 33×33 SCI matrix heatmaps (functional + non-functional) |
| 4 | `figures/supplementary/Supp_Fig_S1.png` | 399 KB | Copy of `data/p0_5plots/sci_per_protein.png` | Composite supplementary figure S1 |
| 5 | `figures/supplementary/Supp_Fig_S2.png` | 4,228 KB | `layer_pair_frequency_map.png` + `layer_pair_mining.png` | Composite supplementary figure S2 (A: AUROC heatmap, B: ranking) |
| 6 | `figures/supplementary/Supp_Fig_S3.png` | 3,438 KB | `sci_heatmap.png` + `best_layer_pair_distribution.png` | Composite supplementary figure S3 (A: 33×33 matrices, B: best-pair distribution) |
| 7 | `figures/supplementary/Supp_Fig_S4.png` | 994 KB | Copy of `fig_B_per_protein_null.png` | Composite supplementary figure S4 (per-protein permutation grid) |
| 8 | `figures/supplementary/Supp_Fig_S5.png` | 432 KB | Copy of `fig_D_qq_significance.png` | Composite supplementary figure S5 (QQ significance plot) |

---

## 2. Files Reused (existing on disk, unchanged)

| # | File | Used for | Source |
|---|------|----------|--------|
| 1 | `data/processed/all_proteins_sci_site_scores_top50.npy` | Generating S1 panels | Existing from p6 pipeline |
| 2 | `data/processed/all_proteins_site_metadata.csv` | Generating S1 panels | Existing from p6 pipeline |
| 3 | `data/processed/all_proteins_sci_site_matrices.dat` | Generating S3 panels | Existing from p6 pipeline |
| 4 | `data/p0_5layer/layer_pair_frequency_map.png` | S2 panel A | Existing from p0_5layer |
| 5 | `data/p0_5layer/layer_pair_mining.png` | S2 panel B | Existing from p0_5layer |
| 6 | `data/p0_5layer/best_layer_pair_distribution.png` | S3 panel B | Existing from p0_5layer |
| 7 | `data/p0_output/p5_permutation/fig_B_per_protein_null.png` | S4 | Existing from p5_permutation pipeline |
| 8 | `data/p0_output/p5_permutation/fig_D_qq_significance.png` | S5 | Existing from p5_permutation pipeline |

---

## 3. Files Modified

| # | File | Change | Lines affected |
|---|------|--------|----------------|
| 1 | `paper/final_draft.md` | Replaced placeholder reference block (`> Placeholder. The reference list will be compiled…`) with complete 30-item numbered reference list from `paper/final_draft_submission.md` | Lines 1412–1444 |

---

## 4. Generation Scripts Created

| # | File | Purpose |
|---|------|---------|
| 1 | `figure_assembly/generate_supp_panels.py` | Generates missing data panels: S1 (per-protein SCI histograms) and S3 (SCI 33×33 heatmap). Can be re-run if data changes. |
| 2 | `figure_assembly/assemble_supp_figures.py` | Assembles all 5 supplementary figure composites from existing + newly generated panels. Can be re-run independently. |

---

## 5. Unresolved Issues

| # | Issue | Status | Action needed |
|---|-------|--------|---------------|
| 1 | `data/p0_5plots/sci_per_protein.png` — Figure legend says "overlaid histograms ... Vertical dashed lines mark the two group means" | ✅ Implemented | None |
| 2 | Main composite figures (Fig1–Fig5) still not assembled as singletons | ❌ Unresolved (out of scope) | Run `figure_assembly/assemble_figures.py` |
| 3 | Supp Fig S2 — the two source panels have different visual styles (different DPI, colour schemes) | ⚠️ Cosmetic | Optional: re-generate both panels at consistent 600 DPI with matching colour palettes |
| 4 | `paper/final_draft.md` — reference block replaced from `final_draft_submission.md`; DOIs and full venue names still need to be added | ⚠️ Bibliographic | Verify each entry against `paper/literature_matrix/literature_matrix.csv` |
| 5 | Table 1 — exists as CSV but not typeset as a formatted journal table | ❌ Unresolved (out of scope) | Typeset at submission time |
| 6 | `paper/supplementary.md` line 58: `(see CSV)` placeholders remain in per-protein T_max z column | ⚠️ Minor | Replace with actual values from `permutation_per_protein_summary.csv` |

---

## 6. Output Summary

```
figures/supplementary/
├── Supp_Fig_S1.png   (399 KB)  # Per-protein SCI distributions [NEW]
├── Supp_Fig_S2.png   (4.2 MB)  # 528-pair AUROC landscape     [NEW]
├── Supp_Fig_S3.png   (3.4 MB)  # SCI matrix + best-pair dist  [NEW]
├── Supp_Fig_S4.png   (994 KB)  # Per-protein permutation grid [NEW]
└── Supp_Fig_S5.png   (432 KB)  # QQ significance summary      [NEW]

data/p0_5plots/
├── sci_per_protein.png   (439 KB)  # Newly generated panel
├── sci_distribution.png  (222 KB)  # Newly generated panel
└── sci_heatmap.png       (371 KB)  # Newly generated panel
```
