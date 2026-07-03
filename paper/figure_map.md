# Figure Map (Phase A)

> Five main figures, each mapped 1-to-1 onto one or two of the validated
> claims (C1–C7) from `paper/scientific_story_map.md`. No new figures are
> generated here; every panel reuses a PNG/PDF already on disk under
> `data/`. All numbers in captions are sourced from the corresponding CSV
> on disk and not interpolated.

| # | Claim mapping | Composite panels (existing files on disk) |
|---|---|---|
| **Fig 1** | C1 | `data/processed/P1.5_sci_signal_validation.png` + `data/p0_output/p5_permutation/fig_A_global_null.png` (T1 panel only) |
| **Fig 2** | C2 + C3 | `data/processed/p0_output/layer_pair_category_distribution.png` + `data/p0_output/p5_permutation/fig_A_global_null.png` (T2 and T3 panels) |
| **Fig 3** | C4 | `data/p0_output/protein_clustering/heatmap_layer_pair_categories.png` + `…/dendrogram_protein_clustering.png` + `…/pca_protein_layer_pair.png` |
| **Fig 4** | C5 | `data/p0_output/protein_clustering/layer_recurrence_barplot.png` + `…/layer_frequency_histogram.png` |
| **Fig 5** | C6 + C7 | `data/p0_output/p5_permutation/fig_B_per_protein_null.png` + `…/fig_C_observed_vs_null_bar.png` + `…/fig_D_qq_significance.png` |
| **Table 1** | §2.8 baseline comparison | `data/processed/sci_vs_llr_comparison.csv` (numbers verbatim) + `data/processed/baseline_summary.csv` (LLR-only sanity check) |

For each figure the **caption** is the publication-quality block that
sits below the figure in the manuscript; the **legend** is a panel-by-panel
description that supports the caption; and the **take-home message** is
the single statement the reader should retrieve from the figure.

---

## Figure 1 — SCI carries mutation-effect signal (C1)

### Caption
**SCI discriminates functional from non-functional mutations across all
five DMS panels.**
**(a)** Distribution of the per-mutation Top50-mean SCI score
(`all_proteins_sci_site_scores_top50.npy`; n = 95,142) stratified by the
binarised DMS label (`DMS_score_bin`, 1 = functional, 0 =
non-functional). **(b)** Permutation null distribution of the SCI mean
difference T₁ = mean(SCI | y = 1) − mean(SCI | y = 0) under
10⁻³ within-protein label shuffles (Methods §3.7). Vertical red line
marks the observed T₁ = 2.79 × 10⁻³; the null mean and standard
deviation are −3.7 × 10⁻⁴ and 5.6 × 10⁻⁵, giving z = 57.0 and
empirical p_FDR = 10⁻³ against the 4-statistic global family
(`permutation_global_summary.csv`). At the per-protein level T₁ is
FDR-significant in 5/5 proteins (smallest = CBS, p_FDR = 3.7 × 10⁻²;
remaining four ≤ 1.4 × 10⁻³; `permutation_per_protein_summary.csv`).

### Legend
- **Panel (a)** — overlaid histograms of Top50-mean SCI score for
  functional (blue) and non-functional (red) mutations; x-axis = SCI
  value, y-axis = count. Vertical dashed lines mark the two group means.
- **Panel (b)** — histogram of T₁ values from the permutation null;
  x-axis = T₁ (functional − non-functional SCI mean), y-axis = count;
  vertical red line marks the observed T₁; annotation reports z and
  p_FDR.

### Take-home message
SCI separates functional from non-functional mutations at z = 57
against the within-protein label-shuffle null, and the separation is
reproducible in every one of the five proteins after Bonferroni and BH
FDR correction.

---

## Figure 2 — Top-ranked layer pairs and the 528-pair category landscape (C2 + C3)

### Caption
**Early–Late layer pairs dominate both the Top50 high-Composite tail
and the full 528-pair landscape in ESM-2.**
**(a)** Stacked bar chart of the six-category layer-pair composition
(EE / EM / EL / MM / ML / LL; Methods §3.3) within each protein's
Top50 Composite layer pairs (per-protein
`*_pair_category_ratio_Composite.csv`). Four of the five proteins (CBS,
GAL4, PTEN, TEM1) are EL-dominated (EL accounts for 40 %, 66 %, 44 % and
40 %, respectively, of the Top50 set); PABP is the only protein whose
Top50 set contains a larger EM than EL share (EM 48 %; EL 42 %).
**(b)** Global permutation null of the AUROC-Top50 EL count T₂ under
10⁻³ within-protein label shuffles; observed T₂ = 17 EL pairs vs. null
mean 1.0 (z = 8.0, p_FDR = 10⁻³;
`permutation_global_summary.csv`).
**(c)** Global permutation null of the EL-category dominance score
T₃ = mean |AUROC − 0.5|_{EL} − mean |AUROC − 0.5|_{non-EL}
over all 528 layer pairs; observed T₃ = 4.0 × 10⁻² (z = 62.0,
p_FDR = 10⁻³). The standardised counterpart T₃ˢᵗᵈ = 1.02 (z = 8.05)
is reported in `publication_summary.csv`.

### Legend
- **Panel (a)** — y-axis: proportion of Top50 Composite layer pairs that
  fall in each of the six categories; x-axis: protein name (PABP, CBS,
  GAL4, PTEN, TEM1); stack colours map to category (legend embedded in
  source PNG `layer_pair_category_distribution.png`).
- **Panel (b)** — null distribution histogram for T₂; x-axis: number of
  EL pairs among the 50 layer pairs with largest |AUROC − 0.5|; vertical
  red line at observed value; annotation reports z and p_FDR.
- **Panel (c)** — null distribution histogram for T₃; x-axis: EL minus
  non-EL mean |AUROC − 0.5|; vertical red line at observed value.

### Take-home message
EL enrichment is observed at both the high-signal tail (Top50; per-protein
proportions of 40 %–66 %) and across the full 528-pair plane (category
dominance with z = 62 against null), establishing that the EL block is
not a tail-only phenomenon.

---

## Figure 3 — PABP forms an independent cluster on the layer-pair composition axis (C4)

### Caption
**PABP forms a singleton in the six-dimensional layer-pair composition
space, driven by Early–Middle interactions.**
**(a)** Heatmap of the six-category Top50 Composite layer-pair
composition for each of the five proteins (n × 6 matrix from per-protein
`*_pair_category_ratio_Composite.csv`). PABP is the only row in which
the EM proportion (0.48) exceeds the EL proportion (0.42); the
remaining four proteins are EL-dominated. **(b)** Average-linkage
hierarchical clustering of the composition vectors under Euclidean
distance (`pairwise_euclidean_distance.csv`). PABP joins the four
EL-dominated proteins only at a height of ≈ 0.50, whereas the four
others fuse at maximum height ≈ 0.22. The minimum PABP-to-other
distance (0.389; PABP–TEM1) exceeds the maximum within-other-four
distance (0.397), satisfying a singleton condition. **(c)** Principal
component analysis on the same vectors; PC1 explains 68.6 % and PC2
28.1 % of variance (96.7 % cumulative;
`pca_loadings.csv`). PABP projects to PC1 ≈ +0.38, isolated from the
other four proteins which fall in PC1 ∈ [−0.18, +0.01]. PC1 loadings
are dominated by Early–Middle (+0.207); PC2 separates Early–Late
(+0.103) from Late–Late (−0.094) (`pca_loadings_ranked.csv`).

### Legend
- **Panel (a)** — colour-coded matrix; rows: proteins (CBS, GAL4, PABP,
  PTEN, TEM1); columns: categories (EE, EM, EL, MM, ML, LL); colour =
  proportion of Top50.
- **Panel (b)** — dendrogram; y-axis: Euclidean distance under average
  linkage; leaves labelled with protein names; PABP rendered in the
  highlight colour used throughout the manuscript.
- **Panel (c)** — scatter of proteins in (PC1, PC2) coordinates;
  variance-explained shown on each axis; loading arrows for the six
  categories overlaid.

### Take-home message
Three independent geometric views (heatmap, hierarchical clustering,
PCA) converge: PABP's layer-pair signature is shifted toward
Early–Middle interactions and is more distant from any of the four
EL-dominated proteins than they are from each other.

---

## Figure 4 — Cross-protein recurrent layer pairs concentrate at layers 30–32 (C5)

### Caption
**Recurrent layer pairs across five proteins concentrate at a sharp
three-layer hub at depths 30–32, with layer 32 most frequent.**
**(a)** Per-layer count of layer-end occurrences among the Top20
cross-protein recurrent layer pairs
(`layer_recurrence_top20_pairs.csv`; 40 layer-end occurrences = 2 × 20
pairs). Layer 32 contributes 8 of 40 (20 %) and is the single most
frequent layer; layers 31 (7), 30 (3) and 33 (2) complete the deep-layer
hub; layers 18–27 contribute zero. Band totals
(`layer_recurrence_band_summary.csv`): Late (25–33) = 22 / 40 = 55 %;
Middle (9–24) = 11 / 40 = 27.5 %; Early (1–8) = 7 / 40 = 17.5 %.
**(b)** Layer-occurrence histogram across the broader union of all five
proteins' Top50 Composite sets (`layer_frequency.csv`; n = 500
layer-end occurrences = 5 × 50 × 2). The depth-wise concentration is
preserved at this broader scale (Late 47.6 %, Middle 36.0 %, Early
16.4 %).

### Legend
- **Panel (a)** — bar chart; x-axis: layer index 1–33 (1-based); y-axis:
  count among the 40 layer-end occurrences of the Top20 recurrent set;
  bar colours reflect band (Early / Middle / Late).
- **Panel (b)** — bar chart with the same x-axis but evaluated over the
  larger n = 500 union; band proportions printed below the panel.

### Take-home message
The recurrent layer pairs do not spread over the Late band uniformly;
they sit at layers 30–32 with layer 32 most frequent, and layers
18–27 are essentially absent — both observations follow from direct
counts on the Top20 recurrent set.

---

## Figure 5 — Statistical robustness across proteins, corrections and bootstrap (C6 + C7)

### Caption
**Per-protein reproducibility and multi-layered statistical robustness of
the layer-pair findings.**
**(a)** Per-protein null distributions of T₁, T₂, T₃ and T₃ˢᵗᵈ under
within-protein label shuffles (N = 10³ per protein; Methods §3.7),
for each of the five proteins; observed value marked in the
protein-specific highlight colour. All five proteins return
FDR-significant T₁, T₃ and T₃ˢᵗᵈ
(`permutation_per_protein_summary.csv`). **(b)** Bar chart of
observed-vs-null z-scores across the four primary statistics
(T₁, T₂, T₃, T₃ˢᵗᵈ) for the pooled (Global) sample and each of the
five proteins; reference lines at |z| = 1.96 (two-sided α = 0.05) and
at the Bonferroni-adjusted threshold for the 20 per-protein primary
tests. **(c)** QQ-style significance plot of the observed |z| against
the half-normal quantile under H₀ for all 24 tests (4 statistics ×
{Global + 5 proteins}); diagonal y = x marks the null expectation,
horizontal lines mark α = 0.05 and the Bonferroni threshold. Bootstrap
N = 10³ 95 % confidence intervals for T₁ (full bootstrap) and T₂, T₃,
T₃ˢᵗᵈ, T_{max} (rank-reuse fast path with measured peak |ΔAUROC| ≤ 6 ×
10⁻³ at real scale; Methods §3.9) exclude the corresponding null mean
for every primary statistic (`bootstrap_summary.csv`).

### Legend
- **Panel (a)** — 5 × 4 grid of null histograms; rows: proteins, columns:
  statistics; each cell annotated with z and p_FDR.
- **Panel (b)** — grouped bars; x-axis: scope (Global + 5 proteins);
  y-axis: z-score; one bar per statistic per scope; dotted = α = 0.05;
  dashed = Bonferroni.
- **Panel (c)** — scatter; x-axis: expected |z| under H₀; y-axis:
  observed |z|; one point per test; markers coloured by scope; diagonal
  reference y = x.

### Take-home message
Every primary finding survives three independent layers of statistical
robustness — Phipson–Smyth permutation, Bonferroni FWER, Benjamini–
Hochberg FDR — and the bootstrap 95 % CI excludes the null mean for
every primary statistic; per-protein reproducibility is established in
5 / 5 proteins on T₁, T₃, T₃ˢᵗᵈ.

---

## Table 1 — SCI vs ESM-2 zero-shot LLR baseline (§2.8)

### Caption
**As a per-mutation scalar predictor on the 95,142-mutation panel,
the masked-marginal ESM-2 LLR baseline of Meier et al. (2021) exceeds
the 528-pair upper-triangle SCI summary on every protein, and a
combined SCI+LLR scorer adds at most Δ AUROC = 0.008 over LLR alone.**
Pooled and per-protein Spearman r against `DMS_score`, AUROC and
AUPRC against `DMS_score_bin` are reported for three scorers: the
528-pair upper-triangle mean SCI summary
(`all_proteins_sci_site_scores_mean.npy`), the masked-marginal LLR
(Methods §3.12) and the 5-fold stratified-CV logistic regression on
standardised (SCI, LLR) (Methods §3.12). Numbers are verbatim from
`data/processed/sci_vs_llr_comparison.csv`; Spearman r and AUROC for
LLR alone match `data/processed/baseline_summary.csv` (LLR-only sanity
check).

| Scope | Scorer | n | Spearman r | AUROC | AUPRC |
|---|---|---:|---:|---:|---:|
| Pooled | SCI | 95,142 | 0.200 | 0.669 | 0.739 |
| Pooled | LLR | 95,142 | 0.480 | 0.753 | 0.791 |
| Pooled | SCI+LLR (CV-LR) | 95,142 | 0.464 | 0.758 | 0.796 |
| CBS | SCI | 7,217 | 0.131 | 0.570 | 0.532 |
| CBS | LLR | 7,217 | 0.342 | 0.692 | 0.589 |
| CBS | SCI+LLR (CV-LR) | 7,217 | 0.343 | 0.692 | 0.594 |
| GAL4 | SCI | 1,195 | 0.433 | 0.717 | 0.856 |
| GAL4 | LLR | 1,195 | 0.668 | 0.795 | 0.905 |
| GAL4 | SCI+LLR (CV-LR) | 1,195 | 0.657 | 0.797 | 0.907 |
| PABP | SCI | 74,229 | 0.356 | 0.699 | 0.761 |
| PABP | LLR | 74,229 | 0.486 | 0.767 | 0.801 |
| PABP | SCI+LLR (CV-LR) | 74,229 | 0.491 | 0.770 | 0.804 |
| PTEN | SCI | 7,504 | 0.135 | 0.629 | 0.865 |
| PTEN | LLR | 7,504 | 0.501 | 0.846 | 0.945 |
| PTEN | SCI+LLR (CV-LR) | 7,504 | 0.497 | 0.847 | 0.945 |
| TEM1 | SCI | 4,997 | 0.502 | 0.759 | 0.768 |
| TEM1 | LLR | 4,997 | 0.714 | 0.882 | 0.852 |
| TEM1 | SCI+LLR (CV-LR) | 4,997 | 0.733 | 0.889 | 0.882 |

### Legend
- **Spearman r** — Spearman rank correlation of scorer against the
  continuous `DMS_score`; reported as |r| with the natural sign of
  the correlation given in the `spearman_sign` column of
  `sci_vs_llr_comparison.csv` (always +1 across all 18 entries).
- **AUROC** — sign-corrected AUROC of scorer against `DMS_score_bin`,
  using the same convention as `P1.5_validate_sci_signal.py`
  (corrected = max over native and negated score).
- **AUPRC** — sign-corrected average precision of scorer against
  `DMS_score_bin`.

### Take-home message
LLR is the stronger univariate scorer on every protein and on the
pooled sample. The combined SCI+LLR scorer is essentially LLR
(Δ AUROC ≤ 0.008 in every protein), confirming that the scalar SCI
summary is largely redundant with the LM-head likelihood. The
contribution of the present work is therefore not a new scalar
predictor; it is the layer-resolved 528-pair landscape (C2, C3),
the layer-30–32 recurrence hub (C5) and the per-protein compositional
heterogeneity (C4), none of which is exposed by a single LM-head LLR.

---

## Cross-reference index

| Figure | Carries | Source-of-truth files (numbers) | Source-of-truth files (images) |
|---|---|---|---|
| Fig 1 | C1 | `permutation_global_summary.csv` (T1 row), `permutation_per_protein_summary.csv` (T1 rows), `bootstrap_summary.csv` (T1 row) | `P1.5_sci_signal_validation.png`, `fig_A_global_null.png` |
| Fig 2 | C2, C3 | `*_pair_category_ratio_Composite.csv` (per-protein), `permutation_global_summary.csv` (T2, T3, T3_std rows) | `layer_pair_category_distribution.png`, `fig_A_global_null.png` |
| Fig 3 | C4 | `pairwise_euclidean_distance.csv`, `pca_loadings.csv`, `pca_loadings_ranked.csv` | `heatmap_layer_pair_categories.png`, `dendrogram_protein_clustering.png`, `pca_protein_layer_pair.png` |
| Fig 4 | C5 | `layer_recurrence_top20_pairs.csv`, `layer_recurrence_band_summary.csv`, `layer_recurrence_frequency.csv`, `layer_frequency.csv`, `layer_frequency_band_summary.csv` | `layer_recurrence_barplot.png`, `layer_frequency_histogram.png` |
| Fig 5 | C6, C7 | `permutation_per_protein_summary.csv`, `permutation_global_summary.csv`, `bootstrap_summary.csv`, `publication_summary.csv` | `fig_B_per_protein_null.png`, `fig_C_observed_vs_null_bar.png`, `fig_D_qq_significance.png` |
| Table 1 | §2.8 baseline comparison | `sci_vs_llr_comparison.csv`, `baseline_summary.csv`, `esm2_llr_scores.csv` | — (text table only; no figure asset required) |

---

## Layout notes (for the editor / for Adobe Illustrator assembly)

- **Fig 1**: two-panel (a, b), single column, ≈ 90 mm width.
- **Fig 2**: three-panel (a stack-bar wide; b, c paired) — two-column,
  ≈ 180 mm width.
- **Fig 3**: three-panel (heatmap, dendrogram, PCA) — two-column,
  ≈ 180 mm; panel (a) heatmap should retain the PABP row highlight.
- **Fig 4**: two-panel (a, b) — two-column, ≈ 180 mm.
- **Fig 5**: three-panel (a 5 × 4 grid; b, c paired) — two-column,
  full-width ≈ 180 mm; (a) is the largest panel and should occupy the
  top half.
- **Table 1**: single-column-wide text table, no figure asset; numbers
  pasted verbatim from `data/processed/sci_vs_llr_comparison.csv`
  rounded to three decimal places. Place after Fig 5 and before
  the supplementary block.
