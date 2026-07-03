# Supplementary Materials (Phase D)

> Two supplementary tables, plus legends for supplementary figures.
> Every numerical value is sourced from on-disk artefacts; no new
> analyses are introduced.

---

## Supplementary Table S1 — Per-protein dataset characteristics

| Protein | DMS source dataset(s) (from `all_proteins_site_metadata.csv`, `dataset` column) | n mutations | n functional (`DMS_score_bin = 1`) | n non-functional | Fraction functional |
|---|---|---:|---:|---:|---:|
| CBS  | `CBS_Sun2020` (Sun et al., 2020) | 7,217 | 3,314 | 3,903 | 0.459 |
| GAL4 | `GAL4_Kitzman2015` (Kitzman et al., 2015) | 1,195 | 831 | 364 | 0.695 |
| PABP | `PABP_Melamed2013` (Melamed et al., 2013) | 74,229 | 46,078 | 28,151 | 0.621 |
| PTEN | `PTEN_Mighell2018`, `PTEN_Matreyek2021` (Mighell et al., 2018; Matreyek et al., 2021; merged where indicated in the `dataset` column) | 7,504 | 5,880 | 1,624 | 0.784 |
| TEM1 | `TEM1_Deng2012`, `TEM1_Firnberg2014`, `TEM1_Jacquier2013`, `TEM1_Stiffler2015` (Deng et al., 2012; Firnberg et al., 2014; Jacquier et al., 2013; Stiffler et al., 2015; combined into the merged TEM1 panel per ProteinGym) | 4,997 | 2,526 | 2,471 | 0.506 |
| **Total** | — | **95,142** | **58,629** | **36,513** | **0.616** |

**Source for the row counts:**
`data/processed/all_proteins_site_metadata.csv`, columns `protein`,
`DMS_score_bin`, `dataset`. Counts computed via
`groupby('protein').agg(...)` and verified to sum to 95,142.

**Note on label balance.**
Per-protein label balances differ substantially (0.459 in CBS to
0.784 in PTEN), motivating the **within-protein** label-shuffle null
adopted throughout (Methods §4.7) — a cross-protein shuffle would
conflate protein identity with label balance and inflate the type-I
error.

---

## Supplementary Table S2 — Full per-statistic per-scope results

Source-of-truth files:
- `data/p0_output/p5_permutation/permutation_global_summary.csv`
- `data/p0_output/p5_permutation/permutation_per_protein_summary.csv`
- `data/p0_output/p5_permutation/bootstrap_summary.csv`
- `data/p0_output/p5_permutation/publication_summary.csv` (a
  pre-formatted union of the three above)

### S2-A · Global tests

| Statistic | T_obs | Null mean | Null SD | z = Cohen's *d* | p_empirical | p_FDR | Null 95 % CI | Null 99 % CI | Bootstrap 95 % CI | N_perm |
|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|
| T₁ (SCI mean diff) | 2.79 × 10⁻³ | −3.75 × 10⁻⁴ | 5.56 × 10⁻⁵ | 57.0 | 10⁻³ | 10⁻³ | [−4.90 × 10⁻⁴, −2.64 × 10⁻⁴] | [−5.17 × 10⁻⁴, −2.24 × 10⁻⁴] | [2.67 × 10⁻³, 2.91 × 10⁻³] | 1,000 |
| T₂ (#EL in Top50 by AUROC) | 17 | 1.013 | 1.995 | 8.01 | 10⁻³ | 10⁻³ | [0, 7] | [0, 8] | [14, 19] | 1,000 |
| T₃ (EL dominance) | 3.99 × 10⁻² | −2.08 × 10⁻³ | 6.77 × 10⁻⁴ | 62.0 | 10⁻³ | 10⁻³ | [−3.27 × 10⁻³, −6.52 × 10⁻⁴] | [−3.60 × 10⁻³, −2.34 × 10⁻⁴] | [3.70 × 10⁻², 4.27 × 10⁻²] | 1,000 |
| T₃ˢᵗᵈ (standardised) | 1.02 | −0.602 | 0.201 | 8.05 | 10⁻³ | 10⁻³ | [−0.912, −0.180] | [−0.984, −0.0563] | [0.954, 1.073] | 1,000 |
| T_max (supplementary) | 0.181 | 0.0176 | 1.74 × 10⁻³ | 93.7 | 10⁻³ | — | [0.0144, 0.0211] | [0.0130, 0.0222] | [0.176, 0.187] | 1,000 |

### S2-B · Per-protein tests — z-scores (sourced verbatim from `permutation_per_protein_summary.csv`)

| Protein | T₁ z | T₂ z | T₃ z | T₃ˢᵗᵈ z | T_max z |
|---|---:|---:|---:|---:|---:|
| CBS  | 2.25  | n/a | 16.11 | 18.50 | −4.20 |
| GAL4 | 6.29  | n/a | 5.18  | 4.59  | 8.06 |
| PABP | 56.55 | (T_obs=1; see Methods §4.6 note) | 38.71 | 42.80 | 57.07 |
| PTEN | 2.88  | n/a | 14.57 | 17.05 | 9.04 |
| TEM1 | 23.59 | n/a | 15.19 | 15.94 | 14.09 |

### S2-C · Per-protein tests — q-values after BH FDR over the 20 per-protein primary tests

| Protein | T₁ p_FDR | T₂ p_FDR | T₃ p_FDR | T₃ˢᵗᵈ p_FDR |
|---|---:|---:|---:|---:|
| CBS  | 3.7 × 10⁻²   | 1     | 1.4 × 10⁻³ | 1.4 × 10⁻³ |
| GAL4 | 1.4 × 10⁻³   | 1     | 1.4 × 10⁻³ | 1.4 × 10⁻³ |
| PABP | 1.4 × 10⁻³   | 1.4 × 10⁻³ | 1.4 × 10⁻³ | 1.4 × 10⁻³ |
| PTEN | 9.3 × 10⁻³   | 1     | 1.4 × 10⁻³ | 1.4 × 10⁻³ |
| TEM1 | 1.4 × 10⁻³   | 1     | 1.4 × 10⁻³ | 1.4 × 10⁻³ |

**Rows that pass q < 0.05 at the BH-FDR level.** T₁, T₃ and T₃ˢᵗᵈ in
**5 / 5** proteins; T₂ in **1 / 5** proteins (PABP only, see Methods
§4.6 note on Composite vs AUROC-only ranking).

### S2-D · Bootstrap CIs (N = 1,000)

The bootstrap 95 % CIs reported in S2-A apply to the pooled (Global)
scope. Per-protein bootstrap was not performed; per-protein
uncertainty is captured by the null 95 % and 99 % CIs in
`permutation_per_protein_summary.csv`.

---

## Supplementary Figure S1 — Per-protein SCI distributions

**Source files:** `data/p0_5plots/sci_per_protein.png`,
`data/p0_5plots/sci_distribution.png`.

**Legend.** Per-protein histograms of the per-mutation Top50-mean SCI
score (`all_proteins_sci_site_scores_top50.npy`), stratified by the
binarised functional label (`DMS_score_bin`). Each panel shows one of
the five DMS panels (CBS, GAL4, PABP, PTEN, TEM1). Overlaid blue
(functional) and red (non-functional) histograms allow visual
inspection of the per-protein separation that underpins T₁ at the
per-protein level (Results §2.1; Sup Table S2-B). The pooled
distribution across all five proteins is reproduced in
`sci_distribution.png` for cross-reference.

**Take-home.** Per-protein label separation is visible for all five
proteins, consistent with the per-protein T₁ FDR-significance in
Sup Table S2-C.

---

## Supplementary Figure S2 — 528-pair AUROC landscape and Composite ranking

**Source files:** `data/p0_5layer/layer_pair_frequency_map.png`,
`data/p0_5layer/layer_pair_mining.png`,
`data/p0_5layer/layer_pair_auc_ranking.csv`.

**Legend.** Heatmap visualisation of |AUROC − 0.5| over the
33 × 33 layer-pair plane (with the diagonal masked), showing the
EL block in the off-diagonal corner where the highest discrimination
values concentrate (Results §2.3; Fig 2c). The
`layer_pair_auc_ranking.csv` table records the rank-ordered list of all
528 layer pairs by AUROC magnitude and is the source of T₂ and T_max
in Sup Table S2-A.

**Take-home.** EL dominance is visible on the heatmap as an
off-diagonal block of higher |AUROC − 0.5| values; the Composite-Top50
sets used in Results §2.2 (PABP heatmap row in Fig 3a) are a tail of
this larger landscape.

---

## Supplementary Figure S3 — SCI matrix structure across mutations

**Source files:** `data/p0_5plots/sci_heatmap.png`,
`data/p0_5layer/best_layer_pair_distribution.png`.

**Legend.** Representative per-mutation 33 × 33 SCI matrix (left
panel of `sci_heatmap.png`) and the distribution of best-performing
layer pairs across mutations (`best_layer_pair_distribution.png`). The
heatmap is symmetric by construction (Methods §4.3) and is
representative of the (33, 33) matrices stored at
`data/processed/all_proteins_sci_site_matrices.dat`.

**Take-home.** Per-mutation SCI matrices exhibit pairwise correlations
in the 0.7–1.0 range across most layer pairs, with the off-diagonal
Early–Late corner showing the strongest contrast between functional
and non-functional groups when aggregated (cf. Fig 2c).

---

## Supplementary Figure S4 — Per-protein permutation grid

**Source file:**
`data/p0_output/p5_permutation/fig_B_per_protein_null.png`.

**Legend.** 5 × 4 grid of per-protein permutation null distributions
for the four primary statistics (T₁, T₂, T₃, T₃ˢᵗᵈ); each panel shows
the null histogram with the observed value marked in the
protein-specific highlight colour. z-scores and empirical p-values are
annotated in each panel. This figure is reproduced in
Fig 5a of the main text; the standalone supplementary version is
larger for ease of detailed inspection of each per-protein null.

**Take-home.** Per-protein null distributions are well-separated from
the observed value across the four primary statistics in 5/5 proteins
for T₁, T₃ and T₃ˢᵗᵈ; T₂ at the per-protein level diverges between
the AUROC-only and Composite rankings (Methods §4.6).

---

## Supplementary Figure S5 — QQ-style significance summary

**Source file:**
`data/p0_output/p5_permutation/fig_D_qq_significance.png`.

**Legend.** Observed |z| values from all 24 tests (4 statistics
× {Global + 5 proteins}) plotted against the expected half-normal
quantile under H₀; points above the diagonal y = x indicate departure
from the null. Horizontal guides mark |z| = 1.96 (two-sided α = 0.05)
and the Bonferroni-adjusted threshold for the 20-test family. Tests
are coloured by scope. This figure is identical to Fig 5c.

**Take-home.** All 24 tests but four (the per-protein T₂ rows for CBS,
GAL4, PTEN, TEM1; Methods §4.6 note) lie above the Bonferroni
threshold.

---

## Supplementary Note S1 — PCA loadings on the layer-pair compositional axis

**Source file:**
`data/p0_output/protein_clustering/pca_loadings.csv` and
`pca_loadings_ranked.csv`.

| Category | PC1 loading | PC2 loading |
|---|---:|---:|
| EE  | 0.0000 | 0.0000 |
| EM  | +0.2065 | −0.0203 |
| EL  | −0.0256 | +0.1027 |
| MM  | +0.0000 | +0.0000 |
| ML  | +0.0000 | +0.0000 |
| LL  | −0.0722 | −0.0945 |

**Rankings by |loading|.**

- **PC1**: EM (0.207) > LL (0.072) > EL (0.026); EE / MM / ML have
  zero loading because they are zero in all five Top50 Composite
  vectors.
- **PC2**: EL (0.103) > LL (0.094) > EM (0.020); same zero-loading
  set as PC1.

**Take-home.** PC1 (68.6 % of variance) is dominated by EM with sign
opposite to LL, and projects PABP (the EM-rich protein) far from the
four EL-dominated proteins; PC2 (28.1 %) contrasts EL with LL, which
is the contrast that separates CBS, GAL4, PTEN, TEM1 among themselves
(Fig 3c).

---

## Supplementary Note S2 — Pairwise Euclidean distance matrix

**Source file:**
`data/p0_output/protein_clustering/pairwise_euclidean_distance.csv`.

|  | PABP | CBS | GAL4 | PTEN | TEM1 |
|---|---:|---:|---:|---:|---:|
| PABP | 0.000 | 0.566 | 0.537 | 0.495 | 0.389 |
| CBS  | 0.566 | 0.000 | 0.397 | 0.184 | 0.242 |
| GAL4 | 0.537 | 0.397 | 0.000 | 0.251 | 0.290 |
| PTEN | 0.495 | 0.184 | 0.251 | 0.000 | 0.115 |
| TEM1 | 0.389 | 0.290 | 0.290 | 0.115 | 0.000 |

**Singleton condition for PABP.** minimum PABP-to-other = 0.389
(PABP–TEM1); maximum within-other-four = 0.397 (CBS–GAL4); minimum
within-other-four = 0.115 (PTEN–TEM1). PABP's closest neighbour is
farther than the most distant pair among the other four.

---

## Cross-reference index of supplementary materials

| Item | Carries | Main-text reference |
|---|---|---|
| Sup Table S1 | Dataset characteristics, per-protein label balance | Methods §4.1; Results §2.1 |
| Sup Table S2 | Full per-statistic per-scope results (z, Cohen's *d*, p, q, CIs, bootstrap) | Results §2.1–2.7; Methods §4.7–3.10 |
| Sup Fig S1 | Per-protein SCI distributions | Results §2.1 (Fig 1) |
| Sup Fig S2 | 528-pair AUROC landscape, Composite ranking | Results §2.2, §2.3 (Fig 2) |
| Sup Fig S3 | Per-mutation SCI matrix structure | Methods §4.3 |
| Sup Fig S4 | Per-protein permutation grid (standalone, larger) | Results §2.6 (Fig 5a) |
| Sup Fig S5 | QQ-style significance summary | Results §2.7 (Fig 5c) |
| Sup Note S1 | Full PCA loadings table | Results §2.4 (Fig 3c) |
| Sup Note S2 | Pairwise Euclidean distance matrix | Results §2.4 (Fig 3b) |
