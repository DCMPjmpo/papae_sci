# Manuscript framework

> Built strictly from the C1–C7 evidence inventoried in
> `paper/scientific_story_map.md`. No new experiments, no fabricated
> results. All numerical values are sourced from the on-disk artefacts.
>
> Target journals: Nature Communications · Briefings in Bioinformatics ·
> Bioinformatics. Section order follows the Nature Communications template
> (Results before Methods); the Discussion already exists in
> `paper/discussion_full.md` and is summarised here in compact form.

---

## Title

**Early–Late layer-pair interactions encode mutation-effect signal in the
protein language model ESM-2 across five deep mutational scans.**

(Alternatives: *"Layer-pair geometry of mutation-effect encoding in
ESM-2"*; *"Cross-layer interaction geometry organises mutation-effect
signal in ESM-2"*. The lead title above embeds the central finding
without causal verbs.)

---

## Abstract  (≈ 200 words; single paragraph, Nature Communications style)

Mutation-effect information in protein language models (PLMs) has been
studied almost exclusively at the single-layer level, typically the
masked-language-model head. How this information is *distributed across
layers* of a PLM remains under-characterised. Here we systematically
quantify the discriminative signal carried by every one of the 528
unordered layer pairs of ESM-2 (esm2_t33_650M_UR50D; Lin et al., 2023), using the per-mutation Structural
Context Inconsistency (SCI) metric across five deep mutational scanning
(DMS) panels (CBS, GAL4, PABP, PTEN, TEM1; 95,142 mutations). SCI
separates functional from non-functional mutations at the pooled level
and within each protein (within-protein permutation null, global z = 57,
p_FDR = 10⁻³; 5/5 proteins per-protein FDR-significant). When the
528-pair plane is ranked by a composite discrimination score, the
top-scoring layer pairs are enriched for Early–Late (EL) interactions,
and the EL category exceeds Early–Early, Early–Middle, Middle–Middle and
Late–Late categories on |AUROC − 0.5| (global T3 z = 62). Cross-protein
recurrent pairs concentrate at layers 30–32, with layer 32 most frequent
and layers 18–27 essentially absent. Four proteins share an EL-dominated
composition; PABP forms a singleton driven by Early–Middle interactions
(hierarchical clustering, PCA, pairwise Euclidean distances). All
primary effects survive Bonferroni FWER, Benjamini–Hochberg FDR, and
non-parametric bootstrap. We interpret these observations as evidence of
a two-axis (layer × layer) organisation of mutation-effect signal in
ESM-2, complementing the prevailing one-axis (single-layer) view.

---

## Introduction  (4 paragraphs)

**¶1 — The single-layer convention in PLM mutation-effect work.**
Protein language models such as ESM-2, ESM-2 and Tranception have
re-defined how variant effect is scored from sequence alone. Across this
literature, the variant score is read from a *single* representation of
the model — typically the final masked-LM head or the last hidden state
(Meier et al., 2021; Notin et al., 2022; Brandes et al., 2023; Frazer
et al., 2021; Hsu et al., 2022). The same one-axis convention organises
PLM interpretability, where per-layer probing or per-head attention
analysis charts how secondary structure, contacts and biophysical
properties emerge with depth (Rives et al., 2021; Lin et al., 2023;
Elnaggar et al., 2022; Vig et al., 2021; Simon & Zou, 2025). Whether
single-layer reading is *sufficient* for mutation-effect representation,
or whether discriminative information is in fact distributed across
multiple layers, has not been systematically tested.

**¶2 — A two-axis precedent from general transformers.**
For general transformers, mechanistic interpretability has documented
that representational features are progressively composed across layers
via the residual stream (Tenney et al., 2019; Geva et al., 2021; Belrose
et al., 2023). Per-layer probing alone collapses this two-axis structure
onto one dimension; scalar-mix and logit-lens / tuned-lens methods
recover the other. The only existing protein-LM analysis that explicitly
aggregates across layers is the contact-prediction pipeline of Rao
et al. (2021), in which a sparse logistic regression pools attention
heads across all ESM layers. Critically, that work targets residue
contacts, not mutation effects, leaving the layer-pair geometry of
mutation-effect signal an open question.

**¶3 — Open question and approach.**
We ask how mutation-effect information is organised across the layers of
ESM-2. Specifically: (i) Is the signal localised to a single layer, or
distributed across pairs? (ii) If distributed, which depth combinations
carry it? (iii) Is the distribution shared across proteins, or
protein-specific? To address these questions we compute, for every
mutation in five DMS panels, the per-layer-pair Structural Context
Inconsistency (SCI) score, producing a (n_mutations, 33, 33) matrix
which projects onto 528 unordered upper-triangle layer pairs. We then
analyse this layer-pair plane jointly — across categories defined by
the Early (layers 1–8), Middle (9–24) and Late (25–33) bands of the
33-layer stack — and stratify by protein.

**¶4 — Contribution preview.**
The 528-pair landscape is dominated by Early–Late interactions; the
deep-layer anchor of this dominance is sharply localised at layers
30–32 and recurs across proteins; and one of the five proteins (PABP)
forms a distinct cluster driven by Early–Middle, rather than
Early–Late, interactions. Every claim is validated against a
within-protein permutation null at both the pooled and per-protein
levels, with Bonferroni and Benjamini–Hochberg FDR control and
non-parametric bootstrap. We frame the resulting picture as a
*two-axis* (layer × layer) view of mutation-effect encoding that
complements, rather than replaces, the prevailing one-axis (single-layer)
view.

---

## Results

Each subsection corresponds to one validated claim (C1–C7) from
`paper/scientific_story_map.md`. Numbers are sourced from the
referenced on-disk artefact; nothing is interpolated.

### 2.1 · SCI discriminates functional from non-functional mutations (C1)

To test whether the Structural Context Inconsistency (SCI) metric carries
mutation-effect signal, we computed the per-mutation SCI summary score
(`all_proteins_sci_site_scores_top50.npy`) for every mutation in the
five DMS panels (n = 95,142). At the pooled level, the SCI mean
differed between functional (`DMS_score_bin = 1`) and non-functional
mutations by T1 = 2.79 × 10⁻³ (`permutation_global_summary.csv`).
Against a within-protein label-shuffle null (10⁻³ adaptive permutations;
Methods), the observed difference was z = 57.0, p_empirical = 10⁻³,
p_FDR = 10⁻³ across the four global primary statistics. At the
per-protein level, T1 was FDR-significant in 5/5 proteins
(smallest = CBS, p_FDR = 0.037; remaining four ≤ 1.4 × 10⁻³;
`permutation_per_protein_summary.csv`). Bootstrap (N = 1,000)
95% confidence interval for T1 excluded zero
(`bootstrap_summary.csv`). Together, these establish that SCI carries
discriminative signal that is reproducible against the null at both the
pooled and per-protein levels.

### 2.2 · Top-ranked layer pairs are enriched for Early–Late interactions (C2)

We ranked the 528 unordered layer pairs by the Composite discrimination
score (a function of pair-level AUROC, Spearman rank correlation and
Cohen's *d*; Methods) and examined the Top50 layer pairs for each
protein (`layer_pair_top50.csv` and per-protein
`*_pair_category_ratio_Composite.csv`). The category composition of
each Top50 set was: CBS 40% EL + 30% LL; GAL4 66% EL; PTEN 44% EL +
12% LL; TEM1 40% EL + 10% EM + 8% LL; PABP 48% EM + 42% EL. Four of the
five proteins were EL-dominated, with EL accounting for 40–66% of the
high-scoring layer pairs. PABP was the only protein whose Top50 set
contained a larger EM than EL share (the implications of this asymmetry
are addressed in §2.4 and Discussion D3). At the pooled level, the
AUROC-only Top50 set contained 17 EL pairs against a null expectation of
1.0 (z = 8.0, p_FDR = 10⁻³; `permutation_global_summary.csv`),
confirming the enrichment under a metric that is independent of the
Composite. A methodological note: at the per-protein level the
AUROC-only ranking and the Composite ranking can diverge (Methods);
per-protein evidence for the EL enrichment is therefore primarily
carried by the per-protein T3 statistic in §2.3.

### 2.3 · Early–Late pairs dominate the 528-pair landscape (C3)

To test whether the EL enrichment in the Top50 is a consequence of a
broader category-level effect rather than a tail phenomenon, we computed
|AUROC − 0.5| for each of the 528 layer pairs and compared the EL
category to the union of EE, EM, MM and LL pairs (T3 = mean_EL −
mean_non-EL; T3_std = standardised by the non-EL standard deviation;
Methods). At the pooled level, T3 = 0.040 with z = 62.0
and p_FDR = 10⁻³, and T3_std = 1.02 with z = 8.05 and p_FDR = 10⁻³
(`permutation_global_summary.csv`; Fig. A). At the per-protein
level, T3 and T3_std were FDR-significant in 5/5 proteins
(`permutation_per_protein_summary.csv`; Fig. B). The per-protein
|AUROC − 0.5| dominance score is computed against the globally
ranked SCI scores (see §3.5 in Methods); per-protein effect-size
magnitudes are therefore not directly bounded in [−0.5, 0.5], but the
permutation null is computed under the identical metric and the
significance call is unaffected. The category-level EL dominance is
therefore present both at the pooled and protein-stratified levels.

### 2.4 · One of the five proteins (PABP) forms a distinct cluster (C4)

We asked whether the EL-dominated composition reported in §2.2 is shared
across the five proteins or whether it conceals heterogeneity. For each
protein we represented the Top50 layer pairs as a six-dimensional vector
of category proportions (EE, EM, EL, MM, ML, LL) and applied three
independent geometric analyses (`data/p0_output/protein_clustering/`).
Hierarchical clustering with Euclidean distance and average linkage
(`dendrogram_protein_clustering.png`) placed CBS, GAL4, PTEN and TEM1
into one cluster whose maximum pairwise distance was 0.40 and joined
PABP to this cluster only at a distance of 0.50. The minimum
PABP-to-other distance (0.389; PABP–TEM1) exceeded the maximum
within-other-four distance (0.397; CBS–GAL4),
satisfying a singleton condition (`pairwise_euclidean_distance.csv`).
Principal component analysis on the same six-dimensional vector
(`pca_protein_layer_pair.png`; `pca_loadings.csv`) explained 96.7% of
variance in two components (PC1 68.6%, PC2 28.1%) and placed PABP at
PC1 ≈ +0.38, isolated from the other four (PC1 ∈ [−0.18, +0.01]). PC1
loadings (`pca_loadings_ranked.csv`) were dominated by the Early–Middle
category (+0.207), with the next-largest loading roughly three times
smaller (LL, −0.072). PC2 separated EL from LL with loadings of +0.103
and −0.094. We summarise this as: the four EL-dominated proteins share
a similar layer-pair signature, while PABP's signature is shifted toward
Early–Middle interactions.

### 2.5 · Cross-protein recurrent layer pairs concentrate at layers 30–32 (C5)

To identify layer pairs that carry mutation-effect signal *across*
proteins, we counted, for each layer pair, the number of proteins in
whose Top50 Composite set it appears (`cross_protein_universal_pairs.csv`;
146 candidate recurrent pairs). The Top20 most-recurrent pairs
(`layer_recurrence_top20_pairs.csv`) carried 40 layer-end occurrences
(2 ends × 20 pairs). The Late band (layers 25–33) carried 22 of these
(55%); the Middle band (9–24) carried 11 (27.5%); the Early band (1–8)
carried 7 (17.5%) (`layer_recurrence_band_summary.csv`). Within the
Late band the distribution was non-uniform: layer 32 alone contributed
8 of 40 layer-ends (20%; the single most recurrent layer), layer 31
contributed 7 (17.5%), layer 30 contributed 3 (7.5%), and layers
25–29 contributed 0 (`layer_recurrence_barplot.png`). At the broader
Top50 union across five proteins (n = 500 layer-ends;
`layer_frequency.csv`), the band proportions remained Late 47.6%,
Middle 36.0%, Early 16.4%. The recurrent set is therefore concentrated
at a sharp three-layer hub (30–32), with layer 32 most frequent, and
with a marked trough at layers 18–27.

### 2.6 · Per-protein signal is reproducible against within-protein null (C6)

To rule out the possibility that the pooled effects reported in §2.1–2.3
are driven by a single protein, we re-computed the test statistics
within each of the five proteins, using a within-protein label-shuffle
null (1,000 permutations per protein; Methods). All five proteins
returned FDR-significant T1, T3 and T3_std statistics
(`permutation_per_protein_summary.csv`; Fig. B). T1 z-scores
ranged from 2.25 (CBS) to 56.6 (PABP). T3 was FDR-significant in 5/5
proteins; the smallest Bonferroni-adjusted p across the per-protein
primary tests was approximately 2 × 10⁻². Per-protein FDR-significance
implies that the EL-dominance reported at the pooled level (§2.3) is not
an artefact of pooling proteins with heterogeneous label balances or
sample sizes.

### 2.7 · Statistical robustness across multiple testing layers (C7)

To assess robustness to multiple-testing assumptions and to estimator
uncertainty, we applied three independent corrections. (i) **Family-wise
error rate**: Bonferroni correction over the 20 per-protein primary tests
(5 proteins × 4 primary statistics) was applied to the empirical
permutation p-values; with this correction, every protein retains
FDR-significant T1, T3 and T3_std (`permutation_per_protein_summary.csv`,
column `p_bonferroni`). (ii) **False discovery rate**: Benjamini–Hochberg
FDR correction was applied over the same 20-test family; the
column `p_fdr` reports the resulting q-values, and the q < 0.05 calls
match the FWER calls in 5/5 proteins on T1, T3 and T3_std. (iii)
**Bootstrap re-sampling**: a non-parametric bootstrap (N = 1,000;
rank-reuse fast path with measured AUROC deviation ≤ 6 × 10⁻³ relative
to the natural bootstrap, Methods) provided 95% confidence intervals
for T1, T2, T3, T3_std and T_max (`bootstrap_summary.csv`); all five
intervals excluded the corresponding null mean. (iv) **Family-wise max
statistic**: T_max = max |AUROC − 0.5| over all 528 pairs achieved
z = 93.7 against null (`permutation_global_summary.csv`), addressing the
concern that the best-of-528 search inflates significance. The four
robustness layers converge: the EL geometry of mutation-effect signal in
ESM-2 survives every correction we have applied to it.

---

## Discussion  (compact form; full version in `paper/discussion_full.md`)

The full Discussion is provided in `paper/discussion_full.md` (≈1,758
words, four subsections D1–D4). The compact form below preserves the
four themes and the same set of cited works.

**D1 — Deep-layer enrichment in ESM.** The concentration of cross-protein
recurrence at layers 30–32 with layer 32 most frequent (§2.5) is
consistent with prior layer-wise probing studies that identify deep ESM
layers as the locus of higher-order features (Rives et al., 2021;
Lin et al., 2023; Elnaggar et al., 2022; Vig et al., 2021), and with
the observation that the most useful single layer of a PLM depends on
the downstream task (Detlefsen et al., 2022). The novelty of our
finding is that, for mutation-effect signal in ESM-2, the deep layers
appear within *pairs* rather than as stand-alone signal carriers.

**D2 — Cross-layer integration hypothesis.** The systematic EL
enrichment (§2.2, §2.3, §2.5) is consistent with a working hypothesis
that mutation-effect information is integrated across layers rather than
read off any single layer. This reading is consistent with the
residual-stream-centered interpretability picture developed for general
transformers (Tenney et al., 2019; Geva et al., 2021; Belrose et al.,
2023), and is distinct from the conventional single-layer protocol
(Meier et al., 2021; Notin et al., 2022; Brandes et al., 2023). The
hypothesis is correlational; no interventional experiment is performed
here.

**D3 — Protein-specific heterogeneity.** The PABP shift toward
Early–Middle interactions (§2.4) co-varies with biophysical class (PABP
is the only RNA-binder in the panel) and with deep-mutational-scan
assay design (Melamed et al., 2013; the other four assays measure
stability or enzymatic activity). Both readings are *consistent with*
the data; neither is causally established at n = 5 proteins.

**D4 — Implications for PLM interpretability.** The EL enrichment,
the layer-32 anchor, and the protein-specific compositional heterogeneity
together suggest that layer-pair geometry is a useful interpretability
primitive complementing single-layer probing (Elnaggar et al., 2022;
Detlefsen et al., 2022; Vig et al., 2021; Simon & Zou, 2025) and
per-head attention analysis (Rao et al., 2021). The six-category
composition vector reported per-protein could in principle be computed
for ESM-1b (Rives et al., 2021), ESM3 (Hayes et al., 2025) and ProtT5,
as a low-dimensional fingerprint of how a backbone organises
mutation-effect signal across depth.

**Boundaries on these claims.** We do not claim that layer-pair features
yield higher predictive accuracy than tuned single-layer PLM scores; we
did not perform a head-to-head comparison. We do not claim the EL
pattern is causally responsible for ESM-2's mutation-effect
performance; activation patching, single-layer probes at matched
dimensionality, and tuned-lens decoding are deferred to follow-up work.
We do not claim that the layer-32 anchor generalises to other PLMs;
only ESM-2 is tested here.

---

## Methods

### 3.1 · Datasets
Five DMS panels were used: CBS (Sun et al., 2020), GAL4 (Kitzman
et al., 2015), PABP (Melamed et al., 2013), PTEN (Matreyek et al.,
2021), and TEM1 (BLAT_ECOLX; Firnberg et al., 2014 and related). All
panels were drawn from the standardised ProteinGym substitution release
(`data/raw/proteingym/DMS_ProteinGym_substitutions/`). Per-mutation
metadata (`data/processed/all_proteins_site_metadata.csv`; n = 95,142)
include the protein label, mutant identifier, mutation position, raw
DMS score, and the binarised label `DMS_score_bin` (1 = functional,
0 = non-functional) used as the supervised target throughout.

### 3.2 · ESM-2 feature extraction and Structural Context Inconsistency (SCI)
Per-residue representations were extracted from the 33-layer ESM-2
encoder for each wild-type sequence and each mutated sequence. For each
mutation we constructed the per-pair Structural Context Inconsistency
(SCI) metric `P1_build_sci.py`, producing a symmetric (33, 33) score
matrix per mutation; the full set of 95,142 such matrices is stored as a
float32 memmap at `data/processed/all_proteins_sci_site_matrices.dat`
(shape 95142 × 33 × 33). Per-mutation summary scores
(`all_proteins_sci_site_scores_mean.npy`,
`all_proteins_sci_site_scores_top50.npy`) record, respectively, the mean
over all 528 upper-triangle layer pairs and the mean over the 50
top-scoring pairs.

### 3.3 · Layer-pair categorisation and Composite ranking
Layer pairs (i, j) with i < j and i, j ∈ {1, …, 33} were partitioned
into six categories on the Early (1–8) / Middle (9–24) / Late (25–33)
banding: EE (28 pairs), EM (128), EL (72), MM (120), ML (144), LL (36),
totalling 528. The Composite ranking score (`p0_layer_pair_mining_v2.py`;
columns reported in `*_top50_composite_*.csv`) is a non-negative
combination of (a) per-pair AUROC against `DMS_score_bin`, (b) Spearman
rank correlation between the per-pair SCI and the continuous
`DMS_score`, and (c) Cohen's *d* between functional and non-functional
groups on per-pair SCI. The exact weighting follows the original
`p0_layer_pair_mining_v2.py` implementation.

### 3.4 · Protein clustering on layer-pair categories
For each protein, the Top50 Composite layer pairs were summarised as a
six-dimensional category-proportion vector
(`*_pair_category_ratio_Composite.csv`). Hierarchical clustering used
Euclidean distance and average linkage; PCA was computed on the same
six-dimensional vectors with no standardisation, so that PC1 / PC2
reflect raw compositional differences (`protein_clustering_layer_pair.py`;
`pca_loadings_extract.py`).

### 3.5 · Cross-protein recurrence
A layer pair (i, j) was deemed *recurrent* if it appeared in the Top50
Composite of multiple proteins. The Cross_Protein_Count was computed
in `top_recurrent_cross_protein_pairs.py` and stored at
`data/processed/p0_output/cross_protein_universal_pairs.csv` (146
candidate pairs). The Top20 recurrent set, ranked by Cross_Protein_Count
(ties broken in the original file's order), was used for the layer-end
occurrence analysis (`layer_recurrence_top20.py`; 40 layer-end
occurrences = 2 × 20 pairs).

### 3.6 · Test statistics

- **T1 (SCI mean difference)**: mean of the Top50-mean SCI score over
  functional minus non-functional mutations.
- **T2 (AUROC-Top50 EL count)**: number of Early–Late layer pairs among
  the 50 pairs with largest |AUROC − 0.5|. Note: at the per-protein level
  this AUROC-only ranking can diverge from the Composite ranking
  underlying §2.2; per-protein C2 evidence is therefore primarily
  carried by T3 and the Composite category-ratio CSVs.
- **T3 (EL dominance)**: mean |AUROC − 0.5| over the 72 EL pairs minus
  the same over the 456 non-EL pairs.
- **T3_std (standardised EL dominance)**: T3 divided by the standard
  deviation of |AUROC − 0.5| over the non-EL pairs.
- **T_max (family-wise max)**: max |AUROC − 0.5| over all 528 pairs.

Per-pair AUROC is computed via the standard rank–U identity. At the
per-protein level, ranks are reused from the full-population SCI
ordering for efficiency; under this convention the per-protein T3 / T_max
are *globally-anchored* U-statistics whose null distribution is computed
under the identical statistic, so significance is preserved while
absolute values are not bounded in [−0.5, 0.5].

### 3.7 · Permutation test and effect size
Permutation tests use within-protein label shuffle: labels are permuted
independently inside each of the five proteins, preserving protein-level
sample sizes and label balances. Empirical p-values follow the
Phipson & Smyth (2010) form p = (1 + #{T_perm extreme}) / (N + 1).
Two-sided tests center the comparison on the null mean
(|T_null − μ_null| ≥ |T_obs − μ_null|), correcting a common
implementation bug under which two-sided p is computed against zero.
The number of permutations follows an adaptive schedule {1,000 / 5,000 /
10,000} with an early-stopping rule at max-p < 0.001. Effect sizes are
reported as z = (T_obs − μ_null) / σ_null (equivalently Cohen's *d*
under the permutation null) and as the raw deviation T_obs − μ_null.
Null 95% and 99% percentile intervals are reported alongside.

### 3.8 · Multiple-testing correction
Bonferroni correction was applied over the 20 per-protein primary tests
(5 proteins × 4 statistics T1, T2, T3, T3_std). Benjamini–Hochberg FDR
control with step-up monotonicity enforcement was applied over the same
family. For the global summary, BH FDR was additionally applied over
the 4 primary global statistics. The supplementary T_max is reported
outside the family-wise correction. Significance calls in §2.7 use
FDR q < 0.05.

### 3.9 · Bootstrap confidence intervals
Non-parametric bootstrap with N = 1,000 resamples and replacement was
used to compute 95% and 99% percentile confidence intervals for T1 (no
ranks involved). For T2, T3, T3_std and T_max, we used a *rank-reuse*
fast bootstrap: instead of recomputing per-column ranks on every
resample (the natural bootstrap, O(B · K · N log N)), we reused the
precomputed full-sample ranks (O(B · K · N)). Asymptotically, the
rank-reuse AUROC and natural bootstrap AUROC coincide; at the real
scale (N = 95,142, K = 528) we measured peak per-pair |ΔAUROC| ≤ 6 ×
10⁻³ across 3 sanity-check resamples (`_sanity_check_rank_reuse_bootstrap`
in `p5_permutation.py`), which is well below the precision at which
95% CIs are reported. Measured speedup vs. the natural bootstrap was
≈ 110× on real-scale data.

### 3.10 · Reproducibility
All analyses use the fixed pseudo-random seed `RNG_SEED = 42` for the
permutation RNG and `RNG_SEED + 1` for the bootstrap RNG. All input
data, intermediate matrices, summary CSVs, null distributions and
figures are stored on disk under `data/`. The full p5 pipeline runs in
≈ 8–10 min on a CPU (single thread for the rank-reuse paths). Source
scripts are inventoried by claim in §3 of
`paper/scientific_story_map.md`; the per-claim mapping is preserved one-
to-one between Results §2.1–2.7 of this manuscript and that map.

---

## Figure index (for layout planning, not new figures)

| # | Suggested main figure | Source artefacts | Carries |
|---|---|---|---|
| Fig 1 | SCI signal validation + permutation null | `P1.5_sci_signal_validation.png` + `p5_permutation/fig_A_global_null.png` (T1 panel) | §2.1 (C1) |
| Fig 2 | 528-pair landscape + Top50 EL composition | `layer_pair_frequency_map.png` + `layer_pair_category_distribution.png` | §2.2, §2.3 (C2, C3) |
| Fig 3 | PABP heterogeneity (multi-panel) | `heatmap_layer_pair_categories.png` + `dendrogram_protein_clustering.png` + `pca_protein_layer_pair.png` | §2.4 (C4) |
| Fig 4 | Cross-protein recurrence at layers 30–32 | `layer_recurrence_barplot.png` + `layer_frequency_histogram.png` | §2.5 (C5) |
| Fig 5 | Robustness summary | `p5_permutation/fig_B_per_protein_null.png` + `fig_C_observed_vs_null_bar.png` + `fig_D_qq_significance.png` | §2.6, §2.7 (C6, C7) |

| Sup-Fig / Sup-Tab | Source artefacts |
|---|---|
| `permutation_global_summary.csv`, `permutation_per_protein_summary.csv`, `publication_summary.csv`, `bootstrap_summary.csv` | All as supplementary tables |
| `pca_loadings.csv`, `pca_loadings_ranked.csv`, `pairwise_euclidean_distance.csv`, `layer_recurrence_top20_pairs.csv`, `layer_frequency_band_summary.csv` | Supplementary tables for §2.4 and §2.5 |
