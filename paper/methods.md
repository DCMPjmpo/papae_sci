# Methods (Phase C)

> Methods text for the manuscript. Section numbering mirrors the
> Results sections in `paper/results.md` and the Methods anchors used
> throughout `paper/figure_map.md` and `paper/manuscript_framework.md`.
> All implementation details are sourced from the actual scripts under
> `p6_scripts/` and the actual on-disk files under `data/`.

---

## 3.1 Datasets

Five deep mutational scanning (DMS) panels covering one yeast
RNA-binding protein (PABP), one fungal transcription factor (GAL4) and
three enzymes (CBS, PTEN, TEM1) were drawn from the ProteinGym
substitution release
(`data/raw/proteingym/DMS_ProteinGym_substitutions/`; Notin et al.,
2023). The mutations were merged into a single panel of 95,142 single
amino-acid substitutions covering the five proteins
(`data/processed/all_proteins_site_metadata.csv`;
columns include `protein`, `mutant`, `parent_mutant`,
`mutation_position`, `n_mutations`, `DMS_score`, `DMS_score_bin`,
`dataset`). The binarised functional label `DMS_score_bin` (1 =
functional, 0 = non-functional) follows the threshold provided by
ProteinGym for each panel and is used throughout as the supervised
target. Per-protein counts of functional and non-functional mutations
are listed in Supplementary Table S1. No additional filtering was
applied; the full set of 95,142 mutations was used in every analysis.

## 3.2 ESM-2 feature extraction

We used the publicly released `esm2_t33_650M_UR50D` model
(33 Transformer blocks, hidden dimension 1280, maximum sequence length
1022; Lin et al., 2023) for all feature extraction
(`p6_scripts/config.py` lines 65–68:
`ESM_MODEL_NAME = "esm2_t33_650M_UR50D"`,
`ESM_N_LAYERS = 33`, `ESM_DIM = 1280`). For each unique wild-type
sequence and for each single-mutation variant we extracted the
per-residue hidden representations at all 33 layers
(`p6_scripts/extract_site.py`). At the mutation position only, we
computed the per-layer difference vector
**Δ**ᵢ = h_mut(layer i, pos) − h_WT(layer i, pos) ∈ ℝ¹²⁸⁰,
stacked across the 33 layers into a per-mutation matrix
Δ ∈ ℝ³³ ˣ ¹²⁸⁰ stored at
`data/processed/all_proteins_delta_site.dat` (memmap float16). No
attention pooling and no averaging across positions is used; SCI is
defined strictly at the per-residue mutation site.

## 3.3 Structural Context Inconsistency (SCI)

The Structural Context Inconsistency (SCI) metric, defined in
`p6_scripts/P1_build_sci.py`, is the Pearson correlation between
per-layer Δ-embedding vectors at the mutation position. For each
mutation we compute the 33 × 33 symmetric matrix

  SCI[i, j] = corr(Δᵢ, Δⱼ)

via `np.corrcoef(Δ)` on the (33, 1280) matrix
(`p6_scripts/P1_build_sci.py`, lines 116–132). Pearson correlation
values are well-defined when both layer-difference vectors have
non-zero variance; samples with zero-variance layers are detected
explicitly (`zero_var_sample_cnt`, line 118), and undefined or
non-finite correlations are mapped to zero by `np.nan_to_num` (line
130), without addition of random noise. The full set of 95,142 SCI
matrices is stored as a float32 memmap at
`data/processed/all_proteins_sci_site_matrices.dat`
(shape 95,142 × 33 × 33; ≈ 414 MB). The 528 unique upper-triangle
entries SCI[i, j] with i < j define the per-pair SCI features for
downstream analysis. Three per-mutation summary scores are also
materialised: the mean of all 528 upper-triangle values
(`all_proteins_sci_site_scores_mean.npy`), the mean of the top-20
(`*_top20.npy`) and the mean of the top-50 values
(`*_top50.npy`); all three are float32 (95,142,). T₁ (Methods §3.6)
uses the Top50 summary score.

## 3.4 Layer-pair categorisation and Composite ranking

Layer pairs (i, j) with i < j and i, j ∈ {1, …, 33} were partitioned
into six categories based on the Early (layers 1–8) / Middle (9–24) /
Late (25–33) banding: Early-Early (EE; 28 pairs), Early-Middle
(EM; 128 pairs), Early-Late (EL; 72 pairs), Middle-Middle (MM; 120
pairs), Middle-Late (ML; 144 pairs) and Late-Late (LL; 36 pairs);
total = 528 (`p6_scripts/p5_permutation.py`, lines 65–86). The
Composite layer-pair score, implemented in
`p6_scripts/p0_layer_pair_mining_v2.py`, is a non-negative combination
of three per-pair quantities computed against the binarised functional
label:

- AUROC of per-pair SCI as predictor of `DMS_score_bin`;
- Spearman rank correlation between per-pair SCI and the continuous
  `DMS_score`;
- Cohen's *d* between functional and non-functional groups on per-pair
  SCI.

The Composite values per protein are stored in
`data/processed/p0_output/{PROTEIN}/top50_composite_{PROTEIN}.csv`
together with the per-component scores; the across-protein union of
Composite-Top50 sets is in `data/p0_5layer/layer_pair_top50.csv`. The
exact arithmetic combination follows the original
`p0_layer_pair_mining_v2.py` implementation and is not redefined here.
The six-category proportions for each protein's Composite-Top50 set are
stored in `*_pair_category_ratio_Composite.csv`.

## 3.5 Protein clustering and cross-protein recurrence

Per-protein Composite-Top50 sets were summarised as
six-dimensional category-proportion vectors. Hierarchical clustering
used Euclidean distance with average linkage, computed in
`p6_scripts/protein_clustering_layer_pair.py`; pairwise distances are
in `pairwise_euclidean_distance.csv`. Principal component analysis
was applied to the same six-dimensional vectors without
standardisation, so that PC1 and PC2 reflect raw compositional
differences rather than rescaled units
(`pca_loadings.csv`, `pca_loadings_ranked.csv`).

For cross-protein recurrence, a layer pair (i, j) is *recurrent* if it
appears in the Composite-Top50 set of multiple proteins; the
Cross_Protein_Count statistic was computed in
`top_recurrent_cross_protein_pairs.py`
(`cross_protein_universal_pairs.csv`; 146 candidate recurrent pairs).
The Top20 recurrent set was extracted by sorting on
Cross_Protein_Count in descending order with ties broken in the
original file's order (`p6_scripts/layer_recurrence_top20.py`).
Layer-end occurrence counts per layer, per-band totals and the
broader Top50 union across the five proteins (n = 500 layer-ends =
5 × 50 × 2) are stored in `layer_recurrence_*.csv` and
`layer_frequency_*.csv`.

## 3.6 Test statistics

Five test statistics were defined to validate findings F1–F5 at both
the pooled and per-protein levels (`p6_scripts/p5_permutation.py`,
lines 174–218):

- **T₁ — SCI mean difference (validates C1).**
  T₁ = mean(SCI_top50 | y = 1) − mean(SCI_top50 | y = 0).
- **T₂ — AUROC-Top50 EL count (validates C2 at pooled level).**
  Number of EL layer pairs among the 50 with largest |AUROC − 0.5|
  per the per-pair AUROC computed against `DMS_score_bin`.
- **T₃ — EL category dominance (validates C3).**
  T₃ = mean |AUROC − 0.5|_EL − mean |AUROC − 0.5|_non-EL.
- **T₃ˢᵗᵈ — standardised EL dominance.**
  T₃ˢᵗᵈ = (mean |AUROC − 0.5|_EL − mean |AUROC − 0.5|_non-EL) /
  σ(|AUROC − 0.5|_non-EL).
- **T_max — family-wise maximum.**
  T_max = max |AUROC − 0.5| over all 528 pairs; reported as a
  supplementary control for the 528-pair search.

Per-pair AUROC is computed via the Mann–Whitney U identity using
average ranks (`scipy.stats.rankdata` with `method = "average"`;
implementation in `aurocs_from_ranks`, lines 164–171 of
`p5_permutation.py`). |AUROC − 0.5| is used to make the discrimination
test direction-agnostic, so that both functional-up and
functional-down layer pairs contribute. At the per-protein level, the
ranks of SCI values are reused from the full-population ordering for
computational efficiency; under this convention the per-protein T₃,
T₃ˢᵗᵈ and T_max are *globally-anchored* U-statistics. Significance
calls are unaffected (the permutation null is computed under the
identical statistic), but per-protein effect-size magnitudes are not
bounded in [−0.5, 0.5] and should be interpreted relative to the
per-protein null mean and standard deviation.

A methodological note on T₂: at the per-protein level the AUROC-only
Top50 set diverges from the Composite-Top50 set used in §2.2 of the
Results, and per-protein T₂ values are uninformative in four of the
five proteins. Per-protein evidence for the C2 enrichment is therefore
primarily carried by T₃ and T₃ˢᵗᵈ.

## 3.7 Permutation test (within-protein label shuffle)

We tested the null hypothesis that SCI carries no information about
mutation functionality using a within-protein label-shuffle
permutation
(`p6_scripts/p5_permutation.py`, `within_protein_shuffle`, lines
224–229). For each permutation, the `DMS_score_bin` labels were
permuted independently within each of the five proteins, holding fixed
the SCI feature matrix, protein identity, and per-protein label
balance. The five test statistics were recomputed under each
permuted label vector.

Empirical p-values follow the Phipson and Smyth (2010) convention
p = (1 + #{T_perm extreme}) / (N + 1), which lower-bounds p away from
zero. Two-sided p-values (used for T₁) centre the comparison on the
null mean:
|T_perm − μ_null| ≥ |T_obs − μ_null|, fixing a common bug under
which two-sided p is computed against zero
(`empirical_p`, lines 252–260). One-sided p-values (used for T₂, T₃,
T₃ˢᵗᵈ, T_max) use T_perm ≥ T_obs.

The number of permutations follows an adaptive schedule with
budgets {1,000 / 5,000 / 10,000} and an early-stopping rule at
max-p < 0.001 across the four primary statistics
(`run_permutation_block`, lines 354–405). Per-protein analyses used a
fixed budget of N = 1,000 permutations.

Effect sizes are reported as the z-score z = (T_obs − μ_null) /
σ_null and as the raw deviation T_obs − μ_null; under the permutation
null the z-score equals Cohen's *d* and is reported under both labels
in the publication summary. Null 95 % and 99 % percentile intervals are
reported alongside (`summarise`, lines 263–284 of `p5_permutation.py`).

All RNG draws use `numpy.random.default_rng(42)` for the permutation
stream and `default_rng(43)` for the bootstrap stream
(`RNG_SEED = 42`; lines 715–716 of `p5_permutation.py`), ensuring full
reproducibility.

## 3.8 Multiple-testing correction

Two independent corrections were applied to the per-protein primary
tests (4 statistics × 5 proteins = 20 tests;
`N_BONFERRONI_TESTS = 20`):

- **Bonferroni FWER** with p_bonferroni = min(1, p_empirical × 20)
  (`bonferroni_adjust`, lines 686–687 of `p5_permutation.py`); the
  `p_bonferroni` column of `permutation_per_protein_summary.csv`.
- **Benjamini–Hochberg FDR** with step-up monotonicity enforcement
  (`bh_fdr`, lines 690–714 of `p5_permutation.py`); the `p_fdr`
  column of `permutation_per_protein_summary.csv`. The same step-up
  procedure is applied independently to the 4 global primary tests in
  `permutation_global_summary.csv`. NaN p-values are passed through
  as NaN and excluded from the BH family size.

T_max is reported outside the family-wise correction as a
supplementary control statistic for the 528-pair search; its
p_empirical and z-score are still reported, but it does not contribute
to the Bonferroni or BH families.

## 3.9 Non-parametric bootstrap

A non-parametric bootstrap with N = 1,000 resamples and replacement
was used to compute 95 % and 99 % confidence intervals for each of
the five test statistics
(`bootstrap_summary.csv`). For T₁ the bootstrap recomputes the SCI
mean difference directly on the resampled (mutation, label) pairs
(`bootstrap_T1`, lines 411–418 of `p5_permutation.py`). For T₂, T₃,
T₃ˢᵗᵈ and T_max we use a rank-reuse fast path
(`bootstrap_pair_stats`, lines 420–518): instead of recomputing
per-column SCI ranks on every bootstrap resample (the natural
non-parametric bootstrap, with complexity O(B · K · N log N)), the
precomputed full-sample ranks are reused under the bootstrap index
vector. This is a U-statistic approximation whose deviation from the
natural bootstrap AUROC was measured at real scale (N = 95,142,
K = 528, three sanity-check resamples in
`_sanity_check_rank_reuse_bootstrap`): the peak per-pair |ΔAUROC|
was 6 × 10⁻³ and the mean was 1.6 × 10⁻³, both well below the
precision at which 95 % CIs are reported. Measured speedup over the
natural bootstrap on real-scale data was ≈ 110×, reducing 1,000
bootstrap resamples from ≈ 320 min to ≈ 2.9 min.

## 3.10 Statistical reporting

All p-values, z-scores, Cohen's *d* values, raw effect sizes, null 95 %
and 99 % percentile intervals, Bonferroni-adjusted p, BH-adjusted q
and bootstrap 95 % and 99 % CIs are reported jointly in three CSV
artefacts produced by the same pipeline run:

- `data/p0_output/p5_permutation/permutation_global_summary.csv`
  — global tests (5 statistics × 1 scope).
- `data/p0_output/p5_permutation/permutation_per_protein_summary.csv`
  — per-protein tests (5 statistics × 5 proteins = 25 rows).
- `data/p0_output/p5_permutation/publication_summary.csv`
  — formatted manuscript-ready table combining the above with
  bootstrap CIs.

Significance calls in the manuscript use q < 0.05 after BH FDR
correction; the corresponding raw `p_empirical` and `p_bonferroni`
columns are listed in Supplementary Table S2.

## 3.11 Code and data availability

All analysis scripts are stored under `p6_scripts/` of the project
repository:
`config.py`, `extract_site.py`, `recover_wt_v2.py`,
`expand_multi_mutations.py`, `merge_data_v2.py`, `P1_build_sci.py`,
`P1.5_validate_sci_signal.py`, `p0_layer_pair_mining_v2.py`,
`p01_all_layer_pair_stat.py`, `layer_pair_category_distribution.py`,
`protein_clustering_layer_pair.py`, `pca_loadings_extract.py`,
`top_recurrent_cross_protein_pairs.py`,
`layer_recurrence_top20.py`, `layer_occurrence_frequency.py`,
`p5_permutation.py`. All intermediate matrices, summary CSVs,
permutation null distributions
(`permutation_null_distributions.npz`) and figures are stored under
`data/processed/` and `data/p0_output/`. The full `p5_permutation.py`
pipeline runs in ≈ 8–10 min single-threaded on a CPU at the
real-scale (95,142 mutations, 528 layer pairs, 1,000 adaptive
permutations + 1,000 bootstrap resamples), with peak resident memory
≈ 410 MB.

## 3.12 ESM-2 zero-shot LLR baseline

For the head-to-head comparison reported in §2.8 we computed the
ESM-2 zero-shot log-likelihood ratio (LLR) baseline of
Meier et al. (2021). Two standard variants are implemented in
`p6_scripts/p7_compute_esm2_llr.py` on the same
`esm2_t33_650M_UR50D` checkpoint and the same WT sequences
(`data/processed/wt_sequences.csv`) used for SCI feature extraction
(§3.2):

- **WT-marginal.** A single forward pass on the WT sequence; the
  per-position log-softmax of the LM-head logits gives
  log p(aa | x_WT) for each of the 20 standard amino acids at every
  residue position. LLR^{WT} = log p(mut | x_WT) − log p(wt | x_WT)
  at the mutation position; total cost is one forward pass per
  protein (5 forwards).
- **Masked-marginal.** For each unique (protein, `mutation_position`)
  pair, a single forward pass is performed with the token at the
  mutated position replaced by the `<mask>` token; the per-residue
  log-softmax of the LM-head logits at the masked position gives
  log p(aa | x_\T). LLR^{MM} = log p(mut | x_\T) − log p(wt | x_\T).
  Across the five proteins this requires one forward pass per unique
  (protein, position), for a total of ≈ 1,300 forwards
  (`p7_compute_esm2_llr.py`, "Phase 2"). The masked-marginal variant
  is the headline LLR scorer in §2.8 and corresponds to the
  ProteinGym ESM-2 zero-shot baseline.

For multi-mutation rows (the 73,042 expanded PABP doubles in
`all_proteins_site_metadata.csv`), each constituent single substitution
is scored independently with single-position masking on the WT
context, following the additive ProteinGym leaderboard convention for
ESM-2 masked-marginal scoring of multi-residue substitutions. Per-row
LLR is written to `data/processed/esm2_llr_scores.csv` together with
the WT-marginal score, the `DMS_score`, `DMS_score_bin` and the
`(wt_aa, mut_aa, mutation_position)` triple parsed from the `mutant`
string (`p7_compute_esm2_llr.py`, "Phase 3"). The implementation
re-uses the existing cached model weights at
`model/hub/checkpoints/esm2_t33_650M_UR50D.pt` and adds no extraction
of the 33-layer hidden states already cached for SCI; the only
additional forward work is the ≈ 1,305 forwards described above
(≈ 70× fewer than the 95,142 forwards used for the SCI Δ-embedding
cache in §3.2).

For the SCI-vs-LLR comparison reported in §2.8, the SCI scorer used is
the 528-pair upper-triangle mean
(`all_proteins_sci_site_scores_mean.npy`), which matches the
per-mutation SCI summary characterised in `P1.5_signal_summary.csv`.
This differs from the Top50 SCI summary used by T₁ in §2.1: T₁ is
a mean-difference test designed for discrimination against a
within-protein permutation null and uses the high-signal tail; the
present comparison uses the all-pairs mean as the simplest scalar SCI
summary, for parity with a scalar zero-shot baseline. The combined
SCI+LLR scorer is a 5-fold stratified cross-validated logistic
regression on the standardised (SCI, LLR) feature pair, implemented in
`p6_scripts/p7_sci_vs_llr_comparison.py` with `StratifiedKFold(
n_splits = 5, shuffle = True, random_state = 42)`, a `StandardScaler`
fit *within* each training fold to prevent leakage, and the held-out
out-of-fold P(class = 1) used as the predictor for each held-out fold.
Per-protein and pooled Spearman, AUROC and AUPRC values, computed with
the same sign-correction convention used in
`P1.5_validate_sci_signal.py` (corrected = max over native and
negated score), are written to
`data/processed/sci_vs_llr_comparison.csv`; the per-protein and
pooled LLR-only summary used as a sanity check is in
`data/processed/baseline_summary.csv`.
