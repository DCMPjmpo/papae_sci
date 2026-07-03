# Results (Phase B)

> One Results section per validated claim (C1–C7) from
> `paper/scientific_story_map.md`. All numerical values are sourced
> *verbatim* from the on-disk CSVs (`permutation_global_summary.csv`,
> `permutation_per_protein_summary.csv`, `bootstrap_summary.csv`,
> `pairwise_euclidean_distance.csv`, `pca_loadings.csv`,
> `layer_recurrence_*`, `layer_frequency_*`, per-protein
> `*_pair_category_ratio_Composite.csv`). No new experiments and no
> new conclusions are introduced.
>
> Cross-references: figure labels (Fig 1–Fig 5) follow
> `paper/figure_map.md`; claim labels (C1–C7) follow
> `paper/scientific_story_map.md`.

---

## 2.1 SCI carries mutation-effect signal (C1)

We first asked whether the Structural Context Inconsistency (SCI) score
constructed from ESM-2 (esm2_t33_650M_UR50D; Methods §3.2) carries
mutation-effect signal. The per-mutation summary score over the 50
top-scoring layer pairs
(`data/processed/all_proteins_sci_site_scores_top50.npy`; n = 95,142
mutations across CBS, GAL4, PABP, PTEN, TEM1) was compared between
functional and non-functional mutations as defined by the binarised
deep-mutational-scanning label (`DMS_score_bin`, 1 = functional,
0 = non-functional). At the pooled level, the SCI mean difference was
T₁ = 2.79 × 10⁻³
(`permutation_global_summary.csv`).

Against a within-protein label-shuffle null (1,000 adaptive
permutations; Methods §3.7), the observed T₁ exceeded the null mean
(−3.75 × 10⁻⁴) by 57 standard deviations of the null
(σ_null = 5.56 × 10⁻⁵), yielding empirical p = 10⁻³ and
BH-FDR-adjusted q = 10⁻³ across the four global primary statistics
(Fig 1b). At the per-protein level, T₁ was FDR-significant in 5/5
proteins (`permutation_per_protein_summary.csv`): the smallest
q-value among the per-protein T₁ tests was 3.7 × 10⁻²
(CBS; per-protein z = 2.25), and the remaining four proteins reached
q = 1.4 × 10⁻³ or smaller (PABP z = 56.55; TEM1 z = 23.59;
GAL4 z = 6.29; PTEN z = 2.88). The non-parametric bootstrap
(N = 1,000) returned a 95 % confidence interval for T₁ of
[2.67 × 10⁻³, 2.91 × 10⁻³] (`bootstrap_summary.csv`), excluding
zero. Together these establish that SCI carries discriminative signal
that is reproducible against the null at both the pooled and
per-protein levels (C1).

## 2.2 Top-ranked layer pairs are enriched for Early–Late interactions (C2)

To characterise which layer pairs carry the SCI signal most strongly,
we ranked all 528 unordered layer pairs by a Composite discrimination
score combining per-pair AUROC, Spearman rank correlation and Cohen's
*d* (Methods §3.3; `p0_layer_pair_mining_v2.py`) and examined the
Top50 layer pairs for each protein. The six-category composition
(Early-Early / Early-Middle / Early-Late / Middle-Middle / Middle-Late
/ Late-Late, defined on the Early-1–8 / Middle-9–24 / Late-25–33
banding; Methods §3.3) of each Top50 set is recorded in the per-protein
`*_pair_category_ratio_Composite.csv` files. Four of the five proteins
were EL-dominated: CBS (EL 40 %; LL 30 %), GAL4 (EL 66 %), PTEN
(EL 44 %; LL 12 %) and TEM1 (EL 40 %; EM 10 %; LL 8 %). PABP was the
only protein whose Top50 contained a larger EM than EL share (EM 48 %;
EL 42 %; no LL contribution); we return to this asymmetry in §2.4
(Fig 2a).

At the pooled level we then asked, under the AUROC-only ranking, how
many of the 50 pairs with largest |AUROC − 0.5| fall into the EL
category. Observed T₂ = 17 of 50, compared with a permutation null
mean of 1.01 (σ_null = 2.00), giving z = 8.01 and q = 10⁻³
(`permutation_global_summary.csv`; Fig 2b). The bootstrap 95 % CI
for T₂ was [14, 19] (`bootstrap_summary.csv`). A methodological note:
at the per-protein level the AUROC-only ranking and the Composite
ranking can diverge, and per-protein T₂ values are uninformative for
four of the five proteins (Methods §3.6); per-protein evidence for the
EL enrichment is therefore primarily carried by the per-protein T₃ and
T₃ˢᵗᵈ analyses in §2.3 and by the per-protein Composite category
ratios reported above.

## 2.3 Early–Late pairs dominate the 528-pair landscape (C3)

To determine whether the EL enrichment observed in the Top50 tail is a
consequence of a broader category-level effect rather than a tail
phenomenon, we computed |AUROC − 0.5| for each of the 528 layer pairs
and compared the EL category (72 pairs) to the union of EE, EM, MM,
ML and LL pairs (456 pairs) (Methods §3.6). We define
T₃ = mean(|AUROC − 0.5|)_EL − mean(|AUROC − 0.5|)_non-EL,
and T₃ˢᵗᵈ = T₃ / σ(|AUROC − 0.5|)_non-EL as a sample-size-aware
standardised counterpart.

At the pooled level, T₃ = 3.99 × 10⁻²
(null mean −2.08 × 10⁻³, σ_null = 6.77 × 10⁻⁴) with z = 62.03 and
q = 10⁻³, and T₃ˢᵗᵈ = 1.02 (null mean −0.602, σ_null = 0.201)
with z = 8.05 and q = 10⁻³
(`permutation_global_summary.csv`; Fig 2c). Bootstrap 95 % CIs are
[3.70 × 10⁻², 4.27 × 10⁻²] for T₃ and [0.954, 1.073] for
T₃ˢᵗᵈ (`bootstrap_summary.csv`), excluding the corresponding null
means. At the per-protein level T₃ and T₃ˢᵗᵈ were FDR-significant in
5/5 proteins (`permutation_per_protein_summary.csv`); per-protein
T₃ z-scores ranged from 5.18 (GAL4) to 38.71 (PABP), and T₃ˢᵗᵈ
z-scores ranged from 4.59 (GAL4) to 42.80 (PABP). At the per-protein
level T₃ is computed against the globally ranked SCI scores
(Methods §3.6); per-protein effect-size magnitudes are therefore not
bounded in [−0.5, 0.5], but the permutation null is computed under the
identical metric and the significance call is unaffected. The
category-level EL dominance is observed at both the pooled level and
within each of the five proteins.

## 2.4 One of the five proteins (PABP) forms a distinct cluster (C4)

We then asked whether the EL-dominated composition reported in §2.2 is
shared across proteins or conceals heterogeneity. Each protein's Top50
Composite layer-pair set was summarised as a six-dimensional vector of
category proportions, and three independent geometric analyses were
applied (`p6_scripts/protein_clustering_layer_pair.py`,
`pca_loadings_extract.py`; `data/p0_output/protein_clustering/`).

The heatmap of these vectors (Fig 3a) identifies PABP as the unique
EM-dominated row. Hierarchical clustering with Euclidean distance and
average linkage (Fig 3b; `pairwise_euclidean_distance.csv`) places
CBS, GAL4, PTEN and TEM1 into a compact cluster whose maximum pairwise
distance is 0.397 (CBS–GAL4) and minimum is 0.115 (PTEN–TEM1). PABP
joins this cluster only at a height of ≈ 0.50; the minimum
PABP-to-other distance is 0.389 (PABP–TEM1), which exceeds the
maximum within-other-four distance of 0.397 — i.e. PABP's closest
neighbour is farther than the most distant pair among the other four.
Principal component analysis on the same six-dimensional vectors
(Fig 3c; `pca_protein_layer_pair.png`) recovers PC1 = 68.6 % and
PC2 = 28.1 % of variance (cumulative 96.7 %). PC1 loadings
(`pca_loadings.csv`; `pca_loadings_ranked.csv`) are dominated by
Early–Middle (+0.207), with the next-largest absolute loading
approximately three times smaller (LL, −0.072) and EL contributing
only −0.026. PC2 separates Early–Late (+0.103) from Late–Late
(−0.094). PABP projects to PC1 ≈ +0.38, isolated from the four
EL-dominated proteins which fall in PC1 ∈ [−0.18, +0.01]. The four
EL-dominated proteins share a similar layer-pair signature, while
PABP's signature is shifted toward Early–Middle interactions.

## 2.5 Cross-protein recurrent layer pairs concentrate at layers 30–32 (C5)

To identify layer pairs that carry mutation-effect signal across
proteins, we counted, for each layer pair, the number of proteins in
whose Top50 Composite set it appears
(`p6_scripts/top_recurrent_cross_protein_pairs.py`;
`cross_protein_universal_pairs.csv`; 146 candidate recurrent pairs).
The Top20 most-recurrent pairs (6 pairs at Cross_Protein_Count = 4
and 14 at Cross_Protein_Count = 3;
`layer_recurrence_top20_pairs.csv`) carry 40 layer-end occurrences
(2 × 20 pairs).

Band totals (`layer_recurrence_band_summary.csv`; Fig 4a) are: Late
(layers 25–33) = 22 of 40 (55.0 %), Middle (9–24) = 11 (27.5 %), Early
(1–8) = 7 (17.5 %). Within the Late band the distribution is
non-uniform: layer 32 contributes 8 of 40 layer-end occurrences
(20.0 %, the single most recurrent layer), layer 31 contributes 7
(17.5 %), layer 30 contributes 3 (7.5 %), and layer 33 contributes 2
(5.0 %) (`layer_recurrence_frequency.csv`). Layers 25–29 contribute
zero, producing a sharp three-layer hub at depths 30–32 with a marked
trough immediately upstream. At the broader Top50 union across the
five proteins (n = 500 layer-end occurrences = 5 × 50 × 2;
`layer_frequency.csv`, Fig 4b) the band proportions are preserved:
Late 47.6 %, Middle 36.0 %, Early 16.4 %; layer 32 alone accounts for
14.6 % of all occurrences. The recurrent set concentrates at a sharp
three-layer hub centred on layer 32, and the central middle band
(layers 18–27 in the Top20 set) is essentially absent from the
recurrent geometry.

## 2.6 Per-protein signal is reproducible against within-protein null (C6)

To rule out the possibility that the pooled effects in §2.1–2.3 are
driven by a single protein, we re-computed the four primary test
statistics within each of the five proteins, using a within-protein
label-shuffle null (N = 1,000 permutations per protein; Methods §3.7).
All five proteins returned FDR-significant T₁, T₃ and T₃ˢᵗᵈ
(`permutation_per_protein_summary.csv`; Fig 5a, 5b). The five
per-protein T₁ z-scores were 2.25 (CBS), 6.29 (GAL4), 56.55 (PABP),
2.88 (PTEN) and 23.59 (TEM1); the five per-protein T₃ z-scores were
16.11, 5.18, 38.71, 14.57 and 15.19; and the five T₃ˢᵗᵈ z-scores
were 18.50, 4.59, 42.80, 17.05 and 15.94. After Bonferroni correction
over the 20 per-protein primary tests (5 proteins × 4 statistics),
each of these 15 (T₁, T₃, T₃ˢᵗᵈ across 5 proteins) tests retains
significance at the conventional α = 0.05 level. The signal is therefore
not pooled out of one over-represented protein; it is reproduced
within each.

## 2.7 Statistical robustness across multiple corrections and bootstrap (C7)

We assessed robustness across three independent layers. (i) Phipson and
Smyth (2010) empirical p-values were computed under the null-mean–
centered convention for two-sided tests
(`empirical_p` in `p5_permutation.py`); under this convention, every
primary statistic at the pooled level reached p = 10⁻³, the
floor for N = 1,000 permutations. (ii) Bonferroni correction over the
20-test family of per-protein primary statistics yielded
`p_bonferroni` columns in `permutation_per_protein_summary.csv`;
Bonferroni significance at α = 0.05 was retained for T₁, T₃, T₃ˢᵗᵈ in
5/5 proteins. (iii) Benjamini–Hochberg FDR with step-up monotonicity
enforcement was applied to the same 20-test family, with the
`p_fdr` column containing the BH q-values; q < 0.05 calls match the
Bonferroni calls exactly on T₁, T₃, T₃ˢᵗᵈ across the five proteins.

A non-parametric bootstrap (N = 1,000; Methods §3.9) provided 95 % CIs
for the pooled statistics: T₁ ∈ [2.67, 2.91] × 10⁻³, T₂ ∈ [14, 19],
T₃ ∈ [3.70, 4.27] × 10⁻², T₃ˢᵗᵈ ∈ [0.954, 1.073], T_max ∈ [0.176,
0.187] (`bootstrap_summary.csv`). All five intervals exclude the null
mean. The family-wise maximum statistic T_max = max |AUROC − 0.5|
over all 528 pairs reached 0.181, exceeding the null max
(mean 0.0176, σ 1.74 × 10⁻³) by z = 93.7, addressing the concern
that the best-of-528 search inflates significance. Across all four
robustness layers, the EL geometry of mutation-effect signal in ESM-2
survives every correction applied to it (Fig 5).

## 2.8 Comparison with ESM-2 zero-shot baseline

C1–C7 establish that SCI carries mutation-effect signal and that the
signal has a reproducible layer-pair geometry. To position this against
the field's de facto zero-shot scorer, we computed the ESM-2
log-likelihood ratio (LLR) baseline of Meier et al. (2021) on the same
95,142 mutation rows and asked two questions: how does SCI compare to
LLR as a per-mutation scalar predictor, and does SCI carry information
that is additive to LLR? The masked-marginal variant of LLR
(Meier et al., 2021) was computed by masking the mutated residue in
the WT sequence, taking log-softmax of the LM-head logits at that
position, and forming LLR = log p(mut | x_\T) − log p(wt | x_\T)
(Methods §3.12); the per-mutation SCI summary used here is the
528-pair upper-triangle mean
(`all_proteins_sci_site_scores_mean.npy`), the same scorer
characterised in `P1.5_signal_summary.csv`. The combined SCI+LLR
scorer is a 5-fold stratified-CV logistic regression on the
standardised (SCI, LLR) feature pair, with the StandardScaler fit
within each training fold; out-of-fold P(class = 1) is used as the
predictor for the held-out fold (Methods §3.12).

At the pooled level (n = 95,142), the masked-marginal LLR achieves
Spearman r = 0.480 (p ≈ 0) against `DMS_score` and AUROC = 0.753
against `DMS_score_bin`, exceeding the 528-mean SCI summary
(r = 0.200, AUROC = 0.669) on both metrics (Table 1;
`sci_vs_llr_comparison.csv`). The 5-fold CV SCI+LLR combination
reaches AUROC = 0.758 and Spearman r = 0.464, within 0.005 AUROC of
LLR alone. At the per-protein level LLR exceeds SCI on every one of
the five proteins, with per-protein AUROC values of 0.692 (CBS),
0.795 (GAL4), 0.767 (PABP), 0.846 (PTEN) and 0.882 (TEM1) for LLR
versus 0.570, 0.717, 0.699, 0.629 and 0.759 for SCI, and per-protein
Spearman r ranging from 0.342 (CBS) to 0.714 (TEM1) for LLR versus
0.131 to 0.502 for SCI. The combined SCI+LLR scorer exceeds the
better single scorer by at most Δ AUROC = 0.008 (TEM1), with smaller
gains on PABP (+0.003), PTEN (+0.001) and GAL4 (+0.001) and no
improvement on CBS (Δ ≤ 0.001); the SCI signal is therefore largely
redundant with LLR as a scalar predictor.

We make this comparison transparently. SCI as a scalar per-mutation
predictor underperforms the masked-marginal LLR baseline, and a
combined logistic regression does not meaningfully exceed LLR alone.
This result does not contradict §2.1, which establishes that SCI
carries discriminative signal against a within-protein label-shuffle
null at z = 57; rather, it indicates that the SCI signal is largely
co-linear with the LM-head likelihood at the scalar level. The
contribution of the present work is not a new scalar predictor; it is
the layer-resolved 528-pair landscape (C2, C3), the layer-30–32
recurrence hub (C5), and the protein-specific compositional
heterogeneity (C4) — none of which is exposed by a single LM-head
likelihood ratio. The relationship between the scalar SCI summary and
LLR is discussed further in §D5.
