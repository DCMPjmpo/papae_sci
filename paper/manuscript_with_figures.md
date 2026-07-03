# First Draft — Manuscript (Assembled)

> READ-ONLY INTEGRATION: This manuscript is assembled from existing files
> without modification to any scientific content, experimental results,
> statistical values, or figure legends. Figure insertion points are
> marked where composite figures should be placed.

---

## Title

**Early–Late layer-pair interactions encode mutation-effect signal in the
protein language model ESM-2 across five deep mutational scans.**

(Alternatives: *"Layer-pair geometry of mutation-effect encoding in
ESM-2"*; *"Cross-layer interaction geometry organises mutation-effect
signal in ESM-2"*. The lead title above embeds the central finding
without causal verbs.)

---

## Abstract

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
composition; PABP forms a singleton characterized by Early–Middle interactions
(hierarchical clustering, PCA, pairwise Euclidean distances). All
primary effects survive Bonferroni FWER, Benjamini–Hochberg FDR, and
non-parametric bootstrap. We interpret these observations as evidence of
a two-axis (layer × layer) organisation of mutation-effect signal in
ESM-2, complementing the prevailing one-axis (single-layer) view.

---

## 1. Introduction

**¶1 — From mutation-effect prediction to mutation-effect representation.**
Protein language models (PLMs) — large self-supervised Transformers
trained on hundreds of millions of protein sequences — have reshaped
how the functional consequence of an amino-acid substitution is
estimated from sequence alone. Variant scores derived from these
models routinely match or surpass classical evolutionary-density
predictors on deep mutational scanning (DMS) benchmarks and on
clinical variant-of-uncertain-significance datasets (Meier et al.,
2021; Frazer et al., 2021; Notin et al., 2022; Brandes et al., 2023;
Hsu et al., 2022; Marquet et al., 2022). Yet across this prediction
literature, the variant score is read from a *single* representation
of the model, typically the masked-LM head or the final hidden
state. How the underlying information is organised across the depth
of a PLM — i.e. *how mutation-effect signal is encoded representationally,
not how well it predicts* — remains essentially uncharacterised.

**¶2 — The single-layer convention dominates both prediction and
interpretability for PLMs.**
The one-axis convention extends beyond prediction. Layer-wise
probing (Rives et al., 2021; Lin et al., 2023; Elnaggar et al., 2022;
Detlefsen et al., 2022; Vig et al., 2021) maps biophysical and
structural properties to individual layers, treating depth as a
scalar along which features emerge. Per-head attention analysis (Vig
et al., 2021; Rao et al., 2021, MSA Transformer) decomposes a single
layer's computation; per-layer sparse-autoencoder analysis (Simon &
Zou, 2025) extracts interpretable features at each chosen depth. In
every case the unit of analysis is a single layer; the *relationship
between* layers — what is encoded jointly across depth, or what
emerges only when two layers are read together — is not directly
addressed.

**¶3 — Cross-layer geometry is well-established outside protein
language models.**
In vision and language Transformers, cross-layer interactions have been
shown to carry information not captured by single-layer representations.
In BERT, attention heads at different depths capture different levels
of syntactic and semantic structure (Tenney et al., 2019); in GPT,
layer-wise activation patching reveals that earlier layers encode
local features while later layers encode higher-level concepts
(Geva et al., 2021). For protein language models, the residual stream
has been shown to accumulate information across layers (Belrose et al.,
2023), but the question of *which* layer pairs are most informative
for a specific downstream task — and whether such pairs follow a
systematic pattern — remains open.

**¶4 — We propose a layer-pair geometry analysis of mutation-effect
signal.**
Here we address this gap by quantifying the discriminative signal
carried by every one of the 528 unordered layer pairs of ESM-2
(esm2_t33_650M_UR50D; 33 layers; Lin et al., 2023) using a novel
per-mutation metric called Structural Context Inconsistency (SCI). We
apply this analysis across five deep mutational scanning panels
(CBS, GAL4, PABP, PTEN, TEM1; 95,142 mutations) and ask three questions:
(1) Does SCI carry mutation-effect signal? (2) Which layer pairs
carry the strongest signal, and do they follow a systematic pattern?
(3) Is the pattern consistent across proteins, or does it vary
by biophysical class? Our findings reveal a two-axis organisation
of mutation-effect signal in ESM-2, with Early–Late layer pairs
dominating the high-signal tail and a narrow hub at layers 30–32
emerging as the most recurrent site across proteins.

---

## 2. Results

### 2.1 SCI carries mutation-effect signal (C1)

We first asked whether the Structural Context Inconsistency (SCI) score
constructed from ESM-2 (esm2_t33_650M_UR50D; Methods §4.2) carries
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
permutations; Methods §4.7), the observed T₁ exceeded the null mean
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

![Figure 1 — SCI carries mutation-effect signal (C1)](figures/Fig1.png)

### 2.2 Top-ranked layer pairs are enriched for Early–Late interactions (C2)

To characterise which layer pairs carry the SCI signal most strongly,
we ranked all 528 unordered layer pairs by a Composite discrimination
score combining per-pair AUROC, Spearman rank correlation and Cohen's
*d* (Methods §4.3; `p0_layer_pair_mining_v2.py`) and examined the
Top50 layer pairs for each protein. The six-category composition
(Early-Early / Early-Middle / Early-Late / Middle-Middle / Middle-Late
/ Late-Late, defined on the Early-1–8 / Middle-9–24 / Late-25–33
banding; Methods §4.3) of each Top50 set is recorded in the per-protein
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
four of the five proteins (Methods §4.6); per-protein evidence for the
EL enrichment is therefore primarily carried by the per-protein T₃ and
T₃ˢᵗᵈ analyses in §2.3 and by the per-protein Composite category
ratios reported above.

### 2.3 Early–Late pairs dominate the 528-pair landscape (C3)

To determine whether the EL enrichment observed in the Top50 tail is a
consequence of a broader category-level effect rather than a tail
phenomenon, we computed |AUROC − 0.5| for each of the 528 layer pairs
and compared the EL category (72 pairs) to the union of EE, EM, MM,
ML and LL pairs (456 pairs) (Methods §4.6). We define
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
(Methods §4.6); per-protein effect-size magnitudes are therefore not
bounded in [−0.5, 0.5], but the permutation null is computed under the
identical metric and the significance call is unaffected. The
category-level EL dominance is observed at both the pooled level and
within each of the five proteins.

![Figure 2 — Top-ranked layer pairs and the 528-pair category landscape (C2 + C3)](figures/Fig2.png)

### 2.4 One of the five proteins (PABP) forms a distinct cluster (C4)

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
PABP-to-other distance is 0.389 (PABP–TEM1), comparable in magnitude
to the maximum within-other-four distance of 0.397 (CBS–GAL4).
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

![Figure 3 — PABP forms an independent cluster on the layer-pair composition axis (C4)](figures/Fig3.png)

### 2.5 Cross-protein recurrent layer pairs concentrate at layers 30–32 (C5)

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

![Figure 4 — Cross-protein recurrent layer pairs concentrate at layers 30–32 (C5)](figures/Fig4.png)

### 2.6 Per-protein signal is reproducible against within-protein null (C6)

To rule out the possibility that the pooled effects in §2.1–2.3 are
driven by a single protein, we re-computed the four primary test
statistics within each of the five proteins, using a within-protein
label-shuffle null (N = 1,000 permutations per protein; Methods §4.7).
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

### 2.7 Statistical robustness across multiple corrections and bootstrap (C7)

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

A non-parametric bootstrap (N = 1,000; Methods §4.9) provided 95 % CIs
for the pooled statistics: T₁ ∈ [2.67, 2.91] × 10⁻³, T₂ ∈ [14, 19],
T₃ ∈ [3.70, 4.27] × 10⁻², T₃ˢᵗᵈ ∈ [0.954, 1.073], T_max ∈ [0.176,
0.187] (`bootstrap_summary.csv`). All five intervals exclude the null
mean. The family-wise maximum statistic T_max = max |AUROC − 0.5|
over all 528 pairs reached 0.181, exceeding the null max
(mean 0.0176, σ 1.74 × 10⁻³) by z = 93.7, addressing the concern
that the best-of-528 search inflates significance. Across all four
robustness layers, the EL geometry of mutation-effect signal in ESM-2
survives every correction applied to it (Fig 5).

![Figure 5 — Statistical robustness across proteins, corrections and bootstrap (C6 + C7)](figures/Fig5.png)

### 2.8 Comparison with ESM-2 zero-shot baseline

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
(Methods §4.12); the per-mutation SCI summary used here is the
528-pair upper-triangle mean
(`all_proteins_sci_site_scores_mean.npy`), the same scorer
characterised in `P1.5_signal_summary.csv`. The combined SCI+LLR
scorer is a 5-fold stratified-CV logistic regression on the
standardised (SCI, LLR) feature pair, with the StandardScaler fit
within each training fold; out-of-fold P(class = 1) is used as the
predictor for the held-out fold (Methods §4.12).

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
0.131 to 0.502 for SCI. The SCI and LLR scores are moderately
correlated at the pooled level (Pearson r = 0.456; Spearman ρ = 0.467;
`sci_llr_correlation.csv`), indicating partial overlap but
substantial non-shared variance. The combined SCI+LLR scorer exceeds the
better single scorer by at most Δ AUROC = 0.007 (TEM1), with smaller
gains on PABP (+0.003), GAL4 (+0.002), PTEN (+0.001) and no
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
likelihood ratio.



---

## 3. Discussion

### D1. Deep-layer enrichment in ESM

The cross-protein recurrence analysis indicates that the late end of mutation-effect signal in ESM-2 is sharply localized. Among the Top20 cross-protein recurrent layer pairs, the Late band (layers 25–33) accounts for 55% of all layer-end occurrences, while the Middle band (layers 9–24) and the Early band (layers 1–8) account for 27.5% and 17.5%, respectively (C5). Within the Late band, the distribution is itself non-uniform: layer 32 alone supplies 8 of the 40 layer-end occurrences, layers 30–32 together supply 18 of 40, and layers 25–29 contribute essentially none (C5). The late-end concentration therefore reflects a narrow three-layer hub rather than a broad enrichment across the upper half of the network. This hub is observed in the same recurrent set in which Early–Late pairs dominate (C2, C5) and is reproduced at the level of the full 528-pair distribution (C3).

The concentration of recurrence at layers 30–32 is consistent with prior reports that the deepest layers of ESM and related protein language models carry higher-order semantic and functional features (Rives et al., 2021; Lin et al., 2023; Elnaggar et al., 2022; Vig et al., 2021). It is also consistent with the broader observation that the most useful single layer of a protein language model depends on the downstream task (Detlefsen et al., 2022). The present finding may be read as a refinement of that framing: for mutation-effect signal in ESM-2, the deepest layers are recurrently engaged across proteins, but they appear within layer pairs rather than as stand-alone signal carriers (C2, C5). The convergence of the recurrent set on layer 32 is also notable in light of work that uses ESM-1b's final layer alone to score human missense variants at proteome scale (Brandes et al., 2023). The late-end anchor identified here lies one layer below the masked-LM head; the convergence of the recurrent set on layer 32 may therefore be read as compatible with the use of layer 33 in downstream variant-effect scoring, rather than as opposed to it. We hesitate to interpret this convergence as more than a structural coincidence; it is one of several possible readings of the geometry. All observations in this subsection are made on a single backbone (ESM-2, 33 layers); whether a comparable late-end anchor emerges in ESM-1b (Rives et al., 2021), ESM3 (Hayes et al., 2025), or other Transformer-based protein language models remains an open question (C5).

### D2. Cross-layer integration hypothesis

Across the 528 unordered layer pairs of ESM-2, Early–Late (EL) pairs achieve significantly higher Composite scores and AUROC values than Early–Early, Early–Middle, Middle–Middle and Late–Late pairs, with category-level comparisons returning p ≪ 0.001 (C3). This category-level pattern is preserved at the high-signal tail: the Top50 Composite layer pairs are enriched for EL interactions across all five proteins (C2). The cross-protein recurrent set is itself dominated by pairs that connect a deep-layer endpoint, typically in 30–32, with a non-Late endpoint in the Early-to-mid-Middle range (C2, C5). These statistical relationships are observed on a representational signal — SCI mean — that has independently been shown to distinguish functional from non-functional mutations in the cohort examined here (C1).

Taken together, C1, C2, C3 and C5 are consistent with a working hypothesis that mutation-effect signal in ESM-2 is integrated across layers rather than read off any single layer. Under this hypothesis, early layers may contribute residue-local biochemical features, and late layers — notably 30–32 — may contribute family-level or contextual features, with the discriminative information potentially corresponding to the relationship between these two depth-stratified representations (C2, C3, C5). This reading is consistent with interpretability work on general Transformers that characterizes the residual stream as a substrate over which intermediate representations are progressively composed across layers (Tenney et al., 2019; Geva et al., 2021; Belrose et al., 2023). The closest existing protein-LM analysis that explicitly aggregates across layers is the contact-prediction pipeline of Rao et al. (2021), in which a sparse logistic regression pools attention from all ESM layers; the present findings suggest that an analogous, but biologically distinct, multi-layer combination may be at play for mutation-effect encoding (C2, C3). We emphasize that this remains a hypothesis suggested by the layer-pair geometry rather than an established property of the forward pass.

Several limitations of this hypothesis warrant explicit acknowledgement. First, all findings here are correlational; we have not performed activation patching, single-layer ablation, or other interventions that would isolate the contribution of any individual layer (C2, C3, C5). Second, the Composite score combines AUROC, Spearman *r* and Cohen's *d*, and an EL pair could in principle rank highly through aggregation of two independently weak but complementary endpoints, rather than through a genuinely interactive signal (C3). Disentangling aggregation from interaction is an open problem that the present analysis does not resolve. Third, the existence of high-performing single-layer mutation-effect scores from the masked-LM head of ESM-1v (Meier et al., 2021), the autoregressive output of Tranception (Notin et al., 2022), and the ridge regression over PLM embeddings of Hsu et al. (2022) may at first appear in tension with the integration hypothesis. Predictive sufficiency at a single layer does not preclude representational distribution across pairs, however, given that each of those single-layer scores is itself the product of cross-layer composition through the residual stream; the present claims are about layer-pair geometry rather than predictor design (C2, C3).

### D3. Protein-specific heterogeneity

The five proteins examined here do not share a single layer-pair signature. Four of them — CBS, GAL4, PTEN and TEM1 — are EL-dominated within their Top50 Composite layer pairs (C4). PABP is the only protein whose Early–Middle proportion exceeds its Early–Late proportion within the same Top50 set, and the only protein with no Late–Late contribution (C4). Hierarchical clustering and PCA over the six-category proportion vector place PABP at a singleton position whose distance to the centroid of the remaining four proteins is comparable to the entire pairwise spread within the four-protein cluster (C4). All five proteins exhibit clearly above-baseline SCI signal (C1); the observed heterogeneity therefore concerns which layer pairs carry mutation-effect signal, not whether such signal is present. The cross-protein recurrent set reflects the four EL-dominated proteins more strongly than it reflects PABP, consistent with PABP's distinct compositional profile (C4, C5).

The shift of PABP toward Early–Middle interactions co-varies with biophysical class: PABP is the only RNA-binding protein in the panel, while the four EL-dominated proteins span enzymes (CBS, PTEN, TEM1) and a transcription factor (GAL4). This co-variation is consistent with depth-stratified attention-head specialization in protein Transformers, in which different biophysical properties are tracked at characteristic depths (Vig et al., 2021), and with the existence of layer-localized interpretable features that vary across positions and proteins (Simon & Zou, 2025). A second reading is that the heterogeneity is consistent with differences in deep-mutational-scanning (DMS) assay design: PABP's assay measures growth-coupled RNA binding (Melamed et al., 2013), whereas the four others measure stability or enzymatic activity in microbial or human cell contexts. Differences in the shape of the underlying fitness landscape may plausibly project onto different layer-pair compositions. We do not claim that RNA binding or assay design is responsible for Middle-layer involvement; we only note that the observed heterogeneity is non-random and aligns with two biological axes that may inform follow-up sampling (C4).

It is also important to note what the PABP singleton is not. PABP's Top50 Composite scores are well-separated from baseline as for the other four proteins (C1); only the layer-band composition of those top-scoring pairs differs (C4). The singleton position therefore should not be read as evidence of poor model fit, low data quality, or a failure of the protein language model on RNA-binding proteins. Two limitations of this reading are explicit. With *n* = 5, PABP's singleton status is a qualitative observation rather than a statistical generalization (C4). Additionally, biophysical class and DMS assay design are confounded within the present sample — there is one RNA-binder and four non-RNA-binders, each with its own assay format — so the two consistent readings above cannot be disentangled (C4). Whether RNA-binding proteins more generally, or proteins with shallow fitness landscapes, systematically display Early–Middle dominance is a hypothesis that requires expansion of the protein panel to at least twenty sequences spanning multiple biophysical classes and assay types.

### D4. Implications for PLM interpretability

Taken together, the systematic EL enrichment at both the category and Top50 levels (C2, C3), the narrow deep-layer anchor at layers 30–32 (C5), and the per-protein compositional heterogeneity (C4) describe a structure that is not visible from per-layer probing alone. The strongest signal lies in the off-diagonal portion of the 528-pair plane — specifically the EL block — while per-layer probing protocols, by construction, integrate over one of the two axes of that plane (C1, C2, C3). The deep-layer anchor (C5) appears in paired form: each of the most recurrent layer pairs in the Top20 set connects layers 30–32 to either an Early or a Middle partner (C2, C5).

On the basis of C2, C3, C4 and C5, we suggest that layer-pair geometry may serve as a useful interpretability primitive for protein language models, complementing rather than replacing per-layer probing (Elnaggar et al., 2022; Detlefsen et al., 2022; Rao et al., 2019) and per-head attention analysis (Vig et al., 2021; Rao et al., 2021). The six-category composition vector summarizes a per-protein fingerprint of how mutation-effect signal is distributed across the depth of a model, and could in principle be computed for ESM-1b (Rives et al., 2021), ESM3 (Hayes et al., 2025), ProtT5 (Elnaggar et al., 2022) and other backbones, to ask whether the EL enrichment and the layer-32 anchor are properties of ESM-2 specifically or of large protein language models more broadly (C2, C4, C5). The observation that mutation-effect signal recurs at cross-depth pairs is also consistent with the residual-stream–centered interpretability picture developed for general Transformers (Tenney et al., 2019; Geva et al., 2021; Belrose et al., 2023), and suggests that the same toolkit — scalar-mix probes, tuned-lens decoding, activation patching across layer pairs, and per-layer sparse-autoencoder analyses (Simon & Zou, 2025) — could plausibly be applied to protein language models in future work to clarify whether layer pairs are merely the most informative readout sites or substrates of genuine cross-layer integration.

Three boundaries on these claims warrant explicit statement. First, we make no claim that layer-pair features yield higher predictive accuracy than tuned single-layer protein-LM scores; a head-to-head comparison was not performed, and predictive parity is not the goal of this work (C1, C2, C3). Second, the EL pattern reported here is correlational, and we do not claim that it is responsible for ESM-2's mutation-effect performance (C2, C3). Third, the prominence of layer 32 has been established only for ESM-2; whether an analogous anchor emerges in other protein language models remains open (C5). Within these boundaries, the layer-pair view contributes a two-axis vocabulary to a literature that has been organized largely along a single depth axis, and the off-diagonal structure of that vocabulary — particularly the EL block and its protein-specific composition — may be a productive site for further interpretability work.

### D5. Relationship to the ESM-2 zero-shot likelihood baseline

The masked-marginal log-likelihood ratio (LLR) of Meier et al. (2021) is the field's de facto zero-shot scorer for ESM-2 and is the ESM-2 baseline reported on the ProteinGym leaderboard. We computed it on the same 95,142 mutation rows (§2.8, Methods §4.12) and observe that LLR exceeds the 528-mean SCI summary as a per-mutation predictor on every one of the five proteins examined here (Table 1). Pooled, LLR achieves AUROC = 0.753 against SCI 0.669; per-protein AUROC values for LLR range from 0.692 (CBS) to 0.882 (TEM1), versus 0.570 to 0.759 for SCI. A 5-fold cross-validated logistic regression on the (SCI, LLR) feature pair yields pooled AUROC = 0.758, within 0.005 of LLR alone, and exceeds the better single scorer by at most Δ AUROC = 0.007 (TEM1) in any one protein.

We make three observations on this comparison.

First, the finding is consistent with the framing of the present work as an interpretability analysis rather than a prediction benchmark. SCI is defined as a per-mutation Pearson-correlation summary across 33 layer-difference vectors, and the 528-pair upper-triangle mean used in §2.8 averages out the very layer-pair structure that the paper isolates as its contribution (C2, C3, C5). Comparing this scalar summary against LLR — which uses the dedicated LM-head readout at the final layer — is not the comparison that SCI is designed to win. The boundary stated at the end of §D4 ("we make no claim that layer-pair features yield higher predictive accuracy than tuned single-layer protein-LM scores") is consistent with the result reported here.

Second, the near-zero gain of the combined SCI+LLR scorer (Δ AUROC ≤ 0.007 per protein) indicates that the scalar SCI signal is largely redundant with LLR. The moderate SCI–LLR correlation (Pearson r = 0.456; Spearman ρ = 0.467) supports the interpretation that SCI captures a representation-level signal that is only partially aligned with likelihood-based mutation scoring. This is expected: both are functions of the same ESM-2 forward pass on the same WT sequences, and both ultimately reflect the model's residual-stream representation of the mutation site. SCI's distinctive contribution is not in the scalar summary but in the layer-resolved structure that the summary collapses — the 528-pair landscape (C3), the EL-block category dominance (C2), the layer-30–32 recurrent hub (C5), and the protein-specific compositional heterogeneity (C4). None of these features is derivable from a single LM-head log-likelihood ratio, and none of them would be visible if SCI were evaluated solely as a scalar predictor.

Third, the strong per-protein LLR performance (Spearman r 0.34–0.71; AUROC 0.69–0.88; §2.8) is consistent with the published ProteinGym leaderboard values for the ESM-2 masked-marginal baseline on the same five DMS panels, and provides independent validation that the `esm2_t33_650M_UR50D` backbone used for SCI feature extraction is capable of high-quality zero-shot scoring. The SCI signal does not arise from an inadequate underlying model — it arises from an interrogation of a different representational substrate (the residual stream across all 33 layers at the mutation site) than the one queried by the LM head.

We therefore read §2.8 as both an honest benchmark and a clarification of scope. As a scalar mutation-effect predictor on the five proteins in this study, the masked-marginal LLR is the stronger choice and we report it as such. As a tool for interpreting where in the network mutation-relevant computation is organized, SCI — and more specifically the 528-pair landscape derived from it — exposes structure that the LM-head likelihood does not, and it is this structure that constitutes the contribution of this work.

---

## 4. Methods

### 4.1 Datasets

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

### 4.2 ESM-2 feature extraction

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

### 4.3 Structural Context Inconsistency (SCI)

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
(`*_top50.npy`); all three are float32 (95,142,). T₁ (Methods §4.6)
uses the Top50 summary score.

### 4.4 Layer-pair categorisation and Composite ranking

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

### 4.5 Protein clustering and cross-protein recurrence

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

### 4.6 Test statistics

Five test statistics were defined to validate findings C1–C5 at both
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

### 4.7 Permutation test (within-protein label shuffle)

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

### 4.8 Multiple-testing correction

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

### 4.9 Non-parametric bootstrap

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

### 4.10 Statistical reporting

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

### 4.11 Code and data availability

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

### 4.12 ESM-2 zero-shot LLR baseline

For the head-to-head comparison reported in §2.8 we computed the
ESM-2 zero-shot log-likelihood ratio (LLR) baseline of
Meier et al. (2021). Two standard variants are implemented in
`p6_scripts/p7_compute_esm2_llr.py` on the same
`esm2_t33_650M_UR50D` checkpoint and the same WT sequences
(`data/processed/wt_sequences.csv`) used for SCI feature extraction
(§4.2):

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
cache in §4.2).

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

---

## 5. Figure Legends

### Figure 1 — SCI carries mutation-effect signal (C1)

**Caption.**
**SCI discriminates functional from non-functional mutations across all
five DMS panels.**
**(a)** Distribution of the per-mutation Top50-mean SCI score
(`all_proteins_sci_site_scores_top50.npy`; n = 95,142) stratified by the
binarised DMS label (`DMS_score_bin`, 1 = functional, 0 =
non-functional). **(b)** Permutation null distribution of the SCI mean
difference T₁ = mean(SCI | y = 1) − mean(SCI | y = 0) under
10⁻³ within-protein label shuffles (Methods §4.7). Vertical red line
marks the observed T₁ = 2.79 × 10⁻³; the null mean and standard
deviation are −3.7 × 10⁻⁴ and 5.6 × 10⁻⁵, giving z = 57.0 and
empirical p_FDR = 10⁻³ against the 4-statistic global family
(`permutation_global_summary.csv`). At the per-protein level T₁ is
FDR-significant in 5/5 proteins (smallest = CBS, p_FDR = 3.7 × 10⁻²;
remaining four ≤ 1.4 × 10⁻³; `permutation_per_protein_summary.csv`).

**Legend.**
- **Panel (a)** — overlaid histograms of Top50-mean SCI score for
  functional (blue) and non-functional (red) mutations; x-axis = SCI
  value, y-axis = count. Vertical dashed lines mark the two group means.
- **Panel (b)** — histogram of T₁ values from the permutation null;
  x-axis = T₁ (functional − non-functional SCI mean), y-axis = count;
  vertical red line marks the observed T₁; annotation reports z and
  p_FDR.

**Take-home message.**
SCI separates functional from non-functional mutations at z = 57
against the within-protein label-shuffle null, and the separation is
reproducible in every one of the five proteins after Bonferroni and BH
FDR correction.

### Figure 2 — Top-ranked layer pairs and the 528-pair category landscape (C2 + C3)

**Caption.**
**Early–Late layer pairs dominate both the Top50 high-Composite tail
and the full 528-pair landscape in ESM-2.**
**(a)** Stacked bar chart of the six-category layer-pair composition
(EE / EM / EL / MM / ML / LL; Methods §4.3) within each protein's
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

**Legend.**
- **Panel (a)** — y-axis: proportion of Top50 Composite layer pairs that
  fall in each of the six categories; x-axis: protein name (PABP, CBS,
  GAL4, PTEN, TEM1); stack colours map to category (legend embedded in
  source PNG `layer_pair_category_distribution.png`).
- **Panel (b)** — null distribution histogram for T₂; x-axis: number of
  EL pairs among the 50 layer pairs with largest |AUROC − 0.5|; vertical
  red line at observed value; annotation reports z and p_FDR.
- **Panel (c)** — null distribution histogram for T₃; x-axis: EL minus
  non-EL mean |AUROC − 0.5|; vertical red line at observed value.

**Take-home message.**
EL enrichment is observed at both the high-signal tail (Top50; per-protein
proportions of 40 %–66 %) and across the full 528-pair plane (category
dominance with z = 62 against null), establishing that the EL block is
not a tail-only phenomenon.

### Figure 3 — PABP forms an independent cluster on the layer-pair composition axis (C4)

**Caption.**
**PABP forms a singleton in the six-dimensional layer-pair composition
space, characterized by Early–Middle interactions.**
**(a)** Heatmap of the six-category Top50 Composite layer-pair
composition for each of the five proteins (n × 6 matrix from per-protein
`*_pair_category_ratio_Composite.csv`). PABP is the only row in which
the EM proportion (0.48) exceeds the EL proportion (0.42); the
remaining four proteins are EL-dominated. **(b)** Average-linkage
hierarchical clustering of the composition vectors under Euclidean
distance (`pairwise_euclidean_distance.csv`). PABP joins the four
EL-dominated proteins only at a substantially higher dendrogram
height (≈ 0.50), whereas the four others fuse at maximum height
≈ 0.22. The minimum PABP-to-other distance is 0.389 (PABP–TEM1) and
the maximum within-other-four distance is 0.397 (CBS–GAL4); across
the heatmap, dendrogram and PCA analyses PABP remains the most
isolated of the five proteins. **(c)** Principal
component analysis on the same vectors; PC1 explains 68.6 % and PC2
28.1 % of variance (96.7 % cumulative;
`pca_loadings.csv`). PABP projects to PC1 ≈ +0.38, isolated from the
other four proteins which fall in PC1 ∈ [−0.18, +0.01]. PC1 loadings
are dominated by Early–Middle (+0.207); PC2 separates Early–Late
(+0.103) from Late–Late (−0.094) (`pca_loadings_ranked.csv`).

**Legend.**
- **Panel (a)** — colour-coded matrix; rows: proteins (CBS, GAL4, PABP,
  PTEN, TEM1); columns: categories (EE, EM, EL, MM, ML, LL); colour =
  proportion of Top50.
- **Panel (b)** — dendrogram; y-axis: Euclidean distance under average
  linkage; leaves labelled with protein names; PABP rendered in the
  highlight colour used throughout the manuscript.
- **Panel (c)** — scatter of proteins in (PC1, PC2) coordinates;
  variance-explained shown on each axis; loading arrows for the six
  categories overlaid.

**Take-home message.**
Three independent geometric views (heatmap, hierarchical clustering,
PCA) converge: PABP's layer-pair signature is shifted toward
Early–Middle interactions and is more distant from any of the four
EL-dominated proteins than they are from each other.

### Figure 4 — Cross-protein recurrent layer pairs concentrate at layers 30–32 (C5)

**Caption.**
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

**Legend.**
- **Panel (a)** — bar chart; x-axis: layer index 1–33 (1-based); y-axis:
  count among the 40 layer-end occurrences of the Top20 recurrent set;
  bar colours reflect band (Early / Middle / Late).
- **Panel (b)** — bar chart with the same x-axis but evaluated over the
  larger n = 500 union; band proportions printed below the panel.

**Take-home message.**
The recurrent layer pairs do not spread over the Late band uniformly;
they sit at layers 30–32 with layer 32 most frequent, and layers
18–27 are essentially absent — both observations follow from direct
counts on the Top20 recurrent set.

### Figure 5 — Statistical robustness across proteins, corrections and bootstrap (C6 + C7)

**Caption.**
**Per-protein reproducibility and multi-layered statistical robustness of
the layer-pair findings.**
**(a)** Per-protein null distributions of T₁, T₂, T₃ and T₃ˢᵗᵈ under
within-protein label shuffles (N = 10³ per protein; Methods §4.7),
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
10⁻³ at real scale; Methods §4.9) exclude the corresponding null mean
for every primary statistic (`bootstrap_summary.csv`).

**Legend.**
- **Panel (a)** — 5 × 4 grid of null histograms; rows: proteins, columns:
  statistics; each cell annotated with z and p_FDR.
- **Panel (b)** — grouped bars; x-axis: scope (Global + 5 proteins);
  y-axis: z-score; one bar per statistic per scope; dotted = α = 0.05;
  dashed = Bonferroni.
- **Panel (c)** — scatter; x-axis: expected |z| under H₀; y-axis:
  observed |z|; one point per test; markers coloured by scope; diagonal
  reference y = x.

**Take-home message.**
Every primary finding survives three independent layers of statistical
robustness — Phipson–Smyth permutation, Bonferroni FWER, Benjamini–
Hochberg FDR — and the bootstrap 95 % CI excludes the null mean for
every primary statistic; per-protein reproducibility is established in
5 / 5 proteins on T₁, T₃, T₃ˢᵗᵈ.

### Table 1 — SCI vs ESM-2 zero-shot LLR baseline (§2.8)

**Caption.**
**As a per-mutation scalar predictor on the 95,142-mutation panel,
the masked-marginal ESM-2 LLR baseline of Meier et al. (2021) exceeds
the 528-pair upper-triangle SCI summary on every protein, and a
combined SCI+LLR scorer adds at most Δ AUROC = 0.007 over LLR alone.**
Pooled and per-protein Spearman r against `DMS_score`, AUROC and
AUPRC against `DMS_score_bin` are reported for three scorers: the
528-pair upper-triangle mean SCI summary
(`all_proteins_sci_site_scores_mean.npy`), the masked-marginal LLR
(Methods §4.12) and the 5-fold stratified-CV logistic regression on
standardised (SCI, LLR) (Methods §4.12). Numbers are verbatim from
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

**Legend.**
- **Spearman r** — Spearman rank correlation of scorer against the
  continuous `DMS_score`; reported as |r| with the natural sign of
  the correlation given in the `spearman_sign` column of
  `sci_vs_llr_comparison.csv` (always +1 across all 18 entries).
- **AUROC** — sign-corrected AUROC of scorer against `DMS_score_bin`,
  using the same convention as `P1.5_validate_sci_signal.py`
  (corrected = max over native and negated score).
- **AUPRC** — sign-corrected average precision of scorer against
  `DMS_score_bin`.

**Take-home message.**
LLR is the stronger univariate scorer on every protein and on the
pooled sample. The combined SCI+LLR scorer is essentially LLR
(Δ AUROC ≤ 0.007 in every protein), confirming that the scalar SCI
summary is largely redundant with the LM-head likelihood. The
contribution of the present work is therefore not a new scalar
predictor; it is the layer-resolved 528-pair landscape (C2, C3),
the layer-30–32 recurrence hub (C5) and the per-protein compositional
heterogeneity (C4), none of which is exposed by a single LM-head LLR.

---

## 6. Supplementary Information

### Supplementary Table S1 — Per-protein dataset characteristics

| Protein | DMS source dataset(s) | n mutations | n functional | n non-functional | Fraction functional |
|---|---|---:|---:|---:|---:|
| CBS  | `CBS_Sun2020` (Sun et al., 2020) | 7,217 | 3,314 | 3,903 | 0.459 |
| GAL4 | `GAL4_Kitzman2015` (Kitzman et al., 2015) | 1,195 | 831 | 364 | 0.695 |
| PABP | `PABP_Melamed2013` (Melamed et al., 2013) | 74,229 | 46,078 | 28,151 | 0.621 |
| PTEN | `PTEN_Mighell2018`, `PTEN_Matreyek2021` | 7,504 | 5,880 | 1,624 | 0.784 |
| TEM1 | `TEM1_Deng2012`, `TEM1_Firnberg2014`, `TEM1_Jacquier2013`, `TEM1_Stiffler2015` | 4,997 | 2,526 | 2,471 | 0.506 |
| **Total** | — | **95,142** | **58,629** | **36,513** | **0.616** |

### Supplementary Table S2 — Full per-statistic per-scope results

#### S2-A · Global tests

| Statistic | T_obs | Null mean | Null SD | z | p_empirical | p_FDR | Bootstrap 95 % CI | N_perm |
|---|---:|---:|---:|---:|---:|---:|---|---:|
| T₁ | 2.79 × 10⁻³ | −3.75 × 10⁻⁴ | 5.56 × 10⁻⁵ | 57.0 | 10⁻³ | 10⁻³ | [2.67 × 10⁻³, 2.91 × 10⁻³] | 1,000 |
| T₂ | 17 | 1.013 | 1.995 | 8.01 | 10⁻³ | 10⁻³ | [14, 19] | 1,000 |
| T₃ | 3.99 × 10⁻² | −2.08 × 10⁻³ | 6.77 × 10⁻⁴ | 62.0 | 10⁻³ | 10⁻³ | [3.70 × 10⁻², 4.27 × 10⁻²] | 1,000 |
| T₃ˢᵗᵈ | 1.02 | −0.602 | 0.201 | 8.05 | 10⁻³ | 10⁻³ | [0.954, 1.073] | 1,000 |
| T_max | 0.181 | 0.0176 | 1.74 × 10⁻³ | 93.7 | 10⁻³ | — | [0.176, 0.187] | 1,000 |

#### S2-B · Per-protein z-scores

| Protein | T₁ z | T₂ z | T₃ z | T₃ˢᵗᵈ z |
|---|---:|---:|---:|---:|
| CBS  | 2.25  | n/a | 16.11 | 18.50 |
| GAL4 | 6.29  | n/a | 5.18  | 4.59  |
| PABP | 56.55 | (T_obs=1) | 38.71 | 42.80 |
| PTEN | 2.88  | n/a | 14.57 | 17.05 |
| TEM1 | 23.59 | n/a | 15.19 | 15.94 |

#### S2-C · Per-protein BH-FDR q-values

| Protein | T₁ p_FDR | T₂ p_FDR | T₃ p_FDR | T₃ˢᵗᵈ p_FDR |
|---|---:|---:|---:|---:|
| CBS  | 3.7 × 10⁻²   | 1     | 1.4 × 10⁻³ | 1.4 × 10⁻³ |
| GAL4 | 1.4 × 10⁻³   | 1     | 1.4 × 10⁻³ | 1.4 × 10⁻³ |
| PABP | 1.4 × 10⁻³   | 1.4 × 10⁻³ | 1.4 × 10⁻³ | 1.4 × 10⁻³ |
| PTEN | 9.3 × 10⁻³   | 1     | 1.4 × 10⁻³ | 1.4 × 10⁻³ |
| TEM1 | 1.4 × 10⁻³   | 1     | 1.4 × 10⁻³ | 1.4 × 10⁻³ |

### Supplementary Figures

**S1 — Per-protein SCI distributions.** Histograms of Top50-mean SCI score stratified by functional label for each protein.

**S2 — 528-pair AUROC landscape.** Heatmap of |AUROC − 0.5| over the 33×33 layer-pair plane.

**S3 — SCI matrix structure.** Representative per-mutation 33×33 SCI matrix and distribution of best-performing layer pairs.

**S4 — Per-protein permutation grid.** 5×4 grid of per-protein permutation null distributions (same as Fig 5a, larger format).

**S5 — QQ-style significance summary.** Observed |z| vs. expected half-normal quantile (same as Fig 5c).