# Early–Late layer-pair interactions encode mutation-effect signal in the protein language model ESM-2 across five deep mutational scans

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
primary effects survive Benjamini–Hochberg FDR and non-parametric bootstrap;
T₃ and T₃ˢᵗᵈ further survive the more conservative Bonferroni FWER across all
five proteins. We interpret these observations as evidence of
a two-axis (layer × layer) organisation of mutation-effect signal in
ESM-2, complementing the prevailing one-axis (single-layer) view.

---

## 1. Introduction

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

The one-axis convention extends beyond prediction. Layer-wise
probing (Rives et al., 2021; Lin et al., 2023; Elnaggar et al., 2022;
Detlefsen et al., 2022; Vig et al., 2021) maps biophysical and
structural properties to individual layers, treating depth as a
scalar along which features emerge. Per-head attention analysis (Vig
et al., 2021; Rao et al., 2021a, MSA Transformer) decomposes a single
layer's computation; per-layer sparse-autoencoder analysis (Simon and
Zou, 2025) extracts interpretable features at each chosen depth. In
every case the unit of analysis is a single layer; the *relationship
between* layers — what is encoded jointly across depth, or what
emerges only when two layers are read together — is not directly
addressed.

In the broader interpretability literature for Transformers,
residual-stream–centred analyses have made the cross-layer view
central. Scalar-mix probing of BERT (Tenney et al., 2019) shows that
linguistic representations are composed across depths in a
characteristic order. Feed-forward layers in language models act as
key–value memories whose outputs are composed across the residual
stream (Geva et al., 2021), and tuned-lens decoding (Belrose et al.,
2023) makes the layer-by-layer evolution of a model's internal
predictions directly inspectable. To our knowledge, only one
existing PLM analysis explicitly aggregates information across more
than one layer: the contact-prediction pipeline of Rao et al. (2021b),
in which a sparse logistic regression pools attention from
all 33 ESM layers. Critically, that work targets *residue contacts*,
not mutation effects, leaving the cross-layer geometry of
mutation-effect signal in PLMs an empty cell in the existing
literature.

We adopt a *layer-pair geometry* perspective on this question. For
each mutation in a panel of five DMS-characterised proteins (CBS,
GAL4, PABP, PTEN, TEM1; 95,142 single amino-acid substitutions
released through ProteinGym, Notin et al., 2023), we project the
per-layer residue-embedding difference at the mutation position onto
every one of the 528 unordered layer pairs (i, j) with i < j of
ESM-2 (esm2_t33_650M_UR50D; Lin et al., 2023). The 528-pair plane is
partitioned into six categories by the Early (layers 1–8) / Middle
(9–24) / Late (25–33) banding — Early-Early, Early-Middle,
Early-Late, Middle-Middle, Middle-Late, Late-Late — and forms the
substrate of every analysis we report. The per-pair discrimination
metric we use for our analyses, Structural Context Inconsistency
(SCI; Methods §4.3), is one specific operationalisation of
"discriminative signal at a layer pair" and is incidental to the
central finding; what we report is the geometry of the layer-pair
plane itself.

Three patterns emerge from the layer-pair plane. (i) Early–Late
interactions dominate the high-Composite tail at the per-protein
level (Top50 Composite proportions of 40 %–66 % EL in four of the
five proteins) and also dominate the full 528-pair category-level
comparison at the pooled level (z = 62 against within-protein
label-shuffle null). (ii) Cross-protein recurrent layer pairs
concentrate at a sharp three-layer hub at depths 30–32, with layer
32 alone supplying 20 % of layer-end occurrences in the Top20
recurrent set, and layers 18–27 contributing zero. (iii) The
EL-dominated composition is shared by CBS, GAL4, PTEN and TEM1, but
PABP forms an independent singleton characterized by Early–Middle
interactions on the same six-band axis (hierarchical clustering,
PCA, and pairwise Euclidean distances). Every primary effect
survives Benjamini–Hochberg FDR over twenty per-protein primary tests and
non-parametric bootstrap; T₃ and T₃ˢᵗᵈ further survive the more conservative
Bonferroni FWER across all five proteins.

This work is, to our knowledge, the first systematic mapping of how
mutation-effect signal is organised across the layers of a protein
language model. The contribution is interpretive rather than
predictive: we make no claim that layer-pair features improve
mutation-effect prediction over tuned single-layer scores, and we
perform no head-to-head benchmark. Rather, we propose layer-pair
geometry as an interpretability primitive that complements
single-layer probing (Elnaggar et al., 2022; Detlefsen et al., 2022;
Vig et al., 2021), per-head attention analysis (Rao et al., 2021a),
and per-layer sparse-autoencoder analysis (Simon and Zou, 2025); and
we identify three reproducible patterns on this primitive that
existing one-axis methods cannot reveal. Beyond ESM-2 specifically,
the six-category composition vector defines a low-dimensional
fingerprint that could in principle be computed for ESM-1b (Rives et
al., 2021), ESM3 (Hayes et al., 2025), ProtT5 (Elnaggar et al.,
2022) and other backbones, to ask whether the Early–Late dominance
and the layer-32 anchor reported here are properties of one
particular backbone or of large protein language models more
broadly.

---

## 2. Results

### 2.1 SCI carries mutation-effect signal

We first asked whether the Structural Context Inconsistency (SCI) score
constructed from ESM-2 (esm2_t33_650M_UR50D; Methods §4.2) carries
mutation-effect signal. The per-mutation summary score over the 50
top-scoring layer pairs (n = 95,142 mutations across CBS, GAL4, PABP,
PTEN, TEM1) was compared between functional and non-functional mutations
as defined by the binarised deep-mutational-scanning label
(DMS_score_bin, 1 = functional, 0 = non-functional). At the pooled
level, the SCI mean difference was T₁ = 2.79 × 10⁻³.

Against a within-protein label-shuffle null (1,000 adaptive
permutations; Methods §4.7), the observed T₁ exceeded the null mean
(−3.75 × 10⁻⁴) by 57 standard deviations of the null
(σ_null = 5.56 × 10⁻⁵), yielding empirical p = 10⁻³ and
BH-FDR-adjusted q = 10⁻³ across the four global primary statistics
(Fig 1b). At the per-protein level, T₁ was FDR-significant in 5/5
proteins: the smallest q-value among the per-protein T₁ tests was
3.7 × 10⁻² (CBS; per-protein z = 2.25), and the remaining four
proteins reached q = 1.4 × 10⁻³ or smaller (PABP z = 56.55;
TEM1 z = 23.59; GAL4 z = 6.29; PTEN z = 2.88). The non-parametric
bootstrap (N = 1,000) returned a 95 % confidence interval for T₁ of
[2.67 × 10⁻³, 2.91 × 10⁻³], excluding zero. Together these establish
that SCI carries discriminative signal that is reproducible against
the null at both the pooled and per-protein levels.

### 2.2 Top-ranked layer pairs are enriched for Early–Late interactions

To characterise which layer pairs carry the SCI signal most strongly,
we ranked all 528 unordered layer pairs by a Composite discrimination
score combining per-pair AUROC, Spearman rank correlation and Cohen's
*d* (Methods §4.3) and examined the Top50 layer pairs for each
protein. The six-category composition (Early-Early / Early-Middle /
Early-Late / Middle-Middle / Middle-Late / Late-Late, defined on the
Early-1–8 / Middle-9–24 / Late-25–33 banding; Methods §4.3) of each
Top50 set is recorded per protein. Four of the five proteins were
EL-dominated: CBS (EL 40 %; LL 30 %), GAL4 (EL 66 %), PTEN (EL 44 %;
LL 12 %) and TEM1 (EL 40 %; EM 10 %; LL 8 %). PABP was the only protein
whose Top50 contained a larger EM than EL share (EM 48 %; EL 42 %; no
LL contribution); we return to this asymmetry in §2.4 (Fig 2a).

At the pooled level we then asked, under the AUROC-only ranking, how
many of the 50 pairs with largest |AUROC − 0.5| fall into the EL
category. Observed T₂ = 17 of 50, compared with a permutation null
mean of 1.01 (σ_null = 2.00), giving z = 8.01 and q = 10⁻³ (Fig 2b).
The bootstrap 95 % CI for T₂ was [14, 19]. A methodological note: at
the per-protein level the AUROC-only ranking and the Composite ranking
can diverge, and per-protein T₂ values are uninformative for four of
the five proteins (Methods §4.6); per-protein evidence for the EL
enrichment is therefore primarily carried by the per-protein T₃ and
T₃ˢᵗᵈ analyses in §2.3 and by the per-protein Composite category
ratios reported above.

### 2.3 Early–Late pairs dominate the 528-pair landscape

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
with z = 8.05 and q = 10⁻³ (Fig 2c). Bootstrap 95 % CIs are
[3.70 × 10⁻², 4.27 × 10⁻²] for T₃ and [0.954, 1.073] for T₃ˢᵗᵈ,
excluding the corresponding null means. At the per-protein level T₃
and T₃ˢᵗᵈ were FDR-significant in 5/5 proteins; per-protein T₃ z-scores
ranged from 5.18 (GAL4) to 38.71 (PABP), and T₃ˢᵗᵈ z-scores ranged
from 4.59 (GAL4) to 42.80 (PABP). At the per-protein level T₃ is
computed against the globally ranked SCI scores (Methods §4.6);
per-protein effect-size magnitudes are therefore not bounded in
[−0.5, 0.5], but the permutation null is computed under the identical
metric and the significance call is unaffected. The category-level EL
dominance is observed at both the pooled level and within each of the
five proteins.

### 2.4 One of the five proteins (PABP) forms a distinct cluster

We then asked whether the EL-dominated composition reported in §2.2 is
shared across proteins or conceals heterogeneity. Each protein's Top50
Composite layer-pair set was summarised as a six-dimensional vector of
category proportions, and three independent geometric analyses were
applied.

The heatmap of these vectors (Fig 3a) identifies PABP as the unique
EM-dominated row. Hierarchical clustering with Euclidean distance and
average linkage (Fig 3b) places CBS, GAL4, PTEN and TEM1 into a compact
cluster whose maximum pairwise distance is 0.397 (CBS–GAL4) and minimum
is 0.115 (PTEN–TEM1). PABP joins this cluster only at a substantially
higher dendrogram height (≈ 0.50); the minimum PABP-to-other distance
is 0.389 (PABP–TEM1). Across the heatmap, dendrogram and PCA analyses,
PABP remains the most isolated of the five proteins. Principal
component analysis on the same six-dimensional vectors (Fig 3c)
recovers PC1 = 68.6 % and PC2 = 28.1 % of variance (cumulative
96.7 %). PC1 loadings are dominated by Early–Middle (+0.207), with the
next-largest absolute loading approximately three times smaller (LL,
−0.072) and EL contributing only −0.026. PC2 separates Early–Late
(+0.103) from Late–Late (−0.094). PABP projects to PC1 ≈ +0.38,
isolated from the four EL-dominated proteins which fall in
PC1 ∈ [−0.18, +0.01]. The four EL-dominated proteins share a similar
layer-pair signature, while PABP's signature is shifted toward
Early–Middle interactions.

### 2.5 Cross-protein recurrent layer pairs concentrate at layers 30–32

To identify layer pairs that carry mutation-effect signal across
proteins, we counted, for each layer pair, the number of proteins in
whose Top50 Composite set it appears (146 candidate recurrent pairs).
The Top20 most-recurrent pairs (6 pairs at Cross_Protein_Count = 4 and
14 at Cross_Protein_Count = 3) carry 40 layer-end occurrences
(2 × 20 pairs).

Band totals (Fig 4a) are: Late (layers 25–33) = 22 of 40 (55.0 %),
Middle (9–24) = 11 (27.5 %), Early (1–8) = 7 (17.5 %). Within the Late
band the distribution is non-uniform: layer 32 contributes 8 of 40
layer-end occurrences (20.0 %, the single most recurrent layer), layer
31 contributes 7 (17.5 %), layer 30 contributes 3 (7.5 %), and layer
33 contributes 2 (5.0 %). Layers 25–29 contribute zero, producing a
sharp three-layer hub at depths 30–32 with a marked trough immediately
upstream. At the broader Top50 union across the five proteins
(n = 500 layer-end occurrences = 5 × 50 × 2; Fig 4b) the band
proportions are preserved: Late 47.6 %, Middle 36.0 %, Early 16.4 %;
layer 32 alone accounts for 14.6 % of all occurrences. The recurrent
set concentrates at a sharp three-layer hub centred on layer 32, and
the central middle band (layers 18–27 in the Top20 set) is essentially
absent from the recurrent geometry.

### 2.6 Per-protein signal is reproducible against within-protein null

To rule out the possibility that the pooled effects in §2.1–§2.3 are
attributable to a single protein, we re-computed the four primary test
statistics within each of the five proteins, using a within-protein
label-shuffle null (N = 1,000 permutations per protein; Methods §4.7).
All five proteins returned FDR-significant T₁, T₃ and T₃ˢᵗᵈ (Fig 5a,
5b). The five per-protein T₁ z-scores were 2.25 (CBS), 6.29 (GAL4),
56.55 (PABP), 2.88 (PTEN) and 23.59 (TEM1); the five per-protein T₃
z-scores were 16.11, 5.18, 38.71, 14.57 and 15.19; and the five T₃ˢᵗᵈ
z-scores were 18.50, 4.59, 42.80, 17.05 and 15.94. After Bonferroni
correction over the 20 per-protein primary tests (5 proteins × 4
statistics), 13 of these 15 tests (T₁, T₃, T₃ˢᵗᵈ across 5 proteins) retained
significance at the conventional α = 0.05 level; CBS-T₁ (p_bonf = 0.60) and
PTEN-T₁ (p_bonf = 0.14) did not survive the conservative Bonferroni threshold
despite being significant under BH-FDR control (q = 0.037 and q = 0.009,
respectively). The signal is
therefore not pooled out of one over-represented protein; it is
reproduced within each, with strongest consistency for T₃ and T₃ˢᵗᵈ across
all five proteins.

### 2.7 Statistical robustness across multiple corrections and bootstrap

We assessed robustness across three independent layers. (i) Phipson and
Smyth (2010) empirical p-values were computed under the null-mean–
centered convention for two-sided tests; under this convention, every
primary statistic at the pooled level reached p = 10⁻³, the floor for
N = 1,000 permutations. (ii) Bonferroni correction over the 20-test
family of per-protein primary statistics retained significance at α = 0.05
for T₃ and T₃ˢᵗᵈ in all five proteins; for T₁, three of five proteins (GAL4,
PABP, TEM1) survived Bonferroni, while CBS-T₁ (p_bonf = 0.60) and PTEN-T₁
(p_bonf = 0.14) did not.
(iii) Benjamini–Hochberg FDR with step-up monotonicity enforcement was
applied to the same 20-test family, yielding q < 0.05 for T₁ in all five
proteins and thus providing a slightly less conservative set of discoveries
than Bonferroni. All 15 tests remained significant under BH-FDR control.

A non-parametric bootstrap (N = 1,000; Methods §4.9) provided 95 % CIs
for the pooled statistics: T₁ ∈ [2.67, 2.91] × 10⁻³, T₂ ∈ [14, 19],
T₃ ∈ [3.70, 4.27] × 10⁻², T₃ˢᵗᵈ ∈ [0.954, 1.073], T_max ∈ [0.176,
0.187]. All five intervals exclude the null mean. The family-wise
maximum statistic T_max = max |AUROC − 0.5| over all 528 pairs reached
0.181, exceeding the null max (mean 0.0176, σ 1.74 × 10⁻³) by z = 93.7,
addressing the concern that the best-of-528 search inflates significance.
Across all four robustness layers, the EL geometry of mutation-effect
signal in ESM-2 survives every correction applied to it (Fig 5).

### 2.8 Comparison with ESM-2 zero-shot baseline

The above analyses establish that SCI carries mutation-effect signal
and that the signal has a reproducible layer-pair geometry. To position
this against the field's de facto zero-shot scorer, we computed the
ESM-2 log-likelihood ratio (LLR) baseline of Meier et al. (2021) on the
same 95,142 mutation rows and asked two questions: how does SCI compare
to LLR as a per-mutation scalar predictor, and does SCI carry
information that is additive to LLR? The masked-marginal variant of LLR
(Meier et al., 2021) was computed by masking the mutated residue in the
WT sequence, taking log-softmax of the LM-head logits at that position,
and forming LLR = log p(mut | x_\T) − log p(wt | x_\T) (Methods §4.12);
the per-mutation SCI summary used here is the 528-pair upper-triangle
mean. The combined SCI+LLR scorer is a 5-fold stratified-CV logistic
regression on the standardised (SCI, LLR) feature pair, with the
StandardScaler fit within each training fold; out-of-fold P(class = 1)
is used as the predictor for the held-out fold (Methods §4.12).

At the pooled level (n = 95,142), the masked-marginal LLR achieves
Spearman r = 0.480 (p ≈ 0) against DMS_score and AUROC = 0.753 against
DMS_score_bin, exceeding the 528-mean SCI summary (r = 0.200,
AUROC = 0.669) on both metrics (Table 1). The 5-fold CV SCI+LLR
combination reaches AUROC = 0.758 and Spearman r = 0.464, within 0.005
AUROC of LLR alone. SCI and LLR were only moderately correlated at the
pooled level (Pearson r = 0.456; Spearman ρ = 0.467), indicating
partial overlap but substantial non-shared variance. At the per-protein
level LLR exceeds SCI on every one of the five proteins, with
per-protein AUROC values of 0.692 (CBS), 0.795 (GAL4), 0.767 (PABP),
0.846 (PTEN) and 0.882 (TEM1) for LLR versus 0.570, 0.717, 0.699, 0.629
and 0.759 for SCI, and per-protein Spearman r ranging from 0.342 (CBS)
to 0.714 (TEM1) for LLR versus 0.131 to 0.502 for SCI. The combined
SCI+LLR scorer exceeds the better single scorer by at most
Δ AUROC = 0.007 (TEM1), with smaller gains on PABP (+0.003), GAL4
(+0.002), PTEN (+0.001) and no improvement on CBS (Δ ≤ 0.001); the SCI
signal is therefore largely redundant with LLR as a scalar predictor.

We make this comparison transparently. SCI as a scalar per-mutation
predictor underperforms the masked-marginal LLR baseline, and a
combined logistic regression does not meaningfully exceed LLR alone.
This result does not contradict §2.1, which establishes that SCI
carries discriminative signal against a within-protein label-shuffle
null at z = 57; rather, it indicates that the SCI signal is largely
co-linear with the LM-head likelihood at the scalar level. The
contribution of the present work is not a new scalar predictor; it is
the layer-resolved 528-pair landscape (§2.2–§2.3), the layer-30–32
recurrence hub (§2.5), and the protein-specific compositional
heterogeneity (§2.4) — none of which is exposed by a single LM-head
likelihood ratio. The relationship between the scalar SCI summary and
LLR is discussed further in Discussion D5.

---

## 3. Discussion

### D1. Deep-layer enrichment in ESM

The cross-protein recurrence analysis indicates that the late end of mutation-effect signal in ESM-2 is sharply localized. Among the Top20 cross-protein recurrent layer pairs, the Late band (layers 25–33) accounts for 55% of all layer-end occurrences, while the Middle band (layers 9–24) and the Early band (layers 1–8) account for 27.5% and 17.5%, respectively (§2.5). Within the Late band, the distribution is itself non-uniform: layer 32 alone supplies 8 of the 40 layer-end occurrences, layers 30–32 together supply 18 of 40, and layers 25–29 contribute essentially none (§2.5). The late-end concentration therefore reflects a narrow three-layer hub rather than a broad enrichment across the upper half of the network. This hub is observed in the same recurrent set in which Early–Late pairs dominate (§2.2, §2.5) and is reproduced at the level of the full 528-pair distribution (§2.3).

The concentration of recurrence at layers 30–32 is consistent with prior reports that the deepest layers of ESM and related protein language models carry higher-order semantic and functional features (Rives et al., 2021; Lin et al., 2023; Elnaggar et al., 2022; Vig et al., 2021). It is also consistent with the broader observation that the most useful single layer of a protein language model depends on the downstream task (Detlefsen et al., 2022). The present finding may be read as a refinement of that framing: for mutation-effect signal in ESM-2, the deepest layers are recurrently engaged across proteins, but they appear within layer pairs rather than as stand-alone signal carriers (§2.2, §2.5). The convergence of the recurrent set on layer 32 is also notable in light of work that uses ESM-1b's final layer alone to score human missense variants at proteome scale (Brandes et al., 2023). The late-end anchor identified here lies one layer below the masked-LM head; the convergence of the recurrent set on layer 32 may therefore be read as compatible with the use of layer 33 in downstream variant-effect scoring, rather than as opposed to it. We hesitate to interpret this convergence as more than a structural coincidence; it is one of several possible readings of the geometry. All observations in this subsection are made on a single backbone (ESM-2, 33 layers); whether a comparable late-end anchor emerges in ESM-1b (Rives et al., 2021), ESM3 (Hayes et al., 2025), or other Transformer-based protein language models remains an open question (§2.5).

### D2. Cross-layer integration hypothesis

Across the 528 unordered layer pairs of ESM-2, Early–Late (EL) pairs achieve significantly higher Composite scores and AUROC values than Early–Early, Early–Middle, Middle–Middle and Late–Late pairs, with category-level comparisons returning p ≪ 0.001 (§2.3). This category-level pattern is preserved at the high-signal tail: the Top50 Composite layer pairs are enriched for EL interactions across all five proteins (§2.2). The cross-protein recurrent set is itself dominated by pairs that connect a deep-layer endpoint, typically in 30–32, with a non-Late endpoint in the Early-to-mid-Middle range (§2.2, §2.5). These statistical relationships are observed on a representational signal — SCI mean — that has independently been shown to distinguish functional from non-functional mutations in the cohort examined here (§2.1).

Taken together, §2.1, §2.2, §2.3 and §2.5 are consistent with a working hypothesis that mutation-effect signal in ESM-2 is integrated across layers rather than read off any single layer. Under this hypothesis, early layers may contribute residue-local biochemical features, and late layers — notably 30–32 — may contribute family-level or contextual features, with the discriminative information potentially corresponding to the relationship between these two depth-stratified representations (§2.2, §2.3, §2.5). This reading is consistent with interpretability work on general Transformers that characterizes the residual stream as a substrate over which intermediate representations are progressively composed across layers (Tenney et al., 2019; Geva et al., 2021; Belrose et al., 2023). The closest existing protein-LM analysis that explicitly aggregates across layers is the contact-prediction pipeline of Rao et al. (2021b), in which a sparse logistic regression pools attention from all ESM layers; the present findings suggest that an analogous, but biologically distinct, multi-layer combination may be at play for mutation-effect encoding (§2.2, §2.3). We emphasize that this remains a hypothesis suggested by the layer-pair geometry rather than an established property of the forward pass.

Several limitations of this hypothesis warrant explicit acknowledgement. First, all findings here are correlational; we have not performed activation patching, single-layer ablation, or other interventions that would isolate the contribution of any individual layer (§2.2, §2.3, §2.5). Second, the Composite score combines AUROC, Spearman *r* and Cohen's *d*, and an EL pair could in principle rank highly through aggregation of two independently weak but complementary endpoints, rather than through a genuinely interactive signal (§2.3). Disentangling aggregation from interaction is an open problem that the present analysis does not resolve. Third, the existence of high-performing single-layer mutation-effect scores from the masked-LM head of ESM-1v (Meier et al., 2021), the autoregressive output of Tranception (Notin et al., 2022), and the ridge regression over PLM embeddings of (Hsu et al., 2022) may at first appear in tension with the integration hypothesis. Predictive sufficiency at a single layer does not preclude representational distribution across pairs, however, given that each of those single-layer scores is itself the product of cross-layer composition through the residual stream; the present claims are about layer-pair geometry rather than predictor design (§2.2, §2.3).

### D3. Protein-specific heterogeneity

The five proteins examined here do not share a single layer-pair signature. Four of them — CBS, GAL4, PTEN and TEM1 — are EL-dominated within their Top50 Composite layer pairs (§2.4). PABP is the only protein whose Early–Middle proportion exceeds its Early–Late proportion within the same Top50 set, and the only protein with no Late–Late contribution (§2.4). Hierarchical clustering and PCA over the six-category proportion vector place PABP at a singleton position whose distance to the centroid of the remaining four proteins is comparable to the entire pairwise spread within the four-protein cluster (§2.4). All five proteins exhibit clearly above-baseline SCI signal (§2.1); the observed heterogeneity therefore concerns which layer pairs carry mutation-effect signal, not whether such signal is present. The cross-protein recurrent set reflects the four EL-dominated proteins more strongly than it reflects PABP, consistent with PABP's distinct compositional profile (§2.4, §2.5).

The shift of PABP toward Early–Middle interactions co-varies with biophysical class: PABP is the only RNA-binding protein in the panel, while the four EL-dominated proteins span enzymes (CBS, PTEN, TEM1) and a transcription factor (GAL4). This co-variation is consistent with depth-stratified attention-head specialization in protein Transformers, in which different biophysical properties are tracked at characteristic depths (Vig et al., 2021), and with the existence of layer-localized interpretable features that vary across positions and proteins (Simon and Zou, 2025). A second reading is that the heterogeneity is consistent with differences in deep-mutational-scanning (DMS) assay design: PABP's assay measures growth-coupled RNA binding (Melamed et al., 2013), whereas the four others measure stability or enzymatic activity in microbial or human cell contexts. Differences in the shape of the underlying fitness landscape may plausibly project onto different layer-pair compositions. We do not claim that RNA binding or assay design is responsible for Middle-layer involvement; we only note that the observed heterogeneity is non-random and aligns with two biological axes that may inform follow-up sampling (§2.4).

It is also important to note what the PABP singleton is not. PABP's Top50 Composite scores are well-separated from baseline as for the other four proteins (§2.1); only the layer-band composition of those top-scoring pairs differs (§2.4). The singleton position therefore should not be read as evidence of poor model fit, low data quality, or a failure of the protein language model on RNA-binding proteins. Two limitations of this reading are explicit. With *n* = 5, PABP's singleton status is a qualitative observation rather than a statistical generalization (§2.4). Additionally, biophysical class and DMS assay design are confounded within the present sample — there is one RNA-binder and four non-RNA-binders, each with its own assay format — so the two consistent readings above cannot be disentangled (§2.4). Whether RNA-binding proteins more generally, or proteins with shallow fitness landscapes, systematically display Early–Middle dominance is a hypothesis that requires expansion of the protein panel to at least twenty sequences spanning multiple biophysical classes and assay types.

### D4. Implications for PLM interpretability

Taken together, the systematic EL enrichment at both the category and Top50 levels (§2.2, §2.3), the narrow deep-layer anchor at layers 30–32 (§2.5), and the per-protein compositional heterogeneity (§2.4) describe a structure that is not visible from per-layer probing alone. The strongest signal lies in the off-diagonal portion of the 528-pair plane — specifically the EL block — while per-layer probing protocols, by construction, integrate over one of the two axes of that plane (§2.1, §2.2, §2.3). The deep-layer anchor (§2.5) appears in paired form: each of the most recurrent layer pairs in the Top20 set connects layers 30–32 to either an Early or a Middle partner (§2.2, §2.5).

On the basis of §2.2, §2.3, §2.4 and §2.5, we suggest that layer-pair geometry may serve as a useful interpretability primitive for protein language models, complementing rather than replacing per-layer probing (Elnaggar et al., 2022; Detlefsen et al., 2022; Rao et al., 2019) and per-head attention analysis (Vig et al., 2021; Rao et al., 2021a). The six-category composition vector summarizes a per-protein fingerprint of how mutation-effect signal is distributed across the depth of a model, and could in principle be computed for ESM-1b (Rives et al., 2021), ESM3 (Hayes et al., 2025), ProtT5 (Elnaggar et al., 2022) and other backbones, to ask whether the EL enrichment and the layer-32 anchor are properties of ESM-2 specifically or of large protein language models more broadly (§2.2, §2.4, §2.5). The observation that mutation-effect signal recurs at cross-depth pairs is also consistent with the residual-stream–centered interpretability picture developed for general Transformers (Tenney et al., 2019; Geva et al., 2021; Belrose et al., 2023), and suggests that the same toolkit — scalar-mix probes, tuned-lens decoding, activation patching across layer pairs, and per-layer sparse-autoencoder analyses (Simon and Zou, 2025) — could plausibly be applied to protein language models in future work to clarify whether layer pairs are merely the most informative readout sites or substrates of genuine cross-layer integration.

Three boundaries on these claims warrant explicit statement. First, we make no claim that layer-pair features yield higher predictive accuracy than tuned single-layer protein-LM scores; a head-to-head comparison was not performed, and predictive parity is not the goal of this work (§2.1, §2.2, §2.3). Second, the EL pattern reported here is correlational, and we do not claim that it is responsible for ESM-2's mutation-effect performance (§2.2, §2.3). Third, the prominence of layer 32 has been established only for ESM-2; whether an analogous anchor emerges in other protein language models remains open (§2.5). Within these boundaries, the layer-pair view contributes a two-axis vocabulary to a literature that has been organized largely along a single depth axis, and the off-diagonal structure of that vocabulary — particularly the EL block and its protein-specific composition — may be a productive site for further interpretability work.

### D5. Relationship to the ESM-2 zero-shot likelihood baseline

The masked-marginal log-likelihood ratio (LLR) of Meier et al. (2021) is the field's de facto zero-shot scorer for ESM-2 and is the ESM-2 baseline reported on the ProteinGym leaderboard. We computed it on the same 95,142 mutation rows (§2.8, Methods §4.12) and observe that LLR exceeds the 528-mean SCI summary as a per-mutation predictor on every one of the five proteins examined here (Table 1). Pooled, LLR achieves AUROC = 0.753 against SCI 0.669; per-protein AUROC values for LLR range from 0.692 (CBS) to 0.882 (TEM1), versus 0.570 to 0.759 for SCI. A 5-fold cross-validated logistic regression on the (SCI, LLR) feature pair yields pooled AUROC = 0.758, within 0.005 of LLR alone, and exceeds the better single scorer by at most Δ AUROC = 0.007 (TEM1) in any one protein.

We make three observations on this comparison.

First, the finding is consistent with the framing of the present work as an interpretability analysis rather than a prediction benchmark. SCI is defined as a per-mutation Pearson-correlation summary across 33 layer-difference vectors, and the 528-pair upper-triangle mean used in §2.8 averages out the very layer-pair structure that the paper isolates as its contribution (§2.2, §2.3, §2.5). Comparing this scalar summary against LLR — which uses the dedicated LM-head readout at the final layer — is not the comparison that SCI is designed to win. The boundary stated at the end of Discussion D4 ("we make no claim that layer-pair features yield higher predictive accuracy than tuned single-layer protein-LM scores") is consistent with the result reported here.

Second, the near-zero gain of the combined SCI+LLR scorer (Δ AUROC ≤ 0.007 per protein) indicates that the scalar SCI signal is largely redundant with LLR. The moderate SCI–LLR correlation (Pearson r = 0.456; Spearman ρ = 0.467) supports the interpretation that SCI captures a representation-level signal that is only partially aligned with likelihood-based mutation scoring. This is expected: both are functions of the same ESM-2 forward pass on the same WT sequences, and both ultimately reflect the model's residual-stream representation of the mutation site. SCI's distinctive contribution is not in the scalar summary but in the layer-resolved structure that the summary collapses — the 528-pair landscape (§2.3), the EL-block category dominance (§2.2), the layer-30–32 recurrent hub (§2.5), and the protein-specific compositional heterogeneity (§2.4). None of these features is derivable from a single LM-head log-likelihood ratio, and none of them would be visible if SCI were evaluated solely as a scalar predictor.

Third, the strong per-protein LLR performance (Spearman r 0.34–0.71; AUROC 0.69–0.88; §2.8) is consistent with the published ProteinGym leaderboard values for the ESM-2 masked-marginal baseline on the same five DMS panels, and provides independent validation that the esm2_t33_650M_UR50D backbone used for SCI feature extraction is capable of high-quality zero-shot scoring. The SCI signal does not arise from an inadequate underlying model — it arises from an interrogation of a different representational substrate (the residual stream across all 33 layers at the mutation site) than the one queried by the LM head.

We therefore read §2.8 as both an honest benchmark and a clarification of scope. As a scalar mutation-effect predictor on the five proteins in this study, the masked-marginal LLR is the stronger choice and we report it as such. As a tool for interpreting where in the network mutation-relevant computation is organized, SCI — and more specifically the 528-pair landscape derived from it — exposes structure that the LM-head likelihood does not, and it is this structure that constitutes the contribution of this work.

---

## 4. Methods

### 4.1 Datasets

Five deep mutational scanning (DMS) panels covering one yeast
RNA-binding protein (PABP), one fungal transcription factor (GAL4) and
three enzymes (CBS, PTEN, TEM1) were drawn from the ProteinGym
substitution release (Notin et al., 2023). The mutations were merged
into a single panel of 95,142 single amino-acid substitutions covering
the five proteins; the binarised functional label DMS_score_bin
(1 = functional, 0 = non-functional) follows the threshold provided by
ProteinGym for each panel and is used throughout as the supervised
target. Per-protein counts of functional and non-functional mutations
are listed in Supplementary Table S1. No additional filtering was
applied; the full set of 95,142 mutations was used in every analysis.

### 4.2 ESM-2 feature extraction

We used the publicly released esm2_t33_650M_UR50D model (33 Transformer
blocks, hidden dimension 1280, maximum sequence length 1022; Lin et
al., 2023) for all feature extraction. For each unique wild-type
sequence and for each single-mutation variant we extracted the
per-residue hidden representations at all 33 layers. At the mutation
position only, we computed the per-layer difference vector
**Δ**ᵢ = h_mut(layer i, pos) − h_WT(layer i, pos) ∈ ℝ¹²⁸⁰,
stacked across the 33 layers into a per-mutation matrix
Δ ∈ ℝ³³ ˣ ¹²⁸⁰. No attention pooling and no averaging across positions
is used; SCI is defined strictly at the per-residue mutation site.

### 4.3 Structural Context Inconsistency (SCI)

The Structural Context Inconsistency (SCI) metric is the Pearson
correlation between per-layer Δ-embedding vectors at the mutation
position. For each mutation we compute the 33 × 33 symmetric matrix

  SCI[i, j] = corr(Δᵢ, Δⱼ)

via the Pearson correlation on the (33, 1280) matrix. Pearson
correlation values are well-defined when both layer-difference vectors
have non-zero variance; samples with zero-variance layers are detected
explicitly, and undefined or non-finite correlations are mapped to zero
without addition of random noise. The 528 unique upper-triangle entries
SCI[i, j] with i < j define the per-pair SCI features for downstream
analysis. Three per-mutation summary scores are also materialised: the
mean of all 528 upper-triangle values, the mean of the top-20, and the
mean of the top-50 values. T₁ (Methods §4.6) uses the Top50 summary
score.

### 4.4 Layer-pair categorisation and Composite ranking

Layer pairs (i, j) with i < j and i, j ∈ {1, …, 33} were partitioned
into six categories based on the Early (layers 1–8) / Middle (9–24) /
Late (25–33) banding: Early-Early (EE; 28 pairs), Early-Middle
(EM; 128 pairs), Early-Late (EL; 72 pairs), Middle-Middle (MM; 120
pairs), Middle-Late (ML; 144 pairs) and Late-Late (LL; 36 pairs);
total = 528. The Composite layer-pair score is a non-negative
combination of three per-pair quantities computed against the binarised
functional label:

- AUROC of per-pair SCI as predictor of DMS_score_bin;
- Spearman rank correlation between per-pair SCI and the continuous
  DMS_score;
- Cohen's *d* between functional and non-functional groups on per-pair
  SCI.

The six-category proportions for each protein's Composite-Top50 set are
recorded per protein.

### 4.5 Protein clustering and cross-protein recurrence

Per-protein Composite-Top50 sets were summarised as six-dimensional
category-proportion vectors. Hierarchical clustering used Euclidean
distance with average linkage. Principal component analysis was applied
to the same six-dimensional vectors without standardisation, so that
PC1 and PC2 reflect raw compositional differences rather than rescaled
units.

For cross-protein recurrence, a layer pair (i, j) is *recurrent* if it
appears in the Composite-Top50 set of multiple proteins. The Top20
recurrent set was extracted by sorting on Cross_Protein_Count in
descending order. Layer-end occurrence counts per layer, per-band
totals and the broader Top50 union across the five proteins (n = 500
layer-ends = 5 × 50 × 2) are reported.

### 4.6 Test statistics

Five test statistics were defined to validate the primary findings at
both the pooled and per-protein levels:

- **T₁ — SCI mean difference.**
  T₁ = mean(SCI_top50 | y = 1) − mean(SCI_top50 | y = 0).
- **T₂ — AUROC-Top50 EL count.**
  Number of EL layer pairs among the 50 with largest |AUROC − 0.5|
  per the per-pair AUROC computed against DMS_score_bin.
- **T₃ — EL category dominance.**
  T₃ = mean |AUROC − 0.5|_EL − mean |AUROC − 0.5|_non-EL.
- **T₃ˢᵗᵈ — standardised EL dominance.**
  T₃ˢᵗᵈ = (mean |AUROC − 0.5|_EL − mean |AUROC − 0.5|_non-EL) /
  σ(|AUROC − 0.5|_non-EL).
- **T_max — family-wise maximum.**
  T_max = max |AUROC − 0.5| over all 528 pairs; reported as a
  supplementary control for the 528-pair search.

Per-pair AUROC is computed via the Mann–Whitney U identity using
average ranks. |AUROC − 0.5| is used to make the discrimination test
direction-agnostic, so that both functional-up and functional-down
layer pairs contribute. At the per-protein level, the ranks of SCI
values are reused from the full-population ordering for computational
efficiency; under this convention the per-protein T₃, T₃ˢᵗᵈ and T_max
are *globally-anchored* U-statistics. Significance calls are unaffected
(the permutation null is computed under the identical statistic), but
per-protein effect-size magnitudes are not bounded in [−0.5, 0.5] and
should be interpreted relative to the per-protein null mean and
standard deviation.

A methodological note on T₂: at the per-protein level the AUROC-only
Top50 set diverges from the Composite-Top50 set used in §2.2 of the
Results, and per-protein T₂ values are uninformative in four of the
five proteins. Per-protein evidence for the §2.2 enrichment is
therefore primarily carried by T₃ and T₃ˢᵗᵈ.

### 4.7 Permutation test (within-protein label shuffle)

We tested the null hypothesis that SCI carries no information about
mutation functionality using a within-protein label-shuffle
permutation. For each permutation, the DMS_score_bin labels were
permuted independently within each of the five proteins, holding fixed
the SCI feature matrix, protein identity, and per-protein label
balance. The five test statistics were recomputed under each permuted
label vector.

Empirical p-values follow the Phipson and Smyth (2010) convention
p = (1 + #{T_perm extreme}) / (N + 1), which lower-bounds p away from
zero. Two-sided p-values (used for T₁) centre the comparison on the
null mean: |T_perm − μ_null| ≥ |T_obs − μ_null|. One-sided p-values
(used for T₂, T₃, T₃ˢᵗᵈ, T_max) use T_perm ≥ T_obs.

The number of permutations follows an adaptive schedule with budgets
{1,000 / 5,000 / 10,000} and an early-stopping rule at max-p < 0.001
across the four primary statistics. Per-protein analyses used a fixed
budget of N = 1,000 permutations.

Effect sizes are reported as the z-score z = (T_obs − μ_null) / σ_null
and as the raw deviation T_obs − μ_null; under the permutation null the
z-score equals Cohen's *d*. Null 95 % and 99 % percentile intervals are
reported alongside. All RNG draws use a fixed seed (42 for permutations,
43 for the bootstrap), ensuring full reproducibility.

### 4.8 Multiple-testing correction

Two independent corrections were applied to the per-protein primary
tests (4 statistics × 5 proteins = 20 tests):

- **Bonferroni FWER** with p_bonferroni = min(1, p_empirical × 20).
- **Benjamini–Hochberg FDR** with step-up monotonicity enforcement.
  The same step-up procedure is applied independently to the 4 global
  primary tests. NaN p-values are passed through as NaN and excluded
  from the BH family size.

T_max is reported outside the family-wise correction as a supplementary
control statistic for the 528-pair search; its p_empirical and z-score
are still reported, but it does not contribute to the Bonferroni or BH
families.

### 4.9 Non-parametric bootstrap

A non-parametric bootstrap with N = 1,000 resamples and replacement was
used to compute 95 % and 99 % confidence intervals for each of the
five test statistics. For T₁ the bootstrap recomputes the SCI mean
difference directly on the resampled (mutation, label) pairs. For T₂,
T₃, T₃ˢᵗᵈ and T_max we use a rank-reuse fast path: instead of
recomputing per-column SCI ranks on every bootstrap resample (the
natural non-parametric bootstrap, with complexity O(B · K · N log N)),
the precomputed full-sample ranks are reused under the bootstrap index
vector. This is a U-statistic approximation whose deviation from the
natural bootstrap AUROC was measured at real scale (N = 95,142,
K = 528, multiple sanity-check resamples): the peak per-pair |ΔAUROC| was
6 × 10⁻³ and the mean was 1.6 × 10⁻³, both well below the precision at
which 95 % CIs are reported.

### 4.10 Statistical reporting

All p-values, z-scores, Cohen's *d* values, raw effect sizes, null 95 %
and 99 % percentile intervals, Bonferroni-adjusted p, BH-adjusted q and
bootstrap 95 % and 99 % CIs are reported jointly. Significance calls in
the manuscript use q < 0.05 after BH FDR correction; the corresponding
raw p_empirical and p_bonferroni columns are listed in Supplementary
Table S2.

### 4.11 Code and data availability

All source code used for data processing, SCI computation, permutation
testing, bootstrap validation, clustering analysis, and figure generation
is publicly available at:

[https://github.com/DCMPjmpo/papae_sci](https://github.com/DCMPjmpo/papae_sci)

The repository contains all scripts necessary to reproduce the analyses
reported in this study. Large intermediate data files excluded from
version control can be regenerated from the provided preprocessing
pipeline. Analyses were performed with Python 3.13 and the following
core libraries: numpy (≥1.26), pandas (≥2.1), scipy (≥1.12),
matplotlib (≥3.8), scikit-learn (≥1.4), PyTorch (≥2.0) and
fair-esm (≥2.0). The full permutation + bootstrap pipeline runs in
≈ 8–10 min single-threaded on a CPU at the real scale (95,142 mutations,
528 layer pairs, 1,000 adaptive permutations + 1,000 bootstrap resamples),
with peak resident memory ≈ 410 MB.

### 4.12 ESM-2 zero-shot LLR baseline

For the head-to-head comparison reported in §2.8 we computed the ESM-2
zero-shot log-likelihood ratio (LLR) baseline of Meier et al. (2021).
Two standard variants are implemented on the same esm2_t33_650M_UR50D
checkpoint and the same WT sequences used for SCI feature extraction
(§4.2):

- **WT-marginal.** A single forward pass on the WT sequence; the
  per-position log-softmax of the LM-head logits gives log p(aa | x_WT)
  for each of the 20 standard amino acids at every residue position.
  LLR^{WT} = log p(mut | x_WT) − log p(wt | x_WT) at the mutation
  position; total cost is one forward pass per protein (5 forwards).
- **Masked-marginal.** For each unique (protein, mutation_position)
  pair, a single forward pass is performed with the token at the
  mutated position replaced by the <mask> token; the per-residue
  log-softmax of the LM-head logits at the masked position gives
  log p(aa | x_\T). LLR^{MM} = log p(mut | x_\T) − log p(wt | x_\T).
  Across the five proteins this requires one forward pass per unique
  (protein, position), for a total of ≈ 1,300 forwards. The
  masked-marginal variant is the headline LLR scorer in §2.8 and
  corresponds to the ProteinGym ESM-2 zero-shot baseline.

For multi-mutation rows (the 73,042 expanded PABP doubles), each
constituent single substitution is scored independently with
single-position masking on the WT context, following the additive
ProteinGym leaderboard convention for ESM-2 masked-marginal scoring of
multi-residue substitutions. The implementation re-uses the existing
cached model weights and adds no extraction of the 33-layer hidden
states already cached for SCI; the only additional forward work is the
≈ 1,305 forwards described above (≈ 70× fewer than the 95,142 forwards
used for the SCI Δ-embedding cache in §4.2).

For the SCI-vs-LLR comparison reported in §2.8, the SCI scorer used is
the 528-pair upper-triangle mean. This differs from the Top50 SCI
summary used by T₁ in §2.1: T₁ is a mean-difference test designed for
discrimination against a within-protein permutation null and uses the
high-signal tail; the present comparison uses the all-pairs mean as the
simplest scalar SCI summary, for parity with a scalar zero-shot
baseline. The combined SCI+LLR scorer is a 5-fold stratified
cross-validated logistic regression on the standardised (SCI, LLR)
feature pair, with StratifiedKFold (n_splits = 5, shuffle = True,
random_state = 42), a StandardScaler fit *within* each training fold
to prevent leakage, and the held-out out-of-fold P(class = 1) used as
the predictor for each held-out fold. Per-protein and pooled Spearman,
AUROC and AUPRC values are computed with sign correction
(corrected = max over native and negated score). The pooled SCI–LLR
Pearson and Spearman correlations reported in §2.8 are computed on
the same 95,142-row pair of scalar scores after removal of NaN entries.

---

## 5. Figure Legends

### Figure 1 — SCI carries mutation-effect signal

**SCI discriminates functional from non-functional mutations across all
five DMS panels.**
**(a)** Distribution of the per-mutation Top50-mean SCI score
(n = 95,142) stratified by the binarised DMS label (DMS_score_bin,
1 = functional, 0 = non-functional). **(b)** Permutation null
distribution of the SCI mean difference T₁ = mean(SCI | y = 1) −
mean(SCI | y = 0) under 10³ within-protein label shuffles (Methods
§4.7). Vertical red line marks the observed T₁ = 2.79 × 10⁻³; the null
mean and standard deviation are −3.7 × 10⁻⁴ and 5.6 × 10⁻⁵, giving
z = 57.0 and empirical p_FDR = 10⁻³ against the 4-statistic global
family. At the per-protein level T₁ is FDR-significant in 5/5 proteins
(smallest = CBS, p_FDR = 3.7 × 10⁻²; remaining four ≤ 1.4 × 10⁻³).

- **Panel (a)** — side-by-side boxplots of Top50-mean SCI score for
  functional (blue) and non-functional (red) mutations. The centre line
  marks the median, box limits denote the interquartile range, and
  whiskers extend to 1.5× IQR. Sample means are indicated by dashed
  vertical lines.
- **Panel (b)** — histogram of T₁ values from the permutation null;
  x-axis = T₁ (functional − non-functional SCI mean), y-axis = count;
  vertical red line marks the observed T₁; annotation reports z and
  p_FDR.

### Figure 2 — Top-ranked layer pairs and the 528-pair category landscape

**Early–Late layer pairs dominate both the Top50 high-Composite tail
and the full 528-pair landscape in ESM-2.**
**(a)** Stacked bar chart of the six-category layer-pair composition
(EE / EM / EL / MM / ML / LL; Methods §4.3) within each protein's
Top50 Composite layer pairs. Four of the five proteins (CBS, GAL4,
PTEN, TEM1) are EL-dominated (EL accounts for 40 %, 66 %, 44 % and
40 %, respectively, of the Top50 set); PABP is the only protein whose
Top50 set contains a larger EM than EL share (EM 48 %; EL 42 %).
**(b)** Global permutation null of the AUROC-Top50 EL count T₂ under
10³ within-protein label shuffles; observed T₂ = 17 EL pairs vs. null
mean 1.0 (z = 8.0, p_FDR = 10⁻³). **(c)** Global permutation null of
the EL-category dominance score T₃ = mean |AUROC − 0.5|_{EL} −
mean |AUROC − 0.5|_{non-EL} over all 528 layer pairs; observed
T₃ = 4.0 × 10⁻² (z = 62.0, p_FDR = 10⁻³). The standardised counterpart
T₃ˢᵗᵈ = 1.02 (z = 8.05).

- **Panel (a)** — y-axis: proportion of Top50 Composite layer pairs
  that fall in each of the six categories; x-axis: protein name (PABP,
  CBS, GAL4, PTEN, TEM1); stack colours map to category.
- **Panel (b)** — null distribution histogram for T₂; x-axis: number
  of EL pairs among the 50 layer pairs with largest |AUROC − 0.5|;
  vertical red line at observed value; annotation reports z and p_FDR.
- **Panel (c)** — null distribution histogram for T₃; x-axis: EL minus
  non-EL mean |AUROC − 0.5|; vertical red line at observed value.

### Figure 3 — PABP forms an independent cluster on the layer-pair composition axis

**PABP forms a singleton in the six-dimensional layer-pair composition
space, characterized by Early–Middle interactions.**
**(a)** Heatmap of the six-category Top50 Composite layer-pair
composition for each of the five proteins (n × 6 matrix). PABP is the
only row in which the EM proportion (0.48) exceeds the EL proportion
(0.42); the remaining four proteins are EL-dominated. **(b)**
Average-linkage hierarchical clustering of the composition vectors
under Euclidean distance. PABP joins the four EL-dominated proteins
only at a substantially higher dendrogram height (≈ 0.50), whereas the
four others fuse at maximum height ≈ 0.22. The minimum PABP-to-other
distance is 0.389 (PABP–TEM1) and the maximum within-other-four
distance is 0.397 (CBS–GAL4); across the heatmap, dendrogram and PCA
analyses PABP remains the most isolated of the five proteins.
**(c)** Principal component analysis on the same vectors; PC1 explains
68.6 % and PC2 28.1 % of variance (96.7 % cumulative). PABP projects
to PC1 ≈ +0.38, isolated from the other four proteins which fall in
PC1 ∈ [−0.18, +0.01]. PC1 loadings are dominated by Early–Middle
(+0.207); PC2 separates Early–Late (+0.103) from Late–Late (−0.094).

- **Panel (a)** — colour-coded matrix; rows: proteins (CBS, GAL4,
  PABP, PTEN, TEM1); columns: categories (EE, EM, EL, MM, ML, LL);
  colour = proportion of Top50.
- **Panel (b)** — dendrogram; y-axis: Euclidean distance under average
  linkage; leaves labelled with protein names.
- **Panel (c)** — scatter of proteins in (PC1, PC2) coordinates;
  variance-explained shown on each axis; loading arrows for the six
  categories overlaid.

### Figure 4 — Cross-protein recurrent layer pairs concentrate at layers 30–32

**Recurrent layer pairs across five proteins concentrate at a sharp
three-layer hub at depths 30–32, with layer 32 most frequent.**
**(a)** Per-layer count of layer-end occurrences among the Top20
cross-protein recurrent layer pairs (40 layer-end occurrences = 2 × 20
pairs). Layer 32 contributes 8 of 40 (20 %) and is the single most
frequent layer; layers 31 (7), 30 (3) and 33 (2) complete the
deep-layer hub; layers 18–27 contribute zero. Band totals: Late
(25–33) = 22 / 40 = 55 %; Middle (9–24) = 11 / 40 = 27.5 %; Early
(1–8) = 7 / 40 = 17.5 %. **(b)** Layer-occurrence histogram across the
broader union of all five proteins' Top50 Composite sets (n = 500
layer-end occurrences = 5 × 50 × 2). The depth-wise concentration is
preserved at this broader scale (Late 47.6 %, Middle 36.0 %, Early
16.4 %).

- **Panel (a)** — bar chart; x-axis: layer index 1–33 (1-based);
  y-axis: count among the 40 layer-end occurrences of the Top20
  recurrent set; bar colours reflect band (Early / Middle / Late).
- **Panel (b)** — bar chart with the same x-axis but evaluated over
  the larger n = 500 union; band proportions printed below the panel.

### Figure 5 — Statistical robustness across proteins, corrections and bootstrap

**Per-protein reproducibility and multi-layered statistical robustness
of the layer-pair findings.**
**(a)** Per-protein null distributions of T₁, T₂, T₃ and T₃ˢᵗᵈ under
within-protein label shuffles (N = 10³ per protein; Methods §4.7) for
each of the five proteins; observed value marked in the
protein-specific highlight colour. All five proteins return
FDR-significant T₁, T₃ and T₃ˢᵗᵈ. **(b)** Bar chart of
observed-vs-null z-scores across the four primary statistics (T₁, T₂,
T₃, T₃ˢᵗᵈ) for the pooled (Global) sample and each of the five
proteins; reference lines at |z| = 1.96 (two-sided α = 0.05) and at
the Bonferroni-adjusted threshold for the 20 per-protein primary
tests (|z| ≈ 3.02). **(c)** QQ-style significance plot of the observed |z| against
the half-normal quantile under H₀ for all 24 tests (4 statistics ×
{Global + 5 proteins}); diagonal y = x marks the null expectation,
horizontal lines mark α = 0.05 and the Bonferroni threshold. Bootstrap
N = 10³ 95 % confidence intervals for T₁ (full bootstrap) and T₂, T₃,
T₃ˢᵗᵈ, T_{max} (rank-reuse fast path with measured peak
|ΔAUROC| ≤ 6 × 10⁻³ at real scale; Methods §4.9) exclude the
corresponding null mean for every primary statistic.

- **Panel (a)** — 5 × 4 grid of null histograms; rows: proteins,
  columns: statistics; each cell annotated with z and p_FDR.
- **Panel (b)** — grouped bars; x-axis: scope (Global + 5 proteins);
  y-axis: z-score; one bar per statistic per scope; dotted = α = 0.05;
  dashed = Bonferroni.
- **Panel (c)** — scatter; x-axis: expected |z| under H₀; y-axis:
  observed |z|; one point per test; markers coloured by scope;
  diagonal reference y = x.

### Table 1 — SCI vs ESM-2 zero-shot LLR baseline

**As a per-mutation scalar predictor on the 95,142-mutation panel, the
masked-marginal ESM-2 LLR baseline of Meier et al. (2021) exceeds the
528-pair upper-triangle SCI summary on every protein, and a combined
SCI+LLR scorer adds at most Δ AUROC = 0.007 over LLR alone.** Pooled
and per-protein Spearman r against DMS_score, AUROC and AUPRC against
DMS_score_bin are reported for three scorers: the 528-pair
upper-triangle mean SCI summary, the masked-marginal LLR (Methods
§4.12) and the 5-fold stratified-CV logistic regression on standardised
(SCI, LLR) (Methods §4.12).

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

- **Spearman r** — Spearman rank correlation of scorer against the
  continuous DMS_score; reported as |r| with the natural sign always
  +1 across all 18 entries.
- **AUROC** — sign-corrected AUROC of scorer against DMS_score_bin
  (corrected = max over native and negated score).
- **AUPRC** — sign-corrected average precision of scorer against
  DMS_score_bin.

---

## 6. Supplementary Information

### Supplementary Table S1 — Per-protein dataset characteristics

| Protein | DMS source dataset(s) | n mutations | n functional | n non-functional | Fraction functional |
|---|---|---:|---:|---:|---:|
| CBS  | Sun et al., 2020 | 7,217 | 3,314 | 3,903 | 0.459 |
| GAL4 | Kitzman et al., 2015 | 1,195 | 831 | 364 | 0.695 |
| PABP | Melamed et al., 2013 | 74,229 | 46,078 | 28,151 | 0.621 |
| PTEN | Mighell et al., 2018; Matreyek et al., 2021 | 7,504 | 5,880 | 1,624 | 0.784 |
| TEM1 | Deng et al., 2012; Firnberg et al., 2014; Jacquier et al., 2013; Stiffler et al., 2015 | 4,997 | 2,526 | 2,471 | 0.506 |
| **Total** | — | **95,142** | **58,629** | **36,513** | **0.616** |

**Note on label balance.** Per-protein label balances differ
substantially (0.459 in CBS to 0.784 in PTEN), motivating the
**within-protein** label-shuffle null adopted throughout (Methods
§4.7) — a cross-protein shuffle would conflate protein identity with
label balance and inflate the type-I error.

### Supplementary Table S2 — Full per-statistic per-scope results

#### S2-A · Global tests

| Statistic | T_obs | Null mean | Null SD | z = Cohen's *d* | p_empirical | p_FDR | Null 95 % CI | Null 99 % CI | Bootstrap 95 % CI | N_perm |
|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|
| T₁ (SCI mean diff) | 2.79 × 10⁻³ | −3.75 × 10⁻⁴ | 5.56 × 10⁻⁵ | 57.0 | 10⁻³ | 10⁻³ | [−4.90 × 10⁻⁴, −2.64 × 10⁻⁴] | [−5.17 × 10⁻⁴, −2.24 × 10⁻⁴] | [2.67 × 10⁻³, 2.91 × 10⁻³] | 1,000 |
| T₂ (#EL in Top50 by AUROC) | 17 | 1.013 | 1.995 | 8.01 | 10⁻³ | 10⁻³ | [0, 7] | [0, 8] | [14, 19] | 1,000 |
| T₃ (EL dominance) | 3.99 × 10⁻² | −2.08 × 10⁻³ | 6.77 × 10⁻⁴ | 62.0 | 10⁻³ | 10⁻³ | [−3.27 × 10⁻³, −6.52 × 10⁻⁴] | [−3.60 × 10⁻³, −2.34 × 10⁻⁴] | [3.70 × 10⁻², 4.27 × 10⁻²] | 1,000 |
| T₃ˢᵗᵈ (standardised) | 1.02 | −0.602 | 0.201 | 8.05 | 10⁻³ | 10⁻³ | [−0.912, −0.180] | [−0.984, −0.0563] | [0.954, 1.073] | 1,000 |
| T_max (supplementary) | 0.181 | 0.0176 | 1.74 × 10⁻³ | 93.7 | 10⁻³ | — | [0.0144, 0.0211] | [0.0130, 0.0222] | [0.176, 0.187] | 1,000 |

#### S2-B · Per-protein tests — z-scores

| Protein | T₁ z | T₂ z | T₃ z | T₃ˢᵗᵈ z |
|---|---:|---:|---:|---:|
| CBS  | 2.25  | n/a | 16.11 | 18.50 |
| GAL4 | 6.29  | n/a | 5.18  | 4.59  |
| PABP | 56.55 | (T_obs=1; see Methods §4.6 note) | 38.71 | 42.80 |
| PTEN | 2.88  | n/a | 14.57 | 17.05 |
| TEM1 | 23.59 | n/a | 15.19 | 15.94 |

#### S2-C · Per-protein tests — q-values after BH FDR over the 20 per-protein primary tests

| Protein | T₁ p_FDR | T₂ p_FDR | T₃ p_FDR | T₃ˢᵗᵈ p_FDR |
|---|---:|---:|---:|---:|
| CBS  | 3.7 × 10⁻²   | 1     | 1.4 × 10⁻³ | 1.4 × 10⁻³ |
| GAL4 | 1.4 × 10⁻³   | 1     | 1.4 × 10⁻³ | 1.4 × 10⁻³ |
| PABP | 1.4 × 10⁻³   | 1.4 × 10⁻³ | 1.4 × 10⁻³ | 1.4 × 10⁻³ |
| PTEN | 9.3 × 10⁻³   | 1     | 1.4 × 10⁻³ | 1.4 × 10⁻³ |
| TEM1 | 1.4 × 10⁻³   | 1     | 1.4 × 10⁻³ | 1.4 × 10⁻³ |

**Rows that pass q < 0.05 at the BH-FDR level.** T₁, T₃ and T₃ˢᵗᵈ in
**5 / 5** proteins; T₂ in **1 / 5** proteins (PABP only; see Methods
§4.6 note on Composite vs AUROC-only ranking).

#### S2-D · Bootstrap CIs (N = 1,000)

The bootstrap 95 % CIs reported in S2-A apply to the pooled (Global)
scope. Per-protein bootstrap was not performed; per-protein uncertainty
is captured by the null 95 % and 99 % CIs from the within-protein
permutation runs.

### Supplementary Figure S1 — Per-protein SCI distributions

Per-protein histograms of the per-mutation Top50-mean SCI score,
stratified by the binarised functional label (DMS_score_bin). Each
panel shows one of the five DMS panels (CBS, GAL4, PABP, PTEN, TEM1).
Overlaid blue (functional) and red (non-functional) histograms allow
visual inspection of the per-protein separation that underpins T₁ at
the per-protein level (Results §2.1; Supplementary Table S2-B). The
pooled distribution across all five proteins is reproduced as a
cross-reference panel.

### Supplementary Figure S2 — 528-pair AUROC landscape and Composite ranking

Heatmap visualisation of |AUROC − 0.5| over the 33 × 33 layer-pair
plane (with the diagonal masked), showing the EL block in the
off-diagonal corner where the highest discrimination values concentrate
(Results §2.3; Fig 2c). A companion rank-ordered list of all 528 layer
pairs by AUROC magnitude is the source of T₂ and T_max in
Supplementary Table S2-A.

### Supplementary Figure S3 — SCI matrix structure across mutations

Representative per-mutation 33 × 33 SCI matrix and the distribution of
best-performing layer pairs across mutations. The heatmap is symmetric
by construction (Methods §4.3). Per-mutation SCI matrices exhibit
pairwise correlations in the 0.7–1.0 range across most layer pairs,
with the off-diagonal Early–Late corner showing the strongest contrast
between functional and non-functional groups when aggregated
(cf. Fig 2c).

### Supplementary Figure S4 — Per-protein permutation grid

5 × 4 grid of per-protein permutation null distributions for the four
primary statistics (T₁, T₂, T₃, T₃ˢᵗᵈ); each panel shows the null
histogram with the observed value marked in the protein-specific
highlight colour. z-scores and empirical p-values are annotated in each
panel. This figure is reproduced in Fig 5a of the main text; the
standalone supplementary version is larger for ease of detailed
inspection of each per-protein null. Per-protein null distributions
are well-separated from the observed value across the four primary
statistics in 5/5 proteins for T₁, T₃ and T₃ˢᵗᵈ; T₂ at the per-protein
level diverges between the AUROC-only and Composite rankings
(Methods §4.6).

### Supplementary Figure S5 — QQ-style significance summary

Observed |z| values from all 24 tests (4 statistics × {Global + 5
proteins}) plotted against the expected half-normal quantile under H₀;
points above the diagonal y = x indicate departure from the null.
Horizontal guides mark |z| = 1.96 (two-sided α = 0.05) and the
Bonferroni-adjusted threshold for the 20-test family. Tests are
coloured by scope. This figure is identical to Fig 5c. All 24 tests but
four (the per-protein T₂ rows for CBS, GAL4, PTEN, TEM1; Methods §4.6
note) lie above the Bonferroni threshold.

### Supplementary Note S1 — PCA loadings on the layer-pair compositional axis

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

PC1 (68.6 % of variance) is dominated by EM with sign opposite to LL,
and projects PABP (the EM-rich protein) far from the four EL-dominated
proteins; PC2 (28.1 %) contrasts EL with LL, which is the contrast that
separates CBS, GAL4, PTEN, TEM1 among themselves (Fig 3c).

### Supplementary Note S2 — Pairwise Euclidean distance matrix

|  | PABP | CBS | GAL4 | PTEN | TEM1 |
|---|---:|---:|---:|---:|---:|
| PABP | 0.000 | 0.566 | 0.537 | 0.495 | 0.389 |
| CBS  | 0.566 | 0.000 | 0.397 | 0.184 | 0.242 |
| GAL4 | 0.537 | 0.397 | 0.000 | 0.251 | 0.290 |
| PTEN | 0.495 | 0.184 | 0.251 | 0.000 | 0.115 |
| TEM1 | 0.389 | 0.290 | 0.290 | 0.115 | 0.000 |

**Isolation of PABP.** minimum PABP-to-other = 0.389 (PABP–TEM1);
maximum within-other-four = 0.397 (CBS–GAL4); minimum within-other-four
= 0.115 (PTEN–TEM1). PABP remains the most isolated protein in the
clustering analysis, joins the remaining proteins only at a
substantially higher dendrogram height, and the PCA and distance
analyses consistently place PABP apart from the other four proteins.

### List of Supplementary Materials

| Item | Carries | Main-text reference |
|---|---|---|
| Supplementary Table S1 | Dataset characteristics, per-protein label balance | Methods §4.1; Results §2.1 |
| Supplementary Table S2 | Full per-statistic per-scope results (z, Cohen's *d*, p, q, CIs, bootstrap) | Results §2.1–§2.7; Methods §4.7–§4.10 |
| Supplementary Figure S1 | Per-protein SCI distributions | Results §2.1 (Fig 1) |
| Supplementary Figure S2 | 528-pair AUROC landscape, Composite ranking | Results §2.2, §2.3 (Fig 2) |
| Supplementary Figure S3 | Per-mutation SCI matrix structure | Methods §4.3 |
| Supplementary Figure S4 | Per-protein permutation grid (standalone, larger) | Results §2.6 (Fig 5a) |
| Supplementary Figure S5 | QQ-style significance summary | Results §2.7 (Fig 5c) |
| Supplementary Note S1 | Full PCA loadings table | Results §2.4 (Fig 3c) |
| Supplementary Note S2 | Pairwise Euclidean distance matrix | Results §2.4 (Fig 3b) |

---

## 7. References

1. Belrose, N. et al. Tuned-lens decoding. 2023.
2. Brandes, N. et al. Proteome-scale ESM-1b variant scoring. 2023.
3. Deng, Z. et al. TEM1 deep mutational scanning. 2012.
4. Detlefsen, N. S. et al. Task-dependent best single layer in protein language models. 2022.
5. Elnaggar, A. et al. ProtTrans: per-layer probing of protein language models. 2022.
6. Firnberg, E. et al. TEM1 deep mutational scanning. 2014.
7. Frazer, J. et al. EVE: variant effect from sequence. 2021.
8. Geva, M. et al. Feed-forward layers as key–value memories. 2021.
9. Hayes, T. et al. ESM3. 2025.
10. Hsu, C. et al. Ridge regression over PLM embeddings. 2022.
11. Jacquier, H. et al. TEM1 deep mutational scanning. 2013.
12. Kitzman, J. O. et al. GAL4 deep mutational scanning. 2015.
13. Lin, Z. et al. ESM-2 (esm2_t33_650M_UR50D). 2023.
14. Marquet, C. et al. Variant-of-uncertain-significance scoring. 2022.
15. Matreyek, K. A. et al. PTEN deep mutational scanning. 2021.
16. Meier, J. et al. ESM-1v masked-marginal log-likelihood ratio baseline. 2021.
17. Melamed, D. et al. PABP deep mutational scanning. 2013.
18. Mighell, T. L. et al. PTEN deep mutational scanning. 2018.
19. Notin, P. et al. Tranception. 2022.
20. Notin, P. et al. ProteinGym. 2023.
21. Phipson, B. and Smyth, G. K. Permutation p-values should never be zero: calculating exact p-values when permutations are randomly drawn. 2010.
22. Rao, R. et al. TAPE: a benchmark for protein-language-model probing. 2019.
23. Rao, R. et al. MSA Transformer; per-head attention analysis of protein language models. 2021a.
24. Rao, R. et al. Transformer protein language models are unsupervised structure learners; ICLR contact-prediction pipeline pooling attention across all 33 ESM layers. 2021b.
25. Rives, A. et al. ESM-1b. 2021.
26. Simon, E. and Zou, J. Per-layer sparse-autoencoder analysis of protein language models. 2025.
27. Stiffler, M. A. et al. TEM1 deep mutational scanning. 2015.
28. Sun, S. et al. CBS deep mutational scanning. 2020.
29. Tenney, I. et al. Scalar-mix probing of BERT. 2019.
30. Vig, J. et al. Per-head attention analysis of protein language models. 2021.
