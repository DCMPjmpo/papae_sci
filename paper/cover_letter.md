# Cover Letter

[Date]

Dear Editors,

We are pleased to submit our manuscript, **"Early–Late layer-pair
interactions encode mutation-effect signal in the protein language
model ESM-2 across five deep mutational scans,"** for consideration
as an Article in *Nature Communications*.

Protein language models (PLMs) have become a default tool for
estimating mutation effects from sequence, yet the field's
interpretive vocabulary remains one-dimensional. Variant scores are
read from a single representation — typically the masked-LM head or
the final hidden state — and probing studies map biophysical
properties onto individual layers one at a time. How
mutation-effect information is organised *across* the depth of a
PLM has not been systematically characterised.

Our manuscript reports the first systematic mapping of cross-layer
organisation of mutation-effect signal in a protein language model.
For each of 95,142 single amino-acid substitutions across five
deep-mutational-scanning panels (CBS, GAL4, PABP, PTEN, TEM1; from
the ProteinGym substitution release), we computed a per-pair
discrimination metric — Structural Context Inconsistency (SCI) —
over every one of the 528 unordered layer pairs of ESM-2
(esm2_t33_650M_UR50D). Three reproducible patterns emerge on this
528-pair plane:

1. **Early–Late dominance.** EL pairs concentrate the
   mutation-effect signal at both the high-Composite tail and the
   full 528-pair category-level comparison (T₃ z = 62 against a
   within-protein label-shuffle null; q = 10⁻³).
2. **A sharp three-layer recurrence hub at depths 30–32.** Layer 32
   alone supplies 20 % of the layer-end occurrences in the Top20
   cross-protein recurrent set, while layers 18–27 contribute zero.
3. **Per-protein heterogeneity.** Four of the five proteins share an
   EL-dominated layer-pair signature, but PABP — the only
   RNA-binding protein in the panel — forms a singleton characterized
   by Early–Middle interactions on hierarchical clustering, PCA, and
   pairwise Euclidean distances.

Every primary effect survives Bonferroni FWER, Benjamini–Hochberg
FDR over twenty per-protein primary tests, and a non-parametric
bootstrap (95 % CI for each statistic excludes the null mean).

We have included an explicit head-to-head comparison against the
field's de facto zero-shot baseline, the ESM-2 masked-marginal
log-likelihood ratio (LLR) of Meier et al., on the same 95,142
mutations. LLR is the stronger univariate scalar predictor
(pooled AUROC = 0.753 vs SCI 0.669), and a cross-validated SCI+LLR
combination reaches AUROC = 0.758 — within 0.005 of LLR alone. SCI
and LLR are only moderately correlated (Pearson r = 0.456;
Spearman ρ = 0.467), consistent with our interpretation that SCI
captures a representation-level signal that is partially independent
of likelihood-based scoring. We make this comparison transparently:
the contribution of this work is *not* a new scalar predictor but a
layer-resolved interpretability primitive — the 528-pair landscape,
the layer-30–32 recurrence hub, and the protein-specific
compositional fingerprint — none of which is exposed by a single
LM-head likelihood ratio.

We believe this work is of broad interest to the *Nature
Communications* readership for three reasons. First, it introduces a
two-axis (layer × layer) vocabulary to a PLM interpretability
literature that has been organised almost exclusively along a single
depth axis. Second, the layer-pair composition vector defines a
low-dimensional fingerprint that can be computed for ESM-1b
(Rives et al., 2021), ESM3 (Hayes et al., 2025), ProtT5
(Elnaggar et al., 2022) and other backbones, providing a route to
test whether the EL block and the layer-32 anchor reported here are
properties of ESM-2 specifically or of large protein language models
more broadly. Third, the result connects naturally to the
residual-stream–centred interpretability picture developed for
general Transformers (Tenney et al., 2019; Geva et al., 2021;
Belrose et al., 2023), opening protein language models to the same
class of mechanistic tools.

This manuscript has not been published elsewhere and is not under
consideration at another journal. All authors have approved the
submission. We declare no competing financial interests. All data
and analysis code are made available through the project repository
listed in the Data and Code Availability statements (see
accompanying `data_availability.md` and `code_availability.md`).

We suggest the following potential reviewers, none of whom have a
conflict of interest with the authors:

- [Suggested reviewer 1 — expertise: PLM interpretability]
- [Suggested reviewer 2 — expertise: deep mutational scanning / variant effect]
- [Suggested reviewer 3 — expertise: Transformer mechanistic interpretability]
- [Suggested reviewer 4 — expertise: ESM family of protein language models]

We thank the editors for considering our work and we look forward to
your assessment.

Sincerely,

[Corresponding author name]
[Affiliation]
[Email]
[Phone]
[ORCID]

On behalf of all co-authors.
