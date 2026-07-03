# Introduction (Phase E)

> Positioning: this work is presented as the **first systematic mapping
> of cross-layer organisation of mutation-effect signal in a protein
> language model** — a *Layer-Pair Geometry Discovery* paper rather
> than a method-introduction paper. SCI is treated as one specific
> operationalisation of "per-pair discriminative signal" and is *not*
> the focal contribution.
>
> Target length ≈ 800 words. Citations refer to the literature matrix
> at `paper/literature_matrix/literature_matrix.md` and the mapping at
> `paper/literature_mapping.md`. All claims are restricted to the
> validated C1–C7 inventory in `paper/scientific_story_map.md`.

---

## Introduction

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
than one layer: the contact-prediction pipeline of Rao et al. (2021,
ICLR), in which a sparse logistic regression pools attention from
all 33 ESM layers. Critically, that work targets *residue contacts*,
not mutation effects, leaving the cross-layer geometry of
mutation-effect signal in PLMs an empty cell in the existing
literature (see also `paper/literature_mapping.md`).

**¶4 — A layer-pair geometry view of mutation-effect encoding.**
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
(SCI; Methods §3.3), is one specific operationalisation of
"discriminative signal at a layer pair" and is incidental to the
central finding; what we report is the geometry of the layer-pair
plane itself.

**¶5 — Three discoveries, all validated against a within-protein
permutation null.**
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
PABP forms an independent singleton driven by Early–Middle
interactions on the same six-band axis (hierarchical clustering,
PCA, and pairwise Euclidean distances). Every primary effect
survives Bonferroni FWER, Benjamini–Hochberg FDR over twenty
per-protein primary tests, and non-parametric bootstrap.

**¶6 — Contribution and positioning.**
This work is, to our knowledge, the first systematic mapping of how
mutation-effect signal is organised across the layers of a protein
language model. The contribution is interpretive rather than
predictive: we make no claim that layer-pair features improve
mutation-effect prediction over tuned single-layer scores, and we
perform no head-to-head benchmark. Rather, we propose layer-pair
geometry as an interpretability primitive that complements
single-layer probing (Elnaggar et al., 2022; Detlefsen et al., 2022;
Vig et al., 2021), per-head attention analysis (Rao et al., 2021),
and per-layer sparse-autoencoder analysis (Simon & Zou, 2025); and
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

## Notes on positioning (for collaborators)

- **"First systematic"** is justified by: (i) no prior PLM study
  computes a per-pair discrimination map across all 528 unordered
  layer pairs; (ii) no prior PLM study reports cross-protein
  recurrence statistics at the layer-pair level; (iii) the only
  near-neighbour (Rao et al., 2021, ICLR) is for *contacts* not
  mutation effects (`paper/literature_mapping.md` §6). The phrase
  *"first systematic mapping of cross-layer organisation of
  mutation-effect signal in a protein language model"* is therefore
  defensible against the literature we have audited.
- **SCI is de-emphasised in this Introduction** — it is mentioned only
  once, in ¶4, as the operationalisation we adopt; the Results
  refer to it neutrally as "the per-pair discrimination metric" and
  the Discussion does not depend on SCI's particular form. Should
  reviewers later ask for an SCI-alternative robustness analysis, that
  would be future work; the Layer-Pair Geometry Discovery framing is
  independent of which specific discrimination score is used.
- **No causal claims are introduced** (`prove`, `demonstrate`,
  `cause`, `drive`, `mechanism`, `outperform` do not appear). The
  three findings are framed throughout as observations on the
  layer-pair plane, not as mechanistic statements about ESM-2's
  internal computation. The cross-layer integration *hypothesis* is
  reserved for Discussion (D2 in `paper/discussion_full.md`).
- **References to ESM-1b** in this Introduction are limited to the
  historical Brandes et al. 2023 work (which used ESM-1b at proteome
  scale) and to the cross-backbone replication target list in ¶6.
  Our model is ESM-2 (esm2_t33_650M_UR50D), cited as Lin et al. 2023.
