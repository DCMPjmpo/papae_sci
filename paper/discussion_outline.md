# Discussion outline

> Cross-layer interaction analysis in ESM-2 for mutation-effect encoding.
>
> Language convention: use *observation*, *suggest*, *indicate*, *consistent with*, *hypothesis*. Avoid causal verbs (*causes*, *drives*, *outperforms*, *proves*).
>
> Findings referenced below (F1–F5):
>
> - **F1.** Early–Late (EL) pairs show significantly higher AUROC and composite scores than EE, EM, MM and LL across all 528 layer pairs.
> - **F2.** Four proteins (CBS, GAL4, PTEN, TEM1) are dominated by Early–Late interactions in their Top50 Composite layer pairs.
> - **F3.** PABP forms a distinct cluster (hierarchical clustering, PCA) and is dominated by Early–Middle interactions.
> - **F4.** Recurrent cross-protein layer pairs overwhelmingly involve layers 30–32 (Late band: 55% of recurrent occurrences).
> - **F5.** Layer 32 is the most recurrent individual layer (8 / 40 occurrences in the Top20 recurrent pairs).

---

## D1 — Relationship to ESM layer hierarchy literature

**Thesis.** Our layer-pair geometry is *consistent with* the previously reported layer hierarchy of protein language models, but indicates that mutation-effect signal is not localized to any single layer of that hierarchy.

- **Anchor in prior layer-wise findings.** Earlier per-layer probing studies of ESM and related PLMs report that biophysical and structural properties emerge at characteristic depths (Rives et al., 2021; Lin et al., 2023; Elnaggar et al., 2022; Vig et al., 2021). Our observation that **layers 30–32 are the dominant late-end anchor (F4, F5)** is consistent with this body of work, which repeatedly identifies the deepest layers of ESM as the locus of higher-order semantic / functional features.
- **Where our observation departs from the "best layer" framing.** Detlefsen, Hauberg & Boomsma (2022) note that the "best layer" of a PLM depends on the downstream task. Our results suggest a further refinement: for mutation-effect signal, *no single layer is privileged*; rather, the strongest discriminative signal is recovered from **pairs of layers** that span the network (F1). This indicates that the canonical "use the last hidden state" or "use the masked-LM head" protocol (e.g. Meier et al., 2021; Brandes et al., 2023) may capture only a marginal projection of the underlying encoding.
- **Hierarchy persists, but in a 2D form.** The Top20 recurrent pairs occupy a bipartite structure — a sharp late-layer hub at 30–32 paired with broader donors in the Early-to-mid-Middle range (F4) — which is consistent with the hierarchical organization documented for both PLMs and general transformers (Tenney et al., 2019; Rogers et al., 2020). Our observation can therefore be read as extending the **one-dimensional layer hierarchy** to a **two-dimensional layer-pair hierarchy**.
- **Limitation.** All observations here are made on a single backbone (ESM-2, 33 layers). Whether the same Early–Late dominance and the Layer 32 anchor generalize to ESM-2 (Lin et al., 2023), ESM3 (Hayes et al., 2025), or non-Transformer PLMs remains an open question.

## D2 — Cross-layer integration hypothesis

**Thesis.** The pattern of cross-protein recurrence is consistent with a working hypothesis that **mutation-effect information is integrated across layers** rather than read off any single layer; we frame this only as a hypothesis suggested by the layer-pair geometry, not as a causal mechanism.

- **The observation that suggests integration.** The EL-dominance across 528 layer pairs (F1) and the absence of any high-scoring Early–Early, Middle–Middle or Late–Late pair within the Top20 recurrent set (F4) together indicate that the discriminative signal is carried *between* depth bands rather than within them. This is consistent with prior mechanistic descriptions of transformer residual streams as cross-layer integrators (Geva et al., 2021; Belrose et al., 2023; Tenney et al., 2019), which characterize the depth axis as a substrate over which intermediate representations are progressively combined.
- **Plausible non-causal reading.** Under this hypothesis, **early layers may contribute residue-local / biochemical features**, while **late layers (notably 30–32) may contribute family-level / contextual features**, and the discriminative signal for mutation effects is then *the relationship between* these two depth-stratified views (F1, F4). We emphasize the conditional language: our data describe a correlation between layer-pair composition and mutation-effect discrimination; they do not isolate a causal contribution of any individual layer.
- **What this is *not*.** This hypothesis is distinct from claims that "deeper PLMs predict mutations better" (Meier et al., 2021; Notin et al., 2022; Brandes et al., 2023) — those works rely on a single chosen layer. The closest existing PLM analysis that explicitly fuses across layers (Rao et al., 2021, ICLR) does so for residue *contacts*. Our observation suggests that **mutation-effect signal may demand a multi-layer combination of similar form but different content** to that used for contact recovery.
- **Falsifiable predictions** (for follow-up work, not claims of this paper).
  1. Ablation or feature-flow interventions targeting **only layer 32** should leave residual mutation-effect signal recoverable from earlier layers, if integration is real.
  2. A linear probe restricted to a single layer band (Early-only, Middle-only, or Late-only) should systematically underperform a probe over an Early–Late pair, *at the same dimensionality*. Without these tests, the integration claim remains a hypothesis suggested by — but not established by — our observations.
- **Limitation.** The Composite score combines AUROC, Spearman *r* and Cohen's *d* and could in principle reward layer pairs that *individually* contain weak but complementary signals. We have not yet decomposed how much of the EL advantage reflects *interaction* versus *aggregation*; this distinction is essential for the integration hypothesis and is reserved for follow-up work.

## D3 — Protein-specific heterogeneity

**Thesis.** The five proteins do not share a single layer-pair signature; PABP's distinct composition (F3) indicates **protein-specific encoding heterogeneity** rather than a deficiency of any single protein or of the model.

- **The observation.** Hierarchical clustering and PCA over the layer-pair category vector place CBS, GAL4, PTEN and TEM1 in one compact group (pairwise Euclidean distances ≤ 0.40) and **PABP at a singleton position**, separated from its nearest neighbour (TEM1) by 0.39 and from the centroid of the four-protein cluster by ~0.5 (F3). PABP's signature is the only one in which **Early–Middle exceeds Early–Late** (48% vs 42%) and the only one with no Late–Late contribution.
- **Consistent reading 1 — biophysical class.** PABP is an RNA-binding protein, whereas CBS / GAL4 / PTEN / TEM1 span enzymes (CBS, PTEN, TEM1) and a DNA-binding transcription factor (GAL4). The observation that PABP's signal is shifted *toward the middle layers* is consistent with the depth-stratified interpretability findings of Vig et al. (2021), which indicate that different biophysical properties are tracked at characteristic depths of protein Transformers. We *do not* claim that "RNA binding causes Middle-layer encoding"; we only note that the heterogeneity we observe **co-varies with a known biological partition** that has previously been associated with depth-dependent representations.
- **Consistent reading 2 — DMS experimental design.** PABP's DMS assay (Melamed et al., 2013) measures growth-coupled binding, whereas the other four are stability/activity assays in microbial or human cell contexts. Differences in fitness landscape topology — i.e. what each DMS assay rewards — could plausibly project onto different layer-pair compositions. This is a non-causal, design-level confound that our analysis cannot rule out; it should be acknowledged explicitly in the paper.
- **What this is NOT.** The PABP singleton is **not** a sign of poor model fit, low data quality, or PLM failure on RNA-binding proteins. ESM-2's Top50 Composite layer pairs are well-separated in score space for PABP as for the other proteins; only their *layer-band composition* differs.
- **Limitation and open question.** With *n* = 5 proteins, PABP's apparent singleton status is a *qualitative observation*, not a statistical generalization. Whether RNA-binding proteins (or, more broadly, proteins with shallow fitness landscapes) systematically display Early–Middle dominance is a hypothesis that requires expansion to ≥ 20 proteins spanning multiple biophysical classes.

## D4 — Implications for PLM interpretability

**Thesis.** Our findings, taken together, suggest that **layer-pair geometry is a useful interpretability primitive for PLMs** — one that complements, rather than replaces, per-layer probing and per-head attention analysis.

- **Adds a 2D axis to existing probing pipelines.** Existing PLM interpretability methods are organized along a single depth axis: layer-by-layer probes (Elnaggar et al., 2022; Detlefsen et al., 2022), per-layer attention (Vig et al., 2021), or per-layer sparse autoencoder features (Simon & Zou, 2025). Our observation suggests these one-axis methods may underrepresent signal that lives in the **off-diagonal** of the layer-pair plane, particularly between Early and Late bands (F1, F4).
- **Reframes the "final layer" convention.** A large body of work computes mutation-effect scores from a single masked-LM head or a single hidden layer (Meier et al., 2021; Brandes et al., 2023; Hsu et al., 2022; Marquet et al., 2022; Frazer et al., 2021). Our observation that Layer 32 is the dominant late-end anchor but appears *paired* with broader donor layers (F4, F5) indicates that this convention captures the *late endpoint* but not the *cross-layer relationship* in which the discriminative signal sits. This reframing is interpretive, not predictive: we make no claim that layer-pair features yield better mutation-effect prediction in absolute terms.
- **Connects PLM interpretability to general transformer mechanistic interpretability.** The notion that information is built up via the residual stream — encoded layer-by-layer and integrated across layers — is a central organizing principle of mechanistic interpretability for language models (Geva et al., 2021; Belrose et al., 2023; Tenney et al., 2019). Our observation that mutation-effect signal recurs at cross-depth pairs is *consistent with* this picture and suggests that **the same residual-stream-centered toolkit could be productively applied to protein language models** in future work, including logit/tuned-lens probes, activation patching across layer pairs, and sparse-feature steering studies.
- **Practical use as a model-comparison probe.** The layer-pair composition vector (six categories × $n$ proteins) can be computed for any PLM with a few dozen layers and used as a low-dimensional fingerprint to compare how different PLMs (ESM-2 vs ESM-2 vs ESM3 vs ProtT5) organize mutation-effect signal across their depth. We propose this as a hypothesis-generating tool rather than as a benchmark.
- **Boundary of these claims.** We do not claim that:
  - layer-pair features yield higher predictive AUROC than tuned PLM scores (we did not test this);
  - the Early–Late pattern is causally responsible for ESM-2's mutation-effect performance (we provide correlational geometry, not interventions);
  - Layer 32's prominence generalizes to other PLMs (we tested only ESM-2).

  These boundaries should be stated explicitly in the manuscript so that subsequent work can target them.

---

## Cross-cutting limitations to surface in the final paragraph

- *n* = 5 proteins, all from existing DMS datasets — narrow sampling of the protein universe.
- Single PLM backbone (ESM-2, 33 layers).
- Top50 / Top20 thresholds are conventional cut-offs; we report robustness across thresholds in Supplementary but the precise numeric proportions are threshold-dependent.
- All claims are correlational. No intervention experiments (ablation, activation patching, layer freezing) were performed; integration remains a hypothesis suggested by the geometry.
