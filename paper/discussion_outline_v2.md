# Discussion outline (v2)

> Cross-layer interaction analysis in ESM-2 for mutation-effect encoding.

## Convention

- **Language.** Use *observation*, *suggest*, *indicate*, *consistent with*, *hypothesis*, *may*, *plausibly*, *co-varies with*. Do **not** use *prove*, *demonstrate*, *show mechanism*, *cause*, *drive*, *outperform*.
- **Tagging.** Every clause that makes a claim must carry an `(F#)` tag identifying the underlying finding. Clauses that synthesize multiple findings carry compound tags `(F#, F#)`.
- **Layout.** Each section is split into four blocks: **Thesis → Observations (factual) → Hypotheses / consistent readings (interpretive) → Limitations → Future directions.** This separation is intentional and should be preserved in the prose draft.

## Validated findings

- **F1.** SCI carries mutation-effect signal.
- **F2.** Top-ranked layer pairs are enriched for Early–Late (EL) interactions.
- **F3.** Across all 528 layer pairs, EL pairs significantly outperform EE, EM, MM and LL categories (p ≪ 0.001).
- **F4.** Four proteins (CBS, GAL4, PTEN, TEM1) are EL-dominated, whereas PABP is EM-shifted.
- **F5.** Cross-protein recurrent layer pairs concentrate at layers 30–32, especially layer 32.

---

## D1 — Deep-layer enrichment in ESM

### Thesis
The late-end anchor of mutation-effect signal in ESM-2 is sharply concentrated at layers 30–32 and is most reliably represented by layer 32, an observation consistent with prior reports of higher-order property encoding at the deepest ESM layers (F5).

### Observations
- Among Top20 cross-protein recurrent layer pairs, the Late band (25–33) accounts for 55% of all layer-end occurrences (F5).
- Layer 32 is the single most recurrent layer (8 / 40 occurrences in the Top20 recurrent set), followed by layers 31 (7) and 30 (3) (F5).
- This deep-layer enrichment is observed in the same recurrent set that also exhibits EL dominance (F2, F5) and is reproducible across the 528-pair analysis (F3).
- The deep-layer enrichment is observed on a representational signal that has already been validated to carry mutation-effect discriminative information (F1).

### Hypotheses / consistent readings (no causal claim)
- The concentration of recurrence at layers 30–32 is *consistent with* prior layer-wise probing studies that identify late ESM layers as the locus of higher-order semantic and functional features (Rives et al., 2021; Lin et al., 2023; Vig et al., 2021; Elnaggar et al., 2022) (F5).
- It is also *consistent with* the broader observation that "best layer" of a PLM depends on the task (Detlefsen et al., 2022); within the mutation-effect task, the deep-layer end is sharp, but it appears in our data *as one end of a pair* rather than as a stand-alone signal locus (F2, F5).
- The fact that the deep-layer hub is narrow (3 layers) but recurrent across 5 proteins *suggests* — but does not establish — that the deep-layer anchor may be a shared property of how ESM-2 organizes high-level evolutionary / functional context (F5).

### Limitations
- Single backbone (ESM-2, 33 layers); the layer 32 anchor cannot be assumed to transfer to ESM-2, ESM3 or ProtT5.
- Recurrence is computed on Top20 cross-protein pairs over 5 proteins, so "concentration at 30–32" is a finite-sample observation (F5), not an asymptotic claim.
- The Late band (25–33) spans 9 layers but recurrence is essentially absent in layers 25–29; the proportion 55% therefore reflects an extremely localized hub, not uniform Late enrichment.

### Future directions
- Replicate the recurrence analysis on ESM-2 and ESM3 to test whether a single late-end anchor layer (analogous to layer 32) is a generic property of large PLMs.
- Stratify by ESM-2 training checkpoints to ask at what stage of pretraining the layer 32 anchor emerges.

---

## D2 — Cross-layer integration hypothesis

### Thesis
The systematic, cross-protein EL enrichment is *consistent with* a working hypothesis that mutation-effect information in ESM-2 is integrated across layers rather than localized to a single layer (F2, F3). We frame this only as a hypothesis suggested by the layer-pair geometry, not as a causal mechanism.

### Observations
- Across all 528 unordered layer pairs, EL pairs achieve significantly higher AUROC and Composite scores than EE, EM, MM and LL pairs (p ≪ 0.001) (F3).
- The Top-ranked subset of layer pairs is enriched for EL interactions, replicating the global pattern at the high-signal tail (F2).
- The recurrent set (F5) is dominated by pairs that connect a deep-layer endpoint (typically in 30–32) to a non-Late endpoint, consistent with cross-depth, rather than within-band, integration (F2, F5).
- These statistical observations are based on a representational signal (SCI) whose ability to discriminate mutation effects has already been validated (F1).

### Hypotheses / consistent readings (no causal claim)
- Under this working hypothesis, **early layers may contribute residue-local biochemical features** and **late layers (notably 30–32) may contribute family-level / contextual features**, and the mutation-effect signal *may correspond to the relationship between* these two depth-stratified representations (F2, F3, F5).
- This reading is *consistent with* mechanistic interpretability work on general transformers that characterizes the residual stream as a substrate for cross-layer integration (Geva et al., 2021; Belrose et al., 2023; Tenney et al., 2019).
- The hypothesis is *distinct from* the conventional protocol of computing mutation-effect scores from a single chosen layer (Meier et al., 2021; Brandes et al., 2023; Hsu et al., 2022) and *distinct from* the only existing PLM analysis that explicitly fuses across layers (Rao et al., 2021), which targets residue contacts rather than mutation effects (F2, F3).

### Limitations
- All findings are correlational. We have not performed activation patching, single-layer ablation, or any intervention that would isolate a causal role for any layer (F2, F3, F5).
- The Composite score combines AUROC, Spearman *r* and Cohen's *d*. In principle, EL pairs could be ranked higher because their two endpoints contain *independently weak but complementary* signals (aggregation) rather than because the signal is genuinely *interactive*. Our current analysis does not separate these two cases (F3).
- "Cross-layer integration" as defined here is a property of the *layer-pair representation*, not of any internal computation; we do not claim ESM-2's forward pass implements integration in any particular way (F2).

### Future directions
- A linear probe restricted to a single layer band (Early-only, Middle-only, Late-only), evaluated at the same dimensionality as an Early–Late paired probe, would test whether EL pairs carry signal *beyond what either band carries alone* (test for F2, F3).
- Activation patching of layers 30–32 should leave residual mutation-effect signal recoverable from earlier layers, if the integration hypothesis holds (test for F5).
- A residual-stream / tuned-lens analysis (Belrose et al., 2023) applied to ESM-2 would map at what depth mutation-effect information first becomes linearly recoverable, and how it evolves layer-by-layer.

---

## D3 — Protein-specific heterogeneity

### Thesis
The five proteins do not share a single layer-pair signature; PABP's EM-shifted composition relative to the four EL-dominated proteins (F4) is best read as an observation of protein-specific encoding heterogeneity rather than as evidence of model failure on PABP.

### Observations
- Four of the five proteins (CBS, GAL4, PTEN, TEM1) are EL-dominated in their Top50 Composite layer pairs (F4).
- PABP is the only protein in which Early–Middle exceeds Early–Late within the Top50 set, and the only one with no Late–Late contribution (F4).
- Hierarchical clustering and PCA over the layer-pair category vector place PABP as a singleton, separated from the centroid of the remaining four proteins by a distance comparable to the entire pairwise spread within that four-protein cluster (F4).
- All five proteins exhibit mutation-effect signal under SCI (F1); the heterogeneity is therefore observed at the level of *which layer pairs* carry the signal, not whether signal is present (F1, F4).
- The cross-protein recurrent set (F5) — which concentrates at layers 30–32 — reflects the four EL-dominated proteins more than it reflects PABP, consistent with PABP's distinct composition (F4, F5).

### Hypotheses / consistent readings (no causal claim)
- The PABP shift toward Early–Middle interactions *co-varies with* biophysical class: PABP is an RNA-binding protein, whereas the four EL-dominated proteins span enzymes (CBS, PTEN, TEM1) and a transcription factor (GAL4). This co-variation is *consistent with* depth-stratified attention-head specialization in protein Transformers (Vig et al., 2021) (F4).
- The shift is also *consistent with* differences in DMS assay design: PABP's DMS measures growth-coupled binding (Melamed et al., 2013), whereas the four others measure stability or enzymatic activity in microbial or human cell contexts. Differences in the shape of the underlying fitness landscape may project onto different layer-pair compositions (F4).
- We do **not** claim that RNA binding or assay design *causes* Middle-layer involvement; we only note that the observed heterogeneity is non-random and aligns with two known biological axes that may inform follow-up sampling (F4).

### What this is NOT
- The PABP singleton is **not** an indicator of poor model fit, low data quality, or PLM failure on RNA-binding proteins. PABP's Top50 Composite scores are well-separated from baseline as for the other proteins (F1); only their *layer-band composition* differs (F4).

### Limitations
- *n* = 5 proteins; PABP's singleton status is a *qualitative* observation, not a statistical generalization (F4).
- Biophysical class and DMS assay design are confounded in this sample (one RNA-binder + four non-RNA-binders, each with its own assay), so the two consistent readings above cannot be disentangled (F4).

### Future directions
- Expand the protein panel to ≥ 20 sequences spanning multiple biophysical classes (enzymes, RNA-binders, scaffold proteins, signaling proteins) and matched assay types, then test whether biophysical class predicts EM vs EL dominance (test for F4).
- Within RNA-binding proteins specifically, compare PABP to other RBPs with available DMS data to ask whether the EM shift is family-specific or a single-protein effect (test for F4).

---

## D4 — Implications for PLM interpretability

### Thesis
The systematic EL enrichment (F2, F3), the narrow deep-layer anchor (F5), and the protein-specific compositional heterogeneity (F4) together *suggest* that layer-pair geometry is a useful interpretability primitive for PLMs — one that complements, rather than replaces, single-layer probing and per-head attention analysis (F1).

### Observations
- Mutation-effect signal carried by SCI is non-trivially structured along the layer-pair axis: the strongest signal is concentrated in the EL off-diagonal block of the 528-pair plane (F1, F3).
- Per-layer probing protocols, by construction, cannot exhibit this off-diagonal structure; they integrate over one axis of a two-axis object (F1, F2, F3).
- The deep-layer anchor at layers 30–32 (F5) appears in *paired* form rather than as a stand-alone signal locus (F2, F5).
- Compositional heterogeneity across proteins (F4) means a single layer-pair fingerprint cannot summarize a PLM's mutation-effect representation; the fingerprint must be reported per protein.

### Hypotheses / consistent readings (no causal claim)
- The layer-pair vector (six categories × *n* proteins) *may* serve as a low-dimensional, model-comparable fingerprint of how a PLM organizes mutation-effect signal across its depth, and could in principle be computed for ESM-2, ESM3, ProtT5 and other backbones (F2, F4).
- The observation that mutation-effect signal recurs at cross-depth pairs is *consistent with* residual-stream-centered mechanistic interpretability for general transformers (Geva et al., 2021; Belrose et al., 2023; Tenney et al., 2019), and *suggests* that the same toolkit — logit/tuned-lens probes, activation patching across layer pairs, sparse-feature steering (Simon & Zou, 2025) — could be productively applied to PLMs in future work (F1, F2, F3).
- A specific working hypothesis (to be tested, not claimed): per-layer sparse-autoencoder features (Simon & Zou, 2025) defined on early and late ESM-2 layers may be linearly *combinable* into a mutation-effect classifier with higher recovery than features from either band alone (F2, F3).

### Boundary of these claims (write into the manuscript as a paragraph)
- We do not claim that layer-pair features yield higher predictive AUROC than tuned single-layer PLM scores (we did not perform a head-to-head comparison) (F1, F2, F3).
- We do not claim that the EL pattern is causally responsible for ESM-2's mutation-effect performance; we report correlational geometry, not interventions (F2, F3).
- We do not claim that layer 32's prominence generalizes to other PLMs; we have tested only ESM-2 (F5).

### Limitations
- The "interpretability primitive" framing is conceptual; we have not benchmarked layer-pair features as a prediction substrate.
- The off-diagonal structure has been characterized via six-category compositions (EE/EM/EL/MM/ML/LL); finer-grained partitions (e.g. 9 or 11 bands) may reveal additional structure that the 6-band view collapses.

### Future directions
- Cross-model fingerprinting: compute the layer-pair composition vector for ESM-2 (33 / 36 layers) and ESM3 and compare to the ESM-2 fingerprint developed here (test for F4 generalization).
- Per-task generalization: compute layer-pair compositions for downstream tasks beyond mutation-effect discrimination (e.g. thermostability, secondary-structure recovery) to ask whether each task has its own off-diagonal signature (test for F2, F3).
- Connect SCI-based layer-pair analysis to sparse-autoencoder features (Simon & Zou, 2025) by asking which per-layer SAE features differentially activate within the Top50 EL pairs vs the Top50 EL pairs of a non-mutation control task.

---

## Cross-cutting limitations (collect in a closing limitation paragraph)

- *n* = 5 proteins; all from existing DMS datasets — narrow sampling of the protein universe (F4).
- Single PLM backbone (ESM-2, 33 layers); no replication on ESM-2, ESM3, ProtT5 within this manuscript (F2, F3, F5).
- Top50 / Top20 thresholds are conventional cut-offs; numeric proportions are threshold-dependent (F2, F4, F5).
- All claims are correlational — no interventional / mechanistic experiments (no ablation, no activation patching, no layer freezing). Integration (D2) remains a hypothesis suggested by, not established by, the geometry (F2, F3, F5).
- SCI itself is one specific operationalization of "mutation-effect signal carried by a layer pair." Alternative scoring functions could in principle re-order layer pairs; we do not claim uniqueness of the EL enrichment under arbitrary scoring (F1).
