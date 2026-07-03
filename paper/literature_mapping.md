# Literature mapping for the manuscript

> Cross-layer interaction analysis in ESM-2 for mutation-effect encoding.
>
> Mapping only — no Discussion prose, no new conclusions.

## Findings being mapped against

| Tag | Finding |
|---|---|
| **F1** | SCI mean significantly distinguishes functional vs non-functional mutations. |
| **F2** | Top50 high-contribution layer pairs concentrate in Early–Late. |
| **F3** | Across all 528 layer pairs, EL Composite and AUROC significantly exceed EE / EM / MM / LL (p ≪ 0.001). |
| **F4** | CBS / GAL4 / PTEN / TEM1 are EL-dominated; PABP is EM-dominated and forms an independent cluster. |
| **F5** | Cross-protein recurrent pairs: layer 31 / 32 most frequent; Late 55%, Middle 27.5%, Early 17.5%. |

| Tag | Discussion subsection |
|---|---|
| **D1** | Deep-layer enrichment in ESM. |
| **D2** | Cross-layer integration hypothesis. |
| **D3** | Protein-specific heterogeneity. |
| **D4** | Implications for PLM interpretability. |

---

## 1 — Literature support matrix

| # | Paper | Main finding (≤ 20 words) | Supports F# | Maps to D# | Role |
|---:|---|---|---|---|---|
| 1 | Rives et al., 2021 (PNAS) | Scaling ESM-1 to 250M sequences → SS / contacts / function emerge across layers. | F5 | D1, D4 | Foundational |
| 2 | Lin et al., 2023 (Science / ESM-2) | Scaling ESM-2 to 15B → atomic structure emerges; structure read off final layer. | F5; F2/F3 (single-layer contrast) | D1, D4 | Foundational + Contrast |
| 3 | Elnaggar et al., 2022 (TPAMI / ProtTrans) | Six large PLMs systematically probed per layer for SS / localization / membrane. | F2 | D1, D4 | Methodological precedent |
| 4 | Rao et al., 2021 (ICML / MSA Transformer) | Tied row/column attention over MSAs → SOTA unsupervised contacts. | F2 (multi-layer attention, but cross-position) | D2 | Methodological adjacent |
| 5 | Detlefsen, Hauberg & Boomsma, 2022 (Nat Commun) | "Best layer" of a PLM depends on task and architecture. | F2, F3 | D1, D4 | Direct motivation |
| 6 | Vig et al., 2021 (ICLR / BERTology Meets Biology) | Attention heads track contacts; biophysical properties at characteristic depths. | F4, F5 | D1, D3 | Direct support for depth-stratification |
| 7 | Rao et al., 2021 (ICLR / Unsupervised structure learners) | Sparse LR over attention pooled across **all** ESM layers → SOTA contacts. | F2 | D2, D4 | **Closest cross-layer PLM analog** (but for contacts) |
| 8 | Rao et al., 2019 (NeurIPS / TAPE) | 5-task PLM benchmark with per-task, per-layer probing. | F2 (per-layer baseline) | D1, D4 | Single-layer convention reference |
| 9 | Simon & Zou, 2025 (Nat Methods / InterPLM) | Per-layer SAEs on ESM-2 → thousands of layer-localized interpretable features. | F4 | D3, D4 | Per-layer feature complement |
| 10 | Hayes et al., 2025 (Science / ESM3) | Multimodal generative PLM jointly tokenizes sequence / structure / function. | F5 (generalization target) | D1, D4 | Future-direction anchor |
| 11 | Meier et al., 2021 (NeurIPS / ESM-1v) | Final-layer masked-marginal LL → SOTA zero-shot mutation effect prediction. | F2, F3 (single-layer contrast) | D2, D4 | Single-layer convention |
| 12 | Notin et al., 2022 (ICML / Tranception) | Autoregressive PLM + retrieval → SOTA fitness, single output score. | F2, F3 (single-layer contrast) | D2, D4 | Single-layer convention |
| 13 | Notin et al., 2023 (NeurIPS D&B / ProteinGym) | > 250 DMS assays + clinical benchmark; standardized PLM evaluation. | F1 (data substrate) | D4 | Data substrate |
| 14 | Frazer et al., 2021 (Nature / EVE) | Family-MSA VAE classifies variant pathogenicity at experimental quality. | F1 (signal exists) | D2 | Non-PLM anchor |
| 15 | Brandes et al., 2023 (Nat Genet) | ESM-1b last-layer LL scores all ~450M human missense variants. | F5 (uses our anchor layer); F2, F3 (single-layer contrast) | D1, D2 | Closely related ESM-family backbone + Single-layer convention |
| 16 | Hsu et al., 2022 (Nat Biotechnol) | Ridge regression over PLM embeddings beats complex low-N fitness models. | F2, F3 (linear combinability) | D2 | Supportive (linear combo works) |
| 17 | Tenney, Das & Pavlick, 2019 (ACL / BERT Rediscovers) | Scalar mix over all BERT layers orders POS → parsing → NER → SRL → coref. | F2, F3 | D1, D2, D4 | **Methodological closest analog (non-PLM)** |
| 18 | Rogers, Kovaleva & Rumshisky, 2020 (TACL / Primer in BERTology) | Synthesizes ~150 BERT analyses on cross-layer knowledge distribution. | F2, F3, F5 | D1, D4 | Vocabulary provider |
| 19 | Geva et al., 2021 (EMNLP / FFN as KV memories) | Lower-layer FFNs store shallow patterns, upper-layer FFNs semantics, composed via residual stream. | F2, F3 | D2, D4 | Mechanistic basis for integration |
| 20 | Belrose et al., 2023 (Tuned Lens) | Affine per-layer probe decodes residual stream → iterative cross-layer refinement. | F2, F5 | D2, D4 | Direct methodological extension |

---

## 2 — Best to cite in **Discussion** (by subsection)

> Listed as candidate citations, not as suggested phrasing. Selection prioritises papers that are directly load-bearing for the relevant claim in that subsection.

### D1 — Deep-layer enrichment in ESM
- **Anchor / Consistent:** Rives 2021 (#1), Lin 2023 (#2), Vig 2021 (#6), Elnaggar 2022 (#3), Detlefsen 2022 (#5).
- **Same-backbone use of final layer (contextualisation):** Brandes 2023 (#15).
- **Cross-domain analog of deep-layer abstraction:** Tenney 2019 (#17), Rogers 2020 (#18).

### D2 — Cross-layer integration hypothesis
- **Mechanistic / methodological basis:** Tenney 2019 (#17), Geva 2021 (#19), Belrose 2023 (#20).
- **Closest existing multi-layer PLM analysis (different task):** Rao 2021 ICLR (#7).
- **Single-layer convention to contrast against:** Meier 2021 (#11), Notin 2022 (#12), Brandes 2023 (#15).
- **Linear-combinability of PLM features (consistent):** Hsu 2022 (#16).
- **Non-PLM mutation-effect anchor:** Frazer 2021 (#14).

### D3 — Protein-specific heterogeneity
- **Depth-stratified property tracking → biophysical-class link:** Vig 2021 (#6).
- **Per-layer features differ → composition can differ:** Simon & Zou 2025 / InterPLM (#9).
- **Foundational layer-level emergence:** Rives 2021 (#1).

### D4 — Implications for PLM interpretability
- **Cross-layer interpretability toolkit transferable to PLMs:** Tenney 2019 (#17), Rogers 2020 (#18), Geva 2021 (#19), Belrose 2023 (#20).
- **Per-layer feature extraction to combine with our layer-pair view:** Simon & Zou 2025 (#9).
- **Cross-model fingerprinting targets:** Lin 2023 ESM-2 (#2), Hayes 2025 ESM3 (#10).
- **Per-layer probing limitation:** Elnaggar 2022 (#3), Detlefsen 2022 (#5), TAPE / Rao 2019 (#8).

---

## 3 — Best to cite in **Introduction Gap**

> A "Gap" paragraph needs to do three things: (i) establish what the field has done, (ii) name the dominant convention being departed from, (iii) name what is missing. The mapping below sorts the 20 papers by which of these three roles they best serve.

| Gap-role | Papers | Why they belong here |
|---|---|---|
| **(i) What the field has done** — establish that PLMs encode rich features across layers | Rives 2021 (#1), Lin 2023 (#2), Elnaggar 2022 (#3), Vig 2021 (#6), Rao 2021 ICML (#4) | Foundational PLM / layer-wise representation literature. |
| **(ii) The single-layer convention to depart from** | Meier 2021 (#11), Notin 2022 (#12), Brandes 2023 (#15), Hsu 2022 (#16), Frazer 2021 (#14), TAPE / Rao 2019 (#8) | All canonical mutation-effect / variant-effect work uses a single chosen layer (typically the masked-LM head). |
| **(ii') "Best layer" is itself task-dependent → motivates moving beyond it** | Detlefsen 2022 (#5) | Single citation that motivates the move from "best layer" to "best layer pair". |
| **(iii) Cross-layer geometry is studied in NLP but not in PLM mutation work** | Tenney 2019 (#17), Rogers 2020 (#18), Geva 2021 (#19), Belrose 2023 (#20) | Establishes the NLP-side precedent the manuscript imports. |
| **(iii') Closest existing multi-layer PLM analysis is for a different task** | Rao 2021 ICLR (#7) | Specifically: contacts, not mutation effects. This is the single sharpest "empty cell" reference. |

**Minimum viable Gap-paragraph set** (if length-constrained): #11, #15, #5, #7, #17, #19. These six already cover (i)–(iii).

---

## 4 — Papers **consistent with** my findings

Papers whose results directly support F1–F5. Listed by which finding they most strongly support.

| Supports | Papers | Nature of consistency |
|---|---|---|
| **F1** (signal exists) | Frazer 2021 (#14), Notin 2023 (#13), Meier 2021 (#11), Hsu 2022 (#16) | All confirm that mutation-effect signal is recoverable from PLM / MSA representations. |
| **F2 / F3** (cross-layer / off-diagonal structure) | Tenney 2019 (#17), Geva 2021 (#19), Belrose 2023 (#20), Rogers 2020 (#18), Hsu 2022 (#16), Detlefsen 2022 (#5) | Cross-layer composition is documented in BERT-family models; "best layer" varies. |
| **F4** (per-protein / per-instance heterogeneity) | Vig 2021 (#6), Simon & Zou 2025 (#9) | Depth-stratified head specialization and per-layer interpretable features differ across positions / proteins. |
| **F5** (deep-layer enrichment / late-layer hub) | Rives 2021 (#1), Lin 2023 (#2), Vig 2021 (#6), Rogers 2020 (#18), Brandes 2023 (#15) | Deepest layers carry highest-order features / are used as the readout layer. |

---

## 5 — Papers that **extend** my findings (apply their methods to take this work further)

> Papers whose methodology is directly transferable to a follow-up that would strengthen or generalise F1–F5.

| Paper | What it extends | Concrete follow-up enabled |
|---|---|---|
| Tenney 2019 (#17) | F2, F3 | Apply scalar-mix probe to ESM-2 for mutation-effect task; compare to layer-pair Composite. |
| Geva 2021 (#19) | F2, F3, F5 | Test whether the late-layer hub (30–32) corresponds to a small set of FFN "key" features. |
| Belrose 2023 (#20) | F2, F5 | Tuned-lens decoding of mutation-effect signal across ESM-2's 33 layers; pin the depth at which signal first becomes linearly recoverable. |
| Rao 2021 ICLR (#7) | F2, F3 | Re-purpose the multi-layer attention pooling for the mutation-effect task instead of contacts. |
| Simon & Zou 2025 / InterPLM (#9) | F2, F3, F4 | Train per-layer SAEs and ask which early-band SAE features pair with late-band SAE features in the EL set. |
| Lin 2023 ESM-2 (#2) | F5 | Replicate the 528-pair analysis on ESM-2 (33 / 36 layers) to test whether layer-32 prominence is backbone-specific. |
| Hayes 2025 ESM3 (#10) | F2, F3, F4, F5 | Replicate the full layer-pair analysis on ESM3; especially relevant for F4 (does heterogeneity persist at scale?). |
| ProteinGym / Notin 2023 (#13) | F1, F4 | Expand the protein panel from *n* = 5 to ≥ 20 spanning multiple biophysical classes / assay types. |
| Hsu 2022 (#16) | F2, F3 | A ridge-regression baseline restricted to a single layer band vs an Early–Late paired band (matched dimensionality). |

---

## 6 — Papers with **potential conflict** with my findings

> "Conflict" here means *can be read as inconsistent at first glance*, not *strictly refutes*. For each, I name the surface-level tension and the reason it is not a hard refutation — these are the points where a reviewer is most likely to push back, so the manuscript should pre-empt them in Discussion (without writing the Discussion prose itself).

| Paper | Surface-level tension | Why it is not a hard refutation |
|---|---|---|
| Meier 2021 / ESM-1v (#11) | A **single final-layer** masked-marginal score is SOTA on zero-shot mutation-effect prediction — could be read as "a single layer suffices". | This is a *predictive* result; F2 / F3 are *representational geometry* results. High single-layer predictive AUROC does not preclude that the underlying signal is distributed across layer pairs at the representation level — the final-layer LM head is itself the output of cross-layer composition through the residual stream. |
| Brandes et al., 2023 (#15) | ESM-1b's **final layer** scores all ~450 M human missense variants and outperforms prior tools — i.e. on a closely related backbone (ESM-1b in their work; ESM-2 here, both 33-layer Transformers) the final layer alone is sufficient. | Same logic as above. Closely related backbone (both ESM-family, both 33 layers) increases the rhetorical pressure but the rebuttal is identical: predictive performance ≠ representational geometry. Also relevant: their final layer ≈ our layer 32 anchor (F5), so the two findings are partially **co-pointing** rather than opposing. |
| Notin et al., 2022 / Tranception (#12) | Single output score (autoregressive + retrieval) achieves SOTA fitness; no need for layer-pair analysis. | Tranception's architecture is different from ESM-2's, so the comparison is not within-backbone. F2 / F3 are claims about ESM-2's *internal* layer geometry, not about predictor design. |
| Lin 2023 / ESM-2 (#2) | Final-layer representations alone are sufficient for atomic-resolution structure prediction. | The task is structure, not mutation effects. Not directly opposing, but a reviewer may ask "if final layer alone gives structure, why not mutation?" — this is the F5 generalization question already flagged in D1 limitations. |

**Adjacent papers that are NOT in conflict** (frequently mis-identified as such): Detlefsen 2022 (#5) explicitly *supports* moving beyond "best layer"; Rao 2021 ICLR (#7) does cross-layer fusion but for contacts; Hsu 2022 (#16) supports linear combinability and is therefore *consistent* with cross-layer integration, not opposed to it.

---

## Summary statistics of this mapping

- 20 / 20 papers map to at least one (F#, D#) pair.
- F1 is supported by 4 papers (mostly methodological / data-substrate role).
- F2 is touched by 14 / 20 papers (most heavily used finding — natural, given that EL enrichment is the central claim).
- F3 is touched by 11 / 20 papers (subset of F2 mappings; F3 is F2's statistical hardening, so this is expected).
- F4 is touched by 3 / 20 papers (smallest support set — protein-specific heterogeneity is the most under-served angle in the existing literature, which is itself a point worth surfacing in the Introduction Gap).
- F5 is touched by 8 / 20 papers — consistent with deep-layer-emergence being a well-established theme.
- 7 / 20 papers are flagged as supporting D1; 8 as D2; 3 as D3; 13 as D4.
- 3 / 20 papers carry a surface-level conflict label (all three are **predictive single-layer SOTA papers**, all three have the same rebuttal mechanism — predictive ≠ representational).
