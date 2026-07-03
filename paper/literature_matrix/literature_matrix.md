# Literature matrix — Cross-layer interaction analysis in PLMs for mutation-effect encoding

**Scope.** 20 papers (2019–2025, heavy 2021–2025) across four discussion categories:

- **D1** — ESM / PLM layer-wise representation studies
- **D2** — Protein language-model interpretability
- **D3** — Mutation-effect prediction using PLMs
- **D4** — General transformer representation analysis (BERTology, cross-layer geometry)

**Y/N criteria.**
- *Uses single-layer representation?* — does the paper's main result rely on a single chosen layer (commonly the final layer) of the model?
- *Studies cross-layer interactions?* — does the paper EXPLICITLY analyze how multiple layers JOINTLY contribute, fuse, or interact (scalar mix, residual-stream composition, multi-layer attention pooling, layer-pair geometry, etc.)? Layer-by-layer probing without joint analysis is **not** counted as cross-layer.

> Out of 20 papers, only **6** explicitly study cross-layer interactions (5 in D4, 1 in D2). No paper to date specifically analyzes **layer-pair geometry of mutation-effect encoding**, which is the gap this manuscript addresses.

---

## D1 — ESM / PLM layer-wise representation studies

| # | Citation | Venue | Yr | Main finding | Single-layer? | Cross-layer? | Relationship to my work |
|---:|---|---|---:|---|:---:|:---:|---|
| 1 | [Rives et al., 2021](https://www.pnas.org/doi/10.1073/pnas.2016239118) | PNAS | 2021 | Scaling masked-LM pretraining (ESM-1) to 250M sequences causes secondary structure, contacts, and function to emerge in learned representations; surveyed layer-by-layer. | Y | N | Establishes that ESM layers carry emergent biological signal — my work asks how these layer-localized signals are **jointly recombined** into mutation-effect representations. |
| 2 | [Lin et al. (ESM-2 / ESMFold), 2023](https://www.science.org/doi/10.1126/science.ade2574) | Science | 2023 | Scaling ESM-2 up to 15B parameters causes atomic-level structural information to emerge, enabling single-sequence structure prediction. | Y | N | Structure recoverable from a final-layer representation alone — my work shows mutation-effect signal is, by contrast, **distributed across pairs** of layers. |
| 3 | [Elnaggar et al. (ProtTrans), 2022](https://arxiv.org/abs/2007.06225) | IEEE TPAMI | 2022 | Trains six large PLMs and probes per-layer embeddings for secondary structure, localization, membrane prediction. | Y | N | Provides a per-layer (not per-pair) probing protocol — my analysis extends this from a one-axis (layer index) to a **two-axis (layer pair) view**. |
| 4 | [Rao et al. (MSA Transformer), 2021](https://proceedings.mlr.press/v139/rao21a.html) | ICML | 2021 | Tied row/column attention over MSAs achieves SOTA unsupervised contact prediction. | Y | N | Cross-**position** (not cross-layer) attention is the inductive bias — my work studies cross-**layer**, not cross-residue, interactions for mutation effects. |
| 5 | [Detlefsen, Hauberg, Boomsma, 2022](https://www.nature.com/articles/s41467-022-29443-w) | Nat Commun | 2022 | Representation geometry across PLM layers strongly affects interpretability; the 'best layer' depends on task and architecture. | N | N | Argues against any single 'best layer' — my work generalizes from 'best layer' to **'best layer pair'**, revealing structure invisible to per-layer probing. |

## D2 — Protein language model interpretability

| # | Citation | Venue | Yr | Main finding | Single-layer? | Cross-layer? | Relationship to my work |
|---:|---|---|---:|---|:---:|:---:|---|
| 6 | [Vig et al. (BERTology Meets Biology), 2021](https://arxiv.org/abs/2006.15222) | ICLR | 2021 | Protein-Transformer attention heads track 3D contacts and binding sites; higher-level biophysical properties are tracked by deeper-layer heads. | N | N | Depth-stratified interpretation of attention — my work extends this depth view to a **layer-pair view** focused on mutation-effect encoding. |
| 7 | [Rao et al. (PLMs are unsupervised structure learners), 2021](https://iclr.cc/virtual/2021/poster/3016) | ICLR | 2021 | Sparse logistic regression over attention heads pooled across **all** ESM layers extracts SOTA unsupervised residue–residue contacts. | N | **Y** | Methodologically the closest existing multi-layer fusion in PLMs — but for **contacts**, not mutation effects; my work shows mutation signal demands a structurally similar yet biologically distinct multi-layer combination. |
| 8 | [Rao et al. (TAPE), 2019](https://proceedings.neurips.cc/paper/2019/hash/37f65c068b7723cd7809ee2d31d7861c-Abstract.html) | NeurIPS | 2019 | Releases a 5-task PLM benchmark with per-task probing of recurrent/conv/Transformer models. | Y | N | Codifies the single-layer probing protocol my work challenges — per-task benchmarks may understate signal carried in **layer pairs**. |
| 9 | [Simon & Zou (InterPLM), 2025](https://www.nature.com/articles/s41592-025-02836-7) | Nat Methods | 2025 | Per-layer sparse autoencoders on ESM-2 expose thousands of interpretable, layer-localized biological features and concept superposition. | N | N | SAE features are extracted **within** each layer separately — my work asks the complementary question of **how features in different layers interact**. |
| 10 | [Hayes et al. (ESM3), 2025](https://www.science.org/doi/10.1126/science.ads0018) | Science | 2025 | Multimodal generative PLM jointly tokenizes sequence, structure, and function for prompted protein design. | Y | N | Scale-up + generative angle; my work focuses on **interpretation** of an earlier, fully analyzable variant (ESM-2). Orthogonal aims. |

## D3 — Mutation-effect prediction using PLMs

| # | Citation | Venue | Yr | Main finding | Single-layer? | Cross-layer? | Relationship to my work |
|---:|---|---|---:|---|:---:|:---:|---|
| 11 | [Meier et al. (ESM-1v), 2021](https://proceedings.neurips.cc/paper/2021/hash/f51338d736f95dd42427296047067694-Abstract.html) | NeurIPS | 2021 | Final-layer masked-marginal log-likelihood gives SOTA zero-shot mutation-effect prediction. | Y | N | The canonical 'final layer = mutation score' approach — my work re-examines this assumption and shows the signal is **distributed across layer pairs**. |
| 12 | [Notin et al. (Tranception), 2022](https://proceedings.mlr.press/v162/notin22a.html) | ICML | 2022 | Autoregressive PLM with multi-scale conv-attention and retrieval-augmented inference achieves SOTA fitness prediction. | Y | N | My work is **not** a competing predictor — it is a representational analysis of where in the network mutation information lives. |
| 13 | [Notin et al. (ProteinGym), 2023](https://openreview.net/forum?id=URoZHqAohf) | NeurIPS D&B | 2023 | >250 DMS assays + clinical benchmark with standardized PLM/MSA evaluation. | Y | N | Provides the DMS datasets used here — my work uses them as **encoding-analysis substrates**, not a leaderboard. |
| 14 | [Frazer et al. (EVE), 2021](https://www.nature.com/articles/s41586-021-04043-8) | Nature | 2021 | Family-MSA VAE classifies variant pathogenicity on par with high-throughput experiments. | Y | N | Non-PLM mutation-effect anchor — my work targets the layer-organization question **specific to deep Transformers**. |
| 15 | [Brandes et al., 2023](https://www.nature.com/articles/s41588-023-01465-0) | Nat Genet | 2023 | Uses ESM-1b last-layer log-likelihood to score all ~450M human missense variants. | Y | N | Closely related ESM-family backbone (ESM-1b there; ESM-2 here) and same DMS evaluations — my work asks not how well its **final layer** scores variants, but **how mutation information is organized across all 33 layers**. |
| 16 | [Hsu et al., 2022](https://www.nature.com/articles/s41587-021-01146-5) | Nat Biotechnol | 2022 | Simple ridge regression over evolutionary density-model scores + assay labels outperforms more complex low-N fitness predictors. | Y | N | Supports the view that mutation-effect signal is geometrically simple and recoverable from **linear combinations of representations** — consistent with our cross-layer integration finding. |

## D4 — General transformer representation analysis

| # | Citation | Venue | Yr | Main finding | Single-layer? | Cross-layer? | Relationship to my work |
|---:|---|---|---:|---|:---:|:---:|---|
| 17 | [Tenney, Das & Pavlick (BERT Rediscovers), 2019](https://aclanthology.org/P19-1452/) | ACL | 2019 | Scalar-mix probe over all BERT layers shows POS → parsing → NER → SRL → coreference localized in order, with cross-layer feedback. | N | **Y** | Methodologically the **closest non-PLM analog** to this work: explicit multi-layer attribution. Inspires the 'where along depth, across which pairs?' framing. |
| 18 | [Rogers, Kovaleva & Rumshisky (Primer in BERTology), 2020](https://aclanthology.org/2020.tacl-1.54/) | TACL | 2020 | Synthesizes ~150 BERT analyses on how knowledge is distributed across layers and heads. | N | **Y** | Provides the conceptual vocabulary (low/mid/high layer roles) imported into PLMs — my project produces the **PLM analog of this layer cartography** for mutation-effect signal. |
| 19 | [Geva et al. (FFN as key-value memories), 2021](https://aclanthology.org/2021.emnlp-main.446/) | EMNLP | 2021 | Transformer FFNs act as key-value memories; lower layers store shallow patterns, upper layers semantic ones; residual-stream composes them across layers. | N | **Y** | Mechanistic basis for treating the residual stream as a cross-layer integrator — directly motivates why **Early–Late layer pairs** could carry mutation-effect signal. |
| 20 | [Belrose et al. (Tuned Lens), 2023](https://arxiv.org/abs/2303.08112) | arXiv | 2023 | Learned affine per-layer probe decodes residual-stream states into vocabulary, exposing iterative refinement across layers. | N | **Y** | Layer-by-layer decoding trajectory — methodological analog to mapping mutation-effect signal across PLM depth; reinforces interpreting mutation encoding as a **cross-layer process**. |

---

## Summary statistics for the matrix

| Property | Count | Proportion |
|---|---:|---:|
| Uses single-layer representation (Y) | 12 / 20 | 60% |
| Studies cross-layer interactions (Y) | 6 / 20 | 30% |
| Cross-layer **AND** PLM (D1+D2+D3, Y) | 1 / 16 PLM papers | 6% |
| Cross-layer **AND** mutation-effect (D3, Y) | 0 / 6 D3 papers | 0% |

## Positioning statement (for Introduction / Discussion)

The literature partitions cleanly into two non-overlapping camps:

- **PLM mutation-effect work (D3)** universally relies on a **single layer** (typically the final masked-LM head) for the variant score, even when the score is downstream of a sophisticated model (Meier 2021; Notin 2022; Brandes 2023; Hsu 2022; Frazer 2021).
- **General transformer interpretability (D4)** does the opposite — it analyzes how layers **jointly compose** representations via residual stream and scalar mixing (Tenney 2019; Geva 2021; Belrose 2023; Rogers 2020), but does so on language, not protein, models, and not for mutation effects.

The single existing protein-LM work that explicitly fuses across layers (Rao et al., 2021, ICLR) does so for **residue contacts**, not mutation effects. The 528-layer-pair analysis presented here therefore occupies an empty cell in this 2D landscape: **cross-layer interaction analysis × mutation-effect encoding × protein language models**.
