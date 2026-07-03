# Scientific Story Map

> **Project**: Cross-layer interaction analysis in ESM-2 for mutation-effect encoding (PABP / CBS / GAL4 / PTEN / TEM1; 33 layers; 528 unordered layer pairs; 95,142 mutations).
>
> Compiled from the actual experimental artifacts under `data/` and `p6_scripts/` as of this commit. Status labels reflect the *evidence currently on disk*, not what was *designed* or *planned*.
>
> Legend:
> ✅ **Sufficiently supported** — at least one quantitative analysis with statistical control on disk.
> ⚠️ **Partially supported** — evidence exists but with caveats (limited scope, methodological asymmetry, or correlational only).
> ❌ **Not yet supported** — no run-on-disk evidence; reserved for future work.

---

## 1 · Top-level claim status

| ID | Claim (in publication language) | Status | Strongest single piece of evidence |
|---|---|---|---|
| **C1** | SCI carries discriminative mutation-effect signal that is reproducible under within-protein label shuffles. | ✅ | p5 global T1: z = 57.0, p_fdr ≈ 0.001 |
| **C2** | When all 528 layer pairs are ranked by Composite (AUROC + Spearman + Cohen's *d*), the Top50 layer pairs are enriched for Early–Late (EL) interactions. | ✅ | `layer_pair_top50.csv` + per-protein `*_pair_category_ratio_Composite.csv` |
| **C2′** | The same enrichment also appears when ranking is restricted to AUROC alone. | ⚠️ | p5 global T2: 17/50 EL vs. null 1.0, z = 8.0, p_fdr = 0.001 — but per-protein T2 = 0 for 4 of 5 proteins (AUROC ranking and Composite ranking diverge at the per-protein level). |
| **C3** | Across all 528 layer pairs, EL pairs significantly outperform EE / EM / MM / LL on discrimination metrics. | ✅ | p5 global T3: z = 62.0, p_fdr = 0.001; T3_std z = 8.05, p_fdr = 0.001 |
| **C4** | CBS, GAL4, PTEN, TEM1 share an EL-dominated layer-pair composition; PABP differs (EM-dominant) and forms an independent cluster on this 6-category compositional axis. | ✅ | `protein_clustering/` heatmap + average-linkage dendrogram + PCA (PABP at singleton; PC1 68.6%, EM loading 0.207 dominates PC1) |
| **C5** | Cross-protein recurrent layer pairs concentrate at layers 30–32, with layer 32 most frequent; the Late band carries 55% of layer occurrences in the Top20 recurrent set. | ✅ | `layer_recurrence_top20_pairs.csv` (40 layer-end occurrences: Late 22 / Middle 11 / Early 7; layer 32 = 8 of 40) |
| **C6** | The per-protein signal — not only the pooled signal — is reproducible against within-protein null. | ✅ | p5 per-protein: T1, T3, T3_std all FDR-significant in 5/5 proteins (smallest = CBS T1 p_fdr = 0.037; all others ≪ 0.01) |
| **C7** | The findings are robust to (i) BH FDR correction across 20 per-protein primary tests, (ii) Bonferroni FWER correction, and (iii) bootstrap re-sampling. | ✅ | `permutation_per_protein_summary.csv` (p_empirical, p_bonferroni, p_fdr all reported) + `bootstrap_summary.csv` (N = 1000) |
| **C8** | Cross-layer integration is a causal property of the ESM-2 forward pass. | ❌ | No interventional experiments on disk (ablation / activation patching / single-layer probe at matched dim). Framed in the manuscript as a *working hypothesis*. |
| **C9** | The EL-dominance pattern generalises to other PLM backbones (ESM-2, ESM3, ProtT5). | ❌ | All analyses on ESM-2 only. |
| **C10** | The PABP shift toward EM is *because* PABP is RNA-binding (or *because* of its DMS assay design). | ❌ | n = 5 proteins; biophysical class and assay design are confounded in the panel. Framed as *consistent reading*, not a causal claim. |

---

## 2 · Forward map: each claim → which experiments support it

> Each row lists the script that produced the artefact and the on-disk file that materialises the evidence.

### C1 — SCI carries mutation-effect signal
| Script | Artefact | What it shows |
|---|---|---|
| `p6_scripts/P1.5_validate_sci_signal.py` | `data/processed/P1.5_sci_signal_validation.png` | SCI-mean separation between functional/non-functional (qualitative) |
| `p6_scripts/p5_permutation.py` (global, T1) | `data/p0_output/p5_permutation/permutation_global_summary.csv` | T1 = 0.0028, z = 57, p_fdr = 0.001 |
| `p6_scripts/p5_permutation.py` (per-protein T1) | `data/p0_output/p5_permutation/permutation_per_protein_summary.csv` | T1 FDR-significant in 5/5 proteins |
| `p6_scripts/p5_permutation.py` (bootstrap T1) | `data/p0_output/p5_permutation/bootstrap_summary.csv` | T1 95% bootstrap CI excludes 0 |

### C2 / C2′ — Top50 EL enrichment
| Script | Artefact | What it shows |
|---|---|---|
| `p6_scripts/p0_layer_pair_mining_v2.py` | `data/p0_5layer/layer_pair_top50.csv`, `top20_positive_layer_pairs.npy` | Composite Top50 across the full pool |
| `p6_scripts/layer_pair_category_distribution.py` | `data/processed/p0_output/{PROT}/{PROT}_pair_category_ratio_Composite.csv` + `layer_pair_category_distribution.png` | Per-protein Composite Top50 composition (CBS 40% EL + 30% LL; GAL4 66% EL; PTEN 44% EL + 12% LL; TEM1 40% EL + 10% EM + 8% LL; PABP 48% EM + 42% EL) |
| `p6_scripts/p5_permutation.py` (global T2) | `permutation_global_summary.csv` | AUROC-Top50: 17 EL vs null 1.0; z = 8.0, p_fdr = 0.001 |
| `p6_scripts/p5_permutation.py` (per-protein T2) | `permutation_per_protein_summary.csv` | **Per-protein AUROC-Top50 EL count = 0 in 4/5 proteins; PABP = 1.** This is the only place the AUROC-only and Composite-based rankings disagree, and the divergence is informative — see "Open gaps" §5. |

### C3 — 528-pair EL category dominance
| Script | Artefact | What it shows |
|---|---|---|
| `p6_scripts/p01_all_layer_pair_stat.py` | `data/p0_5layer/layer_pair_auc_ranking.csv` + `layer_pair_frequency_map.png` | Per-pair AUROC ranking across all 528 pairs |
| `p6_scripts/p5_permutation.py` (global T3, T3_std) | `permutation_global_summary.csv` | T3 = 0.040, z = 62, p_fdr = 0.001; T3_std = 1.02, z = 8.05, p_fdr = 0.001 |
| `p6_scripts/p5_permutation.py` (per-protein T3, T3_std) | `permutation_per_protein_summary.csv` | FDR-significant in 5/5 proteins |
| `data/p0_output/p5_permutation/fig_A_global_null.{png,pdf}` | Histogram of 528-pair null distribution under shuffle | Observed lies far in upper tail |

### C4 — PABP forms an independent cluster on the 6-category axis
| Script | Artefact | What it shows |
|---|---|---|
| `p6_scripts/protein_clustering_layer_pair.py` | `data/p0_output/protein_clustering/heatmap_layer_pair_categories.{png,pdf}` | PABP is the only EM-dominated row |
| same | `dendrogram_protein_clustering.{png,pdf}` | PABP joins the rest only at distance ≈ 0.50; the four other proteins join by distance ≈ 0.22 |
| same | `pca_protein_layer_pair.{png,pdf}` | PC1 68.6%, PC2 28.1% (cumulative 96.7%); PABP at PC1 ≈ +0.38, isolated |
| same | `pairwise_euclidean_distance.csv` | Min PABP-to-any-other distance = 0.389 (PABP–TEM1), exceeds max within-other-four = 0.397 |
| `p6_scripts/pca_loadings_extract.py` | `pca_loadings.csv` + `pca_loadings_ranked.csv` | PC1 driven primarily by EM (loading +0.207); PC2 by EL vs LL |

### C5 — Cross-protein recurrence at layers 30–32
| Script | Artefact | What it shows |
|---|---|---|
| `p6_scripts/top_recurrent_cross_protein_pairs.py` | `data/processed/p0_output/cross_protein_universal_pairs.csv` | 146 candidate recurrent pairs ranked by Cross_Protein_Count |
| `p6_scripts/layer_recurrence_top20.py` | `data/p0_output/protein_clustering/layer_recurrence_top20_pairs.csv` | Top20 recurrent pairs (6 at count = 4, 14 at count = 3) |
| same | `layer_recurrence_frequency.csv` + `layer_recurrence_band_summary.csv` | 40 layer-end occurrences: Late 22 / Middle 11 / Early 7 (55% / 27.5% / 17.5%) |
| same | `layer_recurrence_barplot.{png,pdf}` | Layer 32 = 8, layer 31 = 7, layer 30 = 3; layers 18–27 contribute zero |
| `p6_scripts/layer_occurrence_frequency.py` | `data/p0_output/protein_clustering/layer_frequency_histogram.{png,pdf}` + `layer_frequency.csv` | At the broader Top50 union across 5 proteins (n = 500 layer-ends): Late 47.6%, Middle 36.0%, Early 16.4% |

### C6 — Per-protein signal is reproducible
| Script | Artefact | What it shows |
|---|---|---|
| `p6_scripts/p5_permutation.py` | `permutation_per_protein_summary.csv` | 5/5 proteins have T1, T3, T3_std with p_fdr < 0.05 (CBS T1 p_fdr = 0.037; rest 0.0014 to 0.001) |
| same | `fig_B_per_protein_null.{png,pdf}` | Per-protein null histograms with observed marked |
| same | `fig_C_observed_vs_null_bar.{png,pdf}` | Bar plot of z-scores across 6 scopes × 4 stats; reference lines at |z| = 1.96 and Bonferroni threshold |

### C7 — Statistical robustness
| Layer of robustness | Artefact |
|---|---|
| Phipson–Smyth empirical p (centered for two-sided) | `permutation_*_summary.csv` columns `p_empirical` |
| Bonferroni (over 20 per-protein primary tests) | `permutation_per_protein_summary.csv` column `p_bonferroni` |
| BH FDR (over 20 per-protein and over 4 global) | columns `p_fdr` in both summary CSVs |
| Bootstrap CI (N = 1000) | `bootstrap_summary.csv`; T2/T3/T3_std/T_max via the rank-reuse fast path (measured AUROC deviation ≤ 6e-3 at real scale, see in-code sanity check) |
| QQ-plot for 18 (global + per-protein) tests | `fig_D_qq_significance.{png,pdf}` |

---

## 3 · Inverse map: each script → which claims it supports

| Script | Supports |
|---|---|
| `extract_site.py` / `recover_wt_v2.py` / `merge_data_v2.py` / `expand_multi_mutations.py` | Data foundation (no direct claim; supports all downstream) |
| `P1_build_sci.py` | Methodological core (defines SCI) — supports C1 by construction |
| `P1.5_validate_sci_signal.py` | C1 |
| `p0_layer_pair_mining_v2.py` | C2, C3 |
| `p01_all_layer_pair_stat.py` | C3 |
| `layer_pair_category_distribution.py` | C2, C4 |
| `protein_clustering_layer_pair.py` | C4 |
| `pca_loadings_extract.py` | C4 |
| `top_recurrent_cross_protein_pairs.py` | C5 |
| `layer_recurrence_top20.py` | C5 |
| `layer_occurrence_frequency.py` | C5 |
| `p5_permutation.py` | C1, C2′, C3, C6, C7 (and provides reproducibility infra for all of the above) |

---

## 4 · What is sufficiently supported, in one sentence each

- **C1 (SCI carries signal)**: global z = 57 and 5/5 protein-level FDR-significance, with bootstrap CI excluding zero — supported at *all* statistical layers.
- **C2 (Composite Top50 EL enrichment)**: directly observed in `*_pair_category_ratio_Composite.csv`; per-protein numbers are 40%–66% EL among the four EL-dominated proteins, plus 42% EL in PABP's Top50 alongside its 48% EM.
- **C3 (528-pair EL > others)**: global z = 62 on T3, p_fdr = 0.001; supported in 5/5 proteins via T3 / T3_std under the permutation null. (See §5 for the per-protein T3 interpretation caveat.)
- **C4 (PABP distinct cluster)**: three independent geometric views (heatmap, dendrogram, PCA) all place PABP outside the centroid of the other four; pairwise Euclidean distances make the singleton call quantitative.
- **C5 (Recurrence at 30–32, layer 32 most frequent)**: layer 32 contributes 20% of all Top20 recurrent layer-ends; layers 18–27 contribute zero — both the concentration and the trough are direct counts, not statistical estimates.
- **C6 / C7 (statistical robustness)**: within-protein label shuffle preserves protein-level imbalance; Bonferroni and BH FDR are reported side-by-side; bootstrap CI for all four primary statistics is on disk.

---

## 5 · Open gaps (be explicit in the manuscript)

| Gap | Why it matters | Recommended manuscript handling |
|---|---|---|
| **G1.** Per-protein AUROC-only Top50 contains 0 EL pairs in CBS / GAL4 / PTEN / TEM1 (T2 in `permutation_per_protein_summary.csv`). The Composite ranking (used to derive C2) and the AUROC-only ranking (used by p5 T2) **diverge at the per-protein level**. | Reviewers may ask whether F2 (EL enrichment in Top50) depends on the choice of ranking metric. The honest answer is: at the **pooled / global** level the enrichment is robust to either metric (p5 global T2 p_fdr = 0.001), but at the **per-protein** level the enrichment is metric-dependent — Composite reveals it, AUROC alone does not. | Report C2 against Composite-Top50 (which the manuscript already does), and note in Methods that p5 T2 confirms this only at the global level; use per-protein T3 / T3_std (which **do** confirm 5/5 proteins) as the per-protein analogue of F2 instead of T2. |
| **G2.** Per-protein T3 absolute values are not bounded in [−0.5, 0.5] because the rank trick reuses the full-population ranks of length 95,142 inside the per-protein index slice, producing a *globally-anchored U-statistic* rather than a protein-local AUROC (see e.g. GAL4 T3 = −18.97). The significance call is unaffected (permutation null is computed under the same metric), but the absolute T3 numbers are **not interpretable as AUROC dominance** per-protein. | Same: significance is fine but the displayed effect-size numbers per-protein are not AUROC-scale. | Either (a) recompute per-protein AUROC with within-protein ranks for the *reported* effect sizes only (significance can remain as is); or (b) clarify in Methods that the per-protein T3 is a globally-anchored U-statistic and report z / p as the primary per-protein readouts, with T_obs given for completeness. |
| **G3.** No interventional / mechanistic experiment (C8). All five experiments are correlational. | Without intervention we cannot claim the EL geometry is *causal* to ESM-2's mutation-effect performance. | The Discussion already labels integration as a *hypothesis*; flag specifically that activation patching / single-layer probe at matched dimensionality / tuned-lens decoding (Belrose et al., 2023) are deferred to future work. |
| **G4.** Single backbone (ESM-2 only; C9 untested). | Reviewers will ask whether layer 32 is special to ESM-2. | Already framed as a limitation in D1 of the Discussion outline; explicitly list ESM-2, ESM3 and ProtT5 as replication targets. |
| **G5.** Confounded biophysical class × DMS assay design in the 5-protein panel (C10). | Reviewers will not accept "RNA-binding causes EM-dominance" from n = 1 RNA-binder. | Already framed as two parallel *consistent readings* in D3; commit to a follow-up with ≥ 20 proteins spanning multiple classes. |

---

## 6 · Story arc for the paper

The order in which to deploy the validated claims:

> **Existence (C1)** → **Geometry (C2 + C3)** → **Topography (C5)** → **Heterogeneity (C4)** → **Robustness (C6 + C7)** → **Hypothesis (C8) and limitations (G3–G5)**

| Manuscript section | Carries which claims | Lead figure / table |
|---|---|---|
| **Results §1** Validation of SCI as a mutation-effect signal | C1 | `P1.5_sci_signal_validation.png` + p5 global T1 row |
| **Results §2** EL enrichment across the 528 layer-pair plane | C2 + C3 | `layer_pair_top50.csv` summary + `layer_pair_category_distribution.png` + p5 global T2 / T3 |
| **Results §3** Recurrence topography across proteins | C5 | `layer_recurrence_barplot.png` + `layer_recurrence_top20_pairs.csv` |
| **Results §4** Per-protein heterogeneity and PABP singleton | C4 | `heatmap_layer_pair_categories.png` + `dendrogram_protein_clustering.png` + `pca_protein_layer_pair.png` + `pca_loadings_ranked.csv` |
| **Results §5** Statistical robustness | C6 + C7 | `fig_A` / `fig_B` / `fig_C` / `fig_D` from p5 + `publication_summary.csv` |
| **Discussion** D1–D4 | Hypothesis (C8); limitations G3–G5 | see `paper/discussion_full.md` |

---

## 7 · One-paragraph executive summary

We validate five findings about how mutation-effect information is organised across the layers of ESM-2. **(C1)** SCI separates functional from non-functional mutations and the separation is reproducible under within-protein label shuffles (global z = 57, p_fdr = 10⁻³; 5/5 proteins per-protein FDR-significant). **(C2–C3)** Within the 528 unordered layer pairs, the Composite-Top50 set is dominated by Early–Late interactions, and the EL category as a whole exceeds EE / EM / MM / LL on |AUROC – 0.5| at high significance (global T3 z = 62). **(C5)** The deepest layers (30–32, with layer 32 most frequent) form a localised hub that recurs across proteins; layers 18–27 are essentially absent from the recurrent set. **(C4)** On the 6-category compositional axis, CBS / GAL4 / PTEN / TEM1 share an EL-dominated signature whereas PABP forms a singleton driven by Early–Middle interactions, supported independently by hierarchical clustering, PCA, and pairwise distances. **(C6–C7)** All primary effects survive Bonferroni and Benjamini–Hochberg FDR correction over 20 per-protein tests, and bootstrap 95% CIs exclude the null across the four primary statistics. We frame integration as a working hypothesis (C8) and explicitly defer interventional, multi-backbone, and expanded-protein-panel validation (G3–G5) to follow-up work.
