# Code Availability

All analysis code used in this study is available through the
project repository and will be released to a public archive
(Zenodo / GitHub) with a permanent DOI at the time of acceptance.

## Repository

- **URL** — [https://github.com/<organisation>/<repository>] (to be
  populated before publication; placeholder retained here for
  editorial routing).
- **DOI** — A versioned snapshot of the repository at the time of
  acceptance will be archived on Zenodo and assigned a DOI. The DOI
  will replace this placeholder in the final published version.
- **Licence** — The repository will be released under an
  OSI-approved open-source licence (MIT or Apache-2.0; final choice
  to be confirmed before submission).

## Pipeline scripts

The repository contains the complete pipeline used to reproduce
every numerical value, figure, and table in the manuscript:

| Stage | Purpose |
|---|---|
| Dataset assembly | Loads the five DMS panels from the ProteinGym substitution release; expands multi-mutation rows into single-substitution entries; merges into the 95,142-row site metadata table (Methods §4.1). |
| Wild-type recovery | Reconstructs the wild-type sequence for each protein from the ProteinGym annotations. |
| ESM-2 feature extraction | Runs the ESM-2 esm2_t33_650M_UR50D forward pass on each unique wild-type and per-mutation variant and caches the 33-layer per-residue Δ-embedding at the mutation position (Methods §4.2). |
| SCI computation | Computes the per-mutation 33 × 33 Pearson correlation matrix between per-layer Δ-embedding vectors; materialises the three per-mutation summary scores (528-pair mean, Top20 mean, Top50 mean) (Methods §4.3). |
| Layer-pair mining | Computes the per-pair Composite discrimination score combining per-pair AUROC, Spearman r, and Cohen's *d* against the binarised functional label; ranks the 528 pairs per protein (Methods §4.4). |
| Protein clustering and recurrence | Hierarchical clustering with Euclidean distance and average linkage on the six-category composition vectors; principal component analysis; cross-protein recurrence counting (Methods §4.5). |
| Permutation testing | Within-protein label-shuffle null with adaptive permutation schedule {1,000 / 5,000 / 10,000}; computes T₁, T₂, T₃, T₃ˢᵗᵈ, T_max under each permutation; Phipson and Smyth (2010) empirical p-values; Bonferroni FWER and Benjamini–Hochberg FDR over the 20 per-protein primary tests (Methods §4.6–§4.8). |
| Bootstrap | Non-parametric bootstrap with N = 1,000 resamples; rank-reuse fast path for T₂/T₃/T₃ˢᵗᵈ/T_max (Methods §4.9). |
| ESM-2 LLR baseline | WT-marginal and masked-marginal LLR computation on the same esm2_t33_650M_UR50D checkpoint; one forward pass per protein for WT-marginal, one forward pass per unique (protein, position) for masked-marginal (≈ 1,300 forwards total) (Methods §4.12). |
| SCI vs LLR comparison | 5-fold stratified-CV logistic regression on the standardised (SCI, LLR) feature pair (StratifiedKFold; StandardScaler fit within each training fold; out-of-fold P(class = 1) as the predictor); per-protein and pooled Spearman / AUROC / AUPRC; sign-corrected. Pooled SCI–LLR Pearson and Spearman correlations on the same 95,142-row pair after NaN removal (Methods §4.12). |
| Figure generation | Produces every panel of Figures 1–5 and the rows of Table 1 from the on-disk CSV outputs of the above stages. |

## Reproducibility

- All random number generation uses fixed seeds — `42` for the
  permutation stream and `43` for the bootstrap stream — ensuring
  bit-for-bit reproducibility of the reported z-scores, empirical
  p-values, confidence intervals, and figures.
- The 5-fold stratified-CV logistic regression for the SCI+LLR
  combined scorer uses `StratifiedKFold(n_splits = 5, shuffle = True,
  random_state = 42)`.
- The full permutation + bootstrap pipeline runs in approximately
  8–10 minutes single-threaded on a CPU at the real scale
  (95,142 mutations, 528 layer pairs, 1,000 adaptive permutations
  + 1,000 bootstrap resamples), with peak resident memory ≈ 410 MB.
- The ESM-2 LLR baseline computation requires approximately
  1,305 forward passes through esm2_t33_650M_UR50D. This is ≈ 70×
  fewer than the 95,142 forward passes used for the SCI Δ-embedding
  cache. A single consumer-grade GPU is sufficient.

## Dependencies

The repository documents its complete software environment, including
the exact versions of:

- Python and the scientific Python stack (NumPy, SciPy, pandas,
  scikit-learn).
- `fair-esm` (for the ESM-2 forward pass and the LLR baseline).
- PyTorch (for the underlying model inference).
- Plotting libraries (matplotlib, seaborn).

A pinned `requirements.txt` (or equivalent `environment.yml`) is
included in the repository.

## Issue tracking

Issues, questions, and reproduction reports may be filed on the
repository's public issue tracker. The corresponding author will
respond to substantive reproducibility questions during the period
of active maintenance following publication.

## Correspondence

Requests for code components not yet on the public archive, or for
assistance in reproducing specific analyses, should be addressed to
the corresponding author at [email].
