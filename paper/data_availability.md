# Data Availability

All data analysed in this study are derived from publicly released
deep mutational scanning (DMS) panels and from publicly released
protein-language-model weights; no new experimental data were
generated.

## Source data

The five DMS panels used in this study were drawn from the
**ProteinGym substitution release** (Notin et al., 2023). The
specific source datasets, with author attribution, are:

| Protein | Source dataset | Reference |
|---|---|---|
| CBS  | CBS_HUMAN_Sun_2020 | Sun et al., 2020 |
| GAL4 | GAL4_YEAST_Kitzman_2015 | Kitzman et al., 2015 |
| PABP | PABP_YEAST_Melamed_2013 | Melamed et al., 2013 |
| PTEN | PTEN_HUMAN_Mighell_2018; PTEN_HUMAN_Matreyek_2021 | Mighell et al., 2018; Matreyek et al., 2021 |
| TEM1 (BLAT_ECOLX) | BLAT_ECOLX_Deng_2012; BLAT_ECOLX_Firnberg_2014; BLAT_ECOLX_Jacquier_2013; BLAT_ECOLX_Stiffler_2015 | Deng et al., 2012; Firnberg et al., 2014; Jacquier et al., 2013; Stiffler et al., 2015 |

Per-protein counts, label balance, and merged dataset characteristics
are summarised in Supplementary Table S1 of the manuscript. The full
ProteinGym substitution release is available under its published
licence at the URL provided in Notin et al., 2023.

## Processed data

The following processed datasets generated for this study are
publicly available at
[https://github.com/DCMPjmpo/papae_sci](https://github.com/DCMPjmpo/papae_sci):

- The merged single-substitution table covering all five proteins
  (95,142 rows), including the binarised functional label
  `DMS_score_bin`, the continuous `DMS_score`, the parsed
  `(wt_aa, mut_aa, mutation_position)` triple, and the per-row
  dataset identifier.
- Wild-type sequences for each of the five proteins.
- The per-mutation 33-layer Δ-embedding cache extracted from
  ESM-2 esm2_t33_650M_UR50D at the mutation position
  (95,142 × 33 × 1280 array, ≈ 9.4 GB at float16).
- The per-mutation 33 × 33 SCI correlation matrices
  (95,142 × 33 × 33 array, ≈ 414 MB at float32) and the three
  per-mutation summary vectors (528-pair mean; Top20 mean; Top50 mean).
- The per-protein and pooled SCI-vs-LLR comparison tables
  (Spearman r, AUROC, AUPRC) used to populate Table 1 of the
  manuscript.
- The SCI–LLR Pearson and Spearman correlation tables (pooled and
  per-protein) used in Results §2.8 and Discussion D5.
- The ESM-2 zero-shot LLR scores (WT-marginal and masked-marginal
  variants) for all 95,142 rows.
- The permutation null distributions, per-protein and global
  summary tables, and bootstrap confidence intervals reported in
  Supplementary Table S2.

## Model weights

All analyses use the publicly released **ESM-2 esm2_t33_650M_UR50D**
checkpoint (33 Transformer blocks, hidden dimension 1280; Lin et al.,
2023). The model weights are available from the official ESM
distribution (Meta AI / Facebook Research). No fine-tuning was
performed; the released checkpoint is used as-is.

## Restrictions

No restrictions apply to the use of the data deposited by this study.
The ProteinGym source datasets are subject to the licence of the
ProteinGym release (Notin et al., 2023) and the licences of their
respective primary references. The ESM-2 model weights are
distributed under their original licence.

## Correspondence

Requests for processed data tables, intermediate matrices, or
analysis outputs not yet on the public archive should be addressed
to the corresponding author at [email].
