# P7 — SCI vs LLR correlation summary

> Source numbers: `data/processed/sci_llr_correlation.csv` (Pearson +
> Spearman, Pooled + per-protein) and
> `data/processed/sci_vs_llr_comparison.csv` (AUROC for SCI, LLR,
> SCI+LLR). All values are reported as computed; no rounding beyond
> three decimals is applied.

---

## Results

We next asked how the per-mutation SCI score relates to the canonical
ESM-2 masked-marginal Log-Likelihood Ratio (LLR), using the same
n = 95,142 single-residue rows aligned through
`all_proteins_site_metadata.csv`. At the pooled level, SCI and LLR
were positively but only moderately correlated (Pearson r = 0.456;
Spearman ρ = 0.467), indicating that the two scores share a
substantial fraction of unexplained variance. The per-protein
correlations spanned a wide range: PTEN was the weakest
(Pearson 0.185; Spearman 0.214) and CBS was the second-weakest
(Pearson 0.247; Spearman 0.237), while GAL4 (Pearson 0.561;
Spearman 0.604), PABP (Pearson 0.577; Spearman 0.584) and TEM1
(Pearson 0.464; Spearman 0.493) sat in an intermediate band. None of
the per-protein coefficients reached the magnitude expected if SCI
were a near-linear re-expression of LLR. We then compared the two
scores on the binarised DMS label
(`sci_vs_llr_comparison.csv`): LLR alone reached AUROC = 0.7531 at
the pooled level, SCI alone reached AUROC = 0.6685, and a 5-fold
cross-validated logistic combination of [SCI, LLR] reached
AUROC = 0.7578. The combined model improved over LLR-alone by
ΔAUROC ≈ 0.005, consistent with SCI contributing a small but
non-zero increment on top of the LLR ranking rather than displacing
it. The combined-model gain co-occurs with the partial Pearson /
Spearman overlap reported above, with no per-protein scope showing a
correlation high enough to make the two scores interchangeable.

## Discussion

Two findings constrain the interpretation of SCI relative to LLR.
First, SCI and LLR are only partially correlated (Pooled
Pearson 0.456 / Spearman 0.467), with per-protein coefficients
ranging from 0.185 (PTEN) to 0.604 (GAL4, Spearman). These values
are too low for SCI to be treated as a simple re-expression of LLR:
the two scores share rank information but explain substantively
different residual variance, especially on PTEN and CBS where the
shared signal is weakest. Second, the incremental gain from adding
SCI to LLR in a cross-validated logistic combination is small
(ΔAUROC ≈ 0.005 over LLR-alone). We report this increment without
emphasis, because at the pooled scale and with N = 95,142 rows a
ΔAUROC of this size is best read as "non-zero, not large" rather
than as evidence that SCI is a competitive standalone scorer; LLR
alone (AUROC = 0.7531) remains the stronger univariate baseline,
and SCI alone (AUROC = 0.6685) does not match it. The value of SCI
in this context is therefore not framed as classification accuracy.
SCI is a layer-pair-geometry score: it summarises how the residue's
embedding is repositioned across selected ESM-2 transformer layer
pairs, a representation-internal quantity that the scalar LLR cannot
report by construction. The partial correlation pattern is
consistent with this framing — the two scores agree where the
geometric shift tracks the masked-marginal likelihood and diverge
where it does not. Treating SCI as complementary rather than
competitive is the interpretation supported by these numbers.
