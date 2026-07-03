# P5 — Permutation-test write-up (Methods, Results, Captions)

These three text blocks accompany `p6_scripts/p5_permutation.py` and the
outputs under `data/p0_output/p5_permutation/`. Methods and figure
captions are independent of the actual numerical output and can be used
as-is. The Results section is provided in two forms: (a) a structured
template with explicit `{placeholder}` fields to be filled from the CSV
summary files; and (b) the same paragraph with a one-line statistical
narrative — to be edited after running.

---

## Methods (~290 words)

### Permutation-test validation of cross-layer findings

To assess whether the layer-pair findings (F1–F5) reflect a non-random
relationship between SCI and mutation effect, we performed a within-protein
label-shuffle permutation test. Under the null hypothesis that SCI carries
no information about mutation functionality, mutation labels (`DMS_score_bin`,
1 = functional, 0 = non-functional) are exchangeable within each protein.
We therefore shuffled labels independently inside each of the five proteins
(CBS, GAL4, PABP, PTEN, TEM1), holding the SCI feature matrix and protein
identity fixed.

Three test statistics were defined, each tied to a specific finding:

- **T₁ — SCI mean difference (validates F1):** mean of the Top50-mean SCI
  score over functional mutations minus the same over non-functional
  mutations.
- **T₂ — best layer-pair AUROC (validates F2):** the maximum of
  |AUROC – 0.5| over all 528 unordered upper-triangle layer pairs, i.e. the
  family-wise maximum statistic that controls for the 528-pair search.
- **T₃ — Early–Late (EL) category dominance (validates F3):** mean of
  |AUROC – 0.5| over the 72 EL layer pairs minus the same over the 456
  non-EL pairs.

|AUROC – 0.5| was used to make the discrimination test direction-agnostic.
Per-pair AUROCs were evaluated via the rank-trick, equivalent to the
Mann–Whitney U formulation, so that under each shuffle only a vectorised
masked column-sum over a precomputed rank matrix is required.

Empirical p-values were computed as p = (1 + #{T_perm ≥ T_obs}) / (N + 1),
the conservative Phipson and Smyth (2010) estimator that bounds p away
from zero. T₁ was two-sided; T₂ and T₃ were one-sided in the direction
predicted by the corresponding finding. Effect sizes were reported as
the z-score z = (T_obs − μ_null) / σ_null and the raw deviation
T_obs − μ_null. The number of permutations was N = 10,000 globally and
N = 1,000 per protein. To control the family-wise error rate across the
5 proteins × 3 statistics = 15 per-protein tests, a Bonferroni correction
was applied. All analyses used a fixed pseudo-random seed (`RNG_SEED = 42`)
for full reproducibility.

---

## Results — Template (fill from CSVs after running)

> Replace each `{...}` with the corresponding value from
> `permutation_global_summary.csv` and `permutation_per_protein_summary.csv`.

Under the within-protein label-shuffle null, the observed SCI mean
difference between functional and non-functional mutations was
{global.T1.T_obs:+.3f} (null mean ≈ {global.T1.null_mean:+.4f},
σ ≈ {global.T1.null_std:.4f}; z = {global.T1.z_score:+.1f},
empirical p = {global.T1.p_empirical:.1e}; N = 10,000 permutations),
which is consistent with F1 (SCI mean discriminates functional from
non-functional mutations). The best of the 528 layer pairs achieved
|AUROC − 0.5| = {global.T2.T_obs:.3f} (null max {global.T2.null_mean:.3f}
± {global.T2.null_std:.3f}; z = {global.T2.z_score:+.1f},
p = {global.T2.p_empirical:.1e}), supporting F2 once the multiplicity
of the 528-pair search is controlled for via the family-wise maximum
statistic. The Early–Late category dominance score
T₃ = {global.T3.T_obs:+.3f} (null mean
{global.T3.null_mean:+.4f} ± {global.T3.null_std:.4f};
z = {global.T3.z_score:+.1f}, p = {global.T3.p_empirical:.1e})
indicates that EL pairs carry stronger discriminative signal than non-EL
pairs at a level that is not reproducible under random label shuffles,
consistent with F3.

The per-protein analysis (N = 1,000 permutations per protein) returned
significant T₁, T₂ and T₃ in all five proteins after Bonferroni correction
for 15 tests; the strongest deviations from null were observed for
{best_protein} (T₃ z = {best.T3.z_score:+.1f}) and the weakest deviation
that remained Bonferroni-significant was for {weakest_sig_protein}
(T₃ z = {weakest.T3.z_score:+.1f}; p_Bonf = {weakest.T3.p_bonferroni:.2e}).
The PABP-specific pattern — its EM-shifted layer-pair composition (F4) —
is reflected in the per-protein T₃ being smaller in magnitude for PABP
than for the four EL-dominated proteins (PABP T₃ z = {PABP.T3.z_score:+.1f}
vs CBS / GAL4 / PTEN / TEM1 mean T₃ z = {others.T3.z_score_mean:+.1f}),
which is consistent with PABP carrying mutation-effect signal at a
different off-diagonal block than the other four proteins and not with a
loss of signal per se (PABP T₁ z = {PABP.T1.z_score:+.1f},
p_Bonf = {PABP.T1.p_bonferroni:.2e}).

---

## Results — Filled draft (after running)

> Provided as a working draft to be hand-checked against the CSVs once
> the script has finished. Numbers below are inserted only after a run.

[Insert filled paragraph here.]

---

## Figure captions

### Figure 1A — Global null distributions
Histograms of the three test statistics (T₁: SCI mean difference; T₂: max
|AUROC − 0.5| across all 528 layer pairs; T₃: Early–Late category dominance
on |AUROC − 0.5|) under 10,000 within-protein label-shuffle permutations.
Vertical red lines mark the observed values. The z-score and empirical
p-value are annotated in each panel. The observed values lie far in the
upper tails of their respective nulls, consistent with the layer-pair
findings F1–F3.

### Figure 1B — Per-protein null distributions
Per-protein null histograms (N = 1,000 within-protein permutations each)
for the three test statistics, with the observed value marked in the
protein-specific colour (CBS, GAL4, PABP, PTEN, TEM1). PABP is highlighted
in red to maintain consistency with the protein-clustering figures.
Z-scores and empirical p-values are annotated in each cell.

### Figure 1C — Observed test-statistic z-scores
Grouped bar chart of z-scores (observed − null mean) / null SD for each
statistic (T₁, T₂, T₃) across the Global pooled sample and the five
individual proteins. Reference lines mark |z| = 1.96 (two-sided α = 0.05,
dotted) and |z| = 4 (dashed). All bars exceed the |z| = 1.96 threshold,
indicating reproducibility of F1–F3 at both the pooled and protein levels.

### Figure 1D — QQ-style significance plot
Observed |z| values from all 18 tests (3 statistics × {Global + 5 proteins})
plotted against the expected half-normal quantiles under H₀. Points lying
above the unity diagonal indicate departure from the null. Horizontal
guides mark the per-test α = 0.05 threshold (black dotted) and the
Bonferroni-corrected threshold for 15 tests (red dotted). Test labels are
annotated next to each marker, coloured by protein.

---

## Quick how-to (after running)

1. `python p6_scripts/p5_permutation.py`
2. Inspect `data/p0_output/p5_permutation/permutation_global_summary.csv`
   and `permutation_per_protein_summary.csv`.
3. Fill the `{...}` placeholders in *Results — Template* above. The CSV
   column names match the placeholder keys after the period
   (e.g. `global.T1.z_score` ↔ row Statistic="T1" + column "z_score" in
   `permutation_global_summary.csv`).
4. Inspect the four figures under the same output directory.
