# Reviewer Defense — SCI vs ESM-2 zero-shot LLR baseline

> Companion note to `paper/p7_correlation_summary.md`. All numbers are
> sourced from `data/processed/sci_vs_llr_comparison.csv` and
> `data/processed/sci_llr_correlation.csv`. Nothing in this document
> revises or adds to the headline numerical results — it only
> consolidates the response to one specific anticipated critique.

---

## 1. Reviewer concern

A reviewer might raise the following objection:

> "The proposed SCI score reaches AUROC = 0.6685 at the pooled
> level, while the canonical ESM-2 masked-marginal LLR baseline
> (Meier et al. 2021) reaches AUROC = 0.7531 on the same n = 95,142
> mutations. Adding SCI to LLR yields only ΔAUROC ≈ 0.005. Since
> SCI is the weaker univariate scorer and the gain from combining
> the two is marginal, the practical value of SCI is limited."

The concern is factually correct on the numerical comparison and we
do not contest it. The response below clarifies the scope in which
SCI is presented, and the reasons we still consider the layer-pair
geometry result a substantive contribution rather than a competing
predictor of fitness.

## 2. Response

We acknowledge that, judged purely as a zero-shot fitness scorer,
LLR is the stronger univariate baseline on this benchmark
(AUROC 0.7531 vs SCI's 0.6685) and that the combined logistic model
gains only ΔAUROC ≈ 0.005 over LLR alone. We do not claim that SCI
improves variant-effect prediction beyond what LLR already provides.
SCI is framed in this manuscript as a tool for characterising the
internal layer-pair structure of ESM-2 in the context of
single-residue perturbations, not as a replacement scoring function
for fitness ranking. The numerical observation that SCI and LLR are
only partially correlated (Pooled Pearson 0.456 / Spearman 0.467;
PTEN as low as Pearson 0.185 / Spearman 0.214) supports this
framing: the two quantities share rank information but are not
interchangeable, and the value of SCI is in what it reports about
representational geometry, not in displacing LLR on AUROC.

## 3. Why SCI is not intended as a replacement for LLR

LLR is, by construction, the score the language model itself
assigns to a substitution: log P(mut) − log P(wt) under the masked
position. It directly answers the question "how unlikely is this
substitution under the model." SCI answers a different question:
"how is the residue's embedding repositioned, across selected pairs
of transformer layers, when the substitution is introduced." These
are different observables on the same forward pass, and the
manuscript treats them as such. Replacing LLR with SCI as a fitness
predictor was never the design goal; the n = 95,142-row comparison
is included specifically to make the relationship between the two
quantitatively explicit (Section 2 above; `sci_vs_llr_comparison.csv`).

## 4. Why layer-pair geometry remains valuable

The SCI score is constructed over pairs of ESM-2 layers and is the
substrate on which the layer-pair selection, recurrence, and
Early–Late enrichment analyses in this manuscript are built
(see Results §2.2–§2.5). These analyses report which layer pairs
carry mutation-discriminative geometry, how consistent that
selection is across proteins, and how the selected pairs distribute
over the early/mid/late depth bins. LLR, being a single scalar per
substitution derived from the final-layer LM head, does not expose
any of these layer-resolved quantities. The layer-pair geometry
result is therefore not a re-statement of LLR with extra steps; it
is an analysis layer that LLR neither provides nor competes against.

## 5. What information SCI provides that LLR cannot

By construction, LLR collapses the model's response to a substitution
into one scalar at one position. SCI retains a per-(layer i, layer j)
matrix per residue before the headline summary, which yields three
quantities used in this manuscript:

- **Layer-pair selection** — which (i, j) pairs concentrate the
  mutation-effect signal (Results §2.2; `pca_loadings.csv`,
  `layer_recurrence_top20.csv`).
- **Depth structure** — Early–Late enrichment of the top-ranked
  pairs across the five proteins (Results §2.3).
- **Per-protein geometric consistency** — pair-overlap and
  recurrence patterns across protein boundaries
  (Results §2.4–§2.5).

None of these can be read off the LLR scalar. The partial-correlation
pattern in Section 2 (low on PTEN/CBS, intermediate on TEM1,
moderate on GAL4/PABP) is consistent with SCI and LLR tracking
overlapping but non-coextensive aspects of the same forward pass.

## 6. Limitations

We list the limitations explicitly so the framing above does not
overreach.

- SCI alone is a weaker pooled fitness classifier than LLR on this
  benchmark (AUROC 0.6685 vs 0.7531), and we do not suggest
  otherwise.
- The combined-model gain (ΔAUROC ≈ 0.005) is small relative to
  the sample size (n = 95,142); we do not interpret it as evidence
  that SCI is a competitive standalone scorer.
- The correlation between SCI and LLR is uneven across proteins.
  PTEN in particular shows the lowest agreement (Pearson 0.185;
  Spearman 0.214); we do not present this as a uniform property of
  SCI across proteomes.
- All numerical claims are restricted to the five proteins in this
  manuscript (CBS, GAL4, PABP, PTEN, TEM1) and to ESM-2
  esm2_t33_650M_UR50D; the layer-pair geometry result is not
  asserted to transfer to other PLM backbones without re-validation.
- The analyses in this manuscript do not make causal claims about
  the role of any individual layer pair.

## 7. Future directions

The natural next steps follow from the limitations rather than from
the AUROC gap.

- Repeat the layer-pair geometry analysis on additional PLM
  backbones (e.g. larger ESM-2 variants, ESM-3, ProtTrans family)
  to test whether the depth-enrichment pattern is backbone-specific
  or recurs across architectures.
- Examine whether the (i, j) pairs selected by SCI correspond to
  layer transitions that have been independently characterised in
  the interpretability literature (probing studies of secondary
  structure, contact prediction, evolutionary conservation).
- Extend the per-protein comparison to a wider DMS panel
  (e.g. additional ProteinGym substitutions) to test whether the
  PTEN-style low-correlation regime is associated with specific
  protein-level properties (size, fold class, depth of MSA).
- Investigate per-residue and per-region SCI rather than the single
  scalar headline summary used here, to expose layer-pair structure
  at sub-protein resolution.

None of these directions require SCI to surpass LLR on AUROC; they
ask what the layer-pair representation itself looks like once the
question is no longer "is this substitution rare under the model."
