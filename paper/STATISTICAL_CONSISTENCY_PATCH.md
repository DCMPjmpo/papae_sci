# STATISTICAL CONSISTENCY PATCH

**Date:** 2026-07-02  
**Scope:** Correct 3 statistical description errors across manuscript files  
**Status:** READ-ONLY — awaiting approval before applying  

---

## Key Facts (Verified from CSV)

| Protein | Statistic | p_empirical | p_bonferroni | q (BH-FDR) | Bonferroni sig? | FDR sig? |
|---------|-----------|-------------|--------------|------------|----------------|----------|
| CBS | T1 | 0.0300 | **0.599** | **0.037** | ❌ No | ✅ Yes |
| PTEN | T1 | 0.00699 | **0.140** | **0.009** | ❌ No | ✅ Yes |
| All other 13 | T1/T3/T3std | ≤0.001 | ≤0.02 | ≤0.0014 | ✅ Yes | ✅ Yes |

---

## Error Locations (3 errors, 3 files each = 9 edits)

### Error 1 — §2.6: "all 15 tests retained Bonferroni significance" [WRONG]

**Fact:** CBS T1 (p_bonf=0.599) and PTEN T1 (p_bonf=0.140) do NOT survive Bonferroni.

#### Location 1A: `paper/final_draft_submission.md` lines 276–281

```
After Bonferroni correction over the 20 per-protein primary tests (5 proteins × 4
statistics), each of these 15 (T₁, T₃, T₃ˢᵗᵈ across 5 proteins) tests retains
significance at the conventional α = 0.05 level.
```

**Replace with:**

```
After Bonferroni correction over the 20 per-protein primary tests (5 proteins × 4
statistics), 13 of these 15 tests (T₁, T₃, T₃ˢᵗᵈ across 5 proteins) retained
significance at the conventional α = 0.05 level; CBS-T₁ (p_bonf = 0.60) and
PTEN-T₁ (p_bonf = 0.14) did not survive the conservative Bonferroni threshold
despite being significant under BH-FDR control (q = 0.037 and q = 0.009,
respectively). The signal is therefore not pooled out of one over-represented
protein; it is reproduced within each, with strongest consistency for T₃ and
T₃ˢᵗᵈ across all five proteins.
```

---

#### Location 1B: `paper/final_draft.md` lines 388–393

Same text as 1A. Apply identical replacement.

---

#### Location 1C: `paper/manuscript_first_draft.md` lines 298–303

Same text as 1A. Apply identical replacement.

---

### Error 2 — §2.7: "Bonferroni and BH-FDR calls match exactly" [WRONG]

**Fact:** CBS T1 and PTEN T1 differ (Bonferroni: NS, BH-FDR: sig). 2 of 15 mismatches.

#### Location 2A: `paper/final_draft_submission.md` lines 289–294

```
(ii) Bonferroni correction over the 20-test family of per-protein primary
statistics yielded Bonferroni significance at α = 0.05 for T₁, T₃, T₃ˢᵗᵈ in
5/5 proteins. (iii) Benjamini–Hochberg FDR with step-up monotonicity enforcement
was applied to the same 20-test family; q < 0.05 calls match the Bonferroni
calls exactly on T₁, T₃, T₃ˢᵗᵈ across the five proteins.
```

**Replace with:**

```
(ii) Bonferroni correction over the 20-test family of per-protein primary
statistics retained significance at α = 0.05 for T₃ and T₃ˢᵗᵈ in all five
proteins; for T₁, three of five proteins (GAL4, PABP, TEM1) survived Bonferroni,
while CBS-T₁ and PTEN-T₁ did not (p_bonf = 0.60 and 0.14, respectively).
(iii) Benjamini–Hochberg FDR with step-up monotonicity enforcement was applied
to the same 20-test family, yielding q < 0.05 for T₁ across all five proteins
and thus providing a slightly less conservative set of discoveries than
Bonferroni. All 15 tests remained significant under BH-FDR control.
```

---

#### Location 2B: `paper/final_draft.md` lines 402–409

```
(ii) Bonferroni correction over the 20-test family of per-protein primary
statistics yielded `p_bonferroni` columns in `permutation_per_protein_summary.csv`;
Bonferroni significance at α = 0.05 was retained for T₁, T₃, T₃ˢᵗᵈ in 5/5
proteins. (iii) Benjamini–Hochberg FDR with step-up monotonicity enforcement
was applied to the same 20-test family, with the `p_fdr` column containing the
BH q-values; q < 0.05 calls match the Bonferroni calls exactly on T₁, T₃, T₃ˢᵗᵈ
across the five proteins.
```

**Replace with:**

```
(ii) Bonferroni correction over the 20-test family of per-protein primary
statistics yielded `p_bonferroni` columns in `permutation_per_protein_summary.csv`;
Bonferroni significance at α = 0.05 was retained for T₃ and T₃ˢᵗᵈ in all five
proteins, and for T₁ in three of five proteins (GAL4, PABP, TEM1). CBS-T₁
(p_bonf = 0.60) and PTEN-T₁ (p_bonf = 0.14) did not survive Bonferroni.
(iii) Benjamini–Hochberg FDR with step-up monotonicity enforcement was applied
to the same 20-test family, with the `p_fdr` column containing the BH q-values;
BH-FDR yielded q < 0.05 for T₁ in all five proteins, providing slightly broader
discovery than Bonferroni while still controlling the expected false-discovery
rate among the 20-test family.
```

---

#### Location 2C: `paper/manuscript_first_draft.md` lines 312–319

Same as 2B. Apply identical replacement.

---

### Error 3 — §4.9: "three sanity-check resamples" [INACCURATE]

**Fact:** The code `_sanity_check_rank_reuse_bootstrap` defaults to `n_check=10`, and the separate full validation `bootstrap_validation_exact_vs_fast` runs 100 resamples. Saying "three" is a significant undercount.

#### Location 3A: `paper/final_draft_submission.md` lines 581–584

```
three sanity-check resamples
```

**Replace with:**

```
multiple sanity-check resamples
```

---

#### Location 3B: `paper/final_draft.md` lines 767–768

Same as 3A. Apply identical replacement.

---

#### Location 3C: `paper/manuscript_first_draft.md` lines 679–680

Same as 3A. Apply identical replacement.

---

### Ancillary Locations — Abstract, Introduction, and Take-Home Messages

These make broader claims about "all effects survive Bonferroni" that need qualification.

#### Ancillary-1: Abstract (all 3 files)

Current:
```
All primary effects survive Bonferroni FWER, Benjamini–Hochberg FDR, and
non-parametric bootstrap.
```

Recommended replacement:
```
All primary effects survive Benjamini–Hochberg FDR and non-parametric bootstrap;
T₃ and T₃ˢᵗᵈ also survive the more conservative Bonferroni FWER across all five
proteins.
```

#### Ancillary-2: Introduction ¶5 (all 3 files)

Current:
```
Every primary effect survives Bonferroni FWER, Benjamini–Hochberg FDR over
twenty per-protein primary tests, and non-parametric bootstrap.
```

Recommended replacement:
```
Every primary effect survives Benjamini–Hochberg FDR over twenty per-protein
primary tests and non-parametric bootstrap; T₃ and T₃ˢᵗᵈ further survive the
more conservative Bonferroni FWER across all five proteins.
```

#### Ancillary-3: Figure 1 Take-home (all 3 files)

Current:
```
reproducible in every one of the five proteins after Bonferroni and BH FDR
correction.
```

Recommended replacement:
```
reproducible in every one of the five proteins after BH FDR correction; T₃ and
T₃ˢᵗᵈ further survive Bonferroni correction across all five proteins.
```

#### Ancillary-4: Figure 5 Take-home (all 3 files)

Current:
```
Every primary finding survives three independent layers of statistical robustness
— Phipson–Smyth permutation, Bonferroni FWER, Benjamini–Hochberg FDR — and
the bootstrap 95 % CI excludes the null mean for every primary statistic;
per-protein reproducibility is established in 5 / 5 proteins on T₁, T₃, T₃ˢᵗᵈ.
```

Recommended replacement:
```
Every primary finding survives Phipson–Smyth permutation control and the
bootstrap 95 % CI excludes the null mean for every primary statistic. Per-protein
BH-FDR reproducibility is established in 5/5 proteins on T₁, T₃ and T₃ˢᵗᵈ;
T₃ and T₃ˢᵗᵈ additionally survive Bonferroni FWER across all five proteins.
```

---

## Summary of All Required Edits

| # | File | Section | Original claim | Severity | Fix |
|---|------|---------|----------------|----------|-----|
| 1a | `final_draft_submission.md` | §2.6 | "all 15 tests retain Bonferroni significance" | **Critical** | Qualify: 13/15 pass, CBS+PTEN T1 do not |
| 1b | `final_draft.md` | §2.6 | Same | **Critical** | Same as 1a |
| 1c | `manuscript_first_draft.md` | §2.6 | Same | **Critical** | Same as 1a |
| 2a | `final_draft_submission.md` | §2.7 | "BH-FDR calls match Bonferroni exactly" | **Critical** | Describe the divergence honestly |
| 2b | `final_draft.md` | §2.7 | Same | **Critical** | Same as 2a |
| 2c | `manuscript_first_draft.md` | §2.7 | Same | **Critical** | Same as 2a |
| 3a | `final_draft_submission.md` | §4.9 | "three sanity-check resamples" | Minor | "multiple sanity-check resamples" |
| 3b | `final_draft.md` | §4.9 | Same | Minor | Same as 3a |
| 3c | `manuscript_first_draft.md` | §4.9 | Same | Minor | Same as 3a |
| A1 | Abstract (all 3) | Abstract | "All effects survive Bonferroni" | **Major** | Qualify: only T3/T3std survive Bonferroni |
| A2 | Introduction ¶5 (all 3) | §1 ¶5 | Same | **Major** | Same as A1 |
| A3 | Fig1 take-home (all 3) | Figure legends | "Bonferroni and BH FDR" | Minor | Separate Bonferroni from FDR claims |
| A4 | Fig5 take-home (all 3) | Figure legends | "Bonferroni FWER" in blanket statement | Minor | Narrow to T₃/T₃ˢᵗᵈ |

**Files requiring edits:** 3 (all manuscript files)  
**Total distinct text blocks to replace:** 12 (3 core errors × 3 files + 3 ancillary patterns)  

**No statistical code, data files, figures, or tables are modified.**
