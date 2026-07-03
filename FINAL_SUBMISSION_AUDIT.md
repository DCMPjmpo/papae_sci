# FINAL SUBMISSION AUDIT REPORT

**Date:** 2026-07-02  
**Scope:** Pre-submission audit of ESM-2 layer-pair manuscript  
**Auditor:** Claude Code (audit mode — no files modified)

---

## A. Model Name Consistency

### Verdict: PASS — no incorrect model name in manuscript text or code

| File | Finding | Verdict |
|------|---------|---------|
| `paper/final_draft_submission.md` | All references to "ESM-2 (esm2_t33_650M_UR50D)" are correct as our model | ✅ |
| `paper/final_draft_submission.md` | "ESM-1b" appears 4 times, all correctly referencing prior work (Brandes et al. 2023, Rives et al. 2021) as citations | ✅ Appropriate |
| `paper/final_draft_submission.md` | "ESM-1v" appears 1 time citing Meier et al. 2021 (the LLR baseline) | ✅ Appropriate |
| `paper/final_draft_submission.md` | "ESM3" appears 4 times as future-work replication target | ✅ Appropriate |
| `paper/final_draft_submission.md` | "ProtT5" appears 2 times as future-work replication target | ✅ Appropriate |
| All `.py` scripts | No Python script labels any axis "ESM-1b". All plot labels use "ESM-2 layer index (1-33)" | ✅ |
| Panel PNG images | All x-axis labels say "ESM-2 layer index (1-33)" (verified in `layer_recurrence_top20.py:132`, `layer_occurrence_frequency.py:121`) | ✅ |
| `paper/discussion_outline_v2.md` | Line 38: "ESM-2 vs ESM-2 vs ESM3 vs ProtT5" — duplicate "ESM-2" (should be "ESM-1b vs ESM-2 vs ESM3 vs ProtT5") | ⚠️ Minor (outline file, not submission) |
| `paper/literature_mapping.md` line 30 | "ESM-1" refers to Rives et al. 2021's original ESM-1 work | ✅ Appropriate |

**Summary:** Zero instances of our model being called the wrong name. All "ESM-1b" references correctly cite prior work. No further action required.

---

## B. Figure Consistency

### B1. Main Figures (Fig 1–5)

| Check | Finding | Verdict |
|-------|---------|---------|
| Fig 1 cited in text | Lines 151, 1000 | ✅ |
| Fig 2 cited in text | Lines 174, 179, 202, 1001 | ✅ |
| Fig 3 cited in text | Lines 223, 225, 231, 975, 1005, 1006 | ✅ |
| Fig 4 cited in text | Lines 251, 259 | ✅ |
| Fig 5 cited in text | Lines 272, 304, 1003, 1004 | ✅ |
| Citation order | Fig 1 → Fig 2 → Fig 3 → Fig 4 → Fig 5 | ✅ Correct |
| Panel labels match legend | All (a)/(b)/(c) labels in legends match text descriptions | ✅ |
| Fig 1 panel (a) text reference | Never explicitly called "Fig 1a" in text — only "Fig 1b" appears (line 151) | ⚠️ Minor — panel (a) is described but not cross-referenced by label |
| Figure 2 legend "z = 8.0" (line 695) vs Results "z = 8.01" (line 179) | **Inconsistency:** Fig 2b legend truncates z = 8.01 to 8.0, while Results text uses 8.01. CSV source says 8.01. | ⚠️ Minor — inconsistent rounding |
| Figure 2 legend "mean 1.0" (line 695) vs Results "mean of 1.01" (line 179) | **Inconsistency:** Fig 2b legend rounds 1.01 to 1.0 | ⚠️ Minor — inconsistent rounding |
| Table 1 in manuscript | Numbers from `sci_vs_llr_comparison.csv` — see Section C below | See below |

### B2. Supplementary Figures (S1–S5)

| Figure | File exists | Proper composite | Notes |
|--------|------------|-----------------|-------|
| S1 | ✅ `figures/supplementary/Supp_Fig_S1.png` (400 KB) | ✅ Newly created | Per-protein SCI histograms |
| S2 | ✅ `figures/supplementary/Supp_Fig_S2.png` (4.2 MB) | ✅ Newly created | Two-panel composite |
| S3 | ✅ `figures/supplementary/Supp_Fig_S3.png` (3.4 MB) | ✅ Newly created | Matrix heatmap + best-pair dist |
| S4 | ✅ `figures/supplementary/Supp_Fig_S4.png` (995 KB) | ✅ Identical to Fig 5a | Per-protein null grid |
| S5 | ✅ `figures/supplementary/Supp_Fig_S5.png` (432 KB) | ✅ Identical to Fig 5c | QQ significance plot |

### B3. Figure Legend Cross-Reference Issues in Supplementary

| Location | Reference | Should Be |
|----------|-----------|-----------|
| `paper/supplementary.md` line 136 | `(Methods §3.3)` | `(Methods §4.3)` |
| `paper/supplementary.md` line 163 | `(Methods §3.6)` | `(Methods §4.6)` |
| `paper/supplementary.md` line 180 | `(Methods §3.6 note)` | `(Methods §4.6 note)` |
| `paper/supplementary.md` line 240 | `Methods §3.1` | `Methods §4.1` |
| `paper/supplementary.md` line 241 | `Methods §3.7–3.10` | `Methods §4.7–4.10` |
| `paper/supplementary.md` line 244 | `Methods §3.3` | `Methods §4.3` |
| `paper/supplementary.md` line 28 | `(Methods §3.7)` | `(Methods §4.7)` |
| `paper/supplementary.md` line 59 | `(see Methods §3.6 note)` | `(see Methods §4.6 note)` |
| `paper/supplementary.md` line 75 | `(§3.6 note)` | `(§4.6 note)` |

**🔴 CRITICAL:** The entire `paper/supplementary.md` file uses Methods §3.x numbering from an earlier draft phase, while the main manuscript (`final_draft_submission.md`) uses Methods §4.x. Every cross-reference in the supplementary materials is broken.

---

## C. Numerical Consistency

### C1. Global Statistics: manuscript vs CSV source-of-truth

| Statistic | Manuscript value | CSV value | Match? |
|-----------|-----------------|-----------|--------|
| T₁ observed | 2.79 × 10⁻³ | 0.00279397 | ✅ |
| T₁ z-score | 57.0 | 57.01 | ✅ (rounding) |
| T₁ p_FDR | 10⁻³ | 0.000999 | ✅ |
| T₁ bootstrap 95% CI | [2.67, 2.91] × 10⁻³ | [0.00267356, 0.00290996] | ✅ |
| T₂ observed | 17 | 17 | ✅ |
| T₂ z-score | 8.01 (text) / 8.0 (legend) | 8.01 | ⚠️ |
| T₂ p_FDR | 10⁻³ | 0.000999 | ✅ |
| T₂ bootstrap CI | [14, 19] | [14, 19] | ✅ |
| T₃ observed | 4.0 × 10⁻² | 0.0398858 | ✅ (3.99 × 10⁻²) |
| T₃ z-score | 62.0 | 62.03 | ✅ (rounding) |
| T₃ p_FDR | 10⁻³ | 0.000999 | ✅ |
| T₃ bootstrap 95% CI | [3.70, 4.27] × 10⁻² | [0.0370, 0.0427] | ✅ |
| T₃ˢᵗᵈ observed | 1.02 | 1.01553 | ✅ (rounding) |
| T₃ˢᵗᵈ z-score | 8.05 | 8.05 | ✅ |
| T_max z-score | 93.7 | 93.72 | ✅ |

### C2. Per-Protein z-scores: manuscript vs CSV

| Protein | T₁ z MS | T₁ z CSV | T₃ z MS | T₃ z CSV | T₃ˢᵗᵈ z MS | T₃ˢᵗᵈ z CSV |
|---------|---------|----------|---------|----------|-------------|--------------|
| CBS | 2.25 | 2.25 ✅ | 16.11 | 16.11 ✅ | 18.50 | 18.50 ✅ |
| GAL4 | 6.29 | 6.29 ✅ | 5.18 | 5.18 ✅ | 4.59 | 4.59 ✅ |
| PABP | 56.55 | 56.55 ✅ | 38.71 | 38.71 ✅ | 42.80 | 42.80 ✅ |
| PTEN | 2.88 | 2.88 ✅ | 14.57 | 14.57 ✅ | 17.05 | 17.05 ✅ |
| TEM1 | 23.59 | 23.59 ✅ | 15.19 | 15.19 ✅ | 15.94 | 15.94 ✅ |

### C3. Table 1: Three Rounding Discrepancies (0.001)

| Row | Column | CSV value | Manuscript value | Diff | Severity |
|-----|--------|-----------|------------------|------|----------|
| Pooled SCI | AUROC | **0.668** | **0.669** | 0.001 | ⚠️ Minor — misrounded (0.668 → 0.669) |
| GAL4 SCI | AUROC | **0.718** | **0.717** | 0.001 | ⚠️ Minor — misrounded (0.718 → 0.717) |
| PABP LLR | Spearman r | **0.485** | **0.486** | 0.001 | ⚠️ Minor — misrounded (0.485 → 0.486) |

### C4. Figure Legend rounding inconsistency

| Figure 2b | Results text (line 179) | Figure legend (line 695) | CSV |
|-----------|------------------------|-------------------------|-----|
| T₂ z-score | z = **8.01** | z = **8.0** | 8.01 |
| Null mean | **1.01** | **1.0** | 1.013 |

The figure legend truncates precision compared to the Results text and CSV. The legend should match the text.

**Numerical verdict:** All primary results (z-scores, p-values, CIs) are correct. Three Table 1 values are off by ±0.001 due to rounding direction. One legend vs text inconsistency on z = 8.0/8.01.

---

## D. Table Consistency (Placeholders & Placeholder Text)

### D1. Main Manuscript (`paper/final_draft_submission.md`)

| Line | Content | Issue |
|------|---------|-------|
| 871 | `T₂: n/a` (CBS, GAL4, PTEN, TEM1) | ✅ Acceptable — T₂ is undefined for these proteins (documented in Methods §4.6) |
| 873 | `(T_obs=1; see Methods §4.6 note)` (PABP T₂) | ✅ Acceptable — intentional methodological note |
| — | No TODO, TBD, XXX, TBA, placeholder found | ✅ Clean |

### D2. Supplementary File (`paper/supplementary.md`)

| Line | Content | Issue |
|------|---------|-------|
| 58 | GAL4 T_max z: **(see CSV)** | ❌ Placeholder — needs actual value from CSV |
| 59 | PABP T_max z: **(see CSV)** | ❌ Placeholder — needs actual value from CSV |
| 60 | PTEN T_max z: **(see CSV)** | ❌ Placeholder — needs actual value from CSV |
| 61 | TEM1 T_max z: **(see CSV)** | ❌ Placeholder — needs actual value from CSV |

The T_max z-scores are available in `permutation_per_protein_summary.csv` (column `z_score` for `Statistic == "T_max"`). These should replace `(see CSV)`.

### D3. Reference List (`paper/final_draft.md`)

Previously had placeholder at lines 1412–1444. **Now fixed** — replaced with complete numbered reference list from `final_draft_submission.md`. ✅

### D4. Section Numbering Mismatch

As documented in Section B3 above, `paper/supplementary.md` uses **9 instances** of `§3.x` that should be `§4.x` to match the main manuscript. This affects every cross-reference between supplementary and methods.

---

## E. Submission-Blocking Issues

### 🔴 CRITICAL (Submit-blocking if unresolved)

| # | Issue | Location | Fix |
|---|-------|----------|-----|
| **C1** | `paper/supplementary.md` — all 9 Methods cross-references use §3.x numbering instead of §4.x | Lines 28, 59, 75, 136, 163, 180, 240, 241, 244 | Replace `§3` with `§4` throughout; replace `§3.7–3.10` with `§4.7–4.10` |
| **C2** | `paper/supplementary.md` — 4 T_max values reported as "(see CSV)" | Lines 58, 59, 60, 61 | Fill in actual z-scores from `permutation_per_protein_summary.csv` |
| **C3** | Main composite figures (Fig1–Fig5) still not assembled as single publication files | `figures/` | Run `figure_assembly/assemble_figures.py` to generate `Fig1.tiff`–`Fig5.tiff` |

### ⚠️ MAJOR (Should fix before submission)

| # | Issue | Location | Fix |
|---|-------|----------|-----|
| **M1** | `paper/supplementary.md` cross-reference index says "Methods §3.1" but main manuscript Methods section is §4 | Line 240 | Fix to §4.1 |
| **M2** | Table 1: 3 values off by ±0.001 from source CSV | Lines 807, 815, 819 | Correct rounding: Pooled SCI AUROC = 0.668, GAL4 SCI AUROC = 0.718, PABP LLR Spearman r = 0.485 |
| **M3** | Fig 2b legend reports z = 8.0 and null mean = 1.0, but Results text and CSV say z = 8.01 and null mean = 1.01 | Line 695 vs line 179 | Align legend with text (z = 8.01, mean 1.01) |
| **M4** | `paper/final_draft_submission.md` supplementary section (lines 867–875) cross-references "Methods §4.6 note" correctly, but `paper/supplementary.md` (line 59) uses "Methods §3.6 note" | Both files | Fix supplementary.md to use §4.6 |

### 🟡 MINOR (Fix when convenient)

| # | Issue | Location | Severity |
|---|-------|----------|----------|
| **N1** | Fig 1 panel (a) never explicitly referenced by label in text — only "Fig 1b" appears | Line 151 | Low — many journals accept implicit panel references |
| **N2** | `paper/discussion_outline_v2.md` line 38: "ESM-2 vs ESM-2" should be "ESM-1b vs ESM-2" | Outline file, not submission | Negligible |
| **N3** | `data/p0_5layer/` contains 3 legacy figures from old numbering (figure5, figure6, figure7) that conflict with current Fig 1–5 scheme | `figure5_independent_protein_validation.png`, `figure6_p4_no_leakage.png`, `figure7_p5_statistical_significance.png` | Low — archive these files before submission |
| **N4** | Reference list: DOI and full venue names absent | `final_draft_submission.md` lines 1012–1041 | Add DOIs at typesetting |
| **N5** | `final_draft_submission.md` Figure 1 legend (line 922) says "overlaid histograms" but source panel `P1.5_sci_signal_validation.png` panel A contains boxplots | Legend vs actual image | Minor — either update legend to say "boxplots" or regenerate panel with histograms |

---

## Summary

| Category | Count | Breakdown |
|----------|-------|-----------|
| ✅ Pass | 5 | Model name consistency, Figure citation order, Per-protein z-scores, Bootstrap CIs, Global p-values |
| 🔴 Critical | 3 | Supplementary §3.x→§4.x renumbering (9 instances), T_max placeholders, Main composites not assembled |
| ⚠️ Major | 4 | Table 1 rounding (3 cells), Fig 2b z-score inconsistency, Cross-reference numbering mismatch |
| 🟡 Minor | 5 | Fig 1a label, Outline typo, Legacy files, DOI missing, Legend/image mismatch |

**Overall readiness:** Manuscript text is substantially complete. Supplementary cross-references need full renumbering. Three figures (Fig1–Fig5 composites) still require assembly. Numerically, 0 wrong results — only 3 cells in Table 1 off by ±0.001.
