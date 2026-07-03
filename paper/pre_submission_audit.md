# Pre-Submission Audit — `paper/final_draft.md`

> Scope: read-only audit of the assembled manuscript at
> `paper/final_draft.md` against (a) the C1–C7 claim inventory in
> `paper/scientific_story_map.md`, (b) on-disk numerical artefacts
> (`sci_vs_llr_comparison.csv`, `pairwise_euclidean_distance.csv`,
> `permutation_global_summary.csv`,
> `permutation_per_protein_summary.csv`,
> `bootstrap_summary.csv`, `sci_llr_correlation.csv`),
> (c) the user-stated numeric ground truth
> (SCI AUROC 0.6685, LLR AUROC 0.7531, SCI+LLR AUROC 0.7578,
> SCI–LLR Pearson 0.456, Spearman 0.467) and (d) the user's
> forbidden-vocabulary list (prove / demonstrate / cause / drive /
> mechanism / outperform).
>
> The manuscript was NOT modified. This document only reports findings.
> Line numbers refer to `paper/final_draft.md` as of the audit run.

---

## Critical Issues

### CRIT-1 — PABP "singleton condition" is numerically false

**Locations:** Results §2.4 (line 331–333), Figure 3 caption
(line 982–984), Supplementary Note S2 (line 1373–1375).

**Claim made.**
> "the minimum PABP-to-other distance is 0.389 (PABP–TEM1), which
> exceeds the maximum within-other-four distance of 0.397 — i.e.
> PABP's closest neighbour is farther than the most distant pair
> among the other four."

**Actual numbers (from `pairwise_euclidean_distance.csv` /
Sup Note S2 table on line 1367–1371):**
- min PABP-to-other = 0.389 (PABP–TEM1)
- max within-other-four = 0.397 (CBS–GAL4)

**Verdict.** 0.389 < 0.397 (closer by 0.008). The text claims
"0.389 exceeds 0.397", which is arithmetically false; the so-called
**singleton condition is not satisfied** by these values. C4 still
holds qualitatively (PABP is the outlier by every other measure —
heatmap, PC1, dendrogram), but the specific quantitative inequality
is wrong and is repeated in three independent passages. A reviewer
running the comparison will catch this immediately.

**Suggested resolution path** (do not apply here):
- Soften the wording to "of comparable magnitude to" or "approaches"
  rather than "exceeds";
- Or drop the explicit singleton inequality and lean on the
  hierarchical-clustering height (≈ 0.50 vs ≈ 0.22 fusion height),
  which is a real separation; or on PC1 distance (PABP ≈ +0.38 vs
  others ∈ [−0.18, +0.01]), which is real.

### CRIT-2 — Reference-label collision: "F1–F5" in Discussion vs "Fig 1–Fig 5" in Results / Figure Legends

**Locations:** every `F1`/`F2`/`F3`/`F4`/`F5` occurrence in §D1–§D5
(lines 478–518), Methods §3.6 line 650, throughout Figure 2/3/4/5
captions.

**Mismatch.** The Figure Map and the Results section use **Fig 1 …
Fig 5** to label the five main figures, with the mapping:

| Label | Mapping |
|---|---|
| Fig 1 | C1  |
| Fig 2 | C2 + C3 |
| Fig 3 | C4 |
| Fig 4 | C5 |
| Fig 5 | C6 + C7 |

The Discussion (and Methods §3.6) instead uses **F1 … F5** as
*finding* labels with the mapping `Fn = Cn`. For example, D1 line
478 cites "the full 528-pair distribution (F3)" — but C3 lives in
**Fig 2 panel (c)**, not Fig 3. D3 line 492 cites "EL-dominated
within their Top50 Composite layer pairs (F4)" — but C4 (PABP) is
**Fig 3**, not Fig 4. D4 line 500 cites "the narrow deep-layer
anchor at layers 30–32 (F5)" — but C5 is **Fig 4**, not Fig 5.

**Verdict.** A reader following figure callouts from the Discussion
will land on the wrong figure every time the citation involves F3,
F4 or F5. Nature Communications cross-references must resolve
unambiguously; this collision is a typesetting blocker and will be
flagged by any copy-editor.

**Resolution path:** either (a) rename the Discussion's `F1..F5` to
the matching `Fig 1..Fig 5` mapping above, or (b) rename them to the
`C1..C5` claim labels (which the Discussion is actually pointing at).

---

## Major Issues

### MAJ-1 — `driven` appears 4× in body text after a "drive" ban

**Locations:**
- Line 50, Abstract: "PABP forms a singleton **driven** by Early–Middle interactions"
- Line 153, Introduction ¶5: "PABP forms an independent singleton **driven** by Early–Middle interactions"
- Line 378, Results §2.6: "the pooled effects in §2.1–2.3 are **driven** by a single protein" (used in a "rule-out" frame)
- Line 974, Figure 3 caption: "the six-dimensional layer-pair composition space, **driven** by Early–Middle interactions"

The collaborator-notes block (line 198–200) explicitly asserts that
"`drive`" does not appear in the manuscript. As a past participle of
the banned verb "drive", `driven` clearly violates the user's
forbidden-vocabulary list. The Results §2.6 use is the most
defensible (it negates the claim); the Abstract / Intro / Fig 3
uses describe the actual PABP cluster and should be re-phrased
(suggested: "characterised by" / "shifted toward").

### MAJ-2 — Pearson 0.456 / Spearman 0.467 (SCI–LLR partial-correlation result) is missing from the final draft

**Stated ground truth.**
SCI–LLR Pearson r = 0.456 and Spearman ρ = 0.467 at the pooled
level (per `data/processed/sci_llr_correlation.csv` and the companion
`paper/p7_correlation_summary.md` / `paper/reviewer_defense_baseline.md`).

**What the manuscript reports instead.** Results §2.8 / Discussion
D5 report only the *scorer-vs-DMS_score* Spearman values (Pooled
LLR ρ = 0.480, SCI ρ = 0.200, SCI+LLR ρ = 0.464) and argue
redundancy from ΔAUROC ≤ 0.008.

**Verdict.** Two distinct lines of evidence exist for the SCI–LLR
relationship — the *AUROC additivity* test (in the manuscript) and
the *direct rank/linear correlation* test (in the companion docs)
— but only the first is reported. A reviewer familiar with the
companion docs (and the on-disk `sci_llr_correlation.csv`) will
notice that the "partially correlated / partially independent"
framing from `p7_correlation_summary.md` is absent. This is a
coverage gap, not a contradiction; the two analyses are consistent
with each other.

**Resolution path:** either (a) add the SCI–LLR Pearson 0.456 /
Spearman 0.467 one-liner to §2.8 (and a sentence to D5), or (b)
explicitly state in §D5 that the SCI–LLR correlation is not
reported in the main text and cite the supplementary note.

### MAJ-3 — Numeric mismatches between §2.8 prose and Table 1

§2.8 line 453–456 states the per-protein ΔAUROC of SCI+LLR over the
better single scorer as:
- TEM1 Δ = **0.008**
- PABP Δ = +0.003
- PTEN Δ = +0.001
- GAL4 Δ = +0.001
- CBS Δ ≤ 0.001

Recomputing directly from Table 1 (line 1102–1119):

| Protein | LLR | SCI+LLR | Δ (table) | Δ (prose) | Match? |
|---|---:|---:|---:|---:|---|
| TEM1 | 0.882 | 0.889 | **0.007** | 0.008 | off by 0.001 |
| PABP | 0.767 | 0.770 | 0.003 | 0.003 | ✓ |
| PTEN | 0.846 | 0.847 | 0.001 | 0.001 | ✓ |
| GAL4 | 0.795 | 0.797 | **0.002** | 0.001 | off by 0.001 |
| CBS  | 0.692 | 0.692 | 0.000 | ≤ 0.001 | ✓ |

**Verdict.** Two of five per-protein Δ values reported in the prose
do not match the table to the third decimal. Likely caused by
text-rounding of an upstream four-decimal CSV. Should be reconciled
to the table values before submission.

### MAJ-4 — Top-level numbering: Methods uses "§3.X" subsections under top-level "## 4. Methods"

The assembled file uses the top-level numbering 1. Introduction,
2. Results, 3. Discussion, 4. Methods, 5. Figure Legends,
6. Supplementary Information, 7. References. Inside **§4. Methods**,
however, the subsections retain their source-file labels §3.1
through §3.12. Every cross-reference in the manuscript — e.g.
"Methods §3.2", "Methods §3.7", "Methods §3.12" (lines 229, 241,
432, etc.) — resolves correctly *only because* the prose says
"Methods §3.X", not "§3.X" bare. In a Nature Communications PDF
layout this works, but the mixed numbering (top-level 4 with
subsections 3.X) will read as an editing artefact and most
copy-editors will renumber. Decide before submission whether to
renumber Methods subsections to §4.1–§4.12 or to relabel the
top-level Methods header back to "## 3. Methods" (in which case
Discussion needs a non-numeric anchor).

---

## Minor Issues

### MIN-1 — `causal verbs` / `causal claims` phrasing

Line 27 ("without causal verbs"), line 198 ("No causal claims are
introduced"). These are collaborator-facing positioning notes; the
word `cause` per se does not appear as a verb anywhere in the
manuscript's claims. Acceptable for the audit. Consider removing
the meta-block at lines 184–200 from the submitted version (it
explicitly references collaborator-facing markdown files) — these
notes are useful for the author team but should not ship to
reviewers.

### MIN-2 — Per-protein T_max z-scores partially placeholder in Sup Table S2-B

Sup Table S2-B (line 1230–1235) lists per-protein T_max z as
"(see CSV)" for GAL4, PABP, PTEN, TEM1 (CBS is given as −4.20).
This is a fully traceable placeholder, but a reviewer expects the
actual values. Pull them from `permutation_per_protein_summary.csv`
before submission.

### MIN-3 — "the first systematic mapping" claim

Lines 159–161 and 187–193 assert "to our knowledge, the first
systematic mapping". The collaborator-facing notes block (line
188–194) gives a three-point justification and points to
`paper/literature_mapping.md` §6 as the audit. This is the
strongest "novelty" claim in the manuscript and is the most likely
single sentence for a reviewer to push back on. The justification
is defensible (the audited near-neighbour Rao et al. 2021 ICLR
targets contacts, not mutation effects), but the hedge "to our
knowledge" should be retained verbatim in the camera-ready text.

### MIN-4 — `establish` / `establishes` / `establish that` appears 5×

Lines 253, 422, 461, 689, 692 use "establish(es) that …". This is
not on the forbidden list but is a strong verb. It is consistent
with the manuscript's framing of permutation-validated findings as
"reproducible against the null"; not a blocker. Consider
"indicate" or "show" if a reviewer flags overclaim.

### MIN-5 — "AUROC = 0.669 / 0.753 / 0.758" vs four-decimal ground truth

The manuscript rounds all AUROC values to three decimals. The user
ground-truth values are SCI 0.6685, LLR 0.7531, SCI+LLR 0.7578,
which round to 0.669, 0.753, 0.758. **All three values are
internally consistent at the displayed precision.** No fix needed,
but if Nature Communications house style requires four decimals,
re-export from `sci_vs_llr_comparison.csv` would be needed.

### MIN-6 — `underperforms` (line 459)

§2.8 line 459 uses "underperforms", the morphological inverse of
the banned "outperform". This is not on the forbidden list and is
the honest direction; acceptable. Flagged for awareness in case
the user wants to apply the ban symmetrically.

### MIN-7 — Discussion section header `## 3. Discussion` vs internal labels `D1–D5`

The Discussion subsections are labelled `### D1. … ### D5.` rather
than `### 3.1 … ### 3.5`. This mirrors the source file but creates
two parallel numbering schemes inside the same section
(top-level `3.` and internal `D.`). §D5 is referenced from §2.8 as
"discussed further in §D5" (line 470), which works only because the
reader has been trained to recognise the `D` prefix. Acceptable for
a draft; reviewers usually do not push back, but house style varies.

### MIN-8 — `essentially` and `largely` appear repeatedly

Lines 87, 372, 1039 ("essentially absent / essentially
uncharacterised / essentially absent"); lines 455, 459, 467
("largely redundant / largely co-linear"). Used responsibly to
hedge the strength of "no signal" or "no contribution" claims;
not a flag.

### MIN-9 — All cited works in body text are present in §7 References inventory

Verified: every citation appearing in §1, §2, §3, §4 (Belrose,
Brandes, Deng, Detlefsen, Elnaggar, Firnberg, Frazer, Geva, Hayes,
Hsu, Jacquier, Kitzman, Lin, Marquet, Matreyek, Meier, Melamed,
Mighell, Notin 2022, Notin 2023, Phipson & Smyth, Rao 2019, Rao
2021, Rives, Simon & Zou, Stiffler, Sun, Tenney, Vig) is listed
in §7. No orphan citations and no orphan inventory entries.

---

## Forbidden-vocabulary scan (verbatim)

Body-text occurrences of the user's banned list:

| Word | Body-text count | Locations | Verdict |
|---|---:|---|---|
| prove | 0 | only line 198 (metacommentary listing) | clean |
| demonstrate | 0 | only line 198 (metacommentary listing) | clean |
| cause | 0 as a verb (only `causal verbs` line 27 and `causal claims` line 198, both meta) | — | clean (meta-only; consider stripping the meta-block before submission, see MIN-1) |
| **drive / driven** | **4** | lines 50, 153, 378, 974 | **violates ban — see MAJ-1** |
| mechanism | 0 | only line 199 (metacommentary listing) | clean |
| outperform | 0 | only line 199 (metacommentary listing) | clean |

---

## Claim-to-evidence cross-walk (C1–C7)

| Claim | Where stated | Statistical evidence | Source-of-truth file | Numbers consistent? |
|---|---|---|---|---|
| **C1** SCI carries mutation-effect signal | §2.1 (line 226–255); Fig 1 caption | T₁ = 2.79 × 10⁻³; z = 57; 5/5 proteins FDR-sig; bootstrap CI [2.67, 2.91] × 10⁻³ | `permutation_global_summary.csv`, `permutation_per_protein_summary.csv`, `bootstrap_summary.csv` | **Yes** — values match across §2.1 / Fig 1 / Sup S2-A / Sup S2-B |
| **C2** Top-ranked layer pairs enriched for EL | §2.2 (line 257–285); Fig 2(a–b) caption | T₂ = 17 (z = 8.01, q = 10⁻³); per-protein EL = 40/66/44/40 % (CBS/GAL4/PTEN/TEM1) | `permutation_global_summary.csv`, `*_pair_category_ratio_Composite.csv` | **Yes** |
| **C3** EL dominates 528-pair landscape | §2.3 (line 287–314); Fig 2(c) caption | T₃ = 3.99 × 10⁻² (z = 62.03); T₃ˢᵗᵈ = 1.02 (z = 8.05); 5/5 per-protein FDR-sig | `permutation_global_summary.csv`, `permutation_per_protein_summary.csv` | **Yes** |
| **C4** PABP forms a distinct cluster | §2.4 (line 316–344); Fig 3 caption; Sup Notes S1, S2 | Heatmap, dendrogram (≈0.50 vs ≈0.22 fusion), PCA (PC1 ≈ +0.38 vs ∈ [−0.18, +0.01]); singleton inequality 0.389 vs 0.397 | `pairwise_euclidean_distance.csv`, `pca_loadings.csv` | **Partially** — heatmap, dendrogram, PCA all correct. **Singleton inequality is arithmetically false** (see CRIT-1) |
| **C5** Recurrent layer pairs concentrate at layers 30–32 | §2.5 (line 346–373); Fig 4 caption | Layer 32 = 8/40 (20 %); layers 30–32 = 18/40; layers 18–27 = 0; Top50 union: Late 47.6 %, Middle 36.0 %, Early 16.4 % | `layer_recurrence_band_summary.csv`, `layer_recurrence_frequency.csv`, `layer_frequency.csv` | **Yes** |
| **C6** Per-protein signal is reproducible | §2.6 (line 375–391); Fig 5(a–b) | 5/5 proteins FDR-sig on T₁, T₃, T₃ˢᵗᵈ; per-protein z reported verbatim | `permutation_per_protein_summary.csv` | **Yes** — per-protein z values cross-check with Sup S2-B |
| **C7** Statistical robustness across corrections + bootstrap | §2.7 (line 393–418); Fig 5(c) | Phipson–Smyth p = 10⁻³ floor; Bonferroni 15/15 (T₁/T₃/T₃ˢᵗᵈ × 5); BH q < 0.05 matches; bootstrap CIs exclude null | `permutation_global_summary.csv`, `permutation_per_protein_summary.csv`, `bootstrap_summary.csv` | **Yes** |
| **§2.8** SCI vs ESM-2 LLR baseline | §2.8 (line 420–470); Table 1; §D5 | Pooled SCI AUROC 0.669, LLR AUROC 0.753, SCI+LLR AUROC 0.758; per-protein values verbatim from Table 1 | `sci_vs_llr_comparison.csv`, `baseline_summary.csv` | **Mostly** — Pooled values match user ground truth (0.6685 / 0.7531 / 0.7578) to three decimals; **two per-protein Δ values off by 0.001** (see MAJ-3); **SCI–LLR Pearson/Spearman missing** (see MAJ-2) |

---

## Figure–Results–Discussion cross-reference consistency

| Figure | Results section that introduces the numbers | Discussion section that interprets | Cross-reference matches? |
|---|---|---|---|
| Fig 1 (C1) | §2.1 | D1 (via F1 → C1) | **Yes** |
| Fig 2 (C2 + C3) | §2.2 + §2.3 | D2 (via F2 + F3) | **Label collision** — F2 ≠ Fig 2 collides only with reader expectation; Discussion's F2 / F3 actually point at panels of Fig 2 (correct). The `Fn = Cn` convention is internally consistent but **F3 in D1/D2 reads as "Fig 3" to a casual reader** — see CRIT-2 |
| Fig 3 (C4) | §2.4 | D3 (via F4 → C4) | **Collision** — Discussion D3 cites "F4" for PABP; reader expects Fig 4 (recurrence hub), not Fig 3 (PABP cluster) |
| Fig 4 (C5) | §2.5 | D1 + D4 (via F5 → C5) | **Collision** — Discussion cites "F5" for the layer-30–32 hub; reader expects Fig 5 (robustness), not Fig 4 (recurrence) |
| Fig 5 (C6 + C7) | §2.6 + §2.7 | Not explicitly cited in D1–D5 by `F` label | **Coverage gap** — Discussion does not discuss C6 / C7 explicitly. Justifiable (they are meta-claims), but should be noted |
| Table 1 (§2.8) | §2.8 | D5 | **Yes** — Table 1 cited from both §2.8 ("Table 1") and D5 ("Table 1") |

---

## Ready-to-submit Checklist

| # | Item | Status |
|---|---|---|
| 1 | Title present | ✅ |
| 2 | Abstract present, ≈ 200 words, single paragraph | ✅ |
| 3 | Introduction — all six paragraphs present, citations resolve | ✅ |
| 4 | Results §2.1–§2.8 — every claim sourced to a CSV | ✅ |
| 5 | Discussion D1–D5 — every paragraph sourced to a finding label | ⚠️  see CRIT-2 (label collision) |
| 6 | Methods §3.1–§3.12 — every analysis implementable from the cited script | ✅ |
| 7 | Figure Legends — Fig 1–5 + Table 1 | ✅ |
| 8 | Supplementary Information — S1, S2 (4 sub-tables), S Fig S1–S5, S Note S1–S2 | ⚠️  Sup Table S2-B contains 4 `(see CSV)` placeholders — MIN-2 |
| 9 | References — every body-text citation listed in §7 | ✅ |
| 10 | All headline numbers reproducible from on-disk CSVs | ⚠️  two per-protein ΔAUROC off by 0.001 — MAJ-3 |
| 11 | Forbidden vocabulary (prove / demonstrate / cause / mechanism / outperform) absent from body | ✅ |
| 12 | Forbidden vocabulary (drive / driven) absent from body | ❌ — `driven` in 4 places (Abstract, Intro, Results §2.6, Fig 3 caption) — MAJ-1 |
| 13 | Singleton condition for PABP (Sup Note S2) is arithmetically correct | ❌ — 0.389 < 0.397, inequality fails — CRIT-1 |
| 14 | Discussion `F1–F5` labels resolve unambiguously to Fig 1–Fig 5 | ❌ — collision, F3/F4/F5 point at Fig 2/Fig 3/Fig 4 — CRIT-2 |
| 15 | Top-level section numbering (1–7) internally consistent with subsection numbers (§2.X, §3.X, §D.X) | ⚠️  Methods §4 contains §3.1–§3.12 subsections — MAJ-4 |
| 16 | SCI–LLR Pearson / Spearman correlation (0.456 / 0.467) covered | ❌ — present in `sci_llr_correlation.csv` and companion paper docs but absent from final draft — MAJ-2 |
| 17 | Collaborator-facing meta-blocks stripped before submission (Intro "Notes on positioning", §2 / §3 / §5 / §6 leading quote-blocks) | ❌ — still present; consider stripping for the camera-ready PDF — MIN-1 |

**Blocking issues for submission:** CRIT-1, CRIT-2, MAJ-1, MAJ-2,
MAJ-3.

**Recommended order of fixes:** CRIT-1 (numeric error, fastest to
fix); CRIT-2 (decide on F vs Fig convention, then global
replace); MAJ-1 (4 × `driven` → 4 × `characterised-by` / `shifted-toward`);
MAJ-3 (re-export ΔAUROC from CSV, replace two values in prose);
MAJ-2 (add one Pearson / Spearman sentence to §2.8 + D5 or
explicitly link to supplementary). Minor issues can be addressed
in revision.
