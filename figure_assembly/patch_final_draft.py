#!/usr/bin/env python3
"""Apply remaining fixes to final_draft.md"""
import sys
sys.stdout.reconfigure(encoding='utf-8')

path = "d:/文件/工作室/website/paper/final_draft.md"

with open(path, 'r', encoding='utf-8') as f:
    content = f.read()

changes = 0

# Fix §2.6 Bonferroni claim
old = 'After Bonferroni\ncorrection over the 20 per-protein primary tests (5 proteins \xd7 4\nstatistics), each of these 15 (T₁, T₃, T₃ˢᵗᵈ across 5 proteins) tests\nretains significance at the conventional α = 0.05 level. The signal is\ntherefore not pooled out of one over-represented protein; it is\nreproduced within each.'
new = 'After Bonferroni\ncorrection over the 20 per-protein primary tests (5 proteins \xd7 4\nstatistics), 13 of these 15 tests (T₁, T₃, T₃ˢᵗᵈ across 5 proteins) retained\nsignificance at the conventional α = 0.05 level; CBS-T₁ (p_bonf = 0.60) and\nPTEN-T₁ (p_bonf = 0.14) did not survive the conservative Bonferroni threshold\ndespite being significant under BH-FDR control (q = 0.037 and q = 0.009,\nrespectively). The signal is\ntherefore not pooled out of one over-represented protein; it is\nreproduced within each, with strongest consistency for T₃ and T₃ˢᵗᵈ across\nall five proteins.'

if old in content:
    content = content.replace(old, new, 1)
    changes += 1
    print("[OK] final_draft.md §2.6")
else:
    print("[FAIL] final_draft.md §2.6")

# Fix §2.7 Bonferroni/FDR claim
old2 = 'Bonferroni significance at α = 0.05 was retained for T₁, T₃, T₃ˢᵗᵈ in 5/5\nproteins. (iii) Benjamini–Hochberg FDR with step-up monotonicity enforcement\nwas applied to the same 20-test family, with the\n`p_fdr` column containing the BH q-values; q < 0.05 calls match the\nBonferroni calls exactly on T₁, T₃, T₃ˢᵗᵈ across the five proteins.'
new2 = 'Bonferroni significance at α = 0.05 was retained for T₃ and T₃ˢᵗᵈ in all\nfive proteins, and for T₁ in three of five proteins (GAL4, PABP, TEM1). CBS-T₁\n(p_bonf = 0.60) and PTEN-T₁ (p_bonf = 0.14) did not survive Bonferroni.\n(iii) Benjamini–Hochberg FDR with step-up monotonicity enforcement was\napplied to the same 20-test family, with the\n`p_fdr` column containing the BH q-values;\nBH-FDR yielded q < 0.05 for T₁ in all five proteins, providing slightly broader\ndiscovery than Bonferroni while still controlling the expected false-discovery\nrate among the 20-test family.'

if old2 in content:
    content = content.replace(old2, new2, 1)
    changes += 1
    print("[OK] final_draft.md §2.7")
else:
    print("[FAIL] final_draft.md §2.7")

with open(path, 'w', encoding='utf-8') as f:
    f.write(content)
print(f"Total changes: {changes}")
