"""
Smoke test for the LLR pipeline that does NOT require fair-esm.

Checks:
  1. config + paths resolve, site metadata + WT sequences load.
  2. The mutant-string parsing (wt_aa / mut_aa / position) is sound.
  3. p7_baseline_summary.py and p7_sci_vs_llr_comparison.py run end-to-end on
     a fake esm2_llr_scores.csv whose 'LLR' column is set to (DMS_score + noise)
     so we should see Spearman ~ +0.7 -- proving the wiring is correct.

This does NOT validate the actual ESM-2 forward pass; it validates everything
*around* it. Run p7_compute_esm2_llr.py on the GPU box for the real LLR.
"""

import os
import sys
import shutil
import subprocess
import numpy as np
import pandas as pd

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from config import OUTPUT_DIR, SITE_METADATA, WT_SEQUENCES

STD_AAS = "ACDEFGHIKLMNPQRSTVWY"
LLR_CSV = os.path.join(OUTPUT_DIR, "esm2_llr_scores.csv")
BACKUP = LLR_CSV + ".backup_smoke"


def main():
    print("=" * 72)
    print("Smoke test: p7 pipeline (no fair-esm required)")
    print("=" * 72)

    site = pd.read_csv(SITE_METADATA)
    wt = pd.read_csv(WT_SEQUENCES)
    wt_map = dict(zip(wt["protein"], wt["wt_sequence"]))
    print(f"  site rows : {len(site)}")
    print(f"  proteins  : {list(wt_map.keys())}")

    # Parse mutant
    site["wt_aa"] = site["mutant"].str[0]
    site["mut_aa"] = site["mutant"].str[-1]
    pos_from_str = site["mutant"].map(lambda m: int(m[1:-1]))
    assert (pos_from_str == site["mutation_position"]).all(), "position mismatch"
    bad_aas = (set(site["wt_aa"]) | set(site["mut_aa"])) - set(STD_AAS)
    assert not bad_aas, f"non-std aas: {bad_aas}"
    print("  mutant parsing : OK (all std AAs, positions consistent)")

    # WT-AA sanity vs WT sequence
    n_check = 0; n_bad = 0
    for _, r in site.head(2000).iterrows():
        if wt_map[r["protein"]][r["mutation_position"] - 1] != r["wt_aa"]:
            n_bad += 1
        n_check += 1
    assert n_bad == 0, f"WT-AA mismatch in {n_bad}/{n_check} sampled rows"
    print(f"  WT-AA sanity (first 2000 rows) : OK")

    # Length check vs ESM max
    for p, s in wt_map.items():
        if len(s) > 1022:
            print(f"  WARN: {p} length {len(s)} > 1022 -- masked-marginal pipeline cannot handle without truncation")
        else:
            print(f"  {p:<5} len={len(s):>4}  OK (<=1022)")

    # Fabricate a fake esm2_llr_scores.csv where LLR ~ DMS_score (correlated)
    # so the downstream scripts have something to run on.
    if os.path.exists(LLR_CSV):
        shutil.copy(LLR_CSV, BACKUP)
        print(f"  (saved existing {LLR_CSV} -> {BACKUP})")

    rng = np.random.default_rng(42)
    fake = pd.DataFrame({
        "mutation_id": np.arange(len(site)),
        "protein": site["protein"].values,
        "mutant": site["mutant"].values,
        "parent_mutant": site["parent_mutant"].values,
        "mutation_index": site["mutation_index"].values,
        "n_mutations": site["n_mutations"].values,
        "mutation_position": site["mutation_position"].values,
        "wt_aa": site["wt_aa"].values,
        "mut_aa": site["mut_aa"].values,
        # WT-marginal: weaker signal
        "LLR_wt_marginal": site["DMS_score"].values + rng.normal(0, 1.0, len(site)),
        # masked-marginal (headline LLR): stronger signal
        "LLR": site["DMS_score"].values + rng.normal(0, 0.5, len(site)),
        "DMS_score": site["DMS_score"].values,
        "DMS_score_bin": site["DMS_score_bin"].values,
        "dataset": site["dataset"].values,
    })
    fake.to_csv(LLR_CSV, index=False, encoding="utf-8-sig")
    print(f"  wrote fake LLR csv ({len(fake)} rows) -> {LLR_CSV}")

    # Run the summary + comparison scripts
    here = os.path.dirname(os.path.abspath(__file__))
    for script in ["p7_baseline_summary.py", "p7_sci_vs_llr_comparison.py"]:
        print(f"\n  >>> running {script}")
        rc = subprocess.call([sys.executable, os.path.join(here, script)])
        assert rc == 0, f"{script} exited with {rc}"

    # Restore real LLR file if there was one
    if os.path.exists(BACKUP):
        shutil.move(BACKUP, LLR_CSV)
        print(f"\n  restored {LLR_CSV}")
    else:
        os.remove(LLR_CSV)
        print(f"\n  removed fake {LLR_CSV} (no original existed)")

    print("\nSMOKE TEST PASSED.")
    print("Now run p7_compute_esm2_llr.py on a host with fair-esm installed.")


if __name__ == "__main__":
    main()
