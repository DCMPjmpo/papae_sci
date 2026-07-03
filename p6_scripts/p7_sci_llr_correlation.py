"""
P7 — SCI vs LLR correlation (Pearson + Spearman).

Reads:
    data/processed/esm2_llr_scores.csv
    data/processed/all_proteins_sci_site_scores_mean.npy
    data/processed/all_proteins_site_metadata.csv

Computes correlations between SCI and LLR (masked-marginal, headline LLR column):
    Pooled   : Pearson(SCI, LLR), Spearman(SCI, LLR)
    Per-protein : same two coefficients for each protein

Writes:
    data/processed/sci_llr_correlation.csv
    columns: scope, n, pearson_r, pearson_p, spearman_r, spearman_p
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from config import OUTPUT_DIR, SITE_METADATA, SCI_SCORES_MEAN

LLR_CSV = os.path.join(OUTPUT_DIR, "esm2_llr_scores.csv")
OUT_CSV = os.path.join(OUTPUT_DIR, "sci_llr_correlation.csv")

MIN_N = 5  # minimum sample size to compute a correlation


def correlation_row(scope, sci, llr):
    """Return a dict with Pearson + Spearman for one (sci, llr) pair."""
    mask = ~(np.isnan(sci) | np.isnan(llr))
    sci_v = sci[mask]
    llr_v = llr[mask]
    n = int(len(sci_v))
    if n < MIN_N:
        return {
            "scope": scope, "n": n,
            "pearson_r": np.nan, "pearson_p": np.nan,
            "spearman_r": np.nan, "spearman_p": np.nan,
        }
    pr, pp = pearsonr(sci_v, llr_v)
    sr, sp = spearmanr(sci_v, llr_v)
    return {
        "scope": scope,
        "n": n,
        "pearson_r": float(pr),
        "pearson_p": float(pp),
        "spearman_r": float(sr),
        "spearman_p": float(sp),
    }


def main():
    if not os.path.exists(LLR_CSV):
        sys.exit(f"ERROR: {LLR_CSV} not found. Run p7_compute_esm2_llr.py first.")
    if not os.path.exists(SCI_SCORES_MEAN):
        sys.exit(f"ERROR: {SCI_SCORES_MEAN} not found.")
    if not os.path.exists(SITE_METADATA):
        sys.exit(f"ERROR: {SITE_METADATA} not found.")

    llr_df = pd.read_csv(LLR_CSV).sort_values("mutation_id").reset_index(drop=True)
    sci = np.load(SCI_SCORES_MEAN)
    site = pd.read_csv(SITE_METADATA)

    if not (len(llr_df) == len(sci) == len(site)):
        sys.exit(
            f"ERROR: length mismatch  LLR={len(llr_df)}  SCI={len(sci)}  site={len(site)}"
        )

    # Sanity-check alignment (LLR CSV was built from SITE_METADATA in row order)
    if not (llr_df["protein"].values == site["protein"].values).all() or \
       not (llr_df["mutant"].values == site["mutant"].values).all():
        sys.exit("ERROR: LLR CSV does not align with site metadata. "
                 "Re-run p7_compute_esm2_llr.py.")

    df = pd.DataFrame({
        "protein": llr_df["protein"].values,
        "SCI": sci,
        "LLR": llr_df["LLR"].values,
    })

    print(f"Loaded {len(df)} rows; proteins: {sorted(df['protein'].unique().tolist())}")
    print(df.groupby("protein").size())

    rows = [correlation_row("Pooled", df["SCI"].values, df["LLR"].values)]
    for prot in sorted(df["protein"].unique()):
        sub = df[df["protein"] == prot]
        rows.append(correlation_row(prot, sub["SCI"].values, sub["LLR"].values))

    out = pd.DataFrame(rows)[[
        "scope", "n",
        "pearson_r", "pearson_p",
        "spearman_r", "spearman_p",
    ]]
    out.to_csv(OUT_CSV, index=False, encoding="utf-8-sig")
    print(f"\n[OK] Wrote {OUT_CSV}")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
