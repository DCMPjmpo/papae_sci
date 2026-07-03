"""
P7 — Baseline summary for the ESM-2 LLR scorer.

Reads: data/processed/esm2_llr_scores.csv
Writes: data/processed/baseline_summary.csv

For each scope (Pooled + each of the 5 proteins) computes:
    - n
    - n_pos / n_neg  (DMS_score_bin counts)
    - Spearman(score, DMS_score)        + p-value
    - Spearman(score, DMS_score_z)      + p-value   (z-normalized within protein)
    - AUROC raw      (score, bin)
    - AUROC corrected (sign-aligned to maximize AUROC; matches SCI convention)

Done independently for:
    - LLR_wt_marginal
    - LLR (masked-marginal, canonical Meier 2021)

This mirrors the metrics defined in P1.5_validate_sci_signal.py so SCI and LLR
can be compared on the same axes.
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from config import OUTPUT_DIR

LLR_CSV = os.path.join(OUTPUT_DIR, "esm2_llr_scores.csv")
SUMMARY_CSV = os.path.join(OUTPUT_DIR, "baseline_summary.csv")

SCORERS = ["LLR_wt_marginal", "LLR"]


def _zscore(g):
    s = g.std()
    return (g - g.mean()) / s if s and not np.isnan(s) else g * 0.0


def summarize(df, scope_label, score_col):
    sub = df[[score_col, "DMS_score", "DMS_score_z", "DMS_score_bin"]].dropna()
    n = len(sub)
    if n < 5:
        return {
            "scope": scope_label, "scorer": score_col, "n": n,
            "n_pos": int((sub["DMS_score_bin"] == 1).sum()),
            "n_neg": int((sub["DMS_score_bin"] == 0).sum()),
            "spearman_r_raw": np.nan, "spearman_p_raw": np.nan,
            "spearman_r_z":   np.nan, "spearman_p_z":   np.nan,
            "AUROC_raw":      np.nan, "AUROC_corrected": np.nan,
            "score_sign":     np.nan,
        }
    x = sub[score_col].values
    y_raw = sub["DMS_score"].values
    y_z = sub["DMS_score_z"].values
    y_bin = sub["DMS_score_bin"].values

    sr_raw, sp_raw = spearmanr(x, y_raw)
    sr_z, sp_z = spearmanr(x, y_z)

    n_pos = int((y_bin == 1).sum())
    n_neg = int((y_bin == 0).sum())

    if n_pos > 0 and n_neg > 0:
        auc_raw = roc_auc_score(y_bin, x)
        auc_rev = roc_auc_score(y_bin, -x)
        if auc_rev > auc_raw:
            auc_corrected, sign = auc_rev, -1
        else:
            auc_corrected, sign = auc_raw, +1
    else:
        auc_raw = auc_corrected = sign = np.nan

    return {
        "scope": scope_label, "scorer": score_col, "n": n,
        "n_pos": n_pos, "n_neg": n_neg,
        "spearman_r_raw": round(float(sr_raw), 4), "spearman_p_raw": float(sp_raw),
        "spearman_r_z":   round(float(sr_z), 4),   "spearman_p_z":   float(sp_z),
        "AUROC_raw":      round(float(auc_raw), 4),
        "AUROC_corrected": round(float(auc_corrected), 4),
        "score_sign": sign,
    }


def main():
    if not os.path.exists(LLR_CSV):
        sys.exit(f"ERROR: {LLR_CSV} not found. Run p7_compute_esm2_llr.py first.")

    df = pd.read_csv(LLR_CSV)
    print(f"Loaded {len(df)} rows from {LLR_CSV}")
    print(df.groupby("protein").size())

    # z-normalize DMS within protein (mirrors P1.5_validate_sci_signal.py)
    df["DMS_score_z"] = df.groupby("protein")["DMS_score"].transform(_zscore)

    rows = []
    for scorer in SCORERS:
        rows.append(summarize(df, "Pooled", scorer))
        for prot in sorted(df["protein"].unique()):
            rows.append(summarize(df[df["protein"] == prot], prot, scorer))

    summary = pd.DataFrame(rows)
    summary = summary[[
        "scope", "scorer", "n", "n_pos", "n_neg",
        "spearman_r_raw", "spearman_p_raw",
        "spearman_r_z", "spearman_p_z",
        "AUROC_raw", "AUROC_corrected", "score_sign",
    ]]
    summary.to_csv(SUMMARY_CSV, index=False, encoding="utf-8-sig")
    print(f"\n[OK] Wrote {SUMMARY_CSV}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
