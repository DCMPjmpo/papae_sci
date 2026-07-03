"""
P7 — SCI vs LLR vs SCI+LLR three-way comparison.

Reads:
    data/processed/esm2_llr_scores.csv                  (this branch)
    data/processed/all_proteins_sci_site_scores_mean.npy (existing SCI)
    data/processed/all_proteins_site_metadata.csv        (alignment & labels)

Writes:
    data/processed/sci_vs_llr_comparison.csv

For every scope (Pooled + each protein) and every scorer (SCI, LLR, SCI+LLR):
    n
    Spearman |r| vs DMS_score      ( + sign of natural correlation )
    AUROC_corrected (sign-aligned to maximize, matches SCI script convention)
    AUPRC_corrected

The SCI+LLR combo is a 5-fold CV (stratified) logistic regression on the
standardized [SCI, LLR] feature pair, fit fresh per fold (no leakage). For the
Spearman row of the combo, we use the held-out CV probability as the score.
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_auc_score, average_precision_score

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from config import OUTPUT_DIR, SITE_METADATA, SCI_SCORES_MEAN

LLR_CSV = os.path.join(OUTPUT_DIR, "esm2_llr_scores.csv")
OUT_CSV = os.path.join(OUTPUT_DIR, "sci_vs_llr_comparison.csv")

SEED = 42
N_SPLITS = 5


def signed_metrics(score, dms, bin_):
    """Spearman + AUROC + AUPRC, sign-corrected (matches SCI script convention)."""
    mask = ~np.isnan(score)
    score = score[mask]; dms = dms[mask]; bin_ = bin_[mask]
    n = len(score)
    if n < 5 or len(np.unique(bin_)) < 2:
        return {"n": n,
                "spearman_r": np.nan, "spearman_p": np.nan, "spearman_sign": np.nan,
                "AUROC_raw": np.nan, "AUROC_corrected": np.nan,
                "AUPRC_raw": np.nan, "AUPRC_corrected": np.nan}
    sr, sp = spearmanr(score, dms)
    auc_fwd = roc_auc_score(bin_, score)
    auc_rev = roc_auc_score(bin_, -score)
    auprc_fwd = average_precision_score(bin_, score)
    auprc_rev = average_precision_score(bin_, -score)
    return {
        "n": n,
        "spearman_r": round(abs(float(sr)), 4),
        "spearman_p": float(sp),
        "spearman_sign": int(np.sign(sr)) if not np.isnan(sr) else np.nan,
        "AUROC_raw": round(float(auc_fwd), 4),
        "AUROC_corrected": round(float(max(auc_fwd, auc_rev)), 4),
        "AUPRC_raw": round(float(auprc_fwd), 4),
        "AUPRC_corrected": round(float(max(auprc_fwd, auprc_rev)), 4),
    }


def cv_combo_score(X, y_bin):
    """5-fold stratified CV logistic regression; returns out-of-fold P(class=1)."""
    n = len(X)
    if n < N_SPLITS * 2 or len(np.unique(y_bin)) < 2:
        return np.full(n, np.nan)
    pred = np.full(n, np.nan)
    skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=SEED)
    pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(max_iter=500, solver="lbfgs")),
    ])
    for tr, te in skf.split(X, y_bin):
        if len(np.unique(y_bin[tr])) < 2:
            continue
        pipe.fit(X[tr], y_bin[tr])
        pred[te] = pipe.predict_proba(X[te])[:, 1]
    return pred


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
        sys.exit(f"ERROR: length mismatch  LLR={len(llr_df)}  SCI={len(sci)}  site={len(site)}")

    # Sanity check alignment: the LLR CSV was built from SITE_METADATA in row order
    if not (llr_df["protein"].values == site["protein"].values).all() or \
       not (llr_df["mutant"].values == site["mutant"].values).all():
        sys.exit("ERROR: LLR CSV does not align with site metadata. Re-run p7_compute_esm2_llr.py.")

    df = llr_df.copy()
    df["SCI"] = sci

    print(f"Loaded {len(df)} rows; proteins: {df['protein'].unique().tolist()}")
    print(df.groupby("protein").size())

    rows = []
    scopes = ["Pooled"] + sorted(df["protein"].unique().tolist())

    for scope in scopes:
        sub = df if scope == "Pooled" else df[df["protein"] == scope]
        sub = sub.dropna(subset=["SCI", "LLR"]).reset_index(drop=True)
        if len(sub) < 20:
            continue
        sci_v = sub["SCI"].values
        llr_v = sub["LLR"].values
        dms = sub["DMS_score"].values
        binv = sub["DMS_score_bin"].values

        # SCI alone
        m_sci = signed_metrics(sci_v, dms, binv)
        m_sci.update({"scope": scope, "scorer": "SCI"})
        rows.append(m_sci)

        # LLR alone (masked-marginal, canonical)
        m_llr = signed_metrics(llr_v, dms, binv)
        m_llr.update({"scope": scope, "scorer": "LLR"})
        rows.append(m_llr)

        # SCI + LLR via 5-fold CV logistic regression
        X = np.stack([sci_v, llr_v], axis=1)
        combo_prob = cv_combo_score(X, binv)
        m_combo = signed_metrics(combo_prob, dms, binv)
        m_combo.update({"scope": scope, "scorer": "SCI+LLR (CV-LR)"})
        rows.append(m_combo)

    out = pd.DataFrame(rows)[[
        "scope", "scorer", "n",
        "spearman_r", "spearman_p", "spearman_sign",
        "AUROC_raw", "AUROC_corrected",
        "AUPRC_raw", "AUPRC_corrected",
    ]]
    out.to_csv(OUT_CSV, index=False, encoding="utf-8-sig")
    print(f"\n[OK] Wrote {OUT_CSV}")
    print(out.to_string(index=False))

    # ---------- Delta summary: combo vs the best of the two singles ----------
    pivot = (out
             .pivot_table(index="scope", columns="scorer",
                          values=["AUROC_corrected", "spearman_r"])
             .round(4))
    print("\nPivot (AUROC_corrected | spearman_r):")
    print(pivot.to_string())


if __name__ == "__main__":
    main()
