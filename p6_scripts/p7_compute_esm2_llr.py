"""
P7 — ESM-2 Zero-shot LLR Baseline.

Compute Log-Likelihood Ratio (LLR) = logP(mut) - logP(wt) for every expanded
mutation row, using two standard scorers from Meier et al. 2021:
    * WT-marginal     : softmax over the un-masked WT-context logits
    * Masked-marginal : softmax over the logits at the singly-masked position
                        (the canonical ProteinGym baseline)

DISCLOSURE
----------
The existing forward pass (p0_p5/extract_esm2.py) saved 33-layer hidden states
with return_contacts=False; LM-head logits were NOT cached. LLR strictly
requires logits, so this script runs the minimum-additional forward work:
    5 forwards for WT-marginals  +  ~1,300 forwards for masked-marginals
which is ~70x fewer than the 95,142 forwards used for embedding extraction,
on the same model / weights / sequences.

OUTPUT
------
data/processed/esm2_llr_scores.csv  (95,142 rows)
columns: mutation_id, protein, mutant, parent_mutant, mutation_index,
         n_mutations, mutation_position, wt_aa, mut_aa,
         LLR_wt_marginal, LLR, DMS_score, DMS_score_bin, dataset
(LLR == LLR_masked_marginal, kept as the headline column for downstream code.)
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import torch
from tqdm import tqdm

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from config import (
    OUTPUT_DIR, MODEL_DIR,
    ESM_MODEL_NAME, ESM_MAX_SEQ_LEN,
    SITE_METADATA, WT_SEQUENCES,
)

# Force ESM model loader to use the project's cached weights
os.environ["TORCH_HOME"] = MODEL_DIR
torch.hub.set_dir(MODEL_DIR)

try:
    import esm  # fair-esm
except ImportError:
    sys.stderr.write(
        "\nERROR: fair-esm not installed in this Python environment.\n"
        "       The original pipeline (p0_p5/extract_esm2.py) uses fair-esm — install with:\n"
        "           pip install fair-esm\n"
        "       Then re-run this script (CPU is fine for ~1,300 forwards; GPU recommended).\n\n"
    )
    sys.exit(1)


LLR_OUTPUT = os.path.join(OUTPUT_DIR, "esm2_llr_scores.csv")

# 20 standard amino acids (the only AAs that appear in WT/mut substitutions in the table)
STD_AAS = "ACDEFGHIKLMNPQRSTVWY"


def main():
    print("=" * 72)
    print("P7 — ESM-2 Zero-shot LLR Baseline (WT- + Masked-marginal)")
    print("=" * 72)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device : {device}")
    print(f"Model  : {ESM_MODEL_NAME}")
    print(f"Cache  : {torch.hub.get_dir()}")

    # ---------- Load model ----------
    t0 = time.time()
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.to(device).eval()
    batch_converter = alphabet.get_batch_converter()
    mask_idx = alphabet.mask_idx
    aa_to_tokidx = {aa: alphabet.get_idx(aa) for aa in STD_AAS}
    print(f"Model loaded in {time.time()-t0:.1f}s; mask_idx={mask_idx}; "
          f"{len(aa_to_tokidx)} AA -> token mappings.")

    # ---------- Load data ----------
    site = pd.read_csv(SITE_METADATA)
    wt_df = pd.read_csv(WT_SEQUENCES)
    wt_map = dict(zip(wt_df["protein"], wt_df["wt_sequence"]))
    print(f"Site metadata : {len(site)} rows")
    print(f"WT sequences  : {len(wt_map)} proteins -> {list(wt_map.keys())}")

    missing = [p for p in site["protein"].unique() if p not in wt_map]
    if missing:
        raise KeyError(f"No WT sequence for proteins: {missing}")

    # Parse wt_aa, mut_aa from the 'mutant' string (always single-substitution in expanded table)
    site["wt_aa"] = site["mutant"].str[0]
    site["mut_aa"] = site["mutant"].str[-1]
    pos_match = site["mutant"].map(lambda m: int(m[1:-1])) == site["mutation_position"]
    if not pos_match.all():
        n_bad = int((~pos_match).sum())
        raise ValueError(f"mutation_position vs mutant-string mismatch in {n_bad} rows")

    bad_aas = set(site["wt_aa"]) | set(site["mut_aa"])
    bad_aas -= set(STD_AAS)
    if bad_aas:
        raise ValueError(f"Non-standard AAs in mutant table: {bad_aas}")

    # ---------- Pre-tokenize WT once per protein ----------
    wt_tokens = {}
    wt_len = {}
    for prot, seq in wt_map.items():
        if len(seq) > ESM_MAX_SEQ_LEN:
            raise ValueError(f"{prot} length {len(seq)} > ESM_MAX_SEQ_LEN={ESM_MAX_SEQ_LEN}")
        _, _, tokens = batch_converter([(prot, seq)])
        wt_tokens[prot] = tokens
        wt_len[prot] = len(seq)

    # ---------- Phase 1: WT-marginals (one forward per protein) ----------
    print("\n[Phase 1] WT-marginals: 1 forward per protein (5 total)")
    wt_logp = {}
    for prot, tokens in wt_tokens.items():
        with torch.no_grad():
            out = model(tokens.to(device), repr_layers=[], return_contacts=False)
        L = wt_len[prot]
        log_probs = torch.log_softmax(out["logits"][0, 1:1 + L, :], dim=-1).cpu().numpy()
        wt_logp[prot] = log_probs
        print(f"  {prot:<5}  L={L:>4}  logits {tuple(out['logits'].shape)} "
              f"-> log_probs {log_probs.shape}")

    # ---------- Phase 2: Masked-marginals (one forward per unique (protein, pos)) ----------
    print("\n[Phase 2] Masked-marginals: 1 forward per unique (protein, position)")
    unique_pp = (
        site[["protein", "mutation_position"]]
        .drop_duplicates()
        .sort_values(["protein", "mutation_position"])
        .reset_index(drop=True)
    )
    print(f"  Unique (protein, position) pairs : {len(unique_pp)}")
    print(f"  Per-protein unique positions     : "
          f"{unique_pp.groupby('protein').size().to_dict()}")

    masked_logp = {}  # (protein, pos_1based) -> np.array(20,) log-probs over STD_AAS
    for _, row in tqdm(unique_pp.iterrows(), total=len(unique_pp),
                       desc="masked-marginal forwards"):
        prot = row["protein"]
        pos = int(row["mutation_position"])
        L = wt_len[prot]
        if pos < 1 or pos > L:
            masked_logp[(prot, pos)] = None
            continue
        tok = wt_tokens[prot].clone()
        tok[0, pos] = mask_idx  # token[0]=<cls>, token[i]=residue i for i=1..L
        with torch.no_grad():
            out = model(tok.to(device), repr_layers=[], return_contacts=False)
        log_probs = torch.log_softmax(out["logits"][0, pos, :], dim=-1).cpu().numpy()
        masked_logp[(prot, pos)] = np.array(
            [log_probs[aa_to_tokidx[aa]] for aa in STD_AAS], dtype=np.float64
        )

    # ---------- Phase 3: Compute per-row LLR ----------
    print("\n[Phase 3] Per-row LLR for all 95,142 expanded rows")
    aa_idx_in_std = {aa: i for i, aa in enumerate(STD_AAS)}
    llr_wt = np.full(len(site), np.nan)
    llr_mm = np.full(len(site), np.nan)
    n_mismatch = 0
    n_oor = 0
    for i, row in tqdm(site.iterrows(), total=len(site), desc="rows"):
        prot = row["protein"]
        pos = int(row["mutation_position"])
        wt_aa = row["wt_aa"]
        mut_aa = row["mut_aa"]

        # WT sequence sanity: residue at pos-1 (0-based) must equal wt_aa
        if wt_map[prot][pos - 1] != wt_aa:
            n_mismatch += 1
            continue

        # WT-marginal LLR (always available)
        lp_wt = wt_logp[prot][pos - 1]
        llr_wt[i] = float(lp_wt[aa_to_tokidx[mut_aa]] - lp_wt[aa_to_tokidx[wt_aa]])

        # Masked-marginal LLR
        mm = masked_logp.get((prot, pos))
        if mm is None:
            n_oor += 1
            continue
        llr_mm[i] = float(mm[aa_idx_in_std[mut_aa]] - mm[aa_idx_in_std[wt_aa]])

    if n_mismatch:
        print(f"  WARNING: {n_mismatch} rows had WT-AA mismatch with WT sequence (set to NaN)")
    if n_oor:
        print(f"  WARNING: {n_oor} rows had out-of-range mutation_position (set to NaN)")

    # ---------- Output ----------
    out_df = pd.DataFrame({
        "mutation_id": np.arange(len(site)),
        "protein": site["protein"].values,
        "mutant": site["mutant"].values,
        "parent_mutant": site["parent_mutant"].values,
        "mutation_index": site["mutation_index"].values,
        "n_mutations": site["n_mutations"].values,
        "mutation_position": site["mutation_position"].values,
        "wt_aa": site["wt_aa"].values,
        "mut_aa": site["mut_aa"].values,
        "LLR_wt_marginal": llr_wt,
        "LLR": llr_mm,            # headline = masked-marginal (canonical Meier 2021)
        "DMS_score": site["DMS_score"].values,
        "DMS_score_bin": site["DMS_score_bin"].values,
        "dataset": site["dataset"].values,
    })
    out_df.to_csv(LLR_OUTPUT, index=False, encoding="utf-8-sig")
    print(f"\n[OK] Wrote {len(out_df)} rows -> {LLR_OUTPUT}")
    print(out_df.head().to_string(index=False))


if __name__ == "__main__":
    main()
