"""
Extract PCA loadings from the layer-pair category proportion matrix
and rank categories by |loading| on PC1 and PC2.

This script is intentionally self-contained: it re-loads the exact same
input matrix used by protein_clustering_layer_pair.py and re-fits PCA
with identical settings, so the loadings here correspond one-to-one to
the scatter plot already produced.

Loadings are defined as:
    L = components_ * sqrt(explained_variance_)
i.e. covariance between original (centered) variables and each PC.
This is the standard statistical/chemometrics convention and quantifies
each category's actual contribution to that PC's variance.

Outputs (written to data/p0_output/protein_clustering/):
    - pca_loadings.csv                 (Category, PC1_loading, PC2_loading)
    - pca_loadings_ranked.csv          (sorted by |PC1| and |PC2|, for the paper text)
"""

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

DATA_DIR = Path("D:/文件/工作室/website/data/processed/p0_output")
OUT_DIR = Path("D:/文件/工作室/website/data/p0_output/protein_clustering")
OUT_DIR.mkdir(parents=True, exist_ok=True)

PROTEINS = ["PABP", "CBS", "GAL4", "PTEN", "TEM1"]
CATEGORIES = ["EE", "EM", "EL", "MM", "ML", "LL"]


def load_matrix() -> pd.DataFrame:
    rows = []
    for prot in PROTEINS:
        csv_path = DATA_DIR / prot / f"{prot}_pair_category_ratio_Composite.csv"
        df = pd.read_csv(csv_path)
        row = {"Protein": prot}
        for _, r in df.iterrows():
            row[r["Pair_Type"]] = float(r["Ratio"])
        rows.append(row)
    mat = pd.DataFrame(rows).set_index("Protein")
    for c in CATEGORIES:
        if c not in mat.columns:
            mat[c] = 0.0
    return mat[CATEGORIES]


def main() -> None:
    mat = load_matrix()

    n_comp = min(mat.shape[0], mat.shape[1])
    pca = PCA(n_components=n_comp)
    pca.fit(mat.values)

    # Standard "loadings": component direction scaled by sqrt(eigenvalue)
    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
    pc1, pc2 = loadings[:, 0], loadings[:, 1]

    df_load = pd.DataFrame({
        "Category": CATEGORIES,
        "PC1_loading": pc1,
        "PC2_loading": pc2,
    })
    df_load.to_csv(OUT_DIR / "pca_loadings.csv", index=False, float_format="%.6f")

    # Ranked by absolute magnitude, for the paper text
    df_rank_pc1 = df_load.assign(abs_PC1=df_load["PC1_loading"].abs()) \
                         .sort_values("abs_PC1", ascending=False) \
                         .reset_index(drop=True)
    df_rank_pc2 = df_load.assign(abs_PC2=df_load["PC2_loading"].abs()) \
                         .sort_values("abs_PC2", ascending=False) \
                         .reset_index(drop=True)

    ranked = pd.DataFrame({
        "Rank": np.arange(1, len(CATEGORIES) + 1),
        "PC1_Category": df_rank_pc1["Category"].values,
        "PC1_loading": df_rank_pc1["PC1_loading"].values,
        "PC1_abs": df_rank_pc1["abs_PC1"].values,
        "PC2_Category": df_rank_pc2["Category"].values,
        "PC2_loading": df_rank_pc2["PC2_loading"].values,
        "PC2_abs": df_rank_pc2["abs_PC2"].values,
    })
    ranked.to_csv(OUT_DIR / "pca_loadings_ranked.csv", index=False, float_format="%.6f")

    print("Explained variance ratio:", np.round(pca.explained_variance_ratio_, 4))
    print()
    print("pca_loadings.csv:")
    print(df_load.to_string(index=False, float_format=lambda v: f"{v: .4f}"))
    print()
    print("Ranking by |PC1|:")
    print(df_rank_pc1[["Category", "PC1_loading", "abs_PC1"]]
          .to_string(index=False, float_format=lambda v: f"{v: .4f}"))
    print()
    print("Ranking by |PC2|:")
    print(df_rank_pc2[["Category", "PC2_loading", "abs_PC2"]]
          .to_string(index=False, float_format=lambda v: f"{v: .4f}"))
    print()
    print(f"Wrote: {OUT_DIR / 'pca_loadings.csv'}")
    print(f"Wrote: {OUT_DIR / 'pca_loadings_ranked.csv'}")


if __name__ == "__main__":
    main()
