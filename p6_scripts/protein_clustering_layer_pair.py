"""
Protein-level clustering based on layer-pair category proportions.

Goal (per project narrative):
    Test whether PABP exhibits a distinct cross-layer encoding signature
    compared with CBS, GAL4, PTEN, and TEM1 — using only the proportions
    of layer-pair categories (EE, EM, EL, MM, ML, LL) within each
    protein's Composite Top50 layer pairs.

Outputs (publication-quality, 600 dpi):
    1. heatmap_layer_pair_categories.{png,pdf}
    2. dendrogram_protein_clustering.{png,pdf}
    3. pca_protein_layer_pair.{png,pdf}

Notes:
    - We retain all six categories present in the source CSVs (ML included
      for completeness) but the narrative focus is on the five informative
      Early/Middle/Late pairings.
    - All distance/PCA computations are on the proportion vectors as-is;
      no standardization is applied, so the geometry reflects the actual
      compositional differences across proteins.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
DATA_DIR = Path("D:/文件/工作室/website/data/processed/p0_output")
OUT_DIR = Path("D:/文件/工作室/website/data/p0_output/protein_clustering")
OUT_DIR.mkdir(parents=True, exist_ok=True)

PROTEINS = ["PABP", "CBS", "GAL4", "PTEN", "TEM1"]
CATEGORIES = ["EE", "EM", "EL", "MM", "ML", "LL"]
CATEGORY_LABELS = {
    "EE": "Early–Early",
    "EM": "Early–Middle",
    "EL": "Early–Late",
    "MM": "Middle–Middle",
    "ML": "Middle–Late",
    "LL": "Late–Late",
}

# Highlight PABP consistently across all panels
PABP_COLOR = "#D7263D"
OTHER_COLOR = "#1B4965"
PROTEIN_COLORS = {p: (PABP_COLOR if p == "PABP" else OTHER_COLOR) for p in PROTEINS}

# ---------------------------------------------------------------------------
# Global matplotlib style — publication-quality, neutral
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "axes.linewidth": 1.0,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9,
    "legend.frameon": False,
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def load_proportion_matrix() -> pd.DataFrame:
    """Load each protein's Composite category-ratio CSV into a matrix."""
    rows = []
    for prot in PROTEINS:
        csv_path = DATA_DIR / prot / f"{prot}_pair_category_ratio_Composite.csv"
        if not csv_path.exists():
            raise FileNotFoundError(f"Missing input: {csv_path}")
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


# ---------------------------------------------------------------------------
# Figure 1: Heatmap
# ---------------------------------------------------------------------------
def plot_heatmap(mat: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.2))

    data = mat.values
    im = ax.imshow(
        data,
        aspect="auto",
        cmap="RdBu_r",
        vmin=0.0,
        vmax=float(np.nanmax(data)) if np.nanmax(data) > 0 else 1.0,
    )

    ax.set_xticks(np.arange(len(CATEGORIES)))
    ax.set_xticklabels([CATEGORY_LABELS[c] for c in CATEGORIES], rotation=30, ha="right")
    ax.set_yticks(np.arange(len(mat.index)))
    ax.set_yticklabels(mat.index)

    # Highlight PABP row label
    for tick in ax.get_yticklabels():
        if tick.get_text() == "PABP":
            tick.set_color(PABP_COLOR)
            tick.set_fontweight("bold")

    # Annotate each cell with its proportion
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            v = data[i, j]
            txt_color = "white" if v > 0.35 else "black"
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    color=txt_color, fontsize=9)

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Proportion in Composite Top50 layer pairs", rotation=270, labelpad=14)
    cbar.outline.set_linewidth(0.8)

    ax.set_xlabel("Layer-pair category")
    ax.set_ylabel("Protein")
    ax.set_title("Layer-pair category composition across five proteins",
                 pad=10)

    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

    fig.savefig(OUT_DIR / "heatmap_layer_pair_categories.png")
    fig.savefig(OUT_DIR / "heatmap_layer_pair_categories.pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2: Hierarchical clustering dendrogram
# ---------------------------------------------------------------------------
def plot_dendrogram(mat: pd.DataFrame) -> np.ndarray:
    # Euclidean distance on raw proportion vectors; average linkage is
    # conservative and stable for very small n.
    dist = pdist(mat.values, metric="euclidean")
    Z = linkage(dist, method="average")

    fig, ax = plt.subplots(figsize=(6.4, 4.2))

    ddata = dendrogram(
        Z,
        labels=list(mat.index),
        leaf_font_size=11,
        color_threshold=0,                  # disable scipy's auto-coloring
        above_threshold_color="#444444",
        ax=ax,
    )

    # Color each leaf label by protein identity (PABP highlighted)
    for lbl in ax.get_xmajorticklabels():
        name = lbl.get_text()
        lbl.set_color(PROTEIN_COLORS.get(name, OTHER_COLOR))
        if name == "PABP":
            lbl.set_fontweight("bold")

    ax.set_ylabel("Euclidean distance (average linkage)")
    ax.set_title("Hierarchical clustering of proteins\nby layer-pair category composition",
                 pad=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_linewidth(0.8)

    legend_handles = [
        Patch(facecolor=PABP_COLOR, label="PABP"),
        Patch(facecolor=OTHER_COLOR, label="CBS / GAL4 / PTEN / TEM1"),
    ]
    ax.legend(handles=legend_handles, loc="upper right")

    fig.savefig(OUT_DIR / "dendrogram_protein_clustering.png")
    fig.savefig(OUT_DIR / "dendrogram_protein_clustering.pdf")
    plt.close(fig)

    return squareform(dist)


# ---------------------------------------------------------------------------
# Figure 3: PCA scatter
# ---------------------------------------------------------------------------
def plot_pca(mat: pd.DataFrame) -> None:
    # No standardization: keep the geometry compositional.
    n_comp = min(mat.shape[0], mat.shape[1])
    pca = PCA(n_components=n_comp)
    coords = pca.fit_transform(mat.values)
    ev = pca.explained_variance_ratio_

    fig, ax = plt.subplots(figsize=(5.8, 5.0))

    for i, prot in enumerate(mat.index):
        color = PROTEIN_COLORS[prot]
        ax.scatter(
            coords[i, 0], coords[i, 1],
            s=180 if prot == "PABP" else 130,
            c=color,
            edgecolor="black",
            linewidth=0.8,
            zorder=3,
            label=prot,
        )
        # Offset labels slightly so they don't overlap markers
        ax.annotate(
            prot,
            (coords[i, 0], coords[i, 1]),
            xytext=(7, 6),
            textcoords="offset points",
            fontsize=11,
            fontweight=("bold" if prot == "PABP" else "normal"),
            color=color,
        )

    ax.axhline(0, color="#BBBBBB", linewidth=0.6, zorder=1)
    ax.axvline(0, color="#BBBBBB", linewidth=0.6, zorder=1)

    ax.set_xlabel(f"PC1 ({ev[0] * 100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({ev[1] * 100:.1f}% variance)")
    ax.set_title("PCA of proteins by layer-pair category composition", pad=10)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_linewidth(0.8)

    # Add loading arrows for the six categories on PC1/PC2
    loadings = pca.components_[:2].T * np.sqrt(pca.explained_variance_[:2])
    scale = 0.9 * np.max(np.abs(coords[:, :2])) / (np.max(np.abs(loadings)) + 1e-9)
    for k, cat in enumerate(CATEGORIES):
        lx, ly = loadings[k] * scale
        if np.hypot(lx, ly) < 1e-3:
            continue
        ax.annotate(
            "",
            xy=(lx, ly), xytext=(0, 0),
            arrowprops=dict(arrowstyle="->", color="#888888", lw=0.8),
            zorder=2,
        )
        ax.text(lx * 1.08, ly * 1.08, CATEGORY_LABELS[cat],
                color="#555555", fontsize=8, ha="center", va="center")

    fig.savefig(OUT_DIR / "pca_protein_layer_pair.png")
    fig.savefig(OUT_DIR / "pca_protein_layer_pair.pdf")
    plt.close(fig)

    return ev


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    mat = load_proportion_matrix()
    mat.to_csv(OUT_DIR / "layer_pair_category_matrix.csv")

    print("Input matrix (rows = proteins, cols = categories):")
    print(mat.round(3).to_string())
    print()

    plot_heatmap(mat)
    dist_matrix = plot_dendrogram(mat)
    ev = plot_pca(mat)

    # Pairwise Euclidean distances for the figure caption / Results text
    dist_df = pd.DataFrame(dist_matrix, index=mat.index, columns=mat.index)
    dist_df.to_csv(OUT_DIR / "pairwise_euclidean_distance.csv")
    print("Pairwise Euclidean distances:")
    print(dist_df.round(3).to_string())
    print()

    print("PCA explained variance ratio (first 4 components):")
    print(np.round(ev[:4], 4))
    print()

    print(f"Outputs written to: {OUT_DIR}")


if __name__ == "__main__":
    main()
