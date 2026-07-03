"""
Per-layer occurrence frequency among the Top20 cross-protein recurrent
layer pairs (ranked by Cross_Protein_Count).

Input  : data/processed/p0_output/cross_protein_universal_pairs.csv
         Columns: Layer_Pair_0based, Cross_Protein_Count,
                  layer_i_0, layer_j_0, layer_i_1, layer_j_1
         (1-based columns used here; range 1..33.)

Outputs (data/p0_output/protein_clustering/):
    - layer_recurrence_frequency.csv      (Layer, Count, Proportion, Band)
    - layer_recurrence_top20_pairs.csv    (the Top20 list itself, for traceability)
    - layer_recurrence_band_summary.csv   (Band, Count, Proportion)
    - layer_recurrence_barplot.{png,pdf}  (publication-quality bar chart)
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

INPUT_PATH = Path("D:/文件/工作室/website/data/processed/p0_output/cross_protein_universal_pairs.csv")
OUT_DIR = Path("D:/文件/工作室/website/data/p0_output/protein_clustering")
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_LAYERS = 33
TOP_N = 20

BANDS = {
    "Early (1-8)":   (range(1, 9),   "#1B4965"),
    "Middle (9-24)": (range(9, 25),  "#9AA1A8"),
    "Late (25-33)":  (range(25, 34), "#D7263D"),
}


def layer_band(L: int) -> str:
    for name, (members, _) in BANDS.items():
        if L in members:
            return name
    return "Unknown"


def layer_color(L: int) -> str:
    for _, (members, color) in BANDS.items():
        if L in members:
            return color
    return "#888888"


plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "axes.linewidth": 1.0,
    "xtick.labelsize": 9,
    "ytick.labelsize": 10,
    "legend.fontsize": 9,
    "legend.frameon": False,
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def main() -> None:
    df = pd.read_csv(INPUT_PATH)

    # Rank by Cross_Protein_Count desc; preserve input order within ties.
    df_sorted = df.sort_values(
        "Cross_Protein_Count", ascending=False, kind="mergesort",
    ).reset_index(drop=True)
    top20 = df_sorted.head(TOP_N).copy()

    # 1-based layer indices
    top20["layer_i_1"] = top20["layer_i_1"].astype(int)
    top20["layer_j_1"] = top20["layer_j_1"].astype(int)

    # Persist Top20 list for traceability
    top20_out = top20[["layer_i_1", "layer_j_1", "Cross_Protein_Count"]].rename(
        columns={"layer_i_1": "Layer_i", "layer_j_1": "Layer_j"})
    top20_out.insert(0, "Rank", np.arange(1, len(top20_out) + 1))
    top20_out["Pair_Category"] = [
        f"{layer_band(i).split()[0][0]}{layer_band(j).split()[0][0]}"
        for i, j in zip(top20_out["Layer_i"], top20_out["Layer_j"])
    ]
    top20_out.to_csv(OUT_DIR / "layer_recurrence_top20_pairs.csv", index=False)

    # ---- Count individual layer occurrences in the Top20 ----
    counts = np.zeros(N_LAYERS, dtype=int)
    for col in ("layer_i_1", "layer_j_1"):
        for v in top20[col].values:
            if 1 <= v <= N_LAYERS:
                counts[v - 1] += 1

    total = int(counts.sum())                              # = 2 * TOP_N = 40
    layer_df = pd.DataFrame({
        "Layer": np.arange(1, N_LAYERS + 1),
        "Count": counts,
        "Proportion": counts / total if total else np.zeros(N_LAYERS),
    })
    layer_df["Band"] = layer_df["Layer"].apply(layer_band)
    layer_df.to_csv(OUT_DIR / "layer_recurrence_frequency.csv", index=False)

    # ---- Band summary ----
    band_rows = []
    for name, (members, _) in BANDS.items():
        members = list(members)
        c = int(layer_df.loc[layer_df["Layer"].isin(members), "Count"].sum())
        band_rows.append({
            "Band": name,
            "Layers": f"{min(members)}-{max(members)}",
            "Count": c,
            "Proportion": c / total if total else 0.0,
        })
    band_df = pd.DataFrame(band_rows)
    band_df.to_csv(OUT_DIR / "layer_recurrence_band_summary.csv", index=False)

    # ---- Bar plot ----
    fig, ax = plt.subplots(figsize=(9.0, 4.2))
    bar_colors = [layer_color(L) for L in layer_df["Layer"].values]
    ax.bar(layer_df["Layer"], layer_df["Count"],
           color=bar_colors, edgecolor="black", linewidth=0.4)

    ax.set_xticks(np.arange(1, N_LAYERS + 1))
    ax.set_xticklabels(np.arange(1, N_LAYERS + 1), fontsize=8)
    ax.set_xlim(0.4, N_LAYERS + 0.6)
    ax.set_xlabel("ESM-2 layer index (1-33)")
    ax.set_ylabel(f"Occurrences among Top{TOP_N} recurrent layer pairs")
    ax.set_title(
        f"Per-layer recurrence frequency\n(Top{TOP_N} cross-protein layer pairs)",
        pad=10,
    )

    for x in (8.5, 24.5):
        ax.axvline(x, color="black", linestyle=":", linewidth=0.6, alpha=0.6)

    # Annotate non-zero bars with their count for at-a-glance reading
    for L, c in zip(layer_df["Layer"], layer_df["Count"]):
        if c > 0:
            ax.text(L, c + 0.15, str(c), ha="center", va="bottom",
                    fontsize=8, color="black")

    handles = [Patch(facecolor=col, edgecolor="black", linewidth=0.4, label=name)
               for name, (_, col) in BANDS.items()]
    ax.legend(handles=handles, loc="upper center", ncol=3,
              bbox_to_anchor=(0.5, -0.18))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_linewidth(0.8)
    ax.margins(y=0.12)

    fig.savefig(OUT_DIR / "layer_recurrence_barplot.png")
    fig.savefig(OUT_DIR / "layer_recurrence_barplot.pdf")
    plt.close(fig)

    # ---- Top-10 most recurrent layers ----
    top10 = layer_df.sort_values("Count", ascending=False).head(10).reset_index(drop=True)
    top10.insert(0, "Rank", np.arange(1, len(top10) + 1))

    print(f"Top{TOP_N} recurrent pairs preview:")
    print(top20_out.to_string(index=False))
    print()
    print(f"Total layer-occurrences in Top{TOP_N} pairs: {total}")
    print()
    print("Band proportions:")
    print(band_df.to_string(index=False, float_format=lambda v: f"{v: .4f}"))
    print()
    print("Top-10 most recurrent layers:")
    print(top10[["Rank", "Layer", "Count", "Proportion", "Band"]]
          .to_string(index=False, float_format=lambda v: f"{v: .4f}"))


if __name__ == "__main__":
    main()
