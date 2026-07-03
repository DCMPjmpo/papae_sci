"""
Layer-occurrence frequency across the Top50 Composite layer pairs of
five proteins (CBS, GAL4, PABP, PTEN, TEM1) in ESM-2 (33 layers).

For every protein we read top50_composite_{prot}.csv and accumulate
counts of `layer_i` and `layer_j` (both 1-based, range 1..33). Each
Top50 pair therefore contributes 2 to the running layer-occurrence
total. With 5 proteins × 50 pairs × 2 ends, the grand total is 500.

Outputs (data/p0_output/protein_clustering/):
    - layer_frequency.csv               (Layer, Count, Proportion)
    - layer_frequency_histogram.png     (publication-quality)
    - layer_frequency_band_summary.csv  (proportions for 1-8, 9-24, 25-33)
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DATA_DIR = Path("D:/文件/工作室/website/data/processed/p0_output")
OUT_DIR = Path("D:/文件/工作室/website/data/p0_output/protein_clustering")
OUT_DIR.mkdir(parents=True, exist_ok=True)

PROTEINS = ["CBS", "GAL4", "PABP", "PTEN", "TEM1"]
N_LAYERS = 33

# Band definitions (per user spec)
BANDS = {
    "Early (1-8)":  list(range(1, 9)),
    "Middle (9-24)": list(range(9, 25)),
    "Late (25-33)":  list(range(25, 34)),
}

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


def load_per_protein_pairs() -> dict:
    """Return {protein: DataFrame[Top50 pairs]}."""
    data = {}
    for prot in PROTEINS:
        csv_path = DATA_DIR / prot / f"top50_composite_{prot}.csv"
        if not csv_path.exists():
            raise FileNotFoundError(csv_path)
        df = pd.read_csv(csv_path)
        # Defensive: ensure 1-based columns are integer
        df = df.assign(
            layer_i=df["layer_i"].astype(int),
            layer_j=df["layer_j"].astype(int),
        )
        if len(df) != 50:
            print(f"  warning: {prot} has {len(df)} rows (expected 50)")
        data[prot] = df
    return data


def count_layer_occurrences(pairs_by_prot: dict) -> pd.DataFrame:
    counts = np.zeros(N_LAYERS, dtype=int)              # global
    per_prot = {p: np.zeros(N_LAYERS, dtype=int) for p in pairs_by_prot}

    for prot, df in pairs_by_prot.items():
        for col in ("layer_i", "layer_j"):
            for v in df[col].values:
                if 1 <= v <= N_LAYERS:
                    counts[v - 1] += 1
                    per_prot[prot][v - 1] += 1
                else:
                    print(f"  warning: {prot} {col}={v} outside 1..{N_LAYERS}")

    total = counts.sum()
    layer_df = pd.DataFrame({
        "Layer": np.arange(1, N_LAYERS + 1),
        "Count": counts,
        "Proportion": counts / total if total > 0 else np.zeros(N_LAYERS),
    })
    for prot in pairs_by_prot:
        layer_df[f"Count_{prot}"] = per_prot[prot]
    return layer_df


def plot_histogram(layer_df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(9.0, 4.2))

    layers = layer_df["Layer"].values
    counts = layer_df["Count"].values

    # Color bars by band — keeps the figure interpretable at a glance
    band_colors = {
        "Early (1-8)":   "#1B4965",
        "Middle (9-24)": "#9AA1A8",
        "Late (25-33)":  "#D7263D",
    }
    bar_colors = []
    for L in layers:
        for name, members in BANDS.items():
            if L in members:
                bar_colors.append(band_colors[name])
                break

    ax.bar(layers, counts, color=bar_colors, edgecolor="black", linewidth=0.4)

    ax.set_xticks(layers)
    ax.set_xticklabels(layers, rotation=0, fontsize=8)
    ax.set_xlim(0.4, N_LAYERS + 0.6)
    ax.set_xlabel("ESM-2 layer index (1-33)")
    ax.set_ylabel("Occurrences in Top50 Composite layer pairs\n(summed over 5 proteins)")
    ax.set_title("Layer-occurrence frequency across five proteins", pad=10)

    # Band separators
    for x in (8.5, 24.5):
        ax.axvline(x, color="black", linestyle=":", linewidth=0.6, alpha=0.6)

    # Legend
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=c, edgecolor="black", linewidth=0.4, label=name)
               for name, c in band_colors.items()]
    ax.legend(handles=handles, loc="upper center", ncol=3,
              bbox_to_anchor=(0.5, -0.18))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_linewidth(0.8)

    fig.savefig(OUT_DIR / "layer_frequency_histogram.png")
    fig.savefig(OUT_DIR / "layer_frequency_histogram.pdf")
    plt.close(fig)


def summarize(layer_df: pd.DataFrame) -> tuple:
    total = layer_df["Count"].sum()

    # Top-10 most frequent layers
    top10 = layer_df.sort_values("Count", ascending=False).head(10) \
                    .reset_index(drop=True)
    top10.insert(0, "Rank", np.arange(1, len(top10) + 1))

    band_rows = []
    for name, members in BANDS.items():
        sub = layer_df[layer_df["Layer"].isin(members)]
        c = int(sub["Count"].sum())
        band_rows.append({
            "Band": name,
            "Layers": f"{min(members)}-{max(members)}",
            "Count": c,
            "Proportion": c / total if total > 0 else 0.0,
        })
    band_df = pd.DataFrame(band_rows)

    return top10, band_df, int(total)


def main() -> None:
    pairs = load_per_protein_pairs()
    layer_df = count_layer_occurrences(pairs)
    layer_df.to_csv(OUT_DIR / "layer_frequency.csv", index=False)

    plot_histogram(layer_df)

    top10, band_df, total = summarize(layer_df)
    band_df.to_csv(OUT_DIR / "layer_frequency_band_summary.csv", index=False)

    print(f"Total layer-occurrences (5 proteins x 50 pairs x 2 ends): {total}")
    print()
    print("Top-10 most frequent layers:")
    print(top10[["Rank", "Layer", "Count", "Proportion"]]
          .to_string(index=False, float_format=lambda v: f"{v: .4f}"))
    print()
    print("Band proportions:")
    print(band_df.to_string(index=False, float_format=lambda v: f"{v: .4f}"))
    print()
    print(f"Wrote: {OUT_DIR / 'layer_frequency.csv'}")
    print(f"Wrote: {OUT_DIR / 'layer_frequency_histogram.png'}")
    print(f"Wrote: {OUT_DIR / 'layer_frequency_band_summary.csv'}")


if __name__ == "__main__":
    main()
