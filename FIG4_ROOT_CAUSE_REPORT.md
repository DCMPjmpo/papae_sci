# FIG4_ROOT_CAUSE_REPORT

## Summary

The "ESM-1b layer index (1-33)" label on Fig4_final.png is a **stale-output artifact**: the source panels were generated with an older version of the plotting scripts that incorrectly said "ESM-1b", and the PNG/PDF files were never regenerated after the scripts were corrected.

---

## Step 1: Figure Assembly Mapping

| Figure file | Source panels | Panel A source | Panel B source |
|-------------|---------------|----------------|----------------|
| `figures/Fig4_final.png` (4696×5605 px) | 2 stacked panels | `layer_recurrence_barplot.png` (4539×2851 px) | `layer_frequency_histogram.png` (4696×2739 px) |

The assembly was a simple vertical concatenation. No text was added during assembly — the labels come directly from the source PNGs.

---

## Step 2: Script Analysis

| Script | Line | Current content | Comments |
|--------|------|-----------------|----------|
| `p6_scripts/layer_recurrence_top20.py` | **132** | `ax.set_xlabel("ESM-2 layer index (1-33)")` | Current code says ESM-2 |
| `p6_scripts/layer_occurrence_frequency.py` | **121** | `ax.set_xlabel("ESM-2 layer index (1-33)")` | Current code says ESM-2 |
| `p0_p5/` and `p6_scripts/` | All | No Python file anywhere in the project currently contains the string `"ESM-1b"` in an axis label or `set_xlabel` call. | All current scripts show ESM-2. |

**No script currently writes "ESM-1b" to any axis label.**

---

## Step 3: Source Panel Verification

**PDF text extraction (definitive proof).** The PDF files generated alongside the PNGs contain embedded text:

| File | Extracted text | Timestamp |
|------|---------------|-----------|
| `data/p0_output/protein_clustering/layer_recurrence_barplot.pdf` | `"ESM-1b layer index (1-33)"` | Jun 28 15:05 |
| `data/p0_output/protein_clustering/layer_frequency_histogram.pdf` | `"ESM-1b layer index (1-33)"` | Jun 28 14:45 |

Both PDFs confirm the label reads **"ESM-1b"**, not "ESM-2".

---

## Step 4: Timeline (Root Cause)

```
Jun 28 14:45    layer_frequency_histogram.png + .pdf  GENERATED  ← script said "ESM-1b"
Jun 28 15:05    layer_recurrence_barplot.png + .pdf    GENERATED  ← script said "ESM-1b"

Jun 29 10:13    layer_recurrence_top20.py               MODIFIED   ← "ESM-1b" → "ESM-2" fix
Jun 29 10:13    layer_occurrence_frequency.py           MODIFIED   ← "ESM-1b" → "ESM-2" fix

Jul  1 21:44    figures/Fig4.png                         ASSEMBLED  ← from OLD PNGs (still "ESM-1b")
Jul  1 22:06    figures/Fig4_final.png                   ASSEMBLED  ← from OLD PNGs (still "ESM-1b")
```

**Diagnosis:**
1. The scripts originally contained `ax.set_xlabel("ESM-1b layer index (1-33)")` — a **script error**
2. The PNGs and PDFs were generated with this wrong label on **Jun 28**
3. The scripts were edited to `"ESM-2"` on **Jun 29** — the source-code fix was applied
4. But the **PNGs were never regenerated** after the fix — the stale outputs persisted
5. The figure assembly on **Jul 1** consumed the stale PNGs, propagating the error into the composites

---

## Root Cause Classification

**Primary:** **D. Old-version PNGs used in assembly** (stale-output error)

The scripts were fixed, but the `.png` and `.pdf` output files under `data/p0_output/protein_clustering/` were not re-generated. All downstream consumers (Fig4 composites, the `calibration_crops/`, the `figures/` directory) therefore inherited the incorrect label.

**Secondary:** **A. Script error (historical only)** — the Jun-28 version of both scripts used `"ESM-1b"` instead of `"ESM-2"`.

---

## Affected Files (must be regenerated)

| File | Label shown | Must fix |
|------|------------|----------|
| `data/p0_output/protein_clustering/layer_recurrence_barplot.png` | ESM-1b | ✅ Re-run `layer_recurrence_top20.py` |
| `data/p0_output/protein_clustering/layer_recurrence_barplot.pdf` | ESM-1b | ✅ Re-run `layer_recurrence_top20.py` |
| `data/p0_output/protein_clustering/layer_frequency_histogram.png` | ESM-1b | ✅ Re-run `layer_occurrence_frequency.py` |
| `data/p0_output/protein_clustering/layer_frequency_histogram.pdf` | ESM-1b | ✅ Re-run `layer_occurrence_frequency.py` |
| `figure_assembly/calibration_crops/Fig4_PanelA.png` | ESM-1b (derived) | Will be fixed transitively |
| `figure_assembly/calibration_crops/Fig4_PanelB.png` | ESM-1b (derived) | Will be fixed transitively |
| `figures/Fig4.png` | ESM-1b (derived) | Will be fixed transitively |
| `figures/Fig4_final.png` | ESM-1b (derived) | Will be fixed transitively |

**Fix action:** Re-run the two scripts, then re-assemble Fig4 composites. No code modification needed — the scripts already say `"ESM-2"`.

---

## Other Figures (Not Affected)

| Panel | Script x-axis label | Source PNG status |
|-------|-------------------|-------------------|
| Fig 1a — SCI distribution | No layer-index label, unaffected | Clean |
| Fig 1b — T1 null | No layer-index label, unaffected | Clean |
| Fig 2a — Category distribution | No layer-index label, unaffected | Clean |
| Fig 2b — T2 null | No layer-index label, unaffected | Clean |
| Fig 2c — T3 null | No layer-index label, unaffected | Clean |
| Fig 3a — Heatmap | No layer-index label, unaffected | Clean |
| Fig 3b — Dendrogram | No layer-index label, unaffected | Clean |
| Fig 3c — PCA | No layer-index label, unaffected | Clean |
| **Fig 4a** | **"ESM-1b layer index (1-33)"** | **❌ Affected** |
| **Fig 4b** | **"ESM-1b layer index (1-33)"** | **❌ Affected** |
| Fig 5a — Per-protein null grid | No layer-index label, unaffected | Clean |
| Fig 5b — z-score bar chart | No layer-index label, unaffected | Clean |
| Fig 5c — QQ plot | No layer-index label, unaffected | Clean |

Only Figure 4 is affected.
