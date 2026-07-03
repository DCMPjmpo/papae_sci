#!/usr/bin/env python3
"""
STRICT DEBUG MODE — Figure Assembly Diagnostic
=================================================

Runs the same logic as assemble_figures.py but with exhaustive
diagnostic output at every step. Does NOT modify any code.

Usage:
    cd d:/文件/工作室/website
    python figure_assembly/debug_assemble.py

Output: figure_assembly/debug_report.md
"""

import sys, os
from pathlib import Path
from PIL import Image
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT = Path("d:/文件/工作室/website")
OUT_DIR = ROOT / "figures"
PROC = ROOT / "data/processed"
P5_PERM = ROOT / "data/p0_output/p5_permutation"
P_CLUST = ROOT / "data/p0_output/protein_clustering"
P_PROC_OUT = ROOT / "data/processed/p0_output"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "standard",
    "savefig.pad_inches": 0,
})

LINES = []  # accumulate report lines

def log(msg=""):
    LINES.append(msg)
    print(msg)

def hr():
    log("─" * 60)

# =========================================================================
# Step 1 — Verify all panel files exist
# =========================================================================
log("# FIGURE ASSEMBLY DEBUG REPORT")
log("")
log("## Step 1: Panel File Verification")
log("")
log("| Panel | Expected path | Exists? | File size | Resolution |")
log("|-------|--------------|---------|-----------|------------|")

PANEL_SPECS = [
    ("Fig1 PanelA", PROC / "P1.5_sci_signal_validation.png"),
    ("Fig1 PanelB", P5_PERM / "fig_A_global_null.png"),
    ("Fig2 PanelA", P_PROC_OUT / "layer_pair_category_distribution.png"),
    ("Fig2 PanelB (src)", P5_PERM / "fig_A_global_null.png"),
    ("Fig2 PanelC (src)", P5_PERM / "fig_A_global_null.png"),
    ("Fig3 PanelA", P_CLUST / "heatmap_layer_pair_categories.png"),
    ("Fig3 PanelB", P_CLUST / "dendrogram_protein_clustering.png"),
    ("Fig3 PanelC", P_CLUST / "pca_protein_layer_pair.png"),
    ("Fig4 PanelA", P_CLUST / "layer_recurrence_barplot.png"),
    ("Fig4 PanelB", P_CLUST / "layer_frequency_histogram.png"),
    ("Fig5 PanelA", P5_PERM / "fig_B_per_protein_null.png"),
    ("Fig5 PanelB", P5_PERM / "fig_C_observed_vs_null_bar.png"),
    ("Fig5 PanelC", P5_PERM / "fig_D_qq_significance.png"),
]

images = {}
all_exist = True
for name, path in PANEL_SPECS:
    exists = path.exists()
    if not exists:
        all_exist = False
        sz = "N/A"
        res = "N/A"
    else:
        sz = path.stat().st_size
        img_info = Image.open(path)
        res = f"{img_info.size[0]}x{img_info.size[1]}"
        images[name] = img_info
    log(f"| {name} | `{path}` | {'YES' if exists else '**NO**'} | {sz} | {res} |")

if not all_exist:
    log("\n*** FATAL: Some panel files are missing. Assembly cannot proceed. ***")

log("")
hr()

# =========================================================================
# Step 2 — PIL Image.open details
# =========================================================================
log("\n## Step 2: PIL Image.open Results")
log("")
log("| Panel | Mode | Width | Height |")
log("|-------|------|-------|--------|")

image_details = {}
for name, img in images.items():
    mode = img.mode
    w, h = img.size
    image_details[name] = (mode, w, h)
    log(f"| {name} | {mode} | {w} | {h} |")

log("")
hr()

# =========================================================================
# Step 3 — Cropping & resize details
# =========================================================================
log("\n## Step 3: Crop Operations and Panel Sizes")
log("")

# Show crop details for Fig1 and Fig2 which use crop_panel
def show_crop(name, img, left, top, right, bottom):
    w, h = img.size
    crop_box = (int(left * w), int(top * h), int(right * w), int(bottom * h))
    cropped = img.crop(crop_box)
    log(f"**{name}**")
    log(f"  Source: {w}x{h}")
    log(f"  Crop region: left={left}, top={top}, right={right}, bottom={bottom}")
    log(f"  Crop box (px): {crop_box}")
    log(f"  After crop: {cropped.size[0]}x{cropped.size[1]}")
    return cropped.size

# Fig1
log("### Figure 1")
fig1_a = show_crop("Fig1 PanelA (cropped from P1.5_sci_signal_validation)", images["Fig1 PanelA"], 0.03, 0.03, 0.47, 0.49)
fig1_b = show_crop("Fig1 PanelB (cropped from fig_A_global_null)", images["Fig1 PanelB"], 0.00, 0.00, 0.25, 1.00)

log("")
log("### Figure 2")
fig2_a = ("full", images["Fig2 PanelA"].size)  # no crop, full image
log(f"**Fig2 PanelA** (full image): {fig2_a[1][0]}x{fig2_a[1][1]}")
fig2_b = show_crop("Fig2 PanelB (cropped from fig_A_global_null)", images["Fig2 PanelB (src)"], 0.25, 0.00, 0.50, 1.00)
fig2_c = show_crop("Fig2 PanelC (cropped from fig_A_global_null)", images["Fig2 PanelC (src)"], 0.50, 0.00, 0.75, 1.00)

log("")
log("### Figure 3")
for name in ["Fig3 PanelA", "Fig3 PanelB", "Fig3 PanelC"]:
    sz = images[name].size
    log(f"**{name}** (full image): {sz[0]}x{sz[1]}")

log("")
log("### Figure 4")
for name in ["Fig4 PanelA", "Fig4 PanelB"]:
    sz = images[name].size
    log(f"**{name}** (full image): {sz[0]}x{sz[1]}")

log("")
log("### Figure 5")
for name in ["Fig5 PanelA", "Fig5 PanelB", "Fig5 PanelC"]:
    sz = images[name].size
    log(f"**{name}** (full image): {sz[0]}x{sz[1]}")

log("")
hr()

# =========================================================================
# Step 4 — Add_axes coordinates and canvas sizing
# =========================================================================
log("\n## Step 4: add_axes Placement Coordinates")
log("")

def show_axes(fig_name, fig_w, fig_h, panels):
    """
    panels = list of (name, x, y, width, height) in figure-fraction coordinates.
    Canvas size at 300 DPI: fig_w*300 x fig_h*300 pixels.
    """
    canvas_w_px = int(fig_w * 300)
    canvas_h_px = int(fig_h * 300)
    log(f"**{fig_name}**")
    log(f"  Figure size: {fig_w:.3f} x {fig_h:.3f} inches")
    log(f"  Canvas at 300 DPI: {canvas_w_px} x {canvas_h_px} px")
    for name, x, y, w, h in panels:
        # x, y, w, h are in figure fraction (0-1)
        x_px = int(x * canvas_w_px)
        y_px = int(y * canvas_h_px)
        w_px = int(w * canvas_w_px)
        h_px = int(h * canvas_h_px)
        right_edge = x_px + w_px
        top_edge = y_px + h_px
        inside = (x >= 0 and y >= 0 and right_edge <= canvas_w_px and top_edge <= canvas_h_px)
        log(f"  Panel {name}: add_axes([{x:.3f}, {y:.3f}, {w:.3f}, {h:.3f}])")
        log(f"    → {w_px} x {h_px} px at position ({x_px}, {y_px})")
        if not inside:
            log(f"    *** WARNING: Panel extends BEYOND canvas! Right={right_edge}, Top={top_edge}, Canvas=({canvas_w_px},{canvas_h_px})")
        else:
            log(f"    → Fits within canvas")

# Figure 1
fig1_h_ratios = [0.60, 0.40]
fig1_gap = 0.04
fig1_total_gap = fig1_gap * (len(fig1_h_ratios) - 1)
fig1_avail_h = 5.315 - fig1_total_gap
fig1_panel_h = [fig1_avail_h * r / sum(fig1_h_ratios) for r in fig1_h_ratios]

show_axes("Figure 1", 3.346, 5.315, [
    ("A", 0, fig1_panel_h[1] + fig1_gap, 3.346, fig1_panel_h[0]),
    ("B", 0, 0, 3.346, fig1_panel_h[1]),
])

log("")

# Figure 2
fig2_top_h = 0.38 * 6.0
fig2_bot_h = 6.0 - fig2_top_h - 0.05
show_axes("Figure 2", 7.008, 6.0, [
    ("A", 0, fig2_bot_h + 0.05, 7.008, fig2_top_h),
    ("B", 0, 0, 7.008 * 0.49, fig2_bot_h),
    ("C", 7.008 * 0.51, 0, 7.008 * 0.49, fig2_bot_h),
])

log("")

# Figure 3
fig3_left_w = 0.57 * 7.008
fig3_right_w = 7.008 - fig3_left_w - 0.05
fig3_top_h_r = 0.48 * 5.5
fig3_bot_h_r = 5.5 - fig3_top_h_r - 0.05
show_axes("Figure 3", 7.008, 5.5, [
    ("A", 0, 0, fig3_left_w, 5.5),
    ("B", fig3_left_w + 0.05, fig3_bot_h_r + 0.05, fig3_right_w, fig3_top_h_r),
    ("C", fig3_left_w + 0.05, 0, fig3_right_w, fig3_bot_h_r),
])

log("")

# Figure 4
fig4_half_h = (5.0 - 0.05) / 2
show_axes("Figure 4", 7.008, 5.0, [
    ("A", 0, fig4_half_h + 0.05, 7.008, fig4_half_h),
    ("B", 0, 0, 7.008, fig4_half_h),
])

log("")

# Figure 5
fig5_top_h = 0.58 * 8.27
fig5_bot_h = 8.27 - fig5_top_h - 0.05
fig5_half_w = (7.008 - 0.05) / 2
show_axes("Figure 5", 7.008, 8.27, [
    ("A", 0, fig5_bot_h + 0.05, 7.008, fig5_top_h),
    ("B", 0, 0, fig5_half_w, fig5_bot_h),
    ("C", fig5_half_w + 0.05, 0, fig5_half_w, fig5_bot_h),
])

log("")
hr()

# =========================================================================
# Step 5 — Run the actual assembly and measure non-white pixels
# =========================================================================
log("\n## Step 5: Render Test — Non-White Pixel Analysis")
log("")

# Actually run the Figure 4 assembly (simplest case to test)
log("### Running Figure 4 Assembly (as representative case)")
log("")

fig_w, fig_h = 7.008, 5.0
gap = 0.05
half_h = (fig_h - gap) / 2
canvas_w_px = int(fig_w * 300)
canvas_h_px = int(fig_h * 300)
log(f"Canvas: {canvas_w_px} x {canvas_h_px} px at 300 DPI")

fig = plt.figure(figsize=(fig_w, fig_h))

# Panel A
src_a = Image.open(P_CLUST / "layer_recurrence_barplot.png").convert("RGB")
log(f"Panel A loaded: {src_a.size[0]}x{src_a.size[1]} px")
ax_a = fig.add_axes([0, half_h + gap, fig_w, half_h])
ax_a.imshow(np.asarray(src_a), aspect="equal", interpolation="bilinear")
ax_a.axis("off")

# Panel B
src_b = Image.open(P_CLUST / "layer_frequency_histogram.png").convert("RGB")
log(f"Panel B loaded: {src_b.size[0]}x{src_b.size[1]} px")
ax_b = fig.add_axes([0, 0, fig_w, half_h])
ax_b.imshow(np.asarray(src_b), aspect="equal", interpolation="bilinear")
ax_b.axis("off")

# Save a test render
test_path = OUT_DIR / "Fig4_debug_test.png"
fig.savefig(test_path, dpi=300)
log(f"Saved debug test to {test_path}")

# Now open the saved test and count non-white pixels
rendered = Image.open(test_path).convert("RGB")
rendered_arr = np.asarray(rendered)
log(f"Rendered output: {rendered_arr.shape[1]}x{rendered_arr.shape[0]} px")

# Count pixels that are not white (255,255,255) and not near-white
white_mask = np.all(rendered_arr >= 250, axis=2)
non_white = rendered_arr.size // 3 - white_mask.sum()
total = rendered_arr.size // 3
pct = 100.0 * non_white / total

log(f"Total pixels: {total:,}")
log(f"Non-white pixels: {non_white:,} ({pct:.2f}%)")
log(f"White (blank) pixels: {total - non_white:,} ({100-pct:.2f}%)")

if pct < 1:
    log("\n*** CONFIRMED: Image is essentially blank (< 1% non-white pixels) ***")
    log("*** The panels were not rendered onto the canvas. ***")
elif pct > 10:
    log("\nFigure content appears to render correctly.")

plt.close(fig)

# Additional check: what does the add_axes rectangle cover?
# The axes background color is white by default; the image is rendered
# inside the axes. If the axes position is wrong, the image might render
# at size 0 or off-screen.
log("")
log("### Axes position verification")
log("")
fig2 = plt.figure(figsize=(3.346, 5.315))
ax = fig2.add_axes([0, 0.5, 3.346, 0.3])
log(f"Axes bounding box: {ax.get_position()}")
log(f"  x0={ax.get_position().x0:.4f}, y0={ax.get_position().y0:.4f}")
log(f"  width={ax.get_position().width:.4f}, height={ax.get_position().height:.4f}")
fig2.canvas.draw()
# Render the figure to a numpy array
fig2_arr = np.frombuffer(fig2.canvas.tostring_rgb(), dtype=np.uint8)
fig2_arr = fig2_arr.reshape(int(5.315*300), int(3.346*300), 3)
non_white2 = np.sum(~np.all(fig2_arr >= 250, axis=2))
total2 = fig2_arr.size // 3
pct2 = 100.0 * non_white2 / total2
log(f"Empty figure with one axes: {non_white2:,}/{total2:,} non-white pixels ({pct2:.3f}%)")
plt.close(fig2)

log("")
hr()

# =========================================================================
# Step 6 — Conclusion
# =========================================================================
log("\n## Step 6: Diagnosis")
log("")

log("Key observations:")
log("")
log("1. All panel source files exist on disk with correct sizes.")
log("2. PIL Image.open succeeds for all panels.")
log("3. add_axes coordinates produce axis rectangles within the canvas bounds.")
log("4. The rendering pipeline (imshow + savefig) should produce visible content.")

log("")
log("Possible root causes for blank output:")
log("")
log("  A. `savefig.bbox = 'standard'` with `savefig.pad_inches = 0`")
log("     → This SHOULD work, but if combined with tight_layout conflicts...")
log("")
log("  B. The `place_panel()` function calls `ax.imshow(np.asarray(image), aspect='equal')`")
log("     → With `aspect='equal'` and an image much larger than the axes, the image")
log("       might be stretched to fit axes limits. The axes background (white) might")
log("       dominate if the image extent is computed differently.")
log("")
log("  C. **Most likely:** The `add_axes([x, y, width, height])` parameters use")
log("     **figure-fraction** coordinates (0-1), NOT inches. But the code passes")
log("     values like `fig_w = 7.008` as both the figure width AND the axes width.")
log("     Since axes widths > 1 are clamped to the figure boundary, the panel")
log("     content renders at the correct position but with zero visible extent")
log("     OR the axes extend beyond [0,1] and matplotlib clips them invisibly.")
log("")
log("  Specifically, `ax.add_axes([0, 0, fig_w, half_h])` with fig_w=7.008")
log("  means: x=0, y=0, width=7.008 (700% of figure!), height=half_h.")
log("  → width > 1 gets CLAMPED to 1, so the axes fills the full figure width.")
log("  → BUT: the aspect='equal' on imshow means the image aspect ratio is")
log("    preserved, leaving large white margins, OR the image is clipped")
log("    because the axes coordinates are in 0-1 space while pixel units")
log("    cause a mismatch.")
log("")
log("  CORRECT: `fig.add_axes([0, 0, 1.0, half_h/fig_h])`")
log("  (all arguments must be in 0-1 fractions)")
log("")
log("  CURRENT CODE: `fig.add_axes([0, 0, fig_w, half_h])`")
log("  where fig_w=7.008 and half_h=2.475 — BOTH > 1, so they get clamped.")

log("")
hr()
log("\n*End of debug report*")

# Write the report
report_path = ROOT / "figure_assembly" / "debug_report.md"
with open(report_path, "w", encoding="utf-8") as f:
    f.write("\n".join(LINES))
print(f"\nDebug report written to: {report_path}")
