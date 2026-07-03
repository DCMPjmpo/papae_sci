# DEBUG REPORT — Blank Figure Canvas

## Confirmed Root Cause

**Bug:** `fig.add_axes([x, y, width, height])` receives values in **inches**, but matplotlib expects **figure-fraction coordinates** (0–1).

Every panel in every figure is placed outside the visible canvas.

## Evidence

### Test: add_axes with inch values (current code)

```python
fig = plt.figure(figsize=(7.008, 5.0))
ax = fig.add_axes([0, 2.5, 7.008, 2.475])
print(ax.get_position())
# → Bbox(x0=0.0, y0=2.5, x1=7.008, y1=4.975)
```

- y0 = **2.5** (= 250% of figure height above canvas)
- x1 = **7.008** (= 700% of figure width)
- y1 = **4.975** (= 497% of figure height)

The axes is placed **entirely outside the visible [0,1]×[0,1] figure box**. matplotlib clips it, rendering nothing → **blank white canvas**.

### Test: normalized values (correct)

```python
ax2 = fig.add_axes([0, 2.5/5.0, 1.0, 2.475/5.0])
print(ax2.get_position())
# → Bbox(x0=0.0, y0=0.5, x1=1.0, y1=0.995)
```

- All values within [0,1] → panel renders correctly.

### Non-white pixel count: 0 of 3,153,000 (0%)

The debug test confirms: `Fig4_debug_test.png` is **uniform white** (255,255,255) at every pixel.

---

## Affected Code (all figures)

Every `add_axes()` call in `assemble_figures.py` passes inch values instead of fractions.

### File: `figure_assembly/assemble_figures.py`

| Function | Line | Code (broken) | Correct version |
|----------|------|---------------|-----------------|
| `build_figure_1` | 148 | `[0, panel_heights[1] + gap, fig_w, panel_heights[0]]` | `[0, (panel_heights[1]+gap)/fig_h, 1.0, panel_heights[0]/fig_h]` |
| `build_figure_1` | 156 | `[0, 0, fig_w, panel_heights[1]]` | `[0, 0, 1.0, panel_heights[1]/fig_h]` |
| `build_figure_2` | 181 | `[0, bot_h + gap, fig_w, top_h]` | `[0, (bot_h+gap)/fig_h, 1.0, top_h/fig_h]` |
| `build_figure_2` | 188 | `[0, 0, fig_w*0.49, bot_h]` | `[0, 0, 0.49, bot_h/fig_h]` |
| `build_figure_2` | 194 | `[fig_w*0.51, 0, fig_w*0.49, bot_h]` | `[0.51, 0, 0.49, bot_h/fig_h]` |
| `build_figure_3` | 221 | `[0, 0, left_w, fig_h]` | `[0, 0, left_w/fig_w, 1.0]` |
| `build_figure_3` | 227 | `[left_w+gap, bot_h_right+gap, right_w, top_h_right]` | `[(left_w+gap)/fig_w, (bot_h_right+gap)/fig_h, right_w/fig_w, top_h_right/fig_h]` |
| `build_figure_3` | 233 | `[left_w+gap, 0, right_w, bot_h_right]` | `[(left_w+gap)/fig_w, 0, right_w/fig_w, bot_h_right/fig_h]` |
| `build_figure_4` | 255 | `[0, half_h+gap, fig_w, half_h]` | `[0, (half_h+gap)/fig_h, 1.0, half_h/fig_h]` |
| `build_figure_4` | 261 | `[0, 0, fig_w, half_h]` | `[0, 0, 1.0, half_h/fig_h]` |
| `build_figure_5` | 287 | `[0, bot_h+gap, fig_w, top_h]` | `[0, (bot_h+gap)/fig_h, 1.0, top_h/fig_h]` |
| `build_figure_5` | 293 | `[0, 0, half_w, bot_h]` | `[0, 0, half_w/fig_w, bot_h/fig_h]` |
| `build_figure_5` | 300 | `[half_w+gap, 0, half_w, bot_h]` | `[(half_w+gap)/fig_w, 0, half_w/fig_w, bot_h/fig_h]` |

---

## Why the pre-existing `Fig4_final.png` worked (545 KB)

The earlier `figures/Fig4_final.png` (545 KB, generated Jul 1) was created by a **different assembly method** (PIL paste/copy, not matplotlib `add_axes`). That method directly concatenated pixel arrays and did not use matplotlib at all. The current `assemble_figures.py` was written for publication-quality output but has never successfully run.

---

## Fix

Normalize every `add_axes()` argument by dividing the inch values by the figure dimensions:

```python
# Before:
ax = fig.add_axes([0, half_h + gap, fig_w, half_h])
# After:
ax = fig.add_axes([0, (half_h + gap) / fig_h, 1.0, half_h / fig_h])
```

The `assemble_supp_figures.py` script does NOT use `add_axes` — it uses PIL `Image.open().save()` and `add_axes` with matplotlib. So Supp Fig S1, S4, S5 copies are fine. But Supp Fig S2 and S3 (which also use matplotlib `add_axes`) may have the same bug.

---

## Summary

| Item | Value |
|------|-------|
| **Bug type** | Coordinate system confusion |
| **Root cause file** | `figure_assembly/assemble_figures.py` |
| **Root cause in all 13 add_axes calls** | Passes inches instead of 0–1 fractions |
| **Evidence** | `add_axes([0, 2.5, 7.008, 2.475])` → y0=2.5, x1=7.008 (both >> 1) |
| **Effect** | 100% white canvas, 0 non-white pixels |
| **Impact** | All figures Fig1–Fig5 |
| **Fix** | Divide every `w` and `h` argument by `fig_w`/`fig_h` respectively; set full-width to `1.0` |
