#!/usr/bin/env python3
"""
Remake_Panel_1a_Grant.py
========================
Grant revision of Panel 1a only.

Changes from the original Step05 Panel 1a:
  1. Point size tripled:  s=50  ->  s=150
  2. All points uniform warm red (#ed6a5a) instead of per-region coloring
  3. No background polygon fills (clean white background)
  4. No legend
  5. Y-axis label simplified: "SBS2 Weight" (removed "mutation count")
  6. Horizontal median dashed line removed (diagonal threshold kept)

Everything else is preserved exactly:
  - Broken x-axis logic (find_axis_break, add_break_markers)
  - Mint/teal y-axis line on the main panel
  - Dashed diagonal + horizontal threshold lines
  - Inner axes styling (add_inner_axes)
  - Figure size (16x10 with break, 14x10 without)
  - Font sizes (FONT_SIZE=28, ticks=22, etc.)
  - Point edge colors, alpha, rasterization
  - Axis labels, tick sizes, spine cleanup

Reads:  data/FIG_1/HNSC_A3_SBS2_matched_v3.tsv  (426 DIRECT-only tumors)
Saves:  data/FIG_1/FIGURE_1_PANELS/Panel_1a_A3sum_vs_SBS2_grant.pdf/.png

Usage:
    conda run -n NETWORK python Remake_Panel_1a_Grant.py

Author: Jake Lehle / Claude (2026 NMF Paper, grant revision)
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# =============================================================================
# CONFIGURATION  (matches Step05_Revised exactly)
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
OUTPUT_DIR   = os.path.join(PROJECT_ROOT, "data", "FIG_1")
PANEL_DIR    = os.path.join(OUTPUT_DIR, "FIGURE_1_PANELS")
os.makedirs(PANEL_DIR, exist_ok=True)

INPUT_PATH = os.path.join(OUTPUT_DIR, "HNSC_A3_SBS2_matched_v3.tsv")

# Colors (unchanged from Step05)
COLOR_SBS2_HIGH = "#ed6a5a"   # warm coral/red - now used for ALL points
COLOR_CREAM     = "#f4f1bb"
COLOR_NORMAL    = "#9bc1bc"   # teal/mint
COLOR_DARK_GRAY = "#4D4D4D"
COLOR_BLACK     = "#000000"

FONT_SIZE = 28

plt.rcParams.update({
    'font.size':        FONT_SIZE,
    'axes.titlesize':   FONT_SIZE,
    'axes.labelsize':   FONT_SIZE,
    'xtick.labelsize':  FONT_SIZE - 6,
    'ytick.labelsize':  FONT_SIZE - 6,
    'legend.fontsize':  FONT_SIZE - 8,
    'font.family':      'sans-serif',
    'font.sans-serif':  ['Arial', 'DejaVu Sans'],
})


# =============================================================================
# HELPER FUNCTIONS  (copied verbatim from Step05)
# =============================================================================
def log(msg=""):
    print(msg, flush=True)

def banner(title, char="="):
    log(""); log(char * 80); log(f"  {title}"); log(char * 80)

def save_fig(fig, name):
    for ext in ['pdf', 'png']:
        fig.savefig(os.path.join(PANEL_DIR, f"{name}.{ext}"),
                    dpi=300, bbox_inches='tight')
    plt.close(fig)
    log(f"  Saved: {name}.pdf/.png")

def find_axis_break(values, n_check=10, min_gap_ratio=1.5):
    sv = np.sort(values)[::-1]
    if len(sv) < 3:
        return 0, 0, False
    best_r, best_i = 1.0, -1
    for i in range(min(n_check, len(sv) - 1)):
        r = sv[i] / sv[i + 1] if sv[i + 1] > 0 else 1.0
        if r > best_r:
            best_r, best_i = r, i
    if best_r > min_gap_ratio and best_i >= 0:
        bs = sv[best_i + 1] * 1.10
        be = sv[best_i] * 0.97
        log(f"    Axis break: ratio={best_r:.2f}, {bs:.1f}-{be:.1f}")
        return bs, be, True
    return 0, 0, False

def add_inner_axes(ax, x_lo, x_hi, y_lo, y_hi,
                   y_color=COLOR_BLACK, lw=3, zorder=1):
    if x_lo <= 0 <= x_hi:
        ax.plot([0, 0], [max(0, y_lo), y_hi],
                color=y_color, lw=lw, zorder=zorder, solid_capstyle='butt')
    ax.plot([max(0, x_lo), x_hi], [0, 0],
            color=COLOR_BLACK, lw=lw, zorder=zorder, solid_capstyle='butt')

def add_break_markers(axL, axR, d=0.015):
    for kw in [dict(transform=axL.transAxes)]:
        axL.plot((1 - d, 1 + d), (-d, +d),
                 color=COLOR_BLACK, clip_on=False, lw=1.5, **kw)
        axL.plot((1 - d, 1 + d), (1 - d, 1 + d),
                 color=COLOR_BLACK, clip_on=False, lw=1.5, **kw)
    for kw in [dict(transform=axR.transAxes)]:
        axR.plot((-d, +d), (-d, +d),
                 color=COLOR_BLACK, clip_on=False, lw=1.5, **kw)
        axR.plot((-d, +d), (1 - d, 1 + d),
                 color=COLOR_BLACK, clip_on=False, lw=1.5, **kw)
    axR.spines['left'].set_visible(False)
    axR.tick_params(left=False)
    axL.spines['right'].set_visible(False)

def clean_spines(ax):
    for s in ax.spines.values():
        s.set_visible(False)


# =============================================================================
# LOAD DATA
# =============================================================================
banner("LOAD DATA")
if not os.path.exists(INPUT_PATH):
    log(f"  ERROR: Input file not found: {INPUT_PATH}")
    sys.exit(1)

final = pd.read_csv(INPUT_PATH, sep='\t')
log(f"  Loaded {len(final)} tumors from v3 matched file")

x_all = final['A3A_plus_A3B'].values
y_all = final['SBS2'].values

# Recompute threshold geometry (same logic as Step05)
median_sbs2 = final['SBS2'].median()
lm = (x_all <= 10) & (x_all > 0)
safe_slope = np.max(y_all[lm] / x_all[lm]) * 1.05 if lm.sum() > 0 else np.max(y_all) / 10
X_THR = median_sbs2 / safe_slope if safe_slope > 0 else 10
y_thr = median_sbs2

x_max = x_all.max() * 1.05
y_max = y_all.max() * 1.05

log(f"  median_sbs2={median_sbs2:.0f}, safe_slope={safe_slope:.4f}, X_THR={X_THR:.1f}")
log(f"  x_max={x_max:.1f}, y_max={y_max:.1f}")


# =============================================================================
# PANEL 1a  (GRANT REVISION)
#
# CHANGES from original:
#   - No ax.fill() calls (no background polygons)
#   - All points COLOR_SBS2_HIGH (#ed6a5a), s=100 (was s=50, per-region color)
#   - No legend
# =============================================================================
banner("Panel 1a (Grant Revision)")

bs, be, hb = find_axis_break(x_all)


def draw_1a(ax, xl, xh, is_main=True):
    """Draw Panel 1a on a single axes object.

    Identical to Step05 draw_1a EXCEPT:
      - Background polygon fills removed
      - Scatter uses uniform color + doubled size
      - No legend is added
    """
    # X-axis baseline only (green y-axis removed per PI)
    ax.plot([max(0, xl), xh], [0, 0],
            color=COLOR_BLACK, lw=3, zorder=1, solid_capstyle='butt')

    # Dashed diagonal threshold line only (horizontal median line removed per PI)
    ax.plot([0, X_THR], [0, y_thr], '--',
            color=COLOR_DARK_GRAY, lw=1, alpha=0.7, zorder=2)

    # Scatter: uniform warm red, doubled size
    ax.scatter(x_all, y_all,
               c=COLOR_SBS2_HIGH,       # <-- uniform color (was per-region)
               s=150,                    # <-- even larger (was 100, original 50)
               alpha=0.7,
               edgecolors=COLOR_BLACK,
               linewidths=0.4,
               rasterized=True,
               zorder=3)

    ax.set_xlim(xl, xh)
    ax.set_ylim(-y_max * 0.02, y_max)


if hb:
    # Broken x-axis layout (same gridspec as original)
    fig, (axM, axB) = plt.subplots(
        1, 2, sharey=True, figsize=(16, 10),
        gridspec_kw={'width_ratios': [3, 1], 'wspace': 0.04}
    )
    draw_1a(axM, -x_max * 0.02, bs, True)
    draw_1a(axB, be, x_max, False)
    clean_spines(axM)
    clean_spines(axB)
    add_break_markers(axM, axB)

    axM.set_ylabel('SBS2 Weight', fontsize=FONT_SIZE)
    fig.text(0.5, 0.01, 'A3A + A3B Expression (FPKM-UQ)',
             ha='center', fontsize=FONT_SIZE)

    # NO LEGEND (removed per PI directive)

    plt.tight_layout()
    save_fig(fig, "Panel_1a_A3sum_vs_SBS2_grant")

else:
    # No break needed
    fig, ax = plt.subplots(figsize=(14, 10))
    draw_1a(ax, -x_max * 0.02, x_max, True)
    clean_spines(ax)

    ax.set_xlabel('A3A + A3B Expression (FPKM-UQ)', fontsize=FONT_SIZE)
    ax.set_ylabel('SBS2 Weight', fontsize=FONT_SIZE)

    # NO LEGEND

    plt.tight_layout()
    save_fig(fig, "Panel_1a_A3sum_vs_SBS2_grant")


banner("DONE")
log(f"  Output in: {PANEL_DIR}")
log(f"  Files: Panel_1a_A3sum_vs_SBS2_grant.pdf / .png")
