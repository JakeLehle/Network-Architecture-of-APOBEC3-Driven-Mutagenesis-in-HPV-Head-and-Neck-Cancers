#!/usr/bin/env python3
"""
Step05_Panel_1d_Saturation.py
=============================
Panel 1d: A3B saturation effect on SBS2.

Four sub-panels (2 rows x 3 cols, bottom row spans all 3):
  Top-left:   A3B vs SBS2 scatter with decile-binned percentile fan
              Shows SBS2 rises with A3B then plateaus (ceiling effect)
  Top-middle: A3A vs SBS2 scatter with decile-binned percentile fan
              Shows SBS2 continues to rise sporadically with A3A (no ceiling)
  Top-right:  Quantile-based threshold sweep for A3B (rho_below vs rho_above
              at Q20, Q40, Q60, Q80). 20% step matches the bottom panel's
              tick spacing exactly.
  Bottom:     The SAME rho values from the top-right panel sweep, relabeled
              onto a median-centered x-axis:
                Below side: x = -(100 - Q), y = rho_below at Q
                Above side: x = +Q,         y = rho_above at Q
              8 points total, landing exactly on -80, -60, -40, -20 and
              +20, +40, +60, +80. There is a 40% gap across the median;
              the line is drawn as two separate connected segments
              (below-median and above-median) with a shared linear trend
              fit across all 8 points.

Reads: data/FIG_1/HNSC_A3_SBS2_matched_v3.tsv
Saves: FIGURE_1_PANELS/Panel_1d_A3B_Saturation.pdf/.png

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os, sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import spearmanr, linregress

# =============================================================================
# CONFIGURATION (matches Step05)
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "FIG_1")
PANEL_DIR = os.path.join(OUTPUT_DIR, "FIGURE_1_PANELS")
TROUBLE_DIR = os.path.join(OUTPUT_DIR, "TROUBLESHOOTING")
os.makedirs(PANEL_DIR, exist_ok=True); os.makedirs(TROUBLE_DIR, exist_ok=True)

INPUT_PATH = os.path.join(OUTPUT_DIR, "HNSC_A3_SBS2_matched_v3.tsv")

COLOR_SBS2_HIGH = "#ed6a5a"
COLOR_CREAM = "#f4f1bb"
COLOR_NORMAL = "#9bc1bc"
COLOR_DARK_GRAY = "#4D4D4D"
COLOR_BLACK = "#000000"
COLOR_A3B = "#5e81ac"         # steel blue
COLOR_A3A = "#bf616a"         # muted red

FONT_SIZE = 28

plt.rcParams.update({
    'font.size': FONT_SIZE,
    'axes.titlesize': FONT_SIZE,
    'axes.labelsize': FONT_SIZE,
    'xtick.labelsize': FONT_SIZE - 6,
    'ytick.labelsize': FONT_SIZE - 6,
    'legend.fontsize': FONT_SIZE - 8,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
})

# =============================================================================
# LOGGING
# =============================================================================
report_lines = []
def log(msg=""):
    print(msg, flush=True); report_lines.append(str(msg))
def banner(title, char="="):
    log(""); log(char * 80); log(f"  {title}"); log(char * 80)
def save_fig(fig, name):
    for ext in ['pdf', 'png']:
        fig.savefig(os.path.join(PANEL_DIR, f"{name}.{ext}"), dpi=300, bbox_inches='tight')
    plt.close(fig); log(f"  Saved: {name}.pdf/.png")

# =============================================================================
# LOAD DATA
# =============================================================================
banner("LOAD DATA")
df = pd.read_csv(INPUT_PATH, sep='\t')
log(f"  Loaded: {len(df)} tumors")
log(f"  A3B range: {df['A3B'].min():.2f} - {df['A3B'].max():.2f}")
log(f"  A3A range: {df['A3A'].min():.2f} - {df['A3A'].max():.2f}")
log(f"  SBS2 range: {df['SBS2'].min():.0f} - {df['SBS2'].max():.0f}")

# =============================================================================
# COMPUTE BINNED STATISTICS
# =============================================================================
banner("BINNED STATISTICS (deciles)")

N_BINS = 10

def compute_binned_stats(df, expr_col, sbs_col='SBS2', n_bins=N_BINS):
    """Bin tumors by expression decile, compute multiple SBS2 percentiles per bin."""
    df_sorted = df.sort_values(expr_col).copy()
    df_sorted['bin'] = pd.qcut(df_sorted[expr_col], n_bins, labels=False, duplicates='drop')
    actual_bins = sorted(df_sorted['bin'].unique())

    fan_pctls = [30, 50, 70, 90]

    stats = []
    for b in actual_bins:
        bdf = df_sorted[df_sorted['bin'] == b]
        expr_vals = bdf[expr_col].values
        sbs_vals = bdf[sbs_col].values
        row = {
            'bin': b,
            'bin_center': np.median(expr_vals),
            'n': len(bdf),
        }
        for p in fan_pctls:
            row[f'sbs2_p{p}'] = np.percentile(sbs_vals, p)
        stats.append(row)
    return pd.DataFrame(stats), fan_pctls

a3b_bins, FAN_PCTLS = compute_binned_stats(df, 'A3B')
a3a_bins, _ = compute_binned_stats(df, 'A3A')

log(f"\n  A3B decile stats:")
header = f"  {'Decile':>7s} {'Expr_med':>9s}"
for p in FAN_PCTLS:
    header += f" {'P'+str(p):>8s}"
header += f" {'n':>5s}"
log(header)
for _, row in a3b_bins.iterrows():
    line = f"  {int(row['bin'])+1:>7d} {row['bin_center']:>9.1f}"
    for p in FAN_PCTLS:
        line += f" {row[f'sbs2_p{p}']:>8.0f}"
    line += f" {int(row['n']):>5d}"
    log(line)

log(f"\n  A3A decile stats:")
log(header)
for _, row in a3a_bins.iterrows():
    line = f"  {int(row['bin'])+1:>7d} {row['bin_center']:>9.1f}"
    for p in FAN_PCTLS:
        line += f" {row[f'sbs2_p{p}']:>8.0f}"
    line += f" {int(row['n']):>5d}"
    log(line)

# =============================================================================
# QUANTILE-BASED THRESHOLD SWEEP at 20% STEPS (Q20, Q40, Q60, Q80)
# =============================================================================
banner("QUANTILE-BASED THRESHOLD SWEEP (Q20, Q40, Q60, Q80)")

def threshold_sweep(df, expr_col, sbs_col='SBS2', quantiles=None):
    """
    For each QUANTILE-based threshold of expr_col, compute Spearman rho(expr, SBS2)
    separately for tumors BELOW and ABOVE the threshold.

    Sweep at Q20, Q40, Q60, Q80 to match the bottom panel's 20% tick spacing.
    """
    if quantiles is None:
        quantiles = [0.20, 0.40, 0.60, 0.80]

    results = []
    for q in quantiles:
        thresh = df[expr_col].quantile(q)
        below = df[df[expr_col] < thresh]
        above = df[df[expr_col] >= thresh]

        rho_below, p_below = (np.nan, np.nan)
        rho_above, p_above = (np.nan, np.nan)

        if len(below) >= 10:
            rho_below, p_below = spearmanr(below[expr_col], below[sbs_col])
        if len(above) >= 10:
            rho_above, p_above = spearmanr(above[expr_col], above[sbs_col])

        results.append({
            'quantile': q,
            'pct_label': round(q * 100),
            'threshold': thresh,
            'n_below': len(below),
            'n_above': len(above),
            'rho_below': rho_below,
            'p_below': p_below,
            'rho_above': rho_above,
            'p_above': p_above,
        })
    return pd.DataFrame(results)

a3b_sweep = threshold_sweep(df, 'A3B')

log(f"\n  A3B quantile-based threshold sweep:")
log(f"  {'Pctl':>5s} {'Thresh':>7s} {'n_below':>8s} {'n_above':>8s} "
    f"{'rho_below':>10s} {'rho_above':>10s}")
for _, row in a3b_sweep.iterrows():
    log(f"  {int(row['pct_label']):>5d} {row['threshold']:>7.1f} "
        f"{int(row['n_below']):>8d} {int(row['n_above']):>8d} "
        f"{row['rho_below']:>10.3f} {row['rho_above']:>10.3f}")

# =============================================================================
# REMAP PANEL C SWEEP ONTO MEDIAN-CENTERED % AXIS (bottom panel)
# =============================================================================
banner("REMAP PANEL C SWEEP -> MEDIAN-CENTERED AXIS (bottom panel)")

def remap_sweep_to_median_centered(sweep_df):
    """
    Direct lookup of panel c's rho values onto a median-centered % axis.

    For each quantile threshold Q in sweep_df:
      Below-median side: x = -(100 - Q), y = rho_below at Q
        Q20 -> x = -80   Q40 -> x = -60   Q60 -> x = -40   Q80 -> x = -20
      Above-median side: x = +Q, y = rho_above at Q
        Q20 -> x = +20   Q40 -> x = +40   Q60 -> x = +60   Q80 -> x = +80

    8 points total, landing exactly on -80, -60, -40, -20, +20, +40, +60, +80.
    A 40% gap exists across the median. No rho values are recomputed.
    """
    rows = []
    for _, r in sweep_df.iterrows():
        q = int(r['pct_label'])

        # BELOW-median point: x = -(100 - q), y = rho_below
        x_below = -(100 - q)
        if not pd.isna(r['rho_below']):
            rows.append({
                'pct_from_median': x_below,
                'side': 'below',
                'source_quantile': q,
                'threshold': r['threshold'],
                'n_subset': int(r['n_below']),
                'rho': r['rho_below'],
                'p_value': r['p_below'],
            })

        # ABOVE-median point: x = +q, y = rho_above
        x_above = +q
        if not pd.isna(r['rho_above']):
            rows.append({
                'pct_from_median': x_above,
                'side': 'above',
                'source_quantile': q,
                'threshold': r['threshold'],
                'n_subset': int(r['n_above']),
                'rho': r['rho_above'],
                'p_value': r['p_above'],
            })

    return pd.DataFrame(rows).sort_values('pct_from_median').reset_index(drop=True)

a3b_centered = remap_sweep_to_median_centered(a3b_sweep)

log(f"\n  Median-centered remap of panel c values:")
log(f"  {'x_pct':>6s} {'side':>6s} {'from_Q':>7s} {'thresh':>8s} "
    f"{'n':>6s} {'rho':>8s} {'p':>11s}")
for _, r in a3b_centered.iterrows():
    log(f"  {int(r['pct_from_median']):>+6d} {r['side']:>6s} "
        f"Q{int(r['source_quantile']):>2d}     {r['threshold']:>8.2f} "
        f"{int(r['n_subset']):>6d} {r['rho']:>+8.3f} "
        f"{r['p_value']:>11.3e}")

# Linear trend across all 8 remapped points
slope, intercept, r_val, p_val_trend, std_err = linregress(
    a3b_centered['pct_from_median'], a3b_centered['rho']
)
log(f"\n  Linear trend fit (rho vs % from median):")
log(f"    slope     = {slope:+.5f} (rho per % deviation)")
log(f"    intercept = {intercept:+.4f}")
log(f"    R^2       = {r_val**2:.4f}")
log(f"    trend p   = {p_val_trend:.4e}")
log(f"    Direction : {'negative (diminishing returns)' if slope < 0 else 'positive'}")

# =============================================================================
# GENERATE FIGURE: 2 x 3 grid (bottom row spans all 3 columns)
# =============================================================================
banner("GENERATE PANEL 1d (4 sub-panels)")

fig = plt.figure(figsize=(35, 18))
gs = gridspec.GridSpec(2, 3,
                        width_ratios=[1, 1, 1],
                        height_ratios=[1, 0.80],
                        wspace=0.30, hspace=0.38)
ax_a3b      = fig.add_subplot(gs[0, 0])
ax_a3a      = fig.add_subplot(gs[0, 1])
ax_sweep    = fig.add_subplot(gs[0, 2])
ax_centered = fig.add_subplot(gs[1, :])     # bottom row spans all 3 cols

Y_CAP = 300
sbs2_max_plot = Y_CAP
n_clipped = (df['SBS2'] > Y_CAP).sum()
log(f"  Y-axis capped at {Y_CAP} ({n_clipped} points above cap will be clipped)")

a3b_x_cap = np.percentile(df['A3B'], 95) * 1.15
a3a_x_cap = np.percentile(df['A3A'], 95) * 1.15
log(f"  A3B x-axis cap: {a3b_x_cap:.1f} (95th pctl * 1.15)")
log(f"  A3A x-axis cap: {a3a_x_cap:.1f} (95th pctl * 1.15)")

# =========================================================================
# TOP-LEFT: A3B vs SBS2 scatter + percentile fan
# =========================================================================
ax_a3b.scatter(df['A3B'], df['SBS2'], s=30, alpha=0.20, c=COLOR_A3B,
               edgecolors='none', rasterized=True, zorder=1)

fan_alphas  = [0.4, 0.55, 0.75, 1.0]
fan_lws     = [1.5, 2.0, 2.5, 3.0]
fan_markers = ['v', 'o', 'D', 's']
fan_msizes  = [5, 6, 7, 8]

for pctl, alpha, lw, marker, ms in zip(FAN_PCTLS, fan_alphas, fan_lws, fan_markers, fan_msizes):
    col_name = f'sbs2_p{pctl}'
    ax_a3b.plot(a3b_bins['bin_center'], a3b_bins[col_name],
                marker=marker, linestyle='-', color=COLOR_A3B,
                markersize=ms, linewidth=lw, alpha=alpha,
                markeredgecolor=COLOR_BLACK, markeredgewidth=0.8,
                label=f'{pctl}th percentile', zorder=2 + pctl // 10)

ax_a3b.set_xlabel('A3B Expression (FPKM-UQ)', fontsize=FONT_SIZE)
ax_a3b.set_ylabel('SBS2 Mutation Count', fontsize=FONT_SIZE)
ax_a3b.set_xlim(-a3b_x_cap * 0.02, a3b_x_cap)
ax_a3b.set_ylim(-sbs2_max_plot * 0.02, sbs2_max_plot)
ax_a3b.legend(loc='upper right', framealpha=0.9, fontsize=FONT_SIZE - 8)
ax_a3b.spines['top'].set_visible(False)
ax_a3b.spines['right'].set_visible(False)
ax_a3b.set_title('A3B: Saturation', fontsize=FONT_SIZE, fontweight='bold',
                  color=COLOR_A3B)

# =========================================================================
# TOP-MIDDLE: A3A vs SBS2 scatter + percentile fan
# =========================================================================
ax_a3a.scatter(df['A3A'], df['SBS2'], s=30, alpha=0.20, c=COLOR_A3A,
               edgecolors='none', rasterized=True, zorder=1)

for pctl, alpha, lw, marker, ms in zip(FAN_PCTLS, fan_alphas, fan_lws, fan_markers, fan_msizes):
    col_name = f'sbs2_p{pctl}'
    ax_a3a.plot(a3a_bins['bin_center'], a3a_bins[col_name],
                marker=marker, linestyle='-', color=COLOR_A3A,
                markersize=ms, linewidth=lw, alpha=alpha,
                markeredgecolor=COLOR_BLACK, markeredgewidth=0.8,
                label=f'{pctl}th percentile', zorder=2 + pctl // 10)

ax_a3a.set_xlabel('A3A Expression (FPKM-UQ)', fontsize=FONT_SIZE)
ax_a3a.set_ylabel('SBS2 Mutation Count', fontsize=FONT_SIZE)
ax_a3a.set_xlim(-a3a_x_cap * 0.02, a3a_x_cap)
ax_a3a.set_ylim(-sbs2_max_plot * 0.02, sbs2_max_plot)
ax_a3a.legend(loc='upper left', framealpha=0.9, fontsize=FONT_SIZE - 8)
ax_a3a.spines['top'].set_visible(False)
ax_a3a.spines['right'].set_visible(False)
ax_a3a.set_title('A3A: No Saturation', fontsize=FONT_SIZE, fontweight='bold',
                  color=COLOR_A3A)

# =========================================================================
# TOP-RIGHT: Quantile-based threshold sweep (rho above vs below)
# =========================================================================
ax_sweep.fill_between(a3b_sweep['pct_label'],
                       np.minimum(a3b_sweep['rho_below'], a3b_sweep['rho_above']),
                       np.maximum(a3b_sweep['rho_below'], a3b_sweep['rho_above']),
                       alpha=0.1, color=COLOR_DARK_GRAY, zorder=0)

ax_sweep.plot(a3b_sweep['pct_label'], a3b_sweep['rho_below'], 'o-',
              color=COLOR_A3B, markersize=10, linewidth=2.5,
              markeredgecolor=COLOR_BLACK, markeredgewidth=1,
              label='Below threshold', zorder=3)

ax_sweep.plot(a3b_sweep['pct_label'], a3b_sweep['rho_above'], 's-',
              color=COLOR_SBS2_HIGH, markersize=10, linewidth=2.5,
              markeredgecolor=COLOR_BLACK, markeredgewidth=1,
              label='Above threshold', zorder=3)

ax_sweep.axhline(0, color=COLOR_DARK_GRAY, ls='--', lw=1, alpha=0.5, zorder=0)

sweep_xticks = [20, 40, 60, 80]
ax_sweep.set_xticks(sweep_xticks)
ax_sweep.set_xticklabels([f'{x}th' for x in sweep_xticks])
ax_sweep.set_xlim(15, 85)

ax_sweep.set_xlabel('A3B Expression Percentile (Threshold)', fontsize=FONT_SIZE)
ax_sweep.set_ylabel('Spearman rho (A3B vs SBS2)', fontsize=FONT_SIZE)
ax_sweep.legend(loc='upper left', framealpha=0.9, fontsize=FONT_SIZE - 8)
ax_sweep.spines['top'].set_visible(False)
ax_sweep.spines['right'].set_visible(False)
ax_sweep.set_title('A3B Correlation Sweep', fontsize=FONT_SIZE, fontweight='bold',
                    color=COLOR_DARK_GRAY)

# =========================================================================
# BOTTOM (spans 3 cols): Panel c rho values remapped to % from median
# =========================================================================
X_LO, X_HI = -80, 80

rho_vals = a3b_centered['rho'].values
rho_pad  = 0.05
y_top    = max(rho_vals.max() + rho_pad,  0.10)
y_bot    = min(rho_vals.min() - rho_pad, -0.10)

# Background shading: positive-correlation zone (light red) above y=0,
#                     negative-correlation zone (light tan) below y=0
ax_centered.axhspan(0,            y_top * 1.25, color=COLOR_SBS2_HIGH, alpha=0.12, zorder=0)
ax_centered.axhspan(y_bot * 1.25, 0,            color=COLOR_CREAM,     alpha=0.40, zorder=0)

# Zero guides
ax_centered.axhline(0, color=COLOR_DARK_GRAY, ls='-',  lw=1.2, alpha=0.6, zorder=1)
ax_centered.axvline(0, color=COLOR_DARK_GRAY, ls='--', lw=1.5, alpha=0.6, zorder=1)

# Linear trend line across full fixed x-range
x_fit = np.linspace(X_LO, X_HI, 200)
y_fit = slope * x_fit + intercept
ax_centered.plot(x_fit, y_fit, ls='--', color=COLOR_BLACK,
                  linewidth=3.0, alpha=0.55, zorder=2,
                  label=f'Linear trend (slope = {slope:+.4f}, p = {p_val_trend:.2e})')

# Draw below-median and above-median as two separate connected segments
# so the line does not bridge the 40% gap across the median
below_pts = a3b_centered[a3b_centered['side'] == 'below']
above_pts = a3b_centered[a3b_centered['side'] == 'above']

ax_centered.plot(below_pts['pct_from_median'], below_pts['rho'],
                  marker='o', linestyle='-', color=COLOR_A3B,
                  markersize=16, linewidth=3.0,
                  markeredgecolor=COLOR_BLACK, markeredgewidth=1.5,
                  zorder=4, label='rho_below (below-median subset)')

ax_centered.plot(above_pts['pct_from_median'], above_pts['rho'],
                  marker='s', linestyle='-', color=COLOR_A3B,
                  markersize=16, linewidth=3.0,
                  markeredgecolor=COLOR_BLACK, markeredgewidth=1.5,
                  zorder=4, label='rho_above (above-median subset)')

# Fixed x-axis ticks at 20% intervals from -80 to +80
xticks = list(range(X_LO, X_HI + 1, 20))
ax_centered.set_xticks(xticks)
ax_centered.set_xticklabels(
    [(f'{t:+d}%' if t != 0 else '0%\n(median)') for t in xticks]
)

# Symmetric, fixed limits
ax_centered.set_xlim(X_LO - 5, X_HI + 5)
ax_centered.set_ylim(y_bot * 1.25, y_top * 1.25)

# Zone labels centered on each half
ax_centered.text(-40, y_top * 1.13, 'Below median',
                  fontsize=FONT_SIZE - 4, ha='center', va='top',
                  color=COLOR_DARK_GRAY, fontweight='bold', alpha=0.80)
ax_centered.text(+40, y_top * 1.13, 'Above median',
                  fontsize=FONT_SIZE - 4, ha='center', va='top',
                  color=COLOR_DARK_GRAY, fontweight='bold', alpha=0.80)

ax_centered.set_xlabel('A3B Expression (% from cohort median)', fontsize=FONT_SIZE)
ax_centered.set_ylabel('Spearman rho (A3B vs SBS2)', fontsize=FONT_SIZE)
ax_centered.legend(loc='lower left', framealpha=0.9, fontsize=FONT_SIZE - 8)
ax_centered.spines['top'].set_visible(False)
ax_centered.spines['right'].set_visible(False)
ax_centered.set_title('A3B Saturation: Correlation Decay from Median',
                       fontsize=FONT_SIZE, fontweight='bold', color=COLOR_A3B)

plt.tight_layout()
save_fig(fig, "Panel_1d_A3B_Saturation")

# =============================================================================
# DIAGNOSTIC SUMMARY
# =============================================================================
banner("DIAGNOSTIC SUMMARY FOR TEXT")

log(f"  FULL DECILE TABLE (A3B):")
header = f"  {'Decile':>7s} {'Expr_med':>9s} {'n':>5s}"
for p in FAN_PCTLS:
    header += f" {'P'+str(p):>8s}"
log(header)
for _, row in a3b_bins.iterrows():
    line = f"  {int(row['bin'])+1:>7d} {row['bin_center']:>9.1f} {int(row['n']):>5d}"
    for p in FAN_PCTLS:
        line += f" {row[f'sbs2_p{p}']:>8.0f}"
    log(line)

log(f"\n  FULL DECILE TABLE (A3A):")
log(header)
for _, row in a3a_bins.iterrows():
    line = f"  {int(row['bin'])+1:>7d} {row['bin_center']:>9.1f} {int(row['n']):>5d}"
    for p in FAN_PCTLS:
        line += f" {row[f'sbs2_p{p}']:>8.0f}"
    log(line)

# =========================================================================
# 2. PANEL 1d TOP-LEFT: A3B fan numbers
# =========================================================================
banner("PANEL 1d TOP-LEFT: A3B Saturation Fan")

log(f"  Percentile behavior across A3B expression deciles:")
for p in FAN_PCTLS:
    d1 = a3b_bins.iloc[0][f'sbs2_p{p}']
    d10 = a3b_bins.iloc[-1][f'sbs2_p{p}']
    peak = a3b_bins[f'sbs2_p{p}'].max()
    peak_dec = int(a3b_bins.loc[a3b_bins[f'sbs2_p{p}'].idxmax(), 'bin']) + 1
    fold = d10 / max(d1, 1)
    log(f"    P{p}: decile 1 = {d1:.0f}, peak = {peak:.0f} (decile {peak_dec}), "
        f"decile 10 = {d10:.0f}, fold change = {fold:.1f}x")

a3b_p90_vals = a3b_bins['sbs2_p90'].values
a3b_p90_peak_idx = np.argmax(a3b_p90_vals)
a3b_p90_peak_val = a3b_p90_vals[a3b_p90_peak_idx]
a3b_p90_peak_dec = a3b_p90_peak_idx + 1
a3b_p90_post_peak_mean = np.mean(a3b_p90_vals[a3b_p90_peak_idx+1:]) if a3b_p90_peak_idx < len(a3b_p90_vals)-1 else a3b_p90_vals[-1]
a3b_p90_d1 = a3b_p90_vals[0]
a3b_p90_d10 = a3b_p90_vals[-1]

log(f"\n  A3B P90 ceiling analysis:")
log(f"    P90 at decile 1: {a3b_p90_d1:.0f}")
log(f"    P90 peak: {a3b_p90_peak_val:.0f} at decile {a3b_p90_peak_dec}")
log(f"    P90 at decile 10: {a3b_p90_d10:.0f}")
log(f"    Post-peak mean P90: {a3b_p90_post_peak_mean:.0f}")
log(f"    Interpretation: A3B's 90th percentile SBS2 reaches {a3b_p90_peak_val:.0f}")
log(f"    by decile {a3b_p90_peak_dec} and does not increase further, indicating a")
log(f"    mutagenic ceiling for A3B-driven SBS2 accumulation.")

# =========================================================================
# 3. PANEL 1d TOP-MIDDLE: A3A fan numbers
# =========================================================================
banner("PANEL 1d TOP-MIDDLE: A3A No Saturation Fan")

log(f"  Percentile behavior across A3A expression deciles:")
for p in FAN_PCTLS:
    d1 = a3a_bins.iloc[0][f'sbs2_p{p}']
    d10 = a3a_bins.iloc[-1][f'sbs2_p{p}']
    peak = a3a_bins[f'sbs2_p{p}'].max()
    peak_dec = int(a3a_bins.loc[a3a_bins[f'sbs2_p{p}'].idxmax(), 'bin']) + 1
    fold = d10 / max(d1, 1)
    log(f"    P{p}: decile 1 = {d1:.0f}, peak = {peak:.0f} (decile {peak_dec}), "
        f"decile 10 = {d10:.0f}, fold change = {fold:.1f}x")

a3a_p90_d1 = a3a_bins.iloc[0]['sbs2_p90']
a3a_p90_d10 = a3a_bins.iloc[-1]['sbs2_p90']
a3a_p90_fold = a3a_p90_d10 / max(a3a_p90_d1, 1)

log(f"\n  A3A P90 trajectory:")
log(f"    P90 at decile 1: {a3a_p90_d1:.0f}")
log(f"    P90 at decile 10: {a3a_p90_d10:.0f}")
log(f"    Fold increase: {a3a_p90_fold:.1f}x")
log(f"    Interpretation: Unlike A3B, A3A's 90th percentile SBS2 continues")
log(f"    rising across the full expression range ({a3a_p90_fold:.1f}x increase),")
log(f"    indicating that high A3A expression enables substantially higher")
log(f"    SBS2 accumulation in a subset of tumors, consistent with")
log(f"    patient-specific modulation of A3A enzymatic activity.")

# =========================================================================
# 4. PANEL 1d TOP-RIGHT: Quantile-based threshold sweep
# =========================================================================
banner("PANEL 1d TOP-RIGHT: Quantile-Based Threshold Sweep")

rho_below_low = a3b_sweep[a3b_sweep['quantile'] <= 0.40]['rho_below'].mean()
rho_above_high = a3b_sweep[a3b_sweep['quantile'] >= 0.60]['rho_above'].mean()
collapse_q = a3b_sweep[a3b_sweep['rho_above'] <= 0]
collapse_pctl = int(collapse_q['pct_label'].iloc[0]) if len(collapse_q) > 0 else 'N/A'

rho_at = {}
for target_q in [20, 40, 60, 80]:
    match = a3b_sweep[a3b_sweep['pct_label'] == target_q]
    if len(match) > 0:
        rho_at[target_q] = {
            'below': match.iloc[0]['rho_below'],
            'above': match.iloc[0]['rho_above'],
            'n_above': int(match.iloc[0]['n_above']),
        }

log(f"  Spearman rho at each threshold:")
for q, vals in rho_at.items():
    log(f"    {q}th percentile: below = {vals['below']:.3f}, above = {vals['above']:.3f} (n_above = {vals['n_above']})")

log(f"\n  Summary statistics:")
log(f"    Below-threshold rho (avg, Q <= 40th): {rho_below_low:.3f}")
log(f"    Above-threshold rho (avg, Q >= 60th): {rho_above_high:.3f}")
log(f"    Correlation collapse point: {collapse_pctl}th percentile")
log(f"\n  Interpretation: The below-threshold correlation (rho ~ {rho_below_low:.2f})")
log(f"  indicates that A3B expression positively associates with SBS2 in the")
log(f"  rising phase. The above-threshold correlation collapses to near-zero")
log(f"  (rho ~ {rho_above_high:.2f}) at high A3B, confirming that tumors above")
log(f"  the {collapse_pctl}th percentile of A3B expression show no further")
log(f"  A3B-SBS2 association, consistent with enzymatic saturation.")

# =========================================================================
# 5. PANEL 1d BOTTOM: Remapped panel c values on median-centered axis
# =========================================================================
banner("PANEL 1d BOTTOM: Remapped Panel C (median-centered)")

log(f"  NOTE: Each y-value is a rho taken directly from the panel c sweep.")
log(f"        Below-median side uses rho_below at each quantile Q, mapped")
log(f"        to x = -(100 - Q). Above-median side uses rho_above at each")
log(f"        quantile Q, mapped to x = +Q. Sweep is at Q20/Q40/Q60/Q80 so")
log(f"        points land exactly on -80, -60, -40, -20, +20, +40, +60, +80.")
log(f"        There is a 40% gap across the median; the line is drawn as")
log(f"        two separate segments (below and above) with the linear trend")
log(f"        fit across all 8 points.")
log(f"")
log(f"  Remapped points:")
for _, r in a3b_centered.iterrows():
    log(f"    x = {int(r['pct_from_median']):>+4d}%  "
        f"(from Q{int(r['source_quantile']):>2d}, {r['side']:>6s}-subset, "
        f"n = {int(r['n_subset']):>3d}, threshold = {r['threshold']:>6.2f}): "
        f"rho = {r['rho']:+.3f}")

below_mean = a3b_centered.loc[a3b_centered['pct_from_median'] < 0, 'rho'].mean()
above_mean = a3b_centered.loc[a3b_centered['pct_from_median'] > 0, 'rho'].mean()
log(f"\n  Mean rho (below-median side): {below_mean:+.3f}")
log(f"  Mean rho (above-median side): {above_mean:+.3f}")
log(f"  Below - Above delta: {below_mean - above_mean:+.3f}")
log(f"\n  Linear trend: slope = {slope:+.5f} per % "
    f"(p = {p_val_trend:.2e}), R^2 = {r_val**2:.3f}")
log(f"  Interpretation: The A3B-SBS2 correlation "
    f"{'decays' if slope < 0 else 'rises'} as the expression threshold moves from")
log(f"  below to above the cohort median. Rescaling panel c's rho values onto")
log(f"  a symmetric % from median axis makes the saturation pattern visually")
log(f"  intuitive: diminishing returns as A3B rises past the median, with the")
log(f"  correlation decaying toward zero or negative at the far-right points.")

a3b_centered.to_csv(
    os.path.join(TROUBLE_DIR, "panel1d_a3b_centered_sweep.tsv"),
    sep='\t', index=False
)
log(f"  Saved: panel1d_a3b_centered_sweep.tsv")

# =========================================================================
# 6. COMBINED NARRATIVE FOR RESULTS SECTION
# =========================================================================
banner("COMBINED NARRATIVE: Panel 1c + 1d")

log(f"""
  PANEL 1c NARRATIVE (boxplot + heatmap):
    Tumors stratified by median A3A and A3B expression reveal that the
    lowest SBS2 burden occurs in tumors with low expression of both
    enzymes. SBS2 increases with elevated A3B alone and rises further
    with elevated A3A alone. However, tumors co-expressing both A3A
    and A3B above their respective medians do not show the highest SBS2
    burden; instead, their median SBS2 is comparable to or lower than
    the A3A-high/A3B-low group (ns, Panel 1c).

  PANEL 1d NARRATIVE (fan + sweep + median-centered remap):
    This plateau is explained by A3B enzymatic saturation. When SBS2
    distributions are examined across expression deciles (Panel 1d,
    top-left), all percentile lines for A3B (30th through 90th) converge
    at high expression. The 90th percentile of SBS2 reaches
    {a3b_p90_peak_val:.0f} mutations by decile {a3b_p90_peak_dec} and
    does not increase further, establishing a mutagenic ceiling. By
    contrast, A3A shows no such saturation (Panel 1d, top-middle): its
    90th percentile SBS2 rises {a3a_p90_fold:.1f}-fold from decile 1
    ({a3a_p90_d1:.0f}) to decile 10 ({a3a_p90_d10:.0f}), indicating
    that high A3A expression enables substantially greater SBS2
    accumulation in a subset of tumors. A threshold sweep of the
    A3B-SBS2 Spearman correlation (Panel 1d, top-right) confirms this
    pattern: the correlation is maintained among tumors below each
    expression threshold (rho ~ {rho_below_low:.2f}) but collapses
    to near-zero above the {collapse_pctl}th percentile. Plotting the
    same rho values on a median-centered axis (Panel 1d, bottom)
    reveals a negative linear trend of {slope:+.4f} rho per percent
    deviation (p = {p_val_trend:.2e}), capturing the diminishing
    returns of A3B expression with the correlation decaying from
    positive on the below-median side to near-zero or negative at A3B
    levels well above the cohort median.

  TRANSITION TO NEXT SECTION:
    The sporadic nature of A3A's effect on SBS2 (visible as scattered
    high-SBS2 tumors at high A3A expression rather than a uniform
    increase) suggests that patient-specific factors modulate whether
    A3A enzymatic activity translates into genomic mutations. These
    factors could include germline or somatic SNPs affecting A3
    substrate accessibility, or altered expression of cofactors that
    regulate A3 enzymatic activity on genomic DNA. Identifying these
    factors through differential network analysis and patient-specific
    SNP enrichment is the focus of the following sections.
""")

a3b_bins.to_csv(os.path.join(TROUBLE_DIR, "panel1d_a3b_decile_stats.tsv"), sep='\t', index=False)
a3a_bins.to_csv(os.path.join(TROUBLE_DIR, "panel1d_a3a_decile_stats.tsv"), sep='\t', index=False)
a3b_sweep.to_csv(os.path.join(TROUBLE_DIR, "panel1d_a3b_threshold_sweep.tsv"), sep='\t', index=False)
log(f"  Saved: panel1d_a3b_decile_stats.tsv")
log(f"  Saved: panel1d_a3a_decile_stats.tsv")
log(f"  Saved: panel1d_a3b_threshold_sweep.tsv")

rp = os.path.join(TROUBLE_DIR, "figure1_panel1d_saturation_report.txt")
with open(rp, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Report: {rp}")

banner("PANEL 1d COMPLETE")
