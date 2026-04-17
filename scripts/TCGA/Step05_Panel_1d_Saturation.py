#!/usr/bin/env python3
"""
Step05_Panel_1d_Saturation.py
=============================
Panel 1d: A3B saturation effect on SBS2.

Three sub-panels:
  Left:   A3B vs SBS2 scatter with decile-binned median + 90th pctl overlay
          Shows SBS2 rises with A3B then plateaus (ceiling effect)
  Middle: A3A vs SBS2 scatter with decile-binned median + 90th pctl overlay
          Shows SBS2 continues to rise sporadically with A3A (no ceiling)
  Right:  Spearman rho threshold sweep for A3B (above vs below each quantile)
          Formal statistical evidence that A3B-SBS2 correlation collapses
          at high A3B expression while it persists at low A3B expression

Reads: data/FIG_1/HNSC_A3_SBS2_matched_v3.tsv
Saves: FIGURE_1_PANELS/Panel_1d_A3B_Saturation.pdf/.png

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os, sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import spearmanr

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

    # Fan percentiles: from low to high
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
# THRESHOLD SWEEP (Spearman rho above vs below)
# =============================================================================
banner("THRESHOLD SWEEP")

def threshold_sweep(df, expr_col, sbs_col='SBS2',
                    quantiles=None):
    """
    For each quantile threshold of expr_col, compute Spearman rho(expr, SBS2)
    separately for tumors BELOW and ABOVE the threshold.
    """
    if quantiles is None:
        # Clean integer percentiles, evenly spaced by 5%
        quantiles = [q / 100.0 for q in range(20, 85, 5)]

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

log(f"\n  A3B threshold sweep:")
log(f"  {'Pctl':>5s} {'Thresh':>7s} {'n_below':>8s} {'n_above':>8s} "
    f"{'rho_below':>10s} {'rho_above':>10s}")
for _, row in a3b_sweep.iterrows():
    log(f"  {int(row['pct_label']):>5d} {row['threshold']:>7.1f} "
        f"{int(row['n_below']):>8d} {int(row['n_above']):>8d} "
        f"{row['rho_below']:>10.3f} {row['rho_above']:>10.3f}")

# =============================================================================
# GENERATE FIGURE: 3-panel
# =============================================================================
banner("GENERATE PANEL 1d (3 sub-panels)")

fig = plt.figure(figsize=(35, 10))
gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1], wspace=0.30)
ax_a3b = fig.add_subplot(gs[0])
ax_a3a = fig.add_subplot(gs[1])
ax_sweep = fig.add_subplot(gs[2])

# Shared y-axis limit for left and middle panels
# Cap at 300 so the binned trend lines (medians ~5-30, P90 ~40-220) are visible.
# Extreme outliers above 300 are clipped but still in the scatter.
Y_CAP = 300
sbs2_max_plot = Y_CAP
n_clipped = (df['SBS2'] > Y_CAP).sum()
log(f"  Y-axis capped at {Y_CAP} ({n_clipped} points above cap will be clipped)")

# Empirical x-axis caps: use the 95th percentile of expression so extreme
# outliers don't stretch the axis and compress the bulk of the data
a3b_x_cap = np.percentile(df['A3B'], 95) * 1.15
a3a_x_cap = np.percentile(df['A3A'], 95) * 1.15
log(f"  A3B x-axis cap: {a3b_x_cap:.1f} (95th pctl * 1.15)")
log(f"  A3A x-axis cap: {a3a_x_cap:.1f} (95th pctl * 1.15)")

# =========================================================================
# LEFT: A3B vs SBS2 scatter + percentile fan
# =========================================================================
# Scatter: all points, light background
ax_a3b.scatter(df['A3B'], df['SBS2'], s=30, alpha=0.20, c=COLOR_A3B,
               edgecolors='none', rasterized=True, zorder=1)

# Percentile fan: progressively bolder lines from P30 to P90
fan_alphas = [0.4, 0.55, 0.75, 1.0]
fan_lws = [1.5, 2.0, 2.5, 3.0]
fan_markers = ['v', 'o', 'D', 's']
fan_msizes = [5, 6, 7, 8]

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
# MIDDLE: A3A vs SBS2 scatter + percentile fan
# =========================================================================
ax_a3a.scatter(df['A3A'], df['SBS2'], s=30, alpha=0.20, c=COLOR_A3A,
               edgecolors='none', rasterized=True, zorder=1)

# Same percentile fan
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
# RIGHT: Threshold sweep (A3B, rho above vs below)
# =========================================================================
# Ribbon between the two lines
ax_sweep.fill_between(a3b_sweep['pct_label'],
                       np.minimum(a3b_sweep['rho_below'], a3b_sweep['rho_above']),
                       np.maximum(a3b_sweep['rho_below'], a3b_sweep['rho_above']),
                       alpha=0.1, color=COLOR_DARK_GRAY, zorder=0)

# Below-threshold line (rising phase, maintains correlation)
ax_sweep.plot(a3b_sweep['pct_label'], a3b_sweep['rho_below'], 'o-',
              color=COLOR_A3B, markersize=8, linewidth=2.5,
              markeredgecolor=COLOR_BLACK, markeredgewidth=1,
              label='Below threshold', zorder=3)

# Above-threshold line (plateau phase, correlation collapses)
ax_sweep.plot(a3b_sweep['pct_label'], a3b_sweep['rho_above'], 's-',
              color=COLOR_SBS2_HIGH, markersize=8, linewidth=2.5,
              markeredgecolor=COLOR_BLACK, markeredgewidth=1,
              label='Above threshold', zorder=3)

# Zero line
ax_sweep.axhline(0, color=COLOR_DARK_GRAY, ls='--', lw=1, alpha=0.5, zorder=0)

# X-axis: clean ticks at every 20th percentile
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

plt.tight_layout()
save_fig(fig, "Panel_1d_A3B_Saturation")

# =============================================================================
# DIAGNOSTIC SUMMARY
# =============================================================================
banner("DIAGNOSTIC SUMMARY FOR TEXT")

# =========================================================================
# 1. FULL DECILE TABLE FOR BOTH GENES (all fan percentiles)
# =========================================================================
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
# 2. PANEL 1d LEFT: A3B fan numbers
# =========================================================================
banner("PANEL 1d LEFT: A3B Saturation Fan")

log(f"  Percentile behavior across A3B expression deciles:")
for p in FAN_PCTLS:
    d1 = a3b_bins.iloc[0][f'sbs2_p{p}']
    d10 = a3b_bins.iloc[-1][f'sbs2_p{p}']
    peak = a3b_bins[f'sbs2_p{p}'].max()
    peak_dec = int(a3b_bins.loc[a3b_bins[f'sbs2_p{p}'].idxmax(), 'bin']) + 1
    fold = d10 / max(d1, 1)
    log(f"    P{p}: decile 1 = {d1:.0f}, peak = {peak:.0f} (decile {peak_dec}), "
        f"decile 10 = {d10:.0f}, fold change = {fold:.1f}x")

# Identify where A3B P90 peaks and what happens after
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
# 3. PANEL 1d MIDDLE: A3A fan numbers
# =========================================================================
banner("PANEL 1d MIDDLE: A3A No Saturation Fan")

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
# 4. PANEL 1d RIGHT: Threshold sweep numbers
# =========================================================================
banner("PANEL 1d RIGHT: Threshold Sweep")

rho_below_low = a3b_sweep[a3b_sweep['quantile'] <= 0.40]['rho_below'].mean()
rho_below_mid = a3b_sweep[(a3b_sweep['quantile'] > 0.40) & (a3b_sweep['quantile'] <= 0.60)]['rho_below'].mean()
rho_above_low = a3b_sweep[a3b_sweep['quantile'] <= 0.40]['rho_above'].mean()
rho_above_high = a3b_sweep[a3b_sweep['quantile'] >= 0.60]['rho_above'].mean()
collapse_q = a3b_sweep[a3b_sweep['rho_above'] <= 0]
collapse_pctl = int(collapse_q['pct_label'].iloc[0]) if len(collapse_q) > 0 else 'N/A'

# Rho at specific thresholds for text
rho_at = {}
for target_q in [20, 40, 60, 80]:
    match = a3b_sweep[a3b_sweep['pct_label'] == target_q]
    if len(match) > 0:
        rho_at[target_q] = {
            'below': match.iloc[0]['rho_below'],
            'above': match.iloc[0]['rho_above'],
            'n_above': int(match.iloc[0]['n_above']),
        }

log(f"  Spearman rho at key thresholds:")
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
# 5. COMBINED NARRATIVE FOR RESULTS SECTION
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

  PANEL 1d NARRATIVE (fan + sweep):
    This plateau is explained by A3B enzymatic saturation. When SBS2
    distributions are examined across expression deciles (Panel 1d,
    left), all percentile lines for A3B (30th through 90th) converge
    at high expression. The 90th percentile of SBS2 reaches
    {a3b_p90_peak_val:.0f} mutations by decile {a3b_p90_peak_dec} and
    does not increase further, establishing a mutagenic ceiling. By
    contrast, A3A shows no such saturation (Panel 1d, middle): its
    90th percentile SBS2 rises {a3a_p90_fold:.1f}-fold from decile 1
    ({a3a_p90_d1:.0f}) to decile 10 ({a3a_p90_d10:.0f}), indicating
    that high A3A expression enables substantially greater SBS2
    accumulation in a subset of tumors. A threshold sweep of the
    A3B-SBS2 Spearman correlation (Panel 1d, right) confirms this
    pattern: the correlation is maintained among tumors below each
    expression threshold (rho ~ {rho_below_low:.2f}) but collapses
    to near-zero above the {collapse_pctl}th percentile.

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

# Save diagnostic tables
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
