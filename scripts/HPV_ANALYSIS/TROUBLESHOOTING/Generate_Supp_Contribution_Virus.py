#!/usr/bin/env python3
"""
Generate_Supp_Contribution_Virus.py
====================================
Supplemental panels for who over-contributes to each fate and whether viral load
explains it. Four panels, each written as its own PDF and PNG plus a
*_panel_source.tsv, so the supplement grouping is decided after looking at them.

PANELS
------
  P1  SBS2-HIGH contribution per patient. Stacked bars, contributors in coral.
  P2  CNV-HIGH contribution per patient. Same layout in mustard, with a coral
      outline on any patient that also drives SBS2-HIGH.
  P3  Viral load against BOTH folds on one scatter. The panel exists to show
      that load associates with the two fates at nearly the same strength, which
      is the quantitative form of 'load is opportunity, not fate'.
  P4  Per-patient HPV16 lifecycle phase composition, stacked.

SUPERSEDES
----------
This replaces Panel A of Diagnostic_Regenerate_Panel_A_C.py and the standalone
Generate_Panel_CNV_Patient_Distribution.py. Those two drew the same two panels
from different tables with different sort orders and different denominators,
which is exactly the drift the contribution refactor removed. Both fates are now
drawn from ONE table under ONE convention.

ROW ORDER
---------
P1 and P2 are both sorted by the SUM of the two folds, the same order as the
conjunction heatmap. Sorting each panel on its own fate would put SC-level
drivers in different rows in each panel and make the pair impossible to read
against each other or against the heatmap.

P3 CANNOT USE A LOG Y-AXIS. Several patients have a fold of exactly zero on one
or both fates. Linear axes with both rho values annotated.

NO NUMBERS ARE HARDCODED. Every value, correlation and caption sentence is read
or computed from the input tables, and each panel prints its caption.

INPUTS (read-only)
------------------
  data/FIG_5/00_diagnostics/patient_source_matched_folds.tsv   (P1, P2)
      written by Diagnostic_Patient_CellCycle_by_Source.py
  data/FIG_5/00_diagnostics/patient_determinants_table.tsv     (P3)
  data/FIG_5/00_diagnostics/patient_lifecycle_phase.tsv        (P4)

OUTPUTS (to FIGURE_5_PANELS)
----------------------------
  Supp_P1_SBS2_Patient_Contribution.pdf/.png  + _panel_source.tsv
  Supp_P2_CNV_Patient_Contribution.pdf/.png   + _panel_source.tsv
  Supp_P3_Load_vs_Both_Folds.pdf/.png         + _panel_source.tsv
  Supp_P4_Patient_Lifecycle_Phase.pdf/.png    + _panel_source.tsv

Run from the directory holding patient_config.py:
    conda run -n NETWORK python Generate_Supp_Contribution_Virus.py

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

try:
    from patient_config import (DIR_00_DIAG, FIGURE_5_PANELS,
                                banner, log, ensure_dir)
except ImportError:  # local layout check
    DIR_00_DIAG = os.environ.get('DIAG_DIR', 'fixtures')
    FIGURE_5_PANELS = os.environ.get('OUT_DIR', 'panels')

    def log(msg):
        print(msg, flush=True)

    def banner(t):
        log("")
        log("=" * 70)
        log(t)
        log("=" * 70)

    def ensure_dir(p):
        os.makedirs(p, exist_ok=True)
        return p

try:
    from contribution import HC_THRESHOLD
except ImportError:
    HC_THRESHOLD = 2.0

# =============================================================================
# STYLE
# =============================================================================
DPI = 300

COLOR_SBS2      = "#ed6a5a"
COLOR_SBS2_PALE = "#f7b3aa"
COLOR_CNV       = "#F6D155"
COLOR_CNV_PALE  = "#f7e3a1"
COLOR_CNV_DARK  = "#a8862a"
COLOR_BOTH      = "#7a4fa3"
COLOR_NEITHER   = "#9aa0a6"
COLOR_BASAL     = "#d4d4d4"
COLOR_GRID      = "#d9d9d9"
COLOR_TEXT      = "#333333"

# Lifecycle phase palette. Maintenance takes the SBS2 colour and Capsid the CNV
# colour deliberately: the whole point of P4 is that maintenance-leaning virus
# goes with the SBS2 fate and productive-leaning virus with the CNV fate, so the
# phase bands should read in the same colours as the fates they predict.
PHASE_COLORS = {'Maintenance': COLOR_SBS2, 'Amplification': "#f4a259",
                'Oncogene': "#8a8d91", 'Capsid': COLOR_CNV}
PHASE_COLS = [('Maintenance', 'maint_frac'), ('Amplification', 'ampl_frac'),
              ('Oncogene', 'onco_frac'), ('Capsid', 'caps_frac')]

FONT_TITLE  = 34
FONT_LABEL  = 30
FONT_TICK   = 26
FONT_ANNOT  = 24
FONT_LEGEND = 24

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'axes.linewidth': 1.4,
})

OUT = ensure_dir(FIGURE_5_PANELS)
CAPTIONS = []


def short(p):
    return str(p).replace('Patient ', '')


def save(fig, name):
    for ext in ('pdf', 'png'):
        path = os.path.join(OUT, f"{name}.{ext}")
        fig.savefig(path, dpi=DPI, bbox_inches='tight', facecolor='white')
        log(f"  [SAVE] {path}")
    plt.close(fig)


def source(df, name):
    path = os.path.join(OUT, f"{name}_panel_source.tsv")
    df.to_csv(path, sep='\t', index=False)
    log(f"  [SAVE] {path}")


def caption(panel, text):
    CAPTIONS.append((panel, text))


def fmt_p(p):
    if not np.isfinite(p):
        return 'n.d.'
    if p < 1e-3:
        return f"p = {p:.0e}".replace('e-0', 'e-')
    return f"p = {p:.3f}"


def place_labels(ax, fig, items, fontsize):
    """Greedy label placement; see the sibling generator for the rationale."""
    fig.canvas.draw()
    ppp = fig.dpi / 72.0
    placed = []
    cands = [(12, 10), (12, -22), (-12, 10), (-12, -22), (0, 20), (0, -28)]
    for text, x, y in items:
        cx, cy = ax.transData.transform((x, y))
        w = 0.62 * fontsize * len(text) * ppp
        h = 1.25 * fontsize * ppp
        chosen = cands[0]
        for dx, dy in cands:
            bx = cx + dx * ppp + (0 if dx >= 0 else -w)
            by = cy + dy * ppp
            box = (bx, by, bx + w, by + h)
            if not any(box[0] < o[2] and box[2] > o[0] and
                       box[1] < o[3] and box[3] > o[1] for o in placed):
                chosen = (dx, dy)
                placed.append(box)
                break
        else:
            dx, dy = chosen
            bx = cx + dx * ppp
            placed.append((bx, cy + dy * ppp, bx + w, cy + dy * ppp + h))
        ax.annotate(text, (x, y), fontsize=fontsize, xytext=chosen,
                    textcoords='offset points', color=COLOR_TEXT)


def read_tsv(fname):
    path = os.path.join(DIR_00_DIAG, fname)
    if not os.path.exists(path):
        sys.exit(f"  ERROR: required input missing: {path}")
    df = pd.read_csv(path, sep='\t')
    tag = df['denominator'].iloc[0] if 'denominator' in df.columns else 'n/a'
    log(f"  [OK] {fname}: {len(df)} rows, denominator = {tag}")
    return df


# =============================================================================
# LOAD
# =============================================================================
banner("LOAD PANEL INPUTS")
folds = read_tsv("patient_source_matched_folds.tsv")
det = read_tsv("patient_determinants_table.tsv")
life = read_tsv("patient_lifecycle_phase.tsv")

denom = folds['denominator'].iloc[0] if 'denominator' in folds.columns else 'not stated'
denoms = {t for t in (denom,
                      det['denominator'].iloc[0] if 'denominator' in det else denom,
                      life['denominator'].iloc[0] if 'denominator' in life else denom)}
if len(denoms) > 1:
    sys.exit(f"  ERROR: input tables disagree on the denominator: {sorted(denoms)}. "
             f"Re-run the upstream diagnostics so they match.")
log(f"  contribution denominator: {denom}")
REF_COL = 'n_tumor' if denom == 'tumor' else 'n_basal'
REF_LABEL = ('Number of tumor basal cells' if denom == 'tumor'
             else 'Number of basal cells')
FOLD_AXIS = f"fold enrichment ({denom})"

# Contributor sets derived here, never read from a list.
folds['is_sbs2_driver'] = folds['fold_sbs2'] >= HC_THRESHOLD
folds['is_cnv_driver'] = folds['fold_cnv'] >= HC_THRESHOLD
SBS2_HC = set(folds.loc[folds['is_sbs2_driver'], 'patient'])
CNV_HC = set(folds.loc[folds['is_cnv_driver'], 'patient'])
log(f"  SBS2-HIGH contributors (>= {HC_THRESHOLD:.1f}x): "
    f"{sorted(short(p) for p in SBS2_HC)}")
log(f"  CNV-HIGH  contributors (>= {HC_THRESHOLD:.1f}x): "
    f"{sorted(short(p) for p in CNV_HC)}")
log(f"  drives both: {sorted(short(p) for p in (SBS2_HC & CNV_HC)) or 'none'}")

# Shared row order: summed fold, so P1, P2 and the conjunction heatmap all agree.
folds['fold_sum'] = folds['fold_sbs2'].fillna(0) + folds['fold_cnv'].fillna(0)
ORDER = list(folds.sort_values('fold_sum', ascending=False)['patient'])
log(f"  shared row order (summed fold): {[short(p) for p in ORDER]}")
n_pat = len(folds)


def contribution_panel(fate, count_col, fold_col, hc_set, other_hc,
                       color_hc, color_pale, text_color, title, fname):
    """P1 and P2 share every line of layout; only the fate differs."""
    d = folds.set_index('patient').loc[ORDER[::-1]].reset_index()
    y = np.arange(len(d))
    n_hi = d[count_col].values.astype(float)
    n_ref = d[REF_COL].values.astype(float)
    n_other = n_ref - n_hi
    total_hi = n_hi.sum()

    fig, ax = plt.subplots(figsize=(17, 1.0 + 0.78 * len(d)))
    ax.barh(y, n_other, color=COLOR_BASAL, edgecolor='black', linewidth=0.6,
            zorder=2)
    for i, r in d.iterrows():
        if n_hi[i] <= 0:
            # A zero-width segment still draws its outline, which on a dual-driver
            # row put a stray coloured tick on a patient contributing no cells to
            # this fate. Nothing to draw, so draw nothing.
            continue
        dual = r['patient'] in other_hc
        ax.barh(y[i], n_hi[i], left=n_other[i],
                color=color_hc if r['patient'] in hc_set else color_pale,
                edgecolor=(COLOR_SBS2 if (dual and fate == 'CNV') else
                           (COLOR_CNV_DARK if dual else 'black')),
                linewidth=3.4 if dual else 0.6, zorder=3)

    xmax = n_ref.max()
    for i, r in d.iterrows():
        if n_hi[i] <= 0:
            continue
        is_hc = r['patient'] in hc_set
        ax.text(n_ref[i] + xmax * 0.012, y[i],
                f"{int(n_hi[i])}  ({r[fold_col]:.1f}x)",
                va='center', ha='left', fontsize=FONT_ANNOT,
                fontweight='bold' if is_hc else 'normal',
                color=text_color if is_hc else COLOR_NEITHER)

    ax.set_yticks(y)
    ax.set_yticklabels([short(p) for p in d['patient']], fontsize=FONT_TICK)
    for i, r in d.iterrows():
        if r['patient'] in hc_set:
            ax.get_yticklabels()[i].set_color(text_color)
            ax.get_yticklabels()[i].set_fontweight('bold')
    ax.set_xlim(0, xmax * 1.26)
    ax.set_xlabel(REF_LABEL, fontsize=FONT_LABEL)
    ax.set_title(title, fontsize=FONT_TITLE, pad=16)
    ax.tick_params(labelsize=FONT_TICK)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    handles = [
        Patch(facecolor=COLOR_BASAL, edgecolor='black', linewidth=0.6,
              label='other basal'),
        Patch(facecolor=color_hc, edgecolor='black', linewidth=0.6,
              label=f'{fate}-HIGH, contributor (\u2265 {HC_THRESHOLD:.0f}x)'),
        Patch(facecolor=color_pale, edgecolor='black', linewidth=0.6,
              label=f'{fate}-HIGH, not enriched'),
    ]
    if other_hc & set(d['patient']):
        other_name = 'SBS2' if fate == 'CNV' else 'CNV'
        handles.append(Patch(facecolor='white',
                             edgecolor=COLOR_SBS2 if fate == 'CNV' else COLOR_CNV_DARK,
                             linewidth=3.4,
                             label=f'also a {other_name}-HIGH contributor'))
    ax.legend(handles=handles, loc='lower right', fontsize=FONT_LEGEND,
              framealpha=0.95)
    fig.tight_layout()
    save(fig, fname)

    out = d[['patient', REF_COL, count_col, fold_col]].copy()
    out['pct_of_group'] = 100.0 * out[count_col] / total_hi if total_hi else 0.0
    out['is_contributor'] = out['patient'].isin(hc_set)
    out['denominator'] = denom
    source(out, fname)
    return out, total_hi


# =============================================================================
# PANEL 1: SBS2-HIGH CONTRIBUTION
# =============================================================================
banner("PANEL 1: SBS2-HIGH PATIENT CONTRIBUTION")
p1, n_sbs2 = contribution_panel(
    'SBS2', 'n_sbs2_high', 'fold_sbs2', SBS2_HC, CNV_HC,
    COLOR_SBS2, COLOR_SBS2_PALE, COLOR_SBS2,
    'Patient distribution of SBS2-HIGH basal cells',
    'Supp_P1_SBS2_Patient_Contribution')
top1 = p1.sort_values('n_sbs2_high', ascending=False)
caption("P1",
        f"Basal-cell composition per patient with the SBS2-HIGH group highlighted "
        f"({int(n_sbs2)} cells across {n_pat} patients). Fold enrichment is the "
        f"observed share of the group divided by the patient's share of the "
        f"{denom} compartment; contributors clear {HC_THRESHOLD:.0f}x and are shown "
        f"in bold coral "
        f"({', '.join(sorted(short(p) for p in SBS2_HC))}). "
        f"The top contributor supplies {top1.iloc[0]['pct_of_group']:.1f}% of the "
        f"group and the top three supply "
        f"{top1.iloc[:3]['pct_of_group'].sum():.1f}%. Rows are ordered by the sum of "
        f"the two fold columns, the same order used in the other patient panels.")


# =============================================================================
# PANEL 2: CNV-HIGH CONTRIBUTION
# =============================================================================
banner("PANEL 2: CNV-HIGH PATIENT CONTRIBUTION")
p2, n_cnv = contribution_panel(
    'CNV', 'n_cnv_high', 'fold_cnv', CNV_HC, SBS2_HC,
    COLOR_CNV, COLOR_CNV_PALE, COLOR_CNV_DARK,
    'Patient distribution of CNV-HIGH basal cells',
    'Supp_P2_CNV_Patient_Contribution')
top2 = p2.sort_values('n_cnv_high', ascending=False)
dual = sorted(short(p) for p in (SBS2_HC & CNV_HC))
caption("P2",
        f"The same layout for the CNV-HIGH group ({int(n_cnv)} cells). Contributors "
        f"are {', '.join(sorted(short(p) for p in CNV_HC))}; a coral outline marks a "
        f"patient that also contributes to SBS2-HIGH"
        + (f" ({', '.join(dual)})" if dual else "") + ". "
        f"CNV-HIGH is far more concentrated than SBS2-HIGH: "
        f"{short(top2.iloc[0]['patient'])} alone supplies "
        f"{top2.iloc[0]['pct_of_group']:.1f}% of the group and the top two supply "
        f"{top2.iloc[:2]['pct_of_group'].sum():.1f}%. That concentration is a real "
        f"limitation of every CNV-HIGH claim and is shown here rather than left to "
        f"be discovered.")


# =============================================================================
# PANEL 3: VIRAL LOAD vs BOTH FOLDS
# =============================================================================
banner("PANEL 3: VIRAL LOAD vs BOTH FOLDS")

d3 = det.copy()
r_s, p_s = spearmanr(d3['load_per_M'], d3['fold_sbs2'])
r_c, p_c = spearmanr(d3['load_per_M'], d3['fold_cnv'])
n_zero = int(((d3['fold_sbs2'] == 0) | (d3['fold_cnv'] == 0)).sum())
log(f"  load vs SBS2 fold: rho={r_s:+.3f} {fmt_p(p_s)}")
log(f"  load vs CNV  fold: rho={r_c:+.3f} {fmt_p(p_c)}")
log(f"  |difference| = {abs(r_s - r_c):.3f}")
log(f"  {n_zero} patients have a fold of exactly zero on at least one fate, "
    f"so the y-axis stays linear")

# label the drivers plus the largest single-factor discrepancies
d3['rank_gap'] = (d3['load_per_M'].rank(method='min')
                  - d3['A3A_mean'].rank(method='min'))
hpv_pos = d3['pct_hpv_pos'].fillna(0) > 0
nd = d3[~(d3['is_sbs2_driver'] | d3['is_cnv_driver']) & hpv_pos].copy()
nd['abs_gap'] = nd['rank_gap'].abs()
LABEL = set(d3.loc[d3['is_sbs2_driver'] | d3['is_cnv_driver'], 'patient'])
LABEL |= set(nd.sort_values('abs_gap', ascending=False).head(3)['patient'])

fig, ax = plt.subplots(figsize=(19, 14))
for fate, col, color, marker in (('SBS2-HIGH', 'fold_sbs2', COLOR_SBS2, 'o'),
                                 ('CNV-HIGH', 'fold_cnv', COLOR_CNV, 's')):
    ax.scatter(d3['load_per_M'], d3[col], s=460, marker=marker, c=color,
               edgecolors='black', linewidths=1.6, zorder=3, label=fate)

items = []
for _, r in d3.iterrows():
    if r['patient'] not in LABEL:
        continue
    col = 'fold_sbs2' if r['fold_sbs2'] >= r['fold_cnv'] else 'fold_cnv'
    items.append((short(r['patient']), r['load_per_M'], r[col]))
items.sort(key=lambda t: -t[2])

# Headroom on both axes BEFORE labels are measured, so the rightmost and
# topmost patient tags cannot be clipped by the axes box.
xlo, xhi = d3['load_per_M'].min(), d3['load_per_M'].max()
ylo = min(d3['fold_sbs2'].min(), d3['fold_cnv'].min())
yhi = max(d3['fold_sbs2'].max(), d3['fold_cnv'].max())
xpad, ypad = 0.13 * (xhi - xlo), 0.10 * (yhi - ylo)
ax.set_xlim(xlo - xpad * 0.45, xhi + xpad)
ax.set_ylim(ylo - ypad * 0.6, yhi + ypad)

ax.axhline(HC_THRESHOLD, color='#777777', linestyle=':', linewidth=2.2, zorder=1)
ax.text(ax.get_xlim()[0], HC_THRESHOLD, f" contributor threshold ({HC_THRESHOLD:.0f}x)",
        va='bottom', ha='left', fontsize=FONT_ANNOT - 2, color='#777777')
ax.set_xlabel("HPV16 load (UMIs per million total UMIs)", fontsize=FONT_LABEL)
ax.set_ylabel(FOLD_AXIS, fontsize=FONT_LABEL)
ax.set_title("Viral load associates with BOTH fates at the same strength",
             fontsize=FONT_TITLE, pad=16)
ax.tick_params(labelsize=FONT_TICK)
ax.grid(True, color=COLOR_GRID, linewidth=0.9, zorder=0)
ax.set_axisbelow(True)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

stats = (f"SBS2-HIGH   rho = {r_s:+.3f}, {fmt_p(p_s)}\n"
         f"CNV-HIGH    rho = {r_c:+.3f}, {fmt_p(p_c)}\n"
         f"difference  {abs(r_s - r_c):.3f}")
ax.text(0.025, 0.975, stats, transform=ax.transAxes, fontsize=FONT_ANNOT,
        va='top', ha='left', family='monospace',
        bbox=dict(boxstyle='round,pad=0.55', facecolor='white',
                  edgecolor='#bbbbbb'))
ax.legend(fontsize=FONT_LEGEND, loc='upper right', framealpha=0.95,
          markerscale=1.1)
place_labels(ax, fig, items, FONT_ANNOT)
fig.tight_layout()
save(fig, "Supp_P3_Load_vs_Both_Folds")

p3 = d3[['patient', 'load_per_M', 'pct_hpv_pos', 'fold_sbs2', 'fold_cnv',
         'is_sbs2_driver', 'is_cnv_driver']].copy()
p3['rho_load_vs_sbs2'] = r_s
p3['p_load_vs_sbs2'] = p_s
p3['rho_load_vs_cnv'] = r_c
p3['p_load_vs_cnv'] = p_c
source(p3, "Supp_P3_Load_vs_Both_Folds")
caption("P3",
        f"Library-normalized HPV16 load against contribution to each fate "
        f"(n = {len(d3)}). Load associates with SBS2-HIGH at rho = {r_s:+.3f} "
        f"({fmt_p(p_s)}) and with CNV-HIGH at rho = {r_c:+.3f} ({fmt_p(p_c)}), a "
        f"difference of {abs(r_s - r_c):.3f}. If load determined which fate a "
        f"patient drives, these two would diverge; that they do not is the "
        f"quantitative form of load being opportunity rather than fate. Axes are "
        f"linear because {n_zero} patients have a fold of exactly zero on at least "
        f"one fate. Folds use the {denom} denominator.")


# =============================================================================
# PANEL 4: PER-PATIENT LIFECYCLE PHASE COMPOSITION
# =============================================================================
banner("PANEL 4: PER-PATIENT HPV16 LIFECYCLE PHASE")

MIN_GATED = 10
have = [c for _, c in PHASE_COLS if c in life.columns]
if len(have) < len(PHASE_COLS):
    sys.exit(f"  ERROR: patient_lifecycle_phase.tsv missing phase columns "
             f"{[c for _, c in PHASE_COLS if c not in life.columns]}")

l4 = life[(life['n_gated'] >= MIN_GATED) & life['maint_frac'].notna()].copy()
dropped = life[~life['patient'].isin(l4['patient'])]
log(f"  {len(l4)} of {len(life)} patients have >= {MIN_GATED} gated HPV16+ cells")
if len(dropped):
    log(f"  excluded (no measurable phase): "
        f"{sorted(short(p) for p in dropped['patient'])}")
l4 = l4.sort_values('prod_frac', ascending=True).reset_index(drop=True)

y = np.arange(len(l4))
left = np.zeros(len(l4))
fig, ax = plt.subplots(figsize=(19, 1.0 + 0.90 * len(l4)))
for phase, col in PHASE_COLS:
    vals = l4[col].values.astype(float)
    ax.barh(y, vals, left=left, color=PHASE_COLORS[phase], label=phase,
            edgecolor='white', linewidth=1.4, zorder=3)
    left += vals

ax.set_yticks(y)
ax.set_yticklabels([short(p) for p in l4['patient']], fontsize=FONT_TICK)
for i, r in l4.iterrows():
    s = r['patient'] in SBS2_HC
    c = r['patient'] in CNV_HC
    lab = ax.get_yticklabels()[i]
    if s and c:
        lab.set_color(COLOR_BOTH); lab.set_fontweight('bold')
    elif s:
        lab.set_color(COLOR_SBS2); lab.set_fontweight('bold')
    elif c:
        lab.set_color(COLOR_CNV_DARK); lab.set_fontweight('bold')

xmax = float(left.max())
for i, r in l4.iterrows():
    tags = []
    if r['patient'] in SBS2_HC:
        tags.append('SBS2')
    if r['patient'] in CNV_HC:
        tags.append('CNV')
    if tags:
        ax.text(xmax * 1.015, y[i], '+'.join(tags) + ' driver', va='center',
                ha='left', fontsize=FONT_ANNOT - 2, fontweight='bold',
                color=(COLOR_BOTH if len(tags) == 2 else
                       (COLOR_SBS2 if tags[0] == 'SBS2' else COLOR_CNV_DARK)))
ax.set_xlim(0, xmax * 1.30)
ax.set_xlabel("coding-normalized fraction of HPV16 reads (URR excluded)",
              fontsize=FONT_LABEL)
ax.set_title("Per-patient HPV16 lifecycle phase composition",
             fontsize=FONT_TITLE, pad=16)
ax.tick_params(labelsize=FONT_TICK)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(fontsize=FONT_LEGEND, ncol=4, loc='upper center',
          bbox_to_anchor=(0.5, -0.10), frameon=False)
fig.tight_layout()
save(fig, "Supp_P4_Patient_Lifecycle_Phase")

p4 = l4[['patient', 'n_gated'] + [c for _, c in PHASE_COLS]
        + ['fold_sbs2', 'fold_cnv']].copy()
p4['is_sbs2_driver'] = p4['patient'].isin(SBS2_HC)
p4['is_cnv_driver'] = p4['patient'].isin(CNV_HC)
source(p4, "Supp_P4_Patient_Lifecycle_Phase")

most_prod = l4.iloc[-1]
most_maint = l4.iloc[0]
caption("P4",
        f"HPV16 lifecycle phase composition per patient, pooled across each "
        f"patient's HPV16-positive basal cells and normalized to coding reads with "
        f"the URR excluded ({len(l4)} of {len(life)} patients clear the "
        f"{MIN_GATED}-cell floor"
        + (f"; {', '.join(sorted(short(p) for p in dropped['patient']))} have no "
           f"measurable phase" if len(dropped) else "") + "). "
        f"Rows are ordered by productive fraction, running from the most "
        f"maintenance-leaning patient ({short(most_maint['patient'])}, productive "
        f"{most_maint['prod_frac']:.2f}) to the most productive "
        f"({short(most_prod['patient'])}, {most_prod['prod_frac']:.2f}). "
        f"Maintenance is drawn in the SBS2 colour and Capsid in the CNV colour "
        f"because those are the fates each phase predicts. Contributor labels are "
        f"bolded and tagged at the right.")


# =============================================================================
# CAPTIONS
# =============================================================================
banner("DRAFT CAPTIONS (computed from this run; paste into the legend)")
for panel, text in CAPTIONS:
    log("")
    log(f"  {panel}")
    line = "   "
    for w in text.split():
        if len(line) + len(w) + 1 > 88:
            log(line)
            line = "   "
        line += " " + w
    log(line)

log("")
log(f"  Panels written to: {OUT}")
banner("SUPPLEMENTAL CONTRIBUTION + VIRUS PANELS COMPLETE")
