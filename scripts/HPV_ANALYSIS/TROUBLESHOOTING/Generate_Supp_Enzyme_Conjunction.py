#!/usr/bin/env python3
"""
Generate_Supp_Enzyme_Conjunction.py  (v3)
==========================================
Supplemental panels for the patient-level determinants of fate. Four panels, each
written as its own PDF and PNG plus a *_panel_source.tsv holding the exact values
behind it.

v3 CHANGES (Panel 7 only; P5, P6 and P8 are unchanged)
------------------------------------------------------
The conjunction heatmap now reads in units a reader can check against the stated
rule without consulting Methods.

  - OPPORTUNITY is shown as the percentage of the patient's basal-compartment UMIs
    that are HPV16, rather than as parts per million. 0.35% is legible where 3454
    is not, and a reader can see at a glance which patients fall under 0.1%.
  - DIRECTION is shown as a percentage of viral CODING reads (URR excluded), so
    columns 1 and 2 share a unit.
  - CAPABILITY now carries two numbers: the mean log-normalized expression that the
    threshold acts on, and underneath it the penetrance, written as "1 in N cells".
    The per-patient enzyme mean decomposes as
        mean over all cells = fraction expressing x mean among expressing cells
    and the second term is close to constant across patients, so the mean is
    largely a rescaling of how many cells express the enzyme. Showing both lets a
    reader see that a value clearing the A3A cut may still fail the higher A3B cut.
  - Every column header states the threshold it is being tested against, read from
    the source table rather than hardcoded, so the panel documents its own rule.

All four quantities are derived from columns already written by
Diagnostic_Patient_Determinants_Table.py; the percentage forms are computed here.
If that script is later patched to write them directly, the columns are picked up
in preference to the derived ones.

PANELS
------
  P5  Enzyme double dissociation, 2x2 of scatters.
  P6  A3A / A3B, tumor vs normal-adjacent. CONTROL only.
  P7  Three-axis conjunction, continuous, both fates side by side.
  P8  Threshold sensitivity, 5x5 per fate.

NO NUMBERS ARE HARDCODED. Every value, threshold, correlation and caption sentence
is read or computed from the input tables.

INPUTS (read-only, written by Diagnostic_Patient_Determinants_Table.py)
-----------------------------------------------------------------------
  data/FIG_5/00_diagnostics/patient_determinants_table.tsv
  data/FIG_5/00_diagnostics/patient_conjunction_model.tsv
  data/FIG_5/00_diagnostics/patient_conjunction_sensitivity.tsv
  data/FIG_5/00_diagnostics/patient_enzyme_by_source.tsv

Run from the directory holding patient_config.py:
    conda run -n NETWORK python Generate_Supp_Enzyme_Conjunction.py

Author: Jake Lehle
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D

try:
    from patient_config import DIR_00_DIAG, FIGURE_5_PANELS, banner, log, ensure_dir
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

# =============================================================================
# STYLE
# =============================================================================
DPI = 300

COLOR_SBS2      = "#ed6a5a"
COLOR_CNV       = "#F6D155"
COLOR_CNV_DARK  = "#a8862a"
COLOR_BOTH      = "#7a4fa3"
COLOR_NEITHER   = "#9aa0a6"
COLOR_NULL      = "#8c8c8c"
COLOR_GRID      = "#d9d9d9"
COLOR_TEXT      = "#333333"

FONT_TITLE  = 34
FONT_LABEL  = 30
FONT_TICK   = 26
FONT_ANNOT  = 24
FONT_LEGEND = 24
FONT_CELL   = 22    # primary value inside a heatmap cell
FONT_SUB    = 17    # the penetrance line under a capability value

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


def one_in(pct):
    """
    Penetrance written the way it should be read. Below 4% the reciprocal gets
    unwieldy and misleadingly precise, so it is capped rather than printed.
    """
    if not np.isfinite(pct) or pct <= 0:
        return "none"
    if pct < 4:
        return "1 in >25"
    return f"1 in {int(round(100.0 / pct))}"


def role_of(row):
    s, c = bool(row['is_sbs2_driver']), bool(row['is_cnv_driver'])
    if s and c:
        return 'both'
    if s:
        return 'sbs2'
    if c:
        return 'cnv'
    return 'none'


ROLE_COLOR = {'both': COLOR_BOTH, 'sbs2': COLOR_SBS2,
              'cnv': COLOR_CNV, 'none': COLOR_NEITHER}
ROLE_LABEL = {'both': 'drives both', 'sbs2': 'SBS2-HIGH driver',
              'cnv': 'CNV-HIGH driver', 'none': 'drives neither'}


def place_labels(ax, fig, items, fontsize):
    """Greedy label placement for a small scatter."""
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


def read_tsv(fname, required=True):
    path = os.path.join(DIR_00_DIAG, fname)
    if not os.path.exists(path):
        if required:
            sys.exit(f"  ERROR: required input missing: {path}\n"
                     f"  Run Diagnostic_Patient_Determinants_Table.py first.")
        log(f"  [SKIP] optional input missing: {path}")
        return None
    df = pd.read_csv(path, sep='\t')
    log(f"  [OK] {fname}: {len(df)} rows")
    return df


# =============================================================================
# LOAD
# =============================================================================
banner("LOAD PANEL INPUTS")
det = read_tsv("patient_determinants_table.tsv")
conj = read_tsv("patient_conjunction_model.tsv")
sweep = read_tsv("patient_conjunction_sensitivity.tsv")
enz = read_tsv("patient_enzyme_by_source.tsv", required=False)

denom = det['denominator'].iloc[0] if 'denominator' in det.columns else 'not stated'
log(f"  contribution denominator in these tables: {denom}")
FOLD_AXIS = f"fold enrichment ({denom})"

det['role'] = det.apply(role_of, axis=1)
n_pat = len(det)


# =============================================================================
# PANEL 5: ENZYME DOUBLE DISSOCIATION
# =============================================================================
banner("PANEL 5: ENZYME DOUBLE DISSOCIATION")

det['rank_load'] = det['load_per_M'].rank(method='min')
det['rank_a3a'] = det['A3A_mean'].rank(method='min')
det['rank_gap'] = det['rank_load'] - det['rank_a3a']
hpv_pos = det['pct_hpv_pos'].fillna(0) > 0
nd = det[~(det['is_sbs2_driver'] | det['is_cnv_driver']) & hpv_pos].copy()
nd['abs_gap'] = nd['rank_gap'].abs()
DRIVERS = set(det.loc[det['is_sbs2_driver'] | det['is_cnv_driver'], 'patient'])
LABEL = DRIVERS | set(nd.sort_values('abs_gap', ascending=False).head(3)['patient'])
log(f"  drivers: {sorted(short(p) for p in DRIVERS)}")
log(f"  labelled on prediction panels: {sorted(short(p) for p in LABEL)}")

CELLS = [
    ('A3A_mean', 'fold_sbs2', 'A3A', 'SBS2-HIGH', True,  COLOR_SBS2),
    ('A3A_mean', 'fold_cnv',  'A3A', 'CNV-HIGH',  False, COLOR_NULL),
    ('A3B_mean', 'fold_sbs2', 'A3B', 'SBS2-HIGH', False, COLOR_NULL),
    ('A3B_mean', 'fold_cnv',  'A3B', 'CNV-HIGH',  True,  COLOR_CNV_DARK),
]

fig, axes = plt.subplots(2, 2, figsize=(24, 21))
p5_rows = []
for ax, (xc, yc, ename, fname, is_pred, accent) in zip(axes.flat, CELLS):
    sub = det[[xc, yc]].replace([np.inf, -np.inf], np.nan).dropna()
    rho, pv = spearmanr(sub[xc], sub[yc])
    p5_rows.append({'enzyme': ename, 'fate': fname, 'predicted': is_pred,
                    'rho': rho, 'p': pv, 'n': len(sub)})

    for _, r in det.iterrows():
        ax.scatter(r[xc], r[yc], s=460, zorder=3,
                   c=ROLE_COLOR[r['role']], edgecolors='black', linewidths=1.6)

    labels_here = LABEL if is_pred else DRIVERS
    items = [(short(r['patient']), r[xc], r[yc]) for _, r in det.iterrows()
             if r['patient'] in labels_here and np.isfinite(r[xc]) and np.isfinite(r[yc])]
    items.sort(key=lambda t: -t[2])
    place_labels(ax, fig, items, FONT_ANNOT)

    star = ' *' if pv < 0.05 else ''
    ax.set_title(f"{ename} vs {fname}\nrho = {rho:+.3f}, {fmt_p(pv)}{star}",
                 fontsize=FONT_TITLE - 4, fontweight='bold' if is_pred else 'normal',
                 color=accent if is_pred else COLOR_TEXT, pad=16)
    ax.set_xlabel(f"{ename} mean expression", fontsize=FONT_LABEL - 4)
    ax.set_ylabel(f"{fname} {FOLD_AXIS}", fontsize=FONT_LABEL - 4)
    ax.tick_params(labelsize=FONT_TICK)
    ax.grid(True, color=COLOR_GRID, linewidth=0.9, zorder=0)
    ax.set_axisbelow(True)
    for sp_ in ax.spines.values():
        sp_.set_color(accent if is_pred else '#999999')
        sp_.set_linewidth(4.0 if is_pred else 1.4)

handles = [Line2D([], [], marker='o', linestyle='', markersize=20,
                  markerfacecolor=ROLE_COLOR[k], markeredgecolor='black',
                  label=ROLE_LABEL[k]) for k in ('sbs2', 'cnv', 'both', 'none')]
fig.legend(handles=handles, loc='lower center', ncol=4, fontsize=FONT_LEGEND,
           frameon=False, bbox_to_anchor=(0.5, -0.015))
fig.tight_layout(rect=[0, 0.035, 1, 1])
save(fig, "Supp_P5_Enzyme_Double_Dissociation")

p5 = pd.DataFrame(p5_rows)
source(p5, "Supp_P5_Enzyme_Double_Dissociation")
d_sbs2 = p5[(p5.enzyme == 'A3A') & (p5.fate == 'SBS2-HIGH')].iloc[0]
d_cnv = p5[(p5.enzyme == 'A3B') & (p5.fate == 'CNV-HIGH')].iloc[0]
x_sbs2 = p5[(p5.enzyme == 'A3B') & (p5.fate == 'SBS2-HIGH')].iloc[0]
x_cnv = p5[(p5.enzyme == 'A3A') & (p5.fate == 'CNV-HIGH')].iloc[0]
caption("P5",
        f"Per-patient A3A and A3B mean expression against contribution to each "
        f"fate (n = {n_pat}). A3A tracks the SBS2-HIGH fate (rho = {d_sbs2.rho:+.3f}, "
        f"{fmt_p(d_sbs2.p)}) and A3B the CNV-HIGH fate (rho = {d_cnv.rho:+.3f}, "
        f"{fmt_p(d_cnv.p)}), while both cross-terms are null "
        f"(A3A vs CNV rho = {x_cnv.rho:+.3f}, {fmt_p(x_cnv.p)}; "
        f"A3B vs SBS2 rho = {x_sbs2.rho:+.3f}, {fmt_p(x_sbs2.p)}). Bold frames mark "
        f"the two pre-specified predictions; the off-diagonal panels are controls. "
        f"Folds use the {denom} denominator.")


# =============================================================================
# PANEL 6: A3A / A3B TUMOR vs NORMAL-ADJACENT (control)
# =============================================================================
if enz is not None and len(enz):
    banner("PANEL 6: ENZYME, TUMOR vs NORMAL-ADJACENT (control)")
    MIN_FOLD = 2.0
    e = enz.copy()
    e['n_pass'] = (((e['A3A_fold'] >= MIN_FOLD) & (e['A3A_p'] < 0.05)).astype(int)
                   + ((e['A3B_fold'] >= MIN_FOLD) & (e['A3B_p'] < 0.05)).astype(int))
    e['max_fold'] = e[['A3A_fold', 'A3B_fold']].max(axis=1)
    e = e.sort_values(['n_pass', 'max_fold'],
                      ascending=[False, False]).reset_index(drop=True)
    log("  row order (enzymes passing the gate, then larger fold): "
        + ", ".join(f"{short(r['patient'])} [{int(r['n_pass'])}]"
                    for _, r in e.iterrows()))
    y = np.arange(len(e))
    h = 0.36

    fig, axes = plt.subplots(1, 2, figsize=(24, 4.2 + 2.6 * len(e)), sharey=True)
    for ax, gene, accent in ((axes[0], 'A3A', COLOR_SBS2),
                             (axes[1], 'A3B', COLOR_CNV_DARK)):
        ax.barh(y - h / 2, e[f'{gene}_tumor'], height=h, color=accent,
                edgecolor='black', linewidth=1.2, label='tumor')
        ax.barh(y + h / 2, e[f'{gene}_normal'], height=h, color='white',
                edgecolor=accent, linewidth=2.6, hatch='///',
                label='normal-adjacent')
        xmax = max(e[f'{gene}_tumor'].max(), e[f'{gene}_normal'].max())
        for i, r in e.iterrows():
            passes = (r[f'{gene}_fold'] >= MIN_FOLD) and (r[f'{gene}_p'] < 0.05)
            ax.text(xmax * 1.03, y[i],
                    f"{r[f'{gene}_fold']:.2f}x" + ("  *" if passes else ""),
                    va='center', ha='left', fontsize=FONT_ANNOT,
                    fontweight='bold' if passes else 'normal',
                    color=accent if passes else COLOR_NEITHER)
        ax.set_xlim(0, xmax * 1.30)
        ax.set_xlabel(f"{gene} mean expression", fontsize=FONT_LABEL - 4)
        ax.set_title(gene, fontsize=FONT_TITLE, fontweight='bold',
                     color=accent, pad=14)
        ax.tick_params(labelsize=FONT_TICK)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([short(p) for p in e['patient']], fontsize=FONT_TICK)
    axes[0].invert_yaxis()
    axes[0].legend(fontsize=FONT_LEGEND, loc='lower right', frameon=False)
    fig.text(0.5, -0.02,
             f"*  fold \u2265 {MIN_FOLD:.1f} and p < 0.05. Significance alone is not "
             f"enough at these sample sizes.",
             ha='center', fontsize=FONT_LEGEND, color=COLOR_TEXT)
    fig.tight_layout()
    save(fig, "Supp_P6_Enzyme_Tumor_vs_Normal")
    source(e, "Supp_P6_Enzyme_Tumor_vs_Normal")

    passed = e[(e['A3A_fold'] >= MIN_FOLD) & (e['A3A_p'] < 0.05)]
    caption("P6",
            f"A3A and A3B mean expression in tumor versus normal-adjacent basal "
            f"cells, for the {len(e)} patients contributing both tissue sources. "
            f"A3A clears both a fold change of {MIN_FOLD:.1f} and p < 0.05 in "
            f"{len(passed)} of {len(e)} patients "
            f"({', '.join(sorted(short(p) for p in passed['patient']))}). "
            f"Elsewhere normal-adjacent tissue already carries A3A at close to "
            f"the tumor level, so A3A capability appears partly constitutive at "
            f"the patient level rather than purely tumor-induced. Small n: this "
            f"is a control, not a result.")


# =============================================================================
# PANEL 7: CONJUNCTION HEATMAP
# =============================================================================
banner("PANEL 7: THREE-AXIS CONJUNCTION HEATMAP")

cj = conj.copy()

# ---- derive the reader-facing units -----------------------------------------
# HPV16 UMIs as a PERCENTAGE of the patient's basal-compartment UMIs. load_per_M
# is the same quantity per million, so this is a rescaling and not a new measure.
if 'pct_umi_hpv16' not in cj.columns:
    if 'load_per_M' not in cj.columns:
        sys.exit("  ERROR: patient_conjunction_model.tsv has neither "
                 "'pct_umi_hpv16' nor 'load_per_M'.")
    cj['pct_umi_hpv16'] = cj['load_per_M'] / 1e4
    log("  derived pct_umi_hpv16 from load_per_M (per million / 1e4)")

# Penetrance, so a capability value can be read as a fraction of cells. These
# live in the determinants table rather than the conjunction table.
for col in ('A3A_pct_pos', 'A3B_pct_pos'):
    if col not in cj.columns:
        if col in det.columns:
            cj = cj.merge(det[['patient', col]], on='patient', how='left')
        else:
            cj[col] = np.nan
            log(f"  [WARN] {col} unavailable; capability cells will show the mean only")

# ---- thresholds, read from the table so the panel documents its own rule -----
thr_a = cj['thr_a3a'].iloc[0] if 'thr_a3a' in cj.columns else np.nan
thr_b = cj['thr_a3b'].iloc[0] if 'thr_a3b' in cj.columns else np.nan
if 'thr_opportunity_pct_umi' in cj.columns:
    thr_opp = cj['thr_opportunity_pct_umi'].iloc[0]
elif 'thr_load' in cj.columns:
    thr_opp = cj['thr_load'].iloc[0] / 1e4
else:
    thr_opp = np.nan
if 'thr_direction' in cj.columns:
    thr_dir = cj['thr_direction'].iloc[0]
elif 'direction_min_frac' in cj.columns:
    thr_dir = cj['direction_min_frac'].iloc[0]
else:
    thr_dir = np.nan

dir_rule = cj['direction_rule'].iloc[0] if 'direction_rule' in cj.columns else 'compare'
DIR_HDR = (f"\u2265 {100*thr_dir:.0f}% of coding reads"
           if dir_rule == 'threshold' and np.isfinite(thr_dir)
           else "dominant stage")

log(f"  thresholds read from the table:")
log(f"    opportunity  HPV16 \u2265 {thr_opp:.2f}% of basal-compartment UMIs")
log(f"    direction    {DIR_HDR}")
log(f"    capability   A3A > {thr_a:.2f}, A3B > {thr_b:.2f}")

# ---- penetrance decomposition, reported for the text ------------------------
log("")
log("  Enzyme expression decomposes as: mean over all cells = fraction expressing")
log("  x mean among expressing cells. If the second term is near-constant, the")
log("  per-patient mean is largely a rescaling of how many cells express it.")
for g, thr in (('A3A', thr_a), ('A3B', thr_b)):
    pc = cj[f'{g}_pct_pos'] / 100.0
    ratio = (cj[f'{g}_mean'] / pc).replace([np.inf, -np.inf], np.nan).dropna()
    if len(ratio) and np.isfinite(thr):
        log(f"    {g}: mean among expressing = {ratio.mean():.2f} "
            f"+/- {ratio.std():.2f} (range {ratio.min():.2f} to {ratio.max():.2f}, "
            f"n = {len(ratio)})")
        log(f"        a threshold of {thr:.2f} corresponds to about "
            f"{100*thr/ratio.mean():.0f}% of basal cells expressing {g} "
            f"({one_in(100*thr/ratio.mean())} cells)")

# (rank column, flag column, raw column, header, block, formatter, sub column)
DIR_SHORT = (f"\u2265 {100*thr_dir:.0f}%" if dir_rule == 'threshold'
             and np.isfinite(thr_dir) else "dominant")
AXES_SPEC = [
    ('opp',      'opportunity',   'pct_umi_hpv16',
     f"Opportunity\nHPV16 \u2265 {thr_opp:.2f}%\nof UMIs", 'sbs2', 'pct',  None),
    ('dir_sbs2', 'dir_sbs2_flag', 'maint_frac',
     f"Direction\nmaintenance {DIR_SHORT}\nof coding reads", 'sbs2', 'frac', None),
    ('cap_sbs2', 'cap_sbs2_flag', 'A3A_mean',
     f"Capability\nA3A > {thr_a:.2f}", 'sbs2', 'expr', 'A3A_pct_pos'),
    ('opp',      'opportunity',   'pct_umi_hpv16',
     f"Opportunity\nHPV16 \u2265 {thr_opp:.2f}%\nof UMIs", 'cnv',  'pct',  None),
    ('dir_cnv',  'dir_cnv_flag',  'prod_frac',
     f"Direction\nproductive {DIR_SHORT}\nof coding reads", 'cnv',  'frac', None),
    ('cap_cnv',  'cap_cnv_flag',  'A3B_mean',
     f"Capability\nA3B > {thr_b:.2f}", 'cnv',  'expr', 'A3B_pct_pos'),
]
missing = sorted({c for spec in AXES_SPEC for c in spec[:3] if c not in cj.columns})
if missing:
    sys.exit(f"  ERROR: patient_conjunction_model.tsv missing {missing}. "
             f"Re-run Diagnostic_Patient_Determinants_Table.py.")

# Row order: the SUM of the two folds, so every driver sits at the top and the
# fill fades toward white down the page.
cj['fold_sum'] = cj['fold_sbs2'].fillna(0) + cj['fold_cnv'].fillna(0)
cj = cj.sort_values('fold_sum', ascending=False).reset_index(drop=True)
n_rows, n_cols = len(cj), len(AXES_SPEC)
_n_drv = int((cj['is_sbs2_driver'] | cj['is_cnv_driver']).sum())
_top = cj.head(_n_drv)
log("")
log(f"  rows ordered by summed fold; top {_n_drv}: "
    f"{[short(p) for p in _top['patient']]} "
    f"(all drivers: {bool((_top['is_sbs2_driver'] | _top['is_cnv_driver']).all())})")

cmap_sbs2 = mcolors.LinearSegmentedColormap.from_list(
    'sbs2', ['#ffffff', '#f7a79c', COLOR_SBS2, '#a3352a'])
cmap_cnv = mcolors.LinearSegmentedColormap.from_list(
    'cnv', ['#ffffff', '#f6e2a8', COLOR_CNV, COLOR_CNV_DARK, '#6b5416'])

fig, ax = plt.subplots(figsize=(24, 2.8 + 1.15 * n_rows))

# Pass 1: fills.
for j, (rank_c, flag_c, raw_c, _, block, fmt, sub_c) in enumerate(AXES_SPEC):
    cm = cmap_sbs2 if block == 'sbs2' else cmap_cnv
    for i in range(n_rows):
        v = cj.loc[i, rank_c]
        face = cm(float(v)) if np.isfinite(v) else '#f2f2f2'
        ax.add_patch(Rectangle((j, i), 1, 1, facecolor=face,
                               edgecolor='#c8c8c8', linewidth=1.0, zorder=2))

# Pass 2: values. Capability cells carry a second line giving penetrance.
for j, (rank_c, flag_c, raw_c, _, block, fmt, sub_c) in enumerate(AXES_SPEC):
    for i in range(n_rows):
        v = cj.loc[i, rank_c]
        raw = cj.loc[i, raw_c]
        dark = np.isfinite(v) and v > 0.60
        col = 'white' if dark else COLOR_TEXT
        if not np.isfinite(raw):
            txt = 'n.d.'
        elif fmt == 'pct':
            txt = f"{raw:.2f}%"
        elif fmt == 'frac':
            txt = f"{100*raw:.0f}%"
        else:
            txt = f"{raw:.2f}"
        dy = 0.5 if sub_c is None else 0.38
        ax.text(j + 0.5, i + dy, txt, ha='center', va='center',
                fontsize=FONT_CELL, zorder=4, color=col)
        if sub_c is not None:
            pctpos = cj.loc[i, sub_c] if sub_c in cj.columns else np.nan
            ax.text(j + 0.5, i + 0.70, one_in(pctpos), ha='center', va='center',
                    fontsize=FONT_SUB, zorder=4,
                    color='#f0f0f0' if dark else '#666666')

# Pass 3: threshold borders, complete and on top.
for j, (rank_c, flag_c, raw_c, _, block, fmt, sub_c) in enumerate(AXES_SPEC):
    for i in range(n_rows):
        if bool(cj.loc[i, flag_c]):
            ax.add_patch(Rectangle((j, i), 1, 1, fill=False, edgecolor='black',
                                   linewidth=4.5, joinstyle='miter', zorder=9))

ax.axvline(3, color='black', linewidth=5.0, zorder=10)
ax.text(1.5, -1.05, 'SBS2-HIGH', ha='center', va='bottom',
        fontsize=FONT_TITLE, fontweight='bold', color=COLOR_SBS2)
ax.text(4.5, -1.05, 'CNV-HIGH', ha='center', va='bottom',
        fontsize=FONT_TITLE, fontweight='bold', color=COLOR_CNV_DARK)

ax.set_xlim(0, n_cols + 2.3)
ax.set_ylim(n_rows, -2.2)
ax.set_xticks(np.arange(n_cols) + 0.5)
ax.set_xticklabels([h for _, _, _, h, _, _, _ in AXES_SPEC],
                   fontsize=FONT_TICK - 8, linespacing=1.35)
ax.set_yticks(np.arange(n_rows) + 0.5)
ax.set_yticklabels([short(p) for p in cj['patient']], fontsize=FONT_TICK)
ax.tick_params(length=0)
for sp_ in ax.spines.values():
    sp_.set_visible(False)

for i, r in cj.iterrows():
    lab = ax.get_yticklabels()[i]
    s, c = bool(r['is_sbs2_driver']), bool(r['is_cnv_driver'])
    if s and c:
        lab.set_color(COLOR_BOTH); lab.set_fontweight('bold')
    elif s:
        lab.set_color(COLOR_SBS2); lab.set_fontweight('bold')
    elif c:
        lab.set_color(COLOR_CNV_DARK); lab.set_fontweight('bold')

ax.text(n_cols + 0.55, -0.75, 'SBS2\nfold', ha='center', va='bottom',
        fontsize=FONT_TICK - 2, fontweight='bold', color=COLOR_SBS2)
ax.text(n_cols + 1.65, -0.75, 'CNV\nfold', ha='center', va='bottom',
        fontsize=FONT_TICK - 2, fontweight='bold', color=COLOR_CNV_DARK)
for i, r in cj.iterrows():
    ax.text(n_cols + 0.55, i + 0.5, f"{r['fold_sbs2']:.2f}x", ha='center',
            va='center', fontsize=FONT_CELL,
            fontweight='bold' if r['is_sbs2_driver'] else 'normal',
            color=COLOR_SBS2 if r['is_sbs2_driver'] else COLOR_NEITHER)
    ax.text(n_cols + 1.65, i + 0.5, f"{r['fold_cnv']:.2f}x", ha='center',
            va='center', fontsize=FONT_CELL,
            fontweight='bold' if r['is_cnv_driver'] else 'normal',
            color=COLOR_CNV_DARK if r['is_cnv_driver'] else COLOR_NEITHER)

leg = [Patch(facecolor='white', edgecolor='black', linewidth=4.5,
             label='meets the axis threshold'),
       Patch(facecolor='white', edgecolor='#c8c8c8', linewidth=1.0,
             label='below the axis threshold')]
ax.legend(handles=leg, loc='upper center', bbox_to_anchor=(0.5, -0.105),
          ncol=2, fontsize=FONT_LEGEND, frameon=False)

fig.text(0.5, -0.115,
         "Fill = within-cohort percentile rank on that axis. Opportunity is the same "
         "column in both blocks. Capability cells give the mean expression the "
         "threshold acts on, and beneath it how many cells express the enzyme.",
         ha='center', fontsize=FONT_LEGEND - 2, color=COLOR_TEXT)
fig.tight_layout()
save(fig, "Supp_P7_Conjunction_Heatmap")

keep = (['patient', 'pct_umi_hpv16', 'maint_frac', 'prod_frac',
         'A3A_mean', 'A3A_pct_pos', 'A3B_mean', 'A3B_pct_pos',
         'opportunity', 'dir_sbs2_flag', 'cap_sbs2_flag',
         'dir_cnv_flag', 'cap_cnv_flag', 'all3_sbs2', 'all3_cnv',
         'fold_sbs2', 'fold_cnv', 'is_sbs2_driver', 'is_cnv_driver'])
keep = [c for c in keep if c in cj.columns]
out7 = cj[keep].copy()
out7['thr_opportunity_pct_umi'] = thr_opp
out7['thr_direction'] = thr_dir
out7['thr_a3a'] = thr_a
out7['thr_a3b'] = thr_b
source(out7, "Supp_P7_Conjunction_Heatmap")

all3_s = sorted(short(p) for p in cj.loc[cj['all3_sbs2'], 'patient'])
all3_c = sorted(short(p) for p in cj.loc[cj['all3_cnv'], 'patient'])
miss_c = sorted(short(p) for p in cj.loc[cj['is_cnv_driver'] & ~cj['all3_cnv'], 'patient'])
caption("P7",
        f"Three-axis conjunction for both fates, {n_rows} patients. A heavy border "
        f"marks a patient meeting that axis: HPV16 at or above {thr_opp:.2f}% of "
        f"basal-compartment UMIs, {DIR_HDR} in the matching stage with the URR "
        f"excluded, and mean expression of the matching enzyme above {thr_a:.2f} for "
        f"A3A or {thr_b:.2f} for A3B. Capability cells give the mean expression the "
        f"threshold acts on and, beneath it, how many basal cells express the enzyme; "
        f"the two enzyme thresholds differ, so the same value can clear A3A and fail "
        f"A3B. Opportunity is the same column in both blocks. All three axes are met "
        f"by {', '.join(all3_s)} on the SBS2 side and {', '.join(all3_c)} on the CNV "
        f"side"
        + (f", while {', '.join(miss_c)} drives CNV-HIGH without meeting all three"
           if miss_c else "")
        + f". Cell fill is the within-cohort percentile rank on that axis. Rows are "
        f"ordered by the sum of the two folds, which places every driver at the top. "
        f"Folds use the {denom} denominator.")


# =============================================================================
# PANEL 8: THRESHOLD SENSITIVITY
# =============================================================================
banner("PANEL 8: THRESHOLD SENSITIVITY SWEEP")

loads = sorted(sweep['load_pctl'].unique())
enzs = sorted(sweep['enz_pctl'].unique())
fig, axes = plt.subplots(1, 2, figsize=(26, 13))

for ax, tag, title, accent in ((axes[0], 'sbs2', 'SBS2-HIGH', COLOR_SBS2),
                               (axes[1], 'cnv', 'CNV-HIGH', COLOR_CNV_DARK)):
    err = np.zeros((len(loads), len(enzs)))
    for a, lp in enumerate(loads):
        for b, ep in enumerate(enzs):
            r = sweep[(sweep.load_pctl == lp) & (sweep.enz_pctl == ep)].iloc[0]
            err[a, b] = r[f'{tag}_missed'] + r[f'{tag}_false']
    cm = mcolors.LinearSegmentedColormap.from_list(f'{tag}_err', [accent, '#ffffff'])
    ax.imshow(err, cmap=cm, vmin=0, vmax=max(err.max(), 1), origin='lower',
              aspect='auto')
    for a, lp in enumerate(loads):
        for b, ep in enumerate(enzs):
            r = sweep[(sweep.load_pctl == lp) & (sweep.enz_pctl == ep)].iloc[0]
            perfect = bool(r[f'{tag}_perfect'])
            ax.text(b, a,
                    f"{int(r[f'{tag}_hit'])}/{int(r[f'{tag}_missed'])}/"
                    f"{int(r[f'{tag}_false'])}" + ("\n\u2605" if perfect else ""),
                    ha='center', va='center', fontsize=FONT_CELL,
                    fontweight='bold' if perfect else 'normal', color=COLOR_TEXT)
            if perfect:
                ax.add_patch(Rectangle((b - 0.5, a - 0.5), 1, 1, fill=False,
                                       edgecolor='black', linewidth=4.5))
    n_perf = int(sweep[f'{tag}_perfect'].sum())
    ax.set_title(f"{title}\nperfect separation in {n_perf} of {len(sweep)} "
                 f"threshold combinations",
                 fontsize=FONT_TITLE - 4, fontweight='bold', color=accent, pad=16)
    ax.set_xticks(range(len(enzs)))
    ax.set_xticklabels([f"p{e}" for e in enzs], fontsize=FONT_TICK)
    ax.set_yticks(range(len(loads)))
    ax.set_yticklabels([f"p{l}" for l in loads], fontsize=FONT_TICK)
    ax.set_xlabel("enzyme threshold (percentile)", fontsize=FONT_LABEL - 4)
    ax.set_ylabel("viral load threshold (percentile)", fontsize=FONT_LABEL - 4)

fig.text(0.5, 0.015,
         "Each cell: drivers captured / drivers missed / non-drivers wrongly "
         "flagged.  \u2605 and a heavy border mark perfect separation.",
         ha='center', fontsize=FONT_LEGEND, color=COLOR_TEXT)
fig.tight_layout(rect=[0, 0.045, 1, 1])
save(fig, "Supp_P8_Threshold_Sensitivity")
source(sweep, "Supp_P8_Threshold_Sensitivity")

n_s = int(sweep['sbs2_perfect'].sum())
n_c = int(sweep['cnv_perfect'].sum())
ok = sweep[sweep['sbs2_perfect']]
rng = (f"load p{ok['load_pctl'].min()} to p{ok['load_pctl'].max()} and enzyme "
       f"p{ok['enz_pctl'].min()} to p{ok['enz_pctl'].max()}" if n_s else "no range")
caption("P8",
        f"Sensitivity of the conjunction rule to both thresholds, swept across "
        f"{len(loads)} viral-load and {len(enzs)} enzyme percentiles. The SBS2-HIGH "
        f"rule separates drivers from non-drivers perfectly in {n_s} of "
        f"{len(sweep)} combinations, holding across {rng}. The CNV-HIGH rule "
        f"achieves perfect separation in {n_c} of {len(sweep)} combinations. "
        f"Reporting the range rather than a single cut is what keeps this from "
        f"being a result manufactured at one arbitrary threshold.")


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
banner("SUPPLEMENTAL ENZYME + CONJUNCTION PANELS COMPLETE")
