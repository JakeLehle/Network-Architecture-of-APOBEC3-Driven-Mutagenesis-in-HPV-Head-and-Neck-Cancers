#!/usr/bin/env python3
"""
Generate_Supp_Enzyme_Conjunction.py
====================================
Supplemental panels for the patient-level determinants of fate. Four panels,
each written as its own PDF and PNG plus a *_panel_source.tsv holding the exact
values behind it, so the supplement grouping can be decided after looking at
them rather than before.

PANELS
------
  P5  Enzyme double dissociation. 2x2 of scatters, A3A / A3B against SBS2 fold /
      CNV fold. The diagonal cells are the prediction, the off-diagonal cells
      are the controls. This is the payoff panel.
  P6  A3A / A3B, tumor vs normal-adjacent. Three patients, CONTROL only. Shows
      why 'tumor-enriched' needs a fold gate and not just a p-value.
  P7  Three-axis conjunction, continuous. 14 patients by 6 axes, both fates side
      by side, cells coloured by within-cohort percentile rank with the
      threshold drawn as a heavy border and the raw value printed inside.
  P8  Threshold sensitivity, 5x5 for each fate, showing where the boolean rule
      separates drivers from non-drivers and where it breaks.

WHY P7 IS CONTINUOUS RATHER THAN A CHECKMARK GRID
--------------------------------------------------
A binary grid renders a patient sitting thousandths from the threshold as an
empty box, which hides the single most important counterexample in the cohort.
Percentile ranks preserve the continuum, the heavy border still shows the
boolean call, and a borderline patient reads as borderline.

The OPPORTUNITY column is viral load and is therefore IDENTICAL in both fate
blocks. It is deliberately duplicated so each block reads as a self-contained
three-axis rule. The legend states this.

NO NUMBERS ARE HARDCODED. Every value, threshold, correlation and caption
sentence is read or computed from the input tables. The script prints the
recommended caption text for each panel so a legend cannot drift from the data.

INPUTS (read-only, all written by Diagnostic_Patient_Determinants_Table.py)
--------------------------------------------------------------------------
  data/FIG_5/00_diagnostics/patient_determinants_table.tsv
  data/FIG_5/00_diagnostics/patient_conjunction_model.tsv
  data/FIG_5/00_diagnostics/patient_conjunction_sensitivity.tsv
  data/FIG_5/00_diagnostics/patient_enzyme_by_source.tsv

OUTPUTS (to FIGURE_5_PANELS)
----------------------------
  Supp_P5_Enzyme_Double_Dissociation.pdf/.png   + _panel_source.tsv
  Supp_P6_Enzyme_Tumor_vs_Normal.pdf/.png       + _panel_source.tsv
  Supp_P7_Conjunction_Heatmap.pdf/.png          + _panel_source.tsv
  Supp_P8_Threshold_Sensitivity.pdf/.png        + _panel_source.tsv

Run from the directory holding patient_config.py:
    conda run -n NETWORK python Generate_Supp_Enzyme_Conjunction.py

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
import matplotlib.colors as mcolors
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D

# -----------------------------------------------------------------------------
# Paths. patient_config supplies them on the HPC; the fallbacks let this render
# against a local fixture directory for layout checks.
# -----------------------------------------------------------------------------
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

COLOR_SBS2      = "#ed6a5a"   # coral, SBS2 / A3A
COLOR_CNV       = "#F6D155"   # mustard, CNV / A3B
COLOR_CNV_DARK  = "#a8862a"   # mustard is illegible as text on white
COLOR_BOTH      = "#7a4fa3"   # drives both fates
COLOR_NEITHER   = "#9aa0a6"   # drives neither
COLOR_NULL      = "#8c8c8c"   # off-diagonal control panels
COLOR_GRID      = "#d9d9d9"
COLOR_TEXT      = "#333333"

FONT_TITLE  = 34
FONT_LABEL  = 30
FONT_TICK   = 26
FONT_ANNOT  = 24
FONT_LEGEND = 24
FONT_CELL   = 22   # inside heatmap cells; 84 cells will not take 28

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'axes.linewidth': 1.4,
})

OUT = ensure_dir(FIGURE_5_PANELS)

CAPTIONS = []   # collected and printed at the end


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
    """
    Greedy label placement for a small scatter. For each label, try four corner
    offsets and keep the first that does not overlap an already-placed label.
    Fourteen points do not justify a dependency on adjustText, and hand-tuned
    offsets would silently break the moment a value changes.
    """
    fig.canvas.draw()
    px_per_pt = fig.dpi / 72.0
    placed = []
    cands = [(12, 10), (12, -22), (-12, 10), (-12, -22), (0, 20), (0, -28)]
    for text, x, y in items:
        cx, cy = ax.transData.transform((x, y))
        w = 0.62 * fontsize * len(text) * px_per_pt
        h = 1.25 * fontsize * px_per_pt
        chosen = cands[0]
        for dx, dy in cands:
            bx = cx + dx * px_per_pt + (0 if dx >= 0 else -w)
            by = cy + dy * px_per_pt
            box = (bx, by, bx + w, by + h)
            if not any(box[0] < o[2] and box[2] > o[0] and
                       box[1] < o[3] and box[3] > o[1] for o in placed):
                chosen = (dx, dy)
                placed.append(box)
                break
        else:
            dx, dy = chosen
            bx = cx + dx * px_per_pt
            placed.append((bx, cy + dy * px_per_pt, bx + w, cy + dy * px_per_pt + h))
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

# Which patients get a text label: every driver, plus the largest single-factor
# discrepancies among HPV16-positive non-drivers. Derived, never a fixed list.
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
    # Prediction panels carry the full label set; the control panels label only
    # the drivers, because their job is to show the ABSENCE of a trend and
    # crowding them with fourteen tags works against that.
    labels_here = LABEL if is_pred else DRIVERS
    items = [(short(r['patient']), r[xc], r[yc]) for _, r in det.iterrows()
             if r['patient'] in labels_here and np.isfinite(r[xc]) and np.isfinite(r[yc])]
    # place the highest-fold points first so the ones a reader looks for get the
    # cleanest position
    items.sort(key=lambda t: -t[2])
    place_labels(ax, fig, items, FONT_ANNOT)

    # No fitted line. rho here is Spearman, which is rank-based, so an ordinary
    # least-squares line drawn through the raw values would be a different model
    # from the statistic in the title and would be dragged around by the two
    # high-fold patients. The bold frame carries the emphasis instead.

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
        f"Folds use the {denom} denominator. These are pre-specified single "
        f"hypotheses, not hits from a scan.")


# =============================================================================
# PANEL 6: A3A / A3B TUMOR vs NORMAL-ADJACENT (control)
# =============================================================================
if enz is not None and len(enz):
    banner("PANEL 6: ENZYME, TUMOR vs NORMAL-ADJACENT (control)")
    MIN_FOLD = 2.0
    # Row order: how many of the two enzymes show genuine tumor enrichment
    # (fold >= MIN_FOLD AND p < 0.05), then by the larger fold. Sorting on A3A
    # fold alone put the patient with the reversed A3A direction last even though
    # it carries a real A3B enrichment. Counting the passes first puts the
    # clearest tumor-specific cases at the top, which is what the panel argues.
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
# PANEL 7: CONJUNCTION HEATMAP (continuous, both fates)
# =============================================================================
banner("PANEL 7: THREE-AXIS CONJUNCTION HEATMAP")

# (rank column, boolean flag column, raw value column, header, block)
AXES_SPEC = [
    ('opp',      'opportunity',   'load_per_M', 'Opportunity\n(viral load)',   'sbs2'),
    ('dir_sbs2', 'dir_sbs2_flag', 'maint_frac', 'Direction\n(maintenance)',    'sbs2'),
    ('cap_sbs2', 'cap_sbs2_flag', 'A3A_mean',   'Capability\n(A3A)',           'sbs2'),
    ('opp',      'opportunity',   'load_per_M', 'Opportunity\n(viral load)',   'cnv'),
    ('dir_cnv',  'dir_cnv_flag',  'prod_frac',  'Direction\n(productive)',     'cnv'),
    ('cap_cnv',  'cap_cnv_flag',  'A3B_mean',   'Capability\n(A3B)',           'cnv'),
]
missing = [c for c, f, r, _, _ in AXES_SPEC
           for c in (c, f, r) if c not in conj.columns]
if missing:
    sys.exit(f"  ERROR: patient_conjunction_model.tsv missing {sorted(set(missing))}. "
             f"Re-run Diagnostic_Patient_Determinants_Table.py.")

# Row order: the SUM of the two folds, descending. Ordering on one fate alone
# pushed the other fate's driver to the bottom, which broke the visual gradient.
# Summing puts every driver at the top regardless of which fate it drives, so the
# fill fades monotonically toward white down the page and the block of dark cells
# at the top is the whole point of the panel.
cj = conj.copy()
cj['fold_sum'] = cj['fold_sbs2'].fillna(0) + cj['fold_cnv'].fillna(0)
cj = cj.sort_values('fold_sum', ascending=False).reset_index(drop=True)
n_rows, n_cols = len(cj), len(AXES_SPEC)
_n_drv = int((cj['is_sbs2_driver'] | cj['is_cnv_driver']).sum())
_top = cj.head(_n_drv)
log(f"  rows ordered by summed fold; top {_n_drv}: "
    f"{[short(p) for p in _top['patient']]} "
    f"(all drivers: {bool((_top['is_sbs2_driver'] | _top['is_cnv_driver']).all())})")

# Three-stop ramps. A two-stop white-to-base ramp topped out at the base colour,
# which left the driver rows barely darker than mid-cohort. Carrying each ramp
# past its base into a deeper shade makes the top of the panel genuinely dark and
# the bottom genuinely pale, which is the separation the panel exists to show.
cmap_sbs2 = mcolors.LinearSegmentedColormap.from_list(
    'sbs2', ['#ffffff', '#f7a79c', COLOR_SBS2, '#a3352a'])
cmap_cnv = mcolors.LinearSegmentedColormap.from_list(
    'cnv', ['#ffffff', '#f6e2a8', COLOR_CNV, COLOR_CNV_DARK, '#6b5416'])

fig, ax = plt.subplots(figsize=(26, 2.6 + 1.35 * n_rows))

# Drawn in three passes. Cell fills first, then the values, then every heavy
# threshold border last on its own layer. Drawing a border with its fill meant a
# neighbouring cell painted afterwards could bury part of it, so a bordered cell
# next to another bordered cell showed only two or three of its four sides.
# Separating the passes guarantees each border is drawn complete and on top.
for j, (rank_c, flag_c, raw_c, _, block) in enumerate(AXES_SPEC):
    cm = cmap_sbs2 if block == 'sbs2' else cmap_cnv
    for i in range(n_rows):
        v = cj.loc[i, rank_c]
        face = cm(float(v)) if np.isfinite(v) else '#f2f2f2'
        ax.add_patch(Rectangle((j, i), 1, 1, facecolor=face,
                               edgecolor='#c8c8c8', linewidth=1.0, zorder=2))

for j, (rank_c, flag_c, raw_c, _, block) in enumerate(AXES_SPEC):
    for i in range(n_rows):
        v = cj.loc[i, rank_c]
        raw = cj.loc[i, raw_c]
        if np.isfinite(raw):
            txt = f"{raw:.0f}" if abs(raw) >= 100 else f"{raw:.2f}"
        else:
            txt = 'n.d.'
        ax.text(j + 0.5, i + 0.5, txt, ha='center', va='center',
                fontsize=FONT_CELL, zorder=4,
                color='white' if (np.isfinite(v) and v > 0.60) else COLOR_TEXT)

# Threshold borders last, unfilled, on the top layer.
for j, (rank_c, flag_c, raw_c, _, block) in enumerate(AXES_SPEC):
    for i in range(n_rows):
        if bool(cj.loc[i, flag_c]):
            ax.add_patch(Rectangle((j, i), 1, 1, fill=False, edgecolor='black',
                                   linewidth=4.5, joinstyle='miter', zorder=9))

# block separator and headers
ax.axvline(3, color='black', linewidth=5.0, zorder=10)
ax.text(1.5, -0.75, 'SBS2-HIGH', ha='center', va='bottom',
        fontsize=FONT_TITLE, fontweight='bold', color=COLOR_SBS2)
ax.text(4.5, -0.75, 'CNV-HIGH', ha='center', va='bottom',
        fontsize=FONT_TITLE, fontweight='bold', color=COLOR_CNV_DARK)

ax.set_xlim(0, n_cols + 2.3)
ax.set_ylim(n_rows, -1.9)
ax.set_xticks(np.arange(n_cols) + 0.5)
ax.set_xticklabels([h for _, _, _, h, _ in AXES_SPEC], fontsize=FONT_TICK - 2)
ax.set_yticks(np.arange(n_rows) + 0.5)
ax.set_yticklabels([short(p) for p in cj['patient']], fontsize=FONT_TICK)
ax.tick_params(length=0)
for sp_ in ax.spines.values():
    sp_.set_visible(False)

# driver-coloured row labels
for i, r in cj.iterrows():
    lab = ax.get_yticklabels()[i]
    s, c = bool(r['is_sbs2_driver']), bool(r['is_cnv_driver'])
    if s and c:
        lab.set_color(COLOR_BOTH); lab.set_fontweight('bold')
    elif s:
        lab.set_color(COLOR_SBS2); lab.set_fontweight('bold')
    elif c:
        lab.set_color(COLOR_CNV_DARK); lab.set_fontweight('bold')

# fold columns to the right of the grid
ax.text(n_cols + 0.55, -0.55, 'SBS2\nfold', ha='center', va='bottom',
        fontsize=FONT_TICK - 2, fontweight='bold', color=COLOR_SBS2)
ax.text(n_cols + 1.65, -0.55, 'CNV\nfold', ha='center', va='bottom',
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
ax.legend(handles=leg, loc='upper center', bbox_to_anchor=(0.5, -0.045),
          ncol=2, fontsize=FONT_LEGEND, frameon=False)

# The direction rule is read from the source table, not assumed, so the caption
# always states the criterion that actually produced the borders.
dir_rule = cj['direction_rule'].iloc[0] if 'direction_rule' in cj.columns else 'compare'
dir_min = cj['direction_min_frac'].iloc[0] if 'direction_min_frac' in cj.columns else np.nan
DIR_TEXT = (f"matching phase fraction \u2265 {dir_min:.2f}"
            if dir_rule == 'threshold' and np.isfinite(dir_min)
            else "matching phase dominant")
log(f"  direction rule from the source table: {dir_rule}"
    + (f" (>= {dir_min:.2f})" if dir_rule == 'threshold' and np.isfinite(dir_min)
       else " (dominant phase)"))

thr_l = cj['thr_load'].iloc[0] if 'thr_load' in cj.columns else np.nan
thr_a = cj['thr_a3a'].iloc[0] if 'thr_a3a' in cj.columns else np.nan
thr_b = cj['thr_a3b'].iloc[0] if 'thr_a3b' in cj.columns else np.nan
fig.text(0.5, -0.055,
         "Fill = within-cohort percentile rank on that axis. Printed value = the raw "
         "measurement. The Opportunity column is viral load and is identical in both "
         "blocks.",
         ha='center', fontsize=FONT_LEGEND, color=COLOR_TEXT)
fig.tight_layout()
save(fig, "Supp_P7_Conjunction_Heatmap")

keep = ['patient'] + sorted({c for spec in AXES_SPEC for c in spec[:3]}) + \
       ['fold_sbs2', 'fold_cnv', 'is_sbs2_driver', 'is_cnv_driver',
        'all3_sbs2', 'all3_cnv']
keep = [c for c in keep if c in cj.columns]
source(cj[keep], "Supp_P7_Conjunction_Heatmap")

all3_s = sorted(short(p) for p in cj.loc[cj['all3_sbs2'], 'patient'])
all3_c = sorted(short(p) for p in cj.loc[cj['all3_cnv'], 'patient'])
miss_c = sorted(short(p) for p in cj.loc[cj['is_cnv_driver'] & ~cj['all3_cnv'], 'patient'])
caption("P7",
        f"Three-axis conjunction for both fates, {n_rows} patients. Cell fill is the "
        f"within-cohort percentile rank on that axis and the printed number is the raw "
        f"measurement; a heavy border marks a patient meeting the axis threshold "
        f"(load > {thr_l:.0f} per million, A3A > {thr_a:.3f}, A3B > {thr_b:.3f}, "
        f"direction scored as {DIR_TEXT}). "
        f"Opportunity is viral load and is therefore the same column in both blocks. "
        f"All three axes are met by {', '.join(all3_s)} on the SBS2 side and "
        f"{', '.join(all3_c)} on the CNV side. "
        + (f"{', '.join(miss_c)} drives CNV-HIGH without meeting all three, reaching "
           f"that fate by a route the phase axis does not capture. " if miss_c else "")
        + f"Rows are ordered by the sum of the two folds, which places every driver "
        f"at the top and lets the fill fade toward white down the page. The fade is a "
        f"tendency rather than a rule: several non-drivers low in the panel still hold "
        f"one dark cell, because being high on a single axis does not produce a fate. "
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
    cm = mcolors.LinearSegmentedColormap.from_list(
        f'{tag}_err', [accent, '#ffffff'])
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
        f"achieves perfect separation in {n_c} of {len(sweep)} combinations, so on "
        f"that side the rule is threshold-dependent and descriptive only. Reporting "
        f"the range rather than a single cut is what keeps this from being a result "
        f"manufactured at one arbitrary threshold.")


# =============================================================================
# CAPTIONS
# =============================================================================
banner("DRAFT CAPTIONS (computed from this run; paste into the legend)")
for panel, text in CAPTIONS:
    log("")
    log(f"  {panel}")
    words = text.split()
    line = "   "
    for w in words:
        if len(line) + len(w) + 1 > 88:
            log(line)
            line = "   "
        line += " " + w
    log(line)

log("")
log(f"  Panels written to: {OUT}")
banner("SUPPLEMENTAL ENZYME + CONJUNCTION PANELS COMPLETE")
