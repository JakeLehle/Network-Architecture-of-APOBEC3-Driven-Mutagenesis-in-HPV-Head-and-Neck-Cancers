#!/usr/bin/env python3
"""
Generate_Supp_CellCycle.py
===========================
Supplemental panels for the cell-cycle layer: does APOBEC3 enzyme identity map
onto where a basal cell sits in the cycle? Four panels, each written as its own
PDF and PNG plus a *_panel_source.tsv.

PANELS
------
  P9   Tirosh phase mix within SBS2-HIGH, CNV-HIGH and NORMAL. Descriptive, and
       PARTLY CIRCULAR: CNV-HIGH was selected on CNV score and stemness, which
       track proliferation, so its G2/M enrichment is partly built in. The panel
       says so on its face rather than in a footnote.
  P10  Phase mix by APOBEC3 enzyme, over ALL basal cells. This is the panel that
       carries the claim, because it is not selection-confounded. Three
       contrasts: A3 dominance (the hero), then A3A positivity and A3B
       positivity, which are the literature replication and its extension.
  P11  ccAFv2 confirmation, seven states plus the below-threshold Unknown class.
       An independent classifier reaching the same split, and the panel that
       settles whether the A3A-associated G1 is really quiescence.
  P12  Tumor versus normal-adjacent phase mix. CONTROL only, three patients.

WHY EFFECT SIZE AND NOT p
--------------------------
At tens of thousands of cells every contingency test returns a vanishing
p-value, so every panel here is annotated with the PERCENTAGE-POINT difference
in the target phase. That is the number the legend should quote. Cramer's V is
carried in the source tables as a scale-free cross-check but is deliberately not
drawn: its conventional bands are calibrated for other questions and understate
a shift of this size.

NO NUMBERS ARE HARDCODED. Percentages, sample sizes, test statistics and every
caption sentence are read from the diagnostic outputs.

INPUTS (read-only, written by Diagnostic_Patient_CellCycle_by_Source.py)
------------------------------------------------------------------------
  data/FIG_5/00_diagnostics/cellcycle_phase_mix.tsv
  data/FIG_5/00_diagnostics/cellcycle_ccafv2_mix.tsv
  data/FIG_5/00_diagnostics/cellcycle_a3_tests.tsv
  data/FIG_5/00_diagnostics/patient_cellcycle_by_source.tsv

OUTPUTS (to FIGURE_5_PANELS)
----------------------------
  Supp_P9_CellCycle_by_Group.pdf/.png          + _panel_source.tsv
  Supp_P10_CellCycle_by_A3_Enzyme.pdf/.png     + _panel_source.tsv
  Supp_P11_ccAFv2_Confirmation.pdf/.png        + _panel_source.tsv
  Supp_P12_CellCycle_Tumor_vs_Normal.pdf/.png  + _panel_source.tsv

Run from the directory holding patient_config.py:
    conda run -n NETWORK python Generate_Supp_CellCycle.py

Author: Jake Lehle
"""

import os
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

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

# =============================================================================
# STYLE
# =============================================================================
DPI = 300

PHASE_ORDER = ['G1', 'S', 'G2M']
PHASE_COLORS = {'G1': '#6b8f71', 'S': '#f4a259', 'G2M': '#ed6a5a'}

CCAFV2_ORDER = ['qG0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1',
                'Unknown']
CCAFV2_COLORS = {'qG0': '#3b5c7a', 'G1': '#6b8f71', 'Late G1': '#9bb4c7',
                 'S': '#f4a259', 'S/G2': '#e8823b', 'G2/M': '#ed6a5a',
                 'M/Early G1': '#c04a3b', 'Unknown': '#cccccc'}

COLOR_SBS2     = "#ed6a5a"
COLOR_CNV_DARK = "#a8862a"
COLOR_TEXT     = "#333333"
COLOR_MUTED    = "#888888"

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

PRETTY = {
    'SBS2_HIGH': 'SBS2-HIGH', 'CNV_HIGH': 'CNV-HIGH', 'NORMAL': 'Normal',
    'all_basal': 'all basal', 'A3A_dominant': 'A3A-dominant',
    'A3B_dominant': 'A3B-dominant', 'A3A_pos': 'A3A positive',
    'A3A_neg': 'A3A negative', 'A3B_pos': 'A3B positive',
    'A3B_neg': 'A3B negative', 'tumor': 'tumor',
    'normal_adj': 'normal-adjacent',
}


def pretty(s):
    return PRETTY.get(s, str(s).replace('_', ' '))


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
    if p == 0:
        return 'p < 1e-300'
    if p < 1e-3:
        return f"p = {p:.0e}".replace('e-0', 'e-')
    return f"p = {p:.3f}"


def read_tsv(fname):
    path = os.path.join(DIR_00_DIAG, fname)
    if not os.path.exists(path):
        sys.exit(f"  ERROR: required input missing: {path}\n"
                 f"  Run Diagnostic_Patient_CellCycle_by_Source.py first.")
    df = pd.read_csv(path, sep='\t')
    log(f"  [OK] {fname}: {len(df)} rows")
    return df


def wide(mix, stratum_type, order, states):
    """Long-format mix -> {stratum: {state: pct}} plus per-stratum n."""
    sub = mix[mix['stratum_type'] == stratum_type]
    out, ns = {}, {}
    for s in order:
        rows = sub[sub['stratum'] == s]
        if rows.empty:
            continue
        out[s] = {st: float(rows.loc[rows['state'] == st, 'pct'].sum())
                  for st in states}
        ns[s] = int(rows['n_total'].iloc[0])
    return out, ns


def stacked(ax, labels, data, states, colors, gaps=None, bar_h=0.68):
    """Horizontal stacked bars, top row first. `gaps` inserts blank rows."""
    gaps = gaps or set()
    ypos, ylab, k = [], [], 0
    for i, lab in enumerate(labels):
        if i in gaps:
            k += 0.55
        ypos.append(k)
        ylab.append(lab)
        k += 1
    ypos = [max(ypos) - y for y in ypos]        # first label at the top
    left = np.zeros(len(labels))
    for st in states:
        vals = np.array([data[l][st] for l in labels], dtype=float)
        ax.barh(ypos, vals, left=left, height=bar_h, color=colors[st],
                edgecolor='white', linewidth=1.4, label=st, zorder=3)
        left += vals
    ax.set_yticks(ypos)
    ax.set_yticklabels(ylab, fontsize=FONT_TICK)
    ax.set_xlim(0, 100)
    ax.set_xlabel('% of cells', fontsize=FONT_LABEL)
    ax.tick_params(labelsize=FONT_TICK)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    return ypos


# =============================================================================
# LOAD
# =============================================================================
banner("LOAD PANEL INPUTS")
mix = read_tsv("cellcycle_phase_mix.tsv")
tests = read_tsv("cellcycle_a3_tests.tsv")
src = read_tsv("patient_cellcycle_by_source.tsv")
cc_path = os.path.join(DIR_00_DIAG, "cellcycle_ccafv2_mix.tsv")
ccm = read_tsv("cellcycle_ccafv2_mix.tsv") if os.path.exists(cc_path) else None
if ccm is None:
    log("  [SKIP] cellcycle_ccafv2_mix.tsv absent; P11 will be skipped. "
        "Re-run the cell-cycle diagnostic in an environment where ccAFv2 loads.")

T = tests.set_index('test')


def stat_of(name, field):
    return T.loc[name, field] if name in T.index else np.nan


def note_of(name):
    return str(T.loc[name, 'note']) if name in T.index else ''


# =============================================================================
# PANEL 9: PHASE MIX BY GROUP
# =============================================================================
banner("PANEL 9: CELL-CYCLE MIX BY GROUP")
GORDER = ['SBS2_HIGH', 'CNV_HIGH', 'NORMAL']
gd, gn = wide(mix, 'group', GORDER, PHASE_ORDER)
missing = [g for g in GORDER if g not in gd]
if missing:
    sys.exit(f"  ERROR: cellcycle_phase_mix.tsv missing group strata {missing}")

labels = [f"{pretty(g)}\n(n = {gn[g]:,})" for g in GORDER]
fig, ax = plt.subplots(figsize=(19, 9))
ypos = stacked(ax, labels, {l: gd[g] for l, g in zip(labels, GORDER)},
               PHASE_ORDER, PHASE_COLORS)
for y, g in zip(ypos, GORDER):
    run = 0.0
    for st in PHASE_ORDER:
        v = gd[g][st]
        if v >= 7:
            ax.text(run + v / 2, y, f"{v:.0f}%", ha='center', va='center',
                    fontsize=FONT_ANNOT, fontweight='bold', color='white',
                    zorder=5)
        run += v
d_sc = {st: gd['SBS2_HIGH'][st] - gd['CNV_HIGH'][st] for st in PHASE_ORDER}
mx = max(abs(v) for v in d_sc.values())
ax.set_title("Cell-cycle phase by population", fontsize=FONT_TITLE, pad=16)
ax.legend(fontsize=FONT_LEGEND, ncol=3, loc='upper center',
          bbox_to_anchor=(0.5, -0.16), frameon=False)
fig.text(0.5, -0.10,
         "CNV-HIGH was selected on CNV score and stemness, which track "
         "proliferation.\nIts G2/M enrichment is therefore partly circular; the "
         "clean test is the enzyme panel.",
         ha='center', fontsize=FONT_ANNOT, color=COLOR_MUTED)
fig.tight_layout()
save(fig, "Supp_P9_CellCycle_by_Group")

p9 = pd.DataFrame([{'group': pretty(g), 'n': gn[g],
                    **{st: gd[g][st] for st in PHASE_ORDER}} for g in GORDER])
source(p9, "Supp_P9_CellCycle_by_Group")
caption("P9",
        f"Tirosh cell-cycle phase within each population "
        f"({', '.join(f'{pretty(g)} n = {gn[g]:,}' for g in GORDER)}). "
        f"SBS2-HIGH is G1-weighted at {gd['SBS2_HIGH']['G1']:.1f}% while CNV-HIGH "
        f"is proliferative, with {gd['CNV_HIGH']['G2M']:.1f}% in G2/M and "
        f"{gd['CNV_HIGH']['S'] + gd['CNV_HIGH']['G2M']:.1f}% in S or G2/M. The "
        f"largest single-phase difference between the two tumor states is "
        f"{mx:.1f} percentage points "
        f"({fmt_p(stat_of('group_x_phase', 'p'))}). This comparison is PARTLY "
        f"CIRCULAR: CNV-HIGH was selected on CNV score and stemness, both of which "
        f"track proliferation, so its G2/M enrichment is in part a property of the "
        f"selection. The non-circular test is the enzyme panel.")


# =============================================================================
# PANEL 10: PHASE MIX BY A3 ENZYME (all basal; the clean test)
# =============================================================================
banner("PANEL 10: CELL-CYCLE MIX BY A3 ENZYME")
dom, domn = wide(mix, 'a3_dominance', ['A3A_dominant', 'A3B_dominant'],
                 PHASE_ORDER)
flg, flgn = wide(mix, 'a3_flag', ['A3A_pos', 'A3A_neg', 'A3B_pos', 'A3B_neg'],
                 PHASE_ORDER)
if len(dom) < 2 or len(flg) < 4:
    sys.exit("  ERROR: cellcycle_phase_mix.tsv missing a3_dominance or a3_flag strata")

STRATA = ['A3A_dominant', 'A3B_dominant', 'A3A_pos', 'A3A_neg',
          'A3B_pos', 'A3B_neg']
alld = {**dom, **flg}
alln = {**domn, **flgn}
labels = [f"{pretty(s)}\n(n = {alln[s]:,})" for s in STRATA]
data = {l: alld[s] for l, s in zip(labels, STRATA)}

fig, ax = plt.subplots(figsize=(21, 15))
ypos = stacked(ax, labels, data, PHASE_ORDER, PHASE_COLORS, gaps={2, 4})
for y, s in zip(ypos, STRATA):
    run = 0.0
    for st in PHASE_ORDER:
        v = alld[s][st]
        if v >= 7:
            ax.text(run + v / 2, y, f"{v:.0f}%", ha='center', va='center',
                    fontsize=FONT_ANNOT, fontweight='bold', color='white',
                    zorder=5)
        run += v

# Effect sizes in percentage points, computed from the panel's own numbers.
d_dom = alld['A3A_dominant']['G1'] - alld['A3B_dominant']['G1']
d_a3a = alld['A3A_pos']['G1'] - alld['A3A_neg']['G1']
d_a3b = alld['A3B_pos']['G2M'] - alld['A3B_neg']['G2M']
or_a3a, p_a3a = stat_of('A3A_vs_G1', 'stat'), stat_of('A3A_vs_G1', 'p')
or_a3b, p_a3b = stat_of('A3B_vs_G2M', 'stat'), stat_of('A3B_vs_G2M', 'p')
p_dom = stat_of('dominance_x_phase', 'p')

brackets = [
    (0, 1, f"G1  {d_dom:+.1f} pp", COLOR_SBS2, True),
    (2, 3, f"G1  {d_a3a:+.1f} pp   OR {or_a3a:.2f}", COLOR_SBS2, False),
    (4, 5, f"G2/M  {d_a3b:+.1f} pp   OR {or_a3b:.2f}", COLOR_CNV_DARK, False),
]
for i, j, text, color, hero in brackets:
    ymid = (ypos[i] + ypos[j]) / 2
    ax.annotate('', xy=(103, ypos[i]), xytext=(103, ypos[j]),
                arrowprops=dict(arrowstyle='-', color=color,
                                linewidth=3.5 if hero else 2.2))
    ax.text(105, ymid, text + ('\n(non-circular)' if hero else ''),
            va='center', ha='left', color=color, fontsize=FONT_ANNOT,
            fontweight='bold' if hero else 'normal')
ax.set_xlim(0, 152)
# Percentages only run to 100; ticks past it are an artefact of the room made
# for the bracket annotations, so they are removed rather than left to imply
# a scale that does not exist.
ax.set_xticks([0, 20, 40, 60, 80, 100])
ax.set_title("Cell-cycle phase by APOBEC3 enzyme, all basal cells",
             fontsize=FONT_TITLE, pad=16)
ax.legend(fontsize=FONT_LEGEND, ncol=3, loc='upper center',
          bbox_to_anchor=(0.38, -0.09), frameon=False)
fig.tight_layout()
save(fig, "Supp_P10_CellCycle_by_A3_Enzyme")

p10 = pd.DataFrame([{'stratum': pretty(s), 'n': alln[s],
                     **{st: alld[s][st] for st in PHASE_ORDER}} for s in STRATA])
source(p10, "Supp_P10_CellCycle_by_A3_Enzyme")
caption("P10",
        f"Tirosh cell-cycle phase by APOBEC3 enzyme across ALL basal cells, which "
        f"is not confounded by group selection. Among A3-expressing cells, "
        f"A3A-dominant cells are {d_dom:+.1f} percentage points more likely to be "
        f"in G1 than A3B-dominant cells ({alld['A3A_dominant']['G1']:.1f}% versus "
        f"{alld['A3B_dominant']['G1']:.1f}%, n = {alln['A3A_dominant']:,} and "
        f"{alln['A3B_dominant']:,}, {fmt_p(p_dom)}). Splitting instead on "
        f"positivity, A3B-positive cells are enriched in G2/M "
        f"({d_a3b:+.1f} pp, OR {or_a3b:.2f}, {fmt_p(p_a3b)}), replicating the "
        f"published association, and A3A-positive cells are enriched in G1 "
        f"({d_a3a:+.1f} pp, OR {or_a3a:.2f}, {fmt_p(p_a3a)}), which is the "
        f"extension. Quote the percentage-point differences rather than the "
        f"p-values, which are driven by the tens of thousands of cells involved.")


# =============================================================================
# PANEL 11: ccAFv2 CONFIRMATION
# =============================================================================
if ccm is not None:
    banner("PANEL 11: ccAFv2 CONFIRMATION")
    states = [s for s in CCAFV2_ORDER if s in set(ccm['state'])]
    cg, cgn = wide(ccm, 'group', GORDER, states)
    cd, cdn = wide(ccm, 'a3_dominance', ['A3A_dominant', 'A3B_dominant'], states)
    ca, can = wide(ccm, 'all', ['all_basal'], states)

    fig, axes = plt.subplots(2, 1, figsize=(21, 15),
                             gridspec_kw={'height_ratios': [3, 2]})
    l1 = [f"{pretty(g)}\n(n = {cgn[g]:,})" for g in GORDER]
    stacked(axes[0], l1, {l: cg[g] for l, g in zip(l1, GORDER)}, states,
            CCAFV2_COLORS)
    axes[0].set_title("ccAFv2 state by population", fontsize=FONT_TITLE - 2,
                      pad=14)
    axes[0].set_xlabel('')

    DOM = ['A3A_dominant', 'A3B_dominant']
    l2 = [f"{pretty(d)}\n(n = {cdn[d]:,})" for d in DOM]
    stacked(axes[1], l2, {l: cd[d] for l, d in zip(l2, DOM)}, states,
            CCAFV2_COLORS)
    axes[1].set_title("ccAFv2 state by APOBEC3 dominance",
                      fontsize=FONT_TITLE - 2, pad=14)

    qg0 = {d: cd[d]['qG0'] for d in DOM if 'qG0' in states}
    unk = ca['all_basal'].get('Unknown', np.nan) if ca else np.nan
    handles = [Patch(facecolor=CCAFV2_COLORS[s], edgecolor='white', label=s)
               for s in states]
    fig.legend(handles=handles, loc='lower center', ncol=min(len(states), 8),
               fontsize=FONT_LEGEND - 2, frameon=False,
               bbox_to_anchor=(0.5, -0.035))
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    save(fig, "Supp_P11_ccAFv2_Confirmation")

    p11 = pd.DataFrame(
        [{'stratum_type': 'group', 'stratum': pretty(g), 'n': cgn[g],
          **cg[g]} for g in GORDER]
        + [{'stratum_type': 'a3_dominance', 'stratum': pretty(d), 'n': cdn[d],
            **cd[d]} for d in DOM]
        + [{'stratum_type': 'all', 'stratum': 'all basal', 'n': can['all_basal'],
            **ca['all_basal']}])
    source(p11, "Supp_P11_ccAFv2_Confirmation")

    qtext = (f"qG0 is {qg0['A3A_dominant']:.0f}% of A3A-dominant cells against "
             f"{qg0['A3B_dominant']:.0f}% of A3B-dominant cells, so the claim is "
             f"A3A-with-G1 and NOT A3A-with-quiescence. " if len(qg0) == 2 else "")
    rest_a = sum(cd['A3A_dominant'].get(s, 0) for s in ('qG0', 'G1', 'Late G1'))
    prol_a = sum(cd['A3A_dominant'].get(s, 0) for s in ('S', 'S/G2', 'G2/M'))
    rest_b = sum(cd['A3B_dominant'].get(s, 0) for s in ('qG0', 'G1', 'Late G1'))
    prol_b = sum(cd['A3B_dominant'].get(s, 0) for s in ('S', 'S/G2', 'G2/M'))
    caption("P11",
            f"ccAFv2, an independently trained seven-state classifier, applied to "
            f"the same cells. It reproduces the split: A3A-dominant cells are "
            f"{rest_a:.0f}% quiescent or G1 against {prol_a:.0f}% in S or G2/M, "
            f"while A3B-dominant cells are {rest_b:.0f}% against {prol_b:.0f}%. "
            + qtext +
            f"The Unknown class holds calls falling below the classifier's 0.5 "
            f"probability threshold and is {unk:.1f}% of all basal cells; it is "
            f"drawn rather than dropped so the denominator is honest.")
else:
    banner("PANEL 11: SKIPPED (no ccAFv2 output on disk)")


# =============================================================================
# PANEL 12: TUMOR vs NORMAL-ADJACENT (control)
# =============================================================================
banner("PANEL 12: TUMOR vs NORMAL-ADJACENT (control)")
both = src[src['n_normal_adj'] > 0].copy()
if both.empty:
    log("  no patient contributes both tissue sources; panel skipped")
else:
    both = both.sort_values('n_normal_adj', ascending=False).reset_index(drop=True)
    pooled, pooled_n = wide(mix, 'tissue_pooled', ['tumor', 'normal_adj'],
                            PHASE_ORDER)

    labels, data, gaps = [], {}, set()
    for _, r in both.iterrows():
        pid = str(r['patient']).replace('Patient ', '')
        for src_tag, pre, n in (('tumor', 'tumor', int(r['n_tumor'])),
                                ('norm', 'normal-adj', int(r['n_normal_adj']))):
            lab = f"{pid}  {pre}\n(n = {n:,})"
            key = 'tumor_' if src_tag == 'tumor' else 'norm_'
            labels.append(lab)
            data[lab] = {st: float(r[f"{key}{st}"]) for st in PHASE_ORDER}
    if pooled:
        gaps.add(len(labels))
        for tag in ('tumor', 'normal_adj'):
            lab = f"POOLED  {pretty(tag)}\n(n = {pooled_n[tag]:,})"
            labels.append(lab)
            data[lab] = pooled[tag]

    fig, ax = plt.subplots(figsize=(19, 1.6 + 1.15 * len(labels)))
    ypos = stacked(ax, labels, data, PHASE_ORDER, PHASE_COLORS, gaps=gaps)
    for y, lab in zip(ypos, labels):
        run = 0.0
        for st in PHASE_ORDER:
            v = data[lab][st]
            if v >= 8:
                ax.text(run + v / 2, y, f"{v:.0f}%", ha='center', va='center',
                        fontsize=FONT_ANNOT - 2, fontweight='bold',
                        color='white', zorder=5)
            run += v
    ax.set_title("Cell-cycle phase, tumor versus normal-adjacent basal cells",
                 fontsize=FONT_TITLE, pad=16)
    ax.legend(fontsize=FONT_LEGEND, ncol=3, loc='upper center',
              bbox_to_anchor=(0.5, -0.09), frameon=False)
    fig.text(0.5, -0.055,
             f"CONTROL, not a result: only {len(both)} of {len(src)} patients "
             f"contribute normal-adjacent basal cells.",
             ha='center', fontsize=FONT_ANNOT, color=COLOR_MUTED)
    fig.tight_layout()
    save(fig, "Supp_P12_CellCycle_Tumor_vs_Normal")

    p12 = pd.DataFrame([{'stratum': l, **data[l]} for l in labels])
    source(p12, "Supp_P12_CellCycle_Tumor_vs_Normal")

    pool_txt = ""
    if pooled:
        dg1 = pooled['normal_adj']['G1'] - pooled['tumor']['G1']
        pool_txt = (f" Pooled across those patients, normal-adjacent basal is "
                    f"{dg1:+.1f} percentage points more G1 than tumor basal "
                    f"({pooled['normal_adj']['G1']:.1f}% versus "
                    f"{pooled['tumor']['G1']:.1f}%).")
    caption("P12",
            f"Tirosh cell-cycle phase split by tissue source, for the "
            f"{len(both)} of {len(src)} patients contributing normal-adjacent "
            f"basal cells "
            f"({', '.join(sorted(str(p).replace('Patient ', '') for p in both['patient']))})."
            + pool_txt +
            f" Normal-adjacent epithelium being less proliferative than tumor is "
            f"the expected direction and serves as a sanity check on the phase "
            f"calls. Small n and only a few patients: present this as a control, "
            f"not a result.")


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
banner("SUPPLEMENTAL CELL-CYCLE PANELS COMPLETE")
