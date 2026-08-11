#!/usr/bin/env python3
"""
Generate_Panel_CNV_Patient_Distribution.py  (v2)
=================================================
Companion to the Figure 5 SBS2-HIGH patient-distribution panel
(Diagnostic_Regenerate_Panel_A_C.py, Panel A), rendered for the CNV-HIGH
(productive) group in the CNV-HIGH mustard palette.

v2 CHANGES
----------
  - The local DENOMINATOR switch is gone. The setting now comes from
    patient_config.CONTRIBUTION_DENOMINATOR via contribution.py, so this panel
    and the SBS2 panel cannot drift apart.
  - Folds and counts are read from cnv_high_patient_contribution.tsv, which now
    carries both fold variants and a denominator tag. Nothing is recomputed
    unless that table is missing.
  - The SBS2 contributor set used for the dual-driver outline is derived from
    the SAME table under the SAME denominator, instead of being read from a
    hardcoded list built under a possibly different convention.
  - No result numbers in the docstring. Everything quoted is from this run.

Why this panel exists
---------------------
Figure 5 Panel A shows which patients over-contribute to SBS2-HIGH. Without the
matching CNV-HIGH view, a reader has no way to see that the two groups are
driven by mostly different patients, which is the patient-level evidence for
divergent fate. This panel puts both driver sets on the page, and makes the
concentration of CNV-HIGH in a small number of patients visible rather than
buried in a supplementary table.

Layout (mirrors the SBS2 panel exactly, so the two can sit side by side)
-----------------------------------------------------------------------
Horizontal stacked bars, one per patient, sorted by CNV-HIGH count. Gray segment
is other basal cells; colored segment is CNV-HIGH cells. Annotation to the right
of each bar gives the raw count and the fold enrichment. Patients clearing the
high-contributor threshold get bold mustard y-axis labels. A patient that is
ALSO an SBS2-HIGH high contributor gets a coral bar outline.

Bars show the reference compartment matching the active denominator: all basal
under 'all_basal', tumor basal under 'tumor'. That keeps the bar and the fold
annotation describing the same population.

Inputs (read-only)
------------------
  data/FIG_5/00_diagnostics/cnv_high_patient_contribution.tsv   (preferred)
      written by Diagnostic_CNV_High_Patient_Contributions.py v2
  falls back to computing from adata_final.h5ad + three_group_assignments.tsv

Outputs (to FIGURE_5_PANELS)
----------------------------
  Supp_Panel_CNV_Patient_Distribution.pdf/.png
  cnv_patient_distribution_panel_source.tsv   (exact values behind the figure)

Run from the directory holding patient_config.py and contribution.py:
    conda run -n NETWORK python Generate_Panel_CNV_Patient_Distribution.py

Author: Jake Lehle
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from patient_config import (
    DIR_00_DIAG, FIGURE_5_PANELS,
    PATIENT_COL, CELLTYPE_COL,
    HIGH_CONTRIBUTORS,
    FONT_TITLE, FONT_LABEL, FONT_TICK,
    COLOR_SBS2_HIGH,
    banner, log, ensure_dir, load_adata, load_three_groups,
)

from contribution import (
    CONTRIBUTION_DENOMINATOR, HC_THRESHOLD,
    both_folds, attach_folds, derive_contributors, announce, short,
)

try:
    from patient_config import CNV_HIGH_CONTRIBUTORS
except ImportError:
    CNV_HIGH_CONTRIBUTORS = None

# =============================================================================
# CONFIG
# =============================================================================
DPI = 300
CNV_TABLE_PATH = os.path.join(DIR_00_DIAG, "cnv_high_patient_contribution.tsv")
NORMAL_SOURCE = 'normal tissue adjucent to head and neck squamous cell carcinoma'

# Reference column that the bars depict, matching the active denominator.
REF_COUNT_COL = 'n_tumor' if CONTRIBUTION_DENOMINATOR == 'tumor' else 'n_basal'
REF_LABEL = ('Number of tumor basal cells' if CONTRIBUTION_DENOMINATOR == 'tumor'
             else 'Number of basal cells')

# =============================================================================
# STYLE  (mirrors the SBS2 panel, mustard instead of coral)
# =============================================================================
FONT_ANNOT = FONT_TICK - 2

COLOR_CNV_HIGH = "#F6D155"     # mustard, CNV-HIGH (matches Figure 6)
COLOR_HC_BAR   = COLOR_CNV_HIGH        # enriched patients
COLOR_NON_HC   = "#f7e3a1"             # pale mustard, non-enriched CNV-HIGH cells
COLOR_BASAL    = "#d4d4d4"             # gray, other basal
COLOR_DUAL     = COLOR_SBS2_HIGH       # coral outline = also an SBS2-HIGH driver

# Mustard on white is low-contrast for text; darken the label color only.
COLOR_HC_TEXT = "#a8862a"

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'font.size': FONT_TICK,
    'axes.titlesize': FONT_TITLE,
    'axes.labelsize': FONT_LABEL,
    'xtick.labelsize': FONT_TICK,
    'ytick.labelsize': FONT_TICK,
})

# =============================================================================
# LOAD
# =============================================================================
banner("CNV-HIGH PATIENT DISTRIBUTION PANEL (v2)")

out_dir = ensure_dir(FIGURE_5_PANELS)

if os.path.exists(CNV_TABLE_PATH):
    log(f"  Loading contribution table: {CNV_TABLE_PATH}")
    df = pd.read_csv(CNV_TABLE_PATH, sep='\t')
    log(f"  {len(df)} patients loaded")
    tag = df['denominator'].iloc[0] if 'denominator' in df.columns else 'not stated'
    log(f"  Table denominator tag: {tag}")
    if tag not in ('not stated', CONTRIBUTION_DENOMINATOR):
        raise SystemExit(
            f"  ERROR: table was written under '{tag}' but the active setting is "
            f"'{CONTRIBUTION_DENOMINATOR}'. Re-run "
            f"Diagnostic_CNV_High_Patient_Contributions.py first.")
    needed = {'patient', 'n_basal', 'n_tumor', 'n_cnv_high'}
    if not needed.issubset(df.columns):
        raise SystemExit(f"  ERROR: table missing columns "
                         f"{needed - set(df.columns)}; re-run "
                         f"Diagnostic_CNV_High_Patient_Contributions.py")
else:
    log(f"  Table not found at {CNV_TABLE_PATH}; computing from scratch ...")
    adata = load_adata()
    sbs2_high, cnv_high, normal = load_three_groups()

    basal = adata[adata.obs[CELLTYPE_COL] == 'basal cell'].copy()
    basal.obs['group'] = 'other'
    basal.obs.loc[basal.obs_names.isin(cnv_high), 'group'] = 'CNV_HIGH'
    basal.obs.loc[basal.obs_names.isin(sbs2_high), 'group'] = 'SBS2_HIGH'
    if 'source_name' in basal.obs.columns:
        basal.obs['tissue_grp'] = np.where(
            basal.obs['source_name'].astype(str) == NORMAL_SOURCE, 'normal', 'tumor')
    else:
        basal.obs['tissue_grp'] = 'tumor'

    rows = []
    for patient in sorted(basal.obs[PATIENT_COL].unique()):
        pmask = basal.obs[PATIENT_COL] == patient
        rows.append({
            'patient': patient,
            'n_basal': int(pmask.sum()),
            'n_tumor': int((pmask & (basal.obs['tissue_grp'] == 'tumor')).sum()),
            'n_cnv_high': int((pmask & (basal.obs['group'] == 'CNV_HIGH')).sum()),
            'n_sbs2_high': int((pmask & (basal.obs['group'] == 'SBS2_HIGH')).sum()),
        })
    df = pd.DataFrame(rows)

# =============================================================================
# FOLD ENRICHMENT, BOTH DENOMINATORS
# =============================================================================
banner("FOLD ENRICHMENT")

n_cnv_total   = int(df['n_cnv_high'].sum())
n_basal_total = int(df['n_basal'].sum())
n_tumor_total = int(df['n_tumor'].sum())

announce(log, n_basal_total, n_tumor_total)
log("")

df['pct_cnv_high'] = 100.0 * df['n_cnv_high'] / n_cnv_total

# Recompute both folds locally when the table predates v2, so this script never
# depends on a column that might not be there.
if not {'fold_cnv_all_basal', 'fold_cnv_tumor'}.issubset(df.columns):
    log("  Fold columns absent; recomputing both from the counts in the table.")
    fc = [both_folds(nc, n_cnv_total, nb, n_basal_total, nt, n_tumor_total)
          for nc, nb, nt in zip(df['n_cnv_high'], df['n_basal'], df['n_tumor'])]
    df['fold_cnv_all_basal'] = [f['all_basal'] for f in fc]
    df['fold_cnv_tumor']     = [f['tumor'] for f in fc]

if not {'fold_sbs2_all_basal', 'fold_sbs2_tumor'}.issubset(df.columns):
    if 'n_sbs2_high' in df.columns:
        n_sbs2_total = int(df['n_sbs2_high'].sum())
        fs = [both_folds(ns, n_sbs2_total, nb, n_basal_total, nt, n_tumor_total)
              for ns, nb, nt in zip(df['n_sbs2_high'], df['n_basal'], df['n_tumor'])]
        df['fold_sbs2_all_basal'] = [f['all_basal'] for f in fs]
        df['fold_sbs2_tumor']     = [f['tumor'] for f in fs]
    else:
        log("  WARNING: no SBS2 counts available; falling back to the "
            "patient_config expectation for the dual-driver outline.")
        df['fold_sbs2_all_basal'] = np.nan
        df['fold_sbs2_tumor'] = np.nan

df = attach_folds(df, 'fold_cnv')
df = attach_folds(df, 'fold_sbs2')

log(f"  CNV-HIGH cells: {n_cnv_total} | basal {n_basal_total:,} "
    f"| tumor basal {n_tumor_total:,}\n")
log(f"  {'Patient':<10s} {'N_CNV':>7s} {'%CNV':>7s} {'fold(all)':>10s} "
    f"{'fold(tumor)':>12s}")
log(f"  {'-'*10} {'-'*7} {'-'*7} {'-'*10} {'-'*12}")
for _, r in df.sort_values('n_cnv_high', ascending=False).iterrows():
    log(f"  {short(r['patient']):<10s} {int(r['n_cnv_high']):>7d} "
        f"{r['pct_cnv_high']:>6.1f}% {r['fold_cnv_all_basal']:>9.2f}x "
        f"{r['fold_cnv_tumor']:>11.2f}x")

log("")
cnv_hc = derive_contributors(df, 'fold_cnv', expected=CNV_HIGH_CONTRIBUTORS,
                             label='CNV-HIGH contributors', logger=log)
log("")
if df['fold_sbs2'].notna().any():
    sbs2_hc = derive_contributors(df, 'fold_sbs2', expected=HIGH_CONTRIBUTORS,
                                  label='SBS2-HIGH contributors', logger=log)
else:
    sbs2_hc = set(HIGH_CONTRIBUTORS)
    log(f"  SBS2-HIGH contributors: patient_config expectation in use "
        f"{sorted(short(p) for p in sbs2_hc)}")

dual = sorted(cnv_hc & sbs2_hc)

# Chi-square, patient x (CNV_HIGH, not-CNV_HIGH) on the matching reference set
cont = np.column_stack([df['n_cnv_high'].values,
                        df[REF_COUNT_COL].values - df['n_cnv_high'].values])
chi2, pval, dof, _ = chi2_contingency(cont)
log(f"\n  chi-square (patient x CNV-HIGH membership, reference "
    f"{REF_COUNT_COL}): chi2 = {chi2:.2f}, df = {dof}, p = {pval:.2e}")

top = df.sort_values('n_cnv_high', ascending=False)
log(f"  Concentration: {short(top.iloc[0]['patient'])} alone = "
    f"{top.iloc[0]['pct_cnv_high']:.1f}% of CNV-HIGH; top two = "
    f"{top.iloc[:2]['pct_cnv_high'].sum():.1f}%")
log(f"  Drives BOTH: {[short(p) for p in dual] or 'none'}")

# =============================================================================
# PLOT
# =============================================================================
banner("PLOTTING")

plot_df = df.sort_values('n_cnv_high', ascending=True).reset_index(drop=True)
n_patients = len(plot_df)

fig, ax = plt.subplots(figsize=(16, 10))

y = np.arange(n_patients)
n_high  = plot_df['n_cnv_high'].values
n_ref   = plot_df[REF_COUNT_COL].values
n_other = n_ref - n_high
patients = plot_df['patient'].values
folds = plot_df['fold_cnv'].values

ax.barh(y, n_other, color=COLOR_BASAL, edgecolor='black',
        linewidth=0.5, label='Other basal', zorder=2)

for i in range(n_patients):
    is_hc   = folds[i] >= HC_THRESHOLD
    is_dual = patients[i] in sbs2_hc
    ax.barh(y[i], n_high[i], left=n_other[i],
            color=COLOR_HC_BAR if is_hc else COLOR_NON_HC,
            edgecolor=COLOR_DUAL if is_dual else 'black',
            linewidth=3.0 if is_dual else 0.5, zorder=3)

x_max = n_ref.max()
for i in range(n_patients):
    if n_high[i] <= 0:
        continue
    fold = folds[i]
    if fold >= HC_THRESHOLD:
        fold_color = COLOR_HC_TEXT
    elif fold >= 1.0:
        fold_color = '#333333'
    else:
        fold_color = '#888888'
    ax.text(n_ref[i] + x_max * 0.01, y[i],
            f'{int(n_high[i])}  ({fold:.1f}x)',
            va='center', ha='left', fontsize=FONT_ANNOT,
            fontweight='bold' if fold >= HC_THRESHOLD else 'normal',
            color=fold_color)

ax.set_yticks(y)
ax.set_yticklabels([short(p) for p in patients])
for i, p in enumerate(patients):
    if folds[i] >= HC_THRESHOLD:
        ax.get_yticklabels()[i].set_color(COLOR_HC_TEXT)
        ax.get_yticklabels()[i].set_fontweight('bold')

ax.set_xlabel(REF_LABEL, fontsize=FONT_LABEL)
ax.set_title('Patient Distribution of CNV-HIGH Basal Cells',
             fontsize=FONT_TITLE, pad=15)

legend_handles = [
    Patch(facecolor=COLOR_BASAL, edgecolor='black', linewidth=0.5,
          label='Other basal'),
    Patch(facecolor=COLOR_HC_BAR, edgecolor='black', linewidth=0.5,
          label=f'CNV-HIGH (enriched, fold \u2265 {HC_THRESHOLD:.0f}x)'),
    Patch(facecolor=COLOR_NON_HC, edgecolor='black', linewidth=0.5,
          label='CNV-HIGH (not enriched)'),
    Patch(facecolor='white', edgecolor=COLOR_DUAL, linewidth=3.0,
          label='Also SBS2-HIGH contributor'),
]
ax.legend(handles=legend_handles, loc='lower right',
          fontsize=FONT_ANNOT, framealpha=0.9)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
for ext in ['pdf', 'png']:
    path = os.path.join(out_dir, f"Supp_Panel_CNV_Patient_Distribution.{ext}")
    fig.savefig(path, dpi=DPI, bbox_inches='tight')
    log(f"  [SAVE] {path}")
plt.close()

src_path = os.path.join(out_dir, "cnv_patient_distribution_panel_source.tsv")
df.assign(chi2=chi2, chi2_p=pval, chi2_dof=dof,
          reference_column=REF_COUNT_COL).to_csv(src_path, sep='\t', index=False)
log(f"  [SAVE] {src_path}")

banner("CNV-HIGH PATIENT DISTRIBUTION PANEL COMPLETE (v2)")
