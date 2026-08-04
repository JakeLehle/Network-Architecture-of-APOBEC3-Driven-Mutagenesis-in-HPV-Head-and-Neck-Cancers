#!/usr/bin/env python3
"""
Generate_Panel_CNV_Patient_Distribution.py
===========================================
Companion to the Figure 5 SBS2-HIGH patient-distribution panel
(Diagnostic_Regenerate_Panel_A_C.py, Panel A), rendered for the CNV-HIGH
(productive) group in the CNV-HIGH mustard palette.

Why this panel exists
---------------------
Figure 5 Panel A shows which patients over-contribute to SBS2-HIGH. Without the
matching CNV-HIGH view, a reader has no way to see that the two groups are driven
by mostly DIFFERENT patients, which is the patient-level evidence for divergent
fate. This panel puts both driver sets on the page:

  SBS2-HIGH drivers : SC013 (6.84x), SC029 (2.06x), SC001 (2.49x)
  CNV-HIGH  drivers : SC027 (8.79x), SC001 (8.11x)
  Drives BOTH       : SC001 only

It also makes visible, rather than buried in a supplementary table, that CNV-HIGH
is concentrated in a single patient: SC027 alone contributes 285 of 546 cells
(52.2%), and SC027 plus SC001 account for 82.1%. That concentration is a real
limitation of every CNV-HIGH claim and is better shown than discovered.

Layout (mirrors the SBS2 panel exactly, so the two can sit side by side)
-----------------------------------------------------------------------
Horizontal stacked bars, one per patient, sorted by CNV-HIGH count. Gray segment
is other basal cells; colored segment is CNV-HIGH cells. Annotation to the right
of each bar gives the raw count and the fold enrichment. Patients clearing the
2-fold high-contributor threshold get bold mustard y-axis labels. A patient that
is ALSO an SBS2-HIGH high contributor gets a coral bar outline, so the single
dual driver (SC001) is identifiable at a glance.

Denominator
-----------
DENOMINATOR = 'tumor' (default) matches the unified patient-thread methods:
SBS2-HIGH and CNV-HIGH are selected only from tumor-tissue basal cells, so the
expected fraction is the patient's share of TUMOR basal, not of all basal.
DENOMINATOR = 'all_basal' reproduces the original Figure 5 Panel A convention.
Both are printed in the log. The difference is under 1% in every case (only
SC003/SC005/SC006 contribute any normal-adjacent basal, 554 cells, ~1% of the
compartment, and none of the three is a tumor-group contributor), but the SBS2
panel and this panel must use the SAME setting or the two fold columns are not
comparable. If the SBS2 panel is left on all-basal, set DENOMINATOR accordingly.

Inputs (read-only)
------------------
  data/FIG_5/00_diagnostics/cnv_high_patient_contribution.tsv   (preferred)
      written by Diagnostic_CNV_High_Patient_Contributions.py
  falls back to computing from adata_final.h5ad + three_group_assignments.tsv

Outputs (to FIGURE_5_PANELS)
----------------------------
  Supp_Panel_CNV_Patient_Distribution.pdf/.png
  cnv_patient_distribution_panel_source.tsv   (exact values behind the figure)

Run from scripts/PATIENT_SPECIFIC_EFFECTS/:
    conda run -n NETWORK python Generate_Panel_CNV_Patient_Distribution.py

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from patient_config import (
    DIR_00_DIAG, FIGURE_5_PANELS,
    PATIENT_COL, CELLTYPE_COL,
    HIGH_CONTRIBUTORS,
    FONT_TITLE, FONT_LABEL, FONT_TICK,
    COLOR_SBS2_HIGH,
    banner, log, ensure_dir, load_adata, load_three_groups,
)

# =============================================================================
# CONFIG
# =============================================================================
DENOMINATOR = 'tumor'        # 'tumor' (unified methods) or 'all_basal' (orig Fig 5A)
HC_THRESHOLD = 2.0           # the Figure 5 high-contributor definition
DPI = 300

CNV_TABLE_PATH = os.path.join(DIR_00_DIAG, "cnv_high_patient_contribution.tsv")
NORMAL_SOURCE = 'normal tissue adjucent to head and neck squamous cell carcinoma'

# =============================================================================
# STYLE  (mirrors the SBS2 panel, mustard instead of coral)
# =============================================================================
FONT_ANNOT = FONT_TICK - 2

COLOR_CNV_HIGH = "#F6D155"     # mustard, CNV-HIGH (matches Figure 6)
COLOR_HC_LABEL = COLOR_CNV_HIGH
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

def short(p):
    return str(p).replace('Patient ', '')

# =============================================================================
# LOAD
# =============================================================================
banner("CNV-HIGH PATIENT DISTRIBUTION PANEL")

out_dir = ensure_dir(FIGURE_5_PANELS)

if os.path.exists(CNV_TABLE_PATH):
    log(f"  Loading contribution table: {CNV_TABLE_PATH}")
    df = pd.read_csv(CNV_TABLE_PATH, sep='\t')
    log(f"  {len(df)} patients loaded")
    needed = {'patient', 'n_basal', 'n_tumor', 'n_cnv_high'}
    if not needed.issubset(df.columns):
        raise SystemExit(f"  ERROR: table missing columns {needed - set(df.columns)}; "
                         f"re-run Diagnostic_CNV_High_Patient_Contributions.py")
else:
    log(f"  Table not found at {CNV_TABLE_PATH}; computing from scratch ...")
    adata = load_adata()
    sbs2_high, cnv_high, normal = load_three_groups()

    basal = adata[adata.obs[CELLTYPE_COL] == 'basal cell'].copy()
    basal.obs['group'] = 'other'
    basal.obs.loc[basal.obs_names.isin(cnv_high), 'group'] = 'CNV_HIGH'
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
        })
    df = pd.DataFrame(rows)

# =============================================================================
# FOLD ENRICHMENT, BOTH DENOMINATORS
# =============================================================================
banner("FOLD ENRICHMENT")

n_cnv_total   = int(df['n_cnv_high'].sum())
n_basal_total = int(df['n_basal'].sum())
n_tumor_total = int(df['n_tumor'].sum())

df['pct_cnv_high'] = 100.0 * df['n_cnv_high'] / n_cnv_total
df['fold_all_basal'] = np.where(
    df['n_basal'] > 0,
    (df['n_cnv_high'] / n_cnv_total) / (df['n_basal'] / n_basal_total), 0.0)
df['fold_tumor'] = np.where(
    df['n_tumor'] > 0,
    (df['n_cnv_high'] / n_cnv_total) / (df['n_tumor'] / n_tumor_total), 0.0)
df['fold'] = df['fold_tumor'] if DENOMINATOR == 'tumor' else df['fold_all_basal']

log(f"  CNV-HIGH cells: {n_cnv_total} | basal {n_basal_total:,} "
    f"| tumor basal {n_tumor_total:,}")
log(f"  Denominator in use: {DENOMINATOR}\n")
log(f"  {'Patient':<10s} {'N_CNV':>7s} {'%CNV':>7s} {'fold(all)':>10s} "
    f"{'fold(tumor)':>12s}  role")
log(f"  {'-'*10} {'-'*7} {'-'*7} {'-'*10} {'-'*12}  ----")
for _, r in df.sort_values('n_cnv_high', ascending=False).iterrows():
    role = []
    if r['fold'] >= HC_THRESHOLD:
        role.append('CNV-HC')
    if r['patient'] in HIGH_CONTRIBUTORS:
        role.append('SBS2-HC')
    log(f"  {short(r['patient']):<10s} {int(r['n_cnv_high']):>7d} "
        f"{r['pct_cnv_high']:>6.1f}% {r['fold_all_basal']:>9.2f}x "
        f"{r['fold_tumor']:>11.2f}x  {'+'.join(role)}")

cnv_hc = set(df.loc[df['fold'] >= HC_THRESHOLD, 'patient'])
dual   = sorted(cnv_hc & set(HIGH_CONTRIBUTORS))

# Chi-square, patient x (CNV_HIGH, not-CNV_HIGH) on the matching reference set
ref_col = 'n_tumor' if DENOMINATOR == 'tumor' else 'n_basal'
cont = np.column_stack([df['n_cnv_high'].values,
                        df[ref_col].values - df['n_cnv_high'].values])
chi2, pval, dof, _ = chi2_contingency(cont)
log(f"\n  chi-square (patient x CNV-HIGH membership): chi2 = {chi2:.2f}, "
    f"df = {dof}, p = {pval:.2e}")

# Concentration, the limitation this panel makes visible
top = df.sort_values('n_cnv_high', ascending=False)
log(f"  Concentration: {short(top.iloc[0]['patient'])} alone = "
    f"{top.iloc[0]['pct_cnv_high']:.1f}% of CNV-HIGH; top two = "
    f"{top.iloc[:2]['pct_cnv_high'].sum():.1f}%")
log(f"  CNV-HIGH contributors (>= {HC_THRESHOLD:.0f}x): "
    f"{[short(p) for p in sorted(cnv_hc)]}")
log(f"  SBS2-HIGH contributors: {[short(p) for p in sorted(HIGH_CONTRIBUTORS)]}")
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
n_other = plot_df['n_basal'].values - n_high
patients = plot_df['patient'].values
folds = plot_df['fold'].values

# Other basal (gray)
ax.barh(y, n_other, color=COLOR_BASAL, edgecolor='black',
        linewidth=0.5, label='Other basal', zorder=2)

# CNV-HIGH segment, per-patient color; coral outline marks a dual driver
for i in range(n_patients):
    is_hc   = folds[i] >= HC_THRESHOLD
    is_dual = patients[i] in HIGH_CONTRIBUTORS
    ax.barh(y[i], n_high[i], left=n_other[i],
            color=COLOR_HC_BAR if is_hc else COLOR_NON_HC,
            edgecolor=COLOR_DUAL if is_dual else 'black',
            linewidth=3.0 if is_dual else 0.5, zorder=3)

# Count + fold annotation
x_max = plot_df['n_basal'].max()
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
    ax.text(plot_df['n_basal'].values[i] + x_max * 0.01, y[i],
            f'{int(n_high[i])}  ({fold:.1f}x)',
            va='center', ha='left', fontsize=FONT_ANNOT,
            fontweight='bold' if fold >= HC_THRESHOLD else 'normal',
            color=fold_color)

# Y labels: bold + mustard-dark for CNV high contributors
ax.set_yticks(y)
ax.set_yticklabels([short(p) for p in patients])
for i, p in enumerate(patients):
    if folds[i] >= HC_THRESHOLD:
        ax.get_yticklabels()[i].set_color(COLOR_HC_TEXT)
        ax.get_yticklabels()[i].set_fontweight('bold')

ax.set_xlabel('Number of Basal Cells', fontsize=FONT_LABEL)
ax.set_title('Patient Distribution of CNV-HIGH Basal Cells',
             fontsize=FONT_TITLE, pad=15)

legend_handles = [
    Patch(facecolor=COLOR_BASAL, edgecolor='black', linewidth=0.5,
          label='Other basal'),
    Patch(facecolor=COLOR_HC_BAR, edgecolor='black', linewidth=0.5,
          label=f'CNV-HIGH (enriched, fold > {HC_THRESHOLD:.0f}x)'),
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

# Exact figure source values
src_path = os.path.join(out_dir, "cnv_patient_distribution_panel_source.tsv")
df.assign(denominator=DENOMINATOR, chi2=chi2, chi2_p=pval).to_csv(
    src_path, sep='\t', index=False)
log(f"  [SAVE] {src_path}")

banner("CNV-HIGH PATIENT DISTRIBUTION PANEL COMPLETE")
