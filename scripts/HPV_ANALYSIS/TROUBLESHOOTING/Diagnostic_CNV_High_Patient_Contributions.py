#!/usr/bin/env python3
"""
Diagnostic_CNV_High_Patient_Contributions.py
=============================================
READ-ONLY. Mirror of the Figure 5 SBS2-HIGH patient-contribution analysis, but
for the CNV-HIGH (productive) group.

Figure 5 showed three patients (SC013, SC029, SC001) contribute >2-fold more
SBS2-HIGH cells than expected from their basal share (chi-square p = 4.46e-309).
This asks the same question of CNV-HIGH: are there patients that over-contribute
to the productive group, and are they the SAME patients that drive SBS2-HIGH or
DIFFERENT ones? Different drivers would be direct patient-level support for the
divergent-fate model (some tumors lean maintenance, others productive).

Method (identical to the Figure 5 test, applied to CNV-HIGH)
------------------------------------------------------------
  - Subset to basal cells; tag SBS2_HIGH / CNV_HIGH / NORMAL from
    three_group_assignments.tsv.
  - Per patient: CNV-HIGH count, % of CNV-HIGH, and fold enrichment
    (observed fraction of CNV-HIGH / patient's share of all basal cells).
  - chi2_contingency on patient x (CNV_HIGH, not-CNV_HIGH); the contingency
    expected counts ARE the basal-share-proportional model, so obs/exp equals
    the fold enrichment. Report obs, exp, obs/exp, standardized residual.
  - High-contributor identification at fold thresholds (Fig 5 used >2-fold).
  - Overlap with the SBS2-HIGH high contributors (patient_config.HIGH_CONTRIBUTORS).

Reuses patient_config (same as the SBS2 diagnostics). Run from
scripts/PATIENT_SPECIFIC_EFFECTS/:
    conda run -n NETWORK python Diagnostic_CNV_High_Patient_Contributions.py

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from patient_config import (
    PATIENT_COL, CELLTYPE_COL, TISSUE_COL,
    DIR_00_DIAG, HIGH_CONTRIBUTORS,
    banner, log, ensure_dir, load_adata, load_three_groups,
)

# =============================================================================
# CONFIG
# =============================================================================
OUT_DIR   = ensure_dir(DIR_00_DIAG)
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
FIG5_THRESHOLD = 2.0          # the Fig 5 high-contributor definition
COLOR_CNV   = '#F6D155'       # mustard, CNV-HIGH
COLOR_BASAL = '#d4d4d4'       # other basal
COLOR_SHARED = '#ed6a5a'      # coral, marks a patient that is ALSO an SBS2 HC
DPI = 300

report_lines = []
def rlog(msg=""):
    log(msg)
    report_lines.append(str(msg))

# =============================================================================
# STEP 0: LOAD + TAG
# =============================================================================
banner("STEP 0: LOAD DATA")
adata = load_adata()
sbs2_high, cnv_high, normal = load_three_groups()

basal = adata[adata.obs[CELLTYPE_COL] == 'basal cell'].copy()
basal.obs['group'] = 'other'
basal.obs.loc[basal.obs_names.isin(sbs2_high), 'group'] = 'SBS2_HIGH'
basal.obs.loc[basal.obs_names.isin(cnv_high), 'group'] = 'CNV_HIGH'
basal.obs.loc[basal.obs_names.isin(normal), 'group'] = 'NORMAL'

# tumor vs normal-adjacent source (SBS2/CNV are drawn only from tumor basal per
# Step00B, so tumor-group folds must use the TUMOR basal denominator).
NORMAL_SOURCE = 'normal tissue adjucent to head and neck squamous cell carcinoma'
if 'source_name' in basal.obs.columns:
    basal.obs['tissue_grp'] = np.where(
        basal.obs['source_name'].astype(str) == NORMAL_SOURCE, 'normal', 'tumor')
else:
    basal.obs['tissue_grp'] = 'tumor'
n_tumor = int((basal.obs['tissue_grp'] == 'tumor').sum())

n_cnv = int((basal.obs['group'] == 'CNV_HIGH').sum())
rlog(f"  Total basal cells: {basal.n_obs:,}")
rlog(f"  CNV-HIGH: {n_cnv}  (SBS2-HIGH: {(basal.obs['group']=='SBS2_HIGH').sum()}, "
     f"NORMAL: {(basal.obs['group']=='NORMAL').sum()})")

patients = sorted(basal.obs[PATIENT_COL].unique())
rlog(f"  Patients: {len(patients)}")

# =============================================================================
# STEP 1: PER-PATIENT CNV-HIGH COUNTS + FOLD ENRICHMENT
# =============================================================================
banner("STEP 1: PER-PATIENT CNV-HIGH COUNTS + FOLD ENRICHMENT")

rows = []
for patient in patients:
    pmask = basal.obs[PATIENT_COL] == patient
    n_basal_p = int(pmask.sum())
    n_tumor_p = int((pmask & (basal.obs['tissue_grp'] == 'tumor')).sum())
    n_grp = int((pmask & (basal.obs['group'] == 'CNV_HIGH')).sum())
    expected_frac = n_tumor_p / n_tumor if n_tumor else 0   # tumor-matched
    observed_frac = n_grp / n_cnv if n_cnv else 0
    fold = observed_frac / expected_frac if expected_frac > 0 else 0.0
    tissues = '; '.join(f"{t}={c}" for t, c in
                        basal.obs.loc[pmask, TISSUE_COL].value_counts().items())
    rows.append({'patient': patient, 'n_basal': n_basal_p, 'n_tumor': n_tumor_p,
                 'n_cnv_high': n_grp,
                 'pct_cnv_high': 100.0 * n_grp / n_cnv if n_cnv else 0.0,
                 'fold_cnv_high': fold, 'tissues': tissues})

df = pd.DataFrame(rows).sort_values('n_cnv_high', ascending=False).reset_index(drop=True)
table_path = os.path.join(OUT_DIR, "cnv_high_patient_contribution.tsv")
df.to_csv(table_path, sep='\t', index=False)
rlog(f"  Saved: {table_path}\n")

rlog(f"  {'Patient':<20s} {'N_basal':>8s} {'N_CNV':>7s} {'%_CNV':>7s} {'Fold':>6s}")
rlog(f"  {'-'*20} {'-'*8} {'-'*7} {'-'*7} {'-'*6}")
for _, r in df.iterrows():
    flag = " <<< SBS2 HC" if r['patient'] in HIGH_CONTRIBUTORS else ""
    rlog(f"  {str(r['patient']):<20s} {int(r['n_basal']):>8d} {int(r['n_cnv_high']):>7d} "
         f"{r['pct_cnv_high']:>6.1f}% {r['fold_cnv_high']:>5.1f}x{flag}")

# per-patient source composition (shows why tumor-matching moves the folds so little)
rlog("\n  Per-patient basal by source (tumor / normal-adjacent):")
for p in patients:
    pm = basal.obs[PATIENT_COL] == p
    nt = int((pm & (basal.obs['tissue_grp'] == 'tumor')).sum())
    nn = int((pm & (basal.obs['tissue_grp'] == 'normal')).sum())
    rlog(f"    {str(p).replace('Patient ',''):<8s} tumor={nt:>5d}  normal-adj={nn:>4d}")

# =============================================================================
# STEP 2: CUMULATIVE CONTRIBUTION CURVE
# =============================================================================
banner("STEP 2: CUMULATIVE CONTRIBUTION CURVE (CNV-HIGH)")
rlog(f"  {'Rank':<5s} {'Patient':<20s} {'N_CNV':>7s} {'Cum_N':>7s} {'Cum_%':>7s}")
rlog(f"  {'-'*5} {'-'*20} {'-'*7} {'-'*7} {'-'*7}")
cum = 0
for rank, (_, r) in enumerate(df.iterrows(), 1):
    cum += int(r['n_cnv_high'])
    cum_pct = 100.0 * cum / n_cnv if n_cnv else 0
    rlog(f"  {rank:<5d} {str(r['patient']):<20s} {int(r['n_cnv_high']):>7d} "
         f"{cum:>7d} {cum_pct:>6.1f}%")
    if cum_pct >= 95:
        break

# =============================================================================
# STEP 3: CHI-SQUARE (patient x CNV-HIGH membership)  <-- the requested test
# =============================================================================
banner("STEP 3: CHI-SQUARE TEST (patient x CNV-HIGH vs tumor-not-CNV-HIGH)")
cont = []
for patient in patients:
    pmask = basal.obs[PATIENT_COL] == patient
    n_high = int((pmask & (basal.obs['group'] == 'CNV_HIGH')).sum())
    n_tumor_p = int((pmask & (basal.obs['tissue_grp'] == 'tumor')).sum())
    cont.append({'patient': patient, 'CNV_HIGH': n_high,
                 'not_CNV_HIGH': n_tumor_p - n_high})   # tumor basal reference
cont_df = pd.DataFrame(cont).set_index('patient')
chi2, pval, dof, expected = chi2_contingency(cont_df.values)
rlog(f"  chi2 = {chi2:.2f}, df = {dof}, p = {pval:.2e}")
rlog(f"  {'SIGNIFICANT' if pval < 0.05 else 'Not significant'}: CNV-HIGH cells are "
     f"{'non-uniformly' if pval < 0.05 else 'roughly uniformly'} distributed across patients")
rlog("")
rlog(f"  {'Patient':<20s} {'Obs':>7s} {'Exp':>8s} {'Obs/Exp':>8s} {'Residual':>9s}")
rlog(f"  {'-'*20} {'-'*7} {'-'*8} {'-'*8} {'-'*9}")
for i, patient in enumerate(cont_df.index):
    obs = cont_df.iloc[i, 0]
    exp = expected[i, 0]
    ratio = obs / exp if exp > 0 else 0
    resid = (obs - exp) / np.sqrt(exp) if exp > 0 else 0
    flag = " <<< SBS2 HC" if patient in HIGH_CONTRIBUTORS else ""
    rlog(f"  {str(patient):<20s} {obs:>7.0f} {exp:>8.1f} {ratio:>7.2f}x {resid:>+8.2f}{flag}")

# =============================================================================
# STEP 4: HIGH-CONTRIBUTOR IDENTIFICATION
# =============================================================================
banner("STEP 4: CNV-HIGH HIGH-CONTRIBUTOR IDENTIFICATION")
for thresh in (1.2, 1.5, 2.0, 2.5, 3.0):
    hc = df[df['fold_cnv_high'] >= thresh]
    cum_pct = 100.0 * hc['n_cnv_high'].sum() / n_cnv if n_cnv else 0
    marker = "  <-- Fig 5 threshold" if thresh == FIG5_THRESHOLD else ""
    rlog(f"\n  Fold >= {thresh:.1f}x: {len(hc)} patient(s), "
         f"{int(hc['n_cnv_high'].sum())} cells ({cum_pct:.1f}% of CNV-HIGH){marker}")
    for _, r in hc.iterrows():
        also = " (also SBS2 HC)" if r['patient'] in HIGH_CONTRIBUTORS else ""
        rlog(f"    {str(r['patient']):<20s} fold={r['fold_cnv_high']:.2f}x  "
             f"n_cnv={int(r['n_cnv_high'])}{also}")

cnv_hc = df[df['fold_cnv_high'] >= FIG5_THRESHOLD]['patient'].tolist()

# =============================================================================
# STEP 5: OVERLAP WITH SBS2-HIGH CONTRIBUTORS  (the divergent-fate question)
# =============================================================================
banner("STEP 5: OVERLAP WITH SBS2-HIGH CONTRIBUTORS")
sbs2_hc = list(HIGH_CONTRIBUTORS)
rlog(f"  SBS2-HIGH contributors (Fig 5): {[str(p).replace('Patient ','') for p in sbs2_hc]}")
rlog(f"  CNV-HIGH  contributors (>={FIG5_THRESHOLD}x): "
     f"{[str(p).replace('Patient ','') for p in cnv_hc]}")
shared = sorted(set(sbs2_hc) & set(cnv_hc))
cnv_only = sorted(set(cnv_hc) - set(sbs2_hc))
sbs2_only = sorted(set(sbs2_hc) - set(cnv_hc))
rlog("")
rlog(f"  Drive BOTH groups: {[str(p).replace('Patient ','') for p in shared] or 'none'}")
rlog(f"  CNV-HIGH only:     {[str(p).replace('Patient ','') for p in cnv_only] or 'none'}")
rlog(f"  SBS2-HIGH only:    {[str(p).replace('Patient ','') for p in sbs2_only] or 'none'}")
rlog("")
if not shared:
    rlog("  READ: fully distinct drivers -> strong patient-level support for divergent")
    rlog("        fate (maintenance-leaning vs productive-leaning tumors).")
elif cnv_only or sbs2_only:
    rlog("  READ: partially distinct drivers -> patients lean toward one fate but with")
    rlog("        overlap; consistent with a lifecycle gradient rather than a clean split.")
else:
    rlog("  READ: identical drivers -> the same patients dominate both groups; the")
    rlog("        divergence is within-patient, not between-patient.")

# =============================================================================
# STEP 6: PLOT (diagnostic; promotable to a supplemental panel)
# =============================================================================
banner("STEP 6: PLOT")
dfp = df.sort_values('n_cnv_high', ascending=True).reset_index(drop=True)
y = np.arange(len(dfp))
n_other = dfp['n_basal'].values - dfp['n_cnv_high'].values
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 9))

ax1.barh(y, n_other, color=COLOR_BASAL, edgecolor='black', linewidth=0.4, label='Other basal')
for i, (_, r) in enumerate(dfp.iterrows()):
    ax1.barh(y[i], r['n_cnv_high'], left=n_other[i], color=COLOR_CNV,
             edgecolor='black', linewidth=0.4)
ax1.set_yticks(y)
ax1.set_yticklabels([str(p).replace('Patient ', '') for p in dfp['patient']], fontsize=9)
ax1.set_xlabel('Number of basal cells', fontsize=12)
ax1.set_title('Basal composition per patient\n(mustard = CNV-HIGH)', fontsize=12)
ax1.legend(loc='lower right', fontsize=9)

folds = dfp['fold_cnv_high'].values
bar_colors = [COLOR_SHARED if p in HIGH_CONTRIBUTORS else COLOR_CNV for p in dfp['patient']]
ax2.barh(y, folds, color=bar_colors, edgecolor='black', linewidth=0.4)
ax2.axvline(1.0, color='black', linestyle='--', linewidth=1.2, label='Expected (1.0x)')
ax2.axvline(FIG5_THRESHOLD, color='#333333', linestyle=':', linewidth=1.2,
            label=f'HC threshold ({FIG5_THRESHOLD:.0f}x)')
ax2.set_yticks(y)
ax2.set_yticklabels([str(p).replace('Patient ', '') for p in dfp['patient']], fontsize=9)
ax2.set_xlabel('Fold enrichment (observed / expected)', fontsize=12)
ax2.set_title('CNV-HIGH enrichment per patient\n(coral outline = also SBS2 contributor)',
              fontsize=12)
ax2.legend(loc='lower right', fontsize=9)

plt.suptitle(f'CNV-HIGH patient contribution (chi-square p = {pval:.2e})',
             fontsize=14, y=1.01)
plt.tight_layout()
for ext in ('pdf', 'png'):
    plt.savefig(os.path.join(OUT_DIR, f"cnv_high_patient_contribution.{ext}"),
                dpi=DPI, bbox_inches='tight')
plt.close()
rlog(f"  [SAVE] cnv_high_patient_contribution.pdf/.png")

# =============================================================================
# SAVE REPORT
# =============================================================================
report_path = os.path.join(OUT_DIR, f"cnv_high_contributor_diagnostic_{TIMESTAMP}.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
rlog(f"\n  Report: {report_path}")
banner("CNV-HIGH CONTRIBUTOR DIAGNOSTIC COMPLETE")
