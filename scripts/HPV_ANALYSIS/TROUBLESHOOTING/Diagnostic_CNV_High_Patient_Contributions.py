#!/usr/bin/env python3
"""
Diagnostic_CNV_High_Patient_Contributions.py  (v2 -- denominator toggle)
========================================================================
READ-ONLY. Mirror of the Figure 5 SBS2-HIGH patient-contribution analysis, but
for the CNV-HIGH (productive) group.

Figure 5 showed that a small number of patients contribute more SBS2-HIGH cells
than expected from their share of the basal compartment. This asks the same
question of CNV-HIGH: are there patients that over-contribute to the productive
group, and are they the SAME patients that drive SBS2-HIGH or DIFFERENT ones?
Different drivers would be direct patient-level support for the divergent-fate
model (some tumors lean maintenance, others productive).

v2 CHANGES
----------
  - The denominator is no longer hardcoded. It comes from
    patient_config.CONTRIBUTION_DENOMINATOR via contribution.py. Both folds
    (all-basal and tumor-basal) are always computed, written and printed; the
    setting picks which one is the headline and which reference set the
    chi-square uses.
  - SBS2-HIGH folds are computed here too, under the SAME denominator, so the
    driver-overlap comparison in STEP 5 compares like with like. Previously the
    CNV side was derived from tumor-denominated folds while the SBS2 side was
    read from a hardcoded list built on all-basal folds.
  - Contributor sets on BOTH sides are derived from the fold columns at
    runtime. patient_config lists are treated as expectations and a mismatch is
    logged loudly rather than silently accepted.
  - Result numbers have been removed from the docstring and the READ blocks.
    Everything printed is computed from this run.

Method
------
  - Subset to basal cells; tag SBS2_HIGH / CNV_HIGH / NORMAL from
    three_group_assignments.tsv.
  - Per patient: group counts, % of group, and fold enrichment under BOTH
    denominators (observed share of the group / patient's share of the
    reference set).
  - chi2_contingency on patient x (CNV_HIGH, not-CNV_HIGH) where the reference
    set matches the active denominator, so the expected counts ARE the
    share-proportional model and obs/exp equals the fold enrichment. Reports
    obs, exp, obs/exp and the standardized residual.
  - High-contributor identification swept across fold thresholds.
  - Overlap between the CNV-HIGH and SBS2-HIGH contributor sets.

Inputs (read-only)
------------------
  data/FIG_4/00_input/adata_final.h5ad
  data/FIG_4/01_group_selection/three_group_assignments.tsv

Outputs (to data/FIG_5/00_diagnostics/)
---------------------------------------
  cnv_high_patient_contribution.tsv          per-patient counts + both folds
  cnv_high_patient_contribution.pdf/.png     diagnostic figure
  cnv_high_contributor_diagnostic_<ts>.txt   full console log

Run from the directory holding patient_config.py and contribution.py:
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

from contribution import (
    CONTRIBUTION_DENOMINATOR, HC_THRESHOLD,
    both_folds, reference_count, attach_folds, derive_contributors,
    announce, short,
)

# Optional expectation for the CNV side. Absent from older patient_config
# versions, so it is imported defensively and only used as an expectation.
try:
    from patient_config import CNV_HIGH_CONTRIBUTORS
except ImportError:
    CNV_HIGH_CONTRIBUTORS = None

# =============================================================================
# CONFIG
# =============================================================================
OUT_DIR   = ensure_dir(DIR_00_DIAG)
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")

COLOR_CNV    = '#F6D155'      # mustard, CNV-HIGH
COLOR_BASAL  = '#d4d4d4'      # other basal
COLOR_SHARED = '#ed6a5a'      # coral, marks a patient that is ALSO an SBS2 HC
DPI = 300

SWEEP_THRESHOLDS = (1.2, 1.5, 2.0, 2.5, 3.0)

NORMAL_SOURCE = 'normal tissue adjucent to head and neck squamous cell carcinoma'

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

if 'source_name' in basal.obs.columns:
    basal.obs['tissue_grp'] = np.where(
        basal.obs['source_name'].astype(str) == NORMAL_SOURCE, 'normal', 'tumor')
else:
    rlog("  WARNING: 'source_name' absent; every basal cell treated as tumor.")
    basal.obs['tissue_grp'] = 'tumor'

n_basal_total = int(basal.n_obs)
n_tumor_total = int((basal.obs['tissue_grp'] == 'tumor').sum())
n_cnv  = int((basal.obs['group'] == 'CNV_HIGH').sum())
n_sbs2 = int((basal.obs['group'] == 'SBS2_HIGH').sum())

rlog(f"  Total basal cells: {n_basal_total:,}")
rlog(f"  CNV-HIGH: {n_cnv}  (SBS2-HIGH: {n_sbs2}, "
     f"NORMAL: {int((basal.obs['group']=='NORMAL').sum())})")
rlog("")
announce(rlog, n_basal_total, n_tumor_total)

patients = sorted(basal.obs[PATIENT_COL].unique())
rlog(f"\n  Patients: {len(patients)}")

# =============================================================================
# STEP 1: PER-PATIENT COUNTS + FOLD ENRICHMENT (BOTH DENOMINATORS)
# =============================================================================
banner("STEP 1: PER-PATIENT COUNTS + FOLD ENRICHMENT")

rows = []
for patient in patients:
    pmask = basal.obs[PATIENT_COL] == patient
    n_basal_p = int(pmask.sum())
    n_tumor_p = int((pmask & (basal.obs['tissue_grp'] == 'tumor')).sum())
    n_cnv_p   = int((pmask & (basal.obs['group'] == 'CNV_HIGH')).sum())
    n_sbs2_p  = int((pmask & (basal.obs['group'] == 'SBS2_HIGH')).sum())

    f_cnv = both_folds(n_cnv_p, n_cnv, n_basal_p, n_basal_total,
                       n_tumor_p, n_tumor_total)
    f_sbs2 = both_folds(n_sbs2_p, n_sbs2, n_basal_p, n_basal_total,
                        n_tumor_p, n_tumor_total)

    tissues = '; '.join(f"{t}={c}" for t, c in
                        basal.obs.loc[pmask, TISSUE_COL].value_counts().items())

    rows.append({
        'patient': patient,
        'n_basal': n_basal_p,
        'n_tumor': n_tumor_p,
        'n_normal_adj': n_basal_p - n_tumor_p,
        'n_cnv_high': n_cnv_p,
        'n_sbs2_high': n_sbs2_p,
        'pct_cnv_high': 100.0 * n_cnv_p / n_cnv if n_cnv else 0.0,
        'pct_sbs2_high': 100.0 * n_sbs2_p / n_sbs2 if n_sbs2 else 0.0,
        'fold_cnv_all_basal': f_cnv['all_basal'],
        'fold_cnv_tumor': f_cnv['tumor'],
        'fold_sbs2_all_basal': f_sbs2['all_basal'],
        'fold_sbs2_tumor': f_sbs2['tumor'],
        'tissues': tissues,
    })

df = pd.DataFrame(rows)
df = attach_folds(df, 'fold_cnv')
df = attach_folds(df, 'fold_sbs2')

# Legacy column name kept so downstream readers that predate v2 still work.
df['fold_cnv_high'] = df['fold_cnv']
df = df.sort_values('n_cnv_high', ascending=False).reset_index(drop=True)

rlog(f"  {'Patient':<10s} {'N_basal':>8s} {'N_tum':>7s} {'N_CNV':>7s} {'%_CNV':>7s} "
     f"{'CNV(all)':>9s} {'CNV(tum)':>9s} {'SBS2(all)':>10s} {'SBS2(tum)':>10s}")
rlog(f"  {'-'*10} {'-'*8} {'-'*7} {'-'*7} {'-'*7} {'-'*9} {'-'*9} {'-'*10} {'-'*10}")
for _, r in df.iterrows():
    rlog(f"  {short(r['patient']):<10s} {int(r['n_basal']):>8d} "
         f"{int(r['n_tumor']):>7d} {int(r['n_cnv_high']):>7d} "
         f"{r['pct_cnv_high']:>6.1f}% "
         f"{r['fold_cnv_all_basal']:>8.2f}x {r['fold_cnv_tumor']:>8.2f}x "
         f"{r['fold_sbs2_all_basal']:>9.2f}x {r['fold_sbs2_tumor']:>9.2f}x")

max_gap_cnv = float((df['fold_cnv_all_basal'] - df['fold_cnv_tumor']).abs().max())
max_gap_sbs2 = float((df['fold_sbs2_all_basal'] - df['fold_sbs2_tumor']).abs().max())
rlog(f"\n  Headline column: 'fold_cnv' = fold_cnv_{CONTRIBUTION_DENOMINATOR}")
rlog(f"  Max |all-basal minus tumor| across patients: "
     f"CNV {max_gap_cnv:.4f}x, SBS2 {max_gap_sbs2:.4f}x")

rlog("\n  Per-patient basal by source (this is why the two denominators differ "
     "so little):")
for _, r in df.sort_values('patient').iterrows():
    rlog(f"    {short(r['patient']):<8s} tumor={int(r['n_tumor']):>5d}  "
         f"normal-adj={int(r['n_normal_adj']):>4d}")

table_path = os.path.join(OUT_DIR, "cnv_high_patient_contribution.tsv")
df.to_csv(table_path, sep='\t', index=False)
rlog(f"\n  Saved: {table_path}")

# =============================================================================
# STEP 2: CUMULATIVE CONTRIBUTION CURVE
# =============================================================================
banner("STEP 2: CUMULATIVE CONTRIBUTION CURVE (CNV-HIGH)")
rlog(f"  {'Rank':<5s} {'Patient':<10s} {'N_CNV':>7s} {'Cum_N':>7s} {'Cum_%':>7s}")
rlog(f"  {'-'*5} {'-'*10} {'-'*7} {'-'*7} {'-'*7}")
cum = 0
for rank, (_, r) in enumerate(df.iterrows(), 1):
    cum += int(r['n_cnv_high'])
    cum_pct = 100.0 * cum / n_cnv if n_cnv else 0.0
    rlog(f"  {rank:<5d} {short(r['patient']):<10s} {int(r['n_cnv_high']):>7d} "
         f"{cum:>7d} {cum_pct:>6.1f}%")
    if cum_pct >= 95:
        break

top = df.sort_values('n_cnv_high', ascending=False)
conc_1 = float(top.iloc[0]['pct_cnv_high'])
conc_2 = float(top.iloc[:2]['pct_cnv_high'].sum())
rlog(f"\n  CONCENTRATION (a limitation of every CNV-HIGH claim, stated rather "
     f"than discovered):")
rlog(f"    {short(top.iloc[0]['patient'])} alone contributes {conc_1:.1f}% of "
     f"the CNV-HIGH group.")
rlog(f"    Top two patients together contribute {conc_2:.1f}%.")

# =============================================================================
# STEP 3: CHI-SQUARE (patient x CNV-HIGH membership)
# =============================================================================
banner(f"STEP 3: CHI-SQUARE (patient x CNV-HIGH; reference = "
       f"{CONTRIBUTION_DENOMINATOR})")

cont_rows = []
for patient in patients:
    pmask = basal.obs[PATIENT_COL] == patient
    n_high = int((pmask & (basal.obs['group'] == 'CNV_HIGH')).sum())
    n_basal_p = int(pmask.sum())
    n_tumor_p = int((pmask & (basal.obs['tissue_grp'] == 'tumor')).sum())
    n_ref_p = reference_count(n_basal_p, n_tumor_p)
    cont_rows.append({'patient': patient, 'CNV_HIGH': n_high,
                      'not_CNV_HIGH': n_ref_p - n_high})

cont_df = pd.DataFrame(cont_rows).set_index('patient')
chi2, pval, dof, expected = chi2_contingency(cont_df.values)
rlog(f"  chi2 = {chi2:.2f}, df = {dof}, p = {pval:.2e}")
rlog(f"  {'SIGNIFICANT' if pval < 0.05 else 'Not significant'}: CNV-HIGH cells "
     f"are {'non-uniformly' if pval < 0.05 else 'roughly uniformly'} "
     f"distributed across patients.")
rlog("")
rlog(f"  {'Patient':<10s} {'Obs':>7s} {'Exp':>8s} {'Obs/Exp':>9s} {'Residual':>9s}")
rlog(f"  {'-'*10} {'-'*7} {'-'*8} {'-'*9} {'-'*9}")
resid_rows = []
for i, patient in enumerate(cont_df.index):
    obs = float(cont_df.iloc[i, 0])
    exp = float(expected[i, 0])
    ratio = obs / exp if exp > 0 else 0.0
    resid = (obs - exp) / np.sqrt(exp) if exp > 0 else 0.0
    resid_rows.append({'patient': patient, 'obs': obs, 'exp': exp,
                       'obs_over_exp': ratio, 'std_residual': resid})
    rlog(f"  {short(patient):<10s} {obs:>7.0f} {exp:>8.1f} {ratio:>8.2f}x "
         f"{resid:>+8.2f}")

resid_df = pd.DataFrame(resid_rows)
df = df.merge(resid_df, on='patient', how='left')
df['chi2'] = chi2
df['chi2_p'] = pval
df['chi2_dof'] = dof
df.to_csv(table_path, sep='\t', index=False)

# =============================================================================
# STEP 4: HIGH-CONTRIBUTOR IDENTIFICATION
# =============================================================================
banner("STEP 4: CNV-HIGH HIGH-CONTRIBUTOR IDENTIFICATION")
for thresh in SWEEP_THRESHOLDS:
    hc = df[df['fold_cnv'] >= thresh]
    cum_pct = 100.0 * hc['n_cnv_high'].sum() / n_cnv if n_cnv else 0.0
    marker = "  <-- HC_THRESHOLD" if abs(thresh - HC_THRESHOLD) < 1e-9 else ""
    rlog(f"\n  Fold >= {thresh:.1f}x: {len(hc)} patient(s), "
         f"{int(hc['n_cnv_high'].sum())} cells ({cum_pct:.1f}% of CNV-HIGH){marker}")
    for _, r in hc.sort_values('fold_cnv', ascending=False).iterrows():
        rlog(f"    {short(r['patient']):<8s} fold={r['fold_cnv']:.2f}x  "
             f"(all {r['fold_cnv_all_basal']:.2f}x / tumor {r['fold_cnv_tumor']:.2f}x)"
             f"  n_cnv={int(r['n_cnv_high'])}")

# =============================================================================
# STEP 5: OVERLAP WITH SBS2-HIGH CONTRIBUTORS (the divergent-fate question)
# =============================================================================
banner("STEP 5: OVERLAP WITH SBS2-HIGH CONTRIBUTORS")
rlog("  Both sets are derived from this run's fold columns under the SAME")
rlog("  denominator, so the overlap compares like with like.\n")

cnv_hc = derive_contributors(df, 'fold_cnv', expected=CNV_HIGH_CONTRIBUTORS,
                             label='CNV-HIGH contributors', logger=rlog)
rlog("")
sbs2_hc = derive_contributors(df, 'fold_sbs2', expected=HIGH_CONTRIBUTORS,
                              label='SBS2-HIGH contributors', logger=rlog)

shared    = sorted(cnv_hc & sbs2_hc)
cnv_only  = sorted(cnv_hc - sbs2_hc)
sbs2_only = sorted(sbs2_hc - cnv_hc)

rlog("")
rlog(f"  Drive BOTH groups: {[short(p) for p in shared] or 'none'}")
rlog(f"  CNV-HIGH only:     {[short(p) for p in cnv_only] or 'none'}")
rlog(f"  SBS2-HIGH only:    {[short(p) for p in sbs2_only] or 'none'}")
rlog("")
if not shared:
    rlog("  READ: fully distinct drivers -> strong patient-level support for")
    rlog("        divergent fate (maintenance-leaning vs productive-leaning tumors).")
elif cnv_only or sbs2_only:
    rlog("  READ: partially distinct drivers -> patients lean toward one fate but")
    rlog("        with overlap; consistent with a lifecycle gradient rather than a")
    rlog("        clean split.")
else:
    rlog("  READ: identical drivers -> the same patients dominate both groups; the")
    rlog("        divergence would then be within-patient, not between-patient.")

# Dual drivers, described from the data rather than named in advance.
for p in shared:
    r = df[df['patient'] == p].iloc[0]
    rlog(f"    {short(p)} drives both: SBS2 {r['fold_sbs2']:.2f}x, "
         f"CNV {r['fold_cnv']:.2f}x")

# =============================================================================
# STEP 6: PLOT (diagnostic; promotable to a supplemental panel)
# =============================================================================
banner("STEP 6: PLOT")
dfp = df.sort_values('n_cnv_high', ascending=True).reset_index(drop=True)
y = np.arange(len(dfp))
n_other = dfp['n_basal'].values - dfp['n_cnv_high'].values

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 9))

ax1.barh(y, n_other, color=COLOR_BASAL, edgecolor='black', linewidth=0.4,
         label='Other basal')
for i, (_, r) in enumerate(dfp.iterrows()):
    ax1.barh(y[i], r['n_cnv_high'], left=n_other[i], color=COLOR_CNV,
             edgecolor='black', linewidth=0.4)
ax1.set_yticks(y)
ax1.set_yticklabels([short(p) for p in dfp['patient']], fontsize=9)
ax1.set_xlabel('Number of basal cells', fontsize=12)
ax1.set_title('Basal composition per patient\n(mustard = CNV-HIGH)', fontsize=12)
ax1.legend(loc='lower right', fontsize=9)

bar_colors = [COLOR_SHARED if p in sbs2_hc else COLOR_CNV for p in dfp['patient']]
ax2.barh(y, dfp['fold_cnv'].values, color=bar_colors, edgecolor='black',
         linewidth=0.4)
ax2.axvline(1.0, color='black', linestyle='--', linewidth=1.2,
            label='Expected (1.0x)')
ax2.axvline(HC_THRESHOLD, color='#333333', linestyle=':', linewidth=1.2,
            label=f'HC threshold ({HC_THRESHOLD:.0f}x)')
ax2.set_yticks(y)
ax2.set_yticklabels([short(p) for p in dfp['patient']], fontsize=9)
ax2.set_xlabel(f'Fold enrichment (denominator: {CONTRIBUTION_DENOMINATOR})',
               fontsize=12)
ax2.set_title('CNV-HIGH enrichment per patient\n(coral = also an SBS2 contributor)',
              fontsize=12)
ax2.legend(loc='lower right', fontsize=9)

plt.suptitle(f'CNV-HIGH patient contribution (chi-square p = {pval:.2e}, '
             f'denominator = {CONTRIBUTION_DENOMINATOR})', fontsize=14, y=1.01)
plt.tight_layout()
for ext in ('pdf', 'png'):
    plt.savefig(os.path.join(OUT_DIR, f"cnv_high_patient_contribution.{ext}"),
                dpi=DPI, bbox_inches='tight')
plt.close()
rlog("  [SAVE] cnv_high_patient_contribution.pdf/.png")

# =============================================================================
# SAVE REPORT
# =============================================================================
report_path = os.path.join(OUT_DIR,
                           f"cnv_high_contributor_diagnostic_{TIMESTAMP}.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
rlog(f"\n  Table:  {table_path}")
rlog(f"  Report: {report_path}")
banner("CNV-HIGH CONTRIBUTOR DIAGNOSTIC COMPLETE (v2)")
