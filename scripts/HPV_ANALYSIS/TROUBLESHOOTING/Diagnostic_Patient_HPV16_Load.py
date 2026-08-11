#!/usr/bin/env python3
"""
Diagnostic_Patient_HPV16_Load.py  (v2 -- denominator toggle)
=============================================================
READ-ONLY. Per-patient HPV16 viral load, normalized to library size, as a
control on the patient-contribution result.

Motivation. The contributor tests found that the SBS2-HIGH and CNV-HIGH groups
are driven by mostly different patients, with some overlap. One explanation for
a patient driving BOTH fates is simply that it carries a higher overall HPV16
burden, so it seeds cells of both fates regardless of any fate bias. This
measures per-patient viral load and asks whether load explains the dual
contributors while the fate-specific drivers are not load outliers, which would
leave the divergent-fate reading intact.

v2 CHANGES
----------
  - The contribution denominator comes from patient_config.CONTRIBUTION_DENOMINATOR
    via contribution.py. Both folds (all-basal and tumor-basal) are always
    computed and written; the setting picks the headline column.
  - Contributor sets on both sides are derived from the fold columns at runtime.
    The patient_config lists are expectations only, and a mismatch is logged.
  - The STEP 3 READ block no longer names patients or quotes ranks that were
    written in advance. It classifies the derived driver sets against their
    observed load ranks, so the prose cannot drift from the data.
  - SPOTLIGHT remains a hardcoded list because it is an ANALYSIS CHOICE (which
    patients get a detailed printout), not a result.

What it computes (over basal cells, the group-source population)
---------------------------------------------------------------
  - Library-normalized load: sum(raw_HPV16) / sum(total UMIs) per patient, scaled
    (HPV16 UMIs per million total UMIs; falls back to per-1000-genes if total
    UMIs are unavailable, matching the Phase 3 precedent).
  - Positivity rate: fraction of basal cells with raw_HPV16 >= 8 (L-method).
  - Mean load in positive cells (copy-number-in-infected proxy).
  - Per-patient SBS2-HIGH and CNV-HIGH fold enrichment, both denominators.
  - Spearman of load vs SBS2 fold, CNV fold, and total HIGH contribution.

The HPV-negative control patient should read near-zero load (sanity check).

Inputs (read-only)
------------------
  data/FIG_4/00_input/adata_final.h5ad
  data/FIG_4/01_group_selection/three_group_assignments.tsv
  data/FIG_6/01_raw_hpv16_counts/basal_cell_master_table_with_raw_HPV16.tsv

Outputs (to data/FIG_5/00_diagnostics/)
---------------------------------------
  patient_hpv16_load.tsv                counts, load metrics, both folds
  patient_hpv16_load.pdf/.png
  patient_hpv16_load_<ts>.txt

Run from the directory holding patient_config.py and contribution.py:
    conda run -n NETWORK python Diagnostic_Patient_HPV16_Load.py

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from patient_config import (
    PATIENT_COL, CELLTYPE_COL, DIR_00_DIAG, HIGH_CONTRIBUTORS,
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
MASTER_HPV_PATH = ("/master/jlehle/WORKING/2026_NMF_PAPER/"
                   "data/FIG_6/01_raw_hpv16_counts/"
                   "basal_cell_master_table_with_raw_HPV16.tsv")
OUT_DIR   = ensure_dir(DIR_00_DIAG)
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")

HPV16_POS_THRESHOLD = 8          # L-method positivity call
FIG5_THRESHOLD = HC_THRESHOLD    # imported; never redefine locally

# ANALYSIS CHOICE, not a result: which patients get a detailed printout.
SPOTLIGHT = ['Patient SC001', 'Patient SC027', 'Patient SC013', 'Patient SC029']

NORMAL_SOURCE = 'normal tissue adjucent to head and neck squamous cell carcinoma'
DPI = 300

# Role palette
COLOR_BOTH    = '#7a4fa3'   # purple, drives both fates
COLOR_SBS2    = '#ed6a5a'   # coral
COLOR_CNV     = '#F6D155'   # mustard
COLOR_NEITHER = '#9aa0a6'   # gray

report_lines = []


def rlog(msg=""):
    log(msg)
    report_lines.append(str(msg))


# =============================================================================
# STEP 0: LOAD + TAG + RESOLVE LIBRARY SIZE
# =============================================================================
banner("STEP 0: LOAD DATA")
adata = load_adata()
sbs2_high, cnv_high, normal = load_three_groups()

basal = adata[adata.obs[CELLTYPE_COL] == 'basal cell'].copy()
basal.obs['group'] = 'other'
basal.obs.loc[basal.obs_names.isin(sbs2_high), 'group'] = 'SBS2_HIGH'
basal.obs.loc[basal.obs_names.isin(cnv_high), 'group'] = 'CNV_HIGH'
basal.obs.loc[basal.obs_names.isin(normal), 'group'] = 'NORMAL'
n_sbs2 = int((basal.obs['group'] == 'SBS2_HIGH').sum())
n_cnv  = int((basal.obs['group'] == 'CNV_HIGH').sum())

# tumor vs normal-adjacent source. The tumor groups are seeded only from tumor
# basal (Step00B), which is why the 'tumor' denominator exists as an option.
if 'source_name' in basal.obs.columns:
    basal.obs['tissue_grp'] = np.where(
        basal.obs['source_name'].astype(str) == NORMAL_SOURCE, 'normal', 'tumor')
else:
    rlog("  WARNING: 'source_name' absent; every basal cell treated as tumor.")
    basal.obs['tissue_grp'] = 'tumor'

n_basal_total = int(basal.n_obs)
n_tumor = int((basal.obs['tissue_grp'] == 'tumor').sum())
rlog(f"  Basal cells: {n_basal_total:,}")
rlog("")
announce(rlog, n_basal_total, n_tumor)
rlog("")

# raw HPV16 per cell (merge from the HPV master table)
master = pd.read_csv(MASTER_HPV_PATH, sep='\t', index_col=0)
raw_hpv = master['raw_HPV16'].reindex(basal.obs_names)
n_matched = int(raw_hpv.notna().sum())
basal.obs['raw_HPV16'] = raw_hpv.fillna(0).astype(float).values
rlog(f"  raw_HPV16 matched for {n_matched}/{n_basal_total} basal cells "
     f"({100*n_matched/n_basal_total:.1f}%); unmatched set to 0")

# library size: total UMIs if available, else genes detected (Phase 3 precedent)
LIB_COL = next((c for c in ('total_counts', 'n_counts') if c in basal.obs.columns), None)
if LIB_COL is not None:
    basal.obs['_lib'] = basal.obs[LIB_COL].astype(float).values
    NORM_SCALE, NORM_UNIT = 1e6, "HPV16 UMIs per million total UMIs"
    rlog(f"  library size from adata.obs['{LIB_COL}']")
else:
    ng = master['n_genes'].reindex(basal.obs_names).fillna(0).astype(float)
    basal.obs['_lib'] = ng.values
    NORM_SCALE, NORM_UNIT = 1e3, "HPV16 UMIs per 1000 genes detected (fallback)"
    rlog("  WARNING: no total_counts/n_counts in adata.obs; using genes-detected "
         "(n_genes) as the library-size proxy, matching Phase 3.")

patients = sorted(basal.obs[PATIENT_COL].unique())
rlog(f"  Patients: {len(patients)}")

# =============================================================================
# STEP 1: PER-PATIENT LOAD + CONTRIBUTION TABLE
# =============================================================================
banner("STEP 1: PER-PATIENT HPV16 LOAD + CONTRIBUTION")
rows = []
for p in patients:
    m = basal.obs[PATIENT_COL] == p
    sub = basal.obs.loc[m]
    n_basal_p = int(m.sum())
    raw_sum = float(sub['raw_HPV16'].sum())
    lib_sum = float(sub['_lib'].sum())
    load_norm = NORM_SCALE * raw_sum / lib_sum if lib_sum > 0 else 0.0
    n_pos = int((sub['raw_HPV16'] >= HPV16_POS_THRESHOLD).sum())
    pct_pos = 100.0 * n_pos / n_basal_p if n_basal_p else 0.0
    mean_pos = float(sub.loc[sub['raw_HPV16'] >= HPV16_POS_THRESHOLD,
                             'raw_HPV16'].mean()) if n_pos > 0 else 0.0

    # contributions: BOTH denominators; the headline is set in patient_config
    n_tumor_p = int((m & (basal.obs['tissue_grp'] == 'tumor')).sum())
    n_s = int((m & (basal.obs['group'] == 'SBS2_HIGH')).sum())
    n_c = int((m & (basal.obs['group'] == 'CNV_HIGH')).sum())
    fs = both_folds(n_s, n_sbs2, n_basal_p, n_basal_total, n_tumor_p, n_tumor)
    fc = both_folds(n_c, n_cnv,  n_basal_p, n_basal_total, n_tumor_p, n_tumor)

    rows.append({'patient': p,
                 'n_basal': n_basal_p, 'n_tumor': n_tumor_p,
                 'n_normal_adj': n_basal_p - n_tumor_p,
                 'raw_HPV16_total': raw_sum, 'load_norm': load_norm,
                 'n_pos': n_pos, 'pct_pos': pct_pos, 'mean_load_pos': mean_pos,
                 'n_sbs2': n_s, 'n_cnv': n_c,
                 'fold_sbs2_all_basal': fs['all_basal'],
                 'fold_sbs2_tumor':     fs['tumor'],
                 'fold_cnv_all_basal':  fc['all_basal'],
                 'fold_cnv_tumor':      fc['tumor'],
                 'total_HIGH': n_s + n_c})

df = pd.DataFrame(rows)
df = attach_folds(df, 'fold_sbs2')
df = attach_folds(df, 'fold_cnv')
df = df.sort_values('load_norm', ascending=False).reset_index(drop=True)

# ---- Contributor sets, derived from this run --------------------------------
rlog("")
sbs2_hc = derive_contributors(df, 'fold_sbs2', expected=HIGH_CONTRIBUTORS,
                              label='SBS2-HIGH contributors', logger=rlog)
rlog("")
cnv_hc = derive_contributors(df, 'fold_cnv', expected=CNV_HIGH_CONTRIBUTORS,
                             label='CNV-HIGH contributors', logger=rlog)

df['is_sbs2_driver'] = df['patient'].isin(sbs2_hc)
df['is_cnv_driver']  = df['patient'].isin(cnv_hc)
df['role'] = ['+'.join([t for t, v in (('SBS2-HC', s), ('CNV-HC', c)) if v]) or '-'
              for s, c in zip(df['is_sbs2_driver'], df['is_cnv_driver'])]

df.to_csv(os.path.join(OUT_DIR, "patient_hpv16_load.tsv"), sep='\t', index=False)

rlog(f"\n  normalization: {NORM_UNIT}")
rlog(f"  folds: headline column = fold_*_{CONTRIBUTION_DENOMINATOR}\n")
rlog(f"  {'Patient':<10s} {'Load':>9s} {'%pos':>6s} {'MeanPos':>8s} "
     f"{'SBS2(all)':>10s} {'SBS2(tum)':>10s} {'CNV(all)':>9s} {'CNV(tum)':>9s}  role")
rlog(f"  {'-'*10} {'-'*9} {'-'*6} {'-'*8} {'-'*10} {'-'*10} {'-'*9} {'-'*9}  ----")
for _, r in df.iterrows():
    rlog(f"  {short(r['patient']):<10s} {r['load_norm']:>9.1f} {r['pct_pos']:>5.1f}% "
         f"{r['mean_load_pos']:>8.1f} "
         f"{r['fold_sbs2_all_basal']:>9.2f}x {r['fold_sbs2_tumor']:>9.2f}x "
         f"{r['fold_cnv_all_basal']:>8.2f}x {r['fold_cnv_tumor']:>8.2f}x  {r['role']}")

max_gap_s = float((df['fold_sbs2_all_basal'] - df['fold_sbs2_tumor']).abs().max())
max_gap_c = float((df['fold_cnv_all_basal'] - df['fold_cnv_tumor']).abs().max())
rlog(f"\n  Max |all-basal minus tumor|: SBS2 {max_gap_s:.4f}x, CNV {max_gap_c:.4f}x")

# per-patient source composition, with HPV16 load in each source
# (normal-adjacent load should sit near zero, a built-in control)
rlog("\n  Per-patient basal by source, with normalized HPV16 load in each:")
for p in patients:
    pm = basal.obs[PATIENT_COL] == p
    tm = pm & (basal.obs['tissue_grp'] == 'tumor')
    nm = pm & (basal.obs['tissue_grp'] == 'normal')
    nt, nn = int(tm.sum()), int(nm.sum())
    lt_lib = basal.obs.loc[tm, '_lib'].sum()
    nn_lib = basal.obs.loc[nm, '_lib'].sum()
    lt = NORM_SCALE * basal.obs.loc[tm, 'raw_HPV16'].sum() / lt_lib if lt_lib > 0 else 0.0
    ln = NORM_SCALE * basal.obs.loc[nm, 'raw_HPV16'].sum() / nn_lib if nn_lib > 0 else 0.0
    rlog(f"    {short(p):<8s} tumor={nt:>5d} (load {lt:>7.1f})  "
         f"normal-adj={nn:>4d} (load {ln:>6.1f})")

# =============================================================================
# STEP 2: DOES LOAD PREDICT CONTRIBUTION?  (Spearman)
# =============================================================================
banner(f"STEP 2: LOAD vs CONTRIBUTION (Spearman, n={len(df)})")
rlog(f"  n={len(df)} patients; correlations are low-powered, read as directional.\n")

rho_store = {}


def sp(a, b, label):
    rho, p = spearmanr(df[a], df[b])
    rho_store[(a, b)] = (rho, p)
    rlog(f"  {label:<38s} rho={rho:+.3f}  p={p:.3g}")


sp('load_norm', 'total_HIGH', 'load vs total HIGH contribution')
sp('load_norm', 'fold_sbs2',  'load vs SBS2-HIGH fold')
sp('load_norm', 'fold_cnv',   'load vs CNV-HIGH fold')
sp('load_norm', 'pct_pos',    'load vs % HPV16+ (sanity: should be high)')

r_s = rho_store[('load_norm', 'fold_sbs2')][0]
r_c = rho_store[('load_norm', 'fold_cnv')][0]
rlog("")
rlog(f"  Load associates with the two fates at rho {r_s:+.3f} (SBS2) and "
     f"{r_c:+.3f} (CNV),")
rlog(f"  a difference of {abs(r_s - r_c):.3f}. If load set DIRECTION rather than")
rlog("  OPPORTUNITY these would diverge; similar strength on both sides is the")
rlog("  quantitative form of 'load is opportunity, not fate'.")

# =============================================================================
# STEP 3: CONTRIBUTOR SPOTLIGHT + DERIVED READ
# =============================================================================
banner("STEP 3: CONTRIBUTOR SPOTLIGHT")
load_rank = {p: i + 1 for i, p in enumerate(df['patient'])}   # 1 = highest load
n_pat = len(df)
top_third = max(1, n_pat // 3)

for p in SPOTLIGHT:
    if p not in load_rank:
        continue
    r = df[df['patient'] == p].iloc[0]
    rlog(f"  {short(p):<8s} load={r['load_norm']:.1f} "
         f"(rank {load_rank[p]}/{n_pat}), %pos={r['pct_pos']:.1f}%, "
         f"SBS2 {r['fold_sbs2']:.1f}x / CNV {r['fold_cnv']:.1f}x  [{r['role']}]")

# ---- READ, derived from the sets rather than written in advance -------------
rlog("")
rlog(f"  READ (derived; 'high load' = load rank within the top third, "
     f"rank <= {top_third}/{n_pat}):")

dual_drivers = sorted(sbs2_hc & cnv_hc)
fate_specific = sorted((sbs2_hc | cnv_hc) - (sbs2_hc & cnv_hc))

if dual_drivers:
    for p in dual_drivers:
        rk = load_rank.get(p)
        r = df[df['patient'] == p].iloc[0]
        if rk and rk <= top_third:
            rlog(f"    {short(p)} drives BOTH and is a high-load patient "
                 f"(rank {rk}/{n_pat}, {r['load_norm']:.0f}), so its dual")
            rlog(f"      contribution is consistent with load seeding both fates.")
        else:
            rlog(f"    {short(p)} drives BOTH but is NOT a load outlier "
                 f"(rank {rk}/{n_pat}); load does not explain its")
            rlog(f"      dual contribution, so the dual role is genuine.")
else:
    rlog("    No patient drives both fates at this threshold.")

if fate_specific:
    rlog("")
    hi, lo = [], []
    for p in fate_specific:
        rk = load_rank.get(p)
        which = 'SBS2' if p in sbs2_hc else 'CNV'
        (hi if (rk and rk <= top_third) else lo).append((short(p), which, rk))
    for name, which, rk in lo:
        rlog(f"    {name} ({which}-only driver) is NOT high-load (rank {rk}/{n_pat}); "
             f"its bias is a genuine")
        rlog(f"      fate preference rather than a load artifact.")
    for name, which, rk in hi:
        rlog(f"    {name} ({which}-only driver) is also high-load (rank {rk}/{n_pat}); "
             f"its fate bias sits on")
        rlog(f"      top of load rather than being pure load.")
    n_lo = len(lo)
    rlog("")
    rlog(f"    {n_lo} of {len(fate_specific)} fate-specific drivers are not load "
         f"outliers.")
    if n_lo:
        rlog("    The divergent-fate reading holds: load alone does not select which")
        rlog("    fate a patient drives.")

# High-load patients that drive NEITHER fate: the opportunity-without-fate cases
non_drivers_hi = [short(r['patient']) for _, r in df.iterrows()
                  if (not r['is_sbs2_driver'] and not r['is_cnv_driver']
                      and load_rank[r['patient']] <= top_third)]
if non_drivers_hi:
    rlog("")
    rlog(f"    High-load patients driving NEITHER fate: {sorted(non_drivers_hi)}.")
    rlog("    These are the cases that make load necessary but not sufficient.")

# =============================================================================
# STEP 4: PLOT
# =============================================================================
banner("STEP 4: PLOT")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))


def role_color(p):
    s, c = p in sbs2_hc, p in cnv_hc
    if s and c:
        return COLOR_BOTH
    if s:
        return COLOR_SBS2
    if c:
        return COLOR_CNV
    return COLOR_NEITHER


dfp = df.copy()
colors = [role_color(p) for p in dfp['patient']]
y = np.arange(len(dfp))[::-1]
ax1.barh(y, dfp['load_norm'].values, color=colors, edgecolor='black', linewidth=0.4)
ax1.set_yticks(y)
ax1.set_yticklabels([short(p) for p in dfp['patient']], fontsize=9)
ax1.set_xlabel(NORM_UNIT, fontsize=11)
ax1.set_title('Per-patient HPV16 load\n'
              '(coral=SBS2-HC, mustard=CNV-HC, purple=both)', fontsize=11)

ax2.scatter(df['load_norm'], df['total_HIGH'],
            c=[role_color(p) for p in df['patient']], s=70,
            edgecolor='black', zorder=3)
label_set = set(SPOTLIGHT) | sbs2_hc | cnv_hc
for _, r in df.iterrows():
    if r['patient'] in label_set:
        ax2.annotate(short(r['patient']), (r['load_norm'], r['total_HIGH']),
                     fontsize=8, xytext=(4, 4), textcoords='offset points')
rho, pv = rho_store[('load_norm', 'total_HIGH')]
ax2.set_xlabel(NORM_UNIT, fontsize=11)
ax2.set_ylabel('HIGH cells contributed (SBS2 + CNV)', fontsize=11)
ax2.set_title(f'Load vs total contribution\n'
              f'(Spearman rho={rho:+.2f}, p={pv:.2g}, n={len(df)})', fontsize=11)

plt.tight_layout()
for ext in ('pdf', 'png'):
    plt.savefig(os.path.join(OUT_DIR, f"patient_hpv16_load.{ext}"),
                dpi=DPI, bbox_inches='tight')
plt.close()
rlog("  [SAVE] patient_hpv16_load.pdf/.png")

report_path = os.path.join(OUT_DIR, f"patient_hpv16_load_{TIMESTAMP}.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
rlog(f"\n  Denominator: {CONTRIBUTION_DENOMINATOR}")
rlog(f"  Table:  {os.path.join(OUT_DIR, 'patient_hpv16_load.tsv')}")
rlog(f"  Report: {report_path}")
banner("PATIENT HPV16 LOAD DIAGNOSTIC COMPLETE (v2)")
