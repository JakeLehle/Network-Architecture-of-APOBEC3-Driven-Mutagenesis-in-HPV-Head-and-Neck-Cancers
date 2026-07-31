#!/usr/bin/env python3
"""
Diagnostic_Patient_HPV16_Load.py
================================
READ-ONLY. Per-patient HPV16 viral load, normalized to library size, as a
control on the patient-contribution result.

Motivation. The CNV-HIGH contributor test found SC027 and SC001 drive CNV-HIGH,
SC013 and SC029 drive SBS2-HIGH, and only SC001 drives BOTH. One explanation for
SC001's dual role is simply that it carries a higher overall HPV16 burden, so it
seeds cells of both fates regardless of any fate bias. This measures per-patient
viral load and asks whether load explains the dual contributor (SC001) while the
fate-specific drivers (SC027, SC013, SC029) are not load outliers, which would
leave the divergent-fate reading intact.

What it computes (over basal cells, the group-source population)
---------------------------------------------------------------
  - Library-normalized load: sum(raw_HPV16) / sum(total UMIs) per patient, scaled
    (HPV16 UMIs per million total UMIs; falls back to per-1000-genes if total
    UMIs are unavailable, matching the Phase 3 precedent).
  - Positivity rate: fraction of basal cells with raw_HPV16 >= 8 (L-method).
  - Mean load in positive cells (copy-number-in-infected proxy).
  - Per-patient SBS2-HIGH and CNV-HIGH fold enrichment (recomputed).
  - Spearman (n=14) of load vs SBS2 fold, CNV fold, and total HIGH contribution.
  - A spotlight on SC001 / SC027 / SC013 / SC029.

SC008 is the HPV-negative control patient and should read near-zero load (sanity).

Reuses patient_config. Run from scripts/PATIENT_SPECIFIC_EFFECTS/:
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

# =============================================================================
# CONFIG
# =============================================================================
MASTER_HPV_PATH = ("/master/jlehle/WORKING/2026_NMF_PAPER/"
                   "data/FIG_6/01_raw_hpv16_counts/"
                   "basal_cell_master_table_with_raw_HPV16.tsv")
OUT_DIR   = ensure_dir(DIR_00_DIAG)
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")

HPV16_POS_THRESHOLD = 8          # L-method positivity call
FIG5_THRESHOLD = 2.0
SPOTLIGHT = ['Patient SC001', 'Patient SC027', 'Patient SC013', 'Patient SC029']
DPI = 300

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

# tumor vs normal-adjacent source (SBS2/CNV drawn only from tumor basal per Step00B)
NORMAL_SOURCE = 'normal tissue adjucent to head and neck squamous cell carcinoma'
if 'source_name' in basal.obs.columns:
    basal.obs['tissue_grp'] = np.where(
        basal.obs['source_name'].astype(str) == NORMAL_SOURCE, 'normal', 'tumor')
else:
    basal.obs['tissue_grp'] = 'tumor'
n_tumor = int((basal.obs['tissue_grp'] == 'tumor').sum())
rlog(f"  Basal cells: {basal.n_obs:,}")

# raw HPV16 per cell (merge from the HPV master table)
master = pd.read_csv(MASTER_HPV_PATH, sep='\t', index_col=0)
raw_hpv = master['raw_HPV16'].reindex(basal.obs_names)
n_matched = int(raw_hpv.notna().sum())
basal.obs['raw_HPV16'] = raw_hpv.fillna(0).astype(float).values
rlog(f"  raw_HPV16 matched for {n_matched}/{basal.n_obs} basal cells "
     f"({100*n_matched/basal.n_obs:.1f}%); unmatched set to 0")

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
    mean_pos = float(sub.loc[sub['raw_HPV16'] >= HPV16_POS_THRESHOLD, 'raw_HPV16'].mean()) \
        if n_pos > 0 else 0.0
    # contributions: tumor-matched fold (SBS2/CNV seeded only by tumor basal)
    n_tumor_p = int((m & (basal.obs['tissue_grp'] == 'tumor')).sum())
    exp_frac = n_tumor_p / n_tumor if n_tumor else 0
    n_s = int((m & (basal.obs['group'] == 'SBS2_HIGH')).sum())
    n_c = int((m & (basal.obs['group'] == 'CNV_HIGH')).sum())
    fold_s = (n_s / n_sbs2) / exp_frac if exp_frac > 0 and n_sbs2 else 0.0
    fold_c = (n_c / n_cnv) / exp_frac if exp_frac > 0 and n_cnv else 0.0
    rows.append({'patient': p, 'n_basal': n_basal_p,
                 'raw_HPV16_total': raw_sum, 'load_norm': load_norm,
                 'n_pos': n_pos, 'pct_pos': pct_pos, 'mean_load_pos': mean_pos,
                 'n_sbs2': n_s, 'fold_sbs2': fold_s,
                 'n_cnv': n_c, 'fold_cnv': fold_c,
                 'total_HIGH': n_s + n_c})
df = pd.DataFrame(rows).sort_values('load_norm', ascending=False).reset_index(drop=True)
df.to_csv(os.path.join(OUT_DIR, "patient_hpv16_load.tsv"), sep='\t', index=False)

rlog(f"  normalization: {NORM_UNIT}\n")
rlog(f"  {'Patient':<20s} {'Load':>9s} {'%pos':>6s} {'MeanPos':>8s} "
     f"{'SBS2fold':>9s} {'CNVfold':>8s} {'role'}")
rlog(f"  {'-'*20} {'-'*9} {'-'*6} {'-'*8} {'-'*9} {'-'*8} {'-'*4}")
cnv_hc = set(df[df['fold_cnv'] >= FIG5_THRESHOLD]['patient'])
for _, r in df.iterrows():
    role = []
    if r['patient'] in HIGH_CONTRIBUTORS: role.append('SBS2-HC')
    if r['patient'] in cnv_hc:            role.append('CNV-HC')
    rlog(f"  {str(r['patient']):<20s} {r['load_norm']:>9.1f} {r['pct_pos']:>5.1f}% "
         f"{r['mean_load_pos']:>8.1f} {r['fold_sbs2']:>8.2f}x {r['fold_cnv']:>7.2f}x "
         f"  {'+'.join(role)}")

# per-patient source composition, with HPV16 load in each source
# (normal-adjacent load should sit near zero -- a built-in control)
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
    rlog(f"    {str(p).replace('Patient ',''):<8s} tumor={nt:>5d} (load {lt:>7.1f})  "
         f"normal-adj={nn:>4d} (load {ln:>6.1f})")

# =============================================================================
# STEP 2: DOES LOAD PREDICT CONTRIBUTION?  (Spearman, n=14)
# =============================================================================
banner("STEP 2: LOAD vs CONTRIBUTION (Spearman, n=14)")
rlog("  n=14 patients; correlations are low-powered, read as directional.\n")
def sp(a, b, label):
    rho, p = spearmanr(df[a], df[b])
    rlog(f"  {label:<38s} rho={rho:+.3f}  p={p:.3g}")
sp('load_norm', 'total_HIGH', 'load vs total HIGH contribution')
sp('load_norm', 'fold_sbs2',  'load vs SBS2-HIGH fold')
sp('load_norm', 'fold_cnv',   'load vs CNV-HIGH fold')
sp('load_norm', 'pct_pos',    'load vs % HPV16+ (sanity: should be high)')

# =============================================================================
# STEP 3: SPOTLIGHT ON THE CONTRIBUTORS
# =============================================================================
banner("STEP 3: CONTRIBUTOR SPOTLIGHT")
load_rank = {p: i + 1 for i, p in enumerate(df['patient'])}   # 1 = highest load
n_pat = len(df)
for p in SPOTLIGHT:
    if p not in load_rank:
        continue
    r = df[df['patient'] == p].iloc[0]
    role = []
    if p in HIGH_CONTRIBUTORS: role.append('SBS2-HC')
    if p in cnv_hc:            role.append('CNV-HC')
    rlog(f"  {p.replace('Patient ',''):<8s} load={r['load_norm']:.1f} "
         f"(rank {load_rank[p]}/{n_pat}), %pos={r['pct_pos']:.1f}%, "
         f"SBS2 {r['fold_sbs2']:.1f}x / CNV {r['fold_cnv']:.1f}x  [{'+'.join(role)}]")

rlog("")
sc001_rank = load_rank.get('Patient SC001', None)
sc027_rank = load_rank.get('Patient SC027', None)
top_third = max(1, n_pat // 3)
rlog("  READ:")
if sc001_rank and sc001_rank <= top_third:
    rlog(f"    SC001 is a high-load patient (rank {sc001_rank}/{n_pat}), so its dual "
         f"SBS2+CNV contribution is consistent with load seeding both fates.")
else:
    rlog(f"    SC001 is NOT a load outlier (rank {sc001_rank}/{n_pat}); load does not "
         f"explain its dual contribution, so its dual role is genuine.")
if sc027_rank:
    if sc027_rank <= top_third:
        rlog(f"    SC027 (CNV-only driver) is also high-load (rank {sc027_rank}); its "
             f"productive bias sits on top of load rather than being pure load.")
    else:
        rlog(f"    SC027 (CNV-only driver) is NOT high-load (rank {sc027_rank}); its "
             f"CNV dominance is a genuine fate bias, not a load artifact.")
rlog("    If the fate-specific drivers (SC027 CNV-only; SC013/SC029 SBS2-only) are not "
     "load outliers, the divergent-fate reading holds and load explains only SC001.")

# =============================================================================
# STEP 4: PLOT
# =============================================================================
banner("STEP 4: PLOT")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# panel 1: per-patient load bar (sorted desc), labels colored by role
dfp = df.copy()
def role_color(p):
    s = p in HIGH_CONTRIBUTORS; c = p in cnv_hc
    if s and c: return '#7a4fa3'          # both -> purple
    if s:       return '#ed6a5a'          # SBS2 -> coral
    if c:       return '#F6D155'          # CNV -> mustard
    return '#9aa0a6'                      # neither -> gray
colors = [role_color(p) for p in dfp['patient']]
y = np.arange(len(dfp))[::-1]
ax1.barh(y, dfp['load_norm'].values, color=colors, edgecolor='black', linewidth=0.4)
ax1.set_yticks(y)
ax1.set_yticklabels([str(p).replace('Patient ', '') for p in dfp['patient']], fontsize=9)
ax1.set_xlabel(NORM_UNIT, fontsize=11)
ax1.set_title('Per-patient HPV16 load\n(coral=SBS2-HC, mustard=CNV-HC, purple=both)',
              fontsize=11)

# panel 2: load vs total HIGH contribution
ax2.scatter(df['load_norm'], df['total_HIGH'],
            c=[role_color(p) for p in df['patient']], s=70, edgecolor='black', zorder=3)
for _, r in df.iterrows():
    if r['patient'] in SPOTLIGHT or r['patient'] in cnv_hc or r['patient'] in HIGH_CONTRIBUTORS:
        ax2.annotate(str(r['patient']).replace('Patient ', ''),
                     (r['load_norm'], r['total_HIGH']), fontsize=8,
                     xytext=(4, 4), textcoords='offset points')
rho, pv = spearmanr(df['load_norm'], df['total_HIGH'])
ax2.set_xlabel(NORM_UNIT, fontsize=11)
ax2.set_ylabel('HIGH cells contributed (SBS2 + CNV)', fontsize=11)
ax2.set_title(f'Load vs total contribution\n(Spearman rho={rho:+.2f}, p={pv:.2g}, n=14)',
              fontsize=11)
plt.tight_layout()
for ext in ('pdf', 'png'):
    plt.savefig(os.path.join(OUT_DIR, f"patient_hpv16_load.{ext}"),
                dpi=DPI, bbox_inches='tight')
plt.close()
rlog("  [SAVE] patient_hpv16_load.pdf/.png")

report_path = os.path.join(OUT_DIR, f"patient_hpv16_load_{TIMESTAMP}.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
rlog(f"\n  Report: {report_path}")
banner("PATIENT HPV16 LOAD DIAGNOSTIC COMPLETE")
