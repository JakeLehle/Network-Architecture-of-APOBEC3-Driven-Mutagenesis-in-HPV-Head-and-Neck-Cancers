#!/usr/bin/env python3
"""
Diagnostic_Patient_Lifecycle_Phase.py
=====================================
READ-ONLY. Per-patient HPV16 lifecycle-phase balance, and whether it predicts
which HIGH group a patient drives. Patient-level roll-up of Figure 6 Panel F.

Context. The load control showed viral load is the OPPORTUNITY (necessary, not
sufficient): SC001 (rank-1 load) seeds both fates, but high-load SC010/SC022
contribute to neither. This asks what sets FATE among a patient's infected cells:
the HPV16 lifecycle phase. Prediction, per patient:
  - maintenance-leaning viral genes (E1/E2)         -> drives SBS2-HIGH
  - productive-leaning viral genes (E4/E5, L1/L2)   -> drives CNV-HIGH
  - phase-intermediate (e.g. high-load SC010/SC022) -> drives neither strongly

Method
------
  - Gate to HPV16-positive basal cells (raw_HPV16 >= 8 AND total genome reads > 0),
    the Panel F gate.
  - Phase reads per manuscript: maintenance = E1+E2; productive = E4+E5+L1+L2
    (amplification + capsid); oncogene = E6+E7 (small).
  - Per patient, POOLED coding-normalized balance: sum(phase reads) / sum(coding
    reads) across the patient's gated cells, coding = the 8 ORFs, URR excluded.
    Pooling is depth-weighted (stable per patient); URR exclusion sharpens the
    maintenance-vs-productive contrast. This DIFFERS from Panel F's
    mean-of-per-cell-fractions; both are printed for traceability.
  - Relate per-patient maintenance balance to SBS2-HIGH fold and productive
    balance to CNV-HIGH fold (Spearman, patients with >= 10 gated cells).

Reuses patient_config paths/helpers. Run from scripts/PATIENT_SPECIFIC_EFFECTS/:
    conda run -n NETWORK python Diagnostic_Patient_Lifecycle_Phase.py

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
    THREE_GROUP_PATH, DIR_00_DIAG, HIGH_CONTRIBUTORS,
    banner, log, ensure_dir,
)

# =============================================================================
# CONFIG
# =============================================================================
ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
MASTER_HPV_PATH = os.path.join(ROOT, "data/FIG_6/01_raw_hpv16_counts/"
                                     "basal_cell_master_table_with_raw_HPV16.tsv")
HPV_GENE_PATH   = os.path.join(ROOT, "data/FIG_6/03_hpv16_genome/"
                                     "per_cell_hpv16_gene_counts.tsv")
OUT_DIR   = ensure_dir(DIR_00_DIAG)
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")

HPV16_THRESHOLD = 8
TOTAL_COL = 'total_hpv16_genome_reads'
PATIENT_COL = 'subject id'            # column name in the HPV master table
MIN_GATED = 10                        # min HPV16+ cells for a stable phase estimate
FIG5_THRESHOLD = 2.0

PHASES = {'Maintenance': ['E1', 'E2'], 'Amplification': ['E4', 'E5'],
          'Oncogene': ['E6', 'E7'], 'Capsid': ['L1', 'L2']}
CODING = ['E1', 'E2', 'E4', 'E5', 'E6', 'E7', 'L1', 'L2']
PRODUCTIVE = ['E4', 'E5', 'L1', 'L2']
MAINTENANCE = ['E1', 'E2']

SPOTLIGHT = ['Patient SC013', 'Patient SC029', 'Patient SC027',
             'Patient SC010', 'Patient SC022', 'Patient SC001']
PHASE_COLORS = {'Maintenance': '#ed6a5a', 'Amplification': '#f4a259',
                'Oncogene': '#8a8d91', 'Capsid': '#F6D155'}
DPI = 300

report_lines = []
def rlog(msg=""):
    log(msg)
    report_lines.append(str(msg))

# =============================================================================
# STEP 0: LOAD + MERGE
# =============================================================================
banner("STEP 0: LOAD DATA")
master = pd.read_csv(MASTER_HPV_PATH, sep='\t', index_col=0)
hpv = pd.read_csv(HPV_GENE_PATH, sep='\t', index_col=0)
groups = pd.read_csv(THREE_GROUP_PATH, sep='\t')
rlog(f"  master basal cells: {len(master):,}")
rlog(f"  HPV16 gene table:   {len(hpv):,}")

# spine = basal master (has patient + raw_HPV16); attach phase reads + group
df = master[[PATIENT_COL, 'raw_HPV16']].copy()
for g in CODING + ['URR', 'intergenic', TOTAL_COL]:
    df[g] = hpv[g].reindex(df.index).fillna(0) if g in hpv.columns else 0
grp_map = dict(zip(groups['cell_barcode'], groups['group']))
df['group'] = df.index.map(lambda b: grp_map.get(b, 'other'))

n_sbs2 = int((df['group'] == 'SBS2_HIGH').sum())
n_cnv  = int((df['group'] == 'CNV_HIGH').sum())
rlog(f"  groups on basal: SBS2_HIGH={n_sbs2}, CNV_HIGH={n_cnv}, "
     f"NORMAL={(df['group']=='NORMAL').sum()}")

# tumor vs normal-adjacent source (SBS2/CNV drawn only from tumor basal per Step00B)
NORMAL_SOURCE = 'normal tissue adjucent to head and neck squamous cell carcinoma'
if 'source_name' in master.columns:
    _src = master['source_name'].reindex(df.index)
else:
    import scanpy as sc
    from patient_config import ADATA_PATH
    _src = sc.read_h5ad(ADATA_PATH, backed='r').obs['source_name'].reindex(df.index)
df['tissue_grp'] = np.where(_src.astype(str) == NORMAL_SOURCE, 'normal', 'tumor')
n_tumor = int((df['tissue_grp'] == 'tumor').sum())
rlog(f"  tumor basal: {n_tumor}, normal-adjacent: {len(df) - n_tumor}")

# HPV16-positive gate (Panel F)
df['hpv_pos'] = (df['raw_HPV16'] >= HPV16_THRESHOLD) & (df[TOTAL_COL] > 0)
rlog(f"  HPV16-positive basal cells (raw>=8 & total>0): {int(df['hpv_pos'].sum())}")

patients = sorted(df[PATIENT_COL].dropna().unique())
rlog(f"  Patients: {len(patients)}")

# =============================================================================
# STEP 1: PER-PATIENT PHASE BALANCE + CONTRIBUTION
# =============================================================================
banner("STEP 1: PER-PATIENT LIFECYCLE PHASE BALANCE")
rows = []
for p in patients:
    m = df[PATIENT_COL] == p
    n_basal_p = int(m.sum())
    n_tumor_p = int((m & (df['tissue_grp'] == 'tumor')).sum())
    exp_frac = n_tumor_p / n_tumor if n_tumor else 0   # tumor-matched
    n_s = int((m & (df['group'] == 'SBS2_HIGH')).sum())
    n_c = int((m & (df['group'] == 'CNV_HIGH')).sum())
    fold_s = (n_s / n_sbs2) / exp_frac if exp_frac > 0 and n_sbs2 else 0.0
    fold_c = (n_c / n_cnv) / exp_frac if exp_frac > 0 and n_cnv else 0.0

    gated = df[m & df['hpv_pos']]
    n_gated = len(gated)
    coding_sum = gated[CODING].sum().sum()
    if n_gated >= 1 and coding_sum > 0:
        maint = gated[MAINTENANCE].sum().sum() / coding_sum       # pooled coding-normalized
        prod  = gated[PRODUCTIVE].sum().sum() / coding_sum
        onco  = gated[['E6', 'E7']].sum().sum() / coding_sum
        ampl  = gated[['E4', 'E5']].sum().sum() / coding_sum
        caps  = gated[['L1', 'L2']].sum().sum() / coding_sum
    else:
        maint = prod = onco = ampl = caps = np.nan
    rows.append({'patient': p, 'n_basal': n_basal_p, 'n_gated': n_gated,
                 'maint_frac': maint, 'prod_frac': prod, 'onco_frac': onco,
                 'ampl_frac': ampl, 'caps_frac': caps,
                 'n_sbs2': n_s, 'fold_sbs2': fold_s,
                 'n_cnv': n_c, 'fold_cnv': fold_c})
pt = pd.DataFrame(rows)
pt.to_csv(os.path.join(OUT_DIR, "patient_lifecycle_phase.tsv"), sep='\t', index=False)

cnv_hc = set(pt[pt['fold_cnv'] >= FIG5_THRESHOLD]['patient'])
rlog(f"  pooled coding-normalized phase fractions (URR excluded); "
     f"maintenance=E1/E2, productive=E4/E5+L1/L2\n")
rlog(f"  {'Patient':<20s} {'nGate':>6s} {'Maint':>6s} {'Prod':>6s} "
     f"{'SBS2fold':>9s} {'CNVfold':>8s} {'role'}")
rlog(f"  {'-'*20} {'-'*6} {'-'*6} {'-'*6} {'-'*9} {'-'*8}")
for _, r in pt.sort_values('prod_frac', ascending=False).iterrows():
    role = []
    if r['patient'] in HIGH_CONTRIBUTORS: role.append('SBS2-HC')
    if r['patient'] in cnv_hc:            role.append('CNV-HC')
    mstr = f"{r['maint_frac']:.2f}" if pd.notna(r['maint_frac']) else " N.D."
    pstr = f"{r['prod_frac']:.2f}" if pd.notna(r['prod_frac']) else " N.D."
    flag = "  (n<min)" if r['n_gated'] < MIN_GATED else ""
    rlog(f"  {str(r['patient']):<20s} {int(r['n_gated']):>6d} {mstr:>6s} {pstr:>6s} "
         f"{r['fold_sbs2']:>8.2f}x {r['fold_cnv']:>7.2f}x  {'+'.join(role)}{flag}")

# per-patient source composition, with HPV16-positive counts in each source
rlog("\n  Per-patient basal by source (tumor / normal-adjacent), HPV16+ in each:")
for p in patients:
    pm = df[PATIENT_COL] == p
    tm = pm & (df['tissue_grp'] == 'tumor')
    nm = pm & (df['tissue_grp'] == 'normal')
    nt, nn = int(tm.sum()), int(nm.sum())
    pos_t = int((tm & df['hpv_pos']).sum())
    pos_n = int((nm & df['hpv_pos']).sum())
    rlog(f"    {str(p).replace('Patient ',''):<8s} tumor={nt:>5d} (HPV+ {pos_t:>4d})  "
         f"normal-adj={nn:>4d} (HPV+ {pos_n:>3d})")

# =============================================================================
# STEP 2: LOOP CLOSURE (Spearman, patients with >= MIN_GATED cells)
# =============================================================================
banner("STEP 2: PHASE-TO-CONTRIBUTION LOOP (Spearman)")
ok = pt[pt['n_gated'] >= MIN_GATED].dropna(subset=['maint_frac', 'prod_frac'])
rlog(f"  patients with >= {MIN_GATED} gated cells: {len(ok)} of {len(pt)} "
     f"(n=14 overall; low power, read directionally)\n")
if len(ok) >= 4:
    r1, p1 = spearmanr(ok['maint_frac'], ok['fold_sbs2'])
    r2, p2 = spearmanr(ok['prod_frac'], ok['fold_cnv'])
    r3, p3 = spearmanr(ok['prod_frac'], ok['fold_sbs2'])
    r4, p4 = spearmanr(ok['maint_frac'], ok['fold_cnv'])
    rlog(f"  maintenance balance vs SBS2-HIGH fold : rho={r1:+.3f}  p={p1:.3g}  (predicted +)")
    rlog(f"  productive  balance vs CNV-HIGH  fold : rho={r2:+.3f}  p={p2:.3g}  (predicted +)")
    rlog(f"  cross-checks (predicted weak/negative):")
    rlog(f"    productive  vs SBS2-HIGH fold        : rho={r3:+.3f}  p={p3:.3g}")
    rlog(f"    maintenance vs CNV-HIGH  fold        : rho={r4:+.3f}  p={p4:.3g}")
else:
    rlog("  too few well-covered patients for a correlation.")

# =============================================================================
# STEP 3: SPOTLIGHT
# =============================================================================
banner("STEP 3: SPOTLIGHT")
for p in SPOTLIGHT:
    r = pt[pt['patient'] == p]
    if r.empty:
        continue
    r = r.iloc[0]
    if pd.isna(r['maint_frac']):
        rlog(f"  {p.replace('Patient ',''):<8s} n_gated={int(r['n_gated'])} -> no phase profile")
        continue
    lean = 'maintenance' if r['maint_frac'] > r['prod_frac'] else 'productive'
    role = []
    if p in HIGH_CONTRIBUTORS: role.append('SBS2-HC')
    if p in cnv_hc:            role.append('CNV-HC')
    rlog(f"  {p.replace('Patient ',''):<8s} maint={r['maint_frac']:.2f} prod={r['prod_frac']:.2f} "
         f"-> {lean}-leaning  | SBS2 {r['fold_sbs2']:.1f}x / CNV {r['fold_cnv']:.1f}x "
         f"[{'+'.join(role) or 'none'}] (n={int(r['n_gated'])})")
rlog("\n  READ: SBS2 drivers (SC013/SC029) should lean maintenance; the CNV driver")
rlog("  (SC027) should lean productive; high-load non-contributors (SC010/SC022)")
rlog("  should sit intermediate, explaining why load without a phase bias -> no fate.")

# =============================================================================
# STEP 4: PLOT (3 panels)
# =============================================================================
banner("STEP 4: PLOT")
plot_pt = pt[pt['n_gated'] >= MIN_GATED].dropna(subset=['maint_frac']).copy()
plot_pt = plot_pt.sort_values('prod_frac')
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# panel 1: stacked coding-normalized phase composition per patient
y = np.arange(len(plot_pt))
left = np.zeros(len(plot_pt))
for phase, col in [('Maintenance', 'maint_frac'), ('Amplification', 'ampl_frac'),
                   ('Oncogene', 'onco_frac'), ('Capsid', 'caps_frac')]:
    vals = plot_pt[col].values
    axes[0].barh(y, vals, left=left, color=PHASE_COLORS[phase], label=phase,
                 edgecolor='white', linewidth=0.4)
    left += vals
axes[0].set_yticks(y)
axes[0].set_yticklabels([str(p).replace('Patient ', '') for p in plot_pt['patient']], fontsize=9)
axes[0].set_xlabel('coding-normalized phase fraction', fontsize=11)
axes[0].set_title('Per-patient HPV16 phase composition', fontsize=11)
axes[0].legend(fontsize=8, ncol=2, loc='lower right', frameon=False)

# panel 2: maintenance vs SBS2 fold
def rc(p):
    s = p in HIGH_CONTRIBUTORS; c = p in cnv_hc
    return '#7a4fa3' if s and c else '#ed6a5a' if s else '#F6D155' if c else '#9aa0a6'
for ax, (xcol, ycol, xlab, ylab, ttl) in zip(
        axes[1:], [('maint_frac', 'fold_sbs2', 'maintenance balance', 'SBS2-HIGH fold',
                    'Maintenance -> SBS2'),
                   ('prod_frac', 'fold_cnv', 'productive balance', 'CNV-HIGH fold',
                    'Productive -> CNV')]):
    ax.scatter(plot_pt[xcol], plot_pt[ycol],
               c=[rc(p) for p in plot_pt['patient']], s=70, edgecolor='black', zorder=3)
    for _, r in plot_pt.iterrows():
        ax.annotate(str(r['patient']).replace('Patient ', ''), (r[xcol], r[ycol]),
                    fontsize=8, xytext=(4, 3), textcoords='offset points')
    ax.axhline(FIG5_THRESHOLD, color='#888', linestyle=':', linewidth=1)
    if len(plot_pt) >= 4:
        rho, pv = spearmanr(plot_pt[xcol], plot_pt[ycol])
        ttl = f"{ttl}\n(rho={rho:+.2f}, p={pv:.2g})"
    ax.set_xlabel(xlab, fontsize=11); ax.set_ylabel(ylab, fontsize=11)
    ax.set_title(ttl, fontsize=11)
plt.tight_layout()
for ext in ('pdf', 'png'):
    plt.savefig(os.path.join(OUT_DIR, f"patient_lifecycle_phase.{ext}"),
                dpi=DPI, bbox_inches='tight')
plt.close()
rlog("  [SAVE] patient_lifecycle_phase.pdf/.png")

report_path = os.path.join(OUT_DIR, f"patient_lifecycle_phase_{TIMESTAMP}.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
rlog(f"\n  Report: {report_path}")
banner("PATIENT LIFECYCLE PHASE DIAGNOSTIC COMPLETE")
