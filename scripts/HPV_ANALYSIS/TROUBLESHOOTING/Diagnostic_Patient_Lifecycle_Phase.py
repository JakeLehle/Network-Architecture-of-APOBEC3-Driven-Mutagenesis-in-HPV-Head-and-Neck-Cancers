#!/usr/bin/env python3
"""
Diagnostic_Patient_Lifecycle_Phase.py  (v2 -- denominator toggle)
==================================================================
READ-ONLY. Per-patient HPV16 lifecycle-phase balance, and whether it predicts
which HIGH group a patient drives. Patient-level roll-up of Figure 6 Panel F.

Context. The load control showed viral load is the OPPORTUNITY (necessary, not
sufficient): the rank-1 load patient seeds both fates, while other high-load
patients contribute to neither. This asks what sets FATE among a patient's
infected cells: the HPV16 lifecycle phase. Prediction, per patient:
  - maintenance-leaning viral genes (E1/E2)         -> drives SBS2-HIGH
  - productive-leaning viral genes (E4/E5, L1/L2)   -> drives CNV-HIGH
  - phase-intermediate high-load patients           -> drive neither strongly

v2 CHANGES
----------
  - The contribution denominator comes from patient_config.CONTRIBUTION_DENOMINATOR
    via contribution.py. Both folds are computed and written; the setting picks
    the headline column.
  - Contributor sets are derived from the fold columns at runtime rather than
    read from a hardcoded list, so the role tags printed here agree with every
    other script in the chain.
  - The STEP 3 READ block classifies the DERIVED driver sets by their observed
    phase lean instead of asserting in advance which patient leans which way.
  - SPOTLIGHT stays a hardcoded list because it is an ANALYSIS CHOICE (which
    patients get a detailed printout), not a result.

Note on the denominator. This script's folds are not consumed downstream (the
determinants table takes only n_gated, maint_frac and prod_frac from it), but
its printed role tags would contradict the rest of the pipeline if it kept its
own convention. Phase measurement itself is ALWAYS restricted to HPV16-positive
cells regardless of the fold denominator.

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
    mean-of-per-cell-fractions; both estimators are documented.
  - Relate per-patient maintenance balance to SBS2-HIGH fold and productive
    balance to CNV-HIGH fold (Spearman, patients with >= MIN_GATED gated cells).

Inputs (read-only)
------------------
  data/FIG_6/01_raw_hpv16_counts/basal_cell_master_table_with_raw_HPV16.tsv
  data/FIG_6/03_hpv16_genome/per_cell_hpv16_gene_counts.tsv
  data/FIG_4/01_group_selection/three_group_assignments.tsv

Outputs (to data/FIG_5/00_diagnostics/)
---------------------------------------
  patient_lifecycle_phase.tsv           phase fractions + both folds
  patient_lifecycle_phase.pdf/.png
  patient_lifecycle_phase_<ts>.txt

Run from the directory holding patient_config.py and contribution.py:
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
FIG5_THRESHOLD = HC_THRESHOLD         # imported; never redefine locally

PHASES = {'Maintenance': ['E1', 'E2'], 'Amplification': ['E4', 'E5'],
          'Oncogene': ['E6', 'E7'], 'Capsid': ['L1', 'L2']}
CODING = ['E1', 'E2', 'E4', 'E5', 'E6', 'E7', 'L1', 'L2']
PRODUCTIVE = ['E4', 'E5', 'L1', 'L2']
MAINTENANCE = ['E1', 'E2']

# ANALYSIS CHOICE, not a result: which patients get a detailed printout.
SPOTLIGHT = ['Patient SC013', 'Patient SC029', 'Patient SC027',
             'Patient SC010', 'Patient SC022', 'Patient SC001']

PHASE_COLORS = {'Maintenance': '#ed6a5a', 'Amplification': '#f4a259',
                'Oncogene': '#8a8d91', 'Capsid': '#F6D155'}

# Role palette (matches the load diagnostic)
COLOR_BOTH    = '#7a4fa3'
COLOR_SBS2    = '#ed6a5a'
COLOR_CNV     = '#F6D155'
COLOR_NEITHER = '#9aa0a6'

NORMAL_SOURCE = 'normal tissue adjucent to head and neck squamous cell carcinoma'
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

# tumor vs normal-adjacent source
if 'source_name' in master.columns:
    _src = master['source_name'].reindex(df.index)
else:
    import scanpy as sc
    from patient_config import ADATA_PATH
    _src = sc.read_h5ad(ADATA_PATH, backed='r').obs['source_name'].reindex(df.index)
df['tissue_grp'] = np.where(_src.astype(str) == NORMAL_SOURCE, 'normal', 'tumor')

n_basal_total = int(len(df))
n_tumor = int((df['tissue_grp'] == 'tumor').sum())
rlog(f"  tumor basal: {n_tumor}, normal-adjacent: {n_basal_total - n_tumor}")
rlog("")
announce(rlog, n_basal_total, n_tumor)
rlog("")

# HPV16-positive gate (Panel F)
df['hpv_pos'] = (df['raw_HPV16'] >= HPV16_THRESHOLD) & (df[TOTAL_COL] > 0)
rlog(f"  HPV16-positive basal cells (raw>={HPV16_THRESHOLD} & total>0): "
     f"{int(df['hpv_pos'].sum())}")

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
    n_s = int((m & (df['group'] == 'SBS2_HIGH')).sum())
    n_c = int((m & (df['group'] == 'CNV_HIGH')).sum())

    # contributions: BOTH denominators; the headline is set in patient_config
    fs = both_folds(n_s, n_sbs2, n_basal_p, n_basal_total, n_tumor_p, n_tumor)
    fc = both_folds(n_c, n_cnv,  n_basal_p, n_basal_total, n_tumor_p, n_tumor)

    gated = df[m & df['hpv_pos']]
    n_gated = len(gated)
    coding_sum = gated[CODING].sum().sum()
    if n_gated >= 1 and coding_sum > 0:
        maint = gated[MAINTENANCE].sum().sum() / coding_sum   # pooled coding-normalized
        prod  = gated[PRODUCTIVE].sum().sum() / coding_sum
        onco  = gated[['E6', 'E7']].sum().sum() / coding_sum
        ampl  = gated[['E4', 'E5']].sum().sum() / coding_sum
        caps  = gated[['L1', 'L2']].sum().sum() / coding_sum
    else:
        maint = prod = onco = ampl = caps = np.nan

    rows.append({'patient': p, 'n_basal': n_basal_p, 'n_tumor': n_tumor_p,
                 'n_normal_adj': n_basal_p - n_tumor_p, 'n_gated': n_gated,
                 'maint_frac': maint, 'prod_frac': prod, 'onco_frac': onco,
                 'ampl_frac': ampl, 'caps_frac': caps,
                 'n_sbs2': n_s, 'n_cnv': n_c,
                 'fold_sbs2_all_basal': fs['all_basal'],
                 'fold_sbs2_tumor':     fs['tumor'],
                 'fold_cnv_all_basal':  fc['all_basal'],
                 'fold_cnv_tumor':      fc['tumor']})

pt = pd.DataFrame(rows)
pt = attach_folds(pt, 'fold_sbs2')
pt = attach_folds(pt, 'fold_cnv')

# ---- Contributor sets, derived from this run --------------------------------
rlog("")
sbs2_hc = derive_contributors(pt, 'fold_sbs2', expected=HIGH_CONTRIBUTORS,
                              label='SBS2-HIGH contributors', logger=rlog)
rlog("")
cnv_hc = derive_contributors(pt, 'fold_cnv', expected=CNV_HIGH_CONTRIBUTORS,
                             label='CNV-HIGH contributors', logger=rlog)

pt['is_sbs2_driver'] = pt['patient'].isin(sbs2_hc)
pt['is_cnv_driver']  = pt['patient'].isin(cnv_hc)
pt['role'] = ['+'.join([t for t, v in (('SBS2-HC', s), ('CNV-HC', c)) if v]) or '-'
              for s, c in zip(pt['is_sbs2_driver'], pt['is_cnv_driver'])]

pt.to_csv(os.path.join(OUT_DIR, "patient_lifecycle_phase.tsv"), sep='\t', index=False)

rlog(f"\n  pooled coding-normalized phase fractions (URR excluded); "
     f"maintenance=E1/E2, productive=E4/E5+L1/L2")
rlog(f"  folds: headline column = fold_*_{CONTRIBUTION_DENOMINATOR}\n")
rlog(f"  {'Patient':<10s} {'nGate':>6s} {'Maint':>6s} {'Prod':>6s} "
     f"{'SBS2fold':>9s} {'CNVfold':>8s}  role")
rlog(f"  {'-'*10} {'-'*6} {'-'*6} {'-'*6} {'-'*9} {'-'*8}  ----")
for _, r in pt.sort_values('prod_frac', ascending=False).iterrows():
    mstr = f"{r['maint_frac']:.2f}" if pd.notna(r['maint_frac']) else " N.D."
    pstr = f"{r['prod_frac']:.2f}" if pd.notna(r['prod_frac']) else " N.D."
    flag = "  (n<min)" if r['n_gated'] < MIN_GATED else ""
    rlog(f"  {short(r['patient']):<10s} {int(r['n_gated']):>6d} {mstr:>6s} "
         f"{pstr:>6s} {r['fold_sbs2']:>8.2f}x {r['fold_cnv']:>7.2f}x  "
         f"{r['role']}{flag}")

max_gap_s = float((pt['fold_sbs2_all_basal'] - pt['fold_sbs2_tumor']).abs().max())
max_gap_c = float((pt['fold_cnv_all_basal'] - pt['fold_cnv_tumor']).abs().max())
rlog(f"\n  Max |all-basal minus tumor|: SBS2 {max_gap_s:.4f}x, CNV {max_gap_c:.4f}x")

# per-patient source composition, with HPV16-positive counts in each source
rlog("\n  Per-patient basal by source (tumor / normal-adjacent), HPV16+ in each:")
for p in patients:
    pm = df[PATIENT_COL] == p
    tm = pm & (df['tissue_grp'] == 'tumor')
    nm = pm & (df['tissue_grp'] == 'normal')
    nt, nn = int(tm.sum()), int(nm.sum())
    pos_t = int((tm & df['hpv_pos']).sum())
    pos_n = int((nm & df['hpv_pos']).sum())
    rlog(f"    {short(p):<8s} tumor={nt:>5d} (HPV+ {pos_t:>4d})  "
         f"normal-adj={nn:>4d} (HPV+ {pos_n:>3d})")

# =============================================================================
# STEP 2: LOOP CLOSURE (Spearman, patients with >= MIN_GATED cells)
# =============================================================================
banner("STEP 2: PHASE-TO-CONTRIBUTION LOOP (Spearman)")
ok = pt[pt['n_gated'] >= MIN_GATED].dropna(subset=['maint_frac', 'prod_frac'])
rlog(f"  patients with >= {MIN_GATED} gated cells: {len(ok)} of {len(pt)} "
     f"(low power, read directionally)\n")
if len(ok) >= 4:
    r1, p1 = spearmanr(ok['maint_frac'], ok['fold_sbs2'])
    r2, p2 = spearmanr(ok['prod_frac'], ok['fold_cnv'])
    r3, p3 = spearmanr(ok['prod_frac'], ok['fold_sbs2'])
    r4, p4 = spearmanr(ok['maint_frac'], ok['fold_cnv'])
    rlog(f"  maintenance balance vs SBS2-HIGH fold : rho={r1:+.3f}  p={p1:.3g}  "
         f"(predicted +)")
    rlog(f"  productive  balance vs CNV-HIGH  fold : rho={r2:+.3f}  p={p2:.3g}  "
         f"(predicted +)")
    rlog("  cross-checks (predicted weak/negative):")
    rlog(f"    productive  vs SBS2-HIGH fold        : rho={r3:+.3f}  p={p3:.3g}")
    rlog(f"    maintenance vs CNV-HIGH  fold        : rho={r4:+.3f}  p={p4:.3g}")
    rlog("")
    n_as_predicted = sum([r1 > 0, r2 > 0, r3 < 0, r4 < 0])
    rlog(f"  {n_as_predicted} of 4 signs match the prediction. At n={len(ok)} the")
    rlog("  sign structure is the informative part; the p-values are underpowered")
    rlog("  and should not be presented as the result.")
else:
    rlog("  too few well-covered patients for a correlation.")

# =============================================================================
# STEP 3: SPOTLIGHT + DERIVED READ
# =============================================================================
banner("STEP 3: SPOTLIGHT")
for p in SPOTLIGHT:
    r = pt[pt['patient'] == p]
    if r.empty:
        continue
    r = r.iloc[0]
    if pd.isna(r['maint_frac']):
        rlog(f"  {short(p):<8s} n_gated={int(r['n_gated'])} -> no phase profile")
        continue
    lean = 'maintenance' if r['maint_frac'] > r['prod_frac'] else 'productive'
    rlog(f"  {short(p):<8s} maint={r['maint_frac']:.2f} prod={r['prod_frac']:.2f} "
         f"-> {lean}-leaning  | SBS2 {r['fold_sbs2']:.1f}x / "
         f"CNV {r['fold_cnv']:.1f}x [{r['role']}] (n={int(r['n_gated'])})")

# ---- READ, derived from the driver sets -------------------------------------
rlog("")
rlog("  READ (derived; lean = whichever pooled balance is larger, among patients")
rlog(f"  with >= {MIN_GATED} gated cells):")

measurable = pt[(pt['n_gated'] >= MIN_GATED) & pt['maint_frac'].notna()]


def lean_of(row):
    return 'maintenance' if row['maint_frac'] > row['prod_frac'] else 'productive'


agree = disagree = 0
for _, r in measurable.iterrows():
    if not (r['is_sbs2_driver'] or r['is_cnv_driver']):
        continue
    lean = lean_of(r)
    both = r['is_sbs2_driver'] and r['is_cnv_driver']
    expected = ('maintenance' if r['is_sbs2_driver'] else 'productive')
    if both:
        rlog(f"    {short(r['patient'])} drives BOTH and leans {lean} "
             f"({r['maint_frac']:.2f} vs {r['prod_frac']:.2f}); a phase-based rule")
        rlog(f"      cannot assign it to one fate, so it must reach the other by a "
             f"different route.")
        continue
    if lean == expected:
        agree += 1
        rlog(f"    {short(r['patient'])} drives "
             f"{'SBS2' if r['is_sbs2_driver'] else 'CNV'} and leans {lean} "
             f"({r['maint_frac']:.2f} vs {r['prod_frac']:.2f}) AS PREDICTED.")
    else:
        disagree += 1
        rlog(f"    {short(r['patient'])} drives "
             f"{'SBS2' if r['is_sbs2_driver'] else 'CNV'} but leans {lean} "
             f"({r['maint_frac']:.2f} vs {r['prod_frac']:.2f}) AGAINST prediction.")

rlog("")
rlog(f"    {agree} single-fate drivers lean as predicted, {disagree} against.")

non_driver_leans = [(short(r['patient']), lean_of(r),
                     abs(r['maint_frac'] - r['prod_frac']))
                    for _, r in measurable.iterrows()
                    if not (r['is_sbs2_driver'] or r['is_cnv_driver'])]
if non_driver_leans:
    intermediate = sorted([n for n, _, gap in non_driver_leans if gap < 0.25])
    if intermediate:
        rlog(f"    Phase-intermediate non-drivers (|maint - prod| < 0.25): "
             f"{intermediate}.")
        rlog("    These are the patients where virus is present without a clear")
        rlog("    phase bias, which is why load alone does not produce a fate.")

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
axes[0].set_yticklabels([short(p) for p in plot_pt['patient']], fontsize=9)
axes[0].set_xlabel('coding-normalized phase fraction', fontsize=11)
axes[0].set_title('Per-patient HPV16 phase composition', fontsize=11)
axes[0].legend(fontsize=8, ncol=2, loc='lower right', frameon=False)


def rc(p):
    s, c = p in sbs2_hc, p in cnv_hc
    if s and c:
        return COLOR_BOTH
    if s:
        return COLOR_SBS2
    if c:
        return COLOR_CNV
    return COLOR_NEITHER


for ax, (xcol, ycol, xlab, ylab, ttl) in zip(
        axes[1:], [('maint_frac', 'fold_sbs2', 'maintenance balance',
                    'SBS2-HIGH fold', 'Maintenance -> SBS2'),
                   ('prod_frac', 'fold_cnv', 'productive balance',
                    'CNV-HIGH fold', 'Productive -> CNV')]):
    ax.scatter(plot_pt[xcol], plot_pt[ycol],
               c=[rc(p) for p in plot_pt['patient']], s=70,
               edgecolor='black', zorder=3)
    for _, r in plot_pt.iterrows():
        ax.annotate(short(r['patient']), (r[xcol], r[ycol]),
                    fontsize=8, xytext=(4, 3), textcoords='offset points')
    ax.axhline(HC_THRESHOLD, color='#888', linestyle=':', linewidth=1)
    if len(plot_pt) >= 4:
        rho, pv = spearmanr(plot_pt[xcol], plot_pt[ycol])
        ttl = f"{ttl}\n(rho={rho:+.2f}, p={pv:.2g}, n={len(plot_pt)})"
    ax.set_xlabel(xlab, fontsize=11)
    ax.set_ylabel(f"{ylab} ({CONTRIBUTION_DENOMINATOR})", fontsize=11)
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
rlog(f"\n  Denominator: {CONTRIBUTION_DENOMINATOR}")
rlog(f"  Table:  {os.path.join(OUT_DIR, 'patient_lifecycle_phase.tsv')}")
rlog(f"  Report: {report_path}")
banner("PATIENT LIFECYCLE PHASE DIAGNOSTIC COMPLETE (v2)")
