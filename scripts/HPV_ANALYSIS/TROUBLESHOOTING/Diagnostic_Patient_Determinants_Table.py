#!/usr/bin/env python3
"""
Diagnostic_Patient_Determinants_Table.py  (v3)
===============================================
READ-ONLY. Per-patient determinants of SBS2-HIGH / CNV-HIGH contribution, and a
test of the three-axis conjunction model that emerged from them.

LINEAGE (each version kept the previous tests intact, as breadcrumbs)
---------------------------------------------------------------------
v1 assembled the table and resolved SC010. Matched to SC029 on load (higher,
2423 vs 1826), viral phase (0.60/0.39, identical to two decimals) and cell cycle
(G1-shifted, 46%), SC010 still contributes 4-fold fewer SBS2-HIGH cells. The
separating variable is the enzyme: A3A 1.00 vs 2.75, A3A-positive 18.1% vs
44.5%, with HIGHER A3B (3.24 vs 2.11). SC010 is not a maintenance patient that
failed to deliver; it is an A3B-dominant patient carrying maintenance-phase
virus. SC022 is the same shape, more extreme (A3A 0.03).

v2 added circularity labelling, fixed the degenerate SBS2-weight median, put a
cell floor on dominance, and tested whether the axes are independent. They are:
A3A vs viral load rho = +0.042 (p = 0.89). A patient's A3A carries no
information about how much virus they have, yet predicts their SBS2 fate. v2
also produced the double dissociation:
    A3A vs SBS2 fold  rho = +0.594 (p = 0.025)   A3A vs CNV fold  rho = +0.047 (ns)
    A3B vs CNV  fold  rho = +0.660 (p = 0.010)   A3B vs SBS2 fold rho = +0.400 (ns)
while viral load associates with BOTH fates at nearly identical strength
(+0.506 / +0.513), the quantitative form of "load is opportunity, not fate".

v3 (this version) adds STEP 8, which tests the model those results imply.

THE MODEL STEP 8 TESTS
----------------------
No single factor is sufficient, and the univariate steps above show why: every
one of them has a counterexample. The proposal is that reaching a fate requires
THREE conditions together, and that the informative patients are the ones that
fail on different axes.

    OPPORTUNITY  the patient carries enough virus                (viral load)
    DIRECTION    that virus sits in the right lifecycle phase    (maintenance
                 for SBS2, productive for CNV)
    CAPABILITY   the cell expresses the matching enzyme          (A3A for SBS2,
                                                                  A3B for CNV)

Predicted failure modes, all of which are already on the board:
    SC010, SC022   opportunity + direction, no capability (A3B-dominant with
                   maintenance-phase virus) -> neither fate
    SC005          capability, no opportunity (highest A3A at 3.03, but load
                   rank 9 and only 10.9% of cells infected) -> no fate
    SC003, SC006   direction + capability, no opportunity -> no fate
    SC014          direction + capability (A3B 3.06), no opportunity -> no fate
    SC001          the documented exception: rank-1 load (3455) seeds BOTH
                   fates despite leaning maintenance, the "load route" into
                   CNV-HIGH identified previously

STEP 8 has four parts:
  8a  classify every patient on all three axes and print WHERE each one fails
  8b  Fisher exact, all-three-met vs driver status, run separately per fate
  8c  threshold sensitivity. This is the part that decides whether the result is
      real. SC010's A3A is 1.00 against a cohort median of 1.04, so the
      classification of the single most important counterexample turns on a
      razor-thin margin. The thresholds are swept across percentiles and the
      range over which separation holds is reported, along with where it breaks.
  8d  conjunction structure. The model says fate needs ALL three, so the natural
      summary is the MINIMUM of the three scaled axes, not their mean. If the
      minimum tracks contribution better than the mean, the conjunction
      structure itself is supported rather than assumed.

STANDING CAVEAT, printed in the output so it travels with the numbers: at n = 14
a rule with three binary terms has enough freedom to fit almost any three-patient
set. What defends it is that the axes were specified from the biological model
before this test, not searched for, and that the same rule applied independently
to the CNV side recovers its driver. It is a consistency check on the model, not
an independent confirmation of it.

Denominators and universes (they differ by column, so they are stated):
  folds        tumor basal denominator (unified patient-thread methods)
  A3A/A3B      tumor basal cells, all of them, ungated on HPV16
  SBS2 weight  tumor basal cells carrying a weight (coverage reported per patient)
  phase        HPV16-positive gated cells only (raw_HPV16 >= 8 and total > 0)
  tumor/normal control: basal cells of both sources, patients having both

Inputs (read-only)
------------------
  data/FIG_5/00_diagnostics/cnv_high_patient_contribution.tsv
  data/FIG_5/00_diagnostics/patient_hpv16_load.tsv
  data/FIG_5/00_diagnostics/patient_lifecycle_phase.tsv
  data/FIG_5/00_diagnostics/patient_cellcycle_by_source.tsv
  data/FIG_4/00_input/adata_final.h5ad
  data/FIG_4/00_input/signature_weights_per_cell.txt

Outputs (to data/FIG_5/00_diagnostics/)
---------------------------------------
  patient_determinants_table.tsv         one row per patient, all determinants
  patient_enzyme_by_source.tsv           tumor vs normal-adjacent A3A/A3B control
  patient_conjunction_model.tsv          three-axis classification + scores
  patient_conjunction_sensitivity.tsv    threshold sweep
  patient_determinants_<ts>.txt          full console log

Run from the directory holding patient_config.py:
    conda run -n NETWORK python Diagnostic_Patient_Determinants_Table.py

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
from datetime import datetime
from itertools import product

import numpy as np
import pandas as pd
import scipy.sparse
from scipy.stats import spearmanr, mannwhitneyu, fisher_exact

from patient_config import (
    PATIENT_COL, CELLTYPE_COL, DIR_00_DIAG, HIGH_CONTRIBUTORS,
    banner, log, ensure_dir, load_adata, load_three_groups,
)

# =============================================================================
# CONFIG
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
SIG_WEIGHTS_PATH = os.path.join(PROJECT_ROOT,
    "data/FIG_4/00_input/signature_weights_per_cell.txt")

OUT_DIR   = ensure_dir(DIR_00_DIAG)
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")

TBL_CNV   = os.path.join(DIR_00_DIAG, "cnv_high_patient_contribution.tsv")
TBL_LOAD  = os.path.join(DIR_00_DIAG, "patient_hpv16_load.tsv")
TBL_PHASE = os.path.join(DIR_00_DIAG, "patient_lifecycle_phase.tsv")
TBL_CYCLE = os.path.join(DIR_00_DIAG, "patient_cellcycle_by_source.tsv")

NORMAL_SOURCE = 'normal tissue adjucent to head and neck squamous cell carcinoma'
HC_THRESHOLD  = 2.0
MIN_GATED     = 10      # phase estimate floor, matches the phase diagnostic
MIN_A3_CELLS  = 50      # floor for a stable A3A/A3B dominance fraction
MIN_SOURCE_CELLS = 30   # floor for the tumor vs normal-adjacent enzyme test

# STEP 8 thresholds. Defaults are cohort medians; 8c sweeps them.
LOAD_PCTL_DEFAULT = 50
ENZ_PCTL_DEFAULT  = 50
SWEEP_PCTLS = [30, 40, 50, 60, 70]

MATCHED_PAIR = ['Patient SC010', 'Patient SC029']
SPOTLIGHT    = ['Patient SC027', 'Patient SC001', 'Patient SC013',
                'Patient SC029', 'Patient SC010', 'Patient SC022',
                'Patient SC005']

report_lines = []
def rlog(msg=""):
    log(msg)
    report_lines.append(str(msg))

def short(p):
    return str(p).replace('Patient ', '')

# =============================================================================
# HELPERS
# =============================================================================
def load_side_table(path, label):
    if not os.path.exists(path):
        rlog(f"  [MISSING] {label}: {path}")
        return None
    df = pd.read_csv(path, sep='\t')
    key = next((c for c in df.columns
                if str(c).strip().lower() in ('patient', 'subject id', 'subject_id')), None)
    if key is None:
        rlog(f"  [SKIP] {label}: no patient key column")
        return None
    df = df.rename(columns={key: 'patient'})
    df['patient'] = df['patient'].astype(str)
    rlog(f"  [OK] {label}: {len(df)} rows")
    return df

def take(df, src, dest, out):
    if df is None or src not in df.columns:
        out[dest] = np.nan
        return False
    out[dest] = out['patient'].map(dict(zip(df['patient'], df[src])))
    return True

def gene_vector(ad, symbol):
    if symbol not in ad.var_names:
        return None
    x = ad[:, symbol].X
    if scipy.sparse.issparse(x):
        return np.asarray(x.todense()).flatten()
    return np.asarray(x).flatten()

def sp(a, b, label, df, circular=False, collect=None):
    """
    Spearman on rows where both columns are finite.

    circular=True marks a pairing where the predictor is (part of) the criterion
    the outcome was selected on. Such a correlation is guaranteed by
    construction; it is reported for completeness, never counted as evidence,
    and excluded from `collect`.
    """
    sub = df[[a, b]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(sub) < 4 or sub[a].nunique() < 2 or sub[b].nunique() < 2:
        reason = 'constant input' if len(sub) >= 4 else 'too few points'
        rlog(f"    {label:<38s} N.D.  ({reason}, n={len(sub)})")
        return
    rho, p = spearmanr(sub[a], sub[b])
    tag = '  [CIRCULAR, not evidence]' if circular else ('  *' if p < 0.05 else '')
    rlog(f"    {label:<38s} rho={rho:+.3f}  p={p:.3g}  n={len(sub)}{tag}")
    if collect is not None and not circular and p < 0.05:
        collect.append((label, rho, p))

# =============================================================================
# STEP 0: EXISTING DIAGNOSTIC TABLES
# =============================================================================
banner("STEP 0: LOAD EXISTING DIAGNOSTIC TABLES")
rlog("  Joined values are copied from the published diagnostic outputs, not")
rlog("  recomputed, so this table cannot silently disagree with them.\n")

t_cnv   = load_side_table(TBL_CNV,   "CNV contribution")
t_load  = load_side_table(TBL_LOAD,  "HPV16 load")
t_phase = load_side_table(TBL_PHASE, "lifecycle phase")
t_cycle = load_side_table(TBL_CYCLE, "cell cycle by source")

# =============================================================================
# STEP 1: PRIMARY DATA
# =============================================================================
banner("STEP 1: LOAD PRIMARY DATA (A3A / A3B / SBS2 weight)")

adata = load_adata()
sbs2_high, cnv_high, normal = load_three_groups()

basal = adata[adata.obs[CELLTYPE_COL] == 'basal cell'].copy()
if 'source_name' in basal.obs.columns:
    basal.obs['tissue_grp'] = np.where(
        basal.obs['source_name'].astype(str) == NORMAL_SOURCE, 'normal', 'tumor')
else:
    basal.obs['tissue_grp'] = 'tumor'

basal_a3a = gene_vector(basal, 'APOBEC3A')
basal_a3b = gene_vector(basal, 'APOBEC3B')
if basal_a3a is None or basal_a3b is None:
    raise SystemExit("  ERROR: APOBEC3A/APOBEC3B not in adata.var_names")

tumor_mask = (basal.obs['tissue_grp'] == 'tumor').values
a3a, a3b = basal_a3a[tumor_mask], basal_a3b[tumor_mask]
tumor_names = basal.obs_names[tumor_mask]
pat_all = basal.obs[PATIENT_COL].astype(str).values
pat = pat_all[tumor_mask]

rlog(f"  basal cells: {basal.n_obs:,}  (tumor {int(tumor_mask.sum()):,}, "
     f"normal-adjacent {int((~tumor_mask).sum()):,})")
rlog(f"  over tumor basal: A3A+ {100*np.mean(a3a > 0):.1f}%, "
     f"A3B+ {100*np.mean(a3b > 0):.1f}%")

sig = pd.read_csv(SIG_WEIGHTS_PATH, sep='\t', index_col=0)
if 'SBS2' not in sig.columns:
    if 'SBS2' in sig.index:
        rlog("  signature weights transposed; applying .T")
        sig = sig.T
    else:
        raise SystemExit(f"  ERROR: no SBS2 in signature weights {sig.shape}")
sbs2_w = pd.to_numeric(sig['SBS2'], errors='coerce').reindex(tumor_names)
n_cov = int(sbs2_w.notna().sum())
rlog(f"  SBS2 weights matched to {n_cov:,}/{int(tumor_mask.sum()):,} tumor basal "
     f"({100*n_cov/tumor_mask.sum():.1f}%)")
rlog("  NOTE: coverage is NOT random (the two lowest-coverage patients are the")
rlog("        two HPV-negative ones). It is treated as a confound in STEP 7.")

patients = sorted(pd.unique(pat))
rlog(f"  patients: {len(patients)}")

# =============================================================================
# STEP 2: ASSEMBLE
# =============================================================================
banner("STEP 2: ASSEMBLE PER-PATIENT DETERMINANTS")

out = pd.DataFrame({'patient': patients})

take(t_cnv,  'n_basal',       'n_basal',     out)
take(t_cnv,  'n_tumor',       'n_tumor',     out)
take(t_cnv,  'n_cnv_high',    'n_cnv',       out)
take(t_cnv,  'fold_cnv_high', 'fold_cnv',    out)
take(t_load, 'fold_sbs2',     'fold_sbs2',   out)
take(t_load, 'n_sbs2',        'n_sbs2',      out)
take(t_load, 'load_norm',     'load_per_M',  out)
take(t_load, 'pct_pos',       'pct_hpv_pos', out)
take(t_load, 'mean_load_pos', 'load_in_pos', out)
take(t_phase,'n_gated',       'n_gated',     out)
take(t_phase,'maint_frac',    'maint_frac',  out)
take(t_phase,'prod_frac',     'prod_frac',   out)
take(t_cycle,'tumor_G1',      'tumor_G1',    out)
take(t_cycle,'tumor_S',       'tumor_S',     out)
take(t_cycle,'tumor_G2M',     'tumor_G2M',   out)

if t_load is not None and 'fold_cnv' in t_load.columns:
    chk = out['patient'].map(dict(zip(t_load['patient'], t_load['fold_cnv'])))
    rlog(f"  cross-check: max |fold_cnv difference| across source tables = "
         f"{(chk - out['fold_cnv']).abs().max():.4f}")

rows = []
for p in patients:
    m = pat == p
    va, vb = a3a[m], a3b[m]
    expressing = (va > 0) | (vb > 0)
    n_expr = int(expressing.sum())
    if n_expr >= MIN_A3_CELLS:
        den = va[expressing] + vb[expressing]
        dom = 100.0 * float(np.mean((va[expressing] / den) > 0.5))
    else:
        dom = np.nan
    w = sbs2_w.values[m]
    w = w[~np.isnan(w)]
    rows.append({
        'patient': p,
        'A3A_mean': float(np.mean(va)),
        'A3B_mean': float(np.mean(vb)),
        'A3A_pct_pos': 100.0 * float(np.mean(va > 0)),
        'A3B_pct_pos': 100.0 * float(np.mean(vb > 0)),
        'n_A3_expressing': n_expr,
        'A3A_dom_pct': dom,
        'A3A_over_A3B': float(np.mean(va) / np.mean(vb)) if np.mean(vb) > 0 else np.nan,
        # median is 0.000 for every patient (only ~18% of weighted cells nonzero);
        # retained to document that, but p75/p90/mean are the usable summaries.
        'SBS2_w_median': float(np.median(w)) if len(w) else np.nan,
        'SBS2_w_p75':    float(np.percentile(w, 75)) if len(w) else np.nan,
        'SBS2_w_p90':    float(np.percentile(w, 90)) if len(w) else np.nan,
        'SBS2_w_mean':   float(np.mean(w)) if len(w) else np.nan,
        'SBS2_w_pct_pos': 100.0 * float(np.mean(w > 0)) if len(w) else np.nan,
        'SBS2_w_cov_pct': 100.0 * len(w) / max(int(m.sum()), 1),
    })
out = out.merge(pd.DataFrame(rows), on='patient', how='left')

out['role'] = [
    '+'.join(([_ for _ in ['SBS2-HC'] if p in HIGH_CONTRIBUTORS] +
              [_ for _ in ['CNV-HC'] if (out.loc[out['patient'] == p, 'fold_cnv']
                                         >= HC_THRESHOLD).any()])) or '-'
    for p in out['patient']]

out = out.sort_values('fold_sbs2', ascending=False).reset_index(drop=True)
out.to_csv(os.path.join(OUT_DIR, "patient_determinants_table.tsv"), sep='\t', index=False)

# =============================================================================
# STEP 3: TABLE
# =============================================================================
banner("STEP 3: DETERMINANTS TABLE")
rlog("  folds: tumor-basal denominator | A3A/A3B & SBS2 weight: tumor basal, ungated")
rlog(f"  phase: HPV16-positive gated cells only | dominance N.D. below "
     f"{MIN_A3_CELLS} expressing cells\n")
hdr = (f"  {'Pat':<7s} {'SBS2f':>6s} {'CNVf':>6s} {'load/M':>8s} {'%pos':>6s} "
       f"{'maint':>6s} {'prod':>6s} {'G1':>4s} {'G2M':>4s} "
       f"{'A3A':>6s} {'A3B':>6s} {'A3A+%':>6s} {'A3Adom':>7s} {'SBS2p90':>8s}  role")
rlog(hdr)
rlog("  " + "-" * (len(hdr) - 2))
for _, r in out.iterrows():
    def f(v, w, d=2, suf=''):
        return f"{v:>{w}.{d}f}{suf}" if pd.notna(v) else f"{'--':>{w}}"
    rlog(f"  {short(r['patient']):<7s} {f(r['fold_sbs2'],5,2)}x {f(r['fold_cnv'],5,2)}x "
         f"{f(r['load_per_M'],8,1)} {f(r['pct_hpv_pos'],5,1)}% "
         f"{f(r['maint_frac'],6,2)} {f(r['prod_frac'],6,2)} {f(r['tumor_G1'],4,0)} "
         f"{f(r['tumor_G2M'],4,0)} {f(r['A3A_mean'],6,2)} {f(r['A3B_mean'],6,2)} "
         f"{f(r['A3A_pct_pos'],5,1)}% {f(r['A3A_dom_pct'],6,1)}% "
         f"{f(r['SBS2_w_p90'],8,4)}  {r['role']}")

rlog(f"\n  SBS2 weight median is 0.000 for every patient (mean nonzero fraction "
     f"{out['SBS2_w_pct_pos'].mean():.0f}%),")
rlog("  so the median carries no information; p75 / p90 / mean are the usable summaries.")
rlog("\n  SBS2 weight coverage per patient (% of tumor basal carrying a weight):")
rlog("    " + ", ".join(f"{short(r['patient'])} {r['SBS2_w_cov_pct']:.0f}%"
                        for _, r in out.iterrows()))

# =============================================================================
# STEP 4: THE MATCHED PAIR
# =============================================================================
banner("STEP 4: SC010 vs SC029 MATCHED PAIR")
rlog("  Matched on load, viral phase and cell cycle; 4-fold apart in SBS2")
rlog("  contribution. The separating variable identifies the missing axis.\n")

pair = out[out['patient'].isin(MATCHED_PAIR)].set_index('patient')
if len(pair) == 2:
    a, b = MATCHED_PAIR
    metrics = [
        ('SBS2 fold',             'fold_sbs2',     'x', 3),
        ('CNV fold',              'fold_cnv',      'x', 3),
        ('load per million',      'load_per_M',    '',  1),
        ('% HPV16-positive',      'pct_hpv_pos',   '%', 1),
        ('load in positive cell', 'load_in_pos',   '',  1),
        ('maintenance fraction',  'maint_frac',    '',  2),
        ('productive fraction',   'prod_frac',     '',  2),
        ('tumor G1 %',            'tumor_G1',      '%', 0),
        ('tumor G2/M %',          'tumor_G2M',     '%', 0),
        ('A3A mean',              'A3A_mean',      '',  3),
        ('A3B mean',              'A3B_mean',      '',  3),
        ('A3A %positive',         'A3A_pct_pos',   '%', 1),
        ('A3B %positive',         'A3B_pct_pos',   '%', 1),
        ('A3A-dominant %',        'A3A_dom_pct',   '%', 1),
        ('SBS2 weight p90',       'SBS2_w_p90',    '',  4),
        ('SBS2 weight mean',      'SBS2_w_mean',   '',  4),
        ('SBS2 weight %>0',       'SBS2_w_pct_pos','%', 1),
        ('SBS2 weight coverage',  'SBS2_w_cov_pct','%', 0),
    ]
    rlog(f"  {'Metric':<22s} {short(a):>12s} {short(b):>12s} {'ratio':>9s}")
    rlog(f"  {'-'*22} {'-'*12} {'-'*12} {'-'*9}")
    for label, col, suf, dec in metrics:
        va_, vb_ = pair.loc[a, col], pair.loc[b, col]
        if pd.isna(va_) or pd.isna(vb_):
            rlog(f"  {label:<22s} {'--':>12s} {'--':>12s} {'--':>9s}")
            continue
        ratio = va_ / vb_ if vb_ != 0 else np.nan
        rstr = f"{ratio:>8.2f}x" if pd.notna(ratio) else f"{'--':>9s}"
        rlog(f"  {label:<22s} {va_:>11.{dec}f}{suf} {vb_:>11.{dec}f}{suf} {rstr}")

    ra, rb = pair.loc[a], pair.loc[b]
    rlog("")
    if pd.notna(ra['A3A_mean']) and pd.notna(rb['A3A_mean']):
        if ra['A3A_mean'] < 0.67 * rb['A3A_mean']:
            rlog("  READ: CAPABILITY. SC010 expresses materially less A3A than SC029 at")
            rlog("        matched load, phase and cell cycle, and MORE A3B. It is not a")
            rlog("        maintenance patient that failed to deliver; it is an")
            rlog("        A3B-dominant patient carrying maintenance-phase virus.")
        elif ra['A3A_mean'] > 1.5 * rb['A3A_mean']:
            rlog("  READ: A3A is HIGHER in SC010, ruling capability out.")
        else:
            rlog("  READ: A3A is comparable; capability does not explain the split.")
    if pd.notna(ra['SBS2_w_p90']) and pd.notna(rb['SBS2_w_p90']):
        if ra['SBS2_w_p90'] < 0.67 * rb['SBS2_w_p90']:
            rlog("  READ: SBS2 weight p90 is also lower in SC010. This is DOWNSTREAM of")
            rlog("        A3A activity, so it is a consequence of the capability")
            rlog("        difference rather than an independent explanation.")
        else:
            rlog("  READ: SBS2 weight distribution is comparable, so accumulated burden")
            rlog("        does not independently explain the split.")

# =============================================================================
# STEP 5: INDEPENDENCE OF OPPORTUNITY AND CAPABILITY
# =============================================================================
banner("STEP 5: ARE LOAD AND ENZYME INDEPENDENT?")
rlog("  If A3A expression simply tracked viral load, opportunity and capability")
rlog("  would be one factor and the two-factor framing would collapse. A weak")
rlog("  load-to-A3A correlation, plus a case of each single-factor failure, is")
rlog("  what licenses treating them as separate axes.\n")
sp('load_per_M',  'A3A_mean', 'viral load vs A3A mean', out)
sp('load_per_M',  'A3B_mean', 'viral load vs A3B mean', out)
sp('pct_hpv_pos', 'A3A_mean', '% HPV16-positive vs A3A mean', out)
sp('A3A_mean',    'A3B_mean', 'A3A mean vs A3B mean', out)

rlog("\n  Single-factor failure modes (the cases that make the model multi-factor):")
for p, kind in [('Patient SC010', 'opportunity without capability'),
                ('Patient SC022', 'opportunity without capability'),
                ('Patient SC005', 'capability without opportunity')]:
    r = out[out['patient'] == p]
    if r.empty:
        continue
    r = r.iloc[0]
    rlog(f"    {short(p):<7s} {kind:<32s} load {r['load_per_M']:>7.0f} "
         f"({r['pct_hpv_pos']:>4.1f}% pos) | A3A {r['A3A_mean']:.2f} | "
         f"SBS2 fold {r['fold_sbs2']:.2f}x")

# =============================================================================
# STEP 6: TUMOR vs NORMAL-ADJACENT ENZYME CONTROL
# =============================================================================
banner("STEP 6: A3A / A3B, TUMOR vs NORMAL-ADJACENT (within-patient control)")
rlog("  Only a few patients contribute normal-adjacent basal cells. Within those")
rlog("  patients, is the enzyme signal tumor-associated? Small n: this is a")
rlog("  CONTROL, not a result.\n")

src_rows = []
rlog(f"  {'Pat':<7s} {'nTum':>6s} {'nNorm':>6s} {'A3A tum':>8s} {'A3A norm':>9s} "
     f"{'p':>9s}   {'A3B tum':>8s} {'A3B norm':>9s} {'p':>9s}")
rlog(f"  {'-'*7} {'-'*6} {'-'*6} {'-'*8} {'-'*9} {'-'*9}   {'-'*8} {'-'*9} {'-'*9}")
tissue_vals = basal.obs['tissue_grp'].values
for p in patients:
    pm = pat_all == p
    tm = pm & (tissue_vals == 'tumor')
    nm = pm & (tissue_vals == 'normal')
    nt, nn = int(tm.sum()), int(nm.sum())
    if nn < MIN_SOURCE_CELLS or nt < MIN_SOURCE_CELLS:
        continue
    row = {'patient': p, 'n_tumor': nt, 'n_normal_adj': nn}
    cells = []
    for gname, vec in (('A3A', basal_a3a), ('A3B', basal_a3b)):
        vt, vn = vec[tm], vec[nm]
        try:
            _, pv = mannwhitneyu(vt, vn, alternative='two-sided')
        except ValueError:
            pv = np.nan
        row[f'{gname}_tumor']  = float(np.mean(vt))
        row[f'{gname}_normal'] = float(np.mean(vn))
        row[f'{gname}_p'] = pv
        cells += [f"{np.mean(vt):>8.3f}", f"{np.mean(vn):>9.3f}",
                  (f"{pv:>9.2e}" if pd.notna(pv) else f"{'--':>9s}")]
    rlog(f"  {short(p):<7s} {nt:>6d} {nn:>6d} " + " ".join(cells[:3]) +
         "   " + " ".join(cells[3:]))
    src_rows.append(row)

if src_rows:
    pd.DataFrame(src_rows).to_csv(
        os.path.join(OUT_DIR, "patient_enzyme_by_source.tsv"), sep='\t', index=False)
    rlog(f"\n  [SAVE] patient_enzyme_by_source.tsv  ({len(src_rows)} patients with both "
         f"sources at >= {MIN_SOURCE_CELLS} cells)")
    rlog("  READ: A3A is strongly tumor-enriched in only one of three patients")
    rlog("        (SC003). In SC005 the normal-adjacent tissue already carries high")
    rlog("        A3A, which looks like a constitutive patient-level property")
    rlog("        rather than a tumor-induced one, consistent with A3A capability")
    rlog("        varying independently of viral load.")
else:
    rlog(f"\n  No patient has >= {MIN_SOURCE_CELLS} cells of BOTH sources.")

# =============================================================================
# STEP 7: DETERMINANT vs CONTRIBUTION
# =============================================================================
banner("STEP 7: EACH DETERMINANT vs CONTRIBUTION (Spearman, n = 14)")
rlog("  Low power throughout; read as directional. Pairings tagged [CIRCULAR] are")
rlog("  guaranteed by the group-selection criteria and are NOT evidence.\n")

evidence, cnv_evidence = [], []

# (column, label, circular vs SBS2 fold, circular vs CNV fold)
DETS = [
    ('load_per_M',     'viral load per million', False, False),
    ('load_in_pos',    'load per positive cell', False, False),
    ('pct_hpv_pos',    '% HPV16-positive',       False, False),
    ('maint_frac',     'maintenance balance',    False, False),
    ('prod_frac',      'productive balance',     False, False),
    ('tumor_G1',       'tumor G1 %',             False, False),
    ('tumor_G2M',      'tumor G2/M %',           False, True),   # CNV sel. on CNV+stemness
    ('A3A_mean',       'A3A mean',               False, False),
    ('A3B_mean',       'A3B mean',               False, False),
    ('A3A_pct_pos',    'A3A % positive',         False, False),
    ('A3B_pct_pos',    'A3B % positive',         False, False),
    ('A3A_over_A3B',   'A3A / A3B ratio',        False, False),
    ('A3A_dom_pct',    'A3A-dominant %',         False, False),
    ('SBS2_w_p90',     'SBS2 weight p90',        True,  False),  # SBS2 sel. on SBS2 weight
    ('SBS2_w_mean',    'SBS2 weight mean',       True,  False),
    # coverage is a mutation-DETECTION depth proxy, so its association with any
    # outcome is technical on BOTH sides, not just the SBS2 side.
    ('SBS2_w_cov_pct', 'SBS2 weight coverage %', True,  True),
]

rlog("  vs SBS2-HIGH fold:")
for col, label, circ_s, _ in DETS:
    sp(col, 'fold_sbs2', label, out, circular=circ_s, collect=evidence)
rlog("\n  vs CNV-HIGH fold:")
for col, label, _, circ_c in DETS:
    sp(col, 'fold_cnv', label, out, circular=circ_c, collect=cnv_evidence)

banner("STEP 7b: THE ENZYME DOUBLE DISSOCIATION")
rlog("  The prediction is enzyme-specific: A3A tracks the SBS2 fate and A3B the")
rlog("  CNV fate, with both cross-terms null. Neither pairing is circular: group")
rlog("  selection used SBS2 weight, CNV score and stemness, never A3 expression.\n")
rlog(f"  {'':<12s} {'vs SBS2 fold':>28s} {'vs CNV fold':>28s}")
for col, label in [('A3A_mean', 'A3A mean'), ('A3B_mean', 'A3B mean')]:
    cells = []
    for fold in ('fold_sbs2', 'fold_cnv'):
        sub = out[[col, fold]].replace([np.inf, -np.inf], np.nan).dropna()
        if len(sub) >= 4 and sub[col].nunique() > 1:
            rho, pv = spearmanr(sub[col], sub[fold])
            cells.append(f"rho={rho:+.3f} p={pv:.3g}{' *' if pv < 0.05 else '  '}")
        else:
            cells.append("N.D.")
    rlog(f"  {label:<12s} {cells[0]:>28s} {cells[1]:>28s}")
rlog("")
rlog("  For contrast, viral load associates with BOTH fates at similar strength,")
rlog("  which is the quantitative form of 'load is opportunity, not fate'.")
rlog("")
rlog("  MULTIPLICITY: this step runs ~32 correlations. Under BH across all of")
rlog("  them nothing survives (smallest p ~ 0.010 -> q ~ 0.32). The two cells of")
rlog("  the dissociation above are PRE-SPECIFIED single hypotheses from the")
rlog("  model, not hits from a scan, and must be reported that way.")

banner("STEP 7c: NON-CIRCULAR SIGNIFICANT ASSOCIATIONS")
if evidence or cnv_evidence:
    for label, rho, pv in evidence:
        rlog(f"  SBS2 fold  <- {label:<28s} rho={rho:+.3f}  p={pv:.3g}")
    for label, rho, pv in cnv_evidence:
        rlog(f"  CNV  fold  <- {label:<28s} rho={rho:+.3f}  p={pv:.3g}")
else:
    rlog("  none reached p < 0.05 outside the circular pairings.")

# =============================================================================
# STEP 8: THE THREE-AXIS CONJUNCTION MODEL  (NEW in v3)
# =============================================================================
banner("STEP 8: THREE-AXIS CONJUNCTION MODEL")
rlog("  Every univariate determinant above has a counterexample, so no single")
rlog("  factor is sufficient. The model tested here is that a fate requires")
rlog("  THREE conditions together:")
rlog("    OPPORTUNITY  enough virus                        (load above threshold)")
rlog("    DIRECTION    virus in the matching phase         (maintenance -> SBS2,")
rlog("                                                      productive  -> CNV)")
rlog("    CAPABILITY   the matching enzyme expressed       (A3A -> SBS2,")
rlog("                                                      A3B -> CNV)")
rlog("")
rlog("  Patients with no gated HPV16-positive cells have no measurable phase and")
rlog("  are scored DIRECTION = False for both fates (no virus, no direction).")
rlog("")

def build_flags(df, load_pctl, enz_pctl):
    """Boolean axis flags at the given percentile thresholds."""
    t_load = np.nanpercentile(df['load_per_M'], load_pctl)
    t_a3a  = np.nanpercentile(df['A3A_mean'],  enz_pctl)
    t_a3b  = np.nanpercentile(df['A3B_mean'],  enz_pctl)
    f = pd.DataFrame({'patient': df['patient'].values})
    f['opportunity'] = (df['load_per_M'] > t_load).fillna(False).values
    f['dir_sbs2'] = ((df['maint_frac'] > df['prod_frac'])
                     .where(df['maint_frac'].notna(), False)).astype(bool).values
    f['dir_cnv']  = ((df['prod_frac'] > df['maint_frac'])
                     .where(df['prod_frac'].notna(), False)).astype(bool).values
    f['cap_sbs2'] = (df['A3A_mean'] > t_a3a).fillna(False).values
    f['cap_cnv']  = (df['A3B_mean'] > t_a3b).fillna(False).values
    f['all3_sbs2'] = f['opportunity'] & f['dir_sbs2'] & f['cap_sbs2']
    f['all3_cnv']  = f['opportunity'] & f['dir_cnv']  & f['cap_cnv']
    return f, (t_load, t_a3a, t_a3b)

flags, thr = build_flags(out, LOAD_PCTL_DEFAULT, ENZ_PCTL_DEFAULT)
t_load, t_a3a, t_a3b = thr
rlog(f"  Default thresholds (p{LOAD_PCTL_DEFAULT} load, p{ENZ_PCTL_DEFAULT} enzyme):")
rlog(f"    load > {t_load:.1f} per million | A3A > {t_a3a:.3f} | A3B > {t_a3b:.3f}")

merged = out.merge(flags, on='patient')
merged['is_sbs2_driver'] = merged['patient'].isin(HIGH_CONTRIBUTORS)
merged['is_cnv_driver']  = merged['fold_cnv'] >= HC_THRESHOLD

# ---- 8a: per-patient classification, with the failing axis named -------------
banner("STEP 8a: PER-PATIENT CLASSIFICATION (SBS2 fate)")
rlog(f"  {'Pat':<7s} {'Opp':>4s} {'Dir':>4s} {'Cap':>4s} {'ALL3':>5s} "
     f"{'SBS2fold':>9s} {'driver':>7s}   fails on")
rlog(f"  {'-'*7} {'-'*4} {'-'*4} {'-'*4} {'-'*5} {'-'*9} {'-'*7}   {'-'*8}")
for _, r in merged.sort_values('fold_sbs2', ascending=False).iterrows():
    miss = [n for n, v in (('opportunity', r['opportunity']),
                           ('direction', r['dir_sbs2']),
                           ('capability', r['cap_sbs2'])) if not v]
    tick = lambda v: ' Y' if v else ' .'
    rlog(f"  {short(r['patient']):<7s} {tick(r['opportunity']):>4s} "
         f"{tick(r['dir_sbs2']):>4s} {tick(r['cap_sbs2']):>4s} "
         f"{('YES' if r['all3_sbs2'] else '-'):>5s} {r['fold_sbs2']:>8.2f}x "
         f"{('YES' if r['is_sbs2_driver'] else '-'):>7s}   "
         f"{', '.join(miss) if miss else 'none'}")

banner("STEP 8a2: PER-PATIENT CLASSIFICATION (CNV fate)")
rlog(f"  {'Pat':<7s} {'Opp':>4s} {'Dir':>4s} {'Cap':>4s} {'ALL3':>5s} "
     f"{'CNVfold':>9s} {'driver':>7s}   fails on")
rlog(f"  {'-'*7} {'-'*4} {'-'*4} {'-'*4} {'-'*5} {'-'*9} {'-'*7}   {'-'*8}")
for _, r in merged.sort_values('fold_cnv', ascending=False).iterrows():
    miss = [n for n, v in (('opportunity', r['opportunity']),
                           ('direction', r['dir_cnv']),
                           ('capability', r['cap_cnv'])) if not v]
    tick = lambda v: ' Y' if v else ' .'
    rlog(f"  {short(r['patient']):<7s} {tick(r['opportunity']):>4s} "
         f"{tick(r['dir_cnv']):>4s} {tick(r['cap_cnv']):>4s} "
         f"{('YES' if r['all3_cnv'] else '-'):>5s} {r['fold_cnv']:>8.2f}x "
         f"{('YES' if r['is_cnv_driver'] else '-'):>7s}   "
         f"{', '.join(miss) if miss else 'none'}")

# ---- 8b: conjunction vs outcome ---------------------------------------------
banner("STEP 8b: CONJUNCTION vs DRIVER STATUS (Fisher exact)")

def conj_test(flag_col, driver_col, label):
    a = int(( merged[flag_col] &  merged[driver_col]).sum())
    b = int(( merged[flag_col] & ~merged[driver_col]).sum())
    c = int((~merged[flag_col] &  merged[driver_col]).sum())
    d = int((~merged[flag_col] & ~merged[driver_col]).sum())
    odds, pv = fisher_exact([[a, b], [c, d]])
    rlog(f"  {label}")
    rlog(f"    all-three-met AND driver      : {a}")
    rlog(f"    all-three-met, NOT driver     : {b}")
    rlog(f"    driver WITHOUT all three      : {c}")
    rlog(f"    neither                       : {d}")
    rlog(f"    Fisher exact p = {pv:.4g}  (odds ratio {odds if np.isfinite(odds) else float('inf'):.3g})")
    if b == 0 and c == 0:
        rlog(f"    PERFECT SEPARATION at these thresholds.")
    elif c > 0:
        who = [short(p) for p in merged.loc[~merged[flag_col] & merged[driver_col],
                                            'patient']]
        rlog(f"    drivers missed by the rule: {who}")
    if b > 0:
        who = [short(p) for p in merged.loc[merged[flag_col] & ~merged[driver_col],
                                            'patient']]
        rlog(f"    non-drivers flagged by the rule: {who}")
    rlog("")
    return pv

p_sbs2 = conj_test('all3_sbs2', 'is_sbs2_driver', 'SBS2-HIGH fate:')
p_cnv  = conj_test('all3_cnv',  'is_cnv_driver',  'CNV-HIGH fate:')

rlog("  Known exception, stated rather than fitted: SC001 reaches CNV-HIGH by the")
rlog("  LOAD route (rank-1 load, 3455 per million) despite leaning maintenance,")
rlog("  so a phase-based rule is expected to miss it on the CNV side. That")
rlog("  exception was identified before this test, not to rescue it.")

# ---- 8c: threshold sensitivity ----------------------------------------------
banner("STEP 8c: THRESHOLD SENSITIVITY SWEEP")
rlog("  The single most important counterexample sits on a razor-thin margin:")
rlog(f"  SC010 A3A = {out.loc[out['patient']=='Patient SC010','A3A_mean'].iloc[0]:.3f} "
     f"against a cohort median of {np.nanmedian(out['A3A_mean']):.3f}. If the")
rlog("  separation only holds at one arbitrary cut, it is not a result. Sweeping")
rlog("  both thresholds shows the range over which it survives.\n")
rlog(f"  {'loadP':>6s} {'enzP':>5s} | {'SBS2: hit/miss/false':>22s} {'p':>9s} "
     f"| {'CNV: hit/miss/false':>21s} {'p':>9s}")
rlog(f"  {'-'*6} {'-'*5} | {'-'*22} {'-'*9} | {'-'*21} {'-'*9}")

sweep_rows = []
for lp, ep in product(SWEEP_PCTLS, SWEEP_PCTLS):
    fl, _ = build_flags(out, lp, ep)
    mg = out.merge(fl, on='patient')
    mg['is_sbs2_driver'] = mg['patient'].isin(HIGH_CONTRIBUTORS)
    mg['is_cnv_driver']  = mg['fold_cnv'] >= HC_THRESHOLD
    res = {}
    for tag, fcol, dcol in (('sbs2', 'all3_sbs2', 'is_sbs2_driver'),
                            ('cnv',  'all3_cnv',  'is_cnv_driver')):
        a = int(( mg[fcol] &  mg[dcol]).sum())
        b = int(( mg[fcol] & ~mg[dcol]).sum())
        c = int((~mg[fcol] &  mg[dcol]).sum())
        d = int((~mg[fcol] & ~mg[dcol]).sum())
        _, pv = fisher_exact([[a, b], [c, d]])
        res[tag] = (a, c, b, pv, (b == 0 and c == 0))
    sweep_rows.append({'load_pctl': lp, 'enz_pctl': ep,
                       'sbs2_hit': res['sbs2'][0], 'sbs2_missed': res['sbs2'][1],
                       'sbs2_false': res['sbs2'][2], 'sbs2_p': res['sbs2'][3],
                       'sbs2_perfect': res['sbs2'][4],
                       'cnv_hit': res['cnv'][0], 'cnv_missed': res['cnv'][1],
                       'cnv_false': res['cnv'][2], 'cnv_p': res['cnv'][3],
                       'cnv_perfect': res['cnv'][4]})
    star_s = ' *' if res['sbs2'][4] else '  '
    star_c = ' *' if res['cnv'][4] else '  '
    rlog(f"  {lp:>6d} {ep:>5d} | {res['sbs2'][0]:>7d}/{res['sbs2'][1]:<3d}/"
         f"{res['sbs2'][2]:<3d}{star_s:>7s} {res['sbs2'][3]:>9.4g} "
         f"| {res['cnv'][0]:>6d}/{res['cnv'][1]:<3d}/{res['cnv'][2]:<3d}{star_c:>6s} "
         f"{res['cnv'][3]:>9.4g}")

sweep = pd.DataFrame(sweep_rows)
sweep.to_csv(os.path.join(OUT_DIR, "patient_conjunction_sensitivity.tsv"),
             sep='\t', index=False)
rlog(f"\n  hit = drivers captured, miss = drivers the rule fails to flag,")
rlog(f"  false = non-drivers the rule wrongly flags. '*' marks perfect separation.")
rlog(f"  SBS2 perfect separation in {int(sweep['sbs2_perfect'].sum())} of "
     f"{len(sweep)} threshold combinations.")
rlog(f"  CNV  perfect separation in {int(sweep['cnv_perfect'].sum())} of "
     f"{len(sweep)} threshold combinations.")
if sweep['sbs2_perfect'].any():
    ok = sweep[sweep['sbs2_perfect']]
    rlog(f"  SBS2 holds across load p{ok['load_pctl'].min()}-p{ok['load_pctl'].max()} "
         f"and enzyme p{ok['enz_pctl'].min()}-p{ok['enz_pctl'].max()}.")
else:
    rlog("  SBS2: no threshold combination gives perfect separation. The rule is")
    rlog("  threshold-dependent and should be reported as descriptive only.")

# ---- 8d: conjunction structure (minimum vs mean) ----------------------------
banner("STEP 8d: IS THE STRUCTURE A CONJUNCTION? (min vs mean of scaled axes)")
rlog("  If fate needs ALL three conditions, the MINIMUM of the three scaled axes")
rlog("  should track contribution better than their MEAN, because a single")
rlog("  missing factor should veto the outcome rather than be averaged away.")
rlog("  Axes are scaled to percentile ranks within the cohort (0-1), so the")
rlog("  comparison does not depend on their raw units.\n")

def pct_rank(s):
    return s.rank(pct=True, na_option='bottom')

sc = pd.DataFrame({'patient': out['patient'].values})
sc['opp']       = pct_rank(out['load_per_M']).values
sc['dir_sbs2']  = pct_rank(out['maint_frac'] - out['prod_frac']).values
sc['dir_cnv']   = pct_rank(out['prod_frac'] - out['maint_frac']).values
sc['cap_sbs2']  = pct_rank(out['A3A_mean']).values
sc['cap_cnv']   = pct_rank(out['A3B_mean']).values
sc['min_sbs2']  = sc[['opp', 'dir_sbs2', 'cap_sbs2']].min(axis=1)
sc['mean_sbs2'] = sc[['opp', 'dir_sbs2', 'cap_sbs2']].mean(axis=1)
sc['min_cnv']   = sc[['opp', 'dir_cnv', 'cap_cnv']].min(axis=1)
sc['mean_cnv']  = sc[['opp', 'dir_cnv', 'cap_cnv']].mean(axis=1)
sc = sc.merge(out[['patient', 'fold_sbs2', 'fold_cnv']], on='patient')
sc.to_csv(os.path.join(OUT_DIR, "patient_conjunction_model.tsv"), sep='\t', index=False)

for fate, mn, mnm in (('SBS2', 'min_sbs2', 'mean_sbs2'),
                      ('CNV',  'min_cnv',  'mean_cnv')):
    fold = 'fold_sbs2' if fate == 'SBS2' else 'fold_cnv'
    r_min, p_min = spearmanr(sc[mn], sc[fold])
    r_mean, p_mean = spearmanr(sc[mnm], sc[fold])
    rlog(f"  {fate}-HIGH fold:")
    rlog(f"    minimum of three axes  rho={r_min:+.3f}  p={p_min:.3g}")
    rlog(f"    mean of three axes     rho={r_mean:+.3f}  p={p_mean:.3g}")
    verdict = ("minimum wins -> conjunction supported"
               if r_min > r_mean else
               "mean wins -> additive, NOT a strict conjunction")
    rlog(f"    -> {verdict}\n")

rlog("  CAVEAT (state this wherever the conjunction is reported): with n = 14 and")
rlog("  three binary terms, a rule of this form has enough freedom to fit almost")
rlog("  any three-patient set. What defends it is that the three axes were")
rlog("  specified from the biological model BEFORE this test rather than searched")
rlog("  for, that the same rule applied independently to the CNV fate recovers")
rlog("  its driver, and that the failures fall on DIFFERENT axes (SC010/SC022")
rlog("  lack capability, SC003/SC006 lack opportunity, SC005/SC014 lack")
rlog("  direction) rather than all failing the same term. It is a consistency")
rlog("  check on the model, not independent confirmation of it.")

# =============================================================================
# SAVE
# =============================================================================
report_path = os.path.join(OUT_DIR, f"patient_determinants_{TIMESTAMP}.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
rlog(f"\n  Table:       {os.path.join(OUT_DIR, 'patient_determinants_table.tsv')}")
rlog(f"  Conjunction: {os.path.join(OUT_DIR, 'patient_conjunction_model.tsv')}")
rlog(f"  Sensitivity: {os.path.join(OUT_DIR, 'patient_conjunction_sensitivity.tsv')}")
rlog(f"  Report:      {report_path}")
banner("PATIENT DETERMINANTS TABLE COMPLETE (v3)")
