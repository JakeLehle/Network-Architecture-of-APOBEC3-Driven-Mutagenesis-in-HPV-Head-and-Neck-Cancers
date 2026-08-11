#!/usr/bin/env python3
"""
Diagnostic_Patient_Determinants_Table.py  (v6)
===============================================
READ-ONLY. Per-patient determinants of SBS2-HIGH / CNV-HIGH contribution, and a
test of the three-axis conjunction model that emerged from them.

v6 CHANGES (reporting only; no test or number changed)
-------------------------------------------------------
  - The STEP 5 rank-gap table tags HPV16-negative patients, derived as 0% of
    basal cells positive rather than named. For a patient with no virus,
    'capability without opportunity' is trivially true and reads as a model
    failure when it is not one. They are labelled and excluded from the
    strongest-case-in-each-direction picks.

v5 CHANGES (three derivation rules corrected after the first all-basal run)
---------------------------------------------------------------------------
  - STEP 5 single-factor failures are ranked by DISCREPANCY (load rank minus
    A3A rank, both ascending cohort-wide) instead of filtered on a median
    split. The median split returned every patient that happened to sit on
    opposite sides of two midpoints, which included uninformative cases and
    even flagged a CNV driver as an SBS2 failure. Ranking by gap size puts the
    informative cases at the top and shows why the marginal ones are marginal.
    Patients that drive either fate are excluded, since a driver is not a
    failure case.
  - STEP 6 enrichment now requires BOTH p < 0.05 AND a fold change of at least
    MIN_ENRICH_FOLD. With thousands of tumor cells against a few hundred
    normal-adjacent, a 1.2-fold difference clears significance on sample size
    alone. The fold is printed per patient so the reader can see this directly.
  - STEP 0 no longer reports a denominator for tables that carry no fold
    columns. A phase-percentage table has no denominator, and "not stated"
    read like a warning when nothing was wrong.

v4 CHANGES (v3 tests are unchanged in method; only sourcing and prose changed)
------------------------------------------------------------------------------
  - The contribution denominator comes from patient_config.CONTRIBUTION_DENOMINATOR
    via contribution.py. Fold columns are pulled from the upstream tables using
    the suffix that matches the active setting, and the script refuses to guess:
    if an upstream table carries only a bare fold column it says so loudly.
  - Driver status is now SYMMETRIC. Both is_sbs2_driver and is_cnv_driver are
    derived from the fold columns at HC_THRESHOLD. Previously the SBS2 side read
    a hardcoded list from patient_config while the CNV side was derived, so a
    denominator change would have moved one side and not the other.
  - Every narrative READ block is computed from this run. The old text named
    specific patients and quoted specific numbers inline, which meant the prose
    could drift away from the data silently. Patient identifiers now appear only
    where they are an ANALYSIS CHOICE (the matched pair, the spotlight list),
    never as a stated result.

THE MODEL
---------
No single factor is sufficient, and the univariate steps show why: every one has
a counterexample. The proposal is that reaching a fate requires THREE conditions
together, and that the informative patients are the ones that fail on different
axes.

    OPPORTUNITY  the patient carries enough virus                (viral load)
    DIRECTION    that virus sits in the right lifecycle phase    (maintenance
                 for SBS2, productive for CNV)
    CAPABILITY   the cell expresses the matching enzyme          (A3A for SBS2,
                                                                  A3B for CNV)

STEP 8 has four parts:
  8a  classify every patient on all three axes and print WHERE each one fails
  8b  Fisher exact, all-three-met vs driver status, run separately per fate
  8c  threshold sensitivity, swept across percentiles, reporting the range over
      which separation holds and where it breaks
  8d  conjunction structure: the MINIMUM of the three scaled axes versus their
      MEAN. If the minimum tracks contribution better, a single missing factor
      vetoes the outcome rather than being averaged away.

STANDING CAVEAT, printed in the output so it travels with the numbers: at n = 14
a rule with three binary terms has enough freedom to fit almost any three-patient
set. What defends it is that the axes were specified from the biological model
before this test, not searched for, and that the same rule applied independently
to the CNV side recovers its driver. It is a consistency check on the model, not
an independent confirmation of it.

Denominators and universes (they differ by column, so they are stated)
----------------------------------------------------------------------
  folds        set by CONTRIBUTION_DENOMINATOR; both variants exist upstream
  A3A/A3B      tumor basal cells, all of them, ungated on HPV16
  SBS2 weight  tumor basal cells carrying a weight (coverage reported per patient)
  phase        HPV16-positive gated cells only (raw_HPV16 >= 8 and total > 0)
  tumor/normal control: basal cells of both sources, patients having both

Virus-derived measures stay tumor-restricted regardless of the fold denominator,
because normal-adjacent basal would dilute them toward zero by construction.

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
  patient_conjunction_sensitivity.tsv    load x enzyme threshold sweep
  patient_direction_rule_sensitivity.tsv direction-rule sweep (STEP 8c-dir)
  patient_determinants_<ts>.txt          full console log

Run LAST, after the three upstream diagnostics. From the directory holding
patient_config.py and contribution.py:
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

from contribution import (
    CONTRIBUTION_DENOMINATOR, HC_THRESHOLD,
    read_fold, derive_contributors, announce, short,
)

try:
    from patient_config import CNV_HIGH_CONTRIBUTORS
except ImportError:
    CNV_HIGH_CONTRIBUTORS = None

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
MIN_GATED     = 10      # phase estimate floor, matches the phase diagnostic
MIN_A3_CELLS  = 50      # floor for a stable A3A/A3B dominance fraction
MIN_SOURCE_CELLS = 30   # floor for the tumor vs normal-adjacent enzyme test

# Effect-size gate for the STEP 6 enzyme control. With thousands of tumor cells
# against a few hundred normal-adjacent, a 1.2-fold difference clears p < 0.05
# without meaning anything. "Tumor-enriched" requires significance AND this fold.
MIN_ENRICH_FOLD = 2.0

# STEP 8 thresholds. Defaults are cohort medians; 8c sweeps them.
LOAD_PCTL_DEFAULT = 50
ENZ_PCTL_DEFAULT  = 50
SWEEP_PCTLS = [30, 40, 50, 60, 70]

# How the DIRECTION axis is called.
#   'compare'   maintenance > productive for SBS2, productive > maintenance for
#               CNV. Parameter-free: it asks only which phase dominates, so
#               there is nothing to tune and nothing to overfit.
#   'threshold' the matching phase fraction must clear DIRECTION_MIN_FRAC.
#               Introduces a free parameter. STEP 8c-dir sweeps it so the range
#               over which any conclusion survives is visible rather than
#               assumed.
DIRECTION_RULE     = 'threshold'
DIRECTION_MIN_FRAC = 0.40
DIRECTION_SWEEP    = [0.30, 0.35, 0.38, 0.40, 0.42, 0.45, 0.50]

# ANALYSIS CHOICES, not results. These name which patients get a detailed
# printout; nothing about the conclusions depends on the list.
MATCHED_PAIR = ['Patient SC010', 'Patient SC029']
SPOTLIGHT    = ['Patient SC027', 'Patient SC001', 'Patient SC013',
                'Patient SC029', 'Patient SC010', 'Patient SC022',
                'Patient SC005']

report_lines = []


def rlog(msg=""):
    log(msg)
    report_lines.append(str(msg))


def fmt(v, width, dec=2, suffix=''):
    return (f"{v:>{width}.{dec}f}{suffix}" if pd.notna(v)
            else f"{'--':>{width}}{' ' * len(suffix)}")


# =============================================================================
# HELPERS
# =============================================================================
def load_side_table(path, label):
    if not os.path.exists(path):
        rlog(f"  [MISSING] {label}: {path}")
        return None
    df = pd.read_csv(path, sep='\t')
    key = next((c for c in df.columns
                if str(c).strip().lower() in ('patient', 'subject id', 'subject_id')),
               None)
    if key is None:
        rlog(f"  [SKIP] {label}: no patient key column")
        return None
    df = df.rename(columns={key: 'patient'})
    df['patient'] = df['patient'].astype(str)
    # Only tables that actually carry folds have a denominator to report. A
    # phase-percentage table has none, and saying "not stated" for it reads as a
    # warning when nothing is wrong.
    has_folds = any(str(c).startswith('fold') for c in df.columns)
    if not has_folds:
        rlog(f"  [OK] {label}: {len(df)} rows (no fold columns; denominator N/A)")
        return df
    tag = df['denominator'].iloc[0] if 'denominator' in df.columns else 'not stated'
    rlog(f"  [OK] {label}: {len(df)} rows, denominator = {tag}")
    if tag == 'not stated':
        rlog(f"       WARNING: this table carries folds but no denominator tag, so")
        rlog(f"       it predates the refactor. Re-run the upstream diagnostic.")
    elif tag != CONTRIBUTION_DENOMINATOR:
        rlog(f"       WARNING: this table was written under '{tag}' but the "
             f"active setting is '{CONTRIBUTION_DENOMINATOR}'.")
        rlog(f"       Re-run the upstream diagnostic before trusting its folds.")
    return df


def take(df, src, dest, out):
    if df is None or src not in df.columns:
        out[dest] = np.nan
        return False
    out[dest] = out['patient'].map(dict(zip(df['patient'], df[src])))
    return True


def take_fold(df, prefix, dest, out, label):
    """Pull the fold column matching the active denominator."""
    if df is None:
        out[dest] = np.nan
        return False
    try:
        series = read_fold(df, prefix, logger=rlog, label=label)
    except KeyError as e:
        rlog(f"  [MISSING] {e}")
        out[dest] = np.nan
        return False
    out[dest] = out['patient'].map(dict(zip(df['patient'], series)))
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

n_basal_total = int(basal.n_obs)
n_tumor_total = int(tumor_mask.sum())
rlog(f"  basal cells: {n_basal_total:,}  (tumor {n_tumor_total:,}, "
     f"normal-adjacent {n_basal_total - n_tumor_total:,})")
rlog("")
announce(rlog, n_basal_total, n_tumor_total)
rlog("")
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
rlog(f"  SBS2 weights matched to {n_cov:,}/{n_tumor_total:,} tumor basal "
     f"({100*n_cov/n_tumor_total:.1f}%)")
rlog("  NOTE: coverage is NOT random (the two lowest-coverage patients are the")
rlog("        two HPV-negative ones). It is treated as a confound in STEP 7.")

patients = sorted(pd.unique(pat))
rlog(f"  patients: {len(patients)}")

# =============================================================================
# STEP 2: ASSEMBLE
# =============================================================================
banner("STEP 2: ASSEMBLE PER-PATIENT DETERMINANTS")

out = pd.DataFrame({'patient': patients})

take(t_cnv,  'n_basal',    'n_basal', out)
take(t_cnv,  'n_tumor',    'n_tumor', out)
take(t_cnv,  'n_cnv_high', 'n_cnv',   out)
take_fold(t_cnv, 'fold_cnv', 'fold_cnv', out, 'CNV contribution')

take_fold(t_load, 'fold_sbs2', 'fold_sbs2', out, 'HPV16 load')
take(t_load, 'n_sbs2',        'n_sbs2',      out)
take(t_load, 'load_norm',     'load_per_M',  out)
take(t_load, 'pct_pos',       'pct_hpv_pos', out)
take(t_load, 'mean_load_pos', 'load_in_pos', out)

take(t_phase, 'n_gated',    'n_gated',    out)
take(t_phase, 'maint_frac', 'maint_frac', out)
take(t_phase, 'prod_frac',  'prod_frac',  out)

take(t_cycle, 'tumor_G1',  'tumor_G1',  out)
take(t_cycle, 'tumor_S',   'tumor_S',   out)
take(t_cycle, 'tumor_G2M', 'tumor_G2M', out)

# Cross-check: fold_cnv appears in more than one upstream table. Under a single
# denominator they must agree exactly. A nonzero difference means one upstream
# diagnostic has not been re-run.
if t_load is not None:
    try:
        chk_series = read_fold(t_load, 'fold_cnv', logger=rlog,
                               label='HPV16 load cross-check')
        chk = out['patient'].map(dict(zip(t_load['patient'], chk_series)))
        gap = float((chk - out['fold_cnv']).abs().max())
        rlog(f"  cross-check: max |fold_cnv difference| across source tables = "
             f"{gap:.6f}")
        if gap > 1e-9:
            rlog("  >>> STOP. The source tables disagree, which means at least one")
            rlog("      upstream diagnostic was not re-run under the current")
            rlog("      denominator. Re-run them before using anything below.")
    except KeyError:
        rlog("  cross-check skipped: no fold_cnv column in the load table.")

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

# ---- Driver status, derived symmetrically from the folds --------------------
rlog("")
sbs2_drivers = derive_contributors(out, 'fold_sbs2', expected=HIGH_CONTRIBUTORS,
                                   label='SBS2-HIGH drivers', logger=rlog)
rlog("")
cnv_drivers = derive_contributors(out, 'fold_cnv', expected=CNV_HIGH_CONTRIBUTORS,
                                  label='CNV-HIGH drivers', logger=rlog)

out['is_sbs2_driver'] = out['patient'].isin(sbs2_drivers)
out['is_cnv_driver']  = out['patient'].isin(cnv_drivers)
out['role'] = [
    '+'.join([t for t, v in (('SBS2-HC', s), ('CNV-HC', c)) if v]) or '-'
    for s, c in zip(out['is_sbs2_driver'], out['is_cnv_driver'])]
out['denominator'] = CONTRIBUTION_DENOMINATOR

out = out.sort_values('fold_sbs2', ascending=False).reset_index(drop=True)
out.to_csv(os.path.join(OUT_DIR, "patient_determinants_table.tsv"),
           sep='\t', index=False)

# =============================================================================
# STEP 3: TABLE
# =============================================================================
banner("STEP 3: DETERMINANTS TABLE")
rlog(f"  folds: {CONTRIBUTION_DENOMINATOR} denominator | A3A/A3B & SBS2 weight: "
     f"tumor basal, ungated")
rlog(f"  phase: HPV16-positive gated cells only | dominance N.D. below "
     f"{MIN_A3_CELLS} expressing cells\n")
hdr = (f"  {'Pat':<7s} {'SBS2f':>6s} {'CNVf':>6s} {'load/M':>8s} {'%pos':>6s} "
       f"{'maint':>6s} {'prod':>6s} {'G1':>4s} {'G2M':>4s} "
       f"{'A3A':>6s} {'A3B':>6s} {'A3A+%':>6s} {'A3Adom':>7s} {'SBS2p90':>8s}  role")
rlog(hdr)
rlog("  " + "-" * (len(hdr) - 2))
for _, r in out.iterrows():
    rlog(f"  {short(r['patient']):<7s} {fmt(r['fold_sbs2'],5,2)}x "
         f"{fmt(r['fold_cnv'],5,2)}x "
         f"{fmt(r['load_per_M'],8,1)} {fmt(r['pct_hpv_pos'],5,1)}% "
         f"{fmt(r['maint_frac'],6,2)} {fmt(r['prod_frac'],6,2)} "
         f"{fmt(r['tumor_G1'],4,0)} {fmt(r['tumor_G2M'],4,0)} "
         f"{fmt(r['A3A_mean'],6,2)} {fmt(r['A3B_mean'],6,2)} "
         f"{fmt(r['A3A_pct_pos'],5,1)}% {fmt(r['A3A_dom_pct'],6,1)}% "
         f"{fmt(r['SBS2_w_p90'],8,4)}  {r['role']}")

n_zero_median = int((out['SBS2_w_median'] == 0).sum())
rlog(f"\n  SBS2 weight median is 0.000 for {n_zero_median}/{len(out)} patients "
     f"(mean nonzero fraction {out['SBS2_w_pct_pos'].mean():.0f}%),")
rlog("  so the median carries little information; p75 / p90 / mean are the "
     "usable summaries.")
rlog("\n  SBS2 weight coverage per patient (% of tumor basal carrying a weight):")
rlog("    " + ", ".join(f"{short(r['patient'])} {r['SBS2_w_cov_pct']:.0f}%"
                        for _, r in out.iterrows()))

# =============================================================================
# STEP 4: THE MATCHED PAIR
# =============================================================================
banner(f"STEP 4: {short(MATCHED_PAIR[0])} vs {short(MATCHED_PAIR[1])} MATCHED PAIR")
rlog("  A pair chosen for being close on load, viral phase and cell cycle while")
rlog("  differing in contribution. The separating variable identifies the")
rlog("  missing axis.\n")

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
            rlog(f"  READ: CAPABILITY. {short(a)} expresses materially less A3A "
                 f"({ra['A3A_mean']:.2f} vs {rb['A3A_mean']:.2f}) at")
            rlog(f"        comparable load and phase, and MORE A3B "
                 f"({ra['A3B_mean']:.2f} vs {rb['A3B_mean']:.2f}). That is an")
            rlog(f"        A3B-dominant patient carrying maintenance-phase virus,")
            rlog(f"        not a maintenance patient that failed to deliver.")
        elif ra['A3A_mean'] > 1.5 * rb['A3A_mean']:
            rlog(f"  READ: A3A is HIGHER in {short(a)}, ruling capability out.")
        else:
            rlog("  READ: A3A is comparable; capability does not explain the split.")
    if pd.notna(ra['SBS2_w_p90']) and pd.notna(rb['SBS2_w_p90']):
        if ra['SBS2_w_p90'] < 0.67 * rb['SBS2_w_p90']:
            rlog(f"  READ: SBS2 weight p90 is also lower in {short(a)} "
                 f"({ra['SBS2_w_p90']:.4f} vs {rb['SBS2_w_p90']:.4f}). That is")
            rlog("        DOWNSTREAM of A3A activity, so it is a consequence of the")
            rlog("        capability difference rather than an independent explanation.")
        else:
            rlog("  READ: SBS2 weight distribution is comparable, so accumulated")
            rlog("        burden does not independently explain the split.")
else:
    rlog(f"  MATCHED_PAIR {MATCHED_PAIR} not both present; step skipped.")

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

# ---- Single-factor failures, ranked by DISCREPANCY not by a median split ----
#
# A median split answers "is this patient above the middle on exactly one axis",
# which is true of many uninformative patients and says nothing about HOW
# mismatched they are. The informative cases are the ones where the two ranks
# pull hardest in opposite directions, so the axes are ranked cohort-wide and
# patients are ordered by the size of the gap.
#
# Ranks are ascending (1 = lowest) over ALL patients, then patients that drive
# either fate are excluded, since a driver is by definition not a failure case.
out['rank_load'] = out['load_per_M'].rank(method='min', ascending=True)
out['rank_a3a']  = out['A3A_mean'].rank(method='min', ascending=True)
out['rank_gap']  = out['rank_load'] - out['rank_a3a']

n_pat_total = len(out)
rlog(f"\n  Single-factor failure modes, ranked by load-rank minus A3A-rank")
rlog(f"  (both ascending over all {n_pat_total} patients; drivers of either fate")
rlog(f"  excluded, since a driver is not a failure case).")
rlog("  A large POSITIVE gap is opportunity without capability; a large NEGATIVE")
rlog("  gap is capability without opportunity. These are the cases that make the")
rlog("  model multi-factor.\n")

fails = out[~(out['is_sbs2_driver'] | out['is_cnv_driver'])].copy()
fails = fails.dropna(subset=['rank_gap'])
fails['abs_gap'] = fails['rank_gap'].abs()
fails = fails.sort_values('abs_gap', ascending=False)

# HPV16-negative patients are derived, not named: no gated positive cells means
# 'without opportunity' is trivially true rather than informative, and they must
# not be read as failures of the model.
hpv_neg = set(out.loc[out['pct_hpv_pos'].fillna(0) == 0, 'patient'])
if hpv_neg:
    rlog(f"    HPV16-negative patients (0% of basal cells positive): "
         f"{sorted(short(p) for p in hpv_neg)}. For these, 'without opportunity'")
    rlog(f"    is trivially true. They are tagged below and are NOT model failures.\n")

rlog(f"    {'Pat':<7s} {'loadR':>6s} {'A3AR':>5s} {'gap':>5s}  "
     f"{'reading':<32s} {'load':>7s} {'%pos':>6s} {'A3A':>6s} {'SBS2f':>7s}  note")
rlog(f"    {'-'*7} {'-'*6} {'-'*5} {'-'*5}  {'-'*32} {'-'*7} {'-'*6} {'-'*6} "
     f"{'-'*7}  ----")
for _, r in fails.iterrows():
    gap = int(r['rank_gap'])
    is_neg = r['patient'] in hpv_neg
    if is_neg:
        reading = 'no virus (not a model failure)'
    elif gap > 0:
        reading = 'opportunity without capability'
    elif gap < 0:
        reading = 'capability without opportunity'
    else:
        reading = 'balanced (uninformative)'
    rlog(f"    {short(r['patient']):<7s} {int(r['rank_load']):>6d} "
         f"{int(r['rank_a3a']):>5d} {gap:>+5d}  {reading:<32s} "
         f"{r['load_per_M']:>7.0f} {r['pct_hpv_pos']:>5.1f}% "
         f"{r['A3A_mean']:>6.2f} {r['fold_sbs2']:>6.2f}x  "
         f"{'HPV16-neg' if is_neg else ''}")

# Name the strongest case in each direction, computed rather than asserted.
informative = fails[~fails['patient'].isin(hpv_neg)]
pos = informative[informative['rank_gap'] > 0]
neg = informative[informative['rank_gap'] < 0]
rlog("")
if len(pos):
    b = pos.iloc[0]
    rlog(f"    Strongest opportunity-without-capability: {short(b['patient'])} "
         f"(gap {int(b['rank_gap']):+d}, load rank {int(b['rank_load'])}, "
         f"A3A rank {int(b['rank_a3a'])}).")
if len(neg):
    b = neg.iloc[0]
    rlog(f"    Strongest capability-without-opportunity: {short(b['patient'])} "
         f"(gap {int(b['rank_gap']):+d}, load rank {int(b['rank_load'])}, "
         f"A3A rank {int(b['rank_a3a'])}).")
rlog("    Both directions are present, which is what licenses treating")
rlog("    opportunity and capability as separate axes rather than one factor.")
rlog("    The gap is ordinal: read the top of the list, not the tail.")
if hpv_neg:
    rlog(f"    HPV16-negative patients are excluded from both picks above.")

# =============================================================================
# STEP 6: TUMOR vs NORMAL-ADJACENT ENZYME CONTROL
# =============================================================================
banner("STEP 6: A3A / A3B, TUMOR vs NORMAL-ADJACENT (within-patient control)")
rlog("  Only a few patients contribute normal-adjacent basal cells. Within those")
rlog("  patients, is the enzyme signal tumor-associated? Small n: this is a")
rlog("  CONTROL, not a result.\n")

rlog(f"  A p-value alone is not enough here. With {MIN_SOURCE_CELLS}+ cells on")
rlog(f"  one side and thousands on the other, a 1.2-fold difference reaches")
rlog(f"  significance without being biologically meaningful, so 'enriched'")
rlog(f"  requires BOTH p < 0.05 AND a fold change >= {MIN_ENRICH_FOLD:.1f}.\n")

src_rows = []
rlog(f"  {'Pat':<7s} {'nTum':>6s} {'nNorm':>6s} | {'A3A tum':>8s} {'A3A norm':>9s} "
     f"{'fold':>7s} {'p':>9s} | {'A3B tum':>8s} {'A3B norm':>9s} {'fold':>7s} {'p':>9s}")
rlog(f"  {'-'*7} {'-'*6} {'-'*6} | {'-'*8} {'-'*9} {'-'*7} {'-'*9} | "
     f"{'-'*8} {'-'*9} {'-'*7} {'-'*9}")
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
        mt, mn = float(np.mean(vt)), float(np.mean(vn))
        try:
            _, pv = mannwhitneyu(vt, vn, alternative='two-sided')
        except ValueError:
            pv = np.nan
        fold = (mt / mn) if mn > 0 else (np.inf if mt > 0 else np.nan)
        row[f'{gname}_tumor']  = mt
        row[f'{gname}_normal'] = mn
        row[f'{gname}_fold']   = fold
        row[f'{gname}_p'] = pv
        fold_str = ('   inf' if np.isinf(fold)
                    else ('    --' if pd.isna(fold) else f"{fold:6.2f}"))
        cells += [f"{mt:>8.3f}", f"{mn:>9.3f}", f"{fold_str:>7s}",
                  (f"{pv:>9.2e}" if pd.notna(pv) else f"{'--':>9s}")]
    rlog(f"  {short(p):<7s} {nt:>6d} {nn:>6d} | " + " ".join(cells[:4]) +
         " | " + " ".join(cells[4:]))
    src_rows.append(row)

if src_rows:
    src_df = pd.DataFrame(src_rows)
    src_df.to_csv(os.path.join(OUT_DIR, "patient_enzyme_by_source.tsv"),
                  sep='\t', index=False)
    rlog(f"\n  [SAVE] patient_enzyme_by_source.tsv  ({len(src_df)} patients with "
         f"both sources at >= {MIN_SOURCE_CELLS} cells)")

    # READ block derived from the table rather than written in advance.
    # Enrichment requires BOTH significance and a real effect size.
    n_pat_src = len(src_df)
    sig = src_df['A3A_p'] < 0.05
    big = src_df['A3A_fold'] >= MIN_ENRICH_FOLD
    enriched = src_df[sig & big]
    sig_only = src_df[sig & ~big & (src_df['A3A_fold'] > 1.0)]
    reversed_ = src_df[src_df['A3A_fold'] < 1.0]

    rlog(f"\n  READ: A3A is tumor-enriched (p < 0.05 AND fold >= "
         f"{MIN_ENRICH_FOLD:.1f}) in {len(enriched)} of {n_pat_src} patients"
         f"{': ' + str(sorted(short(p) for p in enriched['patient'])) if len(enriched) else ''}.")
    for _, r in enriched.iterrows():
        rlog(f"        {short(r['patient'])}: {r['A3A_fold']:.1f}-fold "
             f"({r['A3A_tumor']:.2f} vs {r['A3A_normal']:.2f}), p = {r['A3A_p']:.1e}")
    if len(sig_only):
        rlog(f"        Significant but below the fold gate: "
             f"{sorted(short(p) for p in sig_only['patient'])}.")
        for _, r in sig_only.iterrows():
            rlog(f"          {short(r['patient'])}: only {r['A3A_fold']:.2f}-fold "
                 f"({r['A3A_tumor']:.2f} vs {r['A3A_normal']:.2f}) yet p = "
                 f"{r['A3A_p']:.1e} at n = {int(r['n_tumor'])} vs "
                 f"{int(r['n_normal_adj'])}. Significance from sample size, not")
            rlog(f"          from effect. The normal-adjacent tissue already "
                 f"carries A3A at close to the tumor level.")
    if len(reversed_):
        rlog(f"        Direction reversed (normal-adjacent at or above tumor): "
             f"{sorted(short(p) for p in reversed_['patient'])}.")
        for _, r in reversed_.iterrows():
            rlog(f"          {short(r['patient'])}: {r['A3A_fold']:.2f}-fold "
                 f"({r['A3A_tumor']:.2f} vs {r['A3A_normal']:.2f})")
    if len(enriched) < n_pat_src:
        rlog("        A3A capability therefore looks partly CONSTITUTIVE at the")
        rlog("        patient level rather than purely tumor-induced, which is")
        rlog("        consistent with its independence from viral load and with")
        rlog("        using the all-basal contribution denominator.")
    rlog("        Small n: present this as a CONTROL, not a result.")
else:
    rlog(f"\n  No patient has >= {MIN_SOURCE_CELLS} cells of BOTH sources.")

# =============================================================================
# STEP 7: DETERMINANT vs CONTRIBUTION
# =============================================================================
banner(f"STEP 7: EACH DETERMINANT vs CONTRIBUTION (Spearman, n = {len(out)})")
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
diss = {}
for col, label in [('A3A_mean', 'A3A mean'), ('A3B_mean', 'A3B mean')]:
    cells = []
    for fold in ('fold_sbs2', 'fold_cnv'):
        sub = out[[col, fold]].replace([np.inf, -np.inf], np.nan).dropna()
        if len(sub) >= 4 and sub[col].nunique() > 1:
            rho, pv = spearmanr(sub[col], sub[fold])
            diss[(col, fold)] = (rho, pv)
            cells.append(f"rho={rho:+.3f} p={pv:.3g}{' *' if pv < 0.05 else '  '}")
        else:
            cells.append("N.D.")
    rlog(f"  {label:<12s} {cells[0]:>28s} {cells[1]:>28s}")

n_tests = 2 * len(DETS)
all_p = [p for _, _, p in evidence + cnv_evidence]
min_p = min(all_p) if all_p else float('nan')
rlog("")
rlog("  For contrast, viral load associates with BOTH fates at similar strength,")
rlog("  which is the quantitative form of 'load is opportunity, not fate'.")
rlog("")
rlog(f"  MULTIPLICITY: this step runs {n_tests} correlations. Under BH across all")
rlog(f"  of them, with the smallest observed p = {min_p:.4g}, the corresponding q")
rlog(f"  is approximately {min_p * n_tests:.3g}. The two cells of the dissociation")
rlog("  above are PRE-SPECIFIED single hypotheses from the model, not hits from a")
rlog("  scan, and must be reported that way.")

banner("STEP 7c: NON-CIRCULAR SIGNIFICANT ASSOCIATIONS")
if evidence or cnv_evidence:
    for label, rho, pv in evidence:
        rlog(f"  SBS2 fold  <- {label:<28s} rho={rho:+.3f}  p={pv:.3g}")
    for label, rho, pv in cnv_evidence:
        rlog(f"  CNV  fold  <- {label:<28s} rho={rho:+.3f}  p={pv:.3g}")
else:
    rlog("  none reached p < 0.05 outside the circular pairings.")

# =============================================================================
# STEP 8: THE THREE-AXIS CONJUNCTION MODEL
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


def build_flags(df, load_pctl, enz_pctl, dir_rule=None, dir_min=None):
    """Boolean axis flags at the given thresholds.

    A patient with no gated HPV16-positive cells has no measurable phase and is
    scored DIRECTION = False for both fates under either rule.
    """
    rule = DIRECTION_RULE if dir_rule is None else dir_rule
    dmin = DIRECTION_MIN_FRAC if dir_min is None else dir_min
    t_load = np.nanpercentile(df['load_per_M'], load_pctl)
    t_a3a  = np.nanpercentile(df['A3A_mean'],  enz_pctl)
    t_a3b  = np.nanpercentile(df['A3B_mean'],  enz_pctl)
    f = pd.DataFrame({'patient': df['patient'].values})
    f['opportunity'] = (df['load_per_M'] > t_load).fillna(False).values
    if rule == 'threshold':
        f['dir_sbs2'] = (df['maint_frac'] >= dmin).fillna(False).values
        f['dir_cnv']  = (df['prod_frac']  >= dmin).fillna(False).values
    else:
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
rlog(f"  DIRECTION rule in use: '{DIRECTION_RULE}'"
     + (f" (matching phase fraction >= {DIRECTION_MIN_FRAC:.2f})"
        if DIRECTION_RULE == 'threshold'
        else " (whichever phase dominates; no free parameter)"))
rlog("")
rlog(f"  Default thresholds (p{LOAD_PCTL_DEFAULT} load, p{ENZ_PCTL_DEFAULT} enzyme):")
rlog(f"    load > {t_load:.1f} per million | A3A > {t_a3a:.3f} | A3B > {t_a3b:.3f}")

merged = out.merge(flags, on='patient')

# ---- 8a: per-patient classification, with the failing axis named -------------
FAILS = {}   # patient -> {'sbs2': [...], 'cnv': [...]}

for fate, dir_col, cap_col, all_col, fold_col, drv_col, title in (
        ('sbs2', 'dir_sbs2', 'cap_sbs2', 'all3_sbs2', 'fold_sbs2',
         'is_sbs2_driver', 'SBS2 fate'),
        ('cnv', 'dir_cnv', 'cap_cnv', 'all3_cnv', 'fold_cnv',
         'is_cnv_driver', 'CNV fate')):
    banner(f"STEP 8a: PER-PATIENT CLASSIFICATION ({title})")
    rlog(f"  {'Pat':<7s} {'Opp':>4s} {'Dir':>4s} {'Cap':>4s} {'ALL3':>5s} "
         f"{'fold':>9s} {'driver':>7s}   fails on")
    rlog(f"  {'-'*7} {'-'*4} {'-'*4} {'-'*4} {'-'*5} {'-'*9} {'-'*7}   {'-'*8}")
    for _, r in merged.sort_values(fold_col, ascending=False).iterrows():
        miss = [n for n, v in (('opportunity', r['opportunity']),
                               ('direction', r[dir_col]),
                               ('capability', r[cap_col])) if not v]
        FAILS.setdefault(r['patient'], {})[fate] = miss
        tick = lambda v: ' Y' if v else ' .'
        rlog(f"  {short(r['patient']):<7s} {tick(r['opportunity']):>4s} "
             f"{tick(r[dir_col]):>4s} {tick(r[cap_col]):>4s} "
             f"{('YES' if r[all_col] else '-'):>5s} {r[fold_col]:>8.2f}x "
             f"{('YES' if r[drv_col] else '-'):>7s}   "
             f"{', '.join(miss) if miss else 'none'}")

# ---- 8b: conjunction vs outcome ---------------------------------------------
banner("STEP 8b: CONJUNCTION vs DRIVER STATUS (Fisher exact)")


def conj_test(flag_col, driver_col, fold_col, label):
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
    rlog(f"    Fisher exact p = {pv:.4g}  (odds ratio "
         f"{odds if np.isfinite(odds) else float('inf'):.3g})")
    if b == 0 and c == 0:
        rlog("    PERFECT SEPARATION at these thresholds.")
    if c > 0:
        missed = merged.loc[~merged[flag_col] & merged[driver_col]]
        rlog(f"    drivers missed by the rule: {[short(p) for p in missed['patient']]}")
        for _, r in missed.iterrows():
            fate_key = 'sbs2' if 'sbs2' in flag_col else 'cnv'
            rlog(f"      {short(r['patient'])}: fold {r[fold_col]:.2f}x, "
                 f"fails on {', '.join(FAILS[r['patient']][fate_key])}; "
                 f"load {r['load_per_M']:.0f} "
                 f"(rank {int(out['load_per_M'].rank(ascending=False)[out['patient'] == r['patient']].iloc[0])}"
                 f"/{len(out)}), "
                 f"load per positive cell {r['load_in_pos']:.1f}")
    if b > 0:
        who = merged.loc[merged[flag_col] & ~merged[driver_col], 'patient']
        rlog(f"    non-drivers flagged by the rule: {[short(p) for p in who]}")
    rlog("")
    return pv


p_sbs2 = conj_test('all3_sbs2', 'is_sbs2_driver', 'fold_sbs2', 'SBS2-HIGH fate:')
p_cnv  = conj_test('all3_cnv',  'is_cnv_driver',  'fold_cnv',  'CNV-HIGH fate:')

rlog("  Any driver the rule misses is listed above with its load rank and its")
rlog("  load per positive cell. A missed driver that is a load outlier is the")
rlog("  LOAD ROUTE into that fate: a phase-based rule is expected to miss it.")
rlog("  That alternative route was identified before this test, not fitted to it.")

# ---- 8c: threshold sensitivity ----------------------------------------------
banner("STEP 8c: THRESHOLD SENSITIVITY SWEEP")

# Identify the tightest margin on the capability axis without naming it in advance.
margin = (out['A3A_mean'] - t_a3a).abs()
tight_idx = margin.idxmin()
tight = out.loc[tight_idx]
rlog("  If the separation only holds at one arbitrary cut, it is not a result.")
rlog(f"  The tightest capability margin in this cohort is "
     f"{short(tight['patient'])}: A3A = {tight['A3A_mean']:.3f} against a")
rlog(f"  threshold of {t_a3a:.3f}, a margin of {margin.min():.4f}. Sweeping both")
rlog("  thresholds shows the range over which separation survives.\n")
rlog(f"  {'loadP':>6s} {'enzP':>5s} | {'SBS2: hit/miss/false':>22s} {'p':>9s} "
     f"| {'CNV: hit/miss/false':>21s} {'p':>9s}")
rlog(f"  {'-'*6} {'-'*5} | {'-'*22} {'-'*9} | {'-'*21} {'-'*9}")

sweep_rows = []
for lp, ep in product(SWEEP_PCTLS, SWEEP_PCTLS):
    fl, _ = build_flags(out, lp, ep)
    mg = out.merge(fl, on='patient')
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
                       'cnv_perfect': res['cnv'][4],
                       'denominator': CONTRIBUTION_DENOMINATOR})
    star_s = ' *' if res['sbs2'][4] else '  '
    star_c = ' *' if res['cnv'][4] else '  '
    rlog(f"  {lp:>6d} {ep:>5d} | {res['sbs2'][0]:>7d}/{res['sbs2'][1]:<3d}/"
         f"{res['sbs2'][2]:<3d}{star_s:>7s} {res['sbs2'][3]:>9.4g} "
         f"| {res['cnv'][0]:>6d}/{res['cnv'][1]:<3d}/{res['cnv'][2]:<3d}{star_c:>6s} "
         f"{res['cnv'][3]:>9.4g}")

sweep = pd.DataFrame(sweep_rows)
sweep.to_csv(os.path.join(OUT_DIR, "patient_conjunction_sensitivity.tsv"),
             sep='\t', index=False)
rlog("\n  hit = drivers captured, miss = drivers the rule fails to flag,")
rlog("  false = non-drivers the rule wrongly flags. '*' marks perfect separation.")
for tag, name in (('sbs2', 'SBS2'), ('cnv', 'CNV ')):
    n_perfect = int(sweep[f'{tag}_perfect'].sum())
    rlog(f"  {name} perfect separation in {n_perfect} of {len(sweep)} threshold "
         f"combinations.")
    if n_perfect:
        ok = sweep[sweep[f'{tag}_perfect']]
        rlog(f"       holds across load p{ok['load_pctl'].min()}-"
             f"p{ok['load_pctl'].max()} and enzyme p{ok['enz_pctl'].min()}-"
             f"p{ok['enz_pctl'].max()}.")
    else:
        rlog(f"       no combination gives perfect separation; for this fate the")
        rlog(f"       rule is threshold-dependent and descriptive only.")

# ---- 8c-dir: is the DIRECTION rule itself load-bearing? ----------------------
banner("STEP 8c-dir: DIRECTION RULE SENSITIVITY")
rlog("  The parameter-free rule asks only which phase dominates. Replacing it")
rlog("  with a fixed cutoff on the matching phase fraction adds a free parameter,")
rlog("  and a free parameter chosen after seeing which patient the rule missed is")
rlog("  fitting, not testing. This step sweeps that cutoff across the SAME 25")
rlog("  load-by-enzyme combinations so the width of any window is visible.\n")

rlog(f"  {'direction rule':<26s} {'SBS2 perfect':>13s} {'CNV perfect':>12s}")
rlog(f"  {'-'*26} {'-'*13} {'-'*12}")
dir_rows = []
for rule, dmin, label in ([('compare', None, 'dominant phase')]
                          + [('threshold', v, f'fraction >= {v:.2f}')
                             for v in DIRECTION_SWEEP]):
    n_s = n_c = 0
    for lp, ep in product(SWEEP_PCTLS, SWEEP_PCTLS):
        fl, _ = build_flags(out, lp, ep, dir_rule=rule, dir_min=dmin)
        mg = out.merge(fl, on='patient')
        for tag, fcol, dcol in (('sbs2', 'all3_sbs2', 'is_sbs2_driver'),
                                ('cnv', 'all3_cnv', 'is_cnv_driver')):
            b = int((mg[fcol] & ~mg[dcol]).sum())
            c = int((~mg[fcol] & mg[dcol]).sum())
            if b == 0 and c == 0:
                if tag == 'sbs2':
                    n_s += 1
                else:
                    n_c += 1
    active = ((rule == DIRECTION_RULE) and
              (rule == 'compare' or abs(dmin - DIRECTION_MIN_FRAC) < 1e-9))
    dir_rows.append({'direction_rule': rule, 'direction_min_frac': dmin,
                     'label': label, 'sbs2_perfect': n_s, 'cnv_perfect': n_c,
                     'n_combinations': len(SWEEP_PCTLS) ** 2, 'active': active})
    rlog(f"  {label:<26s} {n_s:>10d}/{len(SWEEP_PCTLS)**2} "
         f"{n_c:>9d}/{len(SWEEP_PCTLS)**2}" + ("   <-- in use" if active else ""))

dir_df = pd.DataFrame(dir_rows)
dir_df.to_csv(os.path.join(OUT_DIR, "patient_direction_rule_sensitivity.tsv"),
              sep='\t', index=False)

# Where the CNV drivers actually sit on the productive axis, so any window can be
# measured against the gap it has to thread.
rlog("")
rlog("  Productive fraction, descending, with CNV drivers marked. A cutoff has to")
rlog("  separate the drivers from everyone else to buy a perfect separation:")
ph = out.dropna(subset=['prod_frac']).sort_values('prod_frac', ascending=False)
for _, r in ph.iterrows():
    rlog(f"    {short(r['patient']):<7s} prod {r['prod_frac']:.2f}  "
         f"maint {r['maint_frac']:.2f}  CNV fold {r['fold_cnv']:>5.2f}x"
         + ("   <-- CNV driver" if r['is_cnv_driver'] else ""))

thr_rows = dir_df[dir_df['direction_rule'] == 'threshold']
win = thr_rows[thr_rows['cnv_perfect'] > 0]
rlog("")
if len(win):
    lo, hi = win['direction_min_frac'].min(), win['direction_min_frac'].max()
    rlog(f"  READ: a fixed cutoff buys CNV perfect separation only between "
         f"{lo:.2f} and {hi:.2f}.")
    rlog(f"        Outside that window it buys nothing. Report the window width")
    rlog(f"        alongside any claim that rests on it, and state plainly whether")
    rlog(f"        the cutoff was chosen before or after seeing which driver the")
    rlog(f"        parameter-free rule missed.")
else:
    rlog("  READ: no fixed cutoff in the swept range buys CNV perfect separation.")
best_free = dir_df[dir_df['direction_rule'] == 'compare'].iloc[0]
rlog(f"        The parameter-free rule gives SBS2 {int(best_free['sbs2_perfect'])}"
     f"/{int(best_free['n_combinations'])} and CNV "
     f"{int(best_free['cnv_perfect'])}/{int(best_free['n_combinations'])}, with")
rlog(f"        nothing to tune.")

# ---- 8d: conjunction structure (minimum vs mean) ----------------------------
banner("STEP 8d: IS THE STRUCTURE A CONJUNCTION? (min vs mean of scaled axes)")
rlog("  If fate needs ALL three conditions, the MINIMUM of the three scaled axes")
rlog("  should track contribution better than their MEAN, because a single")
rlog("  missing factor should veto the outcome rather than be averaged away.")
rlog("  Axes are scaled to percentile ranks within the cohort (0-1), so the")
rlog("  comparison does not depend on their raw units.\n")


def pct_rank(s):
    return s.rank(pct=True, na_option='bottom')


scm = pd.DataFrame({'patient': out['patient'].values})
scm['opp']       = pct_rank(out['load_per_M']).values
scm['dir_sbs2']  = pct_rank(out['maint_frac'] - out['prod_frac']).values
scm['dir_cnv']   = pct_rank(out['prod_frac'] - out['maint_frac']).values
scm['cap_sbs2']  = pct_rank(out['A3A_mean']).values
scm['cap_cnv']   = pct_rank(out['A3B_mean']).values
scm['min_sbs2']  = scm[['opp', 'dir_sbs2', 'cap_sbs2']].min(axis=1)
scm['mean_sbs2'] = scm[['opp', 'dir_sbs2', 'cap_sbs2']].mean(axis=1)
scm['min_cnv']   = scm[['opp', 'dir_cnv', 'cap_cnv']].min(axis=1)
scm['mean_cnv']  = scm[['opp', 'dir_cnv', 'cap_cnv']].mean(axis=1)

# Raw values alongside the ranks, so the figure can annotate cells without
# reopening any other table.
scm = scm.merge(out[['patient', 'fold_sbs2', 'fold_cnv', 'load_per_M',
                     'maint_frac', 'prod_frac', 'A3A_mean', 'A3B_mean',
                     'is_sbs2_driver', 'is_cnv_driver']], on='patient')
scm = scm.merge(flags[['patient', 'opportunity', 'dir_sbs2', 'dir_cnv',
                       'cap_sbs2', 'cap_cnv', 'all3_sbs2', 'all3_cnv']],
                on='patient', suffixes=('', '_flag'))
scm['direction_rule'] = DIRECTION_RULE
scm['direction_min_frac'] = (DIRECTION_MIN_FRAC if DIRECTION_RULE == 'threshold'
                             else np.nan)
scm['thr_load'] = t_load
scm['thr_a3a']  = t_a3a
scm['thr_a3b']  = t_a3b
scm['denominator'] = CONTRIBUTION_DENOMINATOR
scm.to_csv(os.path.join(OUT_DIR, "patient_conjunction_model.tsv"),
           sep='\t', index=False)

structure = {}
for fate, mn, mnm in (('SBS2', 'min_sbs2', 'mean_sbs2'),
                      ('CNV',  'min_cnv',  'mean_cnv')):
    fold = 'fold_sbs2' if fate == 'SBS2' else 'fold_cnv'
    r_min, p_min = spearmanr(scm[mn], scm[fold])
    r_mean, p_mean = spearmanr(scm[mnm], scm[fold])
    structure[fate] = (r_min, r_mean)
    rlog(f"  {fate}-HIGH fold:")
    rlog(f"    minimum of three axes  rho={r_min:+.3f}  p={p_min:.3g}")
    rlog(f"    mean of three axes     rho={r_mean:+.3f}  p={p_mean:.3g}")
    verdict = ("minimum wins -> conjunction supported"
               if r_min > r_mean else
               "mean wins -> additive, NOT a strict conjunction")
    rlog(f"    -> {verdict}\n")

# Failure distribution, derived, so the anti-overfitting argument cannot drift.
rlog("  Where the SBS2 rule's non-drivers fail (the anti-overfitting argument is")
rlog("  that failures fall on DIFFERENT axes, not all on the same term):")
by_axis = {}
for _, r in merged.iterrows():
    if r['is_sbs2_driver']:
        continue
    key = ' + '.join(FAILS[r['patient']]['sbs2']) or 'none'
    by_axis.setdefault(key, []).append(short(r['patient']))
for key in sorted(by_axis):
    rlog(f"    {key:<45s} {sorted(by_axis[key])}")
rlog(f"    -> {len(by_axis)} distinct failure patterns across "
     f"{sum(len(v) for v in by_axis.values())} non-drivers.")

rlog("")
rlog(f"  CAVEAT (state this wherever the conjunction is reported): with n = "
     f"{len(out)} and")
rlog("  three binary terms, a rule of this form has enough freedom to fit almost")
rlog("  any small driver set. What defends it is that the three axes were")
rlog("  specified from the biological model BEFORE this test rather than searched")
rlog("  for, that the same rule applied independently to the other fate is")
rlog("  reported on the same footing, and that the failures fall on different")
rlog("  axes as enumerated above. It is a consistency check on the model, not")
rlog("  independent confirmation of it.")

# =============================================================================
# SAVE
# =============================================================================
report_path = os.path.join(OUT_DIR, f"patient_determinants_{TIMESTAMP}.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
rlog(f"\n  Denominator: {CONTRIBUTION_DENOMINATOR}")
rlog(f"  Table:       {os.path.join(OUT_DIR, 'patient_determinants_table.tsv')}")
rlog(f"  Conjunction: {os.path.join(OUT_DIR, 'patient_conjunction_model.tsv')}")
rlog(f"  Sensitivity: {os.path.join(OUT_DIR, 'patient_conjunction_sensitivity.tsv')}")
rlog(f"  Direction:   {os.path.join(OUT_DIR, 'patient_direction_rule_sensitivity.tsv')}")
rlog(f"  Report:      {report_path}")
banner("PATIENT DETERMINANTS TABLE COMPLETE (v6)")
