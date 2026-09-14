#!/usr/bin/env python3
"""
Diagnostic_Fig3_Dominance_Statistics.py   (item 7.5)  -- v3
================================================================================
READ-ONLY. Does A3 dominance actually predict SBS2 burden and CNA score, once
patient structure is accounted for?

v3 CHANGES
----------
  A. NEGATIVE CONTROL NOW JUDGES ON RAW p, NOT q. BH correction on a negative
     control works backwards: it makes it EASIER to declare the control clean,
     which is the opposite of conservative. v2 printed "All null" for
     normal-adjacent while model (c) sat at raw p = 0.0398 (q = 0.064). Any
     normal-adjacent model at raw p < 0.05 is now flagged, and the tumour-vs-
     normal DIRECTION is printed, because an inverted sign is reassuring
     (not a leaked copy of the tumour effect) while a matching sign is fatal.
     GEE rows also now print as "OR=0.675" rather than "+0.6745", which read
     like a signed coefficient.

  B. Non-evaluable patient-level rows print "not evaluable" instead of "+nan".

  C. TWO UNIVERSES, because the figure script uses two.
     Generate_Figure3_A3_Dominance_Slider.py computes its rank z-scores and its
     frac-vs-outcome correlations on EVERY cell expressing either enzyme (ties
     INCLUDED), and only splits by dominance afterwards. v2 dropped ties first
     and correlated on what was left, so the rho anchors checked a universe 1014
     cells smaller than the figure's. The figure script header records
     frac-vs-CNV = -0.437 (tumour), which is what the manuscript rounds to -0.44;
     v2 returned -0.445. Both round the same, so the anchor passed on tolerance,
     but an anchor that is only approximately right is worth little.
     Counts EXCLUDE ties. Correlations and rank z-scores INCLUDE them.
     Anchor tolerance tightened from 0.015 to 0.006 now that it is exact.

  D. THE COMPOSITE AXIS IS NOW TESTED. Figure 3 plots z(SBS2) - z(CNV), and its
     own header notes that axis is carried mainly by CNV. The patient-level test
     ran on SBS2 prevalence and CNA mean SEPARATELY, so overlaying it on the
     figure would invite a reader to check an agreement never measured. The
     composite axis is added as a patient-level outcome, which decides whether
     per-patient marks can go on the plot at all.

v2 changes retained: ties excluded from dominance (Methods 6.2 uses strict
inequalities); CNA patient summary is the mean, not the degenerate median;
verdict reads the mean-based CNA test; p-value underflow prints as "< 1e-300".

WHY THIS EXISTS
---------------
Results 4.1 concludes:

    "A3A dominance is associated with SBS2, whereas A3B dominance is
     associated with CNA."

The A3B-to-CNA half is well supported (A3A fraction vs CNA score, rho = -0.44 in
tumour). The A3A-to-SBS2 half rests on rho = +0.05, reported two sentences
earlier in the same paragraph. A methods-minded Cancer Discovery reviewer will
land there.

Three reasons the raw Spearman is the wrong statistic here:
  1. Pseudoreplication. ~24,000 cells from 14 patients, and the SBS2-HIGH
     population is 74% drawn from three of them.
  2. Zero inflation. SBS2 is zero in ~84% of these cells, which drags a pooled
     correlation toward zero regardless of whether a real difference exists
     among the cells that do carry it.
  3. Continuous A3 expression is a noisy proxy for cumulative activity
     (episodic induction), which is why Figure 3 uses dominance, not level.

WHAT IT RUNS
------------
  HEADLINE  Patient-level paired Wilcoxon signed-rank, n <= 14. Per patient, a
            summary of the outcome in A3A-dominant cells vs the same summary in
            A3B-dominant cells; one pair per patient. No pseudoreplication, and
            it answers the question a reviewer is actually asking: does this
            hold WITHIN tumours, or is it carried by three of them. Reported
            with a Hodges-Lehmann estimate and a bootstrap CI.

            SUMMARY STATISTIC BY OUTCOME, and this matters more than it looks.
            A median suits NEITHER raw outcome: SBS2 is zero in ~84% of cells so
            the median is 0 in both arms for every patient and the signed-rank
            test is undefined; cnv_score is too coarse for per-patient medians to
            separate, so 10 of 13 paired differences come out exactly zero,
            wilcoxon discards them, and the test collapses onto 3 pairs while
            printing a zero-width CI that reads like a null.
            SBS2 -> PREVALENCE (primary) and MEAN (magnitude).
            CNA  -> MEAN. The median is kept, explicitly labelled degenerate,
                    so the failure stays visible rather than being hidden.
            COMPOSITE AXIS -> MEDIAN, which is appropriate there: it is a
                    continuous rank-based score with neither zero inflation nor
                    discretisation.

  SUPPORT   Cell-level models, patient as a random effect:
            (a) CNA        -> LMM, no transformation
            (b) SBS2 rank  -> LMM on the within-tissue rank z-score
            (c) SBS2 > 0   -> GEE logistic, exchangeable, clustered by patient.
                              Population-averaged, cluster-robust SE.
            (d) SBS2 | >0  -> LMM among carriers only.
            All four run in TUMOUR and in NORMAL-ADJACENT (negative control).

DESIGN DECISIONS, and why
-------------------------
  Predictor is BINARY DOMINANCE, the same cut Figure 3 draws, so the test and
  the figure describe one thing. Continuous A3A fraction runs as a sensitivity
  analysis only.

  Covariate is log total UMI: both SBS2 detection (callable sites) and CNA
  scoring depend on sequencing depth.

  CELL-CYCLE PHASE IS DELIBERATELY NOT ADJUSTED FOR. A3B expression peaks at
  G2/M, so phase lies on the pathway under study. It is a mediator, not a
  confounder, and adjusting would be over-control.

  EFFECT SIZES AND CIs BEFORE p-VALUES. At n ~ 24,000 the LMMs return p-values
  near zero for effects of no consequence.

  A random slope for dominance by patient is attempted. With 14 patients it
  frequently converges while emitting singular-covariance or boundary warnings,
  meaning the between-patient slope variance is near zero. Treat convergence as
  "not evidence against" consistency, never as evidence for it.

ANCHORS
-------
Reproduces the Figure 3 / Results 4.1 numbers before running anything new. This
check earned its keep twice: it caught the v1 tie bug (only the A3B arms were
inflated while A3A matched exactly, which pinned the cause immediately) and the
v2 universe mismatch.

INPUT (read-only)
-----------------
  data/FIG_6/01_raw_hpv16_counts/basal_cell_master_table_with_raw_HPV16.tsv

OUTPUT (to data/FIG_3/DIAGNOSTIC_DOMINANCE_STATS/)
--------------------------------------------------
  dominance_statistics_report.txt
  patient_level_medians.tsv
  patient_level_tests.tsv
  cell_level_models.tsv

Env: NETWORK
Usage: conda run -n NETWORK python Diagnostic_Fig3_Dominance_Statistics.py

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
import statsmodels.api as sm

warnings.filterwarnings('ignore')

# =============================================================================
# CONFIG
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
MASTER_TABLE_PATH = os.path.join(PROJECT_ROOT,
    "data/FIG_6/01_raw_hpv16_counts/basal_cell_master_table_with_raw_HPV16.tsv")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_3/DIAGNOSTIC_DOMINANCE_STATS")
os.makedirs(OUTPUT_DIR, exist_ok=True)

DOMINANCE_CUT = 0.5        # A3A/(A3A+A3B); matches Figure 3. Ties are EXCLUDED.
MIN_CELLS_PER_ARM = 10     # a patient needs this many cells in BOTH arms to pair
N_BOOT = 10000
BOOT_SEED = 42
ANCHOR_TOL_RHO = 0.006     # v3: tightened from 0.015 now that the universe is exact

# Figure 3 / Results 4.1 values this script must reproduce before anything else.
ANCHORS = [
    ('tumour A3A-dominant n',        11252, 'n_a3a_tumour',   'count'),
    ('tumour A3B-dominant n',        12438, 'n_a3b_tumour',   'count'),
    ('normal A3A-dominant n',          119, 'n_a3a_normal',   'count'),
    ('normal A3B-dominant n',           57, 'n_a3b_normal',   'count'),
    # Figure script header records -0.437 / +0.070 on the ties-INCLUDED universe.
    ('tumour rho A3Afrac vs CNA',  -0.437, 'rho_frac_cnv_tumour',  'rho'),
    ('normal rho A3Afrac vs CNA',  +0.070, 'rho_frac_cnv_normal',  'rho'),
    ('tumour rho A3Afrac vs SBS2', +0.047, 'rho_frac_sbs2_tumour', 'rho'),
    ('normal rho A3Afrac vs SBS2', +0.037, 'rho_frac_sbs2_normal', 'rho'),
]

COL_CANDIDATES = {
    'patient': ['subject id', 'subject_id', 'subject', 'patient', 'donor'],
    'tissue':  ['tissue type', 'tissue_type', 'tissue', 'sample_type'],
    'sbs2':    ['sig_SBS2', 'SBS2', 'sbs2_weight'],
    'cnv':     ['cnv_score', 'cnv', 'infercnv_score'],
    'a3a':     ['APOBEC3A', 'A3A'],
    'a3b':     ['APOBEC3B', 'A3B'],
    'umi':     ['n_counts', 'total_counts', 'nCount_RNA'],
}
PHASE_CANDIDATES = ['phase', 'ccAFv2', 'cell_cycle_phase', 'ccafv2_phase']

COMPOSITE_NAME = 'composite axis z(SBS2)-z(CNV)'

# =============================================================================
# LOGGING
# =============================================================================
report = []
def log(msg=""):
    print(msg, flush=True)
    report.append(str(msg))

def banner(t, ch="="):
    log(""); log(ch * 88); log(f"  {t}"); log(ch * 88)

def fmt_p(p):
    """A p-value of exactly 0 is floating-point underflow, not a result, and
    must never reach a manuscript as '0'."""
    if p is None:
        return "N.D."
    try:
        p = float(p)
    except (TypeError, ValueError):
        return "N.D."
    if np.isnan(p):
        return "N.D."
    if p == 0:
        return "< 1e-300 (underflow)"
    return f"{p:.3g}"

def _first(cols, cands):
    low = {str(c).strip().lower(): c for c in cols}
    for c in cands:
        if c.strip().lower() in low:
            return low[c.strip().lower()]
    return None

# =============================================================================
# STEP 0: load and resolve columns
# =============================================================================
banner("STEP 0: Load master table and resolve columns")

if not os.path.exists(MASTER_TABLE_PATH):
    sys.exit(f"ERROR: master table not found: {MASTER_TABLE_PATH}")

df = pd.read_csv(MASTER_TABLE_PATH, sep='\t', index_col=0)
log(f"  master table: {df.shape[0]} basal cells x {df.shape[1]} columns")

cols = {}
for key, cands in COL_CANDIDATES.items():
    c = _first(df.columns, cands)
    if c is None:
        log(f"  ERROR: no column for '{key}'. Tried {cands}")
        log(f"  available: {list(df.columns)}")
        sys.exit(1)
    cols[key] = c
    log(f"  {key:<8s} -> '{c}'")

phase_col = _first(df.columns, PHASE_CANDIDATES)
log(f"  phase    -> {phase_col or 'NOT PRESENT (descriptive block skipped)'}")

# --- tissue mapping ---
raw_tissue = df[cols['tissue']].astype(str).str.strip()
log(f"\n  tissue values: {sorted(raw_tissue.unique())}")

def map_tissue(v):
    s = str(v).lower()
    if 'tum' in s or 'canc' in s:
        return 'tumour'
    if 'norm' in s or 'adj' in s or 'health' in s:
        return 'normal'
    return None

df['_tissue'] = raw_tissue.map(map_tissue)
if df['_tissue'].isna().any():
    bad = sorted(raw_tissue[df['_tissue'].isna()].unique())
    log(f"  ERROR: unmapped tissue labels {bad}. Extend map_tissue() and re-run.")
    sys.exit(1)
for t in ('tumour', 'normal'):
    log(f"    {t}: {(df['_tissue'] == t).sum()} cells")

# --- numeric coercion ---
for k in ('sbs2', 'cnv', 'a3a', 'a3b', 'umi'):
    df[f'_{k}'] = pd.to_numeric(df[cols[k]], errors='coerce')
df['_patient'] = df[cols['patient']].astype(str)
log(f"\n  patients: {df['_patient'].nunique()}")

df['_sbs2'] = df['_sbs2'].fillna(0.0)   # no detected SBS2 == zero weight

before = len(df)
df = df[df[['_cnv', '_a3a', '_a3b', '_umi']].notna().all(axis=1)].copy()
if len(df) < before:
    log(f"  dropped {before - len(df)} cells with missing cnv/A3/UMI")

# --- dominance, and the two universes ---------------------------------------
df['_a3sum'] = df['_a3a'] + df['_a3b']

# TWO UNIVERSES, and they are not the same.
# Generate_Figure3_A3_Dominance_Slider.py (process_tissue) computes its rank
# z-scores and its frac-vs-outcome correlations on EVERY cell expressing either
# enzyme, TIES INCLUDED, and only splits by dominance afterwards. So:
#   expr_all -> ties INCLUDED. Used for the rho anchors and for the composite
#               z(SBS2) - z(CNV) axis, so both match the figure exactly.
#   expr     -> ties EXCLUDED. Used for everything keyed on dominance.
expr_all = df[df['_a3sum'] > 0].copy()
expr_all['_a3a_frac'] = expr_all['_a3a'] / expr_all['_a3sum']
log(f"  A3-expressing basal cells (A3A + A3B > 0): {len(expr_all)} of {len(df)}")

# Rank z-scores and the composite axis, WITHIN tissue on the ties-included
# universe, exactly as the figure script's zrank() does (average ranks, ddof=0).
for _t in ('tumour', 'normal'):
    _m = expr_all['_tissue'] == _t
    for _src, _dst in (('_sbs2', '_sbs2_rankz'), ('_cnv', '_cnv_rankz')):
        _r = expr_all.loc[_m, _src].rank()
        _s = _r.std(ddof=0)
        expr_all.loc[_m, _dst] = (_r - _r.mean()) / _s if _s > 0 else 0.0
expr_all['_x'] = expr_all['_sbs2_rankz'] - expr_all['_cnv_rankz']

# Methods 6.2 classifies on STRICT inequalities: above 0.5 is A3A-dominant,
# below 0.5 is A3B-dominant. Cells at EXACTLY 0.5 belong to neither arm.
# This is not a rounding detail -- a cell with one UMI of A3A and one of A3B has
# identical log-normalized values and lands exactly on the cut, so ties are ~4%
# of expressing cells. Sweeping them into A3B (the v1 bug) inflated that arm by
# 1010 cells in tumour and 4 in normal while leaving the A3A arms untouched.
# That asymmetry is what made the bug diagnosable from the anchor table.
expr = expr_all.copy()
_n_tie = int((expr['_a3a_frac'] == DOMINANCE_CUT).sum())
expr = expr[expr['_a3a_frac'] != DOMINANCE_CUT].copy()
expr['_dominant'] = np.where(expr['_a3a_frac'] > DOMINANCE_CUT, 'A3A', 'A3B')
expr['_dom01'] = (expr['_dominant'] == 'A3A').astype(int)
expr['_log_umi'] = np.log1p(expr['_umi'])
expr['_sbs2_pos'] = (expr['_sbs2'] > 0).astype(int)
log(f"  dropped {_n_tie} cells at exactly {DOMINANCE_CUT} (neither arm), "
    f"leaving {len(expr)}")
log(f"  rank z-scores and the composite axis carried over from the")
log(f"  ties-included universe, so they match the figure cell for cell.")

# =============================================================================
# STEP 1: anchors
# =============================================================================
banner("STEP 1: Reproduce the Figure 3 / Results 4.1 numbers")

COMPUTED = {}
for t in ('tumour', 'normal'):
    sub = expr[expr['_tissue'] == t]                   # ties EXCLUDED
    sub_all = expr_all[expr_all['_tissue'] == t]       # ties INCLUDED
    COMPUTED[f'n_a3a_{t}'] = int((sub['_dominant'] == 'A3A').sum())
    COMPUTED[f'n_a3b_{t}'] = int((sub['_dominant'] == 'A3B').sum())
    # Correlations on the ties-INCLUDED universe, matching process_tissue().
    COMPUTED[f'rho_frac_cnv_{t}'] = spearmanr(sub_all['_a3a_frac'], sub_all['_cnv'])[0]
    COMPUTED[f'rho_frac_sbs2_{t}'] = spearmanr(sub_all['_a3a_frac'], sub_all['_sbs2'])[0]

log(f"  {'anchor':<32s} {'stated':>10s} {'computed':>12s}   verdict")
log(f"  {'-'*32} {'-'*10} {'-'*12}   -------")
n_diff = 0
for label, stated, key, kind in ANCHORS:
    got = COMPUTED.get(key)
    if got is None or (isinstance(got, float) and np.isnan(got)):
        v = 'NO VALUE'
    elif kind == 'count':
        v = 'MATCH' if int(got) == int(stated) else 'DIFF'
    else:
        v = 'MATCH' if abs(got - stated) <= ANCHOR_TOL_RHO else 'DIFF'
    if v != 'MATCH':
        n_diff += 1
    shown = f"{got:.3f}" if kind == 'rho' else f"{got}"
    log(f"  {label:<32s} {stated:>10} {shown:>12s}   {v}")

if n_diff:
    log(f"\n  >>> {n_diff} anchor(s) FAILED. The input or the dominance definition")
    log(f"      has drifted from what Figure 3 used. NOTHING BELOW IS TRUSTWORTHY")
    log(f"      until this is resolved.")
    log(f"      Hint 1: if ONLY the A3B arms are inflated while the A3A arms match")
    log(f"              exactly, the tie rule is wrong.")
    log(f"      Hint 2: if the rho anchors move but the counts hold, check which")
    log(f"              universe the correlation runs on. The figure includes ties")
    log(f"              in the correlation and excludes them from the arms.")
    log(f"      Hint 3: if everything moves, check the SBS2 source column -- the")
    log(f"              master table may carry a different refitting than")
    log(f"              signature_refitting_hnscc/ (which needs .T on load).")
else:
    log(f"\n  All anchors reproduce. Proceeding.")

# zero inflation, the reason a pooled Spearman understates SBS2
log("")
for t in ('tumour', 'normal'):
    sub = expr[expr['_tissue'] == t]
    log(f"  {t}: SBS2 > 0 in {sub['_sbs2_pos'].sum()} / {len(sub)} cells "
        f"({100*sub['_sbs2_pos'].mean():.1f}%)")
log("  >>> This is why the pooled Spearman is near zero for SBS2 whether or not")
log("      a real difference exists among carriers. Hence the hurdle split, and")
log("      hence prevalence rather than the median at the patient level.")

# cell-cycle phase, descriptive ONLY
if phase_col:
    banner("Cell-cycle phase by dominance (DESCRIPTIVE, NOT ADJUSTED FOR)", "-")
    log("  Phase is a MEDIATOR (A3B peaks at G2/M), not a confounder.")
    log("  Adjusting for it would be over-control. Shown so the choice is visible.")
    ct = pd.crosstab(expr['_dominant'], expr[phase_col], normalize='index') * 100
    log("\n" + ct.round(1).to_string())

# =============================================================================
# STEP 2: HEADLINE -- patient-level paired Wilcoxon
# =============================================================================
banner("STEP 2: HEADLINE -- patient-level paired Wilcoxon signed-rank")
log("  One pair per patient. No pseudoreplication. This is the number to lead")
log("  with; the cell-level models below are supporting evidence.")

def hodges_lehmann(d):
    """Pseudomedian: median of all Walsh averages. The location estimate that
    goes with a signed-rank test."""
    d = np.asarray(d, dtype=float)
    n = len(d)
    if n == 0:
        return np.nan
    w = [(d[i] + d[j]) / 2.0 for i in range(n) for j in range(i, n)]
    return float(np.median(w))

def boot_ci(d, stat=np.median, n_boot=N_BOOT, seed=BOOT_SEED):
    d = np.asarray(d, dtype=float)
    if len(d) < 3:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    vals = [stat(rng.choice(d, size=len(d), replace=True)) for _ in range(n_boot)]
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)))

def _prevalence(v):
    return float(np.mean(np.asarray(v) > 0))

# The summary statistic has to suit the outcome, and a median suits NEITHER raw
# outcome here. SBS2 is zero in ~84% of cells, so a per-patient median is 0 in
# both arms for every patient and the signed-rank test is undefined. cnv_score
# is too coarse for per-patient medians to separate, so most paired differences
# come out exactly zero, wilcoxon discards them, and the test collapses onto a
# handful of pairs while printing a zero-width CI that reads like a null.
#   composite axis -> median (continuous rank score; neither failure mode applies)
#   SBS2           -> prevalence (primary), mean (magnitude)
#   CNA            -> mean; the median is retained, explicitly labelled, so the
#                     degeneracy stays visible instead of being quietly dropped.
PT_SPECS = [
    # The composite axis Figure 3 actually draws. Tested so that any claim of
    # within-patient consistency ON THE FIGURE is one we have measured, rather
    # than one inferred from the separate SBS2 and CNA tests.
    ('_x',    np.median,   COMPOSITE_NAME),
    ('_sbs2', _prevalence, 'SBS2 prevalence'),
    ('_sbs2', np.mean,     'SBS2 mean weight'),
    ('_cnv',  np.mean,     'CNA mean score'),
    ('_cnv',  np.median,   'CNA median score (degenerate, secondary)'),
]

pt_rows, test_rows = [], []
for t in ('tumour', 'normal'):
    sub = expr[expr['_tissue'] == t]
    for outcome, summ, oname in PT_SPECS:
        paired, used = [], []
        for pid, g in sub.groupby('_patient'):
            a = g.loc[g['_dominant'] == 'A3A', outcome].values
            b = g.loc[g['_dominant'] == 'A3B', outcome].values
            if len(a) >= MIN_CELLS_PER_ARM and len(b) >= MIN_CELLS_PER_ARM:
                sa, sb = float(summ(a)), float(summ(b))
                paired.append(sa - sb)
                used.append(pid)
                pt_rows.append({'tissue': t, 'outcome': oname, 'patient': pid,
                                'n_A3A': len(a), 'n_A3B': len(b),
                                'summary_A3A': sa, 'summary_A3B': sb,
                                'difference': sa - sb})
        n = len(paired)
        n_zero = int(np.sum(np.isclose(paired, 0.0))) if n else 0
        if n < 3 or np.allclose(paired, 0):
            log(f"\n  {t} / {oname}: only {n} pairable patient(s)"
                f"{' (all differences zero)' if n >= 3 else ''}; test not run.")
            test_rows.append({'tissue': t, 'outcome': oname, 'n_patients': n,
                              'n_zero_diffs': n_zero, 'hodges_lehmann': np.nan,
                              'ci_lo': np.nan, 'ci_hi': np.nan, 'raw_p': np.nan,
                              'n_positive': np.nan})
            continue
        stat, p = wilcoxon(paired, alternative='two-sided')
        hl = hodges_lehmann(paired)
        lo, hi = boot_ci(paired)
        n_pos = int(np.sum(np.array(paired) > 0))
        log(f"\n  {t} / {oname}   (n = {n} patients with >= {MIN_CELLS_PER_ARM} "
            f"cells in both arms)")
        log(f"    A3A-dominant higher in {n_pos}/{n} patients")
        log(f"    Hodges-Lehmann shift = {hl:+.4f}  95% CI {lo:+.4f} to {hi:+.4f}")
        log(f"    Wilcoxon signed-rank raw p = {fmt_p(p)}")
        if n_zero:
            log(f"    NOTE: {n_zero}/{n} paired differences are exactly zero and")
            log(f"    are discarded by the signed-rank test, so it effectively")
            log(f"    runs on {n - n_zero} pair(s). Treat with suspicion.")
        test_rows.append({'tissue': t, 'outcome': oname, 'n_patients': n,
                          'n_zero_diffs': n_zero, 'hodges_lehmann': hl,
                          'ci_lo': lo, 'ci_hi': hi, 'raw_p': p,
                          'n_positive': n_pos})

td = pd.DataFrame(test_rows)
ok = td['raw_p'].notna()
td['bh_q'] = np.nan
if ok.sum():
    td.loc[ok, 'bh_q'] = multipletests(td.loc[ok, 'raw_p'], method='fdr_bh')[1]
log("\n  BH-corrected across the patient-level family:")
for _, r in td.iterrows():
    if pd.isna(r['hodges_lehmann']):
        log(f"    {r['tissue']:<7s} {r['outcome']:<42s} n={r['n_patients']:<3.0f} "
            f"not evaluable")
        continue
    log(f"    {r['tissue']:<7s} {r['outcome']:<42s} n={r['n_patients']:<3.0f} "
        f"HL={r['hodges_lehmann']:+.4f}  q={fmt_p(r['bh_q'])}")
pd.DataFrame(pt_rows).to_csv(os.path.join(OUTPUT_DIR, "patient_level_medians.tsv"),
                             sep='\t', index=False)
td.to_csv(os.path.join(OUTPUT_DIR, "patient_level_tests.tsv"), sep='\t', index=False)

# =============================================================================
# STEP 3: SUPPORT -- cell-level mixed models
# =============================================================================
banner("STEP 3: SUPPORT -- cell-level models, patient as a random effect")
log("  Read the coefficient and its CI first. With ~24,000 cells the p-values")
log("  will be tiny for effects of no consequence.")
log("  Coefficient = effect of A3A dominance (reference is A3B-dominant),")
log("  adjusted for log total UMI, NOT adjusted for cell-cycle phase.")

model_rows = []

def run_lmm(data, formula, label, tissue, try_slope=True):
    try:
        m = smf.mixedlm(formula, data, groups=data['_patient']).fit(reml=True)
        beta = m.params.get('_dom01', np.nan)
        se = m.bse.get('_dom01', np.nan)
        p = m.pvalues.get('_dom01', np.nan)
        lo, hi = beta - 1.96 * se, beta + 1.96 * se
        slope_note = 'intercept only'
        if try_slope:
            try:
                ms = smf.mixedlm(formula, data, groups=data['_patient'],
                                 re_formula="~_dom01").fit(reml=True)
                # A converged slope here often still carries singular-covariance
                # or boundary warnings, meaning the between-patient slope
                # variance is near zero. Convergence is "not evidence against"
                # consistency, never evidence for it.
                slope_note = ('random slope converged (check for boundary warnings)'
                              if ms.converged else 'random slope did not converge')
            except Exception as e:
                slope_note = f'random slope failed ({type(e).__name__})'
        log(f"\n  {tissue} / {label}   n = {len(data)}")
        log(f"    beta(A3A dominance) = {beta:+.4f}  95% CI {lo:+.4f} to {hi:+.4f}"
            f"   raw p = {fmt_p(p)}")
        log(f"    {slope_note}")
        model_rows.append({'tissue': tissue, 'model': label, 'n': len(data),
                           'beta': beta, 'ci_lo': lo, 'ci_hi': hi, 'raw_p': p,
                           'note': slope_note})
    except Exception as e:
        log(f"\n  {tissue} / {label}: FAILED ({type(e).__name__}: {e})")
        model_rows.append({'tissue': tissue, 'model': label, 'n': len(data),
                           'beta': np.nan, 'ci_lo': np.nan, 'ci_hi': np.nan,
                           'raw_p': np.nan, 'note': f'failed: {e}'})

def run_gee(data, label, tissue):
    """Logistic GEE, exchangeable, clustered by patient. Population-averaged
    odds ratio with a cluster-robust SE."""
    try:
        m = smf.gee("_sbs2_pos ~ _dom01 + _log_umi", groups="_patient",
                    data=data, family=sm.families.Binomial(),
                    cov_struct=sm.cov_struct.Exchangeable()).fit()
        beta = m.params.get('_dom01', np.nan)
        se = m.bse.get('_dom01', np.nan)
        p = m.pvalues.get('_dom01', np.nan)
        orv = float(np.exp(beta))
        lo, hi = float(np.exp(beta - 1.96 * se)), float(np.exp(beta + 1.96 * se))
        log(f"\n  {tissue} / {label}   n = {len(data)}")
        log(f"    OR(A3A dominance) = {orv:.3f}  95% CI {lo:.3f} to {hi:.3f}"
            f"   raw p = {fmt_p(p)}")
        log(f"    population-averaged, cluster-robust by patient")
        model_rows.append({'tissue': tissue, 'model': label, 'n': len(data),
                           'beta': orv, 'ci_lo': lo, 'ci_hi': hi, 'raw_p': p,
                           'note': 'GEE odds ratio, not a beta'})
    except Exception as e:
        log(f"\n  {tissue} / {label}: FAILED ({type(e).__name__}: {e})")
        model_rows.append({'tissue': tissue, 'model': label, 'n': len(data),
                           'beta': np.nan, 'ci_lo': np.nan, 'ci_hi': np.nan,
                           'raw_p': np.nan, 'note': f'failed: {e}'})

for t in ('tumour', 'normal'):
    sub = expr[expr['_tissue'] == t].copy()
    if sub['_patient'].nunique() < 3 or len(sub) < 50:
        log(f"\n  {t}: only {len(sub)} cells across {sub['_patient'].nunique()} "
            f"patients; models not run.")
        continue
    run_lmm(sub, "_cnv ~ _dom01 + _log_umi", "(a) CNA, LMM", t)
    run_lmm(sub, "_sbs2_rankz ~ _dom01 + _log_umi", "(b) SBS2 rank-z, LMM", t)
    run_gee(sub, "(c) SBS2 > 0, GEE logistic", t)
    carriers = sub[sub['_sbs2_pos'] == 1].copy()
    if len(carriers) >= 50 and carriers['_patient'].nunique() >= 3:
        run_lmm(carriers, "_sbs2 ~ _dom01 + _log_umi",
                "(d) SBS2 | >0, LMM", t, try_slope=False)
    else:
        log(f"\n  {t} / (d) SBS2 | >0: only {len(carriers)} carriers; not run.")

md = pd.DataFrame(model_rows)
ok = md['raw_p'].notna()
md['bh_q'] = np.nan
if ok.sum():
    md.loc[ok, 'bh_q'] = multipletests(md.loc[ok, 'raw_p'], method='fdr_bh')[1]
md.to_csv(os.path.join(OUTPUT_DIR, "cell_level_models.tsv"), sep='\t', index=False)

# --- sensitivity: continuous A3A fraction instead of the binary cut ---
banner("SENSITIVITY: continuous A3A fraction instead of binary dominance", "-")
log("  Reported so the binary cut is a visible choice rather than a hidden one.")
for t in ('tumour', 'normal'):
    sub = expr[expr['_tissue'] == t].copy()
    if len(sub) < 50:
        continue
    for f, nm in (("_cnv ~ _a3a_frac + _log_umi", "CNA"),
                  ("_sbs2_rankz ~ _a3a_frac + _log_umi", "SBS2 rank-z")):
        try:
            m = smf.mixedlm(f, sub, groups=sub['_patient']).fit(reml=True)
            b, se = m.params.get('_a3a_frac', np.nan), m.bse.get('_a3a_frac', np.nan)
            log(f"  {t:<7s} {nm:<12s} beta = {b:+.4f}  95% CI "
                f"{b-1.96*se:+.4f} to {b+1.96*se:+.4f}   p = "
                f"{fmt_p(m.pvalues.get('_a3a_frac', np.nan))}")
        except Exception:
            log(f"  {t:<7s} {nm:<12s} FAILED")

# =============================================================================
# STEP 4: verdict
# =============================================================================
banner("STEP 4: VERDICT on the Results 4.1 sentence")

def verdict_for(outcome):
    r = td[(td['tissue'] == 'tumour') & (td['outcome'] == outcome)]
    if r.empty or pd.isna(r.iloc[0]['bh_q']):
        return None, None, None
    r = r.iloc[0]
    return r['hodges_lehmann'], r['bh_q'], (r['ci_lo'], r['ci_hi'])

hl_s, q_s, ci_s = verdict_for('SBS2 prevalence')
hl_sm, q_sm, ci_sm = verdict_for('SBS2 mean weight')
hl_c, q_c, ci_c = verdict_for('CNA mean score')
hl_x, q_x, ci_x = verdict_for(COMPOSITE_NAME)

log('  Sentence under test: "A3A dominance is associated with SBS2, whereas')
log('  A3B dominance is associated with CNA."')
log("")
if q_c is not None:
    ok_c = (q_c < 0.05) and (hl_c < 0)
    log(f"  A3B-to-CNA half:  HL = {hl_c:+.4f} (95% CI {ci_c[0]:+.4f} to "
        f"{ci_c[1]:+.4f}), q = {fmt_p(q_c)}  -> "
        f"{'SUPPORTED' if ok_c else 'NOT SUPPORTED'}")
    log(f"    (expected NEGATIVE: A3A-dominant cells should carry LESS CNA)")
else:
    log("  A3B-to-CNA half: patient-level test not evaluable; read model (a).")

if q_s is not None:
    ok_s = (q_s < 0.05) and (hl_s > 0)
    log(f"  A3A-to-SBS2 half: HL = {hl_s:+.4f} (95% CI {ci_s[0]:+.4f} to "
        f"{ci_s[1]:+.4f}), q = {fmt_p(q_s)}  -> "
        f"{'SUPPORTED' if ok_s else 'NOT SUPPORTED'}")
    log(f"    (expected POSITIVE: A3A-dominant cells should carry MORE SBS2)")
    if q_sm is not None:
        log(f"  A3A-to-SBS2, magnitude: HL = {hl_sm:+.4f} "
            f"(95% CI {ci_sm[0]:+.4f} to {ci_sm[1]:+.4f}), q = {fmt_p(q_sm)}")
        log(f"    Prevalence is the primary reading. SBS2 is zero in most cells,")
        log(f"    so prevalence is where a dominance effect should appear first;")
        log(f"    the mean is reported for magnitude.")

    # --- can Figure 3 carry this? ---
    if q_x is not None:
        _rx = td[(td['tissue'] == 'tumour') & (td['outcome'] == COMPOSITE_NAME)].iloc[0]
        log("")
        log(f"  FIGURE 3 AXIS: HL = {hl_x:+.4f} (95% CI {ci_x[0]:+.4f} to "
            f"{ci_x[1]:+.4f}), q = {fmt_p(q_x)}; A3A-dominant sits further toward")
        log(f"    the SBS2 pole in {_rx['n_positive']:.0f}/{_rx['n_patients']:.0f} patients")
        if (q_x < 0.05) and (hl_x > 0):
            log("    The composite axis separates within patients too, so")
            log("    per-patient median ticks could be overlaid on the two tumour")
            log("    rows of Figure 3 without adding a panel.")
        else:
            log("    The composite axis does NOT separate within patients, even")
            log("    though SBS2 and CNA do separately. Do NOT overlay per-patient")
            log("    marks on Figure 3; the figure cannot carry this claim.")
            log("    Report the patient-level result in text only.")

    log("")
    if ok_s:
        log("  ACTION: the sentence stands, and it can be strengthened. Replace")
        log("  the rho = +0.05 framing with the patient-level test, which is a")
        log("  within-tumour result rather than a pooled correlation. Add the")
        log("  Methods sentence covering the test and the prevalence choice.")
        log("  Whether a Figure 3 overlay goes in now or waits for revision is a")
        log("  separate call; the TEXT should change either way.")
    else:
        log("  ACTION: the sentence outruns the evidence and should be softened")
        log("  BEFORE submission, not after a reviewer finds it. Something like:")
        log("    'A3B dominance was associated with CNA burden, whereas the")
        log("     association between A3A dominance and SBS2 was weak.'")
        log("  This is the outcome that makes running this now worthwhile: it is")
        log("  cheaper to concede the limit than to defend it in review.")

log("")
log("  NORMAL-ADJACENT is the negative control. Dominance should NOT predict")
log("  burden there. If it does, the association is not tumour-specific and the")
log("  whole framing needs rethinking.")
log("  Read from the CELL-LEVEL models: normal-adjacent has too few patients")
log("  with both arms represented for a patient-level paired test to run.")
_nm = md[md['tissue'] == 'normal']
if _nm.empty:
    log("    no normal-adjacent models were fitted")
else:
    for _, r in _nm.iterrows():
        # GEE rows carry an ODDS RATIO, not a signed coefficient. Printing a
        # ratio of 0.675 as '+0.6745' invites it to be read as a beta.
        is_or = 'GEE' in str(r['model'])
        est = (f"OR={r['beta']:.3f} [{r['ci_lo']:.3f}, {r['ci_hi']:.3f}]" if is_or
               else f"est={r['beta']:+.4f} [{r['ci_lo']:+.4f}, {r['ci_hi']:+.4f}]")
        log(f"    normal / {r['model']:<28s} {est}  "
            f"raw p={fmt_p(r['raw_p'])}  q={fmt_p(r['bh_q'])}")

    # Judge on RAW p, not q. A negative control is the one place where being
    # generous about what counts as a signal is the CONSERVATIVE choice: BH
    # correction here makes it EASIER to declare the control clean, which is
    # backwards. raw p < 0.05 gets flagged regardless of q.
    _flag = _nm['raw_p'].notna() & (_nm['raw_p'] < 0.05)
    _hard = _nm['bh_q'].notna() & (_nm['bh_q'] < 0.05)
    if _hard.any():
        log("    >>> FAILS. Dominance predicts burden in NORMAL-ADJACENT tissue")
        log("        after correction. The association is not tumour-specific.")
    elif _flag.any():
        log("    >>> MARGINAL, do not claim a clean negative control.")
        for _, r in _nm[_flag].iterrows():
            log(f"        {r['model']}: raw p = {fmt_p(r['raw_p'])}, "
                f"q = {fmt_p(r['bh_q'])}")
        log("        Check DIRECTION against tumour. An INVERTED sign means this")
        log("        is not a leaked copy of the tumour effect and the framing")
        log("        survives; say 'no evidence of association in normal-adjacent")
        log("        tissue' and note the small number of A3-expressing normal")
        log("        cells. A MATCHING sign means stop.")
        for _, rn in _nm[_flag].iterrows():
            rt = md[(md['tissue'] == 'tumour') & (md['model'] == rn['model'])]
            if not rt.empty and pd.notna(rt.iloc[0]['beta']):
                ref = 1.0 if 'GEE' in str(rn['model']) else 0.0
                same = ((rn['beta'] - ref) * (rt.iloc[0]['beta'] - ref)) > 0
                log(f"        {rn['model']}: tumour {rt.iloc[0]['beta']:+.4f} vs "
                    f"normal {rn['beta']:+.4f} -> direction "
                    f"{'MATCHES tumour (bad)' if same else 'INVERTED (reassuring)'}")
    else:
        log("    >>> No association detected in normal-adjacent tissue at raw")
        log("        p < 0.05. Note the small number of A3-expressing normal")
        log("        cells; this is weak reassurance, not a clean control.")

for _, r in td[td['tissue'] == 'normal'].iterrows():
    log(f"    (patient-level normal / {r['outcome']}: n={r['n_patients']:.0f}, "
        f"not evaluable)")

banner("COMPLETE")
rp = os.path.join(OUTPUT_DIR, "dominance_statistics_report.txt")
with open(rp, 'w') as f:
    f.write("\n".join(report))
log(f"  [SAVE] {rp}")
log(f"  [SAVE] patient_level_medians.tsv / patient_level_tests.tsv / cell_level_models.tsv")
log(f"  Output: {OUTPUT_DIR}")
