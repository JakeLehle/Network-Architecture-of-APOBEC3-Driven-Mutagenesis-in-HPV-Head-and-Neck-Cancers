#!/usr/bin/env python3
"""
Diagnostic_Patient_CellCycle_by_Source.py  (v4)
===============================================
READ-ONLY. Cell-cycle distribution of basal cells, resolved by group, by
patient, and by tissue source (tumor vs normal-adjacent), plus the A3A/A3B
cell-cycle association and the contribution folds under both denominators.

v4 CHANGES (reporting only; no test, number or output file changed)
--------------------------------------------------------------------
  - Cramer's V no longer carries a verbal grade above 'NEGLIGIBLE'. Cohen's
    bands called the A3 dominance split 'small' despite a 25.7 percentage
    point shift in G1 fraction, which is an adjective a reviewer would quote
    back at us. Every chi-square now LEADS with the per-phase percentage-point
    difference; V follows as a scale-free cross-check, never as a verdict.
  - The hero dominance block prints the exact legend wording to use, so the
    figure caption cannot drift from the numbers.
  - Percentage-point differences are written into cellcycle_a3_tests.tsv, so a
    figure script reading that file gets the effect size, not just a p-value.

v2 CHANGES
----------
  - Denominator comes from patient_config.CONTRIBUTION_DENOMINATOR via
    contribution.py. Both folds are always computed and written.
  - patient_source_matched_folds.tsv now also carries the raw counts
    (n_basal, n_tumor, n_sbs2_high, n_cnv_high) so it is a complete standalone
    source for the two contribution figure panels. No figure script needs to
    reopen the AnnData.
  - THREE NEW TSVs so the supplemental figure generators read from disk instead
    of rescoring Tirosh and re-running ccAFv2:
        cellcycle_phase_mix.tsv     Tirosh G1/S/G2M by every stratum
        cellcycle_ccafv2_mix.tsv    ccAFv2 states by every stratum, incl. Unknown
        cellcycle_a3_tests.tsv      every enrichment / distribution test statistic
  - Contributor sets are derived from the fold columns at runtime; the
    patient_config list is an expectation only.
  - The ccAFv2 block is live (the old docstring still called it a stub).

v3 CHANGES
----------
  - Every chi-square is now reported with Cramer's V and the per-phase
    percentage-point differences. At tens of thousands of cells a chi-square
    p-value is close to meaningless on its own, and the contributor-versus-other
    comparison is the case in point: it clears p < 0.05 while the phase mix is
    practically identical. That comparison is KEPT rather than dropped, because
    a null that says patient-level cell-cycle composition does not explain
    contribution is worth seeing, but it is now labelled as a null.
  - Cramer's V and the effect-size labels are written into cellcycle_a3_tests.tsv
    so a figure legend can quote them without recomputation.

Motivation
----------
1. Baseline: cell-cycle phase mix within NORMAL / SBS2-HIGH / CNV-HIGH.
2. A3 test (the clean, non-selection-confounded one, over ALL basal cells):
   replicate the published A3B -> G2/M result with the same Tirosh scoring, and
   test whether A3A -> G1 in our data (the extension that paper left on the table).
3. Per-patient basal composition split by source, and whether the high
   contributors are skewed toward a cell-cycle phase.
4. Contribution folds under both denominators, so the choice is auditable.

Scoring
-------
Primary: scanpy sc.tl.score_genes_cell_cycle with the Tirosh/Regev S and G2/M
gene sets (G1/S/G2M), the same method as the A3B G2/M reference, so the A3
tests are comparable with the literature.
Confirmatory: ccAFv2 (seven states including quiescent qG0). Two environment
fixes are applied at import: the shipped model is Keras-2 so it is loaded with
compile=False, and _prep_predict_data is reimplemented with columns=list(missing)
because pandas 2.x rejects a set.

Outputs (to data/FIG_5/00_diagnostics/)
---------------------------------------
  patient_cellcycle_by_source.tsv        per-patient phase mix by source
  patient_source_matched_folds.tsv       per-patient counts + both folds
  cellcycle_phase_mix.tsv                NEW: Tirosh mix, long format
  cellcycle_ccafv2_mix.tsv               NEW: ccAFv2 mix, long format
  cellcycle_a3_tests.tsv                 NEW: test statistics
  patient_cellcycle_by_source.pdf/.png
  patient_cellcycle_ccAFv2.pdf/.png
  patient_cellcycle_by_source_<ts>.txt

Run from the directory holding patient_config.py and contribution.py:
    conda run -n NETWORK python Diagnostic_Patient_CellCycle_by_Source.py

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from scipy.stats import chi2_contingency, fisher_exact

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ccAFv2 (comprehensive, G0-aware) confirmation tool. The shipped model is Keras-2;
# Keras 3 rejects its saved loss reduction='auto' when rebuilding the TRAINING
# config at import. We only do inference, so force compile=False on model load.
CCAFV2_OK = False
try:
    import keras, keras.models
    _orig_load = keras.models.load_model
    keras.models.load_model = lambda *a, **k: _orig_load(*a, **{**k, 'compile': False})
    try:
        import keras.saving as _ks
        if hasattr(_ks, 'load_model'):
            _ks.load_model = keras.models.load_model
    except Exception:
        pass
    import ccAFv2
    import ccAFv2.ccAFv2 as _ccmod

    # ccAFv2's _prep_predict_data pads missing genes via pd.DataFrame(columns=<set>),
    # which pandas 2.x rejects. Faithful reimplementation, only change: list(missing).
    def _prep_fixed(data, genes):
        print('    Preparing data for classification...')
        data.var_names_make_unique()
        _ccmod.sc.pp.filter_genes(data, min_cells=1)
        in_both = list(set(genes).intersection(data.var_names))
        if len(in_both) > 0:
            print('    Marker genes present in this dataset: ' + str(len(in_both)))
            print('    Missing marker genes in this dataset: ' + str(len(set(genes)) - len(in_both)))
            data2 = data[:, in_both]
            if _ccmod.isspmatrix(data.X):
                data2 = _ccmod.pd.DataFrame(data2.X.todense(), index=data2.obs_names, columns=data2.var_names)
            else:
                data2 = _ccmod.pd.DataFrame(data2.X, index=data2.obs_names, columns=data2.var_names)
            data3 = _ccmod.pd.DataFrame(_ccmod._scale(data2), index=data2.index, columns=data2.columns)
            missing = set(genes).difference(data3.columns)
            if len(missing) > 0:
                data4 = _ccmod.pd.concat(
                    [data3, _ccmod.pd.DataFrame(data3.values.min(), index=data3.index,
                                                columns=list(missing))], axis=1)
                return data4[list(genes)]
            return data3
        raise RuntimeError('No overlap between input genes and classifier genes (check species/gene_id).')

    _ccmod._prep_predict_data = _prep_fixed
    CCAFV2_OK = True
except Exception as _e:
    print(f"[warn] ccAFv2 unavailable ({type(_e).__name__}); running Tirosh-only.",
          flush=True)

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
OUT_DIR   = ensure_dir(DIR_00_DIAG)
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")

SOURCE_COL    = 'source_name'
NORMAL_SOURCE = 'normal tissue adjucent to head and neck squamous cell carcinoma'
PHASE_ORDER   = ['G1', 'S', 'G2M']
POP_ORDER     = ['SBS2_HIGH', 'CNV_HIGH', 'NORMAL']
POP_LABELS    = {'SBS2_HIGH': 'SBS2-HIGH', 'CNV_HIGH': 'CNV-HIGH', 'NORMAL': 'Normal'}

# Phase palette. G1 deliberately does NOT reuse the NORMAL population blue
# (#4682b4); a phase band and a population band must not read as the same thing.
PHASE_COLORS = {'G1': '#6b8f71', 'S': '#f4a259', 'G2M': '#ed6a5a'}

CCAFV2_ORDER  = ['qG0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1']
CCAFV2_COLORS = {'qG0': '#3b5c7a', 'G1': '#6b8f71', 'Late G1': '#9bb4c7',
                 'S': '#f4a259', 'S/G2': '#e8823b', 'G2/M': '#ed6a5a',
                 'M/Early G1': '#c04a3b'}
DPI = 300

# Tirosh / Regev S and G2/M gene sets (scanpy skips any absent from var_names)
S_GENES = ['MCM5','PCNA','TYMS','FEN1','MCM2','MCM4','RRM1','UNG','GINS2','MCM6',
           'CDCA7','DTL','PRIM1','UHRF1','MLF1IP','HELLS','RFC2','RPA2','NASP',
           'RAD51AP1','GMNN','WDR76','SLBP','CCNE2','UBR7','POLD3','MSH2','ATAD2',
           'RAD51','RRM2','CDC45','CDC6','EXO1','TIPIN','DSCC1','BLM','CASP8AP2',
           'USP1','CLSPN','POLA1','CHAF1B','BRIP1','E2F8']
G2M_GENES = ['HMGB2','CDK1','NUSAP1','UBE2C','BIRC5','TPX2','TOP2A','NDC80','CKS2',
             'NUF2','CKS1B','MKI67','TMPO','CENPF','TACC3','FAM64A','SMC4','CCNB2',
             'CKAP2L','CKAP2','AURKB','BUB1','KIF11','ANP32E','TUBB4B','GTSE1',
             'KIF20B','HJURP','CDCA3','HN1','CDC20','TTK','CDC25C','KIF2C','RANGAP1',
             'NCAPD2','DLGAP5','CDCA2','CDCA8','ECT2','KIF23','HMMR','AURKA','PSRC1',
             'ANLN','LBR','CKAP5','CENPE','CTCF','NEK2','G2E3','GAS2L3','CBX5','CENPA']

report_lines = []


def rlog(msg=""):
    log(msg)
    report_lines.append(str(msg))


# Accumulators for the three new TSVs
phase_mix_rows  = []
ccafv2_mix_rows = []
a3_test_rows    = []


def record_phase_mix(stratum_type, stratum, mask):
    """Append the Tirosh phase mix for a boolean mask to the long-format table."""
    counts, pct, n = phase_pcts(mask)
    for ph in PHASE_ORDER:
        phase_mix_rows.append({
            'scorer': 'tirosh', 'stratum_type': stratum_type, 'stratum': stratum,
            'n_total': n, 'state': ph, 'n': counts[ph], 'pct': pct[ph]})
    return counts, pct, n


def record_ccafv2_mix(stratum_type, stratum, mask):
    """Append the ccAFv2 state mix for a boolean mask. Includes every observed
    label, Unknown included, because the below-threshold fraction has to be
    quotable in the figure legend."""
    if 'ccAFv2' not in basal.obs.columns:
        return None, 0
    s = basal.obs.loc[mask, 'ccAFv2']
    n = len(s)
    observed = list(CCAFV2_ORDER) + [x for x in pd.unique(s) if x not in CCAFV2_ORDER]
    d = {}
    for st in observed:
        cnt = int((s == st).sum())
        pct = 100.0 * cnt / n if n else 0.0
        d[st] = pct
        ccafv2_mix_rows.append({
            'scorer': 'ccAFv2', 'stratum_type': stratum_type, 'stratum': stratum,
            'n_total': n, 'state': st, 'n': cnt, 'pct': pct})
    return d, n


def record_test(test, comparison, stat_type, stat, p, n1=np.nan, n2=np.nan,
                note=''):
    a3_test_rows.append({'test': test, 'comparison': comparison,
                         'stat_type': stat_type, 'stat': stat, 'p': p,
                         'n1': n1, 'n2': n2, 'note': note})


def get_expr(ad, gene):
    if gene not in ad.var_names:
        return np.zeros(ad.n_obs)
    x = ad[:, gene].X
    return (np.asarray(x.todense()).flatten() if scipy.sparse.issparse(x)
            else np.asarray(x).flatten())


def phase_pcts(mask):
    """Return (counts, pct, n) for cells under a boolean mask."""
    ph = basal.obs.loc[mask, 'phase']
    n = len(ph)
    counts = {p: int((ph == p).sum()) for p in PHASE_ORDER}
    pct = {p: (100.0 * counts[p] / n if n else 0.0) for p in PHASE_ORDER}
    return counts, pct, n


def cramers_v(chi2_stat, table):
    """
    Effect size for a contingency table. At tens of thousands of cells a
    chi-square p-value is close to meaningless on its own, so every chi-square
    in this script is reported with V alongside it.

    Conventional reading for these table shapes: V < 0.10 negligible,
    0.10-0.30 small, 0.30-0.50 moderate, > 0.50 large.
    """
    n = int(np.asarray(table).sum())
    k = min(np.asarray(table).shape) - 1
    if n <= 0 or k <= 0:
        return float('nan')
    return float(np.sqrt(chi2_stat / (n * k)))


def v_label(v):
    """
    Deliberately ONE-SIDED. The only verdict this function issues is
    'NEGLIGIBLE', which is the case worth flagging because a large n can hand
    you a tiny p-value on a difference of a percentage point or two.

    Above that threshold it returns an empty string on purpose. Cohen's
    conventional bands (0.10 small, 0.30 moderate, 0.50 large) are calibrated
    for questions unlike this one, and applying them literally here labels the
    A3 dominance split 'small' despite a 25.7 percentage point shift in G1
    fraction. A reviewer would rightly quote that adjective back at us. The
    percentage-point difference is the honest effect measure for a phase mix,
    so that is what gets reported and what belongs in a figure legend; V is
    kept alongside it as a scale-free cross-check, not as a verdict.
    """
    if not np.isfinite(v):
        return 'undefined'
    if v < 0.10:
        return 'NEGLIGIBLE'
    return ''


def v_note(v):
    """Formatted V with the negligible flag appended only when it applies."""
    if not np.isfinite(v):
        return "Cramer's V = undefined"
    lab = v_label(v)
    return f"Cramer's V = {v:.3f}" + (f" ({lab})" if lab else "")


def pp_spread(pct_a, pct_b):
    """
    Per-phase difference in percentage points, and the largest absolute one.
    This is the effect measure that should be quoted in prose and legends.
    """
    deltas = {ph: pct_a[ph] - pct_b[ph] for ph in PHASE_ORDER}
    return deltas, max(abs(d) for d in deltas.values())


# =============================================================================
# STEP 0: LOAD, TAG, SCORE
# =============================================================================
banner("STEP 0: LOAD + SCORE CELL CYCLE")
adata = load_adata()
sbs2_high, cnv_high, normal = load_three_groups()
basal = adata[adata.obs[CELLTYPE_COL] == 'basal cell'].copy()

basal.obs['group'] = 'other'
basal.obs.loc[basal.obs_names.isin(sbs2_high), 'group'] = 'SBS2_HIGH'
basal.obs.loc[basal.obs_names.isin(cnv_high), 'group'] = 'CNV_HIGH'
basal.obs.loc[basal.obs_names.isin(normal), 'group'] = 'NORMAL'

if SOURCE_COL in basal.obs.columns:
    basal.obs['tissue_grp'] = np.where(
        basal.obs[SOURCE_COL].astype(str) == NORMAL_SOURCE, 'normal', 'tumor')
else:
    raise SystemExit(f"ERROR: '{SOURCE_COL}' not in adata.obs; cannot split source.")

n_basal_total = int(basal.n_obs)
n_tumor       = int((basal.obs['tissue_grp'] == 'tumor').sum())
n_normal_t    = n_basal_total - n_tumor
rlog(f"  Basal: {n_basal_total:,}  (tumor {n_tumor:,}, "
     f"normal-adjacent {n_normal_t:,})")
rlog("")
announce(rlog, n_basal_total, n_tumor)
rlog("")

basal.obs['A3A'] = get_expr(basal, 'APOBEC3A')
basal.obs['A3B'] = get_expr(basal, 'APOBEC3B')

s_in  = [g for g in S_GENES if g in basal.var_names]
g2_in = [g for g in G2M_GENES if g in basal.var_names]
rlog(f"  Tirosh gene sets present: S {len(s_in)}/{len(S_GENES)}, "
     f"G2M {len(g2_in)}/{len(G2M_GENES)}")
sc.tl.score_genes_cell_cycle(basal, s_genes=s_in, g2m_genes=g2_in)
basal.obs['phase'] = basal.obs['phase'].astype(str)
overall = basal.obs['phase'].value_counts()
rlog("  Overall basal phase mix: "
     + ", ".join(f"{p} {100*overall.get(p,0)/n_basal_total:.1f}%"
                 for p in PHASE_ORDER))

if CCAFV2_OK:
    try:
        # pass a copy: _prep_predict_data mutates its input (var_names_make_unique,
        # filter_genes). Labels align to basal.obs (cells are not dropped).
        res = ccAFv2.predict_labels(basal.copy(), species='human', gene_id='symbol')
        labels = res[0] if isinstance(res, tuple) else res   # returns (labels, probs)
        basal.obs['ccAFv2'] = np.asarray(labels).astype(str)
        ccc = basal.obs['ccAFv2'].value_counts()
        rlog("  ccAFv2 states: " + ", ".join(
            f"{s} {100*ccc.get(s,0)/n_basal_total:.1f}%"
            for s in CCAFV2_ORDER if s in ccc.index))
        extra = [s for s in ccc.index if s not in CCAFV2_ORDER]
        if extra:
            rlog("  (other labels incl. below-threshold 'Unknown': "
                 + ", ".join(f"{s} {100*ccc.get(s,0)/n_basal_total:.1f}%"
                             for s in extra) + ")")
            rlog("  These are written to cellcycle_ccafv2_mix.tsv and must be "
                 "stated in any ccAFv2 figure legend.")
    except Exception as e:
        rlog(f"  [warn] ccAFv2 failed: {type(e).__name__}: {str(e)[:200]}")
        CCAFV2_OK = False
else:
    rlog("  ccAFv2 not available; Tirosh-only.")

# =============================================================================
# STEP A: PER-GROUP BASELINE
# =============================================================================
banner("STEP A: CELL-CYCLE MIX WITHIN NORMAL / SBS2-HIGH / CNV-HIGH")
rlog("  (CNV-HIGH was selected on CNV+stemness ~ proliferation, so its G2/M")
rlog("   enrichment is partly circular; the clean test is STEP B.)\n")
rlog(f"  {'Group':<12s} {'n':>6s} {'G1':>7s} {'S':>7s} {'G2M':>7s}")
for p in POP_ORDER:
    _, pct, n = record_phase_mix('group', p, basal.obs['group'] == p)
    rlog(f"  {POP_LABELS[p]:<12s} {n:>6d} {pct['G1']:>6.1f}% {pct['S']:>6.1f}% "
         f"{pct['G2M']:>6.1f}%")
# whole basal compartment, for context in the figure
record_phase_mix('all', 'all_basal', pd.Series(True, index=basal.obs.index))

ct = np.array([[int(((basal.obs['group'] == p) & (basal.obs['phase'] == ph)).sum())
                for ph in PHASE_ORDER] for p in POP_ORDER])
chi2, pv, dof, _ = chi2_contingency(ct)
v = cramers_v(chi2, ct)
_pct = {p: phase_pcts(basal.obs['group'] == p)[1] for p in POP_ORDER}
_d_sc, _max_sc = pp_spread(_pct['SBS2_HIGH'], _pct['CNV_HIGH'])
rlog(f"\n  group x phase chi-square: chi2={chi2:.1f}, dof={dof}, p={pv:.2e}")
rlog(f"  EFFECT (SBS2-HIGH minus CNV-HIGH, percentage points): "
     + ", ".join(f"{ph} {_d_sc[ph]:+.1f}" for ph in PHASE_ORDER))
rlog(f"  largest single-phase difference: {_max_sc:.1f} percentage points")
rlog(f"  {v_note(v)} over n={int(ct.sum()):,} cells")
rlog(f"  Quote the percentage points in prose and legends, not V and not the p.")
record_test('group_x_phase', 'SBS2-HIGH vs CNV-HIGH vs NORMAL', 'chi2', chi2, pv,
            n1=int(ct.sum()),
            note=f"dof={dof}; cramers_v={v:.4f}; max phase delta "
                 f"{_max_sc:.1f} pp (SBS2 vs CNV); "
                 f"CNV-HIGH selection is partly circular here")

# =============================================================================
# STEP B: A3 CELL-CYCLE ASSOCIATION (all basal; clean)
# =============================================================================
banner("STEP B: A3A / A3B CELL-CYCLE ASSOCIATION (all basal)")


def a3_test(flag_name, pos_mask, target_phase):
    _, pct_pos, n_pos = record_phase_mix('a3_flag', f'{flag_name}_pos', pos_mask)
    _, pct_neg, n_neg = record_phase_mix('a3_flag', f'{flag_name}_neg', ~pos_mask)
    a = int((pos_mask & (basal.obs['phase'] == target_phase)).sum())
    b = n_pos - a
    c = int(((~pos_mask) & (basal.obs['phase'] == target_phase)).sum())
    d = n_neg - c
    orr, pf = fisher_exact([[a, b], [c, d]])
    rlog(f"  {flag_name}+ (n={n_pos}) vs {flag_name}- (n={n_neg}), "
         f"target {target_phase}:")
    rlog(f"    {flag_name}+ phase mix: "
         + ", ".join(f"{p} {pct_pos[p]:.1f}%" for p in PHASE_ORDER))
    rlog(f"    {flag_name}- phase mix: "
         + ", ".join(f"{p} {pct_neg[p]:.1f}%" for p in PHASE_ORDER))
    rlog(f"    {target_phase} enrichment: OR={orr:.2f}, Fisher p={pf:.2e}  "
         f"({abs(pct_pos[target_phase] - pct_neg[target_phase]):.1f} percentage "
         f"points)")
    rlog("")
    record_test(f'{flag_name}_vs_{target_phase}',
                f'{flag_name}+ vs {flag_name}-', 'odds_ratio', orr, pf,
                n1=n_pos, n2=n_neg,
                note=f'all basal cells, not selection-confounded; '
                     f'{abs(pct_pos[target_phase] - pct_neg[target_phase]):.1f} pp '
                     f'difference in {target_phase}')


rlog("  Replication target (A3B -> G2/M) and extension (A3A -> G1):\n")
a3_test('A3B', basal.obs['A3B'] > 0, 'G2M')
a3_test('A3A', basal.obs['A3A'] > 0, 'G1')

# dominance version among A3-expressing cells (cleanest enzyme axis)
expr = (basal.obs['A3A'] + basal.obs['A3B']) > 0
frac = basal.obs['A3A'] / (basal.obs['A3A'] + basal.obs['A3B']).replace(0, np.nan)
a3a_dom = expr & (frac > 0.5)
a3b_dom = expr & (frac <= 0.5)
_, pa, na = record_phase_mix('a3_dominance', 'A3A_dominant', a3a_dom)
_, pb, nb = record_phase_mix('a3_dominance', 'A3B_dominant', a3b_dom)
rlog("  Among A3-expressing basal cells (dominance split):")
rlog(f"    A3A-dominant (n={na}): "
     + ", ".join(f"{p} {pa[p]:.1f}%" for p in PHASE_ORDER))
rlog(f"    A3B-dominant (n={nb}): "
     + ", ".join(f"{p} {pb[p]:.1f}%" for p in PHASE_ORDER))
ctd = np.array([[int((a3a_dom & (basal.obs['phase'] == ph)).sum()) for ph in PHASE_ORDER],
                [int((a3b_dom & (basal.obs['phase'] == ph)).sum()) for ph in PHASE_ORDER]])
if ctd.sum() > 0 and (ctd.sum(axis=1) > 0).all():
    c2, p2, dof2, _ = chi2_contingency(ctd)
    v2 = cramers_v(c2, ctd)
    d_dom, max_dom = pp_spread(pa, pb)
    rlog(f"    EFFECT (A3A-dominant minus A3B-dominant, percentage points): "
         + ", ".join(f"{ph} {d_dom[ph]:+.1f}" for ph in PHASE_ORDER))
    rlog(f"    largest single-phase difference: {max_dom:.1f} percentage points")
    rlog(f"    dominance x phase chi-square: chi2={c2:.1f}, p={p2:.2e}; "
         f"{v_note(v2)} over n={na + nb:,} A3-expressing cells")
    rlog(f"    LEGEND WORDING for this panel: quote the percentage-point shift,")
    rlog(f"    e.g. 'A3A-dominant cells are {d_dom['G1']:+.1f} percentage points "
         f"more likely to be")
    rlog(f"    in G1 ({pa['G1']:.1f}% vs {pb['G1']:.1f}%)'. Do NOT quote "
         f"Cramer's V: its conventional")
    rlog(f"    bands are calibrated for other questions and would understate "
         f"a shift this size.")
    record_test('dominance_x_phase', 'A3A-dominant vs A3B-dominant', 'chi2',
                c2, p2, n1=na, n2=nb,
                note=f"dof={dof2}; cramers_v={v2:.4f}; "
                     f"G1 delta {d_dom['G1']:+.1f} pp; max delta {max_dom:.1f} pp; "
                     f"the non-circular hero comparison")

# =============================================================================
# STEP B2: ccAFv2 CONFIRMATION (7-state, G0-aware)
# =============================================================================
da = db = None
if CCAFV2_OK and 'ccAFv2' in basal.obs.columns:
    banner("STEP B2: ccAFv2 CONFIRMATION (comprehensive, quiescence-resolved)")

    rlog("  Per-group ccAFv2 state mix:")
    for p in POP_ORDER:
        d, n = record_ccafv2_mix('group', p, basal.obs['group'] == p)
        rlog(f"    {POP_LABELS[p]:<10s} (n={n}): "
             + ", ".join(f"{st} {d[st]:.0f}%" for st in d if d[st] >= 1))
    record_ccafv2_mix('all', 'all_basal', pd.Series(True, index=basal.obs.index))

    da, na_c = record_ccafv2_mix('a3_dominance', 'A3A_dominant', a3a_dom)
    db, nb_c = record_ccafv2_mix('a3_dominance', 'A3B_dominant', a3b_dom)
    rlog("\n  By A3 dominance (does A3A carry qG0 quiescence Tirosh lumped into G1?):")
    rlog(f"    A3A-dominant (n={na_c}): "
         + ", ".join(f"{st} {da[st]:.0f}%" for st in da if da[st] >= 1))
    rlog(f"    A3B-dominant (n={nb_c}): "
         + ", ".join(f"{st} {db[st]:.0f}%" for st in db if db[st] >= 1))
    qg1 = lambda d: d.get('qG0', 0) + d.get('G1', 0) + d.get('Late G1', 0)
    sg2 = lambda d: d.get('S', 0) + d.get('S/G2', 0) + d.get('G2/M', 0)
    rlog(f"    A3A-dom quiescent+G1 = {qg1(da):.0f}% vs S+G2/M = {sg2(da):.0f}%")
    rlog(f"    A3B-dom quiescent+G1 = {qg1(db):.0f}% vs S+G2/M = {sg2(db):.0f}%")
    rlog(f"    qG0 specifically: A3A-dom {da.get('qG0',0):.0f}% vs "
         f"A3B-dom {db.get('qG0',0):.0f}%")
    if abs(da.get('qG0', 0) - db.get('qG0', 0)) < 2.0:
        rlog("    READ: qG0 is comparable between the two, so the claim is")
        rlog("          A3A-with-G1, NOT A3A-with-quiescence.")

# =============================================================================
# STEP C: PER-PATIENT BASAL BY SOURCE + PHASE
# =============================================================================
banner("STEP C: PER-PATIENT BASAL COMPOSITION BY SOURCE")
patients = sorted(basal.obs[PATIENT_COL].unique())
rows = []
for p in patients:
    pm = basal.obs[PATIENT_COL] == p
    tmask = pm & (basal.obs['tissue_grp'] == 'tumor')
    nmask = pm & (basal.obs['tissue_grp'] == 'normal')
    _, pct_t, nt = phase_pcts(tmask)
    _, pct_n, nn = phase_pcts(nmask)
    rows.append({'patient': p, 'n_tumor': nt, 'n_normal_adj': nn,
                 'tumor_G1': pct_t['G1'], 'tumor_S': pct_t['S'],
                 'tumor_G2M': pct_t['G2M'],
                 'norm_G1': pct_n['G1'], 'norm_S': pct_n['S'],
                 'norm_G2M': pct_n['G2M']})
    # tissue strata for the tumor-vs-normal control panel
    if nt > 0:
        record_phase_mix('tissue', f'{short(p)}_tumor', tmask)
    if nn > 0:
        record_phase_mix('tissue', f'{short(p)}_normal_adj', nmask)

pc = pd.DataFrame(rows)
rlog(f"  {'Patient':<10s} {'nTum':>6s} {'nNorm':>6s} | tumor G1/S/G2M "
     f"| normal-adj G1/S/G2M")
for _, r in pc.iterrows():
    rlog(f"  {short(r['patient']):<10s} {int(r['n_tumor']):>6d} "
         f"{int(r['n_normal_adj']):>6d} | "
         f"{r['tumor_G1']:.0f}/{r['tumor_S']:.0f}/{r['tumor_G2M']:.0f} | "
         f"{r['norm_G1']:.0f}/{r['norm_S']:.0f}/{r['norm_G2M']:.0f}")

# pooled tumor vs normal-adjacent, over patients that have both
both_src = set(pc.loc[(pc['n_tumor'] > 0) & (pc['n_normal_adj'] > 0), 'patient'])
if both_src:
    bm = basal.obs[PATIENT_COL].isin(both_src)
    record_phase_mix('tissue_pooled', 'tumor',
                     bm & (basal.obs['tissue_grp'] == 'tumor'))
    record_phase_mix('tissue_pooled', 'normal_adj',
                     bm & (basal.obs['tissue_grp'] == 'normal'))
    rlog(f"\n  Patients with BOTH sources ({len(both_src)}): "
         f"{sorted(short(p) for p in both_src)}")
    rlog("  Pooled tumor vs normal-adjacent written as stratum_type "
         "'tissue_pooled' (CONTROL, not a result).")

# =============================================================================
# STEP D: CONTRIBUTION FOLDS (BOTH DENOMINATORS) + HC PHASE SKEW
# =============================================================================
banner("STEP D: CONTRIBUTION FOLDS + HIGH-CONTRIBUTOR PHASE SKEW")

n_s_total = int((basal.obs['group'] == 'SBS2_HIGH').sum())
n_c_total = int((basal.obs['group'] == 'CNV_HIGH').sum())

frows = []
for p in patients:
    pm = basal.obs[PATIENT_COL] == p
    n_basal_p = int(pm.sum())
    n_tum_p   = int((pm & (basal.obs['tissue_grp'] == 'tumor')).sum())
    n_sbs2_p  = int((pm & (basal.obs['group'] == 'SBS2_HIGH')).sum())
    n_cnv_p   = int((pm & (basal.obs['group'] == 'CNV_HIGH')).sum())

    fs = both_folds(n_sbs2_p, n_s_total, n_basal_p, n_basal_total,
                    n_tum_p, n_tumor)
    fc = both_folds(n_cnv_p, n_c_total, n_basal_p, n_basal_total,
                    n_tum_p, n_tumor)

    frows.append({
        'patient': p,
        'n_basal': n_basal_p,
        'n_tumor': n_tum_p,
        'n_normal_adj': n_basal_p - n_tum_p,
        'n_sbs2_high': n_sbs2_p,
        'n_cnv_high': n_cnv_p,
        'fold_sbs2_all_basal': fs['all_basal'], 'fold_sbs2_tumor': fs['tumor'],
        'fold_cnv_all_basal':  fc['all_basal'], 'fold_cnv_tumor':  fc['tumor'],
    })

fd = pd.DataFrame(frows)
fd = attach_folds(fd, 'fold_sbs2')
fd = attach_folds(fd, 'fold_cnv')

rlog(f"  {'Patient':<10s} {'SBS2 all -> tumor':>20s}   {'CNV all -> tumor':>20s}")
rlog(f"  {'-'*10} {'-'*20}   {'-'*20}")
for _, r in fd.sort_values('fold_cnv', ascending=False).iterrows():
    rlog(f"  {short(r['patient']):<10s} "
         f"{r['fold_sbs2_all_basal']:>8.2f}x ->{r['fold_sbs2_tumor']:>7.2f}x   "
         f"{r['fold_cnv_all_basal']:>8.2f}x ->{r['fold_cnv_tumor']:>7.2f}x")

rlog("")
sbs2_hc = derive_contributors(fd, 'fold_sbs2', expected=HIGH_CONTRIBUTORS,
                              label='SBS2-HIGH contributors', logger=rlog)
rlog("")
cnv_hc = derive_contributors(fd, 'fold_cnv', expected=CNV_HIGH_CONTRIBUTORS,
                             label='CNV-HIGH contributors', logger=rlog)

# phase skew of HC tumor basal cells vs other tumor basal cells
hc_all = sbs2_hc | cnv_hc
hc_mask = basal.obs[PATIENT_COL].isin(hc_all) & (basal.obs['tissue_grp'] == 'tumor')
other_mask = (~basal.obs[PATIENT_COL].isin(hc_all)) & (basal.obs['tissue_grp'] == 'tumor')
_, pct_hc, nhc = record_phase_mix('contributor', 'high_contributor_tumor', hc_mask)
_, pct_ot, notr = record_phase_mix('contributor', 'other_patient_tumor', other_mask)
rlog("\n  Tumor basal phase mix, HIGH contributors vs other patients:")
rlog(f"    HC patients    (n={nhc}): "
     + ", ".join(f"{p} {pct_hc[p]:.1f}%" for p in PHASE_ORDER))
rlog(f"    other patients (n={notr}): "
     + ", ".join(f"{p} {pct_ot[p]:.1f}%" for p in PHASE_ORDER))

cth = np.array([[int((hc_mask & (basal.obs['phase'] == ph)).sum()) for ph in PHASE_ORDER],
                [int((other_mask & (basal.obs['phase'] == ph)).sum()) for ph in PHASE_ORDER]])
if cth.sum() > 0 and (cth.sum(axis=1) > 0).all():
    c3, p3, dof3, _ = chi2_contingency(cth)
    v3 = cramers_v(c3, cth)
    deltas, max_delta = pp_spread(pct_hc, pct_ot)
    rlog(f"    EFFECT (HC minus other, percentage points): "
         + ", ".join(f"{ph} {deltas[ph]:+.1f}" for ph in PHASE_ORDER))
    rlog(f"    largest single-phase difference: {max_delta:.1f} percentage points")
    rlog(f"    contributor x phase chi-square: chi2={c3:.1f}, p={p3:.2e}; "
         f"{v_note(v3)} over n={nhc + notr:,} tumor basal cells")
    if v3 < 0.10:
        rlog("    READ: this is a NULL. The p-value is small only because n is")
        rlog("          large; the effect size is negligible and the phase mix is")
        rlog("          practically identical between contributors and everyone")
        rlog("          else. Report the percentages, never the p-value alone.")
        rlog("          A useful null: high contribution is NOT explained by the")
        rlog("          patient's overall cell-cycle composition, which is what")
        rlog("          the per-patient cell-cycle correlations also found.")
    record_test('contributor_x_phase', 'high contributors vs other patients',
                'chi2', c3, p3, n1=nhc, n2=notr,
                note=f"dof={dof3}; cramers_v={v3:.4f}; "
                     f"max phase delta {max_delta:.1f} pp; tumor basal only; "
                     f"{'NULL by effect size' if v3 < 0.10 else 'see effect size'}")

# =============================================================================
# STEP E: WRITE TABLES + PLOTS
# =============================================================================
banner("STEP E: WRITE TABLES + PLOTS")

pc.to_csv(os.path.join(OUT_DIR, "patient_cellcycle_by_source.tsv"),
          sep='\t', index=False)
fd.to_csv(os.path.join(OUT_DIR, "patient_source_matched_folds.tsv"),
          sep='\t', index=False)

pd.DataFrame(phase_mix_rows).to_csv(
    os.path.join(OUT_DIR, "cellcycle_phase_mix.tsv"), sep='\t', index=False)
pd.DataFrame(a3_test_rows).to_csv(
    os.path.join(OUT_DIR, "cellcycle_a3_tests.tsv"), sep='\t', index=False)
if ccafv2_mix_rows:
    pd.DataFrame(ccafv2_mix_rows).to_csv(
        os.path.join(OUT_DIR, "cellcycle_ccafv2_mix.tsv"), sep='\t', index=False)

rlog(f"  [SAVE] patient_cellcycle_by_source.tsv   ({len(pc)} rows)")
rlog(f"  [SAVE] patient_source_matched_folds.tsv  ({len(fd)} rows, counts + both folds)")
rlog(f"  [SAVE] cellcycle_phase_mix.tsv           ({len(phase_mix_rows)} rows)")
rlog(f"  [SAVE] cellcycle_a3_tests.tsv            ({len(a3_test_rows)} rows)")
if ccafv2_mix_rows:
    rlog(f"  [SAVE] cellcycle_ccafv2_mix.tsv          ({len(ccafv2_mix_rows)} rows)")
else:
    rlog("  [SKIP] cellcycle_ccafv2_mix.tsv (ccAFv2 unavailable this run); the")
    rlog("         ccAFv2 supplemental panel cannot be built until this succeeds.")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
x = np.arange(len(POP_ORDER)); left = np.zeros(len(POP_ORDER))
for ph in PHASE_ORDER:
    vals = np.array([phase_pcts(basal.obs['group'] == p)[1][ph] for p in POP_ORDER])
    ax1.barh(x, vals, left=left, color=PHASE_COLORS[ph], label=ph, edgecolor='white')
    left += vals
ax1.set_yticks(x); ax1.set_yticklabels([POP_LABELS[p] for p in POP_ORDER])
ax1.set_xlabel('% of cells'); ax1.set_title('Cell-cycle mix by group')
ax1.legend(fontsize=9)

left2 = np.zeros(2)
for ph in PHASE_ORDER:
    vals = np.array([pa[ph], pb[ph]])
    ax2.barh([0, 1], vals, left=left2, color=PHASE_COLORS[ph], label=ph,
             edgecolor='white')
    left2 += vals
ax2.set_yticks([0, 1]); ax2.set_yticklabels(['A3A-dominant', 'A3B-dominant'])
ax2.set_xlabel('% of cells'); ax2.set_title('Cell-cycle mix by A3 dominance')
ax2.legend(fontsize=9)
plt.tight_layout()
for ext in ('pdf', 'png'):
    plt.savefig(os.path.join(OUT_DIR, f"patient_cellcycle_by_source.{ext}"),
                dpi=DPI, bbox_inches='tight')
plt.close()
rlog("  [SAVE] patient_cellcycle_by_source.pdf/.png")

if CCAFV2_OK and 'ccAFv2' in basal.obs.columns and da is not None:
    figc, (axc1, axc2) = plt.subplots(1, 2, figsize=(14, 5))
    xg = np.arange(len(POP_ORDER)); leftg = np.zeros(len(POP_ORDER))
    for st in CCAFV2_ORDER:
        vals = np.array([
            100.0 * ((basal.obs['group'] == p) & (basal.obs['ccAFv2'] == st)).sum()
            / max((basal.obs['group'] == p).sum(), 1) for p in POP_ORDER])
        axc1.barh(xg, vals, left=leftg, color=CCAFV2_COLORS[st], label=st,
                  edgecolor='white')
        leftg += vals
    axc1.set_yticks(xg); axc1.set_yticklabels([POP_LABELS[p] for p in POP_ORDER])
    axc1.set_xlabel('% of cells'); axc1.set_title('ccAFv2 states by group')
    axc1.legend(fontsize=7, ncol=2, loc='lower right')
    leftd = np.zeros(2)
    for st in CCAFV2_ORDER:
        vals = np.array([da.get(st, 0), db.get(st, 0)])
        axc2.barh([0, 1], vals, left=leftd, color=CCAFV2_COLORS[st], label=st,
                  edgecolor='white')
        leftd += vals
    axc2.set_yticks([0, 1]); axc2.set_yticklabels(['A3A-dominant', 'A3B-dominant'])
    axc2.set_xlabel('% of cells'); axc2.set_title('ccAFv2 states by A3 dominance')
    axc2.legend(fontsize=7, ncol=2, loc='lower right')
    plt.tight_layout()
    for ext in ('pdf', 'png'):
        plt.savefig(os.path.join(OUT_DIR, f"patient_cellcycle_ccAFv2.{ext}"),
                    dpi=DPI, bbox_inches='tight')
    plt.close()
    rlog("  [SAVE] patient_cellcycle_ccAFv2.pdf/.png")

report_path = os.path.join(OUT_DIR, f"patient_cellcycle_by_source_{TIMESTAMP}.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
rlog(f"\n  Report: {report_path}")
banner("CELL-CYCLE BY SOURCE DIAGNOSTIC COMPLETE (v4)")
