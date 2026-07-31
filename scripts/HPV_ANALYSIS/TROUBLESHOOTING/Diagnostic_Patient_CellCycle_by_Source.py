#!/usr/bin/env python3
"""
Diagnostic_Patient_CellCycle_by_Source.py
=========================================
READ-ONLY. Cell-cycle distribution of basal cells, resolved by group, by
patient, and by tissue source (tumor vs normal-adjacent), plus the A3A/A3B
cell-cycle association and the source-matched contribution folds.

Motivation
----------
1. Baseline: cell-cycle phase mix within NORMAL / SBS2-HIGH / CNV-HIGH.
2. A3 test (the clean, non-selection-confounded one, over ALL basal cells):
   replicate the published A3B -> G2/M result with the same Tirosh scoring, and
   test whether A3A -> G1 in our data (the extension that paper left on the table).
3. Per-patient basal composition split by source, and whether the high
   contributors are skewed toward a cell-cycle phase.
4. Source-matched contribution folds: the tumor groups can only be seeded by
   tumor basal cells (Step00B), so tumor-group fold should divide by TUMOR basal
   share, not all-basal share. Both are reported so the change is visible.

Scoring
-------
Primary: scanpy sc.tl.score_genes_cell_cycle with the Tirosh/Regev S and G2/M
gene sets (G1/S/G2M) -- the same method as the A3B G2/M reference, so the A3
tests are apples-to-apples with the literature.
[HOOK] ccAFv2 (comprehensive, adds a quiescent G0 state) to be wired in once its
2.0.6 API is confirmed; the block is stubbed and skipped for now.

Reuses patient_config. Run from scripts/PATIENT_SPECIFIC_EFFECTS/:
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
    print(f"[warn] ccAFv2 unavailable ({type(_e).__name__}); running Tirosh-only.", flush=True)

from patient_config import (
    PATIENT_COL, CELLTYPE_COL, DIR_00_DIAG, HIGH_CONTRIBUTORS,
    banner, log, ensure_dir, load_adata, load_three_groups,
)

# =============================================================================
# CONFIG
# =============================================================================
OUT_DIR   = ensure_dir(DIR_00_DIAG)
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")

SOURCE_COL   = 'source_name'
NORMAL_SOURCE = 'normal tissue adjucent to head and neck squamous cell carcinoma'
PHASE_ORDER  = ['G1', 'S', 'G2M']
POP_ORDER    = ['SBS2_HIGH', 'CNV_HIGH', 'NORMAL']
POP_LABELS   = {'SBS2_HIGH': 'SBS2-HIGH', 'CNV_HIGH': 'CNV-HIGH', 'NORMAL': 'Normal'}
FIG5_THRESHOLD = 2.0
PHASE_COLORS = {'G1': '#5B7C99', 'S': '#f4a259', 'G2M': '#ed6a5a'}
CCAFV2_ORDER = ['qG0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1']
CCAFV2_COLORS = {'qG0': '#3b5c7a', 'G1': '#5B7C99', 'Late G1': '#9bb4c7',
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
    log(msg); report_lines.append(str(msg))

def get_expr(ad, gene):
    if gene not in ad.var_names:
        return np.zeros(ad.n_obs)
    x = ad[:, gene].X
    return np.asarray(x.todense()).flatten() if scipy.sparse.issparse(x) else np.asarray(x).flatten()

def phase_pcts(mask):
    """Return dict phase->pct and the counts, for cells under boolean mask."""
    ph = basal.obs.loc[mask, 'phase']
    n = len(ph)
    counts = {p: int((ph == p).sum()) for p in PHASE_ORDER}
    pct = {p: (100.0 * counts[p] / n if n else 0.0) for p in PHASE_ORDER}
    return counts, pct, n

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
n_tumor = int((basal.obs['tissue_grp'] == 'tumor').sum())
n_normal_t = int((basal.obs['tissue_grp'] == 'normal').sum())
rlog(f"  Basal: {basal.n_obs:,}  (tumor {n_tumor:,}, normal-adjacent {n_normal_t:,})")

basal.obs['A3A'] = get_expr(basal, 'APOBEC3A')
basal.obs['A3B'] = get_expr(basal, 'APOBEC3B')

s_in  = [g for g in S_GENES if g in basal.var_names]
g2_in = [g for g in G2M_GENES if g in basal.var_names]
rlog(f"  Tirosh gene sets present: S {len(s_in)}/{len(S_GENES)}, "
     f"G2M {len(g2_in)}/{len(G2M_GENES)}")
sc.tl.score_genes_cell_cycle(basal, s_genes=s_in, g2m_genes=g2_in)
basal.obs['phase'] = basal.obs['phase'].astype(str)
overall = basal.obs['phase'].value_counts()
rlog(f"  Overall basal phase mix: "
     + ", ".join(f"{p} {100*overall.get(p,0)/basal.n_obs:.1f}%" for p in PHASE_ORDER))

# ccAFv2 comprehensive scoring (7 states incl. quiescent qG0) -- confirmation tool.
# Our var_names are gene symbols, so gene_id='symbol'. ccAFv2's missing-gene padding
# builds a DataFrame with set() columns, which pandas 2.x rejects ("columns cannot
# be a set"); we pre-pad any absent marker genes with zeros so that path never runs.
if CCAFV2_OK:
    try:
        # pass a copy: _prep_predict_data mutates its input (var_names_make_unique,
        # filter_genes). Labels align to basal.obs (cells are not dropped).
        res = ccAFv2.predict_labels(basal.copy(), species='human', gene_id='symbol')
        labels = res[0] if isinstance(res, tuple) else res      # returns (labels, probabilities)
        basal.obs['ccAFv2'] = np.asarray(labels).astype(str)
        ccc = basal.obs['ccAFv2'].value_counts()
        rlog("  ccAFv2 states: " + ", ".join(
            f"{s} {100*ccc.get(s,0)/basal.n_obs:.1f}%"
            for s in CCAFV2_ORDER if s in ccc.index))
        extra = [s for s in ccc.index if s not in CCAFV2_ORDER]
        if extra:
            rlog("  (other labels incl. below-threshold 'Unknown': "
                 + ", ".join(f"{s} {100*ccc.get(s,0)/basal.n_obs:.1f}%" for s in extra) + ")")
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
    _, pct, n = phase_pcts(basal.obs['group'] == p)
    rlog(f"  {POP_LABELS[p]:<12s} {n:>6d} {pct['G1']:>6.1f}% {pct['S']:>6.1f}% {pct['G2M']:>6.1f}%")
# chi-square on group x phase
ct = np.array([[int(((basal.obs['group'] == p) & (basal.obs['phase'] == ph)).sum())
                for ph in PHASE_ORDER] for p in POP_ORDER])
chi2, pv, dof, _ = chi2_contingency(ct)
rlog(f"\n  group x phase chi-square: chi2={chi2:.1f}, dof={dof}, p={pv:.2e}")

# =============================================================================
# STEP B: A3 CELL-CYCLE ASSOCIATION (all basal; clean)
# =============================================================================
banner("STEP B: A3A / A3B CELL-CYCLE ASSOCIATION (all basal)")

def a3_test(flag_name, pos_mask, target_phase):
    _, pct_pos, n_pos = phase_pcts(pos_mask)
    _, pct_neg, n_neg = phase_pcts(~pos_mask)
    a = int(((pos_mask) & (basal.obs['phase'] == target_phase)).sum())
    b = n_pos - a
    c = int(((~pos_mask) & (basal.obs['phase'] == target_phase)).sum())
    d = n_neg - c
    orr, pf = fisher_exact([[a, b], [c, d]])
    rlog(f"  {flag_name}+ (n={n_pos}) vs {flag_name}- (n={n_neg}), target {target_phase}:")
    rlog(f"    {flag_name}+ phase mix: " + ", ".join(f"{p} {pct_pos[p]:.1f}%" for p in PHASE_ORDER))
    rlog(f"    {flag_name}- phase mix: " + ", ".join(f"{p} {pct_neg[p]:.1f}%" for p in PHASE_ORDER))
    rlog(f"    {target_phase} enrichment: OR={orr:.2f}, Fisher p={pf:.2e}\n")

rlog("  Replication target (A3B -> G2/M) and extension (A3A -> G1):\n")
a3_test('A3B', basal.obs['A3B'] > 0, 'G2M')
a3_test('A3A', basal.obs['A3A'] > 0, 'G1')

# dominance version among A3-expressing cells (cleanest enzyme axis)
expr = (basal.obs['A3A'] + basal.obs['A3B']) > 0
frac = basal.obs['A3A'] / (basal.obs['A3A'] + basal.obs['A3B']).replace(0, np.nan)
a3a_dom = expr & (frac > 0.5)
a3b_dom = expr & (frac <= 0.5)
_, pa, na = phase_pcts(a3a_dom)
_, pb, nb = phase_pcts(a3b_dom)
rlog(f"  Among A3-expressing basal cells (dominance split):")
rlog(f"    A3A-dominant (n={na}): " + ", ".join(f"{p} {pa[p]:.1f}%" for p in PHASE_ORDER))
rlog(f"    A3B-dominant (n={nb}): " + ", ".join(f"{p} {pb[p]:.1f}%" for p in PHASE_ORDER))
ctd = np.array([[int((a3a_dom & (basal.obs['phase'] == ph)).sum()) for ph in PHASE_ORDER],
                [int((a3b_dom & (basal.obs['phase'] == ph)).sum()) for ph in PHASE_ORDER]])
if ctd.sum() > 0 and (ctd.sum(axis=1) > 0).all():
    c2, p2, _, _ = chi2_contingency(ctd)
    rlog(f"    dominance x phase chi-square: chi2={c2:.1f}, p={p2:.2e}")

# =============================================================================
# STEP B2: ccAFv2 CONFIRMATION (7-state, G0-aware)
# =============================================================================
if CCAFV2_OK and 'ccAFv2' in basal.obs.columns:
    banner("STEP B2: ccAFv2 CONFIRMATION (comprehensive, quiescence-resolved)")

    def cc_pcts(mask):
        s = basal.obs.loc[mask, 'ccAFv2']
        n = len(s)
        return {st: (100.0 * (s == st).sum() / n if n else 0.0) for st in CCAFV2_ORDER}, n

    rlog("  Per-group ccAFv2 state mix:")
    for p in POP_ORDER:
        d, n = cc_pcts(basal.obs['group'] == p)
        rlog(f"    {POP_LABELS[p]:<10s} (n={n}): "
             + ", ".join(f"{st} {d[st]:.0f}%" for st in CCAFV2_ORDER if d[st] >= 1))

    da, na = cc_pcts(a3a_dom)
    db, nb = cc_pcts(a3b_dom)
    rlog("\n  By A3 dominance (does A3A carry qG0 quiescence Tirosh lumped into G1?):")
    rlog(f"    A3A-dominant (n={na}): "
         + ", ".join(f"{st} {da[st]:.0f}%" for st in CCAFV2_ORDER if da[st] >= 1))
    rlog(f"    A3B-dominant (n={nb}): "
         + ", ".join(f"{st} {db[st]:.0f}%" for st in CCAFV2_ORDER if db[st] >= 1))
    qg1 = lambda d: d['qG0'] + d['G1'] + d['Late G1']
    sg2 = lambda d: d['S'] + d['S/G2'] + d['G2/M']
    rlog(f"    A3A-dom quiescent+G1 = {qg1(da):.0f}% vs S+G2/M = {sg2(da):.0f}%")
    rlog(f"    A3B-dom quiescent+G1 = {qg1(db):.0f}% vs S+G2/M = {sg2(db):.0f}%")
    rlog(f"    qG0 specifically: A3A-dom {da['qG0']:.0f}% vs A3B-dom {db['qG0']:.0f}%")
else:
    da = db = None

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
                 'tumor_G1': pct_t['G1'], 'tumor_S': pct_t['S'], 'tumor_G2M': pct_t['G2M'],
                 'norm_G1': pct_n['G1'], 'norm_S': pct_n['S'], 'norm_G2M': pct_n['G2M']})
pc = pd.DataFrame(rows)
rlog(f"  {'Patient':<10s} {'nTum':>6s} {'nNorm':>6s} | tumor G1/S/G2M | normal-adj G1/S/G2M")
for _, r in pc.iterrows():
    rlog(f"  {str(r['patient']).replace('Patient ',''):<10s} {int(r['n_tumor']):>6d} "
         f"{int(r['n_normal_adj']):>6d} | {r['tumor_G1']:.0f}/{r['tumor_S']:.0f}/{r['tumor_G2M']:.0f} "
         f"| {r['norm_G1']:.0f}/{r['norm_S']:.0f}/{r['norm_G2M']:.0f}")

# =============================================================================
# STEP D: SOURCE-MATCHED FOLDS + HIGH-CONTRIBUTOR PHASE SKEW
# =============================================================================
banner("STEP D: SOURCE-MATCHED CONTRIBUTION FOLDS + HC PHASE SKEW")
total_basal = basal.n_obs
n_s, n_c = 546, 546
frows = []
for p in patients:
    pm = basal.obs[PATIENT_COL] == p
    n_basal_p = int(pm.sum())
    n_tum_p = int((pm & (basal.obs['tissue_grp'] == 'tumor')).sum())
    n_sbs2_p = int((pm & (basal.obs['group'] == 'SBS2_HIGH')).sum())
    n_cnv_p  = int((pm & (basal.obs['group'] == 'CNV_HIGH')).sum())
    # old (all-basal) vs new (tumor-matched) folds
    fold_s_all = (n_sbs2_p/n_s) / (n_basal_p/total_basal) if n_basal_p else 0
    fold_c_all = (n_cnv_p/n_c)  / (n_basal_p/total_basal) if n_basal_p else 0
    fold_s_tm  = (n_sbs2_p/n_s) / (n_tum_p/n_tumor) if n_tum_p else 0
    fold_c_tm  = (n_cnv_p/n_c)  / (n_tum_p/n_tumor) if n_tum_p else 0
    frows.append({'patient': p, 'fold_sbs2_all': fold_s_all, 'fold_sbs2_tumor': fold_s_tm,
                  'fold_cnv_all': fold_c_all, 'fold_cnv_tumor': fold_c_tm})
fd = pd.DataFrame(frows)
rlog(f"  {'Patient':<10s} {'SBS2 all->tumor':>18s}   {'CNV all->tumor':>18s}")
for _, r in fd.sort_values('fold_cnv_tumor', ascending=False).iterrows():
    rlog(f"  {str(r['patient']).replace('Patient ',''):<10s} "
         f"{r['fold_sbs2_all']:>7.2f}x ->{r['fold_sbs2_tumor']:>6.2f}x   "
         f"{r['fold_cnv_all']:>7.2f}x ->{r['fold_cnv_tumor']:>6.2f}x")

sbs2_hc_tm = set(fd[fd['fold_sbs2_tumor'] >= FIG5_THRESHOLD]['patient'])
cnv_hc_tm  = set(fd[fd['fold_cnv_tumor']  >= FIG5_THRESHOLD]['patient'])
rlog(f"\n  Source-matched high contributors (tumor-denominator, >= {FIG5_THRESHOLD}x):")
rlog(f"    SBS2: {[str(p).replace('Patient ','') for p in sbs2_hc_tm]}")
rlog(f"    CNV:  {[str(p).replace('Patient ','') for p in cnv_hc_tm]}")
rlog(f"    (prior all-basal SBS2 HC were: {[str(p).replace('Patient ','') for p in HIGH_CONTRIBUTORS]})")

# phase skew of HC tumor basal cells vs other tumor basal cells
hc_all = sbs2_hc_tm | cnv_hc_tm
hc_mask = basal.obs[PATIENT_COL].isin(hc_all) & (basal.obs['tissue_grp'] == 'tumor')
other_mask = (~basal.obs[PATIENT_COL].isin(hc_all)) & (basal.obs['tissue_grp'] == 'tumor')
_, pct_hc, nhc = phase_pcts(hc_mask)
_, pct_ot, notr = phase_pcts(other_mask)
rlog(f"\n  Tumor basal phase mix, HIGH contributors vs other patients:")
rlog(f"    HC patients   (n={nhc}): " + ", ".join(f"{p} {pct_hc[p]:.1f}%" for p in PHASE_ORDER))
rlog(f"    other patients(n={notr}): " + ", ".join(f"{p} {pct_ot[p]:.1f}%" for p in PHASE_ORDER))

# =============================================================================
# STEP E: PLOTS + SAVE
# =============================================================================
banner("STEP E: PLOTS")
pc.to_csv(os.path.join(OUT_DIR, "patient_cellcycle_by_source.tsv"), sep='\t', index=False)
fd.to_csv(os.path.join(OUT_DIR, "patient_source_matched_folds.tsv"), sep='\t', index=False)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
# group baseline
x = np.arange(len(POP_ORDER)); left = np.zeros(len(POP_ORDER))
for ph in PHASE_ORDER:
    vals = np.array([phase_pcts(basal.obs['group'] == p)[1][ph] for p in POP_ORDER])
    ax1.barh(x, vals, left=left, color=PHASE_COLORS[ph], label=ph, edgecolor='white')
    left += vals
ax1.set_yticks(x); ax1.set_yticklabels([POP_LABELS[p] for p in POP_ORDER])
ax1.set_xlabel('% of cells'); ax1.set_title('Cell-cycle mix by group'); ax1.legend(fontsize=9)
# A3 dominance
labels = ['A3A-dominant', 'A3B-dominant']
left2 = np.zeros(2)
for ph in PHASE_ORDER:
    vals = np.array([pa[ph], pb[ph]])
    ax2.barh([0, 1], vals, left=left2, color=PHASE_COLORS[ph], label=ph, edgecolor='white')
    left2 += vals
ax2.set_yticks([0, 1]); ax2.set_yticklabels(labels)
ax2.set_xlabel('% of cells'); ax2.set_title('Cell-cycle mix by A3 dominance'); ax2.legend(fontsize=9)
plt.tight_layout()
for ext in ('pdf', 'png'):
    plt.savefig(os.path.join(OUT_DIR, f"patient_cellcycle_by_source.{ext}"), dpi=DPI, bbox_inches='tight')
plt.close()
rlog("  [SAVE] patient_cellcycle_by_source.pdf/.png + 2 TSVs")

if CCAFV2_OK and 'ccAFv2' in basal.obs.columns and da is not None:
    figc, (axc1, axc2) = plt.subplots(1, 2, figsize=(14, 5))
    xg = np.arange(len(POP_ORDER)); leftg = np.zeros(len(POP_ORDER))
    for st in CCAFV2_ORDER:
        vals = np.array([100.0 * ((basal.obs['group'] == p) & (basal.obs['ccAFv2'] == st)).sum()
                         / max((basal.obs['group'] == p).sum(), 1) for p in POP_ORDER])
        axc1.barh(xg, vals, left=leftg, color=CCAFV2_COLORS[st], label=st, edgecolor='white')
        leftg += vals
    axc1.set_yticks(xg); axc1.set_yticklabels([POP_LABELS[p] for p in POP_ORDER])
    axc1.set_xlabel('% of cells'); axc1.set_title('ccAFv2 states by group')
    axc1.legend(fontsize=7, ncol=2, loc='lower right')
    leftd = np.zeros(2)
    for st in CCAFV2_ORDER:
        vals = np.array([da[st], db[st]])
        axc2.barh([0, 1], vals, left=leftd, color=CCAFV2_COLORS[st], label=st, edgecolor='white')
        leftd += vals
    axc2.set_yticks([0, 1]); axc2.set_yticklabels(['A3A-dominant', 'A3B-dominant'])
    axc2.set_xlabel('% of cells'); axc2.set_title('ccAFv2 states by A3 dominance')
    axc2.legend(fontsize=7, ncol=2, loc='lower right')
    plt.tight_layout()
    for ext in ('pdf', 'png'):
        plt.savefig(os.path.join(OUT_DIR, f"patient_cellcycle_ccAFv2.{ext}"), dpi=DPI, bbox_inches='tight')
    plt.close()
    rlog("  [SAVE] patient_cellcycle_ccAFv2.pdf/.png")

report_path = os.path.join(OUT_DIR, f"patient_cellcycle_by_source_{TIMESTAMP}.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
rlog(f"\n  Report: {report_path}")
banner("CELL-CYCLE BY SOURCE DIAGNOSTIC COMPLETE")
