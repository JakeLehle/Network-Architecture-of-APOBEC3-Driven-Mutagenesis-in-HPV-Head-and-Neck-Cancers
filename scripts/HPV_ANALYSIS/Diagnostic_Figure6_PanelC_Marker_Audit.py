#!/usr/bin/env python3
"""
Diagnostic_Figure6_PanelC_Marker_Audit.py
==========================================
READ-ONLY audit of the Figure 6 Panel C host-marker dot plot.

Purpose
-------
Decide, for every gene currently in Panel C plus a set of candidate additions
drawn from the BRD4 / HPV16 lifecycle review, whether each gene (a) pulls its
weight in the dot plot and (b) sits in the right category. Nothing here edits
the figure generator. It writes a table + report we use to revise Panel C by
hand.

Universe
--------
The three canonical basal populations from Figure 4 (SBS2_HIGH, CNV_HIGH,
NORMAL), 546 cells each, UNGATED. Host expression is not gated to HPV16+; only
Panel F is. Group membership is read from three_group_assignments.tsv and mapped
onto adata exactly as Diagnostic_Figure6_HostMarkers_and_IntegrationProxy.py.

Weight flags (thresholds are first-pass; revisit after seeing output)
--------------------------------------------------------------------
  clutter_flag : expressed in < 10% of cells in ALL three populations
  flat_flag    : BH q >= 0.05 in ALL three pairwise contrasts
  verdict      : DROP (both flags), REVIEW (one flag), KEEP (neither)

Category fit
------------
  1. Peak-population + expected-peak direction check per gene.
  2. Data-driven clustering of genes by their z-scored [SBS2, CNV, NORMAL]
     expression pattern (ward, scipy), cross-tabbed against assigned category.
     A gene whose cluster majority category differs from its assigned category
     is a re-assignment candidate.

Stats
-----
  Mann-Whitney U per pairwise contrast (matches Panels B / D), Benjamini-Hochberg
  across the whole panel (all found genes x 3 contrasts) as one family.

Inputs (read-only)
------------------
  data/FIG_4/00_input/adata_final.h5ad
  data/FIG_4/01_group_selection/three_group_assignments.tsv

Outputs (to data/FIG_6/DIAGNOSTIC_PANELC_AUDIT/)
------------------------------------------------
  panelC_marker_audit.tsv           one row per gene, all metrics + verdict
  panelC_category_coherence.tsv     per-category KEEP/REVIEW/DROP + cluster spread
  panelC_proposed_panel.txt         proposed revised gene list with justification
  panelC_marker_audit_report.txt    full console log
  panelC_pattern_heatmap.pdf/.png   z-scored pattern heatmap, clustered

Run in the NETWORK conda env, from scripts/HPV_ANALYSIS/:
    conda run -n NETWORK python Diagnostic_Figure6_PanelC_Marker_Audit.py

Author: Jake Lehle / Claude (2026 NMF Paper)
Texas Biomedical Research Institute
"""

import os
from collections import OrderedDict, Counter
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from scipy.stats import mannwhitneyu
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT     = "/master/jlehle/WORKING/2026_NMF_PAPER"
ADATA_PATH       = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/adata_final.h5ad")
THREE_GROUP_PATH = os.path.join(PROJECT_ROOT,
                                "data/FIG_4/01_group_selection/three_group_assignments.tsv")
OUTPUT_DIR       = os.path.join(PROJECT_ROOT, "data/FIG_6/DIAGNOSTIC_PANELC_AUDIT")
os.makedirs(OUTPUT_DIR, exist_ok=True)

POP_ORDER  = ['SBS2_HIGH', 'CNV_HIGH', 'NORMAL']
POP_LABELS = {'SBS2_HIGH': 'SBS2-HIGH', 'CNV_HIGH': 'CNV-HIGH', 'NORMAL': 'Normal'}
EXPECTED_N = 546   # per population; anything else => stale / wrong file

# weight-flag thresholds (first pass; revisit after output)
CLUTTER_MAX_FRAC = 0.10   # < 10% expressing in ALL three populations
FLAT_Q           = 0.05   # q >= this in ALL three contrasts

# clustering
N_CLUSTERS = 6

DPI = 300

# figure hexes where the three groups appear
POP_COLORS = {'SBS2_HIGH': '#ed6a5a', 'CNV_HIGH': '#F6D155', 'NORMAL': '#5B7C99'}
VERDICT_COLORS = {'KEEP': '#3a9d5d', 'REVIEW': '#e6a020', 'DROP': '#cc3b3b'}

# =============================================================================
# GENE PANELS
# =============================================================================
# Current Panel C, exact figure order (Generate_Figure6_Lifecycle_Panels.py)
CURRENT_PANEL = OrderedDict([
    ('APOBEC/Innate Immune',
        ['APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D',
         'APOBEC3F', 'APOBEC3G', 'APOBEC3H', 'CGAS', 'STING1']),
    ('Immune Signaling',
        ['STAT1', 'HLA-A', 'HLA-B', 'HLA-C', 'IRF1', 'TAP1', 'B2M']),
    ('Differentiation',
        ['KRT5', 'KRT14', 'KRT1', 'KRT10', 'CDH1', 'IVL']),
    ('ATM/DNA Damage',
        ['ATM', 'CHEK2', 'CHEK1', 'BRCA1', 'MRE11', 'RAD50', 'NBN',
         'H2AX', 'STAT5A', 'STAT5B']),
    ('Transformation/Proliferation',
        ['CDKN2A', 'MCM7', 'CCNE1', 'MKI67', 'TOP2A', 'PCNA',
         'BRD4', 'MED1', 'E2F1', 'E2F2']),
    ('p53/Rb Pathway',
        ['CDKN1A', 'BAX', 'MDM2', 'RB1', 'TP53']),
    ('G2/M Arrest',
        ['CDC25A', 'CDC25C', 'CDK1', 'CCNB1']),
])

# Candidate additions, grouped by the category we would slot them into.
# PCNA is already in Transformation, so it is intentionally NOT repeated in NER.
CANDIDATE_ADDITIONS = OrderedDict([
    # dropped from the figure once (present in the diagnostic, absent from Panel C)
    ('Caspase (re-add?)',
        ['CASP3', 'CASP7']),
    # BRD4-recruited DDR arm; predicted to track CNV-HIGH with BRCA1/H2AX/CHEK2
    ('ATM/DNA Damage (+BRD4 arm)',
        ['BARD1', 'TP53BP1', 'RIF1', 'NSD2']),
    # exploratory set A: BRD4-L NER network
    ('NER (exploratory)',
        ['XPC', 'RFC1', 'CUL4A', 'XAB2']),
    # exploratory set A: chromatin / remodeling
    ('Chromatin remodeling (exploratory)',
        ['INO80', 'MCRS1', 'ACTR8']),
    # exploratory set B: SWI/SNF host-silencing (ties to Section 4.5 escape)
    ('SWI/SNF silencing (exploratory)',
        ['SMARCA4', 'SMARCB1', 'ARID1A', 'BICRA']),
])

# Expected peak population per category, from the manuscript's Question 4 findings.
# None = no single directional expectation (report the observed peak instead).
EXPECTED_PEAK = {
    'APOBEC/Innate Immune':               None,        # A3A->SBS2, A3B->CNV; mixed
    'Immune Signaling':                   'SBS2_HIGH', # SBS2-HIGH is immune-visible
    'Differentiation':                    None,        # basal KRT5/14->CNV, terminal IVL->SBS2
    'ATM/DNA Damage':                     'CNV_HIGH',  # ATM-dependent DDR in productive
    'Transformation/Proliferation':       'CNV_HIGH',  # proliferation higher in CNV-HIGH
    'p53/Rb Pathway':                     None,        # mixed
    'G2/M Arrest':                        'CNV_HIGH',  # G2/M arrest in CNV-HIGH
    'Caspase (re-add?)':                  None,
    'ATM/DNA Damage (+BRD4 arm)':         'CNV_HIGH',  # prediction under test
    'NER (exploratory)':                  None,
    'Chromatin remodeling (exploratory)': None,
    'SWI/SNF silencing (exploratory)':    'CNV_HIGH',  # soft prediction under test
}

# Symbol aliases (primary first). resolve_gene tries each in turn.
ALIASES = {
    'TP53BP1': ['TP53BP1', '53BP1'],
    'NSD2':    ['NSD2', 'WHSC1', 'MMSET'],
    'BICRA':   ['BICRA', 'GLTSCR1'],
    'H2AX':    ['H2AX', 'H2AFX'],
    'MRE11':   ['MRE11', 'MRE11A'],
    'CGAS':    ['CGAS', 'MB21D1', 'C6orf150'],
    'STING1':  ['STING1', 'TMEM173'],
    'SMARCA4': ['SMARCA4', 'BRG1'],
    'SMARCB1': ['SMARCB1', 'SNF5', 'INI1', 'BAF47'],
    'ARID1A':  ['ARID1A', 'BAF250A'],
    'XAB2':    ['XAB2', 'SYF1', 'HCNP'],
    'MCRS1':   ['MCRS1', 'MSP58'],
    'NBN':     ['NBN', 'NBS1'],
}

# =============================================================================
# LOGGING
# =============================================================================
report_lines = []

def log(msg=""):
    print(msg, flush=True)
    report_lines.append(str(msg))

def banner(title):
    log("")
    log("=" * 80)
    log(f"  {title}")
    log("=" * 80)

# =============================================================================
# HELPERS
# =============================================================================
def get_expression(adata, symbol):
    """Extract expression vector for a gene symbol (var_names, then gene_symbol col)."""
    if symbol in adata.var_names:
        idx = adata.var_names.get_loc(symbol)
        x = adata.X[:, idx]
        return (np.asarray(x.todense()).flatten()
                if scipy.sparse.issparse(x) else np.asarray(x).flatten())
    if 'gene_symbol' in adata.var.columns:
        mask = (adata.var['gene_symbol'] == symbol).values
        if mask.any():
            idx = int(np.where(mask)[0][0])
            x = adata.X[:, idx]
            return (np.asarray(x.todense()).flatten()
                    if scipy.sparse.issparse(x) else np.asarray(x).flatten())
    return None

def resolve_gene(adata, symbol):
    """Return (matched_symbol, values), trying the primary symbol then any aliases."""
    candidates = ALIASES.get(symbol, [symbol])
    if symbol not in candidates:
        candidates = [symbol] + candidates
    for s in candidates:
        v = get_expression(adata, s)
        if v is not None:
            return s, v
    return None, None

def bh_adjust(pvals):
    """
    Benjamini-Hochberg FDR, NaN-aware, values returned in original order.

    Corrected form (guards against the reversal bug): rank ascending, compute
    p*n/rank, then enforce monotone non-decreasing q along ascending p via a
    cumulative minimum taken from the LARGEST p downward, and scatter the result
    back to the original positions. The value array is never left reversed.
    """
    p = np.asarray(pvals, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    n = int(ok.sum())
    if n == 0:
        return q
    idx = np.where(ok)[0]
    pv = p[idx]
    order = np.argsort(pv)                       # ascending
    ranked = pv[order]
    raw = ranked * n / np.arange(1, n + 1)
    q_sorted = np.minimum.accumulate(raw[::-1])[::-1]   # monotone from the top down
    q_sorted = np.clip(q_sorted, 0, 1)
    q_back = np.empty(n)
    q_back[order] = q_sorted
    q[idx] = q_back
    return q

def pairwise_mwu(vbp):
    """Three pairwise MWU raw p-values: (SBS2,CNV), (CNV,NORM), (SBS2,NORM)."""
    pairs = [('SBS2_HIGH', 'CNV_HIGH'), ('CNV_HIGH', 'NORMAL'), ('SBS2_HIGH', 'NORMAL')]
    out = []
    for a, b in pairs:
        va, vb = vbp[a], vbp[b]
        if len(va) > 5 and len(vb) > 5:
            try:
                _, pp = mannwhitneyu(va, vb, alternative='two-sided')
            except ValueError:            # e.g. all values identical
                pp = np.nan
            out.append(pp)
        else:
            out.append(np.nan)
    return out

def fmt(x, nd=3):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return 'N.D.'
    return f"{x:.{nd}f}"

def fmt_q(x):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return 'N.D.'
    if x < 1e-3:
        return f"{x:.1e}"
    return f"{x:.3f}"

# =============================================================================
# STEP 0: LOAD (READ-ONLY)
# =============================================================================
banner("STEP 0: Load groups + adata (READ-ONLY)")

groups = pd.read_csv(THREE_GROUP_PATH, sep='\t')
counts = groups['group'].value_counts().to_dict()
log(f"  three_group_assignments: {len(groups)} cells")
for p in POP_ORDER:
    n = counts.get(p, 0)
    flag = "" if n == EXPECTED_N else f"   <-- expected {EXPECTED_N}"
    log(f"    {p:<10s}: {n}{flag}")
if any(counts.get(p, 0) != EXPECTED_N for p in POP_ORDER):
    log("  WARNING: population sizes deviate from the canonical 546/546/546. "
        "Confirm this is the current three_group_assignments.tsv, not a stale copy.")

log("\n  Loading adata_final.h5ad ...")
adata = sc.read_h5ad(ADATA_PATH)
log(f"  adata: {adata.n_obs} cells x {adata.n_vars} genes")

# sanity: confirm .X is log-normalized (dot-plot color is mean log-norm expression)
xmax = adata.X.max()
try:
    xmax = float(xmax)
except Exception:
    xmax = float(np.asarray(xmax).flatten()[0])
verdict_x = "looks log-normalized" if xmax < 20 else "LOOKS LIKE RAW COUNTS -- CHECK"
log(f"  adata.X max value = {xmax:.2f}  ({verdict_x}; means assume log-normalized)")

# map populations exactly as the host-marker diagnostic does
adata.obs['population'] = 'other'
for p in POP_ORDER:
    cells = set(groups.loc[groups['group'] == p, 'cell_barcode'])
    adata.obs.loc[adata.obs_names.isin(cells), 'population'] = p
adata_pop = adata[adata.obs['population'].isin(POP_ORDER)].copy()
log(f"\n  cells mapped into the three populations: {adata_pop.n_obs}")
for p in POP_ORDER:
    log(f"    {POP_LABELS[p]:<10s}: {(adata_pop.obs['population'] == p).sum()}")

pop_masks = {p: (adata_pop.obs['population'] == p).values for p in POP_ORDER}

# =============================================================================
# STEP 1: RESOLVE GENES + PER-POPULATION MEAN / FRACTION
# =============================================================================
banner("STEP 1: Resolve genes + per-population mean / fraction expressing")

assigned = OrderedDict()   # gene -> category
source   = {}              # gene -> 'current' | 'candidate'
for cat, glist in CURRENT_PANEL.items():
    for g in glist:
        assigned[g] = cat
        source[g] = 'current'
for cat, glist in CANDIDATE_ADDITIONS.items():
    for g in glist:
        if g in assigned:
            log(f"  NOTE: candidate {g} already in current panel "
                f"({assigned[g]}); not duplicated")
            continue
        assigned[g] = cat
        source[g] = 'candidate'

expr_means = {}     # gene -> {pop: mean}
expr_fracs = {}     # gene -> {pop: frac expressing}
vbp_by_gene = {}    # gene -> {pop: values}
matched_symbol = {}
not_found = []

for g in assigned:
    msym, vals = resolve_gene(adata_pop, g)
    if vals is None:
        not_found.append(g)
        continue
    matched_symbol[g] = msym
    means, fracs, vbp = {}, {}, {}
    for p in POP_ORDER:
        vv = vals[pop_masks[p]]
        vbp[p] = vv
        means[p] = float(np.mean(vv)) if len(vv) else np.nan
        fracs[p] = float(np.mean(vv > 0)) if len(vv) else np.nan
    expr_means[g] = means
    expr_fracs[g] = fracs
    vbp_by_gene[g] = vbp
    note = "" if msym == g else f"  (matched via alias '{msym}')"
    log(f"    {g:<10s} [{source[g]:<9s}] mean SBS2/CNV/NORM = "
        f"{fmt(means['SBS2_HIGH'])}/{fmt(means['CNV_HIGH'])}/{fmt(means['NORMAL'])}"
        f"{note}")

if not_found:
    log(f"\n  Genes NOT found in adata (excluded from stats + clustering): {not_found}")
    log("  Any wanted addition that is missing here is a data-availability issue to note.")

found_genes = [g for g in assigned if g in vbp_by_gene]
log(f"\n  Genes carried into stats: {len(found_genes)}  "
    f"(current: {sum(source[g]=='current' for g in found_genes)}, "
    f"candidate: {sum(source[g]=='candidate' for g in found_genes)})")

# =============================================================================
# STEP 2: PAIRWISE MANN-WHITNEY + PANEL-WIDE BH + WEIGHT FLAGS
# =============================================================================
banner("STEP 2: Pairwise Mann-Whitney + panel-wide BH + weight flags")

raw_p = []
for g in found_genes:
    raw_p.extend(pairwise_mwu(vbp_by_gene[g]))
q_all = bh_adjust(raw_p)
n_valid = int(np.sum(~np.isnan(raw_p)))
log(f"  {len(found_genes)} genes x 3 contrasts = {len(raw_p)} tests; "
    f"{n_valid} valid, BH-corrected as one family")

audit_rows = []
for i, g in enumerate(found_genes):
    means = expr_means[g]
    fracs = expr_fracs[g]
    qs = q_all[3 * i: 3 * i + 3]                       # SBS2vCNV, CNVvNORM, SBS2vNORM
    max_frac = float(np.nanmax([fracs[p] for p in POP_ORDER]))
    valid_q = [x for x in qs if not np.isnan(x)]
    min_q = min(valid_q) if valid_q else np.nan
    clutter = bool(max_frac < CLUTTER_MAX_FRAC)
    flat = bool(len(valid_q) == 0 or all(x >= FLAT_Q for x in valid_q))
    verdict = 'DROP' if (clutter and flat) else ('REVIEW' if (clutter or flat) else 'KEEP')
    peak = max(POP_ORDER,
               key=lambda p: (means[p] if not np.isnan(means[p]) else -np.inf))
    exp_peak = EXPECTED_PEAK.get(assigned[g])
    direction_ok = None if exp_peak is None else bool(peak == exp_peak)
    audit_rows.append({
        'gene': g,
        'matched_symbol': matched_symbol[g],
        'source': source[g],
        'assigned_category': assigned[g],
        'mean_SBS2': means['SBS2_HIGH'], 'mean_CNV': means['CNV_HIGH'], 'mean_NORM': means['NORMAL'],
        'frac_SBS2': fracs['SBS2_HIGH'], 'frac_CNV': fracs['CNV_HIGH'], 'frac_NORM': fracs['NORMAL'],
        'max_frac': max_frac,
        'q_SBS2vCNV': qs[0], 'q_CNVvNORM': qs[1], 'q_SBS2vNORM': qs[2],
        'min_q': min_q,
        'peak_population': peak,
        'expected_peak': exp_peak if exp_peak is not None else '',
        'direction_ok': direction_ok if direction_ok is not None else '',
        'clutter_flag': clutter,
        'flat_flag': flat,
        'verdict': verdict,
    })

audit = pd.DataFrame(audit_rows)
vc = audit['verdict'].value_counts().to_dict()
log(f"  Verdicts: KEEP={vc.get('KEEP',0)}, REVIEW={vc.get('REVIEW',0)}, DROP={vc.get('DROP',0)}")

# =============================================================================
# STEP 3: DATA-DRIVEN PATTERN CLUSTERING
# =============================================================================
banner("STEP 3: Data-driven pattern clustering (z-scored [SBS2, CNV, NORMAL])")

M = audit[['mean_SBS2', 'mean_CNV', 'mean_NORM']].values.astype(float)
Z = np.zeros_like(M)
for i in range(M.shape[0]):
    row = M[i]
    mu = np.nanmean(row)
    sd = np.nanstd(row)
    Z[i] = (row - mu) / sd if sd > 1e-9 else 0.0

k = min(N_CLUSTERS, max(2, len(Z) - 1))
if len(Z) > k:
    link = linkage(Z, method='ward')
    cl = fcluster(link, t=k, criterion='maxclust')
    leaf_order = leaves_list(link)
else:
    link = None
    cl = np.ones(len(Z), dtype=int)
    leaf_order = np.arange(len(Z))

audit['cluster'] = cl

# describe each cluster by its mean pattern + majority assigned category
cluster_major, cluster_peak = {}, {}
for c in np.unique(cl):
    sub = audit[audit['cluster'] == c]
    cluster_major[c] = Counter(sub['assigned_category']).most_common(1)[0][0]
    zmean = sub[['mean_SBS2', 'mean_CNV', 'mean_NORM']].mean().values
    cluster_peak[c] = POP_ORDER[int(np.argmax(zmean))]
audit['cluster_majority_category'] = audit['cluster'].map(cluster_major)
audit['cluster_peak_pop'] = audit['cluster'].map(cluster_peak)
audit['category_mismatch'] = audit['assigned_category'] != audit['cluster_majority_category']

for c in np.unique(cl):
    sub = audit[audit['cluster'] == c]
    log(f"\n  Cluster {c} (n={len(sub)}, peaks at {POP_LABELS[cluster_peak[c]]}, "
        f"majority '{cluster_major[c]}'):")
    for _, r in sub.iterrows():
        tag = "  <-- off-category" if r['category_mismatch'] else ""
        log(f"    {r['gene']:<10s}  [{r['assigned_category']}]{tag}")

# =============================================================================
# STEP 4: PER-CATEGORY COHERENCE
# =============================================================================
banner("STEP 4: Per-category coherence")

cat_order = list(CURRENT_PANEL.keys()) + \
            [c for c in CANDIDATE_ADDITIONS if c not in CURRENT_PANEL]
cat_rows = []
for cat in cat_order:
    sub = audit[audit['assigned_category'] == cat]
    if len(sub) == 0:
        continue
    n = len(sub)
    keep = int((sub['verdict'] == 'KEEP').sum())
    review = int((sub['verdict'] == 'REVIEW').sum())
    drop = int((sub['verdict'] == 'DROP').sum())
    exp_peak = EXPECTED_PEAK.get(cat)
    if exp_peak is not None:
        dir_ok = int((sub['direction_ok'] == True).sum())
        dir_str = f"{dir_ok}/{n} peak at {POP_LABELS[exp_peak]} (expected)"
    else:
        modal = Counter(sub['peak_population']).most_common(1)[0]
        dir_str = f"no fixed expectation; modal peak {POP_LABELS[modal[0]]} ({modal[1]}/{n})"
    n_span = int(sub['cluster'].nunique())
    mism = int(sub['category_mismatch'].sum())
    cat_rows.append({
        'category': cat, 'n_genes': n,
        'keep': keep, 'review': review, 'drop': drop,
        'direction': dir_str,
        'clusters_spanned': n_span, 'off_category_genes': mism,
    })
    log(f"\n  {cat}  (n={n})")
    log(f"    verdicts     : KEEP {keep} / REVIEW {review} / DROP {drop}")
    log(f"    direction    : {dir_str}")
    log(f"    cluster spread: spans {n_span} cluster(s), {mism} off-category gene(s)")
    if drop or review:
        weak = sub[sub['verdict'] != 'KEEP']
        for _, r in weak.iterrows():
            reasons = []
            if r['clutter_flag']:
                reasons.append(f"clutter (max frac {r['max_frac']:.2f})")
            if r['flat_flag']:
                reasons.append("flat (no significant contrast)")
            log(f"      - {r['gene']}: {r['verdict']} :: {', '.join(reasons)}")

cat_summary = pd.DataFrame(cat_rows)

# =============================================================================
# STEP 5: PROPOSED REVISED PANEL
# =============================================================================
banner("STEP 5: Proposed revised Panel C")

prop_lines = []
def pline(s=""):
    prop_lines.append(s)

pline("PROPOSED PANEL C REVISION")
pline("=" * 60)
pline(f"Generated {datetime.now():%Y-%m-%d %H:%M:%S}")
pline("Thresholds: clutter < %.0f%% expressing in all pops; flat q >= %.2f in all contrasts"
      % (100 * CLUTTER_MAX_FRAC, FLAT_Q))
pline("")

cur = audit[audit['source'] == 'current']
cand = audit[audit['source'] == 'candidate']

pline("CURRENT GENES TO DROP (both flags):")
drop_cur = cur[cur['verdict'] == 'DROP']
if len(drop_cur):
    for _, r in drop_cur.iterrows():
        pline(f"  - {r['gene']} [{r['assigned_category']}]  "
              f"max frac {r['max_frac']:.2f}, min q {fmt_q(r['min_q'])}")
else:
    pline("  (none)")

pline("")
pline("CURRENT GENES TO REVIEW (one flag):")
rev_cur = cur[cur['verdict'] == 'REVIEW']
if len(rev_cur):
    for _, r in rev_cur.iterrows():
        why = "clutter" if r['clutter_flag'] else "flat"
        pline(f"  - {r['gene']} [{r['assigned_category']}]  ({why}); "
              f"peaks {POP_LABELS[r['peak_population']]}")
else:
    pline("  (none)")

pline("")
pline("CANDIDATE ADDITIONS THAT EARN A PLACE (verdict KEEP):")
add_ok = cand[cand['verdict'] == 'KEEP']
if len(add_ok):
    for _, r in add_ok.iterrows():
        tgt = (r['assigned_category'] if not r['category_mismatch']
               else f"{r['assigned_category']} -> consider '{r['cluster_majority_category']}' "
                    f"(clusters with {POP_LABELS[r['cluster_peak_pop']]} peak)")
        pline(f"  + {r['gene']}  into {tgt};  "
              f"peaks {POP_LABELS[r['peak_population']]}, min q {fmt_q(r['min_q'])}")
else:
    pline("  (none clear enough to add)")

pline("")
pline("CANDIDATE ADDITIONS THAT DO NOT EARN A PLACE:")
add_no = cand[cand['verdict'] != 'KEEP']
if len(add_no):
    for _, r in add_no.iterrows():
        reasons = []
        if r['clutter_flag']:
            reasons.append(f"clutter (max frac {r['max_frac']:.2f})")
        if r['flat_flag']:
            reasons.append("flat")
        pline(f"  x {r['gene']} [{r['assigned_category']}]  {r['verdict']} :: {', '.join(reasons)}")
else:
    pline("  (none)")

pline("")
pline("RE-ASSIGNMENT CANDIDATES (gene clusters away from its labelmates):")
mism_all = audit[audit['category_mismatch']]
if len(mism_all):
    for _, r in mism_all.iterrows():
        pline(f"  ~ {r['gene']}: assigned '{r['assigned_category']}' but clusters with "
              f"'{r['cluster_majority_category']}' ({POP_LABELS[r['cluster_peak_pop']]} peak)")
else:
    pline("  (none)")

for s in prop_lines:
    log(s)

# =============================================================================
# STEP 6: WRITE TABLES + REPORT
# =============================================================================
banner("STEP 6: Write outputs")

col_order = [
    'gene', 'matched_symbol', 'source', 'assigned_category',
    'mean_SBS2', 'mean_CNV', 'mean_NORM',
    'frac_SBS2', 'frac_CNV', 'frac_NORM', 'max_frac',
    'q_SBS2vCNV', 'q_CNVvNORM', 'q_SBS2vNORM', 'min_q',
    'peak_population', 'expected_peak', 'direction_ok',
    'clutter_flag', 'flat_flag', 'verdict',
    'cluster', 'cluster_peak_pop', 'cluster_majority_category', 'category_mismatch',
]
audit_out = audit[col_order].copy()
audit_path = os.path.join(OUTPUT_DIR, "panelC_marker_audit.tsv")
audit_out.to_csv(audit_path, sep='\t', index=False)
log(f"  [SAVE] {audit_path}")

cat_path = os.path.join(OUTPUT_DIR, "panelC_category_coherence.tsv")
cat_summary.to_csv(cat_path, sep='\t', index=False)
log(f"  [SAVE] {cat_path}")

prop_path = os.path.join(OUTPUT_DIR, "panelC_proposed_panel.txt")
with open(prop_path, 'w') as f:
    f.write("\n".join(prop_lines))
log(f"  [SAVE] {prop_path}")

# =============================================================================
# STEP 7: PATTERN HEATMAP (diagnostic style; not a paper figure)
# =============================================================================
banner("STEP 7: Pattern heatmap")

order = list(leaf_order)
Z_ord = Z[order]
labels = [audit.iloc[j]['gene'] for j in order]
label_syms = [audit.iloc[j]['matched_symbol'] for j in order]
label_cats = [audit.iloc[j]['assigned_category'] for j in order]
label_verd = [audit.iloc[j]['verdict'] for j in order]

# category -> color
uniq_cats = list(dict.fromkeys(assigned.values()))
cmap_cat = plt.get_cmap('tab20', max(20, len(uniq_cats)))
cat_color = {c: cmap_cat(i % 20) for i, c in enumerate(uniq_cats)}

n = len(order)
fig_h = max(6.0, 0.24 * n)
fig, ax = plt.subplots(figsize=(6.4, fig_h))
im = ax.imshow(Z_ord, aspect='auto', cmap='RdBu_r', vmin=-1.5, vmax=1.5)
ax.set_xticks(range(len(POP_ORDER)))
ax.set_xticklabels([POP_LABELS[p] for p in POP_ORDER], rotation=30, ha='right', fontsize=9)
ax.set_yticks(range(n))
ylabels = []
for sym, cat, v in zip(label_syms, label_cats, label_verd):
    suffix = '' if v == 'KEEP' else (' (rev)' if v == 'REVIEW' else ' (DROP)')
    ylabels.append(f"{sym}{suffix}")
ax.set_yticklabels(ylabels, fontsize=7)
for tick, cat, v in zip(ax.get_yticklabels(), label_cats, label_verd):
    tick.set_color(cat_color[cat])
    if v == 'DROP':
        tick.set_fontweight('bold')
ax.set_title("Panel C markers: z-scored expression pattern\n"
             "(row color = assigned category; RdBu = high/low across the 3 populations)",
             fontsize=9)
cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
cbar.set_label("z-score (per gene, across populations)", fontsize=8)
cbar.ax.tick_params(labelsize=7)

legend_handles = [Patch(facecolor=cat_color[c], label=c) for c in uniq_cats]
ax.legend(handles=legend_handles, loc='upper left', bbox_to_anchor=(1.15, 1.0),
          fontsize=6, frameon=False, title='category', title_fontsize=7)

plt.tight_layout()
for ext in ('pdf', 'png'):
    hp = os.path.join(OUTPUT_DIR, f"panelC_pattern_heatmap.{ext}")
    fig.savefig(hp, dpi=DPI, bbox_inches='tight')
log(f"  [SAVE] panelC_pattern_heatmap.pdf/.png")
plt.close(fig)

# =============================================================================
# WRITE REPORT
# =============================================================================
report_path = os.path.join(OUTPUT_DIR, "panelC_marker_audit_report.txt")
with open(report_path, 'w') as f:
    f.write("\n".join(report_lines))
print(f"\n[SAVE] {report_path}")
print("PANEL C MARKER AUDIT COMPLETE")
