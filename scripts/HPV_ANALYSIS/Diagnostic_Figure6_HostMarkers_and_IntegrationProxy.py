#!/usr/bin/env python3
"""
Diagnostic_Figure6_HostMarkers_and_IntegrationProxy.py  (v3 -- text-audit harness)
===================================================================================
Figure 6 diagnostic + Section 4.4 text-number verification.

v3 adds three sections on top of v2, so that every quantitative claim in
Section 4.4 can be checked from one script against the same source files and
the same methodology the figure uses:

  SECTION 1 (NEW): LIFECYCLE FRACTIONS. Mirrors Generate_Figure6_Lifecycle_Panels.py
    Panel F EXACTLY: gated HPV16-positive set (raw_HPV16 >= 8 AND total > 0),
    per-cell gene fractions = gene / total (no pseudocount), permutation test on
    the difference of means (10,000 perms, seed 42), BH-FDR within the 8-gene
    family and separately within the 4-phase family. Reproduces the figure's
    "PANEL F MEAN FRACTIONS" block so the Para-3 numbers are confirmable here.

  SECTION 2 (NEW): READ-CLASS / URR BREAKDOWN. Reports, per group on the gated
    set, the URR / ORF / intergenic read fractions BOTH ways: pooled
    (sum reads / sum total, the estimator behind the prose "two-thirds of reads
    in the URR") and per-cell mean (the estimator in the figure's internal URR
    log). These differ for CNV-HIGH (~63.5% pooled vs ~67.1% per-cell mean);
    the prose should cite pooled.

  SECTION 3 (NEW): TEXT NUMBER AUDIT. Diffs the current Section 4.4 prose
    (hardcoded below from the manuscript draft) against freshly computed values
    and prints MATCH / DIFF / OUT-OF-SCOPE per claim. q-values compared on a
    log10 tolerance to absorb 1-2 sig-fig rounding; means/fractions on relative
    or absolute tolerance. Para-1 numbers (basal enrichment, tier counts, Fisher)
    are marked OUT-OF-SCOPE with their correct source (Phase3 L-method).

v2 carried over unchanged:
  - DIAGNOSTIC A: integration proxy on the gated >=8 set, split pseudocount,
    floor of 10 (NORMAL -> N.D.), BH across the proxy family.
  - VIRAL LOAD SUMMARY: raw_HPV16 (all cells) vs total reads (gated set).
  - DIAGNOSTIC B: host-marker panel, ungated 546/546/546, BH per contrast.

INPUTS (identical to the figure script):
  - data/FIG_4/01_group_selection/three_group_assignments.tsv
  - data/FIG_4/00_input/adata_final.h5ad
  - data/FIG_6/01_raw_hpv16_counts/basal_cell_master_table_with_raw_HPV16.tsv
  - data/FIG_6/03_hpv16_genome/per_cell_hpv16_gene_counts.tsv

OUTPUTS (to data/FIG_6/DIAGNOSTIC_LIFECYCLE_MARKERS/):
  - diagnostic_figure6_report.txt
  - integration_proxy_metrics.tsv
  - viral_load_summary.tsv
  - host_marker_expression_summary.tsv
  - host_marker_per_cell_values.tsv
  - lifecycle_fractions_panelF_mirror.tsv
  - readclass_urr_breakdown.tsv
  - section4_4_text_audit.tsv

Env: NETWORK
Usage: conda run -n NETWORK python Diagnostic_Figure6_HostMarkers_and_IntegrationProxy.py

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from scipy.stats import mannwhitneyu, kruskal
from statsmodels.stats.multitest import multipletests
from collections import OrderedDict
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION  (paths and constants copied from the figure script)
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"

THREE_GROUP_PATH = os.path.join(PROJECT_ROOT,
    "data/FIG_4/01_group_selection/three_group_assignments.tsv")
ADATA_PATH = os.path.join(PROJECT_ROOT,
    "data/FIG_4/00_input/adata_final.h5ad")
MASTER_TABLE_PATH = os.path.join(PROJECT_ROOT,
    "data/FIG_6/01_raw_hpv16_counts/basal_cell_master_table_with_raw_HPV16.tsv")
HPV_GENE_PATH = os.path.join(PROJECT_ROOT,
    "data/FIG_6/03_hpv16_genome/per_cell_hpv16_gene_counts.tsv")

OUTPUT_DIR = os.path.join(PROJECT_ROOT,
    "data/FIG_6/DIAGNOSTIC_LIFECYCLE_MARKERS")
os.makedirs(OUTPUT_DIR, exist_ok=True)

HPV16_THRESHOLD = 8
TOTAL_COL = 'total_hpv16_genome_reads'
MIN_CELLS_FOR_STATS = 10        # figure floor; NORMAL (n=8) -> N.D. on the gated set
N_PERM = 10000
PERM_SEED = 42

POP_ORDER = ['SBS2_HIGH', 'CNV_HIGH', 'NORMAL']
POP_LABELS = {'SBS2_HIGH': 'SBS2-HIGH', 'CNV_HIGH': 'CNV-HIGH', 'NORMAL': 'Normal'}

# Lifecycle phases (figure order)
HPV16_PHASES = OrderedDict([
    ('Maintenance',   ['E1', 'E2']),
    ('Amplification', ['E4', 'E5']),
    ('Oncogene',      ['E6', 'E7']),
    ('Capsid',        ['L1', 'L2']),
])
ALL_HPV_GENES = [g for genes in HPV16_PHASES.values() for g in genes]

# =============================================================================
# HOST MARKER GENE PANEL (53 genes)
# =============================================================================
MARKER_GENES = OrderedDict([
    ('Transformation', ['CDKN2A', 'MCM7', 'CCNE1', 'MKI67', 'TOP2A',
                        'PCNA', 'BRD4', 'MED1', 'E2F1', 'E2F2']),
    ('p53_Rb_pathway', ['CDKN1A', 'BAX', 'MDM2', 'RB1', 'TP53']),
    ('ATM_DNA_damage', ['ATM', 'CHEK2', 'BRCA1', 'NBN', 'MRE11',
                        'RAD50', 'H2AX', 'CHEK1', 'STAT5A', 'STAT5B']),
    ('G2M_arrest',     ['CDC25A', 'CDC25C', 'CDK1', 'CCNB1']),
    ('Caspase',        ['CASP3', 'CASP7']),
    ('Immune',         ['STAT1', 'HLA-A', 'HLA-B', 'HLA-C', 'IRF1', 'TAP1', 'B2M']),
    ('Differentiation', ['KRT5', 'KRT14', 'KRT1', 'KRT10', 'CDH1', 'IVL']),
    ('Innate_APOBEC',  ['APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D',
                        'APOBEC3F', 'APOBEC3G', 'APOBEC3H', 'CGAS', 'STING1']),
])

# =============================================================================
# LOGGING
# =============================================================================
report_lines = []

def log(msg=""):
    print(msg, flush=True)
    report_lines.append(str(msg))

def banner(title, char="="):
    log("")
    log(char * 90)
    log(f"  {title}")
    log(char * 90)

# =============================================================================
# HELPERS  (copied from the figure script so methodology is byte-identical)
# =============================================================================
def get_expression(adata, gene_symbol):
    if gene_symbol in adata.var_names:
        idx = adata.var_names.get_loc(gene_symbol)
        x = adata.X[:, idx]
        if scipy.sparse.issparse(x):
            return np.asarray(x.todense()).flatten()
        return np.asarray(x).flatten()
    if 'gene_symbol' in adata.var.columns:
        mask = adata.var['gene_symbol'] == gene_symbol
        if mask.any():
            idx = np.where(mask)[0][0]
            x = adata.X[:, idx]
            if scipy.sparse.issparse(x):
                return np.asarray(x.todense()).flatten()
            return np.asarray(x).flatten()
    return None

def mw_with_floor(v1, v2):
    v1 = np.asarray(v1, dtype=float)
    v2 = np.asarray(v2, dtype=float)
    if (len(v1) >= MIN_CELLS_FOR_STATS and len(v2) >= MIN_CELLS_FOR_STATS
            and len(v1) > 5 and len(v2) > 5):
        _, p = mannwhitneyu(v1, v2, alternative='two-sided')
        return p
    return np.nan

def permutation_test_means(v1, v2, n_perm=N_PERM, seed=PERM_SEED):
    """Two-sided permutation test on |difference of means|. Figure-identical."""
    v1 = np.asarray(v1, dtype=float)
    v2 = np.asarray(v2, dtype=float)
    if len(v1) < MIN_CELLS_FOR_STATS or len(v2) < MIN_CELLS_FOR_STATS:
        return np.nan
    obs = abs(v1.mean() - v2.mean())
    pooled = np.concatenate([v1, v2])
    n1 = len(v1)
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(pooled)
        if abs(perm[:n1].mean() - perm[n1:].mean()) >= obs:
            count += 1
    return (count + 1) / (n_perm + 1)

def compute_pairwise_perm(data_dict):
    pairs = [(0, 1), (1, 2), (0, 2)]
    out = []
    for i, j in pairs:
        v1 = data_dict[POP_ORDER[i]]
        v2 = data_dict[POP_ORDER[j]]
        if len(v1) >= MIN_CELLS_FOR_STATS and len(v2) >= MIN_CELLS_FOR_STATS:
            out.append(permutation_test_means(v1, v2))
        else:
            out.append(np.nan)
    return out

def bh(raw_list):
    pvals = np.array(raw_list, dtype=float)
    out = np.full_like(pvals, np.nan)
    valid = ~np.isnan(pvals)
    if valid.sum() == 0:
        return out.tolist()
    _, adj, _, _ = multipletests(pvals[valid], method='fdr_bh')
    out[valid] = adj
    return out.tolist()

def stars(q):
    if q is None or (isinstance(q, float) and np.isnan(q)):
        return 'N.D.'
    if q < 1e-4: return '****'
    if q < 1e-3: return '***'
    if q < 0.01: return '**'
    if q < 0.05: return '*'
    return 'ns'

def fmt_p(p):
    return 'N.D.' if (p is None or (isinstance(p, float) and np.isnan(p))) else f"{p:.2e}"

# =============================================================================
# STEP 0: LOAD DATA  (mirrors the figure script's STEP 0)
# =============================================================================
banner("STEP 0: Load data (mirroring figure script)")

groups = pd.read_csv(THREE_GROUP_PATH, sep='\t')
sbs2_cells   = set(groups.loc[groups['group'] == 'SBS2_HIGH', 'cell_barcode'])
cnv_cells    = set(groups.loc[groups['group'] == 'CNV_HIGH',  'cell_barcode'])
normal_cells = set(groups.loc[groups['group'] == 'NORMAL',   'cell_barcode'])
cell_to_group = dict(zip(groups['cell_barcode'], groups['group']))
log(f"  Populations: {len(sbs2_cells)} SBS2-HIGH, {len(cnv_cells)} CNV-HIGH, "
    f"{len(normal_cells)} Normal")

log("  Loading adata_final.h5ad ...")
adata = sc.read_h5ad(ADATA_PATH)
log(f"  adata: {adata.shape[0]} cells x {adata.shape[1]} genes")
adata.obs['population'] = 'other'
adata.obs.loc[adata.obs_names.isin(sbs2_cells), 'population'] = 'SBS2_HIGH'
adata.obs.loc[adata.obs_names.isin(cnv_cells), 'population'] = 'CNV_HIGH'
adata.obs.loc[adata.obs_names.isin(normal_cells), 'population'] = 'NORMAL'
adata_pop = adata[adata.obs['population'].isin(POP_ORDER)].copy()
log(f"  Cells in three populations (host-marker set, ungated): {adata_pop.shape[0]}")

master = pd.read_csv(MASTER_TABLE_PATH, sep='\t', index_col=0)
master['group'] = master.index.map(lambda x: cell_to_group.get(x, 'other'))
master_pop = master[master['group'].isin(POP_ORDER)].copy()
log(f"  Master table rows (all basal cells): {len(master)}")
log(f"  Master table cells in three populations: {len(master_pop)}")

hpv_genes = pd.read_csv(HPV_GENE_PATH, sep='\t', index_col=0)
hpv_genes['population'] = 'other'
hpv_genes.loc[hpv_genes.index.isin(sbs2_cells), 'population'] = 'SBS2_HIGH'
hpv_genes.loc[hpv_genes.index.isin(cnv_cells), 'population'] = 'CNV_HIGH'
hpv_genes.loc[hpv_genes.index.isin(normal_cells), 'population'] = 'NORMAL'
log(f"  HPV16 gene counts: {hpv_genes.shape}")

# Gate to the figure's Panel F positive set
hpv_pos = hpv_genes[hpv_genes['population'].isin(POP_ORDER)].copy()
hpv_pos['raw_HPV16'] = hpv_pos.index.map(master['raw_HPV16'])
hpv_pos = hpv_pos[(hpv_pos['raw_HPV16'] >= HPV16_THRESHOLD) &
                  (hpv_pos[TOTAL_COL] > 0)].copy()
for col in ALL_HPV_GENES + ['URR', 'intergenic']:
    if col not in hpv_pos.columns:
        hpv_pos[col] = 0.0
    hpv_pos[col] = hpv_pos[col].fillna(0.0)
gated_counts = {p: int((hpv_pos['population'] == p).sum()) for p in POP_ORDER}
log(f"\n  Gated HPV16-positive set (raw_HPV16 >= {HPV16_THRESHOLD} AND {TOTAL_COL} > 0):")
for p in POP_ORDER:
    log(f"    {POP_LABELS[p]}: {gated_counts[p]}")
log(f"  >>> MUST match Panel F (expected 197 / 446 / 8).")

# Store computed values for the audit at the end
AUDIT = {}   # key -> computed value
AUDIT['F_count_SBS2'] = gated_counts['SBS2_HIGH']
AUDIT['F_count_CNV']  = gated_counts['CNV_HIGH']
AUDIT['F_count_NORM'] = gated_counts['NORMAL']


# =============================================================================
# SECTION 1: LIFECYCLE FRACTIONS  (mirror of figure Panel F; Para-3 source)
# =============================================================================
banner("SECTION 1: Lifecycle fractions (Panel F mirror; per-cell frac = gene/total)")

# Per-cell gene fractions (bare total, exactly like the figure)
for g in ALL_HPV_GENES:
    hpv_pos[f'{g}_frac'] = hpv_pos[g] / hpv_pos[TOTAL_COL]

gene_frac = {}
for g in ALL_HPV_GENES:
    gene_frac[g] = {p: hpv_pos.loc[hpv_pos['population'] == p, f'{g}_frac'].values.astype(float)
                    for p in POP_ORDER}

# Per-gene permutation, BH within the 8-gene family
gene_raw = []
for g in ALL_HPV_GENES:
    gene_raw.extend(compute_pairwise_perm(gene_frac[g]))
gene_q = bh(gene_raw)
gene_qvals = {g: gene_q[i*3:(i+1)*3] for i, g in enumerate(ALL_HPV_GENES)}

# Per-phase fractions, permutation, BH within the 4-phase family (separate)
phase_frac = OrderedDict()
for phase, genes in HPV16_PHASES.items():
    hpv_pos[f'{phase}_frac'] = hpv_pos[[f'{g}_frac' for g in genes]].sum(axis=1)
    phase_frac[phase] = {p: hpv_pos.loc[hpv_pos['population'] == p, f'{phase}_frac'].values.astype(float)
                         for p in POP_ORDER}
phase_raw = []
for phase in HPV16_PHASES:
    phase_raw.extend(compute_pairwise_perm(phase_frac[phase]))
phase_q = bh(phase_raw)
phase_qvals = {ph: phase_q[i*3:(i+1)*3] for i, ph in enumerate(HPV16_PHASES)}

# Print identical to the figure's "PANEL F MEAN FRACTIONS" block
log(f"\n  {'Item':<14s}  {'SBS2-HIGH':>10s}  {'CNV-HIGH':>10s}  {'Normal':>10s}  {'q SBS2vCNV':>13s}")
log(f"  {'-'*14}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*13}")
lc_rows = []
for g in ALL_HPV_GENES:
    m = {p: (100*np.mean(gene_frac[g][p]) if len(gene_frac[g][p]) else np.nan) for p in POP_ORDER}
    q01 = gene_qvals[g][0]
    log(f"  {g:<14s}  {m['SBS2_HIGH']:>9.2f}%  {m['CNV_HIGH']:>9.2f}%  {m['NORMAL']:>9.2f}%  {q01:>13.4e}")
    lc_rows.append({'item': g, 'kind': 'gene',
                    'SBS2_pct': m['SBS2_HIGH'], 'CNV_pct': m['CNV_HIGH'], 'NORM_pct': m['NORMAL'],
                    'q_SBS2vCNV': q01})
log(f"  {'-'*14}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*13}")
for phase in HPV16_PHASES:
    m = {p: (100*np.mean(phase_frac[phase][p]) if len(phase_frac[phase][p]) else np.nan) for p in POP_ORDER}
    q01 = phase_qvals[phase][0]
    log(f"  {phase:<14s}  {m['SBS2_HIGH']:>9.2f}%  {m['CNV_HIGH']:>9.2f}%  {m['NORMAL']:>9.2f}%  {q01:>13.4e}")
    lc_rows.append({'item': phase, 'kind': 'phase',
                    'SBS2_pct': m['SBS2_HIGH'], 'CNV_pct': m['CNV_HIGH'], 'NORM_pct': m['NORMAL'],
                    'q_SBS2vCNV': q01})
pd.DataFrame(lc_rows).to_csv(os.path.join(OUTPUT_DIR, "lifecycle_fractions_panelF_mirror.tsv"),
                             sep='\t', index=False)

# Stash for audit (SBS2-vs-CNV q on pair index 0)
AUDIT['q_E1']   = gene_qvals['E1'][0]
AUDIT['q_E2']   = gene_qvals['E2'][0]
AUDIT['q_E5']   = gene_qvals['E5'][0]
AUDIT['q_L1']   = gene_qvals['L1'][0]
AUDIT['q_L2']   = gene_qvals['L2'][0]
AUDIT['q_Oncogene'] = phase_qvals['Oncogene'][0]
AUDIT['Oncogene_SBS2_pct'] = 100*np.mean(phase_frac['Oncogene']['SBS2_HIGH'])
AUDIT['Oncogene_CNV_pct']  = 100*np.mean(phase_frac['Oncogene']['CNV_HIGH'])


# =============================================================================
# SECTION 2: READ-CLASS / URR BREAKDOWN  (pooled vs per-cell mean)
# =============================================================================
banner("SECTION 2: Read-class / URR breakdown (pooled vs per-cell mean)")

hpv_pos['ORF_sum'] = hpv_pos[ALL_HPV_GENES].sum(axis=1)

log(f"\n  {'Group':<12s} {'n':>5s}  {'URR pooled':>11s} {'URR percell':>12s}  "
    f"{'ORF pooled':>11s} {'intergenic pooled':>18s}")
log(f"  {'-'*12} {'-'*5}  {'-'*11} {'-'*12}  {'-'*11} {'-'*18}")
rc_rows = []
for p in POP_ORDER:
    sub = hpv_pos[hpv_pos['population'] == p]
    n = len(sub)
    if n == 0:
        continue
    pooled_urr   = 100 * sub['URR'].sum() / sub[TOTAL_COL].sum()
    pooled_orf   = 100 * sub['ORF_sum'].sum() / sub[TOTAL_COL].sum()
    pooled_int   = 100 * sub['intergenic'].sum() / sub[TOTAL_COL].sum()
    percell_urr  = 100 * (sub['URR'] / sub[TOTAL_COL]).mean()
    log(f"  {POP_LABELS[p]:<12s} {n:>5d}  {pooled_urr:>10.1f}% {percell_urr:>11.1f}%  "
        f"{pooled_orf:>10.1f}% {pooled_int:>17.1f}%")
    rc_rows.append({'group': POP_LABELS[p], 'n': n,
                    'URR_pooled_pct': pooled_urr, 'URR_percell_mean_pct': percell_urr,
                    'ORF_pooled_pct': pooled_orf, 'intergenic_pooled_pct': pooled_int})
    if p == 'SBS2_HIGH': AUDIT['URR_pooled_SBS2'] = pooled_urr
    if p == 'CNV_HIGH':  AUDIT['URR_pooled_CNV']  = pooled_urr
    if p == 'NORMAL':    AUDIT['URR_pooled_NORM'] = pooled_urr
pd.DataFrame(rc_rows).to_csv(os.path.join(OUTPUT_DIR, "readclass_urr_breakdown.tsv"),
                             sep='\t', index=False)
log("\n  NOTE: prose 'two-thirds of reads in the URR' should cite POOLED URR.")
log("  The figure's internal URR log uses the per-cell mean (differs for CNV-HIGH).")


# =============================================================================
# DIAGNOSTIC A: INTEGRATION PROXY  (gated >=8 set, figure-matched)  [v2]
# =============================================================================
banner("DIAGNOSTIC A: Integration proxy (gated HPV16+ set, n = 197/446/8)")

hpv_pos['E6E7_sum']    = hpv_pos['E6'] + hpv_pos['E7']
hpv_pos['total_early'] = hpv_pos[['E1', 'E2', 'E4', 'E5', 'E6', 'E7']].sum(axis=1)
hpv_pos['E2_to_E6E7']          = hpv_pos['E2'] / (hpv_pos['E6E7_sum'] + 0.5)
hpv_pos['E2_fraction_of_early'] = hpv_pos['E2'] / (hpv_pos['total_early'] + 0.5)
hpv_pos['E6E7_frac_of_total'] = hpv_pos['E6E7_sum'] / hpv_pos[TOTAL_COL]
hpv_pos['E2_frac_of_total']   = hpv_pos['E2'] / hpv_pos[TOTAL_COL]
hpv_pos['L1L2_frac_of_total'] = (hpv_pos['L1'] + hpv_pos['L2']) / hpv_pos[TOTAL_COL]

PROXY_METRICS = ['E2_to_E6E7', 'E2_fraction_of_early', 'E6E7_sum', 'E2', TOTAL_COL,
                 'E6E7_frac_of_total', 'E2_frac_of_total', 'L1L2_frac_of_total']
if 'early_late_ratio' in hpv_pos.columns:
    PROXY_METRICS.append('early_late_ratio')
PSEUDOCOUNT_FLAG = {'E2_to_E6E7': '(+0.5)', 'E2_fraction_of_early': '(+0.5)'}

proxy_rows = []
hvc_raw = []
for metric in PROXY_METRICS:
    vals = {p: hpv_pos.loc[hpv_pos['population'] == p, metric].dropna().values for p in POP_ORDER}
    means = {p: (np.mean(vals[p]) if len(vals[p]) else np.nan) for p in POP_ORDER}
    p_hvc = mw_with_floor(vals['SBS2_HIGH'], vals['CNV_HIGH'])
    hvc_raw.append(p_hvc)
    direction = ('SBS2 > CNV' if means['SBS2_HIGH'] > means['CNV_HIGH'] else 'CNV > SBS2')
    proxy_rows.append({'metric': metric, **{f'mean_{p}': means[p] for p in POP_ORDER},
                       'hvc_raw_p': p_hvc, 'direction': direction})
hvc_q = bh(hvc_raw)
for row, q in zip(proxy_rows, hvc_q):
    row['hvc_bh_q'] = q

log(f"\n  {'Metric':<24s} {'SBS2-HIGH':>12s} {'CNV-HIGH':>12s} {'NORMAL*':>12s}  "
    f"{'HvC raw p':>11s} {'HvC BH q':>11s}  {'Direction':>12s}")
log(f"  {'-'*24} {'-'*12} {'-'*12} {'-'*12}  {'-'*11} {'-'*11}  {'-'*12}")
for row in proxy_rows:
    flag = PSEUDOCOUNT_FLAG.get(row['metric'], '')
    name = f"{row['metric']}{(' ' + flag) if flag else ''}"
    log(f"  {name:<24s} {row['mean_SBS2_HIGH']:>12.4f} {row['mean_CNV_HIGH']:>12.4f} "
        f"{row['mean_NORMAL']:>12.4f}  {fmt_p(row['hvc_raw_p']):>11s} "
        f"{fmt_p(row['hvc_bh_q']):>11s}  {row['direction']:>12s} {stars(row['hvc_bh_q'])}")
log("  * NORMAL (n=8) below the 10-cell floor; descriptive only.")
pd.DataFrame(proxy_rows).to_csv(os.path.join(OUTPUT_DIR, "integration_proxy_metrics.tsv"),
                                sep='\t', index=False)


# =============================================================================
# VIRAL LOAD SUMMARY  [v2]
# =============================================================================
banner("VIRAL LOAD SUMMARY (pick one; label its cell set and measure in text)")

load_a = {p: np.mean(master_pop.loc[master_pop['group'] == p, 'raw_HPV16'].values.astype(float))
          for p in POP_ORDER}
fold_a = load_a['CNV_HIGH'] / load_a['SBS2_HIGH'] if load_a['SBS2_HIGH'] > 0 else np.nan
load_b = {p: np.mean(hpv_pos.loc[hpv_pos['population'] == p, TOTAL_COL].values.astype(float))
          if (hpv_pos['population'] == p).sum() else np.nan for p in POP_ORDER}
fold_b = load_b['CNV_HIGH'] / load_b['SBS2_HIGH']
load_b_q    = next(r['hvc_bh_q'] for r in proxy_rows if r['metric'] == TOTAL_COL)
load_b_rawp = next(r['hvc_raw_p'] for r in proxy_rows if r['metric'] == TOTAL_COL)

log(f"\n  (a) raw_HPV16 UMI, ALL cells per group (n=546)  [Panel D measure]")
log(f"      SBS2 {load_a['SBS2_HIGH']:.1f} | CNV {load_a['CNV_HIGH']:.1f} | "
    f"NORM {load_a['NORMAL']:.1f}  -> CNV/SBS2 = {fold_a:.2f}x   (Panel D q: 2.32e-73)")
log(f"  (b) {TOTAL_COL}, gated HPV16+ set (n=197/446/8)  [Panel F cell set]")
log(f"      SBS2 {load_b['SBS2_HIGH']:.1f} | CNV {load_b['CNV_HIGH']:.1f} | "
    f"NORM {load_b['NORMAL']:.1f} (desc)  -> CNV/SBS2 = {fold_b:.2f}x   "
    f"HvC raw p = {fmt_p(load_b_rawp)}, BH q = {fmt_p(load_b_q)}")
log(f"  RECOMMENDATION: cite (b) for a 'per HPV16-positive cell' load sentence.")
pd.DataFrame([
    {'measure': 'raw_HPV16_all_cells', 'cell_set': 'all_546',
     **{f'mean_{p}': load_a[p] for p in POP_ORDER}, 'fold_CNV_over_SBS2': fold_a},
    {'measure': TOTAL_COL, 'cell_set': 'gated_pos_197_446_8',
     **{f'mean_{p}': load_b[p] for p in POP_ORDER}, 'fold_CNV_over_SBS2': fold_b},
]).to_csv(os.path.join(OUTPUT_DIR, "viral_load_summary.tsv"), sep='\t', index=False)

AUDIT['load_SBS2'] = load_b['SBS2_HIGH']
AUDIT['load_CNV']  = load_b['CNV_HIGH']
AUDIT['load_fold'] = fold_b
AUDIT['load_q']    = load_b_q


# =============================================================================
# DIAGNOSTIC B: HOST MARKER PANEL  (ungated, all 1,638; BH per contrast)  [v2]
# =============================================================================
banner("DIAGNOSTIC B: Host marker panel (ungated, 546/546/546; BH per contrast)")

flat_genes = [(cat, g) for cat, genes in MARKER_GENES.items() for g in genes]
records = OrderedDict()
raw_hvc, raw_hvn, raw_cvn, order_genes = [], [], [], []
per_cell_rows = []
missing = []

for cat, gene in flat_genes:
    expr = get_expression(adata_pop, gene)
    if expr is None:
        missing.append(gene)
        continue
    vals = {p: expr[(adata_pop.obs['population'] == p).values] for p in POP_ORDER}
    means = {p: float(np.mean(vals[p])) for p in POP_ORDER}
    pcts  = {p: 100.0 * np.sum(vals[p] > 0) / max(len(vals[p]), 1) for p in POP_ORDER}
    try:
        _, kw_p = kruskal(vals['SBS2_HIGH'], vals['CNV_HIGH'], vals['NORMAL'])
    except Exception:
        kw_p = np.nan
    p_hvc = mw_with_floor(vals['SBS2_HIGH'], vals['CNV_HIGH'])
    p_hvn = mw_with_floor(vals['SBS2_HIGH'], vals['NORMAL'])
    p_cvn = mw_with_floor(vals['CNV_HIGH'],  vals['NORMAL'])
    records[gene] = {'category': cat, 'gene': gene,
                     **{f'mean_{p}': means[p] for p in POP_ORDER},
                     **{f'pct_{p}': pcts[p] for p in POP_ORDER},
                     'kw_p': kw_p, 'hvc_raw_p': p_hvc, 'hvn_raw_p': p_hvn, 'cvn_raw_p': p_cvn,
                     'hvc_dir': 'SBS2 > CNV' if means['SBS2_HIGH'] > means['CNV_HIGH'] else 'CNV > SBS2'}
    raw_hvc.append(p_hvc); raw_hvn.append(p_hvn); raw_cvn.append(p_cvn)
    order_genes.append(gene)
    for p in POP_ORDER:
        m = (adata_pop.obs['population'] == p).values
        for bc, v in zip(adata_pop.obs_names[m], vals[p]):
            per_cell_rows.append({'cell_barcode': bc, 'population': p,
                                  'gene': gene, 'category': cat, 'expression': float(v)})

q_hvc = bh(raw_hvc); q_hvn = bh(raw_hvn); q_cvn = bh(raw_cvn)
for gene, qh, qn, qc in zip(order_genes, q_hvc, q_hvn, q_cvn):
    records[gene]['hvc_bh_q'] = qh
    records[gene]['hvn_bh_q'] = qn
    records[gene]['cvn_bh_q'] = qc
if missing:
    log(f"  WARNING: {len(missing)} marker gene(s) not found: {missing}")

log(f"\n  {'Category':<16s} {'Gene':<10s} {'SBS2':>8s} {'CNV':>8s} {'NORM':>8s}  "
    f"{'KW p':>9s} {'HvC q':>9s} {'Dir':>11s}")
log(f"  {'-'*16} {'-'*10} {'-'*8} {'-'*8} {'-'*8}  {'-'*9} {'-'*9} {'-'*11}")
for cat, genes in MARKER_GENES.items():
    for gene in genes:
        if gene not in records:
            continue
        r = records[gene]
        log(f"  {cat:<16s} {gene:<10s} {r['mean_SBS2_HIGH']:>8.3f} {r['mean_CNV_HIGH']:>8.3f} "
            f"{r['mean_NORMAL']:>8.3f}  {fmt_p(r['kw_p']):>9s} {fmt_p(r['hvc_bh_q']):>9s} "
            f"{r['hvc_dir']:>11s} {stars(r['hvc_bh_q'])}")
pd.DataFrame([records[g] for g in order_genes]).to_csv(
    os.path.join(OUTPUT_DIR, "host_marker_expression_summary.tsv"), sep='\t', index=False)
pd.DataFrame(per_cell_rows).to_csv(
    os.path.join(OUTPUT_DIR, "host_marker_per_cell_values.tsv"), sep='\t', index=False)

log(f"\n  Panel B cross-check (must equal figure Panel B):")
for gene in ['APOBEC3A', 'APOBEC3B']:
    r = records[gene]
    log(f"    {gene}: SBS2 {r['mean_SBS2_HIGH']:.4f}  CNV {r['mean_CNV_HIGH']:.4f}  "
        f"NORM {r['mean_NORMAL']:.4f}")

# Stash host-marker values for the audit
def stash(gene, key):
    if gene in records:
        AUDIT[f'{key}_SBS2'] = records[gene]['mean_SBS2_HIGH']
        AUDIT[f'{key}_CNV']  = records[gene]['mean_CNV_HIGH']
        AUDIT[f'q_{key}']    = records[gene]['hvc_bh_q']
for gene, key in [('B2M','B2M'),('HLA-A','HLAA'),('HLA-B','HLAB'),('TAP1','TAP1'),
                  ('IVL','IVL'),('APOBEC3A','A3A'),('APOBEC3B','A3B'),
                  ('BRCA1','BRCA1'),('H2AX','H2AX'),('CHEK2','CHEK2'),
                  ('CDK1','CDK1'),('CCNB1','CCNB1'),('CDC25C','CDC25C'),
                  ('KRT14','KRT14'),('MKI67','MKI67'),('KRT5','KRT5')]:
    stash(gene, key)


# =============================================================================
# SECTION 3: TEXT NUMBER AUDIT  (diff current Section 4.4 prose vs computed)
# =============================================================================
banner("SECTION 3: Section 4.4 text-number audit")

# Claims hardcoded from the current manuscript draft of Section 4.4.
# kind: 'q' (log10 tol), 'mean' (rel tol), 'pct' (abs tol), 'fold' (abs tol),
#       'count' (exact), 'scope' (out-of-scope, source noted)
CLAIMS = [
    # (label, claimed, computed_key, kind)
    ('Gated count SBS2 = 197',        197,    'F_count_SBS2', 'count'),
    ('Gated count CNV = 446',         446,    'F_count_CNV',  'count'),
    ('Gated count NORMAL = 8',        8,      'F_count_NORM', 'count'),
    ('URR SBS2 63.5% (pooled)',       63.5,   'URR_pooled_SBS2', 'pct'),
    ('URR CNV 63.5% (pooled)',        63.5,   'URR_pooled_CNV',  'pct'),
    ('URR NORMAL 64.9% (pooled)',     64.9,   'URR_pooled_NORM', 'pct'),
    ('Load SBS2 90.1',                90.1,   'load_SBS2',    'mean'),
    ('Load CNV 235.1',                235.1,  'load_CNV',     'mean'),
    ('Load fold 2.6x',                2.6,    'load_fold',    'fold'),
    ('Load q 1.7e-14',                1.7e-14,'load_q',       'q'),
    ('E1 q 2.0e-4',                   2.0e-4, 'q_E1',         'q'),
    ('L1 q 2.0e-4',                   2.0e-4, 'q_L1',         'q'),
    ('L2 q 2.0e-4',                   2.0e-4, 'q_L2',         'q'),
    ('E5 q 2.0e-4',                   2.0e-4, 'q_E5',         'q'),
    ('Oncogene q 0.10',               0.10,   'q_Oncogene',   'q'),
    ('Oncogene SBS2 <1% (0.53)',      0.53,   'Oncogene_SBS2_pct', 'pct'),
    ('Oncogene CNV <1% (0.79)',       0.79,   'Oncogene_CNV_pct',  'pct'),
    ('E2 q 0.11',                     0.11,   'q_E2',         'q'),
    ('B2M q 4.0e-100',                4.0e-100,'q_B2M',       'q'),
    ('HLA-A q 5.8e-42',               5.8e-42,'q_HLAA',       'q'),
    ('HLA-B q 1.1e-15',               1.1e-15,'q_HLAB',       'q'),
    ('TAP1 q 1.2e-4',                 1.2e-4, 'q_TAP1',       'q'),
    ('IVL SBS2 2.68',                 2.68,   'IVL_SBS2',     'mean'),
    ('IVL CNV 0.09',                  0.09,   'IVL_CNV',      'mean'),
    ('IVL q 3.1e-70',                 3.1e-70,'q_IVL',        'q'),
    ('A3A SBS2 6.46',                 6.46,   'A3A_SBS2',     'mean'),
    ('A3A CNV 2.08',                  2.08,   'A3A_CNV',      'mean'),
    ('A3A q 5.5e-134',                5.5e-134,'q_A3A',       'q'),
    ('A3B SBS2 2.21',                 2.21,   'A3B_SBS2',     'mean'),
    ('A3B CNV 4.95',                  4.95,   'A3B_CNV',      'mean'),
    ('A3B q 8.3e-66',                 8.3e-66,'q_A3B',        'q'),
    ('BRCA1 q 1.2e-8',                1.2e-8, 'q_BRCA1',      'q'),
    ('H2AX q 4.6e-11',                4.6e-11,'q_H2AX',       'q'),
    ('CHEK2 q 1.6e-5',                1.6e-5, 'q_CHEK2',      'q'),
    ('CDK1 q 4.9e-18',                4.9e-18,'q_CDK1',       'q'),
    ('CCNB1 q 9.9e-27',               9.9e-27,'q_CCNB1',      'q'),
    ('CDC25C q 3.3e-12',              3.3e-12,'q_CDC25C',     'q'),
    ('KRT14 q 4.8e-57',               4.8e-57,'q_KRT14',      'q'),
    ('MKI67 q 2.6e-19',               2.6e-19,'q_MKI67',      'q'),
]

def verdict(claimed, computed, kind):
    if computed is None or (isinstance(computed, float) and np.isnan(computed)):
        return 'NO VALUE'
    if kind == 'q':
        if computed <= 0:
            return 'DIFF'
        return 'MATCH' if abs(np.log10(computed) - np.log10(claimed)) < 0.06 else 'DIFF'
    if kind == 'mean':
        return 'MATCH' if abs(computed - claimed) <= max(0.03, 0.02*abs(claimed)) else 'DIFF'
    if kind == 'pct':
        return 'MATCH' if abs(computed - claimed) <= 0.2 else 'DIFF'
    if kind == 'fold':
        return 'MATCH' if abs(computed - claimed) <= 0.1 else 'DIFF'
    if kind == 'count':
        return 'MATCH' if int(round(computed)) == int(claimed) else 'DIFF'
    return '?'

log(f"\n  {'Claim':<28s} {'claimed':>12s} {'computed':>14s}   Verdict")
log(f"  {'-'*28} {'-'*12} {'-'*14}   -------")
audit_rows = []
n_match = n_diff = n_novalue = 0
for label, claimed, key, kind in CLAIMS:
    computed = AUDIT.get(key)
    v = verdict(claimed, computed, kind)
    if v == 'MATCH': n_match += 1
    elif v == 'NO VALUE': n_novalue += 1
    else: n_diff += 1
    comp_str = ('--' if computed is None else
                (f"{computed:.3g}" if kind in ('q',) else f"{computed:.4g}"))
    cl_str = f"{claimed:.3g}" if kind == 'q' else f"{claimed:g}"
    log(f"  {label:<28s} {cl_str:>12s} {comp_str:>14s}   {v}")
    audit_rows.append({'claim': label, 'claimed': claimed, 'computed': computed, 'verdict': v})

log(f"\n  MATCH: {n_match}   DIFF: {n_diff}   NO VALUE: {n_novalue}")
if n_diff or n_novalue:
    log("  >>> Inspect any DIFF / NO VALUE rows before the text is finalized.")

# Out-of-scope Para 1 numbers: source is the Phase3 L-method population step.
log(f"\n  OUT OF SCOPE for this diagnostic (source = Phase3 L-method / population step):")
log(f"    - 94.6% of HPV16+ cells are basal  (needs non-basal HPV+ counts; this")
log(f"      master table is basal-only, cannot be reconstructed here)")
log(f"    - Tier counts 22,153 / 14,046 / 15,927 (needs the ambiguous-band thresholds)")
log(f"    - Fisher OR = 1.01, p = 0.91 (needs the positivity-vs-SBS2-HIGH contrast set)")
# Partial anchor: positive count at threshold 8 over ALL basal master rows
n_pos_allbasal = int((master['raw_HPV16'] >= HPV16_THRESHOLD).sum())
log(f"    Partial anchor: raw_HPV16 >= {HPV16_THRESHOLD} over all {len(master)} basal "
    f"cells = {n_pos_allbasal}  (compare to the tier 'positive' = 15,927)")

pd.DataFrame(audit_rows).to_csv(os.path.join(OUTPUT_DIR, "section4_4_text_audit.tsv"),
                                sep='\t', index=False)


# =============================================================================
# SAVE REPORT
# =============================================================================
banner("COMPLETE")
report_path = os.path.join(OUTPUT_DIR, "diagnostic_figure6_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Report saved: {report_path}")
log(f"  Output directory: {OUTPUT_DIR}")
