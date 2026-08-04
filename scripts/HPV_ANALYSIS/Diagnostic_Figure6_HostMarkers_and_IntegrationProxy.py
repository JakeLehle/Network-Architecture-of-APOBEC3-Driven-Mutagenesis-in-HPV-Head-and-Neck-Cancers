#!/usr/bin/env python3
"""
Diagnostic_Figure6_HostMarkers_and_IntegrationProxy.py  (v4 -- 57-gene panel)
===================================================================================
Figure 6 diagnostic + Section 4.4 text-number verification.

v4 changes (Panel C alignment; Sections 1/2 and Diagnostic A unchanged):
  - MARKER_GENES replaced with the locked 57-gene / 7-tier Panel C structure,
    mirroring DOTPLOT_CATEGORIES in Generate_Figure6_Lifecycle_Panels.py v6.2.
    The BH family here now equals the panel drawn in Figure 6c, so every host
    q-value quoted in Section 4.4 is corrected across exactly the genes shown.
  - Gene aliases are resolved at load time by renaming adata.var_names. DDX58 is
    stored as 'RIGI' in this transcriptome (2024-A / GENCODE v44); without this
    it drops silently and the family becomes 56, shifting every q.
  - APOBEC3A / APOBEC3B are REMOVED from the BH family. They are Panel B genes,
    not Panel C genes, and their q-values come from the figure script's Panel B
    family of 18. Their means are still computed here, outside the family, as
    the cross-check that this script and the figure script agree.
  - Genes evaluated and dropped from Panel C (CASP3, KRT1, CGAS, STING1, the
    remaining A3 family, IFITM1, BST2, SMC5/6, NSMCE2) are computed outside the
    family as an audit trail for Figure6_PanelC_Tier_Reference.md, without
    contaminating the family.
  - CLAIMS expanded from 39 to a full lock: all 57 panel q-values, including the
    negative (ns) results. Auditing the ns genes is deliberate; the tier
    reference doc previously recorded BARD1 as strong at q=1e-41, which was the
    three-group Kruskal-Wallis p and not the SBS2-vs-CNV contrast (ns, q=0.076).
    Nothing in the old harness could catch that because ns genes were unaudited.

v3 sections retained:
  SECTION 1: LIFECYCLE FRACTIONS. Mirrors Generate_Figure6_Lifecycle_Panels.py
    Panel F EXACTLY: gated HPV16-positive set (raw_HPV16 >= 8 AND total > 0),
    per-cell gene fractions = gene / total (no pseudocount), permutation test on
    the difference of means (10,000 perms, seed 42), BH-FDR within the 8-gene
    family and separately within the 4-phase family.

  SECTION 2: READ-CLASS / URR BREAKDOWN. Per group on the gated set, the
    URR / ORF / intergenic read fractions BOTH ways: pooled (sum reads / sum
    total, the estimator behind the prose "two-thirds of reads in the URR") and
    per-cell mean (the estimator in the figure's internal URR log). These differ
    for CNV-HIGH (~63.5% pooled vs ~67.1% per-cell mean); the prose cites pooled.

  SECTION 3: TEXT NUMBER AUDIT. Diffs the Section 4.4 prose (hardcoded below)
    against freshly computed values and prints MATCH / DIFF / NO VALUE per claim.
    q-values compared on a log10 tolerance to absorb 1-2 sig-fig rounding;
    means/fractions on relative or absolute tolerance. Para-1 numbers and the
    Panel B q-values are marked OUT-OF-SCOPE with their correct source.

  DIAGNOSTIC A: integration proxy on the gated >=8 set, split pseudocount,
    floor of 10 (NORMAL -> N.D.), BH across the proxy family.
  VIRAL LOAD SUMMARY: raw_HPV16 (all cells) vs total reads (gated set).
  DIAGNOSTIC B: host-marker panel, ungated 546/546/546, BH per contrast.

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
  - host_marker_outside_family.tsv
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
# HOST MARKER GENE PANEL (57 genes, 7 tiers)
# Mirrors DOTPLOT_CATEGORIES in Generate_Figure6_Lifecycle_Panels.py v6.2
# exactly, so the BH family here equals the panel drawn in Figure 6c.
# Locked structure: Figure6_PanelC_Tier_Reference.md
# =============================================================================
MARKER_GENES = OrderedDict([
    ('MHCI_AgPres_IFN', ['HLA-A', 'HLA-B', 'HLA-C', 'B2M', 'TAP1',
                         'STAT1', 'IRF1', 'STAT2', 'DDX58']),
    ('IFN_effectors',   ['IFI27', 'ISG15', 'IRF9', 'MX1',
                         'OAS1', 'RSAD2', 'IFI44L', 'IFIT1']),
    ('Differentiation', ['KRT5', 'KRT14', 'IVL', 'KRT10', 'CDH1']),
    ('DDR_ATM_ATR',     ['CHEK2', 'BRCA1', 'NBN', 'H2AX',
                         'BARD1', 'TP53BP1', 'RIF1',
                         'ATM', 'MRE11', 'RAD50',
                         'TOPBP1', 'CHEK1', 'STAT5A', 'STAT5B',
                         'CASP7', 'NSD2']),
    ('CellCycle_Prolif', ['MKI67', 'TOP2A', 'MCM7', 'PCNA', 'CCNE1',
                          'CDKN2A', 'E2F1', 'E2F2', 'BRD4', 'MED1']),
    ('p53_Rb_pathway',  ['CDKN1A', 'MDM2', 'BAX', 'TP53', 'RB1']),
    ('G2M_arrest',      ['CDC25A', 'CDC25C', 'CDK1', 'CCNB1']),
])

EXPECTED_PANELC_GENES = 57   # must equal the figure script's constant

# Panel B genes. Computed for the means cross-check ONLY. They are NOT in the
# Panel C BH family, because they are not in Panel C: their q-values come from
# the figure script's Panel B family. Keeping them here would silently change
# every Panel C q.
PANELB_CROSSCHECK = ['APOBEC3A', 'APOBEC3B']

# Reference values from Generate_Figure6_Lifecycle_Panels.py v6.2, Panel B family
# (18 Mann-Whitney tests across Panels B and D, BH-corrected together).
PANELB_REFERENCE_Q = {'APOBEC3A': 2.6644e-135, 'APOBEC3B': 8.0616e-67}

# Genes evaluated and dropped from Panel C. Means are printed for the audit
# trail in Figure6_PanelC_Tier_Reference.md, but they are EXCLUDED from the BH
# family so the family matches the rendered figure.
CONTEXT_GENES_DROPPED = ['CASP3', 'KRT1', 'CGAS', 'STING1',
                         'APOBEC3C', 'APOBEC3D', 'APOBEC3F',
                         'APOBEC3G', 'APOBEC3H',
                         'IFITM1', 'BST2', 'SMC5', 'SMC6', 'NSMCE2']

# DDX58 is stored as 'RIGI' in this transcriptome (2024-A / GENCODE v44).
# Without this the gene drops silently and the family becomes 56.
GENE_ALIASES = {
    'DDX58':   ['RIGI'],
    'H2AX':    ['H2AFX'],
    'MRE11':   ['MRE11A'],
    'NBN':     ['NBS1'],
    'NSD2':    ['WHSC1', 'MMSET'],
    'TP53BP1': ['TP53BP'],
}


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

# Alias resolution: DDX58 is stored as 'RIGI' in this transcriptome
# (2024-A / GENCODE v44). Rename in place so every downstream lookup and every
# printed table uses the canonical symbol and the BH family stays at 57.
_rename = {}
for _canon, _aliases in GENE_ALIASES.items():
    if _canon not in adata.var_names:
        for _a in _aliases:
            if _a in adata.var_names:
                _rename[_a] = _canon
                break
if _rename:
    adata.var_names = pd.Index([_rename.get(v, v) for v in adata.var_names])
    for _a, _c in _rename.items():
        log(f"  alias resolved: {_c} <- '{_a}'")

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
AUDIT['Maintenance_SBS2_pct'] = 100*np.mean(phase_frac['Maintenance']['SBS2_HIGH'])
AUDIT['Maintenance_CNV_pct']  = 100*np.mean(phase_frac['Maintenance']['CNV_HIGH'])
AUDIT['Capsid_SBS2_pct']      = 100*np.mean(phase_frac['Capsid']['SBS2_HIGH'])
AUDIT['Capsid_CNV_pct']       = 100*np.mean(phase_frac['Capsid']['CNV_HIGH'])
AUDIT['q_Maintenance'] = phase_qvals['Maintenance'][0]
AUDIT['q_Capsid']      = phase_qvals['Capsid'][0]


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

# NOTE for the text: E6E7_frac_of_total here (Mann-Whitney, q = 5.6e-05) and the
# Panel F 'Oncogene' phase fraction (permutation on the difference of means,
# q = 0.10) are the SAME quantity on the SAME cells, tested two ways. The
# permutation test compares means; Mann-Whitney tests distributional shift and
# is far more sensitive to a small consistent offset at n = 197 vs 446. The
# claim that survives either test is the effect size: E6/E7 is under 1% of viral
# reads in both populations. Prose should not assert 'no difference'.
log("\n  NOTE: E6E7_frac_of_total (MW) and the Panel F Oncogene phase fraction")
log("  (permutation) are the same quantity tested two ways and disagree on")
log("  significance. Cite the effect size (<1% of viral reads in both groups),")
log("  not 'no difference'.")


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
# DIAGNOSTIC B: HOST MARKER PANEL  (ungated, all 1,638; BH per contrast)
#   The BH family is exactly the 57 genes rendered in Figure 6c. Panel B genes
#   and dropped-candidate genes are computed separately, outside the family.
# =============================================================================
banner("DIAGNOSTIC B: Host marker panel (ungated, 546/546/546; BH per contrast)")

flat_genes = [(cat, g) for cat, genes in MARKER_GENES.items() for g in genes]
log(f"  Panel C BH family: {len(flat_genes)} genes requested "
    f"(expected {EXPECTED_PANELC_GENES})")
if len(flat_genes) != EXPECTED_PANELC_GENES:
    log(f"  ERROR: MARKER_GENES holds {len(flat_genes)} genes, expected "
        f"{EXPECTED_PANELC_GENES}. Fix before trusting any q-value.")

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
                     'hvc_dir': 'SBS2 > CNV' if means['SBS2_HIGH'] > means['CNV_HIGH'] else 'CNV > SBS2',
                     'peak_pop': max(POP_ORDER, key=lambda p: means[p])}
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
    log(f"  >>> The BH family is now {len(order_genes)}, NOT {EXPECTED_PANELC_GENES}. "
        f"Every q below is wrong until this is resolved.")
else:
    log(f"  All {len(order_genes)} genes resolved; BH family matches Figure 6c.")

log(f"\n  {'Category':<16s} {'Gene':<10s} {'SBS2':>8s} {'CNV':>8s} {'NORM':>8s}  "
    f"{'KW p':>9s} {'HvC q':>9s} {'Dir':>11s} {'Peak':>10s}")
log(f"  {'-'*16} {'-'*10} {'-'*8} {'-'*8} {'-'*8}  {'-'*9} {'-'*9} {'-'*11} {'-'*10}")
for cat, genes in MARKER_GENES.items():
    for gene in genes:
        if gene not in records:
            continue
        r = records[gene]
        log(f"  {cat:<16s} {gene:<10s} {r['mean_SBS2_HIGH']:>8.3f} {r['mean_CNV_HIGH']:>8.3f} "
            f"{r['mean_NORMAL']:>8.3f}  {fmt_p(r['kw_p']):>9s} {fmt_p(r['hvc_bh_q']):>9s} "
            f"{r['hvc_dir']:>11s} {POP_LABELS[r['peak_pop']]:>10s} {stars(r['hvc_bh_q'])}")
pd.DataFrame([records[g] for g in order_genes]).to_csv(
    os.path.join(OUTPUT_DIR, "host_marker_expression_summary.tsv"), sep='\t', index=False)
pd.DataFrame(per_cell_rows).to_csv(
    os.path.join(OUTPUT_DIR, "host_marker_per_cell_values.tsv"), sep='\t', index=False)

# Per-tier peak-direction summary (the structural claim in the Results text)
log(f"\n  Tier peak-direction summary (how many genes peak where, and how many")
log(f"  reach significance in the SBS2-vs-CNV contrast):")
log(f"    {'Tier':<18s} {'n':>3s}  {'peak SBS2':>10s} {'peak CNV':>9s} {'peak NORM':>10s}  {'sig':>4s}")
log(f"    {'-'*18} {'-'*3}  {'-'*10} {'-'*9} {'-'*10}  {'-'*4}")
for cat, genes in MARKER_GENES.items():
    present = [g for g in genes if g in records]
    n_s = sum(records[g]['peak_pop'] == 'SBS2_HIGH' for g in present)
    n_c = sum(records[g]['peak_pop'] == 'CNV_HIGH' for g in present)
    n_n = sum(records[g]['peak_pop'] == 'NORMAL' for g in present)
    n_sig = sum((records[g]['hvc_bh_q'] is not None)
                and (not np.isnan(records[g]['hvc_bh_q']))
                and (records[g]['hvc_bh_q'] < 0.05) for g in present)
    log(f"    {cat:<18s} {len(present):>3d}  {n_s:>10d} {n_c:>9d} {n_n:>10d}  "
        f"{n_sig:>2d}/{len(present):<2d}")

# -----------------------------------------------------------------------------
# OUTSIDE THE BH FAMILY: Panel B cross-check + dropped-candidate audit trail.
# These genes are deliberately excluded from the family above so that the family
# equals the rendered panel. Means are still computed for traceability.
# -----------------------------------------------------------------------------
banner("OUTSIDE THE PANEL C FAMILY: Panel B cross-check + dropped candidates", char="-")

def means_outside_family(gene):
    """Per-population means for a gene NOT in the BH family. Returns dict or None."""
    expr = get_expression(adata_pop, gene)
    if expr is None:
        return None
    vals = {p: expr[(adata_pop.obs['population'] == p).values] for p in POP_ORDER}
    return {'gene': gene,
            **{f'mean_{p}': float(np.mean(vals[p])) for p in POP_ORDER},
            **{f'pct_{p}': 100.0 * np.sum(vals[p] > 0) / max(len(vals[p]), 1)
               for p in POP_ORDER}}

outside_rows = []

log(f"\n  Panel B cross-check (means MUST equal the figure script's Panel B):")
for gene in PANELB_CROSSCHECK:
    m = means_outside_family(gene)
    if m is None:
        log(f"    {gene}: NOT FOUND in adata.var_names")
        continue
    log(f"    {gene}: SBS2 {m['mean_SBS2_HIGH']:.4f}  CNV {m['mean_CNV_HIGH']:.4f}  "
        f"NORM {m['mean_NORMAL']:.4f}   [outside Panel C family; "
        f"Panel B q = {PANELB_REFERENCE_Q.get(gene, float('nan')):.2e}]")
    m['role'] = 'panelB_crosscheck'
    outside_rows.append(m)
    key = 'A3A' if gene == 'APOBEC3A' else 'A3B'
    AUDIT[f'{key}_SBS2'] = m['mean_SBS2_HIGH']
    AUDIT[f'{key}_CNV']  = m['mean_CNV_HIGH']

log(f"\n  Dropped Panel C candidates (audit trail for the tier reference doc;")
log(f"  NOT in the BH family, so no q is reported):")
log(f"    {'Gene':<10s} {'SBS2':>8s} {'CNV':>8s} {'NORM':>8s}")
log(f"    {'-'*10} {'-'*8} {'-'*8} {'-'*8}")
for gene in CONTEXT_GENES_DROPPED:
    m = means_outside_family(gene)
    if m is None:
        log(f"    {gene:<10s} {'--':>8s} {'--':>8s} {'--':>8s}   (not found)")
        continue
    log(f"    {gene:<10s} {m['mean_SBS2_HIGH']:>8.3f} {m['mean_CNV_HIGH']:>8.3f} "
        f"{m['mean_NORMAL']:>8.3f}")
    m['role'] = 'dropped_candidate'
    outside_rows.append(m)

if outside_rows:
    pd.DataFrame(outside_rows).to_csv(
        os.path.join(OUTPUT_DIR, "host_marker_outside_family.tsv"), sep='\t', index=False)

# -----------------------------------------------------------------------------
# Stash host-marker values for the audit. Covers all 57 panel genes, including
# the ns results: an ns gene that later drifts significant, or a q misattributed
# from the three-group Kruskal-Wallis, is only catchable if it is audited.
# -----------------------------------------------------------------------------
def stash(gene, key):
    if gene in records:
        AUDIT[f'{key}_SBS2'] = records[gene]['mean_SBS2_HIGH']
        AUDIT[f'{key}_CNV']  = records[gene]['mean_CNV_HIGH']
        AUDIT[f'{key}_NORM'] = records[gene]['mean_NORMAL']
        AUDIT[f'q_{key}']    = records[gene]['hvc_bh_q']

for gene, key in [
        # MHC-I antigen presentation + IFN signaling/sensing
        ('HLA-A','HLAA'), ('HLA-B','HLAB'), ('HLA-C','HLAC'), ('B2M','B2M'),
        ('TAP1','TAP1'), ('STAT1','STAT1'), ('IRF1','IRF1'), ('STAT2','STAT2'),
        ('DDX58','DDX58'),
        # Type I interferon effectors
        ('IFI27','IFI27'), ('ISG15','ISG15'), ('IRF9','IRF9'), ('MX1','MX1'),
        ('OAS1','OAS1'), ('RSAD2','RSAD2'), ('IFI44L','IFI44L'), ('IFIT1','IFIT1'),
        # Keratinocyte differentiation
        ('KRT5','KRT5'), ('KRT14','KRT14'), ('IVL','IVL'), ('KRT10','KRT10'),
        ('CDH1','CDH1'),
        # DDR, ATM arm
        ('CHEK2','CHEK2'), ('BRCA1','BRCA1'), ('NBN','NBN'), ('H2AX','H2AX'),
        ('BARD1','BARD1'), ('TP53BP1','TP53BP1'), ('RIF1','RIF1'),
        # DDR, post-translationally regulated members
        ('ATM','ATM'), ('MRE11','MRE11'), ('RAD50','RAD50'),
        # DDR, ATR arm + E1-cleavage node + chromatin
        ('TOPBP1','TOPBP1'), ('CHEK1','CHEK1'), ('STAT5A','STAT5A'),
        ('STAT5B','STAT5B'), ('CASP7','CASP7'), ('NSD2','NSD2'),
        # Cell-cycle re-entry / proliferation
        ('MKI67','MKI67'), ('TOP2A','TOP2A'), ('MCM7','MCM7'), ('PCNA','PCNA'),
        ('CCNE1','CCNE1'), ('CDKN2A','CDKN2A'), ('E2F1','E2F1'), ('E2F2','E2F2'),
        ('BRD4','BRD4'), ('MED1','MED1'),
        # p53/Rb
        ('CDKN1A','CDKN1A'), ('MDM2','MDM2'), ('BAX','BAX'), ('TP53','TP53'),
        ('RB1','RB1'),
        # G2/M arrest
        ('CDC25A','CDC25A'), ('CDC25C','CDC25C'), ('CDK1','CDK1'), ('CCNB1','CCNB1'),
]:
    stash(gene, key)


# =============================================================================
# SECTION 3: TEXT NUMBER AUDIT  (diff current Section 4.4 prose vs computed)
# =============================================================================
banner("SECTION 3: Section 4.4 text-number audit")

# Claims hardcoded from the manuscript draft of Section 4.4, updated to the
# 57-gene BH family. Every q in the panel is audited, including ns results.
# kind: 'q' (log10 tol), 'mean' (rel tol), 'pct' (abs tol), 'fold' (abs tol),
#       'count' (exact)
CLAIMS = [
    # ---- Panel F cell set + lifecycle fractions -----------------------------
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
    ('E2 q 0.11',                     0.11,   'q_E2',         'q'),
    ('Oncogene q 0.10 (perm)',        0.10,   'q_Oncogene',   'q'),
    ('Oncogene SBS2 <1% (0.53)',      0.53,   'Oncogene_SBS2_pct', 'pct'),
    ('Oncogene CNV <1% (0.79)',       0.79,   'Oncogene_CNV_pct',  'pct'),
    ('Maintenance SBS2 25.9%',        25.89,  'Maintenance_SBS2_pct', 'pct'),
    ('Maintenance CNV 13.2%',         13.23,  'Maintenance_CNV_pct',  'pct'),
    ('Maintenance q 1.3e-4',          1.3332e-4, 'q_Maintenance', 'q'),
    ('Capsid SBS2 10.3%',             10.29,  'Capsid_SBS2_pct', 'pct'),
    ('Capsid CNV 17.8%',              17.80,  'Capsid_CNV_pct',  'pct'),
    ('Capsid q 1.3e-4',               1.3332e-4, 'q_Capsid',    'q'),

    # ---- Panel B means (outside the Panel C family; q is out-of-scope) ------
    ('A3A SBS2 6.46',                 6.46,   'A3A_SBS2',     'mean'),
    ('A3A CNV 2.08',                  2.08,   'A3A_CNV',      'mean'),
    ('A3B SBS2 2.21',                 2.21,   'A3B_SBS2',     'mean'),
    ('A3B CNV 4.95',                  4.95,   'A3B_CNV',      'mean'),

    # ---- Tier 1: MHC-I antigen presentation + IFN signaling ----------------
    ('B2M q 8.7e-100',                8.67e-100,'q_B2M',      'q'),
    ('HLA-A q 7.3e-42',               7.27e-42,'q_HLAA',      'q'),
    ('HLA-B q 1.2e-15',               1.22e-15,'q_HLAB',      'q'),
    ('HLA-C q 6.1e-8',                6.10e-8, 'q_HLAC',      'q'),
    ('TAP1 q 1.1e-4',                 1.06e-4, 'q_TAP1',      'q'),
    ('STAT1 ns (q 0.31)',             3.07e-1, 'q_STAT1',     'q'),
    ('IRF1 ns (q 0.62)',              6.21e-1, 'q_IRF1',      'q'),
    ('STAT2 ns (q 0.14)',             1.35e-1, 'q_STAT2',     'q'),
    ('DDX58 ns (q 0.094)',            9.35e-2, 'q_DDX58',     'q'),

    # ---- Tier 2: Type I interferon effectors -------------------------------
    ('IFI27 q 2.8e-48',               2.77e-48,'q_IFI27',     'q'),
    ('ISG15 ns (q 0.20)',             1.99e-1, 'q_ISG15',     'q'),
    ('IRF9 q 1.3e-15',                1.25e-15,'q_IRF9',      'q'),
    ('MX1 q 1.4e-13',                 1.37e-13,'q_MX1',       'q'),
    ('OAS1 q 2.3e-30',                2.25e-30,'q_OAS1',      'q'),
    ('RSAD2 q 2.9e-28',               2.94e-28,'q_RSAD2',     'q'),
    ('IFI44L q 4.6e-7',               4.55e-7, 'q_IFI44L',    'q'),
    ('IFIT1 q 7.9e-3',                7.89e-3, 'q_IFIT1',     'q'),

    # ---- Tier 3: Keratinocyte differentiation ------------------------------
    ('KRT5 q 4.6e-7',                 4.55e-7, 'q_KRT5',      'q'),
    ('KRT14 q 8.6e-57',               8.56e-57,'q_KRT14',     'q'),
    ('IVL SBS2 2.68',                 2.68,   'IVL_SBS2',     'mean'),
    ('IVL CNV 0.09',                  0.09,   'IVL_CNV',      'mean'),
    ('IVL q 5.1e-70',                 5.05e-70,'q_IVL',       'q'),
    ('KRT10 q 2.2e-4',                2.15e-4, 'q_KRT10',     'q'),
    ('CDH1 q 2.5e-11',                2.52e-11,'q_CDH1',      'q'),

    # ---- Tier 4: HPV-activated DNA damage response -------------------------
    ('CHEK2 q 1.4e-5',                1.38e-5, 'q_CHEK2',     'q'),
    ('BRCA1 q 1.1e-8',                1.10e-8, 'q_BRCA1',     'q'),
    ('NBN q 1.3e-3',                  1.28e-3, 'q_NBN',       'q'),
    ('H2AX q 4.3e-11',                4.28e-11,'q_H2AX',      'q'),
    ('BARD1 ns (q 0.076)',            7.62e-2, 'q_BARD1',     'q'),
    ('TP53BP1 q 2.4e-6',              2.42e-6, 'q_TP53BP1',   'q'),
    ('RIF1 q 6.7e-11',                6.67e-11,'q_RIF1',      'q'),
    ('ATM q 6.7e-3',                  6.70e-3, 'q_ATM',       'q'),
    ('ATM peaks NORMAL 0.536',        0.536,  'ATM_NORM',     'mean'),
    ('MRE11 q 2.1e-3',                2.05e-3, 'q_MRE11',     'q'),
    ('RAD50 ns (q 0.44)',             4.37e-1, 'q_RAD50',     'q'),
    ('RAD50 peaks NORMAL 1.573',      1.573,  'RAD50_NORM',   'mean'),
    ('TOPBP1 q 8.6e-14',              8.55e-14,'q_TOPBP1',    'q'),
    ('CHEK1 q 7.2e-11',               7.19e-11,'q_CHEK1',     'q'),
    ('STAT5A ns (q 0.31)',            3.07e-1, 'q_STAT5A',    'q'),
    ('STAT5B ns (q 0.34)',            3.41e-1, 'q_STAT5B',    'q'),
    ('CASP7 q 5.7e-3',                5.72e-3, 'q_CASP7',     'q'),
    ('NSD2 q 8.7e-5',                 8.71e-5, 'q_NSD2',      'q'),

    # ---- Tier 5: Cell-cycle re-entry / proliferation -----------------------
    ('MKI67 q 2.5e-19',               2.52e-19,'q_MKI67',     'q'),
    ('TOP2A q 1.2e-25',               1.19e-25,'q_TOP2A',     'q'),
    ('MCM7 q 1.7e-47',                1.67e-47,'q_MCM7',      'q'),
    ('PCNA q 3.9e-13',                3.89e-13,'q_PCNA',      'q'),
    ('CCNE1 q 1.4e-4',                1.40e-4, 'q_CCNE1',     'q'),
    ('CDKN2A ns (q 0.19)',            1.87e-1, 'q_CDKN2A',    'q'),
    ('E2F1 q 1.1e-5',                 1.11e-5, 'q_E2F1',      'q'),
    ('E2F2 q 5.3e-5',                 5.25e-5, 'q_E2F2',      'q'),
    ('BRD4 q 1.3e-11',                1.34e-11,'q_BRD4',      'q'),
    ('MED1 q 1.8e-8',                 1.84e-8, 'q_MED1',      'q'),

    # ---- Tier 6: p53/Rb pathway --------------------------------------------
    ('CDKN1A ns (q 0.19)',            1.87e-1, 'q_CDKN1A',    'q'),
    ('MDM2 q 1.5e-7',                 1.51e-7, 'q_MDM2',      'q'),
    ('BAX q 9.1e-19',                 9.12e-19,'q_BAX',       'q'),
    ('TP53 q 2.6e-11',                2.62e-11,'q_TP53',      'q'),
    ('RB1 ns (q 0.53)',               5.33e-1, 'q_RB1',       'q'),

    # ---- Tier 7: G2/M arrest ------------------------------------------------
    ('CDC25A q 2.5e-16',              2.53e-16,'q_CDC25A',    'q'),
    ('CDC25C q 3.0e-12',              3.02e-12,'q_CDC25C',    'q'),
    ('CDK1 q 4.9e-18',                4.90e-18,'q_CDK1',      'q'),
    ('CCNB1 q 9.5e-27',               9.46e-27,'q_CCNB1',     'q'),
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

log(f"\n  {'Claim':<30s} {'claimed':>12s} {'computed':>14s}   Verdict")
log(f"  {'-'*30} {'-'*12} {'-'*14}   -------")
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
    log(f"  {label:<30s} {cl_str:>12s} {comp_str:>14s}   {v}")
    audit_rows.append({'claim': label, 'claimed': claimed, 'computed': computed, 'verdict': v})

log(f"\n  MATCH: {n_match}   DIFF: {n_diff}   NO VALUE: {n_novalue}   "
    f"(total {len(CLAIMS)})")
if n_diff or n_novalue:
    log("  >>> Inspect any DIFF / NO VALUE rows before the text is finalized.")
else:
    log("  ALL CLAIMS VERIFIED against this run. Section 4.4 numbers are locked to")
    log("  the 57-gene Panel C family drawn in Figure 6c.")

# Out-of-scope numbers, with their correct source.
log(f"\n  OUT OF SCOPE for this diagnostic:")
log(f"    Source = Phase3 L-method / population step:")
log(f"    - 94.6% of HPV16+ cells are basal  (needs non-basal HPV+ counts; this")
log(f"      master table is basal-only, cannot be reconstructed here)")
log(f"    - Tier counts 22,153 / 14,046 / 15,927 (needs the ambiguous-band thresholds)")
log(f"    - Fisher OR = 1.01, p = 0.91 (needs the positivity-vs-SBS2-HIGH contrast set)")
log(f"    Source = Generate_Figure6_Lifecycle_Panels.py, Panel B family of 18:")
log(f"    - A3A q = {PANELB_REFERENCE_Q['APOBEC3A']:.2e} and "
    f"A3B q = {PANELB_REFERENCE_Q['APOBEC3B']:.2e}. These are Panel B genes and")
log(f"      are NOT in the 57-gene Panel C family, so they are not corrected or")
log(f"      audited here. Their means are cross-checked above and must match.")
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
