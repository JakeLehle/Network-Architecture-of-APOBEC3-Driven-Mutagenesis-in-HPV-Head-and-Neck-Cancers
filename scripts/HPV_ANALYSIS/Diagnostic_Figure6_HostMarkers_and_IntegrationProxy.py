#!/usr/bin/env python3
"""
Diagnostic_Figure6_HostMarkers_and_IntegrationProxy.py  (v5 -- 59-gene panel)
===================================================================================
Figure 6 diagnostic + Section 4.3 text-number verification.

v5 changes
----------
  1. SECTION RENUMBER. Everything previously labelled "Section 4.4" is now
     "Section 4.3". The host-marker / lifecycle material is Results 4.3 in the
     current manuscript draft; the old label was audit-against-the-wrong-section
     waiting to happen. Output file renamed section4_3_text_audit.tsv.

  2. NEW TIER 8: BET_BRD4_axis = [BRD4, BRD3, BRD2].
     Rationale (Wu et al. Mol Cell 2024, PMID 38103559): pBRD4 recruits the DDR
     factors 53BP1 (TP53BP1) and BARD1 to the HPV origin of replication in a
     BRD4-L / BRD4-S isoform-specific manner, and BRD3 -- but explicitly NOT
     BRD2 -- is required for differentiation-associated HPV genome amplification
     (their Fig 6F/6G). BRD2 is therefore carried as an internal negative
     control, which is only legible if the three sit in their own tier.
     BRD4 MOVED OUT of CellCycle_Prolif (10 -> 9 genes) into this tier, because
     Wu et al. place BRD4 upstream of the DDR factors at the viral ori rather
     than as a proliferation marker.
     Panel C family: 57 -> 59. EXPECTED_PANELC_GENES updated accordingly.

     >>> Results 4.3 prose consequences, all of which need re-checking:
         - "57 host genes across seven functional tiers" -> 59 / eight
         - "45 differed significantly" -> read the new count off the tier
           summary printed by DIAGNOSTIC B
         - "7 of the 10 genes in that tier" (cell-cycle re-entry) -> that tier
           is now 9 genes; re-read the peak-direction summary
         - EVERY q-value in the section shifts, because the BH family grew from
           57 to 59 (roughly a 3.5% inflation). See item 4.

  3. NEW SECTION 4: BRD4 ISOFORM RATIO.
     Reads the per-SRR table written by Diagnostic_BRD4_Isoform_Population_Ratio.py
     rather than touching BAMs, so this harness stays fast and single-purpose.
     Computes, for BRD4-S (CCDS46004 / O60885-2) versus BRD4-L (CCDS12328 /
     O60885-1):
       - pooled S fraction per population
       - the SBS2-vs-CNV risk ratio with a 95% CI (the number the manuscript
         should quote; the p-value alone is not informative for a null)
       - Fisher exact on the pooled 2x2
       - Cochran-Mantel-Haenszel stratified by patient (or by SRR if no patient
         column is available), because the pooled test treats ~800 molecules
         from 36 samples across 14 patients as independent and SBS2-HIGH is
         74% drawn from three patients. The pooled CI is therefore optimistically
         narrow; CMH is the honest version.
     S(b) (CCDS82307) is reported but EXCLUDED from the ratio: at 794 aa it is
     UniProt O60885-3, a third isoform, not a short variant of BRD4-S.
     Identity confirmed by Diagnostic_BRD4_CCDS_Isoform_Identity.py.

  4. CLAIMS EMITTER. Growing the BH family invalidates all 57 previously locked
     q-values at once, so the audit would report DIFF on everything and tell you
     nothing. After the audit, this script now prints a copy-pasteable CLAIMS
     block built from the freshly computed values.
     WORKFLOW, run twice:
       run 1: expect a wall of DIFF. Copy the emitted CLAIMS block.
       paste:  replace the CLAIMS list below with the emitted block.
       run 2:  expect ALL MATCH. Only now are the numbers locked, and only now
               should the manuscript text be updated from them.
     Do NOT hand-transcribe 59 q-values.

v4 sections retained unchanged
------------------------------
  SECTION 1: LIFECYCLE FRACTIONS. Mirrors Generate_Figure6_Lifecycle_Panels.py
    Panel F EXACTLY: gated HPV16-positive set (raw_HPV16 >= 8 AND total > 0),
    per-cell gene fractions = gene / total (no pseudocount), permutation test on
    the difference of means (10,000 perms, seed 42), BH-FDR within the 8-gene
    family and separately within the 4-phase family.

  SECTION 2: READ-CLASS / URR BREAKDOWN. Per group on the gated set, the
    URR / ORF / intergenic read fractions BOTH ways: pooled and per-cell mean.
    These differ for CNV-HIGH (~63.5% pooled vs ~67.1% per-cell mean); the
    prose cites pooled.

  SECTION 3: TEXT NUMBER AUDIT. Diffs the Section 4.3 prose (hardcoded below)
    against freshly computed values and prints MATCH / DIFF / NO VALUE.

  DIAGNOSTIC A: integration proxy on the gated >=8 set.
  VIRAL LOAD SUMMARY: raw_HPV16 (all cells) vs total reads (gated set).
  DIAGNOSTIC B: host-marker panel, ungated 546/546/546, BH per contrast.

  Gene aliases are still resolved at load time. DDX58 is stored as 'RIGI' in
  this transcriptome (2024-A / GENCODE v44); without this it drops silently and
  the family shrinks, shifting every q. The same guard now protects BRD3/BRD2.

  APOBEC3A / APOBEC3B remain OUTSIDE the Panel C BH family (they are Panel B
  genes); their means are computed here only as a cross-check.

INPUTS
  - data/FIG_4/01_group_selection/three_group_assignments.tsv
  - data/FIG_4/00_input/adata_final.h5ad
  - data/FIG_6/01_raw_hpv16_counts/basal_cell_master_table_with_raw_HPV16.tsv
  - data/FIG_6/03_hpv16_genome/per_cell_hpv16_gene_counts.tsv
  - data/FIG_6/DIAGNOSTIC_BRD4_ISOFORM/brd4_population_ratio_by_srr.tsv   [v5]

OUTPUTS (to data/FIG_6/DIAGNOSTIC_LIFECYCLE_MARKERS/)
  - diagnostic_figure6_report.txt
  - integration_proxy_metrics.tsv
  - viral_load_summary.tsv
  - host_marker_expression_summary.tsv
  - host_marker_per_cell_values.tsv
  - host_marker_outside_family.tsv
  - lifecycle_fractions_panelF_mirror.tsv
  - readclass_urr_breakdown.tsv
  - brd4_isoform_ratio_summary.tsv          [v5]
  - section4_3_text_audit.tsv               [renamed in v5]
  - emitted_claims_block.py                 [v5]

Env: NETWORK
Usage: conda run -n NETWORK python Diagnostic_Figure6_HostMarkers_and_IntegrationProxy.py

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
import re
from scipy.stats import mannwhitneyu, kruskal, fisher_exact
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.contingency_tables import StratifiedTable
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

# v5: written by Diagnostic_BRD4_Isoform_Population_Ratio.py
BRD4_ISOFORM_TSV = os.path.join(PROJECT_ROOT,
    "data/FIG_6/DIAGNOSTIC_BRD4_ISOFORM/brd4_population_ratio_by_srr.tsv")

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
# HOST MARKER GENE PANEL (59 genes, 8 tiers)   [v5]
# Must mirror DOTPLOT_CATEGORIES in Generate_Figure6_Lifecycle_Panels.py
# exactly, so the BH family here equals the panel drawn in Figure 6.
# Locked structure: Figure6_PanelC_Tier_Reference.md  (UPDATE THAT DOC TOO)
#
# Tier order below sets the row-group order in the dot plot. BET_BRD4_axis is
# placed immediately after the DDR tier so the BRD4 -> 53BP1/BARD1 relationship
# reads top to bottom. Move it above DDR_ATM_ATR if you prefer strict
# mechanistic ordering (BRD4 is upstream of the factors it recruits).
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
    # v5 NEW TIER. BRD4 relocated here from CellCycle_Prolif.
    # BRD2 is the internal negative control (Wu et al. Fig 6F: BRD3 knockdown
    # suppresses HPV genome amplification, BRD2 knockdown does not).
    ('BET_BRD4_axis',   ['BRD4', 'BRD3', 'BRD2']),
    ('CellCycle_Prolif', ['MKI67', 'TOP2A', 'MCM7', 'PCNA', 'CCNE1',
                          'CDKN2A', 'E2F1', 'E2F2', 'MED1']),
    ('p53_Rb_pathway',  ['CDKN1A', 'MDM2', 'BAX', 'TP53', 'RB1']),
    ('G2M_arrest',      ['CDC25A', 'CDC25C', 'CDK1', 'CCNB1']),
])

EXPECTED_PANELC_GENES = 59   # v5: was 57. Must equal the figure script's constant.

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
                         'IFITM1', 'BST2', 'SMC5', 'SMC6', 'NSMCE2',
                         # v5: Wu et al. BRD4-L DSB interactors deliberately NOT
                         # added to the family. RAD21 and NIPBL are cohesin and
                         # will track proliferation for reasons unrelated to
                         # BRD4, so a significant result would be uninterpretable.
                         'RAD21', 'NIPBL']

# DDX58 is stored as 'RIGI' in this transcriptome (2024-A / GENCODE v44).
# Without this the gene drops silently and the family shrinks.
GENE_ALIASES = {
    'DDX58':   ['RIGI'],
    'H2AX':    ['H2AFX'],
    'MRE11':   ['MRE11A'],
    'NBN':     ['NBS1'],
    'NSD2':    ['WHSC1', 'MMSET'],
    'TP53BP1': ['TP53BP'],
    # v5 additions. No alias is expected to be needed for the BET genes in
    # GENCODE v44, but if either fails to resolve the missing-gene guard in
    # DIAGNOSTIC B will fire and every q below it is invalid until fixed.
    'BRD3':    ['RING3L'],
    'BRD2':    ['RING3', 'FSRG1'],
}

# --- v5: BRD4 isoform constants ---------------------------------------------
# Identities confirmed by Diagnostic_BRD4_CCDS_Isoform_Identity.py against the
# pipeline GTF and the UniProt O60885 CCDS cross-references:
#   CCDS12328 -> O60885-1, 1362 aa, BRD4-L
#   CCDS46004 -> O60885-2,  722 aa, BRD4-S   <- the Wu et al. short isoform
#   CCDS82307 -> O60885-3,  794 aa           <- a THIRD isoform, not BRD4-S
BRD4_ISO_LABELS = {'L': 'BRD4-L (CCDS12328 / O60885-1)',
                   'Sa': 'BRD4-S (CCDS46004 / O60885-2)',
                   'Sb': 'O60885-3 (CCDS82307), excluded from the ratio'}


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

def _first(cols, candidates):
    """Return the first candidate column name present in cols, else None."""
    lower = {str(c).strip().lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None

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

# Alias resolution. Rename in place so every downstream lookup and every printed
# table uses the canonical symbol and the BH family stays at EXPECTED_PANELC_GENES.
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
# SECTION 1: LIFECYCLE FRACTIONS  (mirror of figure Panel F)
# =============================================================================
banner("SECTION 1: Lifecycle fractions (Panel F mirror; per-cell frac = gene/total)")

for g in ALL_HPV_GENES:
    hpv_pos[f'{g}_frac'] = hpv_pos[g] / hpv_pos[TOTAL_COL]

gene_frac = {}
for g in ALL_HPV_GENES:
    gene_frac[g] = {p: hpv_pos.loc[hpv_pos['population'] == p, f'{g}_frac'].values.astype(float)
                    for p in POP_ORDER}

gene_raw = []
for g in ALL_HPV_GENES:
    gene_raw.extend(compute_pairwise_perm(gene_frac[g]))
gene_q = bh(gene_raw)
gene_qvals = {g: gene_q[i*3:(i+1)*3] for i, g in enumerate(ALL_HPV_GENES)}

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
log("  Prose currently quotes 63.5 / 63.5 / 64.9. The two tumour groups agreeing")
log("  to one decimal is a coincidence worth re-confirming, not a copy error --")
log("  check the pooled column above matches before the text is finalised.")


# =============================================================================
# DIAGNOSTIC A: INTEGRATION PROXY  (gated >=8 set, figure-matched)
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

log("\n  NOTE: E6E7_frac_of_total (MW) and the Panel F Oncogene phase fraction")
log("  (permutation) are the same quantity tested two ways and disagree on")
log("  significance. Cite the effect size (<1% of viral reads in both groups),")
log("  not 'no difference'.")


# =============================================================================
# VIRAL LOAD SUMMARY
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
#   The BH family is exactly the EXPECTED_PANELC_GENES rendered in Figure 6.
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
    log(f"  All {len(order_genes)} genes resolved; BH family matches Figure 6.")

log(f"\n  {'Category':<18s} {'Gene':<10s} {'SBS2':>8s} {'CNV':>8s} {'NORM':>8s}  "
    f"{'KW p':>9s} {'HvC q':>9s} {'Dir':>11s} {'Peak':>10s}")
log(f"  {'-'*18} {'-'*10} {'-'*8} {'-'*8} {'-'*8}  {'-'*9} {'-'*9} {'-'*11} {'-'*10}")
for cat, genes in MARKER_GENES.items():
    for gene in genes:
        if gene not in records:
            continue
        r = records[gene]
        log(f"  {cat:<18s} {gene:<10s} {r['mean_SBS2_HIGH']:>8.3f} {r['mean_CNV_HIGH']:>8.3f} "
            f"{r['mean_NORMAL']:>8.3f}  {fmt_p(r['kw_p']):>9s} {fmt_p(r['hvc_bh_q']):>9s} "
            f"{r['hvc_dir']:>11s} {POP_LABELS[r['peak_pop']]:>10s} {stars(r['hvc_bh_q'])}")
pd.DataFrame([records[g] for g in order_genes]).to_csv(
    os.path.join(OUTPUT_DIR, "host_marker_expression_summary.tsv"), sep='\t', index=False)
pd.DataFrame(per_cell_rows).to_csv(
    os.path.join(OUTPUT_DIR, "host_marker_per_cell_values.tsv"), sep='\t', index=False)

# Per-tier peak-direction summary (the structural claim in the Results text)
log(f"\n  Tier peak-direction summary. The Results 4.3 sentences 'N of {EXPECTED_PANELC_GENES}")
log(f"  differed significantly' and 'X of the Y genes in that tier' come from here.")
log(f"    {'Tier':<18s} {'n':>3s}  {'peak SBS2':>10s} {'peak CNV':>9s} {'peak NORM':>10s}  {'sig':>6s}")
log(f"    {'-'*18} {'-'*3}  {'-'*10} {'-'*9} {'-'*10}  {'-'*6}")
total_sig = 0
total_n = 0
for cat, genes in MARKER_GENES.items():
    present = [g for g in genes if g in records]
    n_s = sum(records[g]['peak_pop'] == 'SBS2_HIGH' for g in present)
    n_c = sum(records[g]['peak_pop'] == 'CNV_HIGH' for g in present)
    n_n = sum(records[g]['peak_pop'] == 'NORMAL' for g in present)
    n_sig = sum((records[g]['hvc_bh_q'] is not None)
                and (not np.isnan(records[g]['hvc_bh_q']))
                and (records[g]['hvc_bh_q'] < 0.05) for g in present)
    total_sig += n_sig
    total_n += len(present)
    log(f"    {cat:<18s} {len(present):>3d}  {n_s:>10d} {n_c:>9d} {n_n:>10d}  "
        f"{n_sig:>2d}/{len(present):<3d}")
log(f"    {'-'*18} {'-'*3}  {'-'*10} {'-'*9} {'-'*10}  {'-'*6}")
log(f"    {'TOTAL':<18s} {total_n:>3d}  {'':>10s} {'':>9s} {'':>10s}  {total_sig:>2d}/{total_n:<3d}")
log(f"\n  >>> Results 4.3 should read '{total_sig} differed significantly' out of "
    f"{total_n} host genes.")
AUDIT['panelC_family_n'] = total_n
AUDIT['panelC_n_sig'] = total_sig

# -----------------------------------------------------------------------------
# OUTSIDE THE BH FAMILY: Panel B cross-check + dropped-candidate audit trail.
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
# Stash host-marker values for the audit. Covers ALL panel genes including ns
# results: an ns gene that later drifts significant, or a q misattributed from
# the three-group Kruskal-Wallis, is only catchable if it is audited.
# (BARD1 was previously recorded at q=1e-41 in the tier reference doc, which was
#  the three-group KW p, not the SBS2-vs-CNV contrast. That is why ns genes are
#  audited.)
# -----------------------------------------------------------------------------
def stash(gene, key):
    if gene in records:
        AUDIT[f'{key}_SBS2'] = records[gene]['mean_SBS2_HIGH']
        AUDIT[f'{key}_CNV']  = records[gene]['mean_CNV_HIGH']
        AUDIT[f'{key}_NORM'] = records[gene]['mean_NORMAL']
        AUDIT[f'q_{key}']    = records[gene]['hvc_bh_q']

# key is the gene symbol with non-identifier characters stripped
STASH_KEY = {'HLA-A': 'HLAA', 'HLA-B': 'HLAB', 'HLA-C': 'HLAC'}
for _cat, _genes in MARKER_GENES.items():
    for _g in _genes:
        stash(_g, STASH_KEY.get(_g, _g))


# =============================================================================
# SECTION 4: BRD4 ISOFORM RATIO  [v5]
#   BRD4-S (CCDS46004 / O60885-2) vs BRD4-L (CCDS12328 / O60885-1).
#   Reads the per-SRR table from Diagnostic_BRD4_Isoform_Population_Ratio.py.
# =============================================================================
banner("SECTION 4: BRD4-S vs BRD4-L isoform ratio (from the isoform script's TSV)")

AUDIT['brd4_iso_available'] = False

if not os.path.exists(BRD4_ISOFORM_TSV):
    log(f"  [SKIP] not found: {BRD4_ISOFORM_TSV}")
    log( "  Run Diagnostic_BRD4_Isoform_Population_Ratio.py first. Note that it")
    log( "  currently WARNS and continues when a BAM index is missing; index")
    log( "  every sample with the conda samtools before quoting these numbers,")
    log( "  or the value depends on which .bai files happened to exist at run time.")
else:
    iso = pd.read_csv(BRD4_ISOFORM_TSV, sep='\t')
    log(f"  loaded {len(iso)} rows from {os.path.basename(BRD4_ISOFORM_TSV)}")
    log(f"  columns: {list(iso.columns)}")

    c_pop = _first(iso.columns, ['population', 'group', 'pop'])
    c_srr = _first(iso.columns, ['srr', 'sample', 'sample_id', 'srr_id', 'run'])
    c_pat = _first(iso.columns, ['patient', 'donor', 'patient_id', 'donor_id', 'subject'])
    c_L   = _first(iso.columns, ['L', 'L_umis', 'L_count', 'n_L', 'umis_L'])
    c_Sa  = _first(iso.columns, ['Sa', 'S_a', 'Sa_umis', 'Sa_count', 'n_Sa', 'umis_Sa', 'S(a)'])
    c_Sb  = _first(iso.columns, ['Sb', 'S_b', 'Sb_umis', 'Sb_count', 'n_Sb', 'umis_Sb', 'S(b)'])

    if not all([c_pop, c_L, c_Sa]):
        log("  ERROR: could not resolve the population / L / S(a) columns.")
        log("  Expected something like: population, srr, L, Sa, Sb")
        log("  Rename the columns in the isoform script or extend the _first() lists.")
    else:
        strat_col = c_pat or c_srr
        strat_kind = 'patient' if c_pat else ('SRR' if c_srr else None)
        if strat_kind == 'SRR':
            log("  NOTE: no patient column found; stratifying by SRR instead.")
            log("  Samples nest within patients, so this is conservative but not")
            log("  identical to patient-level stratification. Add a patient column")
            log("  to the isoform script's output if you want the exact version.")
        elif strat_kind is None:
            log("  NOTE: no patient or SRR column; CMH cannot be computed.")

        iso[c_L]  = pd.to_numeric(iso[c_L],  errors='coerce').fillna(0)
        iso[c_Sa] = pd.to_numeric(iso[c_Sa], errors='coerce').fillna(0)
        if c_Sb:
            iso[c_Sb] = pd.to_numeric(iso[c_Sb], errors='coerce').fillna(0)

        # normalise population labels to POP_ORDER
        norm = {'sbs2-high': 'SBS2_HIGH', 'sbs2_high': 'SBS2_HIGH',
                'cnv-high': 'CNV_HIGH', 'cnv_high': 'CNV_HIGH',
                'normal': 'NORMAL'}
        iso['_pop'] = iso[c_pop].astype(str).str.strip().str.lower().map(norm)
        unmapped = iso.loc[iso['_pop'].isna(), c_pop].unique()
        if len(unmapped):
            log(f"  WARNING: unmapped population labels dropped: {list(unmapped)}")
        iso = iso[iso['_pop'].notna()].copy()

        # ---- pooled per population ----
        log(f"\n  {'Population':<12s} {'L':>7s} {'S(a)':>7s} {'S(b)':>7s} "
            f"{'S(a) frac':>10s}  {'L/S(a)':>8s}")
        log(f"  {'-'*12} {'-'*7} {'-'*7} {'-'*7} {'-'*10}  {'-'*8}")
        pooled = {}
        iso_rows = []
        for p in POP_ORDER:
            sub = iso[iso['_pop'] == p]
            nL  = int(sub[c_L].sum())
            nSa = int(sub[c_Sa].sum())
            nSb = int(sub[c_Sb].sum()) if c_Sb else 0
            denom = nL + nSa
            frac = nSa / denom if denom else np.nan
            ratio = (nL / nSa) if nSa else np.nan
            pooled[p] = {'L': nL, 'Sa': nSa, 'Sb': nSb, 'frac': frac}
            log(f"  {POP_LABELS[p]:<12s} {nL:>7d} {nSa:>7d} {nSb:>7d} "
                f"{100*frac:>9.1f}%  {ratio:>8.2f}")
            iso_rows.append({'population': POP_LABELS[p], 'L': nL, 'Sa': nSa, 'Sb': nSb,
                             'Sa_fraction': frac, 'L_over_Sa': ratio})

        log(f"\n  {BRD4_ISO_LABELS['Sb']}")
        log( "  S(b) counts are reported for completeness only. At 794 aa this is a")
        log( "  third isoform, not a short variant of BRD4-S, so it is excluded")
        log( "  from the ratio and must not be pooled with S(a).")
        log( "\n  The L : S(a) ratio is NOT an abundance ratio. Unique-interval")
        log( "  lengths differ (L 4658 bp, S(a) 2244 bp) under 3'-biased chemistry,")
        log( "  so capture geometry contributes. The across-population comparison")
        log( "  is valid because the same intervals are used in all three groups.")
        log( "  Never place this number next to the Wu et al. protein-level '<3%'.")

        # ---- SBS2 vs CNV risk ratio with 95% CI (THE number for the text) ----
        a1, n1 = pooled['SBS2_HIGH']['Sa'], pooled['SBS2_HIGH']['Sa'] + pooled['SBS2_HIGH']['L']
        a2, n2 = pooled['CNV_HIGH']['Sa'],  pooled['CNV_HIGH']['Sa']  + pooled['CNV_HIGH']['L']
        if a1 > 0 and a2 > 0 and n1 > 0 and n2 > 0:
            p1, p2 = a1 / n1, a2 / n2
            rr = p1 / p2
            se_log = np.sqrt((1 - p1) / (n1 * p1) + (1 - p2) / (n2 * p2))
            lo = float(np.exp(np.log(rr) - 1.96 * se_log))
            hi = float(np.exp(np.log(rr) + 1.96 * se_log))
        else:
            rr = lo = hi = np.nan

        odds, fisher_p = fisher_exact([[a1, n1 - a1], [a2, n2 - a2]])

        log(f"\n  SBS2-HIGH vs CNV-HIGH, BRD4-S fraction:")
        log(f"    SBS2-HIGH  {a1}/{n1} = {100*a1/n1:.1f}%")
        log(f"    CNV-HIGH   {a2}/{n2} = {100*a2/n2:.1f}%")
        log(f"    risk ratio = {rr:.3f}   95% CI {lo:.3f} to {hi:.3f}")
        log(f"    Fisher exact p = {fisher_p:.3f}")
        log(f"\n  >>> QUOTE THE RATIO AND CI, NOT THE p-VALUE. A p-value cannot")
        log(f"      express what a null excludes; the CI can. The manuscript")
        log(f"      sentence is 'did not differ (RR {rr:.2f}, 95% CI {lo:.2f} to {hi:.2f})',")
        log(f"      plus the exclusion bound this implies.")
        log(f"  >>> The correct verb is 'does not argue against' the BRD4 axis.")
        log(f"      NOT 'supports' and NOT 'confirms'. A null is equally consistent")
        log(f"      with the recruitment-level model, with insufficient power, and")
        log(f"      with a modest transcriptional difference below detection.")

        # ---- Cochran-Mantel-Haenszel stratified by patient (or SRR) ----
        cmh_or = cmh_lo = cmh_hi = cmh_p = np.nan
        n_strata_used = 0
        if strat_kind is not None:
            tables = []
            strata_detail = []
            for key, sub in iso.groupby(strat_col):
                s = sub[sub['_pop'].isin(['SBS2_HIGH', 'CNV_HIGH'])]
                if s.empty:
                    continue
                g1 = s[s['_pop'] == 'SBS2_HIGH']
                g2 = s[s['_pop'] == 'CNV_HIGH']
                t = np.array([[g1[c_Sa].sum(), g1[c_L].sum()],
                              [g2[c_Sa].sum(), g2[c_L].sum()]], dtype=float)
                # a stratum contributes nothing unless both arms have molecules
                if t[0].sum() > 0 and t[1].sum() > 0 and t[:, 0].sum() > 0 and t[:, 1].sum() > 0:
                    tables.append(t)
                    strata_detail.append({
                        'stratum': key,
                        'SBS2_Sa': int(t[0, 0]), 'SBS2_L': int(t[0, 1]),
                        'CNV_Sa':  int(t[1, 0]), 'CNV_L':  int(t[1, 1]),
                    })
            n_strata_used = len(tables)

            # WHICH strata carry the CMH, and how much of the data they hold.
            # If these collapse onto one patient, the CMH is a within-patient
            # estimate and cannot be reported as a population comparison.
            if strata_detail:
                sd = pd.DataFrame(strata_detail)
                held = sd[['SBS2_Sa', 'SBS2_L', 'CNV_Sa', 'CNV_L']].to_numpy().sum()
                total_mol = (pooled['SBS2_HIGH']['L'] + pooled['SBS2_HIGH']['Sa'] +
                             pooled['CNV_HIGH']['L'] + pooled['CNV_HIGH']['Sa'])
                log(f"\n  INFORMATIVE STRATA ({n_strata_used}), stratified by {strat_kind}:")
                log(f"    {'stratum':<16s} {'SBS2 Sa/L':>12s} {'CNV Sa/L':>12s} "
                    f"{'SBS2 Sa%':>9s} {'CNV Sa%':>9s}")
                log(f"    {'-'*16} {'-'*12} {'-'*12} {'-'*9} {'-'*9}")
                for r in strata_detail:
                    s_tot = r['SBS2_Sa'] + r['SBS2_L']
                    c_tot = r['CNV_Sa'] + r['CNV_L']
                    log(f"    {str(r['stratum']):<16s} "
                        f"{r['SBS2_Sa']:>5d}/{r['SBS2_L']:<6d} "
                        f"{r['CNV_Sa']:>5d}/{r['CNV_L']:<6d} "
                        f"{100*r['SBS2_Sa']/s_tot if s_tot else float('nan'):>8.1f}% "
                        f"{100*r['CNV_Sa']/c_tot if c_tot else float('nan'):>8.1f}%")
                log(f"    These strata hold {held} of {total_mol} tumour molecules "
                    f"({100*held/total_mol:.1f}%).")
                log(f"    >>> The other {100 - 100*held/total_mol:.1f}% contributes")
                log(f"        NOTHING to the CMH. If these strata map to one patient,")
                log(f"        the CMH is a within-patient estimate and MUST NOT be")
                log(f"        reported as an SBS2-vs-CNV population comparison.")
                sd.to_csv(os.path.join(OUTPUT_DIR, "brd4_cmh_strata_detail.tsv"),
                          sep='\t', index=False)            
            if n_strata_used >= 2:
                try:
                    st = StratifiedTable(tables)
                    cmh_or = float(st.oddsratio_pooled)
                    cmh_lo, cmh_hi = [float(v) for v in st.oddsratio_pooled_confint()]
                    cmh_p = float(st.test_null_odds().pvalue)
                except Exception as e:
                    log(f"    CMH failed ({e})")
            log(f"\n  Cochran-Mantel-Haenszel, stratified by {strat_kind} "
                f"({n_strata_used} informative strata):")
            if np.isfinite(cmh_or):
                log(f"    pooled OR = {cmh_or:.3f}   95% CI {cmh_lo:.3f} to {cmh_hi:.3f}   "
                    f"p = {cmh_p:.3f}")
                log(f"    Compare against the unstratified CI above. The pooled test")
                log(f"    treats every molecule as independent even though SBS2-HIGH is")
                log(f"    74% drawn from three patients, so the unstratified interval is")
                log(f"    optimistically narrow. IF THE CMH INTERVAL IS WIDER, QUOTE IT.")
            else:
                log(f"    not enough informative strata to compute a pooled estimate")

        AUDIT['brd4_iso_available'] = True
        AUDIT['brd4_Sa_frac_SBS2'] = 100 * a1 / n1 if n1 else np.nan
        AUDIT['brd4_Sa_frac_CNV']  = 100 * a2 / n2 if n2 else np.nan
        AUDIT['brd4_rr']    = rr
        AUDIT['brd4_rr_lo'] = lo
        AUDIT['brd4_rr_hi'] = hi
        AUDIT['brd4_fisher_p'] = fisher_p
        AUDIT['brd4_cmh_or'] = cmh_or
        AUDIT['brd4_cmh_p']  = cmh_p

        iso_rows.append({'population': 'SBS2_vs_CNV', 'L': np.nan, 'Sa': np.nan, 'Sb': np.nan,
                         'Sa_fraction': np.nan, 'L_over_Sa': np.nan,
                         'risk_ratio': rr, 'rr_ci_lo': lo, 'rr_ci_hi': hi,
                         'fisher_p': fisher_p, 'cmh_or': cmh_or,
                         'cmh_ci_lo': cmh_lo, 'cmh_ci_hi': cmh_hi, 'cmh_p': cmh_p,
                         'cmh_strata': n_strata_used, 'cmh_stratified_by': strat_kind})
        pd.DataFrame(iso_rows).to_csv(
            os.path.join(OUTPUT_DIR, "brd4_isoform_ratio_summary.tsv"), sep='\t', index=False)


# =============================================================================
# SECTION 3: TEXT NUMBER AUDIT  (diff current Section 4.3 prose vs computed)
# =============================================================================
banner("SECTION 3: Section 4.3 text-number audit")

log("  !! v5 WARNING: the CLAIMS list below is the v4 lock, built on a 57-gene")
log("     BH family. This run uses a 59-gene family, so EVERY host q-value has")
log("     shifted (roughly 3.5% inflation) and will report DIFF. That is")
log("     expected on the first v5 run. Copy the emitted CLAIMS block printed")
log("     at the end of this script, paste it over the list below, and re-run.")
log("     Only after the second run reports ALL MATCH are the numbers locked")
log("     and safe to copy into the manuscript.")

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

    # ---- Panel C family size + significant count (v5: was 45 of 57) --------
    ('Panel C family n = 59',         59,     'panelC_family_n', 'count'),

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

    # ---- Tier 5 (v5 NEW): BET / BRD4 axis ----------------------------------
    # BRD4 relocated from the cell-cycle tier; BRD3 and BRD2 are new. No prior
    # lock exists for BRD3/BRD2, so these are placeholders that WILL report
    # DIFF or NO VALUE on run 1. The emitted block supplies the real values.
    ('BRD4 q 1.3e-11 (was cellcycle tier)', 1.34e-11, 'q_BRD4', 'q'),
    ('BRD3 q PLACEHOLDER',            1.0,    'q_BRD3',       'q'),
    ('BRD2 q PLACEHOLDER',            1.0,    'q_BRD2',       'q'),

    # ---- Tier 6: Cell-cycle re-entry / proliferation (v5: BRD4 removed) ----
    ('MKI67 q 2.5e-19',               2.52e-19,'q_MKI67',     'q'),
    ('TOP2A q 1.2e-25',               1.19e-25,'q_TOP2A',     'q'),
    ('MCM7 q 1.7e-47',                1.67e-47,'q_MCM7',      'q'),
    ('PCNA q 3.9e-13',                3.89e-13,'q_PCNA',      'q'),
    ('CCNE1 q 1.4e-4',                1.40e-4, 'q_CCNE1',     'q'),
    ('CDKN2A ns (q 0.19)',            1.87e-1, 'q_CDKN2A',    'q'),
    ('E2F1 q 1.1e-5',                 1.11e-5, 'q_E2F1',      'q'),
    ('E2F2 q 5.3e-5',                 5.25e-5, 'q_E2F2',      'q'),
    ('MED1 q 1.8e-8',                 1.84e-8, 'q_MED1',      'q'),

    # ---- Tier 7: p53/Rb pathway --------------------------------------------
    ('CDKN1A ns (q 0.19)',            1.87e-1, 'q_CDKN1A',    'q'),
    ('MDM2 q 1.5e-7',                 1.51e-7, 'q_MDM2',      'q'),
    ('BAX q 9.1e-19',                 9.12e-19,'q_BAX',       'q'),
    ('TP53 q 2.6e-11',                2.62e-11,'q_TP53',      'q'),
    ('RB1 ns (q 0.53)',               5.33e-1, 'q_RB1',       'q'),

    # ---- Tier 8: G2/M arrest ------------------------------------------------
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
    if kind == 'ratio':
        return 'MATCH' if abs(computed - claimed) <= 0.02 else 'DIFF'
    return '?'

log(f"\n  {'Claim':<38s} {'claimed':>12s} {'computed':>14s}   Verdict")
log(f"  {'-'*38} {'-'*12} {'-'*14}   -------")
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
    log(f"  {label:<38s} {cl_str:>12s} {comp_str:>14s}   {v}")
    audit_rows.append({'claim': label, 'claimed': claimed, 'computed': computed, 'verdict': v})

log(f"\n  MATCH: {n_match}   DIFF: {n_diff}   NO VALUE: {n_novalue}   "
    f"(total {len(CLAIMS)})")
if n_diff or n_novalue:
    log("  >>> Expected on the FIRST v5 run. Paste the emitted block below and re-run.")
else:
    log("  ALL CLAIMS VERIFIED against this run. Section 4.3 numbers are locked to")
    log(f"  the {EXPECTED_PANELC_GENES}-gene Panel C family drawn in Figure 6.")

log(f"\n  OUT OF SCOPE for this diagnostic:")
log(f"    Source = Phase3 L-method / population step:")
log(f"    - 94.6% of HPV16+ cells are basal")
log(f"    - Tier counts 22,153 / 14,046 / 15,927")
log(f"    - Fisher OR = 1.01, p = 0.91")
log(f"    Source = Generate_Figure6_Lifecycle_Panels.py, Panel B family of 18:")
log(f"    - A3A q = {PANELB_REFERENCE_Q['APOBEC3A']:.2e} and "
    f"A3B q = {PANELB_REFERENCE_Q['APOBEC3B']:.2e}")
n_pos_allbasal = int((master['raw_HPV16'] >= HPV16_THRESHOLD).sum())
log(f"    Partial anchor: raw_HPV16 >= {HPV16_THRESHOLD} over all {len(master)} basal "
    f"cells = {n_pos_allbasal}  (compare to the tier 'positive' = 15,927)")

pd.DataFrame(audit_rows).to_csv(os.path.join(OUTPUT_DIR, "section4_3_text_audit.tsv"),
                                sep='\t', index=False)


# =============================================================================
# CLAIMS EMITTER  [v5]
#   Rebuilds the CLAIMS list from this run's computed values so the 59-gene
#   family can be re-locked without hand-transcribing every q.
# =============================================================================
banner("CLAIMS EMITTER: copy the block below over the CLAIMS list, then re-run")

def _emit_val(key, kind):
    v = AUDIT.get(key)
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return None
    if kind == 'q':
        return f"{v:.4g}"
    if kind == 'count':
        return f"{int(round(v))}"
    return f"{v:.6g}"

# A digit immediately preceded by a letter belongs to a gene symbol (BRD3, SBS2,
# B2M, CDC25A, E2F1), not to a value. Without the lookbehind, 'BRD3 q' rewrites
# to 'BRD1.8e-25 q'.
_NUM_RE = re.compile(r'(?<![A-Za-z])[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?')

def _refresh_label(label, val):
    """Rewrite the number embedded in a claim label so it matches the new value.
    Replaces the LAST numeric token; appends if the label carries none."""
    lab = label.replace(' PLACEHOLDER', '').rstrip()
    hits = list(_NUM_RE.finditer(lab))
    if hits:
        last = hits[-1]
        return lab[:last.start()] + str(val) + lab[last.end():]
    return f"{lab} {val}"

emit = []
emit.append("CLAIMS = [")
skipped = []
for label, claimed, key, kind in CLAIMS:
    val = _emit_val(key, kind)
    if val is None:
        skipped.append((label, key))
        emit.append(f"    # NO VALUE this run, kept the previous claim:")
        emit.append(f"    ({label!r}, {claimed!r}, {key!r}, {kind!r}),")
        continue
    new_label = _refresh_label(label, val)
    emit.append(f"    ({new_label!r}, {val}, {key!r}, {kind!r}),")

# v5 isoform claims, only if the isoform TSV was present
if AUDIT.get('brd4_iso_available'):
    for lbl, key, kind in [
            ('BRD4-S frac SBS2 (%)',  'brd4_Sa_frac_SBS2', 'pct'),
            ('BRD4-S frac CNV (%)',   'brd4_Sa_frac_CNV',  'pct'),
            ('BRD4-S risk ratio',     'brd4_rr',           'ratio'),
            ('BRD4-S RR CI low',      'brd4_rr_lo',        'ratio'),
            ('BRD4-S RR CI high',     'brd4_rr_hi',        'ratio'),
    ]:
        val = _emit_val(key, kind)
        if val is not None:
            emit.append(f"    ({_refresh_label(lbl, val)!r}, {val}, {key!r}, {kind!r}),")
emit.append("]")

for line in emit:
    log("  " + line)

if skipped:
    log(f"\n  {len(skipped)} claim(s) had NO VALUE and kept their previous entry:")
    for lbl, key in skipped:
        log(f"    {lbl}  (AUDIT key '{key}' was never set)")
    log("  Investigate each one; a NO VALUE usually means a gene failed to resolve")
    log("  or an upstream section was skipped, not that the number is fine.")

emit_path = os.path.join(OUTPUT_DIR, "emitted_claims_block.py")
with open(emit_path, 'w') as f:
    f.write("# Emitted by Diagnostic_Figure6_HostMarkers_and_IntegrationProxy.py v5\n")
    f.write("# Paste over the CLAIMS list in that script, then re-run to confirm ALL MATCH.\n")
    f.write("\n".join(emit) + "\n")
log(f"\n  [SAVE] {emit_path}")


# =============================================================================
# SAVE REPORT
# =============================================================================
banner("COMPLETE")
log("  NEXT STEPS")
log("   1. If run 1: paste emitted_claims_block.py over CLAIMS and re-run.")
log("   2. Confirm ALL MATCH on run 2.")
log("   3. Update Generate_Figure6_Lifecycle_Panels.py DOTPLOT_CATEGORIES to the")
log("      same 8-tier / 59-gene structure, or the figure and the q-values diverge.")
log("   4. Update Figure6_PanelC_Tier_Reference.md.")
log("   5. Only then update Results 4.3 from the locked values.")
report_path = os.path.join(OUTPUT_DIR, "diagnostic_figure6_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"\n  Report saved: {report_path}")
log(f"  Output directory: {OUTPUT_DIR}")
