#!/usr/bin/env python3
"""
diagnostic_section6_verify.py
=============================

Independent verification that Generate_Figure6_Lifecycle_Panels.py produced
the correct Panel B, D, and F statistics.

This does NOT trust any remembered numbers. It re-derives every B/D/F value
from the same source files the figure script reads, using the figure script's
EXACT methodology, and then diffs that independent recompute against the
numbers the figure actually wrote into figure6_lifecycle_report.txt.

Methodology mirrored from the figure script (must match for the check to be valid):
  - Panel F cells = master raw_HPV16 >= 8 (the L-method elbow) AND
    total_hpv16_genome_reads > 0   (the 197 / 383 / 8 = 588 set)
  - Panel B + D: Mann-Whitney U, BH-FDR corrected together as ONE family (18 tests)
  - Panel F per-gene: permutation test on difference of means (10,000 perms,
    seed 42), BH-FDR within the 8 genes (24 entries, 16 N.D., 8 valid)
  - Panel F per-phase: same permutation test, BH-FDR within the 4 phases
    (12 entries, 8 N.D., 4 valid)  -- separate family from the genes
  - SBS2-weight NORMAL group uses only cells present in the signature-weights
    file (the n=447 quirk), reproduced here.

Output:
  - Independent recompute of every B/D/F/phase q-value, mean fraction, and URR
    fraction, printed in the figure's format.
  - A PASS/FAIL diff against the parsed figure report (best-effort parse; if a
    value cannot be found in the report it is reported as such rather than
    silently passing).
  - data/FIG_6/TROUBLESHOOTING/diagnostic_section6_verify_report.txt

Env: NETWORK
Usage: conda run -n NETWORK python diagnostic_section6_verify.py

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import re
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from collections import OrderedDict
from datetime import datetime
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
SIG_WEIGHTS_PATH = os.path.join(PROJECT_ROOT,
    "data/FIG_4/00_input/signature_weights_per_cell.txt")
MASTER_TABLE_PATH = os.path.join(PROJECT_ROOT,
    "data/FIG_6/01_raw_hpv16_counts/basal_cell_master_table_with_raw_HPV16.tsv")
HPV_GENE_PATH = os.path.join(PROJECT_ROOT,
    "data/FIG_6/03_hpv16_genome/per_cell_hpv16_gene_counts.tsv")

REPORT_PATH = os.path.join(PROJECT_ROOT,
    "data/FIG_6/FIGURE_6_PANELS/figure6_lifecycle_report.txt")

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_6/TROUBLESHOOTING")
os.makedirs(OUTPUT_DIR, exist_ok=True)

HPV16_THRESHOLD = 8
TOTAL_COL = 'total_hpv16_genome_reads'
MIN_CELLS_FOR_STATS = 10

POP_ORDER = ['SBS2_HIGH', 'CNV_HIGH', 'NORMAL']
POP_LABELS = {'SBS2_HIGH': 'SBS2-HIGH', 'CNV_HIGH': 'CNV-HIGH', 'NORMAL': 'Normal'}

HPV16_PHASES = OrderedDict([
    ('Maintenance',   ['E1', 'E2']),
    ('Amplification', ['E4', 'E5']),
    ('Oncogene',      ['E6', 'E7']),
    ('Capsid',        ['L1', 'L2']),
])
ALL_HPV_GENES = [g for genes in HPV16_PHASES.values() for g in genes]

# Canonical pair order, matching the figure's pair_names
PAIR_NAMES = ['SBS2-HIGH vs CNV-HIGH', 'CNV-HIGH vs Normal', 'SBS2-HIGH vs Normal']

# Comparison tolerances
QTOL_RTOL = 1e-3      # q-values printed to ~4 sig figs in the report
PCT_ATOL = 0.02       # mean fractions printed to 2 decimals (percentage points)

# =============================================================================
# LOGGING
# =============================================================================
report_lines = []
def log(msg=""):
    print(msg, flush=True)
    report_lines.append(str(msg))

def banner(title, char="="):
    log("")
    log(char * 80)
    log(f"  {title}")
    log(char * 80)

# =============================================================================
# HELPERS  (copied verbatim from the figure script so methodology is identical)
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

def compute_pairwise_pvals(data_dict, pop_order):
    pairs = [(0, 1), (1, 2), (0, 2)]
    out = []
    for i, j in pairs:
        v1 = data_dict[pop_order[i]]
        v2 = data_dict[pop_order[j]]
        if (len(v1) >= MIN_CELLS_FOR_STATS and len(v2) >= MIN_CELLS_FOR_STATS
                and len(v1) > 5 and len(v2) > 5):
            _, p = mannwhitneyu(v1, v2, alternative='two-sided')
            out.append(p)
        else:
            out.append(np.nan)
    return out

def permutation_test_means(v1, v2, n_perm=10000, seed=42):
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

def compute_pairwise_perm(data_dict, pop_order, n_perm=10000, seed=42):
    pairs = [(0, 1), (1, 2), (0, 2)]
    out = []
    for i, j in pairs:
        v1 = data_dict[pop_order[i]]
        v2 = data_dict[pop_order[j]]
        if len(v1) >= MIN_CELLS_FOR_STATS and len(v2) >= MIN_CELLS_FOR_STATS:
            out.append(permutation_test_means(v1, v2, n_perm=n_perm, seed=seed))
        else:
            out.append(np.nan)
    return out

def apply_bh_correction(raw_flat):
    pvals = np.array(raw_flat, dtype=float)
    adjusted = np.full_like(pvals, np.nan)
    valid = ~np.isnan(pvals)
    if valid.sum() == 0:
        return adjusted.tolist()
    _, adj_p, _, _ = multipletests(pvals[valid], method='fdr_bh')
    adjusted[valid] = adj_p
    return adjusted.tolist()

def stars(q):
    if np.isnan(q): return 'N.D.'
    if q < 1e-4: return '****'
    if q < 1e-3: return '***'
    if q < 0.01: return '**'
    if q < 0.05: return '*'
    return 'ns'

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
adata.obs['population'] = 'other'
adata.obs.loc[adata.obs_names.isin(sbs2_cells), 'population'] = 'SBS2_HIGH'
adata.obs.loc[adata.obs_names.isin(cnv_cells), 'population'] = 'CNV_HIGH'
adata.obs.loc[adata.obs_names.isin(normal_cells), 'population'] = 'NORMAL'
adata_pop = adata[adata.obs['population'].isin(POP_ORDER)].copy()
log(f"  Cells in three populations: {adata_pop.shape[0]}")

# Signature weights (with the same transpose detection as the figure)
sig_weights = pd.read_csv(SIG_WEIGHTS_PATH, sep='\t', index_col=0)
sbs2_col = None
for c in sig_weights.columns:
    if str(c).strip() == 'SBS2':
        sbs2_col = c; break
if sbs2_col is None:
    for c in sig_weights.columns:
        if 'SBS2' in str(c):
            sbs2_col = c; break
if sbs2_col is None and 'SBS2' in sig_weights.index:
    sig_weights = sig_weights.T
    sbs2_col = 'SBS2'
log(f"  SBS2 column: '{sbs2_col}'  (sig_weights shape {sig_weights.shape})")

# Master table
master = pd.read_csv(MASTER_TABLE_PATH, sep='\t', index_col=0)
master['group'] = master.index.map(lambda x: cell_to_group.get(x, 'other'))
master_pop = master[master['group'].isin(POP_ORDER)].copy()

# HPV gene counts
hpv_genes = pd.read_csv(HPV_GENE_PATH, sep='\t', index_col=0)
hpv_genes['population'] = 'other'
hpv_genes.loc[hpv_genes.index.isin(sbs2_cells), 'population'] = 'SBS2_HIGH'
hpv_genes.loc[hpv_genes.index.isin(cnv_cells), 'population'] = 'CNV_HIGH'
hpv_genes.loc[hpv_genes.index.isin(normal_cells), 'population'] = 'NORMAL'

# =============================================================================
# STEP 1: BUILD DATA DICTS  (mirrors the figure script's STEP 1)
# =============================================================================
banner("STEP 1: Recompute data dicts")

all_violin_data = OrderedDict()

# Panel B
sbs2_data = {}
for pop in POP_ORDER:
    cells = sbs2_cells if pop == 'SBS2_HIGH' else (cnv_cells if pop == 'CNV_HIGH' else normal_cells)
    overlap = cells & set(sig_weights.index)
    sbs2_data[pop] = (sig_weights.loc[list(overlap), sbs2_col].values.astype(float)
                      if (sbs2_col and overlap) else np.array([0.0]))
all_violin_data['B_SBS2'] = sbs2_data

for col, label in [('cnv_score', 'B_CNV'), ('CytoTRACE2_Score', 'B_CytoTRACE2')]:
    d = {}
    for pop in POP_ORDER:
        m = adata_pop.obs['population'] == pop
        d[pop] = adata_pop.obs.loc[m, col].values.astype(float)
    all_violin_data[label] = d

for gene, label in [('APOBEC3A', 'B_A3A'), ('APOBEC3B', 'B_A3B')]:
    expr = get_expression(adata_pop, gene)
    d = {}
    for pop in POP_ORDER:
        m = (adata_pop.obs['population'] == pop).values
        d[pop] = expr[m] if expr is not None else np.array([0.0])
    all_violin_data[label] = d

# Panel D
hpv_dist = {}
for pop in POP_ORDER:
    m = master_pop['group'] == pop
    hpv_dist[pop] = np.log1p(master_pop.loc[m, 'raw_HPV16'].values.astype(float))
all_violin_data['D_HPV16_dist'] = hpv_dist

for label in all_violin_data:
    ns = ", ".join(f"{POP_LABELS[p]} n={len(all_violin_data[label][p])}" for p in POP_ORDER)
    log(f"  {label:<14s}: {ns}")

# Panel F cell set: EXACT figure logic (raw_HPV16 >= 8 from master AND total > 0)
hpv_pop = hpv_genes[hpv_genes['population'].isin(POP_ORDER)].copy()
hpv_pop['raw_HPV16'] = hpv_pop.index.map(master['raw_HPV16'])
hpv_pop = hpv_pop[(hpv_pop['raw_HPV16'] >= HPV16_THRESHOLD) &
                  (hpv_pop[TOTAL_COL] > 0)].copy()
recomp_F_counts = {pop: int((hpv_pop['population'] == pop).sum()) for pop in POP_ORDER}
log(f"\n  Panel F cell set (raw_HPV16 >= {HPV16_THRESHOLD} AND {TOTAL_COL} > 0): "
    f"{len(hpv_pop)} total")
for pop in POP_ORDER:
    log(f"    {POP_LABELS[pop]}: {recomp_F_counts[pop]}")

for gene in ALL_HPV_GENES:
    hpv_pop[f'{gene}_frac'] = (hpv_pop[gene] / hpv_pop[TOTAL_COL]
                               if gene in hpv_pop.columns else 0.0)

for gene in ALL_HPV_GENES:
    d = {}
    for pop in POP_ORDER:
        m = hpv_pop['population'] == pop
        d[pop] = hpv_pop.loc[m, f'{gene}_frac'].values.astype(float)
    all_violin_data[f'F_{gene}'] = d

phase_frac_data = OrderedDict()
for phase, genes in HPV16_PHASES.items():
    cols = [f'{g}_frac' for g in genes]
    hpv_pop[f'{phase}_frac'] = hpv_pop[cols].sum(axis=1)
    d = {}
    for pop in POP_ORDER:
        m = hpv_pop['population'] == pop
        d[pop] = hpv_pop.loc[m, f'{phase}_frac'].values.astype(float)
    phase_frac_data[phase] = d

# =============================================================================
# STEP 2: STATISTICS  (exact family partition from the figure script)
# =============================================================================
banner("STEP 2: Recompute statistics (B+D one MW family; F genes; F phases)")

mw_labels = ['B_SBS2', 'B_CNV', 'B_CytoTRACE2', 'B_A3A', 'B_A3B', 'D_HPV16_dist']
mw_raw, mw_slices = [], {}
for label in mw_labels:
    s = len(mw_raw)
    mw_raw.extend(compute_pairwise_pvals(all_violin_data[label], POP_ORDER))
    mw_slices[label] = (s, len(mw_raw))
mw_adj = apply_bh_correction(mw_raw)

perm_raw, perm_slices = [], {}
for gene in ALL_HPV_GENES:
    s = len(perm_raw)
    perm_raw.extend(compute_pairwise_perm(all_violin_data[f'F_{gene}'], POP_ORDER))
    perm_slices[f'F_{gene}'] = (s, len(perm_raw))
perm_adj = apply_bh_correction(perm_raw)

phase_raw, phase_slices = [], {}
for phase in HPV16_PHASES:
    s = len(phase_raw)
    phase_raw.extend(compute_pairwise_perm(phase_frac_data[phase], POP_ORDER))
    phase_slices[phase] = (s, len(phase_raw))
phase_adj = apply_bh_correction(phase_raw)

qvals = {}
for label in mw_labels:
    a, b = mw_slices[label]; qvals[label] = mw_adj[a:b]
for gene in ALL_HPV_GENES:
    a, b = perm_slices[f'F_{gene}']; qvals[f'F_{gene}'] = perm_adj[a:b]
phase_qvals = {ph: phase_adj[slice(*phase_slices[ph])] for ph in HPV16_PHASES}

n_mw_valid = int(np.sum(~np.isnan(mw_raw)))
n_perm_valid = int(np.sum(~np.isnan(perm_raw)))
n_phase_valid = int(np.sum(~np.isnan(phase_raw)))
log(f"  Mann-Whitney family (B+D): {n_mw_valid} valid / {len(mw_raw)} "
    f"(expected 18 / 18)")
log(f"  Permutation family (F genes): {n_perm_valid} valid / {len(perm_raw)} "
    f"(expected 8 / 24)")
log(f"  Permutation family (F phases): {n_phase_valid} valid / {len(phase_raw)} "
    f"(expected 4 / 12)")

# Mean fractions + URR
mean_frac = {}
for gene in ALL_HPV_GENES:
    vd = all_violin_data[f'F_{gene}']
    mean_frac[gene] = {p: (100 * np.mean(vd[p]) if len(vd[p]) else np.nan) for p in POP_ORDER}
for phase in HPV16_PHASES:
    pd_ = phase_frac_data[phase]
    mean_frac[phase] = {p: (100 * np.mean(pd_[p]) if len(pd_[p]) else np.nan) for p in POP_ORDER}
urr_frac = {}
for pop in POP_ORDER:
    sub = hpv_pop[hpv_pop['population'] == pop]
    urr_frac[pop] = 100 * (sub['URR'] / sub[TOTAL_COL]).mean() if len(sub) else np.nan

# =============================================================================
# STEP 3: PRINT INDEPENDENT RECOMPUTE
# =============================================================================
banner("STEP 3: Independent recompute (authoritative, from source files)")

log(f"\n  {'Label':<16s}  {'Pair':<22s}  {'Test':>5s}  {'BH q':>13s}  Stars")
log(f"  {'-'*16}  {'-'*22}  {'-'*5}  {'-'*13}  -----")
for label in mw_labels:
    for k in range(3):
        q = qvals[label][k]
        qs = 'N.D.' if np.isnan(q) else f"{q:.4e}"
        log(f"  {label:<16s}  {PAIR_NAMES[k]:<22s}  {'MW':>5s}  {qs:>13s}  {stars(q)}")
for gene in ALL_HPV_GENES:
    for k in range(3):
        q = qvals[f'F_{gene}'][k]
        qs = 'N.D.' if np.isnan(q) else f"{q:.4e}"
        log(f"  {'F_'+gene:<16s}  {PAIR_NAMES[k]:<22s}  {'perm':>5s}  {qs:>13s}  {stars(q)}")

log(f"\n  PHASE q (SBS2 vs CNV) and MEAN FRACTIONS")
log(f"  {'Item':<14s}  {'SBS2-HIGH':>10s}  {'CNV-HIGH':>10s}  {'Normal':>10s}  {'q SBS2vCNV':>13s}")
log(f"  {'-'*14}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*13}")
for gene in ALL_HPV_GENES:
    m = mean_frac[gene]; q = qvals[f'F_{gene}'][0]
    log(f"  {gene:<14s}  {m['SBS2_HIGH']:>9.2f}%  {m['CNV_HIGH']:>9.2f}%  "
        f"{m['NORMAL']:>9.2f}%  {q:>13.4e}")
log(f"  {'-'*14}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*13}")
for phase in HPV16_PHASES:
    m = mean_frac[phase]; q = phase_qvals[phase][0]
    log(f"  {phase:<14s}  {m['SBS2_HIGH']:>9.2f}%  {m['CNV_HIGH']:>9.2f}%  "
        f"{m['NORMAL']:>9.2f}%  {q:>13.4e}")
log("")
for pop in POP_ORDER:
    log(f"  URR mean per-cell fraction {POP_LABELS[pop]}: {urr_frac[pop]:.1f}%")

# =============================================================================
# STEP 4: PARSE FIGURE REPORT
# =============================================================================
banner("STEP 4: Parse figure report for comparison")

report_q = {}        # (label, pair_idx) -> q
report_meanfrac = {} # item -> (sbs2, cnv, normal, q)
report_urr = {}      # pop -> pct
report_F_counts = {} # pop -> int
report_found = os.path.exists(REPORT_PATH)

def _to_float(tok):
    tok = tok.strip().rstrip('%')
    if tok in ('N.D.', 'nan', ''):
        return np.nan
    try:
        return float(tok)
    except ValueError:
        return None

if not report_found:
    log(f"  [WARNING] Report not found at {REPORT_PATH}")
    log(f"  Skipping diff; the recompute above is the authoritative result.")
else:
    log(f"  Reading {REPORT_PATH}")
    lines = open(REPORT_PATH).read().splitlines()

    pair_to_idx = {p: i for i, p in enumerate(PAIR_NAMES)}

    in_qtable = False
    in_meanfrac = False
    for line in lines:
        s = line.strip()

        # Panel F cell-count lines under "HPV16-positive cells (>= 8 UMI) ..."
        for pop in POP_ORDER:
            mobj = re.match(rf'^{POP_LABELS[pop]}:\s+(\d+)\s*$', s)
            if mobj:
                report_F_counts.setdefault(pop, int(mobj.group(1)))

        # URR lines
        mobj = re.match(r'^URR mean per-cell fraction (.+?):\s+([\d.]+)%$', s)
        if mobj:
            lbl = mobj.group(1)
            for pop in POP_ORDER:
                if POP_LABELS[pop] == lbl:
                    report_urr[pop] = float(mobj.group(2))

        # Detect tables
        if 'Item' in s and 'q SBS2vCNV' in s:
            in_meanfrac = True; in_qtable = False; continue
        if 'Panel' in s and 'Pair' in s and 'Raw p' in s and 'BH q' in s:
            in_qtable = True; in_meanfrac = False; continue
        if set(s) <= set('- '):  # separator / blank
            continue

        if in_meanfrac:
            parts = s.split()
            # gene/phase rows: Item  pct% pct% pct%  q
            if len(parts) >= 5 and parts[0] in (ALL_HPV_GENES + list(HPV16_PHASES)):
                vals = [_to_float(parts[i]) for i in (1, 2, 3, 4)]
                if all(v is not None for v in vals):
                    report_meanfrac[parts[0]] = tuple(vals)
                continue
            else:
                in_meanfrac = False  # table ended

        if in_qtable:
            parts = s.split()
            # label  popA vs popB  test  raw  q  stars
            if len(parts) >= 7 and parts[2] == 'vs':
                label = parts[0]
                pair = f"{parts[1]} vs {parts[3]}"
                q = _to_float(parts[6])
                if pair in pair_to_idx and q is not None:
                    report_q[(label, pair_to_idx[pair])] = q
                continue
            else:
                in_qtable = False

    log(f"  Parsed: {len(report_q)} q-table entries, "
        f"{len(report_meanfrac)} mean-fraction rows, "
        f"{len(report_urr)} URR values, {len(report_F_counts)} F-count lines")

# =============================================================================
# STEP 5: DIFF
# =============================================================================
banner("STEP 5: Recompute vs figure report")

n_pass = n_fail = n_missing = 0

def check(name, recomp, reported, rtol=QTOL_RTOL, atol=0.0):
    global n_pass, n_fail, n_missing
    if reported is None:
        n_missing += 1
        log(f"  [MISSING] {name}: recompute={recomp}, not found in report")
        return
    # both NaN -> match
    if (isinstance(recomp, float) and np.isnan(recomp)) and \
       (isinstance(reported, float) and np.isnan(reported)):
        n_pass += 1
        log(f"  [PASS] {name}: N.D. == N.D.")
        return
    if (isinstance(recomp, float) and np.isnan(recomp)) != \
       (isinstance(reported, float) and np.isnan(reported)):
        n_fail += 1
        log(f"  [FAIL] {name}: recompute={recomp}, report={reported}")
        return
    if np.isclose(recomp, reported, rtol=rtol, atol=atol):
        n_pass += 1
        log(f"  [PASS] {name}: {recomp:.4g} ~ {reported:.4g}")
    else:
        n_fail += 1
        log(f"  [FAIL] {name}: recompute={recomp:.6g}, report={reported:.6g}")

if report_found:
    # Panel F cell counts (the elbow set)
    log("\n  -- Panel F cell counts (raw_HPV16 >= 8 set) --")
    for pop in POP_ORDER:
        check(f"F count {POP_LABELS[pop]}", float(recomp_F_counts[pop]),
              float(report_F_counts[pop]) if pop in report_F_counts else None,
              rtol=0, atol=0)

    # B + D q-values (all three pairs)
    log("\n  -- Panel B + D BH q-values --")
    for label in mw_labels:
        for k in range(3):
            check(f"{label} [{PAIR_NAMES[k]}]", qvals[label][k],
                  report_q.get((label, k)))

    # F per-gene q-values (all three pairs)
    log("\n  -- Panel F per-gene BH q-values --")
    for gene in ALL_HPV_GENES:
        for k in range(3):
            check(f"F_{gene} [{PAIR_NAMES[k]}]", qvals[f'F_{gene}'][k],
                  report_q.get((f'F_{gene}', k)))

    # F per-phase q (SBS2 vs CNV) + mean fractions
    log("\n  -- Panel F phase q (SBS2 vs CNV) + mean fractions --")
    for phase in HPV16_PHASES:
        rep = report_meanfrac.get(phase)
        check(f"{phase} q SBS2vCNV", phase_qvals[phase][0],
              rep[3] if rep else None)
        if rep:
            check(f"{phase} mean SBS2", mean_frac[phase]['SBS2_HIGH'], rep[0],
                  rtol=0, atol=PCT_ATOL)
            check(f"{phase} mean CNV", mean_frac[phase]['CNV_HIGH'], rep[1],
                  rtol=0, atol=PCT_ATOL)
        else:
            check(f"{phase} mean SBS2", mean_frac[phase]['SBS2_HIGH'], None)

    # Per-gene mean fractions (SBS2 + CNV)
    log("\n  -- Panel F per-gene mean fractions --")
    for gene in ALL_HPV_GENES:
        rep = report_meanfrac.get(gene)
        if rep:
            check(f"{gene} mean SBS2", mean_frac[gene]['SBS2_HIGH'], rep[0],
                  rtol=0, atol=PCT_ATOL)
            check(f"{gene} mean CNV", mean_frac[gene]['CNV_HIGH'], rep[1],
                  rtol=0, atol=PCT_ATOL)
        else:
            check(f"{gene} mean SBS2", mean_frac[gene]['SBS2_HIGH'], None)

    # URR
    log("\n  -- URR per-cell fractions --")
    for pop in POP_ORDER:
        check(f"URR {POP_LABELS[pop]}", urr_frac[pop],
              report_urr.get(pop), rtol=0, atol=0.1)

    banner("SUMMARY")
    log(f"  PASS:    {n_pass}")
    log(f"  FAIL:    {n_fail}")
    log(f"  MISSING: {n_missing}  (value not located in report)")
    if n_fail == 0 and n_missing == 0:
        log("\n  VERDICT: figure report matches independent recompute. CONFIRMED.")
    elif n_fail == 0:
        log("\n  VERDICT: no mismatches, but some values were not found in the "
            "report (check parsing or report completeness).")
    else:
        log("\n  VERDICT: MISMATCHES FOUND. Inspect [FAIL] lines above; the figure "
            "report disagrees with a fresh recompute from source.")

# =============================================================================
# SAVE
# =============================================================================
out_path = os.path.join(OUTPUT_DIR, "diagnostic_section6_verify_report.txt")
with open(out_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"\n  Report saved: {out_path}")
banner("VERIFICATION COMPLETE")
