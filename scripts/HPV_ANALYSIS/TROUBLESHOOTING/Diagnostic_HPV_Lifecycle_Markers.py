#!/usr/bin/env python3
"""
Diagnostic_HPV_Lifecycle_Markers.py
=====================================
Figure 6 Diagnostic: HPV16 Lifecycle Phase and Host Marker Gene Profiling

PURPOSE:
  Comprehensive diagnostic to inform Figure 6 panel design. Three analyses:

  DIAGNOSTIC 1 — HPV16 gene expression by four functional phases
    Maintenance (E1, E2), Amplification (E4, E5), Oncogene (E6, E7),
    Capsid (L1, L2) across the three populations.

  DIAGNOSTIC 2 — Integration proxy metrics
    E2:(E6+E7) ratio, total HPV16 reads per cell, per-gene read fractions.
    Tests the prediction: SBS2-HIGH = integrated (low E2:E6), CNV-HIGH = episomal (high E2:E6).

  DIAGNOSTIC 3 — Host marker gene expression panel
    Canonical HPV field markers organized by biological process:
    transformation, p53/Rb pathway, ATM/DNA damage, cell cycle arrest,
    caspase, immune signaling, differentiation, innate immune/APOBEC regulation.

POPULATIONS (from Figure 4 Step00B, n=546 each):
  SBS2-HIGH:  High SBS2 weight, A3A-dominant, low CNV
  CNV-HIGH:   SBS2=0, high CNV, A3-matched, late gene fraction
  NORMAL:     Normal adjacent tissue basal cells

INPUTS:
  - data/FIG_4/01_group_selection/three_group_assignments.tsv
  - data/FIG_4/00_input/adata_final.h5ad
  - data/FIG_6/03_hpv16_genome/per_cell_hpv16_gene_counts.tsv

OUTPUTS (to data/FIG_6/DIAGNOSTIC_LIFECYCLE_MARKERS/):
  - diagnostic_lifecycle_report.txt
  - hpv16_gene_by_phase_and_population.tsv
  - integration_proxy_metrics.tsv
  - host_marker_expression_summary.tsv
  - host_marker_pct_expressing.tsv
  - host_marker_per_cell_values.tsv

Env: NETWORK
Usage: conda run -n NETWORK python Diagnostic_HPV_Lifecycle_Markers.py

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from scipy.stats import mannwhitneyu, kruskal, spearmanr
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"

# Population definitions (from Figure 4)
THREE_GROUP_PATH = os.path.join(PROJECT_ROOT,
    "data/FIG_4/01_group_selection/three_group_assignments.tsv")

# Single-cell expression data
ADATA_PATH = os.path.join(PROJECT_ROOT,
    "data/FIG_4/00_input/adata_final.h5ad")

# HPV16 genome-level read counts (from Phase 4)
HPV_GENE_PATH = os.path.join(PROJECT_ROOT,
    "data/FIG_6/03_hpv16_genome/per_cell_hpv16_gene_counts.tsv")

# Output
OUTPUT_DIR = os.path.join(PROJECT_ROOT,
    "data/FIG_6/DIAGNOSTIC_LIFECYCLE_MARKERS")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================================================
# HPV16 GENE FUNCTIONAL GROUPINGS (revised four-phase model)
# =============================================================================
HPV16_PHASES = {
    'Maintenance':   ['E1', 'E2'],       # Episome replication, E6/E7 suppression
    'Amplification': ['E4', 'E5'],       # Genome amplification, immune evasion
    'Oncogene':      ['E6', 'E7'],       # Cell cycle deregulation, p53/Rb inactivation
    'Capsid':        ['L1', 'L2'],       # Virion assembly (productive infection)
}

# =============================================================================
# HOST MARKER GENE PANEL (organized by biological process)
# Sources: Doorbar (lifecycle), Laimins (ATM/DNA damage), McBride (integration)
# =============================================================================
MARKER_GENES = {
    # --- Transformation/integration indicators ---
    # Predict HIGH in SBS2-HIGH (integrated, deregulated E6/E7)
    'Transformation': [
        'CDKN2A',   # p16, THE clinical marker for E7-driven Rb inactivation (Doorbar, field standard)
        'MCM7',     # Diagnostic panel marker alongside p16 (Doorbar)
        'CCNE1',    # Cyclin E, direct E2F target from Rb inactivation (Laimins: E7->Rb->E2F)
        'MKI67',    # Ki-67, universal proliferation marker (field standard)
        'TOP2A',    # Proliferation marker
        'PCNA',     # Proliferation / replication
        'BRD4',     # Super-enhancer / transcriptional activation at integration sites (McBride)
        'MED1',     # Mediator complex, super-enhancer marker (McBride)
        'E2F1',     # Rb pathway target, interacts with MRN complex (Laimins)
        'E2F2',     # E7-induced, directs MRN to HPV origins (Laimins 2009)
    ],

    # --- p53/Rb pathway readouts ---
    # Predict LOW in SBS2-HIGH (E6 degrades p53, E7 inactivates Rb)
    'p53_Rb_pathway': [
        'CDKN1A',   # p21, p53 target — predict LOW if E6 degrades p53 (Laimins/McBride)
        'BAX',      # p53 apoptosis target
        'MDM2',     # p53 regulator (complicated by E6 interaction)
        'RB1',      # Direct E7 target
        'TP53',     # Direct E6 target (mRNA may not change, protein does)
    ],

    # --- ATM / DNA damage cascade ---
    # Active in BOTH populations but different downstream consequences
    # Laimins 2009: ATM activated by E7, required for productive replication
    'ATM_DNA_damage': [
        'ATM',      # Master kinase, binds E7 via LXCXE domain (Laimins 2009)
        'CHEK2',    # CHK2, phosphorylated by ATM, required for genome amplification (Laimins 2009)
        'BRCA1',    # Phosphorylated by ATM in differentiating HPV+ cells (Laimins 2009)
        'NBN',      # NBS1, MRN complex — phosphorylation concurrent with amplification (Laimins 2009)
        'MRE11',    # MRE11, MRN complex component (Laimins 2009) [HUGO: MRE11]
        'RAD50',    # MRN complex component (Laimins 2009)
        'H2AX',     # gamma-H2AX, DNA double-strand break marker (Laimins 2009) [HUGO: H2AX]
        'CHEK1',    # ATR target, phosphorylated independently of ATM (Laimins 2009)
        'STAT5A',   # ATM/DNA damage signaling (field: upregulated while STAT1 down)
        'STAT5B',   # ATM/DNA damage signaling
    ],

    # --- G2/M arrest machinery ---
    # Predict enriched in CNV-HIGH (productive replication requires G2/M arrest)
    # Laimins 2009: CHK2 -> CDC25C phosphorylation -> G2/M arrest -> genome amplification
    'G2M_arrest': [
        'CDC25A',   # CHK2 target phosphatase (Laimins 2009)
        'CDC25C',   # CHK2 target, G2/M arrest mediator (Laimins 2009)
        'CDK1',     # Cdc2, G2/M entry kinase — inhibitory phosphorylation in HPV+ (Laimins 2009)
        'CCNB1',    # Cyclin B1, CDK1 partner for M-phase entry
    ],

    # --- Caspases ---
    # Laimins 2009: low-level caspase activation cleaves E1 for E2-independent replication
    'Caspase': [
        'CASP3',    # Caspase 3, cleaves E1 (Laimins 2009)
        'CASP7',    # Caspase 7, consensus cleavage motifs in E1 (Laimins 2009)
    ],

    # --- Immune signaling ---
    # E5/E7 downregulate MHC class I; STAT1 down while STAT5 up
    'Immune': [
        'STAT1',    # Predict DOWN in infected populations (immune evasion)
        'HLA-A',    # MHC class I, downregulated by E5/E7
        'HLA-B',    # MHC class I
        'HLA-C',    # MHC class I
        'IRF1',     # Interferon regulatory factor
        'TAP1',     # Antigen processing
        'B2M',      # Beta-2-microglobulin, MHC class I component
    ],

    # --- Differentiation markers ---
    # Productive infection tied to differentiation; relevant to cell naming
    'Differentiation': [
        'KRT5',     # Basal keratin
        'KRT14',    # Basal keratin
        'KRT1',     # Suprabasal differentiation
        'KRT10',    # Suprabasal differentiation
        'CDH1',     # E-cadherin
        'IVL',      # Involucrin, late differentiation
    ],

    # --- Innate immune / APOBEC regulation ---
    # A3A is interferon-responsive; cGAS-STING may sense integration
    'Innate_APOBEC': [
        'APOBEC3A',
        'APOBEC3B',
        'APOBEC3C',
        'APOBEC3D',
        'APOBEC3F',
        'APOBEC3G',
        'APOBEC3H',
        'CGAS',      # cGAS [HUGO: CGAS, previously MB21D1]
        'STING1',   # STING [HUGO: STING1, previously TMEM173]
    ],
}

# Population display order and colors (for consistency with Figure 4)
POP_ORDER = ['SBS2_HIGH', 'CNV_HIGH', 'NORMAL']
POP_LABELS = {'SBS2_HIGH': 'SBS2-HIGH', 'CNV_HIGH': 'CNV-HIGH', 'NORMAL': 'Normal'}

# =============================================================================
# LOGGING
# =============================================================================
report_lines = []

def log(msg=""):
    print(msg, flush=True)
    report_lines.append(msg)

def log_sep(title=""):
    log("")
    log("=" * 90)
    if title:
        log(f"  {title}")
        log("=" * 90)

# =============================================================================
# HELPER: Extract gene expression from adata
# =============================================================================
def get_expression(adata, gene_symbol):
    """Extract expression vector for a gene. Returns None if not found."""
    if gene_symbol in adata.var_names:
        idx = adata.var_names.get_loc(gene_symbol)
        x = adata.X[:, idx]
        if scipy.sparse.issparse(x):
            return np.asarray(x.todense()).flatten()
        return np.asarray(x).flatten()
    # Try gene_symbol column if var_names are Ensembl IDs
    if 'gene_symbol' in adata.var.columns:
        mask = adata.var['gene_symbol'] == gene_symbol
        if mask.any():
            idx = np.where(mask)[0][0]
            x = adata.X[:, idx]
            if scipy.sparse.issparse(x):
                return np.asarray(x.todense()).flatten()
            return np.asarray(x).flatten()
    return None

# =============================================================================
# HELPER: Mann-Whitney with direction
# =============================================================================
def mwu_test(vals_a, vals_b, label_a, label_b):
    """Run Mann-Whitney U test, return dict with stats."""
    vals_a = vals_a.dropna() if hasattr(vals_a, 'dropna') else vals_a[~np.isnan(vals_a)]
    vals_b = vals_b.dropna() if hasattr(vals_b, 'dropna') else vals_b[~np.isnan(vals_b)]
    if len(vals_a) < 5 or len(vals_b) < 5:
        return None
    stat, p = mannwhitneyu(vals_a, vals_b, alternative='two-sided')
    direction = f"{label_a} > {label_b}" if np.mean(vals_a) > np.mean(vals_b) else f"{label_b} > {label_a}"
    return {'U': stat, 'p': p, 'direction': direction,
            f'mean_{label_a}': np.mean(vals_a), f'mean_{label_b}': np.mean(vals_b)}


# =============================================================================
# STEP 0: LOAD DATA
# =============================================================================
log_sep("STEP 0: Load data")

# --- Population assignments ---
groups = pd.read_csv(THREE_GROUP_PATH, sep='\t')
log(f"  Three-group assignments: {groups.shape}")
log(f"  Columns: {list(groups.columns)}")
log(f"  Group counts:")
for grp, n in groups['group'].value_counts().items():
    log(f"    {grp:15s}: {n}")

sbs2_cells = set(groups.loc[groups['group'] == 'SBS2_HIGH', 'cell_barcode'])
cnv_cells  = set(groups.loc[groups['group'] == 'CNV_HIGH',  'cell_barcode'])
normal_cells = set(groups.loc[groups['group'] == 'NORMAL',  'cell_barcode'])

log(f"\n  SBS2-HIGH: {len(sbs2_cells)}")
log(f"  CNV-HIGH:  {len(cnv_cells)}")
log(f"  NORMAL:    {len(normal_cells)}")

# --- adata ---
log(f"\n  Loading adata_final.h5ad...")
adata = sc.read_h5ad(ADATA_PATH)
log(f"  adata: {adata.shape[0]} cells x {adata.shape[1]} genes")

# Verify overlap
all_pop_cells = sbs2_cells | cnv_cells | normal_cells
overlap = all_pop_cells & set(adata.obs_names)
log(f"  Population cells found in adata: {len(overlap)} / {len(all_pop_cells)}")

# Create population label in adata.obs
adata.obs['population'] = 'other'
adata.obs.loc[adata.obs_names.isin(sbs2_cells), 'population'] = 'SBS2_HIGH'
adata.obs.loc[adata.obs_names.isin(cnv_cells), 'population'] = 'CNV_HIGH'
adata.obs.loc[adata.obs_names.isin(normal_cells), 'population'] = 'NORMAL'

pop_check = adata.obs['population'].value_counts()
log(f"  Population labels in adata.obs:")
for pop, n in pop_check.items():
    log(f"    {pop:15s}: {n}")

# --- HPV16 gene counts (optional, may not exist yet) ---
hpv_available = os.path.exists(HPV_GENE_PATH)
if hpv_available:
    hpv_genes = pd.read_csv(HPV_GENE_PATH, sep='\t', index_col=0)
    log(f"\n  HPV16 gene counts: {hpv_genes.shape}")
    log(f"  Columns: {list(hpv_genes.columns)}")

    # Overlap with populations
    hpv_in_sbs2 = len(set(hpv_genes.index) & sbs2_cells)
    hpv_in_cnv  = len(set(hpv_genes.index) & cnv_cells)
    hpv_in_norm = len(set(hpv_genes.index) & normal_cells)
    log(f"  HPV16 genome data overlap:")
    log(f"    SBS2-HIGH: {hpv_in_sbs2}")
    log(f"    CNV-HIGH:  {hpv_in_cnv}")
    log(f"    NORMAL:    {hpv_in_norm}")
else:
    log(f"\n  WARNING: HPV16 gene count file not found: {HPV_GENE_PATH}")
    log(f"  Diagnostics 1 and 2 will be skipped.")
    hpv_genes = None


# =============================================================================
# DIAGNOSTIC 1: HPV16 GENE EXPRESSION BY FOUR-PHASE GROUPING
# =============================================================================
log_sep("DIAGNOSTIC 1: HPV16 gene expression by lifecycle phase")

if hpv_genes is not None:
    # Identify HPV16 gene columns present in the data
    all_hpv_gene_names = []
    for phase, genes in HPV16_PHASES.items():
        for g in genes:
            if g in hpv_genes.columns:
                all_hpv_gene_names.append(g)
            else:
                log(f"  WARNING: {g} not found in HPV16 gene count columns")

    log(f"  HPV16 genes found: {all_hpv_gene_names}")

    # Tag each HPV cell with its population
    hpv_genes['population'] = 'other'
    hpv_genes.loc[hpv_genes.index.isin(sbs2_cells), 'population'] = 'SBS2_HIGH'
    hpv_genes.loc[hpv_genes.index.isin(cnv_cells), 'population'] = 'CNV_HIGH'
    hpv_genes.loc[hpv_genes.index.isin(normal_cells), 'population'] = 'NORMAL'

    # --- Per-gene, per-population summary ---
    log(f"\n  Mean reads per cell by HPV16 gene and population:")
    log(f"  {'Phase':<15s} {'Gene':<6s} {'SBS2-HIGH':>12s} {'CNV-HIGH':>12s} {'NORMAL':>12s}  "
        f"{'KW p':>10s} {'HIGH vs CNV p':>14s} {'Direction':>18s}")
    log(f"  {'-'*15} {'-'*6} {'-'*12} {'-'*12} {'-'*12}  {'-'*10} {'-'*14} {'-'*18}")

    phase_summary_rows = []

    for phase, genes in HPV16_PHASES.items():
        for gene in genes:
            if gene not in hpv_genes.columns:
                continue

            row = {'phase': phase, 'gene': gene}
            means = {}
            for pop in POP_ORDER:
                vals = hpv_genes.loc[hpv_genes['population'] == pop, gene]
                means[pop] = vals.mean()
                row[f'mean_{pop}'] = vals.mean()
                row[f'median_{pop}'] = vals.median()
                row[f'pct_nonzero_{pop}'] = 100 * (vals > 0).sum() / max(len(vals), 1)
                row[f'n_{pop}'] = len(vals)

            # Kruskal-Wallis across all three
            groups_vals = [hpv_genes.loc[hpv_genes['population'] == pop, gene].values
                          for pop in POP_ORDER
                          if (hpv_genes['population'] == pop).sum() > 0]
            if len(groups_vals) == 3 and all(len(g) > 5 for g in groups_vals):
                kw_stat, kw_p = kruskal(*groups_vals)
                row['kw_p'] = kw_p
            else:
                kw_p = np.nan
                row['kw_p'] = np.nan

            # Pairwise: SBS2-HIGH vs CNV-HIGH
            v_sbs2 = hpv_genes.loc[hpv_genes['population'] == 'SBS2_HIGH', gene].values
            v_cnv  = hpv_genes.loc[hpv_genes['population'] == 'CNV_HIGH', gene].values
            if len(v_sbs2) > 5 and len(v_cnv) > 5:
                u, p = mannwhitneyu(v_sbs2, v_cnv, alternative='two-sided')
                direction = "SBS2 > CNV" if means['SBS2_HIGH'] > means['CNV_HIGH'] else "CNV > SBS2"
                row['hvc_p'] = p
                row['hvc_direction'] = direction
            else:
                p = np.nan
                direction = "N/A"
                row['hvc_p'] = np.nan
                row['hvc_direction'] = "N/A"

            log(f"  {phase:<15s} {gene:<6s} {means.get('SBS2_HIGH', 0):12.4f} "
                f"{means.get('CNV_HIGH', 0):12.4f} {means.get('NORMAL', 0):12.4f}  "
                f"{kw_p:10.2e} {p:14.2e} {direction:>18s}")

            phase_summary_rows.append(row)

    phase_df = pd.DataFrame(phase_summary_rows)
    phase_path = os.path.join(OUTPUT_DIR, "hpv16_gene_by_phase_and_population.tsv")
    phase_df.to_csv(phase_path, sep='\t', index=False)
    log(f"\n  Saved: {phase_path}")

    # --- Phase-level aggregates ---
    log(f"\n  Aggregate reads by lifecycle phase:")
    log(f"  {'Phase':<15s} {'SBS2-HIGH':>12s} {'CNV-HIGH':>12s} {'NORMAL':>12s}  {'KW p':>10s}")
    log(f"  {'-'*15} {'-'*12} {'-'*12} {'-'*12}  {'-'*10}")

    for phase, genes in HPV16_PHASES.items():
        available_genes = [g for g in genes if g in hpv_genes.columns]
        if not available_genes:
            continue

        for pop in POP_ORDER:
            mask = hpv_genes['population'] == pop
            hpv_genes.loc[mask, f'phase_{phase}'] = hpv_genes.loc[mask, available_genes].sum(axis=1)

        means = {}
        groups_vals = []
        for pop in POP_ORDER:
            vals = hpv_genes.loc[hpv_genes['population'] == pop, f'phase_{phase}'].values
            means[pop] = np.mean(vals)
            groups_vals.append(vals)

        if all(len(g) > 5 for g in groups_vals):
            _, kw_p = kruskal(*groups_vals)
        else:
            kw_p = np.nan

        log(f"  {phase:<15s} {means.get('SBS2_HIGH', 0):12.4f} "
            f"{means.get('CNV_HIGH', 0):12.4f} {means.get('NORMAL', 0):12.4f}  {kw_p:10.2e}")

    # =========================================================================
    # DIAGNOSTIC 1B: NORMALIZED HPV16 GENE FRACTIONS
    # =========================================================================
    log(f"\n  --- NORMALIZED: Gene reads as FRACTION of total HPV16 reads per cell ---")
    log(f"  (Controls for 2x higher viral load in CNV-HIGH)")

    # Compute per-gene fractions (only for cells with any HPV reads)
    total_col = 'total_hpv16_genome_reads' if 'total_hpv16_genome_reads' in hpv_genes.columns else None
    if total_col:
        for gene in all_hpv_gene_names:
            hpv_genes[f'{gene}_frac'] = hpv_genes[gene] / (hpv_genes[total_col] + 0.5)

        # Also compute phase-level fractions
        for phase, genes in HPV16_PHASES.items():
            available_genes = [g for g in genes if g in hpv_genes.columns]
            if available_genes:
                hpv_genes[f'phase_{phase}_frac'] = hpv_genes[available_genes].sum(axis=1) / (hpv_genes[total_col] + 0.5)

        # Report normalized per-gene fractions
        log(f"\n  {'Phase':<15s} {'Gene':<6s} {'SBS2-HIGH':>12s} {'CNV-HIGH':>12s} {'NORMAL':>12s}  "
            f"{'KW p':>10s} {'HvC p':>14s} {'Direction':>18s}")
        log(f"  {'-'*15} {'-'*6} {'-'*12} {'-'*12} {'-'*12}  {'-'*10} {'-'*14} {'-'*18}")

        norm_rows = []
        for phase, genes in HPV16_PHASES.items():
            for gene in genes:
                frac_col = f'{gene}_frac'
                if frac_col not in hpv_genes.columns:
                    continue

                means = {}
                groups_vals = []
                for pop in POP_ORDER:
                    mask = hpv_genes['population'] == pop
                    vals = hpv_genes.loc[mask, frac_col].dropna()
                    means[pop] = vals.mean() if len(vals) > 0 else 0
                    groups_vals.append(vals.values)

                valid = [g for g in groups_vals if len(g) > 5]
                kw_p = np.nan
                if len(valid) == 3:
                    _, kw_p = kruskal(*valid)

                hvc_p = np.nan
                direction = "N/A"
                if len(groups_vals[0]) > 5 and len(groups_vals[1]) > 5:
                    _, hvc_p = mannwhitneyu(groups_vals[0], groups_vals[1], alternative='two-sided')
                    direction = "SBS2 > CNV" if means['SBS2_HIGH'] > means['CNV_HIGH'] else "CNV > SBS2"

                log(f"  {phase:<15s} {gene:<6s} {means.get('SBS2_HIGH', 0):12.6f} "
                    f"{means.get('CNV_HIGH', 0):12.6f} {means.get('NORMAL', 0):12.6f}  "
                    f"{kw_p:10.2e} {hvc_p:14.2e} {direction:>18s}")

                norm_rows.append({
                    'phase': phase, 'gene': gene, 'metric': 'fraction_of_total',
                    'mean_SBS2_HIGH': means.get('SBS2_HIGH', np.nan),
                    'mean_CNV_HIGH': means.get('CNV_HIGH', np.nan),
                    'mean_NORMAL': means.get('NORMAL', np.nan),
                    'kw_p': kw_p, 'hvc_p': hvc_p, 'direction': direction,
                })

        # Phase-level normalized fractions
        log(f"\n  Normalized phase fractions:")
        log(f"  {'Phase':<15s} {'SBS2-HIGH':>12s} {'CNV-HIGH':>12s} {'NORMAL':>12s}  {'KW p':>10s} {'HvC p':>14s} {'Direction':>18s}")
        log(f"  {'-'*15} {'-'*12} {'-'*12} {'-'*12}  {'-'*10} {'-'*14} {'-'*18}")

        for phase in HPV16_PHASES:
            frac_col = f'phase_{phase}_frac'
            if frac_col not in hpv_genes.columns:
                continue

            means = {}
            groups_vals = []
            for pop in POP_ORDER:
                mask = hpv_genes['population'] == pop
                vals = hpv_genes.loc[mask, frac_col].dropna().values
                means[pop] = np.mean(vals) if len(vals) > 0 else 0
                groups_vals.append(vals)

            valid = [g for g in groups_vals if len(g) > 5]
            kw_p = np.nan
            if len(valid) == 3:
                _, kw_p = kruskal(*valid)

            hvc_p = np.nan
            direction = "N/A"
            if len(groups_vals[0]) > 5 and len(groups_vals[1]) > 5:
                _, hvc_p = mannwhitneyu(groups_vals[0], groups_vals[1], alternative='two-sided')
                direction = "SBS2 > CNV" if means['SBS2_HIGH'] > means['CNV_HIGH'] else "CNV > SBS2"

            log(f"  {phase:<15s} {means.get('SBS2_HIGH', 0):12.6f} "
                f"{means.get('CNV_HIGH', 0):12.6f} {means.get('NORMAL', 0):12.6f}  "
                f"{kw_p:10.2e} {hvc_p:14.2e} {direction:>18s}")

        # Save normalized fractions
        if norm_rows:
            norm_df = pd.DataFrame(norm_rows)
            norm_path = os.path.join(OUTPUT_DIR, "hpv16_gene_fractions_normalized.tsv")
            norm_df.to_csv(norm_path, sep='\t', index=False)
            log(f"\n  Saved: {norm_path}")

    else:
        log(f"  WARNING: total_hpv16_genome_reads column not found, cannot normalize")


# =============================================================================
# DIAGNOSTIC 2: INTEGRATION PROXY METRICS
# =============================================================================
log_sep("DIAGNOSTIC 2: Integration proxy — E2:(E6+E7) ratio")

if hpv_genes is not None:
    # Compute integration proxy metrics
    has_e2 = 'E2' in hpv_genes.columns
    has_e6 = 'E6' in hpv_genes.columns
    has_e7 = 'E7' in hpv_genes.columns

    if has_e2 and has_e6 and has_e7:
        # E2:(E6+E7) ratio — add pseudocount to avoid division by zero
        hpv_genes['E6E7_sum'] = hpv_genes['E6'] + hpv_genes['E7']
        hpv_genes['E2_to_E6E7'] = hpv_genes['E2'] / (hpv_genes['E6E7_sum'] + 0.5)

        # Also compute E2 fraction of total early reads
        early_cols = [c for c in ['E1', 'E2', 'E4', 'E5', 'E6', 'E7'] if c in hpv_genes.columns]
        hpv_genes['total_early'] = hpv_genes[early_cols].sum(axis=1)
        hpv_genes['E2_fraction_of_early'] = hpv_genes['E2'] / (hpv_genes['total_early'] + 0.5)

        # Total HPV16 reads (episomal copy number proxy)
        total_col = 'total_hpv16_genome_reads' if 'total_hpv16_genome_reads' in hpv_genes.columns else None

        log(f"\n  Integration proxy metrics (all cells with HPV16 alignment data):")
        log(f"\n  {'Metric':<30s} {'SBS2-HIGH':>12s} {'CNV-HIGH':>12s} {'NORMAL':>12s}  "
            f"{'KW p':>10s} {'HvC p':>10s} {'Direction':>15s}")
        log(f"  {'-'*30} {'-'*12} {'-'*12} {'-'*12}  {'-'*10} {'-'*10} {'-'*15}")

        proxy_metrics = ['E2_to_E6E7', 'E2_fraction_of_early', 'E6E7_sum', 'E2']
        if total_col:
            proxy_metrics.append(total_col)
            # Add normalized fractions (gene reads / total HPV reads)
            hpv_genes['E6E7_frac_of_total'] = hpv_genes['E6E7_sum'] / (hpv_genes[total_col] + 0.5)
            hpv_genes['E2_frac_of_total'] = hpv_genes['E2'] / (hpv_genes[total_col] + 0.5)
            hpv_genes['L1L2_frac_of_total'] = (hpv_genes['L1'].fillna(0) + hpv_genes['L2'].fillna(0)) / (hpv_genes[total_col] + 0.5)
            proxy_metrics.extend(['E6E7_frac_of_total', 'E2_frac_of_total', 'L1L2_frac_of_total'])
        if 'early_late_ratio' in hpv_genes.columns:
            proxy_metrics.append('early_late_ratio')

        proxy_rows = []
        for metric in proxy_metrics:
            if metric not in hpv_genes.columns:
                continue

            means = {}
            groups_vals = []
            for pop in POP_ORDER:
                vals = hpv_genes.loc[hpv_genes['population'] == pop, metric].dropna()
                means[pop] = vals.mean() if len(vals) > 0 else np.nan
                groups_vals.append(vals.values)

            # Kruskal-Wallis
            valid_groups = [g for g in groups_vals if len(g) > 5]
            kw_p = np.nan
            if len(valid_groups) == 3:
                _, kw_p = kruskal(*valid_groups)

            # SBS2-HIGH vs CNV-HIGH
            hvc_p = np.nan
            direction = "N/A"
            if len(groups_vals[0]) > 5 and len(groups_vals[1]) > 5:
                _, hvc_p = mannwhitneyu(groups_vals[0], groups_vals[1], alternative='two-sided')
                direction = "SBS2 > CNV" if means['SBS2_HIGH'] > means['CNV_HIGH'] else "CNV > SBS2"

            log(f"  {metric:<30s} {means.get('SBS2_HIGH', 0):12.4f} "
                f"{means.get('CNV_HIGH', 0):12.4f} {means.get('NORMAL', 0):12.4f}  "
                f"{kw_p:10.2e} {hvc_p:10.2e} {direction:>15s}")

            proxy_rows.append({
                'metric': metric,
                'mean_SBS2_HIGH': means.get('SBS2_HIGH', np.nan),
                'mean_CNV_HIGH': means.get('CNV_HIGH', np.nan),
                'mean_NORMAL': means.get('NORMAL', np.nan),
                'kw_p': kw_p, 'hvc_p': hvc_p, 'direction': direction,
            })

        proxy_df = pd.DataFrame(proxy_rows)
        proxy_path = os.path.join(OUTPUT_DIR, "integration_proxy_metrics.tsv")
        proxy_df.to_csv(proxy_path, sep='\t', index=False)
        log(f"\n  Saved: {proxy_path}")

        # --- Detailed distribution stats for E2:(E6+E7) ---
        log(f"\n  E2:(E6+E7) ratio distribution detail:")
        for pop in POP_ORDER:
            vals = hpv_genes.loc[hpv_genes['population'] == pop, 'E2_to_E6E7'].dropna()
            if len(vals) > 0:
                log(f"    {POP_LABELS[pop]:12s}: n={len(vals):4d}, "
                    f"mean={vals.mean():.4f}, median={vals.median():.4f}, "
                    f"std={vals.std():.4f}, "
                    f"Q25={vals.quantile(0.25):.4f}, Q75={vals.quantile(0.75):.4f}")

        # --- Cells with ANY HPV16 reads by population ---
        log(f"\n  Cells with any HPV16 genome alignment reads:")
        if total_col:
            for pop in POP_ORDER:
                mask = hpv_genes['population'] == pop
                n_total = mask.sum()
                n_with_reads = (hpv_genes.loc[mask, total_col] > 0).sum()
                log(f"    {POP_LABELS[pop]:12s}: {n_with_reads} / {n_total} "
                    f"({100*n_with_reads/max(n_total,1):.1f}%)")

    else:
        log(f"  WARNING: Missing E2/E6/E7 columns in HPV gene data. Cannot compute integration proxy.")


# =============================================================================
# DIAGNOSTIC 3: HOST MARKER GENE EXPRESSION PANEL
# =============================================================================
log_sep("DIAGNOSTIC 3: Host marker gene expression panel")

# Subset adata to our three populations
pop_mask = adata.obs['population'].isin(POP_ORDER)
adata_pop = adata[pop_mask].copy()
log(f"  Cells in three populations: {adata_pop.shape[0]}")

# Scan which genes are available
all_marker_genes = []
for category, genes in MARKER_GENES.items():
    for gene in genes:
        all_marker_genes.append((category, gene))

log(f"\n  Checking {len(all_marker_genes)} marker genes in adata...")

found_genes = []
missing_genes = []
for category, gene in all_marker_genes:
    expr = get_expression(adata_pop, gene)
    if expr is not None:
        found_genes.append((category, gene))
    else:
        missing_genes.append((category, gene))

log(f"  Found: {len(found_genes)}")
log(f"  Missing: {len(missing_genes)}")
if missing_genes:
    for cat, gene in missing_genes:
        log(f"    NOT FOUND: {gene} ({cat})")

# --- Compute expression summary per population ---
log(f"\n  Computing expression statistics per population...")

summary_rows = []
per_cell_rows = []

log(f"\n  {'Category':<20s} {'Gene':<12s} "
    f"{'SBS2 mean':>10s} {'CNV mean':>10s} {'NORM mean':>10s}  "
    f"{'SBS2 %+':>8s} {'CNV %+':>8s} {'NORM %+':>8s}  "
    f"{'KW p':>10s} {'HvC p':>10s} {'Direction':>15s}")
log(f"  {'-'*20} {'-'*12} "
    f"{'-'*10} {'-'*10} {'-'*10}  "
    f"{'-'*8} {'-'*8} {'-'*8}  "
    f"{'-'*10} {'-'*10} {'-'*15}")

for category, gene in found_genes:
    expr_all = get_expression(adata_pop, gene)
    if expr_all is None:
        continue

    row = {'category': category, 'gene': gene}
    means = {}
    pct_pos = {}
    groups_vals = []

    for pop in POP_ORDER:
        mask = adata_pop.obs['population'] == pop
        vals = expr_all[mask.values]
        means[pop] = np.mean(vals)
        pct_pos[pop] = 100 * np.sum(vals > 0) / max(len(vals), 1)
        groups_vals.append(vals)

        row[f'mean_{pop}'] = np.mean(vals)
        row[f'median_{pop}'] = np.median(vals)
        row[f'pct_nonzero_{pop}'] = pct_pos[pop]
        row[f'n_{pop}'] = len(vals)

        # Store per-cell values for downstream use
        cell_barcodes = adata_pop.obs_names[mask.values]
        for bc, v in zip(cell_barcodes, vals):
            per_cell_rows.append({
                'cell_barcode': bc, 'population': pop,
                'gene': gene, 'category': category, 'expression': v
            })

    # Kruskal-Wallis
    valid = [g for g in groups_vals if len(g) > 5]
    kw_p = np.nan
    if len(valid) == 3:
        _, kw_p = kruskal(*valid)
    row['kw_p'] = kw_p

    # SBS2-HIGH vs CNV-HIGH
    hvc_p = np.nan
    direction = "N/A"
    if len(groups_vals[0]) > 5 and len(groups_vals[1]) > 5:
        _, hvc_p = mannwhitneyu(groups_vals[0], groups_vals[1], alternative='two-sided')
        direction = "SBS2 > CNV" if means['SBS2_HIGH'] > means['CNV_HIGH'] else "CNV > SBS2"
    row['hvc_p'] = hvc_p
    row['hvc_direction'] = direction

    # SBS2-HIGH vs NORMAL
    hvn_p = np.nan
    hvn_dir = "N/A"
    if len(groups_vals[0]) > 5 and len(groups_vals[2]) > 5:
        _, hvn_p = mannwhitneyu(groups_vals[0], groups_vals[2], alternative='two-sided')
        hvn_dir = "SBS2 > NORM" if means['SBS2_HIGH'] > means['NORMAL'] else "NORM > SBS2"
    row['hvn_p'] = hvn_p
    row['hvn_direction'] = hvn_dir

    # CNV-HIGH vs NORMAL
    cvn_p = np.nan
    cvn_dir = "N/A"
    if len(groups_vals[1]) > 5 and len(groups_vals[2]) > 5:
        _, cvn_p = mannwhitneyu(groups_vals[1], groups_vals[2], alternative='two-sided')
        cvn_dir = "CNV > NORM" if means['CNV_HIGH'] > means['NORMAL'] else "NORM > CNV"
    row['cvn_p'] = cvn_p
    row['cvn_direction'] = cvn_dir

    summary_rows.append(row)

    log(f"  {category:<20s} {gene:<12s} "
        f"{means['SBS2_HIGH']:10.4f} {means['CNV_HIGH']:10.4f} {means['NORMAL']:10.4f}  "
        f"{pct_pos['SBS2_HIGH']:7.1f}% {pct_pos['CNV_HIGH']:7.1f}% {pct_pos['NORMAL']:7.1f}%  "
        f"{kw_p:10.2e} {hvc_p:10.2e} {direction:>15s}")

# --- Save summary ---
summary_df = pd.DataFrame(summary_rows)
summary_path = os.path.join(OUTPUT_DIR, "host_marker_expression_summary.tsv")
summary_df.to_csv(summary_path, sep='\t', index=False)
log(f"\n  Saved: {summary_path}")

# --- Save per-cell values (for downstream plotting) ---
if per_cell_rows:
    per_cell_df = pd.DataFrame(per_cell_rows)
    per_cell_path = os.path.join(OUTPUT_DIR, "host_marker_per_cell_values.tsv")
    per_cell_df.to_csv(per_cell_path, sep='\t', index=False)
    log(f"  Saved: {per_cell_path} ({len(per_cell_df):,} rows)")

# =============================================================================
# SUMMARY: PREDICTION SCORECARD
# =============================================================================
log_sep("PREDICTION SCORECARD")

log(f"""
  Based on our HPV lifecycle model, we predicted:

  INTEGRATION INDICATORS (SBS2-HIGH > CNV-HIGH):
    - Low E2:(E6+E7) ratio in SBS2-HIGH (E2 disrupted by integration)
    - High CDKN2A (p16) in SBS2-HIGH (E7 -> Rb inactivation -> p16 feedback)
    - High CCNE1, E2F1, E2F2 in SBS2-HIGH (E2F targets from Rb inactivation)
    - Low CDKN1A (p21) in SBS2-HIGH (E6 -> p53 degradation)
    - High BRD4, MED1 in SBS2-HIGH (super-enhancer at integration site)
    - Lower total HPV16 reads in SBS2-HIGH (integrated, not amplified episome)

  PRODUCTIVE INFECTION INDICATORS (CNV-HIGH > SBS2-HIGH):
    - High L1, L2 in CNV-HIGH (capsid assembly, productive lifecycle)
    - Higher E2:(E6+E7) ratio in CNV-HIGH (E2 intact, episomal)
    - Higher total HPV16 reads in CNV-HIGH (episomal amplification)
    - G2/M arrest genes enriched in CNV-HIGH (CDC25C, CDK1, CCNB1)
    - ATM/MRN cascade active in CNV-HIGH (productive replication machinery)
    - Differentiation markers in CNV-HIGH (KRT1, KRT10)

  Check the summary tables above against these predictions.
""")

# Check predictions from summary_df if available
if len(summary_rows) > 0:
    log(f"  Checking key predictions against observed data:")
    log(f"  {'Gene':<12s} {'Predicted':>20s} {'Observed':>20s} {'p-value':>12s} {'Match':>6s}")
    log(f"  {'-'*12} {'-'*20} {'-'*20} {'-'*12} {'-'*6}")

    predictions = [
        # Integration/transformation (predicted SBS2 > CNV)
        ('CDKN2A',  'SBS2 > CNV'),
        ('CCNE1',   'SBS2 > CNV'),
        ('MKI67',   'SBS2 > CNV'),
        ('BRD4',    'SBS2 > CNV'),
        ('H2AX',    'SBS2 > CNV'),   # DNA damage
        # p53/Rb pathway
        ('CDKN1A',  'CNV > SBS2'),   # p21 LOW in SBS2 means CNV > SBS2
        ('MDM2',    'CNV > SBS2'),   # p53 regulator
        # Productive infection (predicted CNV > SBS2)
        ('CDC25C',  'CNV > SBS2'),
        ('CDK1',    'CNV > SBS2'),
        ('CCNB1',   'CNV > SBS2'),
        ('CHEK2',   'CNV > SBS2'),
        ('BRCA1',   'CNV > SBS2'),
        # Differentiation
        ('KRT1',    'CNV > SBS2'),
        ('KRT10',   'CNV > SBS2'),
        ('IVL',     'CNV > SBS2'),   # Differentiation marker
        # Immune (predicted DOWN in infected via E5/E7 evasion)
        ('HLA-A',   'NORM > SBS2'),
        ('HLA-B',   'NORM > SBS2'),
        ('B2M',     'NORM > SBS2'),
        ('STAT1',   'NORM > SBS2'),
        # APOBEC
        ('APOBEC3A','SBS2 > CNV'),
        ('APOBEC3B','CNV > SBS2'),
    ]

    for gene, predicted_dir in predictions:
        match_rows = [r for r in summary_rows if r['gene'] == gene]
        if match_rows:
            r = match_rows[0]
            observed = r['hvc_direction']
            p = r['hvc_p']
            match = "YES" if predicted_dir == observed else "NO"
            log(f"  {gene:<12s} {predicted_dir:>20s} {observed:>20s} {p:12.2e} {match:>6s}")
        else:
            log(f"  {gene:<12s} {predicted_dir:>20s} {'NOT FOUND':>20s} {'N/A':>12s} {'---':>6s}")


# =============================================================================
# SAVE REPORT
# =============================================================================
log_sep("COMPLETE")

report_path = os.path.join(OUTPUT_DIR, "diagnostic_lifecycle_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Report saved: {report_path}")
log(f"  Output directory: {OUTPUT_DIR}")
