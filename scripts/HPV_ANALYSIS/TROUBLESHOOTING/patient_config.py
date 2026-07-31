#!/usr/bin/env python3
"""
patient_config.py
=================

Shared configuration for Figure 5: Patient-Specific Effects analysis.
All scripts in PATIENT_SPECIFIC_EFFECTS/ import from this file.

Updated May 2026 to align with Figure 4 V4 three-group pipeline:
  - Groups loaded from three_group_assignments.tsv (SBS2_HIGH, CNV_HIGH, NORMAL)
  - LOPO targets SBS2_VS_NORMAL network (activating chain analysis)
  - Network parameters match V4 auto-selected values
  - Expression matrices read from NETWORK_SBS2_VS_NORMAL/ subdirectory

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime

# =============================================================================
# PROJECT PATHS
# =============================================================================
BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"

# Input data (from Figure 4 pipeline)
ADATA_PATH        = os.path.join(BASE_DIR, "data", "FIG_4", "00_input", "adata_final.h5ad")
WEIGHTS_PATH      = os.path.join(BASE_DIR, "data", "FIG_4", "00_input", "signature_weights_per_cell.txt")

# Three-group assignments (V4: SBS2_HIGH, CNV_HIGH, NORMAL -- 546 each)
THREE_GROUP_PATH  = os.path.join(BASE_DIR, "data", "FIG_4", "01_group_selection", "three_group_assignments.tsv")

# Per-network paths for SBS2_VS_NORMAL (the LOPO target network)
NETWORK_DIR       = os.path.join(BASE_DIR, "data", "FIG_4", "NETWORK_SBS2_VS_NORMAL")
DE_GENES_DIR      = os.path.join(NETWORK_DIR, "02_differential_expression")
CORR_DIR          = os.path.join(NETWORK_DIR, "03_correlation_networks")
COMMUNITIES_DIR   = os.path.join(NETWORK_DIR, "04_communities")
CENTRALITY_DIR    = os.path.join(NETWORK_DIR, "05_centrality_metrics")

# Expression matrices (from Step00B, SBS2_VS_NORMAL comparison)
EXPR_DIR          = os.path.join(BASE_DIR, "data", "FIG_4", "01_group_selection", "NETWORK_SBS2_VS_NORMAL")
HIGH_EXPR_PATH    = os.path.join(EXPR_DIR, "SC_Basal_SBS2_HIGH_expression.tsv")
NORMAL_EXPR_PATH  = os.path.join(EXPR_DIR, "SC_Basal_SBS2_LOW_expression.tsv")

# Legacy alias for scripts that still reference LOW_EXPR_PATH
LOW_EXPR_PATH     = NORMAL_EXPR_PATH

# SComatic mutation calls (from ClusterCatcher pipeline, interactive run)
SCOMATIC_PATH = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/all_samples.single_cell_genotype.filtered.tsv"

# SC network pipeline scripts (for reference; LOPO re-implements V4 internally)
SC_NETWORK_SCRIPTS = os.path.join(BASE_DIR, "scripts", "NETWORK_SINGLE_CELL")

# Harris A3 interactor lists
HARRIS_ALL_PATH  = os.path.join(BASE_DIR, "data", "FIG_4", "00_input", "Harris_A3_interactors.txt")
HARRIS_A3B_PATH  = os.path.join(BASE_DIR, "data", "FIG_4", "00_input", "Harris_A3_interactors_A3B_only.txt")

# ENSG-to-symbol mapping
ENSG_TO_SYMBOL_PATH = os.path.join(BASE_DIR, "data", "FIG_2", "01_cleaned_expression", "ensg_to_symbol.json")

# Full-analysis network results (for LOPO comparison baseline)
FULL_ANALYSIS_PARAMS = os.path.join(COMMUNITIES_DIR, "SC_selected_parameters.txt")
FULL_ANALYSIS_PARTITION = os.path.join(COMMUNITIES_DIR, "SC_best_partition.csv")

# =============================================================================
# OUTPUT DIRECTORIES
# =============================================================================
FIG5_ROOT          = os.path.join(BASE_DIR, "data", "FIG_5")
DIR_00_DIAG        = os.path.join(FIG5_ROOT, "00_diagnostics")
DIR_01_EXPRESSION  = os.path.join(FIG5_ROOT, "01_patient_expression")
DIR_02_SNP         = os.path.join(FIG5_ROOT, "02_snp_haplotype")
DIR_03_SENSITIVITY = os.path.join(FIG5_ROOT, "03_sensitivity")
FIGURE_5_PANELS    = os.path.join(FIG5_ROOT, "FIGURE_5_PANELS")

# =============================================================================
# ADATA COLUMN NAMES
# =============================================================================
PATIENT_COL      = "subject id"
TISSUE_COL       = "tissue type"
CELLTYPE_COL     = "final_annotation"
RUN_COL          = "run_accession"

# =============================================================================
# PATIENT DEFINITIONS (from v2 diagnostic)
# =============================================================================
HIGH_CONTRIBUTORS = ["Patient SC029", "Patient SC013", "Patient SC001"]
ALL_PATIENTS = [
    "Patient SC001", "Patient SC003", "Patient SC005", "Patient SC006",
    "Patient SC008", "Patient SC010", "Patient SC013", "Patient SC014",
    "Patient SC019", "Patient SC022", "Patient SC025", "Patient SC026",
    "Patient SC027", "Patient SC029",
]

# =============================================================================
# A3 GENE DEFINITIONS
# =============================================================================
A3_GENES_SYMBOLS = [
    "APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
    "APOBEC3F", "APOBEC3G", "APOBEC3H",
]
A3_LOCUS_CHR   = "chr22"
A3_LOCUS_START = 38950000
A3_LOCUS_END   = 39100000
# Legacy alias (kept for backward compatibility with Tier2B)
C0_INTERACTORS = ["HNRNPA2B1", "HSPD1", "RPL5", "TIMM8B"]

# Figure 4 concordance analysis gene sets (SBS2_VS_NORMAL)
A3_INTERACTOR_ANCHORS = ["RALY", "HNRNPA2B1"]
ACTIVATING_CHAIN_GENES = [
    "RALY", "HNRNPA2B1",           # Harris A3 interactor anchors
    "CCL20", "KRT24", "LCN2",      # Inflammatory/epithelial
    "LINC00278", "RRAD", "SMOX",   # Oxidative stress/signaling
    "UTY",                          # Chromatin
]
INHIBITING_CHAIN_ANCHORS = ["SNHG3", "THYN1", "ZNG1A"]  # A3B boundary genes

# =============================================================================
# NETWORK PARAMETERS (V4: matching SBS2_VS_NORMAL full analysis)
# =============================================================================
RANDOM_SEED       = 42

# DE parameters (V4: scanpy rank_genes_groups, FDR < 0.05)
DE_FDR_THRESHOLD  = 0.05
FORCE_KEEP_A3     = False    # A3 passes FDR naturally in SC

# DIFF threshold (V4 auto-selected for SBS2_VS_NORMAL: 0.40)
# LOPO uses the same threshold as the full analysis for direct comparison.
# If auto-selection is enabled, each LOPO run can re-select independently.
DIFF_THRESHOLD       = 0.40
DIFF_THRESHOLD_AUTO  = True
SWEEP_THRESHOLDS     = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65,
                        0.70, 0.75, 0.80, 0.85, 0.90]

# Leiden parameters (V4: full-network, component-aware merge)
LEIDEN_RESOLUTION       = 0.70
COMMUNITY_RESOLUTIONS   = [0.2, 0.4, 0.6, 0.8, 1.0]
RUNS_PER_RESOLUTION     = 15
COMMUNITY_BASE_SEED     = 42
USE_LARGEST_COMPONENT   = False   # Full-network Leiden (V4)
MIN_COMMUNITY_SIZE      = 10      # Satellites preserved below this
TARGET_BIG_COMMUNITIES  = 14

# GSEA parameters (Tier 1 analyses)
GSEA_MIN_SIZE     = 10
GSEA_MAX_SIZE     = 500
GSEA_PERMUTATIONS = 1000

# =============================================================================
# STYLE SETTINGS (consistent with Figure 4)
# =============================================================================
FONT_TITLE = 34
FONT_LABEL = 30
FONT_TICK  = 28

COLOR_SBS2_HIGH = "#ed6a5a"
COLOR_CNV_HIGH  = "#F6D155"
COLOR_NORMAL    = "#4682b4"
COLOR_OTHER     = "#cccccc"

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================
def log(msg):
    timestamp = datetime.now().strftime('%H:%M:%S')
    print(f"[{timestamp}] {msg}", flush=True)

def banner(title):
    log("")
    log("=" * 70)
    log(title)
    log("=" * 70)

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path

def load_adata():
    """Load ClusterCatcher AnnData with signature weights."""
    import scanpy as sc
    log(f"Loading: {ADATA_PATH}")
    adata = sc.read_h5ad(ADATA_PATH)
    log(f"  Total cells: {adata.n_obs:,} | Genes: {adata.n_vars:,}")
    if 'SBS2' not in adata.obs.columns:
        if os.path.exists(WEIGHTS_PATH):
            log(f"  Loading signature weights...")
            weights = pd.read_csv(WEIGHTS_PATH, sep='\t', index_col=0)
            for sig in weights.index:
                adata.obs[sig] = adata.obs.index.map(weights.loc[sig]).fillna(0)
            log(f"  Added {len(weights.index)} signatures")
        else:
            log("  ERROR: SBS2 not in adata.obs and weights file not found")
            sys.exit(1)
    return adata

def load_three_groups():
    """Load all three cell populations from the V4 three-group assignments.

    Returns:
        sbs2_high (set): SBS2-HIGH cell barcodes (n=546)
        cnv_high (set):  CNV-HIGH cell barcodes (n=546)
        normal (set):    NORMAL cell barcodes (n=546)
    """
    groups = pd.read_csv(THREE_GROUP_PATH, sep='\t')
    sbs2_high = set(groups.loc[groups['group'] == 'SBS2_HIGH', 'cell_barcode'])
    cnv_high  = set(groups.loc[groups['group'] == 'CNV_HIGH',  'cell_barcode'])
    normal    = set(groups.loc[groups['group'] == 'NORMAL',    'cell_barcode'])
    log(f"  Three groups loaded: {len(sbs2_high)} SBS2-HIGH / "
        f"{len(cnv_high)} CNV-HIGH / {len(normal)} NORMAL")
    return sbs2_high, cnv_high, normal

def load_groups():
    """Backward-compatible wrapper: returns (sbs2_high, normal) for
    SBS2_VS_NORMAL analyses. Scripts that need all three groups
    should call load_three_groups() directly."""
    sbs2_high, _cnv_high, normal = load_three_groups()
    log(f"  SBS2_VS_NORMAL groups: {len(sbs2_high)} HIGH / {len(normal)} NORMAL")
    return sbs2_high, normal

def get_gene_expression(adata, gene_symbol):
    """Extract expression vector for a single gene by symbol."""
    import scipy.sparse
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

def load_harris_interactors():
    """Load Harris A3 interactor gene lists."""
    harris_all = set()
    harris_a3b = set()
    if os.path.exists(HARRIS_ALL_PATH):
        with open(HARRIS_ALL_PATH) as f:
            harris_all = {line.strip() for line in f if line.strip()}
    if os.path.exists(HARRIS_A3B_PATH):
        with open(HARRIS_A3B_PATH) as f:
            harris_a3b = {line.strip() for line in f if line.strip()}
    log(f"  Harris interactors: {len(harris_all)} all, {len(harris_a3b)} A3B-specific")
    return harris_all, harris_a3b

def load_full_analysis_params():
    """Load the V4 full-analysis network parameters for LOPO comparison."""
    params = {}
    if os.path.exists(FULL_ANALYSIS_PARAMS):
        with open(FULL_ANALYSIS_PARAMS) as f:
            for line in f:
                line = line.strip()
                if '=' in line:
                    key, val = line.split('=', 1)
                    params[key.strip()] = val.strip()
        log(f"  Full analysis params: threshold={params.get('DIFF_THRESHOLD')}, "
            f"resolution={params.get('LEIDEN_RESOLUTION')}, "
            f"communities={params.get('N_COMMUNITIES')}, "
            f"genes={params.get('N_GENES')}")
    else:
        log(f"  WARNING: Full analysis params not found at {FULL_ANALYSIS_PARAMS}")
    return params
