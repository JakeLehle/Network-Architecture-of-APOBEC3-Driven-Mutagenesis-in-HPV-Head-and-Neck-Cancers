#!/usr/bin/env python3
"""
patient_config.py
=================

Shared configuration for Figure 5: Patient-Specific Effects analysis.
All scripts in PATIENT_SPECIFIC_EFFECTS/ import from this file.

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
GROUP_ASSIGN_PATH = os.path.join(BASE_DIR, "data", "FIG_4", "01_group_selection", "SC_Basal_group_assignments.tsv")
DE_GENES_DIR      = os.path.join(BASE_DIR, "data", "FIG_4", "02_differential_expression")
COMMUNITIES_DIR   = os.path.join(BASE_DIR, "data", "FIG_4", "04_communities")
CENTRALITY_DIR    = os.path.join(BASE_DIR, "data", "FIG_4", "05_centrality_metrics")

# Expression matrices (from Figure 4 Step 00)
HIGH_EXPR_PATH = os.path.join(BASE_DIR, "data", "FIG_4", "01_group_selection", "SC_Basal_SBS2_HIGH_expression.tsv")
LOW_EXPR_PATH  = os.path.join(BASE_DIR, "data", "FIG_4", "01_group_selection", "SC_Basal_SBS2_LOW_expression.tsv")

# SComatic mutation calls (from ClusterCatcher pipeline, interactive run)
SCOMATIC_PATH = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/all_samples.single_cell_genotype.filtered.tsv"

# SC network pipeline (for Tier 3 re-runs)
SC_NETWORK_SCRIPTS = os.path.join(BASE_DIR, "scripts", "NETWORK_SINGLE_CELL")

# Harris A3 interactor lists
HARRIS_ALL_PATH  = os.path.join(BASE_DIR, "data", "FIG_4", "00_input", "Harris_A3_interactors.txt")
HARRIS_A3B_PATH  = os.path.join(BASE_DIR, "data", "FIG_4", "00_input", "Harris_A3_interactors_A3B_only.txt")

# ENSG-to-symbol mapping
ENSG_TO_SYMBOL_PATH = os.path.join(BASE_DIR, "data", "FIG_2", "01_cleaned_expression", "ensg_to_symbol.json")

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
C0_INTERACTORS = ["HNRNPA2B1", "HSPD1", "RPL5", "TIMM8B"]

# =============================================================================
# ANALYSIS PARAMETERS
# =============================================================================
RANDOM_SEED       = 42
DIFF_THRESHOLD    = 0.40
LEIDEN_RESOLUTION = 1.0
GSEA_MIN_SIZE     = 10
GSEA_MAX_SIZE     = 500
GSEA_PERMUTATIONS = 1000

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
            log("  ERROR: SBS2 not in adata.obs and weights file not found"); sys.exit(1)
    return adata

def load_groups():
    groups = pd.read_csv(GROUP_ASSIGN_PATH, sep='\t')
    high_cells = set(groups.loc[groups['group'] == 'HIGH', 'cell_barcode'])
    low_cells  = set(groups.loc[groups['group'] == 'LOW', 'cell_barcode'])
    log(f"  Groups loaded: {len(high_cells)} HIGH / {len(low_cells)} LOW")
    return high_cells, low_cells

def get_gene_expression(adata, gene_symbol):
    import scipy.sparse
    if gene_symbol in adata.var_names:
        idx = adata.var_names.get_loc(gene_symbol)
        x = adata.X[:, idx]
        if scipy.sparse.issparse(x): return np.asarray(x.todense()).flatten()
        return np.asarray(x).flatten()
    if 'gene_symbol' in adata.var.columns:
        mask = adata.var['gene_symbol'] == gene_symbol
        if mask.any():
            idx = np.where(mask)[0][0]
            x = adata.X[:, idx]
            if scipy.sparse.issparse(x): return np.asarray(x.todense()).flatten()
            return np.asarray(x).flatten()
    return None
