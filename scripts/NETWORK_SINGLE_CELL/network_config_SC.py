#!/usr/bin/env python3
"""
network_config_SC.py
====================

Centralized configuration for the Figure 4 Single-Cell Network Analysis.
Mirrors network_config.py from Figure 2 (TCGA bulk) but adapted for
single-cell basal epithelial cell expression data.

Key differences from Figure 2 (TCGA bulk):
  - Input is single-cell expression (cells = samples)
  - No TCGA clinical metadata -- groups defined by SBS2 weight
  - Gene IDs are HYBRID: gene symbols for annotated genes (APOBEC3B, KRT5),
    ENSG IDs only for unannotated transcripts (NaN-replaced in adata.var)
  - The SC data from Cell Ranger already contains ~20K protein-coding genes,
    so NO protein-coding filter is needed (unlike TCGA bulk which starts at ~60K)
  - DIFF threshold auto-selected using unified criterion (A3 connectivity +
    component peak), identical to TCGA bulk
  - Leiden resolution auto-selected using composite score (modularity x ARI x
    evenness), identical to TCGA bulk Step05
  - Output goes to data/FIG_4/

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os

# =============================================================================
# DIRECTORY PATHS
# =============================================================================
BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"

SC_DATA_DIR = os.path.join(BASE_DIR, "data", "FIG_4", "00_input")

FIG4_ROOT = os.path.join(BASE_DIR, "data", "FIG_4")

DIR_00_INPUT       = os.path.join(FIG4_ROOT, "00_input")
DIR_01_GROUPS      = os.path.join(FIG4_ROOT, "01_group_selection")
DIR_02_DE          = os.path.join(FIG4_ROOT, "02_differential_expression")
DIR_03_NETWORKS    = os.path.join(FIG4_ROOT, "03_correlation_networks")
DIR_04_COMMUNITIES = os.path.join(FIG4_ROOT, "04_communities")
DIR_05_CENTRALITY  = os.path.join(FIG4_ROOT, "05_centrality_metrics")
DIR_06_OVERLAP     = os.path.join(FIG4_ROOT, "06_overlap_analysis")
FIGURE_4_PANELS    = os.path.join(FIG4_ROOT, "FIGURE_4_PANELS")

# Figure 2 data (for cross-reference in overlap analysis)
FIG2_ROOT          = os.path.join(BASE_DIR, "data", "FIG_2")
FIG2_COMMUNITIES   = os.path.join(FIG2_ROOT, "05_communities")
FIG2_CLEANED       = os.path.join(FIG2_ROOT, "01_cleaned_expression")

ENSG_TO_SYMBOL_PATH = os.path.join(FIG2_CLEANED, "ensg_to_symbol.json")
ENSG_TO_BIOTYPE_PATH = os.path.join(FIG2_CLEANED, "ensg_to_biotype.json")

HARRIS_ALL_PATH    = os.path.join(DIR_00_INPUT, "Harris_A3_interactors.txt")
HARRIS_A3B_PATH    = os.path.join(DIR_00_INPUT, "Harris_A3_interactors_A3B_only.txt")

# =============================================================================
# INPUT FILES
# =============================================================================
ADATA_FINAL_PATH = os.path.join(DIR_00_INPUT, "adata_final.h5ad")
WEIGHTS_PATH     = os.path.join(DIR_00_INPUT, "signature_weights_per_cell.txt")

HIGH_EXPR_PATH = os.path.join(DIR_01_GROUPS, "SC_Basal_SBS2_HIGH_expression.tsv")
LOW_EXPR_PATH  = os.path.join(DIR_01_GROUPS, "SC_Basal_SBS2_LOW_expression.tsv")
GROUP_ASSIGN_PATH = os.path.join(DIR_01_GROUPS, "SC_Basal_group_assignments.tsv")

# =============================================================================
# STEP 00 -- CELL SELECTION PARAMETERS
# =============================================================================
TARGET_CELL_TYPE = "basal cell"
SBS2_HIGH_PERCENTILE = 0.80
N_CONTROLS_MATCH     = True
SBS2_LOW_MAX_WEIGHT  = 0.01

# =============================================================================
# STEP 01 -- DIFFERENTIAL EXPRESSION PARAMETERS
# =============================================================================
MIN_CELLS_DETECTED  = 10
RAW_P_THRESHOLD     = 0.05
LOGFC_THRESHOLD     = 0
FORCE_KEEP_A3       = False      # Let A3 genes pass/fail DE on their own merit
FILTER_PROTEIN_CODING = False    # SC data already ~20K protein-coding genes

# =============================================================================
# STEP 02 -- CORRELATION NETWORK PARAMETERS
# =============================================================================
CORRELATION_METHOD = "spearman"

CORR_THRESHOLD     = 0.80     # |rho| threshold for HIGH/LOW network edges

# DIFF threshold: auto-selected using unified criterion.
# Among thresholds where A3A and/or A3B (those passing DE) have degree >= 1,
# select the threshold with the highest connected component count.
# Ties broken by the lower (more inclusive) threshold.
# This fallback value is used only if auto-selection fails.
DIFF_THRESHOLD     = 0.40

# Threshold sweep range for auto-selection
SWEEP_THRESHOLDS = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70,
                    0.75, 0.80, 0.85, 0.90]

# =============================================================================
# STEP 03 -- COMMUNITY DETECTION PARAMETERS
# =============================================================================
COMMUNITY_METHOD    = "leiden"

# Resolution sweep: 0.1 to 0.8 in 0.1 steps (matches TCGA Step05)
# Selection: composite score = modularity x ARI x evenness
COMMUNITY_RESOLUTIONS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
RUNS_PER_RESOLUTION = 15
COMMUNITY_BASE_SEED = 42
USE_LARGEST_COMPONENT = True
TARGET_BIG_COMMUNITIES = 14
MIN_COMMUNITY_SIZE     = 10

# =============================================================================
# STEP 05 -- OVERLAP ANALYSIS PARAMETERS
# =============================================================================
OVERLAP_FDR_METHOD  = "BH"
OVERLAP_P_THRESHOLD = 0.05
JACCARD_DISPLAY     = True

# =============================================================================
# APOBEC3 GENE DEFINITIONS
# =============================================================================

# PRIMARY: Gene symbols (matches SC expression matrix index)
A3_GENES_SYMBOLS = [
    "APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
    "APOBEC3F", "APOBEC3G", "APOBEC3H",
]

A3_SYMBOL_TO_ALIAS = {
    "APOBEC3A": "A3A", "APOBEC3B": "A3B", "APOBEC3C": "A3C",
    "APOBEC3D": "A3D", "APOBEC3F": "A3F", "APOBEC3G": "A3G",
    "APOBEC3H": "A3H",
}

A3_ALIAS_TO_SYMBOL = {v: k for k, v in A3_SYMBOL_TO_ALIAS.items()}

# SECONDARY: ENSG IDs (for TCGA cross-reference)
A3_GENES_ENSG = [
    "ENSG00000128383", "ENSG00000179750", "ENSG00000244509",
    "ENSG00000243811", "ENSG00000128394", "ENSG00000239713",
    "ENSG00000100298",
]

A3_ENSG_TO_ALIAS = {
    "ENSG00000128383": "A3A", "ENSG00000179750": "A3B",
    "ENSG00000244509": "A3C", "ENSG00000243811": "A3D",
    "ENSG00000128394": "A3F", "ENSG00000239713": "A3G",
    "ENSG00000100298": "A3H",
}

A3_ENSG_TO_SYMBOL = {
    "ENSG00000128383": "APOBEC3A", "ENSG00000179750": "APOBEC3B",
    "ENSG00000244509": "APOBEC3C", "ENSG00000243811": "APOBEC3D",
    "ENSG00000128394": "APOBEC3F", "ENSG00000239713": "APOBEC3G",
    "ENSG00000100298": "APOBEC3H",
}

# Backward compatibility aliases
A3_GENES = A3_GENES_SYMBOLS
A3_ID_TO_ALIAS = A3_SYMBOL_TO_ALIAS

# =============================================================================
# BIOMARKER GENES -- gene symbols
# =============================================================================
BIOMARKERS_SYMBOLS = {
    'CD4+ T cell':              ['IL7R', 'CD3E'],
    'B cell':                   ['CD79A', 'MS4A1'],
    'fibroblast':               ['DCN', 'COL1A1'],
    'CD8+ T cell':              ['CD8A', 'CCL5'],
    'basal cell':               ['KRT5', 'TACSTD2'],
    'macrophage':               ['CD68', 'SPI1'],
    'endothelial cell':         ['PECAM1', 'VWF'],
    'smooth muscle cell':       ['TAGLN', 'RGS5'],
    'regulatory T cell':        ['TIGIT', 'CTLA4'],
    'mast cell':                ['TPSAB1', 'CPA3'],
    'plasmacytoid dendritic':   ['IL3RA', 'LILRA4'],
    'myeloid dendritic':        ['LAMP3', 'CCR7'],
}

BIOMARKERS = BIOMARKERS_SYMBOLS

# =============================================================================
# GENERAL
# =============================================================================
VERBOSE = True
RANDOM_SEED = 42


# =============================================================================
# HELPERS
# =============================================================================

def banner(title, char="=", width=100):
    """Print a visible section banner to stdout."""
    line = char * width
    print(f"\n{line}\n{title}\n{line}", flush=True)


def log(msg):
    """Print a log message (respects VERBOSE flag)."""
    if VERBOSE:
        print(msg, flush=True)


def ensure_dir(path):
    """Create directory if it doesn't exist."""
    os.makedirs(path, exist_ok=True)
    return path


def load_ensg_to_symbol():
    """Load ENSG-to-symbol mapping from Figure 2 reference."""
    import json
    if os.path.exists(ENSG_TO_SYMBOL_PATH):
        with open(ENSG_TO_SYMBOL_PATH) as f:
            return json.load(f)
    return {}


def convert_tcga_genes_to_symbols(gene_list, ensg_to_symbol):
    """Convert ENSG IDs to gene symbols where possible."""
    converted = []
    for g in gene_list:
        if g.startswith("ENSG") and g in ensg_to_symbol:
            converted.append(ensg_to_symbol[g])
        else:
            converted.append(g)
    return converted


def load_harris_interactors(path=None, a3b_only=False):
    """
    Load Harris et al. A3 interactor gene list.

    Parameters
    ----------
    path : str or None
        Path to interactor file. If None, uses HARRIS_ALL_PATH or HARRIS_A3B_PATH.
    a3b_only : bool
        If True and path is None, loads A3B-specific interactors (52 genes).
        If False and path is None, loads all A3 interactors (175 genes).

    Returns
    -------
    list of str
        Gene symbols of Harris interactors.
    """
    import pandas as pd
    if path is None:
        path = HARRIS_A3B_PATH if a3b_only else HARRIS_ALL_PATH

    if not os.path.exists(path):
        log(f"  WARNING: Harris interactors file not found: {path}")
        return []

    df = pd.read_csv(path, sep="\t")

    # Identify gene column
    gene_col = None
    for col in ["gene_symbol", "gene", "Gene", "symbol", "Gene_Symbol"]:
        if col in df.columns:
            gene_col = col
            break

    if gene_col is None:
        # Try first column as fallback
        gene_col = df.columns[0]
        log(f"  WARNING: No standard gene column found, using '{gene_col}'")

    genes = df[gene_col].dropna().astype(str).tolist()
    log(f"  Loaded {len(genes)} Harris interactors from {os.path.basename(path)}")
    return genes


def load_tcga_bulk_communities(cancer_type="TCGA-HNSC"):
    """
    Load TCGA bulk community assignments from Figure 2 for cross-reference.

    Parameters
    ----------
    cancer_type : str
        Cancer type label (default "TCGA-HNSC").

    Returns
    -------
    dict
        {gene_ensg: community_id} mapping, or empty dict if file not found.
    """
    import pandas as pd
    comm_dir = os.path.join(FIG2_COMMUNITIES, cancer_type)
    comm_path = os.path.join(comm_dir, f"{cancer_type}_best_partition.csv")

    if not os.path.exists(comm_path):
        log(f"  WARNING: TCGA bulk partition not found: {comm_path}")
        return {}

    df = pd.read_csv(comm_path)

    # Identify gene and community columns
    gene_col = None
    for col in ["gene", "gene_symbol", "Gene"]:
        if col in df.columns:
            gene_col = col
            break
    if gene_col is None:
        gene_col = df.columns[0]

    comm_col = None
    for col in ["community", "Community", "comm"]:
        if col in df.columns:
            comm_col = col
            break
    if comm_col is None:
        comm_col = df.columns[-1]

    result = dict(zip(df[gene_col].astype(str), df[comm_col].astype(int)))
    log(f"  Loaded {len(result)} TCGA bulk community assignments from {cancer_type}")
    return result
