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
  - DIFF threshold auto-selected per network using component peak + LCC criterion
  - Leiden resolution auto-selected using modularity jump detection
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
FORCE_KEEP_A3       = True
FILTER_PROTEIN_CODING = False   # SC data already ~20K protein-coding genes

# =============================================================================
# STEP 02 -- CORRELATION NETWORK PARAMETERS
# =============================================================================
CORRELATION_METHOD = "spearman"

CORR_THRESHOLD     = 0.80     # |rho| threshold for HIGH/LOW network edges

# DIFF threshold: auto-selected per network comparison.
# The auto-selector (Auto_Select_Network_Parameters.py) finds the threshold
# where connected components peak and LCC < DIFF_THRESHOLD_MAX_LCC.
# This produces pathway-scale communities that are biologically interpretable.
#
# Typical auto-selected thresholds:
#   SBS2_VS_CNV:    0.65 (narrowest biological contrast)
#   SBS2_VS_NORMAL: 0.65 (cancer vs normal)
#   CNV_VS_NORMAL:  0.70 (widest biological contrast)
#
# Fallback value used if auto-selection is disabled or sweep data unavailable.
DIFF_THRESHOLD          = 0.40     # Fallback (only used if auto disabled)
DIFF_THRESHOLD_AUTO     = True     # Enable auto-selection from sweep data
DIFF_THRESHOLD_MAX_LCC  = 300      # Max LCC size for auto-selection

# Finer granularity in the 0.50-0.70 range where thresholds typically land
SWEEP_THRESHOLDS = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70,
                    0.75, 0.80, 0.85, 0.90]

# =============================================================================
# STEP 03 -- COMMUNITY DETECTION PARAMETERS
# =============================================================================
COMMUNITY_METHOD    = "leiden"
COMMUNITY_RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0]
RUNS_PER_RESOLUTION = 15
COMMUNITY_BASE_SEED = 42
USE_LARGEST_COMPONENT = True
TARGET_BIG_COMMUNITIES = 14
MIN_COMMUNITY_SIZE     = 10

# Resolution auto-selection: modularity jump detection with stability floor.
# Finds the resolution where modularity increases most (structural jump)
# and verifies ARI stability is above the floor.
RESOLUTION_AUTO     = True     # Enable auto-selection from resolution sweep
RESOLUTION_MIN_ARI  = 0.65    # Stability floor (TCGA bulk used ARI ~0.63)

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

BIOMARKERS = [gene for genes in BIOMARKERS_SYMBOLS.values() for gene in genes]

# =============================================================================
# GENERAL
# =============================================================================
VERBOSE = True
RANDOM_SEED = 42

# =============================================================================
# HELPERS
# =============================================================================
def banner(title, char="=", width=100):
    line = char * width
    print(f"\n{line}\n{title}\n{line}", flush=True)

def log(msg):
    if VERBOSE:
        print(msg, flush=True)

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path

def load_ensg_to_symbol():
    import json
    if os.path.exists(ENSG_TO_SYMBOL_PATH):
        with open(ENSG_TO_SYMBOL_PATH) as f:
            return json.load(f)
    else:
        log(f"WARNING: ENSG->symbol mapping not found at {ENSG_TO_SYMBOL_PATH}")
        return {}

def load_ensg_to_biotype():
    import json
    if os.path.exists(ENSG_TO_BIOTYPE_PATH):
        with open(ENSG_TO_BIOTYPE_PATH) as f:
            return json.load(f)
    else:
        log(f"WARNING: ENSG->biotype mapping not found at {ENSG_TO_BIOTYPE_PATH}")
        return {}

def load_harris_interactors():
    """
    Load Harris A3 interactors from TSV files.
    Files have columns: gene_symbol, A3_baits, source, R_loop_associated,
    confirmed_coIP, A3B_interactor.
    Returns (harris_all, harris_a3b) as sets of gene symbols.
    """
    import pandas as pd

    harris_all = set()
    harris_a3b = set()

    if os.path.exists(HARRIS_ALL_PATH):
        try:
            df = pd.read_csv(HARRIS_ALL_PATH, sep='\t')
            harris_all = set(df['gene_symbol'].dropna().values)
            log(f"  Harris interactors (all): {len(harris_all)} genes")
        except Exception as e:
            log(f"  WARNING: Could not load Harris interactors: {e}")
    else:
        log(f"  WARNING: Harris interactors not found: {HARRIS_ALL_PATH}")

    if os.path.exists(HARRIS_A3B_PATH):
        try:
            df = pd.read_csv(HARRIS_A3B_PATH, sep='\t')
            harris_a3b = set(df['gene_symbol'].dropna().values)
            log(f"  Harris interactors (A3B-specific): {len(harris_a3b)} genes")
        except Exception as e:
            log(f"  WARNING: Could not load A3B interactors: {e}")
    else:
        log(f"  WARNING: A3B interactors not found: {HARRIS_A3B_PATH}")

    return harris_all, harris_a3b

def load_tcga_bulk_communities():
    """
    Load TCGA bulk community assignments for overlap analysis.
    Converts ENSG IDs to gene symbols using ensg_to_symbol.json.
    Returns dict of {community_id: set of gene symbols}.
    """
    import pandas as pd

    tcga_comms = {}
    cancer_type = "TCGA-HNSC"
    part_path = os.path.join(FIG2_COMMUNITIES, cancer_type,
                              f"{cancer_type}_best_partition.csv")

    if not os.path.exists(part_path):
        log(f"  WARNING: TCGA partition not found: {part_path}")
        return tcga_comms

    ensg_to_sym = load_ensg_to_symbol()
    df = pd.read_csv(part_path)

    for _, row in df.iterrows():
        c = int(row['community'])
        gene = row['gene']
        # Convert ENSG to symbol
        if gene.startswith('ENSG'):
            symbol = ensg_to_sym.get(gene)
            if symbol is None:
                symbol = A3_ENSG_TO_SYMBOL.get(gene, gene)
        else:
            symbol = gene
        tcga_comms.setdefault(c, set()).add(symbol)

    log(f"  TCGA bulk communities: {len(tcga_comms)} communities, "
        f"{sum(len(v) for v in tcga_comms.values())} genes")

    return tcga_comms

def convert_tcga_genes_to_symbols(gene_set, ensg_to_symbol=None):
    """Convert TCGA ENSG IDs to gene symbols for overlap analysis."""
    if ensg_to_symbol is None:
        ensg_to_symbol = load_ensg_to_symbol()
    converted = set()
    for gene in gene_set:
        if gene.startswith("ENSG"):
            symbol = ensg_to_symbol.get(gene)
            if symbol:
                converted.add(symbol)
            else:
                a3_sym = A3_ENSG_TO_SYMBOL.get(gene)
                converted.add(a3_sym if a3_sym else gene)
        else:
            converted.add(gene)
    return converted
