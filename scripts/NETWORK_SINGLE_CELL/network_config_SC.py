#!/usr/bin/env python3
"""
network_config_SC.py
====================

Centralized configuration for the Figure 4 Single-Cell Network Analysis.
Mirrors network_config.py from Figure 2 (TCGA bulk) but adapted for
single-cell basal epithelial cell expression data.

Key differences from Figure 2 (TCGA bulk):
  - Input is single-cell expression (cells = samples)
  - No TCGA clinical metadata — groups defined by SBS2 weight
  - Gene IDs are gene symbols (not ENSG IDs)
  - Correlation thresholds may need relaxation for sparser data
  - Output goes to data/FIG_4/

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os

# =============================================================================
# DIRECTORY PATHS
# =============================================================================
BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"

# Input: single-cell data from ClusterCatcher
SC_DATA_DIR = os.path.join(BASE_DIR, "data", "FIG_4", "00_input")

# Output directories (mirror Figure 2 structure)
FIG4_ROOT = os.path.join(BASE_DIR, "data", "FIG_4")

DIR_00_INPUT       = os.path.join(FIG4_ROOT, "00_input")
DIR_01_GROUPS      = os.path.join(FIG4_ROOT, "01_group_selection")
DIR_02_DE          = os.path.join(FIG4_ROOT, "02_differential_expression")
DIR_03_NETWORKS    = os.path.join(FIG4_ROOT, "03_correlation_networks")
DIR_04_COMMUNITIES = os.path.join(FIG4_ROOT, "04_communities")
DIR_05_CENTRALITY  = os.path.join(FIG4_ROOT, "05_centrality_metrics")
DIR_06_OVERLAP     = os.path.join(FIG4_ROOT, "06_overlap_analysis")
FIGURE_4_PANELS    = os.path.join(FIG4_ROOT, "FIGURE_4_PANELS")

# Figure 2 community data (for cross-reference)
FIG2_COMMUNITIES = os.path.join(BASE_DIR, "data", "FIG_2", "05_communities")

# Harris A3 interactor list
HARRIS_INTERACTORS = os.path.join(DIR_00_INPUT, "harris_A3B_interactors.txt")

# =============================================================================
# INPUT FILES (user must place these in data/FIG_4/00_input/)
# =============================================================================
ADATA_FINAL_PATH = os.path.join(DIR_00_INPUT, "adata_final.h5ad")
WEIGHTS_PATH     = os.path.join(DIR_00_INPUT, "signature_weights_per_cell.txt")

# =============================================================================
# STEP 00 — CELL SELECTION PARAMETERS
# =============================================================================
# Cell type to subset
TARGET_CELL_TYPE = "basal cell"

# SBS2 grouping: top percentile = HIGH, matched controls = LOW
SBS2_HIGH_PERCENTILE = 0.80    # Top 20% SBS2 basal cells → HIGH group
N_CONTROLS_MATCH     = True    # Match control count to HIGH count
SBS2_LOW_MAX_WEIGHT  = 0.01   # Maximum SBS2 weight for control eligibility

# =============================================================================
# STEP 01 — DIFFERENTIAL EXPRESSION PARAMETERS
# =============================================================================
# Wilcoxon rank-sum test between HIGH and LOW groups
MIN_CELLS_DETECTED  = 10       # Gene must be detected in >= N cells
RAW_P_THRESHOLD     = 0.05     # Raw p-value cutoff for gene selection
LOGFC_THRESHOLD     = 0        # No fold-change filter
FORCE_KEEP_A3       = True     # Always retain A3 genes in the network

# =============================================================================
# STEP 02 — CORRELATION NETWORK PARAMETERS
# =============================================================================
CORRELATION_METHOD = "spearman"

# NOTE: These may need adjustment for single-cell data sparsity.
# Start with same thresholds as TCGA; relax if network is too sparse.
CORR_THRESHOLD     = 0.80     # |rho| threshold for HIGH/LOW network edges
DIFF_THRESHOLD     = 0.70     # |delta-rho| threshold for DIFF network edges

# =============================================================================
# STEP 03 — COMMUNITY DETECTION PARAMETERS
# =============================================================================
COMMUNITY_METHOD    = "leiden"
COMMUNITY_RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0]
RUNS_PER_RESOLUTION = 15
COMMUNITY_BASE_SEED = 42
USE_LARGEST_COMPONENT = True
TARGET_BIG_COMMUNITIES = 14
MIN_COMMUNITY_SIZE     = 10

# =============================================================================
# APOBEC3 GENE DEFINITIONS (gene symbols for single-cell data)
# =============================================================================
A3_GENES_SYMBOLS = [
    "APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
    "APOBEC3F", "APOBEC3G", "APOBEC3H",
]

A3_SYMBOL_TO_ALIAS = {
    "APOBEC3A": "A3A",
    "APOBEC3B": "A3B",
    "APOBEC3C": "A3C",
    "APOBEC3D": "A3D",
    "APOBEC3F": "A3F",
    "APOBEC3G": "A3G",
    "APOBEC3H": "A3H",
}

# =============================================================================
# BIOMARKER GENES (gene symbols — for annotation cross-reference)
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
    """Print a log message."""
    if VERBOSE:
        print(msg, flush=True)


def ensure_dir(path):
    """Create directory if it doesn't exist."""
    os.makedirs(path, exist_ok=True)
    return path
