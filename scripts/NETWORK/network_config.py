"""
network_config.py

Shared configuration for the APOBEC3-SBS2 Network Analysis Pipeline (Figure 2).

All file paths, analysis parameters, and gene definitions live here.
Every Step script imports this file — edit paths/thresholds in one place.

Usage:
    from network_config import *
"""

import os
import sys

# =============================================================================
# PROJECT ROOT
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"

# =============================================================================
# INPUT PATHS (from Figure 1 pipeline)
# =============================================================================
# Pan-cancer expression matrix (rows 0-2 = annotations, rows 3+ = samples)
#   Row 0: ENSG IDs (versioned)
#   Row 1: Gene symbols
#   Row 2: Gene biotypes
#   Columns 0-4: Project_ID, Tissue_Type, Case_ID, File_ID, Entity_ID
TCGA_EXPRESSION_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_1", "TCGA_master_FPKM_UQ.tsv")

# COSMIC SBS mutation signature weights per TCGA sample
#   Key column: TCGA_Gene_Expression_Entity_ID (full 28-char TCGA barcode)
#   Contains SBS1, SBS2, SBS3, ... columns with signature weights
MUTATION_SIGNATURE_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_1", "Mutation_Table_Tumors_TCGA.tsv")

# Column in the mutation signature file that holds the TCGA barcode
SIGNATURE_SAMPLE_COL = "TCGA_Gene_Expression_Entity_ID"

# =============================================================================
# OUTPUT ROOT (Figure 2)
# =============================================================================
FIG2_ROOT = os.path.join(PROJECT_ROOT, "data", "FIG_2")

# Per-step output directories (created automatically by each script)
DIR_01_CLEANED    = os.path.join(FIG2_ROOT, "01_cleaned_expression")
DIR_02_MERGED     = os.path.join(FIG2_ROOT, "02_merged_with_SBS")
DIR_03_DIFFEXPR   = os.path.join(FIG2_ROOT, "03_differential_expression")
DIR_04_GROUPS     = os.path.join(FIG2_ROOT, "04_TOP_BOTTOM_groups")
DIR_05_NETWORKS   = os.path.join(FIG2_ROOT, "05_correlation_networks")
DIR_06_COMMUNITIES = os.path.join(FIG2_ROOT, "06_communities")
DIR_07_CENTRALITY = os.path.join(FIG2_ROOT, "07_centrality_metrics")

# =============================================================================
# CANCER TYPE(S) TO ANALYZE
# =============================================================================
CANCER_TYPES = ["TCGA-HNSC"]   # extend list for pan-cancer runs

# =============================================================================
# CLINICAL METADATA COLUMNS (first 5 columns of the master TSV)
# =============================================================================
CLINICAL_COLS = ["Project_ID", "Tissue_Type", "Case_ID", "File_ID", "Entity_ID"]
N_CLINICAL = len(CLINICAL_COLS)   # = 5

# =============================================================================
# APOBEC3 GENE DEFINITIONS (ENSG IDs, cleaned — no version suffix)
# =============================================================================
A3_GENES = [
    "ENSG00000128383",  # A3A
    "ENSG00000179750",  # A3B
    "ENSG00000244509",  # A3C
    "ENSG00000243811",  # A3D
    "ENSG00000128394",  # A3F
    "ENSG00000239713",  # A3G
    "ENSG00000100298",  # A3H
]

A3_ID_TO_ALIAS = {
    "ENSG00000128383": "A3A",
    "ENSG00000179750": "A3B",
    "ENSG00000244509": "A3C",
    "ENSG00000243811": "A3D",
    "ENSG00000128394": "A3F",
    "ENSG00000239713": "A3G",
    "ENSG00000100298": "A3H",
}

A3_ALIAS_TO_ID = {v: k for k, v in A3_ID_TO_ALIAS.items()}

# =============================================================================
# BIOMARKER GENES (cell-type markers — used for annotation, NOT force-kept)
# =============================================================================
BIOMARKERS = [
    # CD4+ T cell
    "ENSG00000168685",  # IL7R
    "ENSG00000198851",  # CD3E
    # B cell
    "ENSG00000007314",  # CD79A
    "ENSG00000156738",  # MS4A1 (CD20)
    # Fibroblast
    "ENSG00000011465",  # DCN
    "ENSG00000108821",  # COL1A1
    # CD8+ T cell
    "ENSG00000153563",  # CD8A
    "ENSG00000271503",  # CCL5
    # Basal cell
    "ENSG00000186081",  # KRT5
    "ENSG00000166886",  # TACSTD2
    # Macrophage
    "ENSG00000129226",  # CD68
    "ENSG00000066336",  # SPI1
    # Endothelial
    "ENSG00000261371",  # PECAM1
    "ENSG00000110799",  # VWF
    # Smooth muscle
    "ENSG00000149591",  # TAGLN
    "ENSG00000143248",  # RGS5
    # Treg
    "ENSG00000181847",  # TIGIT
    "ENSG00000163599",  # CTLA4
    # Mast cell
    "ENSG00000172236",  # TPSAB1
    "ENSG00000163736",  # CPA3
    # pDC
    "ENSG00000185291",  # IL3RA
    "ENSG00000162692",  # LILRA4
    # mDC
    "ENSG00000110395",  # LAMP3
    "ENSG00000126353",  # CCR7
]

# =============================================================================
# STEP 02 — MERGE PARAMETERS
# =============================================================================
# How to shorten TCGA barcodes for merge (first N dash-separated parts)
BARCODE_PARTS_FOR_MERGE = 4   # TCGA-XX-YYYY-ZZZZ

# =============================================================================
# STEP 03 — DIFFERENTIAL EXPRESSION PARAMETERS
# =============================================================================
# SBS2 percentile thresholds for high vs low groups (for t-test)
SBS2_HIGH_PERCENTILE = 0.80
SBS2_LOW_PERCENTILE  = 0.20

# p-value threshold (raw, before FDR)
P_VALUE_THRESHOLD = 0.025

# Always retain A3 genes even if not significant
FORCE_KEEP_A3 = True

# =============================================================================
# STEP 04 — GROUP DEFINITION PARAMETERS
# =============================================================================
# A3 score uses A3A + A3B (normalized within cohort)
# High-A3 cutoff: median of A3 score
A3_SCORE_PERCENTILE = 0.50   # median

# Within high-A3 tumors, split SBS2:
GROUP_SBS2_HIGH_PERCENTILE = 0.80   # TOP group = top 20%
GROUP_SBS2_LOW_PERCENTILE  = 0.20   # BOTTOM group = bottom 20%

# Minimum group size to proceed
MIN_GROUP_SIZE = 8

# =============================================================================
# STEP 05 — CORRELATION NETWORK PARAMETERS
# =============================================================================
CORRELATION_METHOD = "spearman"
CORR_THRESHOLD     = 0.80    # |rho| threshold for TOP/BOTTOM network edges
DIFF_THRESHOLD     = 0.40    # |diff| threshold for DIFF network edges

# =============================================================================
# STEP 06 — COMMUNITY DETECTION PARAMETERS
# =============================================================================
COMMUNITY_METHOD    = "leiden"
COMMUNITY_RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4]
RUNS_PER_RESOLUTION = 15
COMMUNITY_BASE_SEED = 42
USE_LARGEST_COMPONENT = True
TARGET_BIG_COMMUNITIES = 8    # merge to top K + Other
MIN_COMMUNITY_SIZE     = 10   # smaller communities get merged

# =============================================================================
# GENERAL
# =============================================================================
VERBOSE = True
RANDOM_SEED = 42


# =============================================================================
# HELPER: banner for console logging
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
