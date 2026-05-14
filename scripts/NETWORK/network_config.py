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

# NEW: SigProfiler v3.4 counts (COSMIC v3.4, 86 signatures, WES-barcode-indexed)
NEW_COUNTS_PATH = os.path.join("/master/jlehle/SHARED/TCGA/VCF/SigProfiler_output",
                                "TCGA_SBS_signature_counts.tsv")

# Crosswalk: old file used ONLY for RNA -> WES barcode mapping
CROSSWALK_PATH = MUTATION_SIGNATURE_PATH  # same file, different role
CROSSWALK_RNA_COL = "TCGA_Gene_Expression_Entity_ID"
CROSSWALK_WES_COL = "Mutation_Signature__File_Orginal_Entity_ID"

# Column in the mutation signature file that holds the TCGA barcode
SIGNATURE_SAMPLE_COL = "TCGA_Gene_Expression_Entity_ID"

# =============================================================================
# OUTPUT ROOT (Figure 2)
# =============================================================================
FIG2_ROOT = os.path.join(PROJECT_ROOT, "data", "FIG_2")

# Per-step output directories (created automatically by each script)
DIR_01_CLEANED     = os.path.join(FIG2_ROOT, "01_cleaned_expression")
DIR_02_MERGED      = os.path.join(FIG2_ROOT, "02_merged_with_SBS")
DIR_03_DIFFEXPR    = os.path.join(FIG2_ROOT, "03_differential_expression")
# NOTE: Step03 now also outputs group definitions (HIGH/LOW) and filtered gene lists.
#       The old DIR_04_GROUPS has been eliminated.
DIR_04_NETWORKS    = os.path.join(FIG2_ROOT, "04_correlation_networks")
DIR_05_COMMUNITIES = os.path.join(FIG2_ROOT, "05_communities")
DIR_06_CENTRALITY  = os.path.join(FIG2_ROOT, "06_centrality_metrics")

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
    "ENSG00000105369",  # CD79A
    "ENSG00000156738",  # MS4A1 (CD20)
    # Fibroblast
    "ENSG00000011465",  # DCN
    "ENSG00000108821",  # COL1A1
    # CD8+ T cell
    "ENSG00000153563",  # CD8A
    "ENSG00000271503",  # CCL5
    # Basal cell
    "ENSG00000186081",  # KRT5
    "ENSG00000184292",  # TACSTD2
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
    "ENSG00000163751",  # CPA3
    # pDC
    "ENSG00000185291",  # IL3RA
    "ENSG00000239961",  # LILRA4
    # mDC
    "ENSG00000078081",  # LAMP3
    "ENSG00000126353",  # CCR7
]

# =============================================================================
# STEP 02 — MERGE PARAMETERS
# =============================================================================
# (No barcode shortening needed — direct Entity_ID merge)

# =============================================================================
# STEP 03 — DIFFERENTIAL EXPRESSION PARAMETERS
# =============================================================================
# Gene filtering (aligned with single-cell workflow)
MIN_SAMPLES_DETECTED = 20     # gene must have FPKM-UQ > 0 in at least this many samples
                               # (analog of sc.pp.filter_genes(min_cells=20))

# Group definition (A3-controlled, within Step03)
# 1. Sum A3A + A3B raw FPKM-UQ per tumor
# 2. Keep tumors above median of summed A3A+A3B
# 3. Within high-A3: rank by SBS2, take top/bottom quartiles
A3_SUM_PERCENTILE = 0.50      # median — minimum summed A3A+A3B to include
SBS2_HIGH_PERCENTILE = 0.75   # top 25% of SBS2 within high-A3 = SBS2-HIGH group
SBS2_LOW_PERCENTILE  = 0.25   # bottom 25% of SBS2 within high-A3 = SBS2-LOW group

# Statistical test: Wilcoxon rank-sum on log1p(FPKM-UQ)
# (matches sc.tl.rank_genes_groups method='wilcoxon')

# Selection thresholds (using raw p-value — FDR too conservative for bulk TCGA)
# NOTE: When adapting the single-cell pipeline, switch to the same raw p approach
RAW_P_THRESHOLD = 0.05            # raw Wilcoxon p-value
LOGFC_THRESHOLD = 0               # No fold-change filter (bulk TCGA effect sizes too small)
                                   # NOTE: Match this in the single-cell pipeline

# Legacy FDR threshold (kept for reference, not used in selection)
FDR_THRESHOLD  = 0.05

# Always retain A3 genes even if not significant
FORCE_KEEP_A3 = True

# Legacy (kept for reference, no longer used in Step03)
P_VALUE_THRESHOLD = 0.025

# =============================================================================
# STEP 04 — GROUP DEFINITION PARAMETERS
# =============================================================================
# Step04 now INHERITS groups from Step03 (same TOP/BOTTOM groups).
# These legacy parameters are kept for reference but are not used.
A3_SCORE_PERCENTILE = 0.50
GROUP_SBS2_HIGH_PERCENTILE = 0.75
GROUP_SBS2_LOW_PERCENTILE  = 0.25

# Minimum group size to proceed
MIN_GROUP_SIZE = 8

# =============================================================================
# STEP 04 — CORRELATION NETWORK PARAMETERS
# =============================================================================
CORRELATION_METHOD = "spearman"
CORR_THRESHOLD     = 0.80    # |rho| threshold for TOP/BOTTOM network edges
DIFF_THRESHOLD     = 0.70    # |diff| threshold for DIFF network edges
                              # Selected via sweep: 2195 nodes, 8364 edges,
                              # avg degree 7.6, modularity 0.44 at res=1.0
SWEEP_THRESHOLDS = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90]

# =============================================================================
# STEP 05 — COMMUNITY DETECTION PARAMETERS
# =============================================================================
COMMUNITY_METHOD    = "leiden"
COMMUNITY_RESOLUTIONS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]   # Single resolution (selected via sweep)
RUNS_PER_RESOLUTION = 15
COMMUNITY_BASE_SEED = 42
USE_LARGEST_COMPONENT = True
TARGET_BIG_COMMUNITIES = 14   # merge to top K + Other
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
