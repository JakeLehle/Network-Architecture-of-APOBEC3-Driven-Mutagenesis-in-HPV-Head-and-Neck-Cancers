#!/usr/bin/env python3
"""
Step02_Supplemental_Marker_Validation.py
=========================================

Supplemental Figure: Marker gene expression UMAPs validating popV cell type
annotations from ClusterCatcher.

Pipeline:
  1. Load adata_final.h5ad (ClusterCatcher output with popV annotations)
  2. Run Wilcoxon rank-sum differential expression to identify top 20 marker
     genes per cell type (transparency: shows how curated markers were selected)
  3. Save full marker gene table to TSV
  4. Generate 6x4 UMAP grid showing 2 curated marker genes per cell type
  5. Save supplemental figure as PDF and PNG (300 DPI)

Input:
  - ClusterCatcher adata_final.h5ad with 'final_annotation' in .obs

Output (all to data/FIG_3/):
  - Supplemental_Marker_Genes_Top20_Per_CellType.tsv
  - selected_marker_genes.tsv
  - Supplemental_Figure_Marker_Validation.pdf
  - Supplemental_Figure_Marker_Validation.png

Usage:
  conda run -n NETWORK python Step02_Supplemental_Marker_Validation.py

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from datetime import datetime


# =============================================================================
# CONFIGURATION
# =============================================================================

# Input paths (ClusterCatcher output)
CLUSTERCATCHER_OUTPUT = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/ClusterCatcher_results"
ADATA_PATH = os.path.join(CLUSTERCATCHER_OUTPUT, "signatures", "adata_final.h5ad")

# Output paths (paper figure directory)
OUTPUT_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_3"

# Plot parameters
DPI = 300
POINT_SIZE = 10

# Custom colormap: slate gray to scarlet red
SLATE_GRAY = '#708090'
BRIGHT_RED = '#FF2400'
CMAP_MARKERS = LinearSegmentedColormap.from_list(
    'slate_to_red',
    [SLATE_GRAY, '#A0522D', '#CD5C5C', '#FF4500', BRIGHT_RED]
)

# Curated marker genes: 2 well-established markers per cell type
# Selected from Wilcoxon rank-sum results and validated against literature
CELL_TYPE_MARKERS = {
    'CD4-Positive alpha-beta T cell': ['IL7R', 'CD3E'],
    'B Cells': ['CD79A', 'MS4A1'],
    'fibroblasts': ['DCN', 'COL1A1'],
    'CD8-Positive, alpha-beta T cell': ['CD8A', 'CCL5'],
    'basal cell': ['KRT5', 'TACSTD2'],
    'macrophage': ['CD68', 'SPI1'],
    'endothelial cell': ['PECAM1', 'VWF'],
    'smooth muscle cell': ['TAGLN', 'RGS5'],
    'regulatory T cell': ['TIGIT', 'CTLA4'],
    'mast cell': ['TPSAB1', 'CPA3'],
    'plasmacytoid dendritic cell': ['IL3RA', 'LILRA4'],
    'myeloid dendritic cell': ['LAMP3', 'CCR7'],
}

# Column layout order: 2 cell types per column, 6 columns total
CELL_TYPE_ORDER = [
    'CD4-Positive alpha-beta T cell', 'B Cells',                    # Column 0
    'CD8-Positive, alpha-beta T cell', 'regulatory T cell',         # Column 1
    'macrophage', 'myeloid dendritic cell',                         # Column 2
    'plasmacytoid dendritic cell', 'mast cell',                     # Column 3
    'fibroblasts', 'smooth muscle cell',                            # Column 4
    'basal cell', 'endothelial cell',                               # Column 5
]


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def log(msg):
    """Print timestamped log message."""
    timestamp = datetime.now().strftime('%H:%M:%S')
    print(f"[{timestamp}] {msg}", flush=True)


# =============================================================================
# STEP 1: LOAD DATA
# =============================================================================

def load_data():
    """Load adata_final.h5ad and verify popV annotations are present."""
    log("STEP 1: LOADING DATA")
    log(f"  Input: {ADATA_PATH}")

    if not os.path.exists(ADATA_PATH):
        log(f"  ERROR: File not found: {ADATA_PATH}")
        sys.exit(1)

    adata = sc.read_h5ad(ADATA_PATH)
    log(f"  Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # Verify popV annotations exist
    if 'final_annotation' not in adata.obs.columns:
        log("  ERROR: 'final_annotation' column not found in adata.obs")
        log(f"  Available columns: {list(adata.obs.columns)}")
        sys.exit(1)

    # Log cell type counts
    ct_counts = adata.obs['final_annotation'].value_counts()
    log(f"  Cell types found: {len(ct_counts)}")
    for ct, n in ct_counts.items():
        log(f"    {ct}: {n:,}")

    return adata


# =============================================================================
# STEP 2: DIFFERENTIAL EXPRESSION (MARKER GENE IDENTIFICATION)
# =============================================================================

def compute_marker_genes(adata):
    """
    Run Wilcoxon rank-sum test to identify top 20 marker genes per cell type.
    This step is included for full transparency of how curated markers were
    selected and validated.
    """
    log("\n" + "=" * 80)
    log("STEP 2: DIFFERENTIAL EXPRESSION (MARKER GENE IDENTIFICATION)")
    log("=" * 80)
    log("  Method: Wilcoxon rank-sum test")
    log("  Groupby: final_annotation")
    log("  Top N genes: 20 per cell type")

    sc.tl.rank_genes_groups(
        adata,
        groupby='final_annotation',
        method='wilcoxon',
        n_genes=20,
    )

    # Extract top 20 markers per cell type into a DataFrame
    cell_types = adata.obs['final_annotation'].unique()
    marker_dict = {}

    for ct in cell_types:
        genes = sc.get.rank_genes_groups_df(adata, group=ct)['names'].head(20).tolist()
        marker_dict[ct] = genes

    marker_df = pd.DataFrame(marker_dict)
    marker_df.index = range(1, 21)
    marker_df.index.name = 'rank'

    # Save full marker table
    full_path = os.path.join(OUTPUT_DIR, "Supplemental_Marker_Genes_Top20_Per_CellType.tsv")
    marker_df.to_csv(full_path, sep='\t')
    log(f"  Saved full marker table: {full_path}")

    # Log where curated markers rank
    log("\n  Curated marker rankings:")
    for ct, markers in CELL_TYPE_MARKERS.items():
        if ct in marker_dict:
            for gene in markers:
                if gene in marker_dict[ct]:
                    rank = marker_dict[ct].index(gene) + 1
                    log(f"    {ct}: {gene} = rank {rank}")
                else:
                    log(f"    {ct}: {gene} = not in top 20")

    return marker_df


# =============================================================================
# STEP 3: SAVE CURATED MARKER TABLE
# =============================================================================

def save_curated_markers():
    """Save the curated 2-gene-per-cell-type marker table."""
    log("\n" + "=" * 80)
    log("STEP 3: SAVING CURATED MARKER GENE TABLE")
    log("=" * 80)

    curated_df = pd.DataFrame(CELL_TYPE_MARKERS)
    curated_df.index = ['Marker_1', 'Marker_2']
    curated_df.index.name = 'marker_rank'

    curated_path = os.path.join(OUTPUT_DIR, "selected_marker_genes.tsv")
    curated_df.to_csv(curated_path, sep='\t')
    log(f"  Saved curated marker table: {curated_path}")
    log(f"  Cell types: {len(CELL_TYPE_MARKERS)}")
    log(f"  Markers per cell type: 2")

    return curated_df


# =============================================================================
# STEP 4: GENERATE SUPPLEMENTAL UMAP FIGURE
# =============================================================================

def plot_marker_validation_umaps(adata):
    """
    Generate 6-column x 4-row UMAP grid.
    Each column shows one pair of cell types (2 cell types x 2 markers = 4 rows).
    """
    log("\n" + "=" * 80)
    log("STEP 4: GENERATING SUPPLEMENTAL MARKER VALIDATION UMAPS")
    log("=" * 80)

    fig, axes = plt.subplots(4, 6, figsize=(24, 16))
    plt.subplots_adjust(wspace=0.1, hspace=0.35)

    for col_idx in range(6):
        # Two cell types per column
        ct1 = CELL_TYPE_ORDER[col_idx * 2]
        ct2 = CELL_TYPE_ORDER[col_idx * 2 + 1]

        # Cell type 1: rows 0 and 1
        for gene_idx, gene in enumerate(CELL_TYPE_MARKERS[ct1]):
            ax = axes[gene_idx, col_idx]
            sc.pl.umap(
                adata,
                color=gene,
                cmap=CMAP_MARKERS,
                ax=ax,
                show=False,
                title=gene,
                colorbar_loc='right' if col_idx == 5 else None,
                frameon=False,
                size=POINT_SIZE,
                na_color=SLATE_GRAY,
            )
            ax.set_title(gene, fontsize=30, fontweight='bold')
            ax.set_xlabel('')
            ax.set_ylabel('')

            # Log expression stats
            if gene in adata.var_names:
                expr = adata[:, gene].X
                if hasattr(expr, 'toarray'):
                    expr = expr.toarray().flatten()
                else:
                    expr = np.asarray(expr).flatten()
                nonzero = (expr > 0).sum()
                log(f"  {ct1} | {gene}: {nonzero:,}/{adata.n_obs:,} cells expressing "
                    f"({100 * nonzero / adata.n_obs:.1f}%)")

        # Cell type 2: rows 2 and 3
        for gene_idx, gene in enumerate(CELL_TYPE_MARKERS[ct2]):
            ax = axes[gene_idx + 2, col_idx]
            sc.pl.umap(
                adata,
                color=gene,
                cmap=CMAP_MARKERS,
                ax=ax,
                show=False,
                title=gene,
                colorbar_loc='right' if col_idx == 5 else None,
                frameon=False,
                size=POINT_SIZE,
                na_color=SLATE_GRAY,
            )
            ax.set_title(gene, fontsize=30, fontweight='bold')
            ax.set_xlabel('')
            ax.set_ylabel('')

            if gene in adata.var_names:
                expr = adata[:, gene].X
                if hasattr(expr, 'toarray'):
                    expr = expr.toarray().flatten()
                else:
                    expr = np.asarray(expr).flatten()
                nonzero = (expr > 0).sum()
                log(f"  {ct2} | {gene}: {nonzero:,}/{adata.n_obs:,} cells expressing "
                    f"({100 * nonzero / adata.n_obs:.1f}%)")

    # Overall title
    fig.suptitle(
        'Marker Gene Expression Validating popV Cell Type Annotations',
        fontsize=28, fontweight='bold', y=0.98,
    )

    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # Save
    pdf_path = os.path.join(OUTPUT_DIR, "Supplemental_Figure_Marker_Validation.pdf")
    png_path = os.path.join(OUTPUT_DIR, "Supplemental_Figure_Marker_Validation.png")
    plt.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    plt.savefig(png_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close()

    log(f"\n  Saved: {pdf_path}")
    log(f"  Saved: {png_path}")


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def main():
    """Run the complete supplemental marker validation pipeline."""

    log("=" * 80)
    log("SUPPLEMENTAL MARKER VALIDATION PIPELINE")
    log("=" * 80)

    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    log(f"Output directory: {OUTPUT_DIR}")
    log("")

    # Step 1: Load data
    adata = load_data()

    # Step 2: Compute marker genes (Wilcoxon rank-sum)
    compute_marker_genes(adata)

    # Step 3: Save curated marker table
    save_curated_markers()

    # Step 4: Generate UMAP figure
    plot_marker_validation_umaps(adata)

    # Summary
    log("\n" + "=" * 80)
    log("SUPPLEMENTAL MARKER VALIDATION COMPLETE")
    log("=" * 80)
    log(f"\nOutput files in {OUTPUT_DIR}:")
    for f in sorted(os.listdir(OUTPUT_DIR)):
        fpath = os.path.join(OUTPUT_DIR, f)
        if os.path.isfile(fpath):
            size = os.path.getsize(fpath)
            log(f"  {f} ({size / 1024:.1f} KB)")

    log("\nSupplemental figure mapping:")
    log("  Supp. Table  : Supplemental_Marker_Genes_Top20_Per_CellType.tsv")
    log("  Supp. Table  : selected_marker_genes.tsv")
    log("  Supp. Figure : Supplemental_Figure_Marker_Validation.pdf/.png")
    log("")


if __name__ == "__main__":
    main()
