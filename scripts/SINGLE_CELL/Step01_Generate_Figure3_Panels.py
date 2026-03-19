#!/usr/bin/env python3
"""
Step01_Generate_Figure3_Panels.py
=================================

Generates the publication-ready panels for Figure 3 that are NOT produced
by the ClusterCatcher pipeline directly. ClusterCatcher produces:
  - Figure 3a: Stacked bar + cluster UMAP + final annotation UMAP
  - Supplemental Figure 1: Per-cell popV prediction UMAP

This script produces:
  - Figure 3b: 4 UMAPs (normalized mutations, SBS2, A3A, A3B)
  - Figure 3c: Aggregated 96-context signature comparison vs COSMIC SBS2
  - Supplemental Figure 2: All A3 family member expression UMAPs

Input:
  - adata_final.h5ad from ClusterCatcher (contains annotations, mutations,
    signature weights, cancer status, and normalized expression)
  - COSMIC_v3.4_SBS_GRCh38.txt (COSMIC reference signatures)
  - signature_weights_per_cell.txt (per-cell signature weights)
  - Mutation matrix from SComatic (for aggregated profile)

Output (→ data/FIG_3/):
  - Panel_3b_Normalized_Mutations_SBS2_A3A_A3B.pdf
  - Panel_3c_Aggregated_SBS2_vs_COSMIC.pdf
  - Supplemental_Figure_2_A3_Family_Expression.pdf

Usage:
  Run interactively or via SLURM after ClusterCatcher pipeline completes.

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
import matplotlib.patheffects as pe
from matplotlib.patches import Patch, Rectangle
from pathlib import Path

# =============================================================================
# CONFIGURATION
# =============================================================================

# Input paths — ClusterCatcher outputs
CLUSTERCATCHER_OUTPUT = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/ClusterCatcher_results"
ADATA_PATH = os.path.join(CLUSTERCATCHER_OUTPUT, "signatures", "adata_final.h5ad")
WEIGHTS_PATH = os.path.join(CLUSTERCATCHER_OUTPUT, "signatures", "signature_weights_per_cell.txt")
MUTATION_MATRIX_PATH = os.path.join(CLUSTERCATCHER_OUTPUT, "signatures", "mutation_matrix_96context.txt")
COSMIC_FILE = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/COSMIC_v3.4_SBS_GRCh38.txt"

# Output paths — Paper figure directory
OUTPUT_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_3"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Plot parameters
POINT_SIZE = 15
DPI = 300
CMAP_EXPRESSION = 'plasma'
CMAP_MUTATIONS = 'plasma'

# COSMIC signature colors (standard 6-class mutation type colors)
COSMIC_COLORS = {
    'C>A': '#1EBFF0',
    'C>G': '#050708',
    'C>T': '#E62725',
    'T>A': '#CBCACB',
    'T>C': '#A1CE63',
    'T>G': '#EDB6C2'
}

# A3 family members to plot
A3_FAMILY = ['APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D',
             'APOBEC3F', 'APOBEC3G', 'APOBEC3H']
A3_SHORT = ['A3A', 'A3B', 'A3C', 'A3D', 'A3F', 'A3G', 'A3H']


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def log(msg):
    """Print timestamped log message."""
    from datetime import datetime
    timestamp = datetime.now().strftime('%H:%M:%S')
    print(f"[{timestamp}] {msg}", flush=True)


def get_gene_expression(adata, gene_symbol):
    """
    Extract expression values for a gene by symbol.
    Handles both gene_symbol column and direct var_names lookup.
    """
    # Try gene_symbol column first
    if 'gene_symbol' in adata.var.columns:
        mask = adata.var['gene_symbol'] == gene_symbol
        if mask.any():
            idx = adata.var.index[mask][0]
            values = adata[:, idx].X
            if hasattr(values, 'toarray'):
                values = values.toarray().flatten()
            return values.flatten()

    # Try feature_name column
    if 'feature_name' in adata.var.columns:
        mask = adata.var['feature_name'] == gene_symbol
        if mask.any():
            idx = adata.var.index[mask][0]
            values = adata[:, idx].X
            if hasattr(values, 'toarray'):
                values = values.toarray().flatten()
            return values.flatten()

    # Try direct var_names
    if gene_symbol in adata.var_names:
        values = adata[:, gene_symbol].X
        if hasattr(values, 'toarray'):
            values = values.toarray().flatten()
        return values.flatten()

    log(f"  WARNING: Gene {gene_symbol} not found in adata")
    return None


# =============================================================================
# STEP 1: LOAD DATA
# =============================================================================

def load_data():
    """Load ClusterCatcher outputs."""
    log("=" * 80)
    log("STEP 1: LOADING CLUSTERCATCHER OUTPUTS")
    log("=" * 80)

    log(f"Loading AnnData: {ADATA_PATH}")
    adata = sc.read_h5ad(ADATA_PATH)
    log(f"  Shape: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    log(f"  Cell types: {adata.obs['final_annotation'].nunique()}")

    # Verify required columns exist
    required_obs = ['final_annotation', 'normalized_total_mutations', 'SBS2']
    for col in required_obs:
        if col in adata.obs.columns:
            log(f"  ✓ {col} found in adata.obs")
        else:
            log(f"  ✗ {col} NOT found — check ClusterCatcher output")

    # Load signature weights if not already in obs
    if 'SBS2' not in adata.obs.columns and os.path.exists(WEIGHTS_PATH):
        log(f"\nLoading signature weights: {WEIGHTS_PATH}")
        weights = pd.read_csv(WEIGHTS_PATH, sep='\t', index_col=0)
        for sig in weights.index:
            adata.obs[sig] = adata.obs.index.map(weights.loc[sig]).fillna(0)
        log(f"  Added {len(weights.index)} signatures to adata.obs")

    log(f"\nCell type distribution:")
    for ct, count in adata.obs['final_annotation'].value_counts().items():
        log(f"  {ct}: {count:,}")

    return adata


# =============================================================================
# STEP 2: FIGURE 3b — MUTATION AND EXPRESSION UMAPS
# =============================================================================

def plot_figure_3b(adata):
    """
    Generate Figure 3b: 4-panel UMAP showing normalized somatic mutations,
    SBS2 signature weight, A3A expression, and A3B expression.
    """
    log("\n" + "=" * 80)
    log("STEP 2: GENERATING FIGURE 3b — MUTATION AND EXPRESSION UMAPS")
    log("=" * 80)

    # Add A3A and A3B expression to obs for plotting
    for gene_full, gene_short in [('APOBEC3A', 'A3A'), ('APOBEC3B', 'A3B')]:
        if gene_short not in adata.obs.columns:
            expr = get_gene_expression(adata, gene_full)
            if expr is not None:
                adata.obs[gene_short] = expr
                log(f"  Added {gene_short} expression to adata.obs")
                log(f"    Non-zero cells: {(expr > 0).sum():,} / {len(expr):,}")
            else:
                log(f"  WARNING: Could not extract {gene_full} expression")
                return

    # Define the 4 panels
    panels = [
        ('normalized_total_mutations', 'Normalized Somatic Mutations', CMAP_MUTATIONS),
        ('SBS2', 'SBS2 Signature Weight', CMAP_MUTATIONS),
        ('A3A', 'APOBEC3A Expression', CMAP_EXPRESSION),
        ('A3B', 'APOBEC3B Expression', CMAP_EXPRESSION),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(28, 6))

    for ax, (col, title, cmap) in zip(axes, panels):
        sc.pl.umap(
            adata,
            color=col,
            cmap=cmap,
            ax=ax,
            show=False,
            title=title,
            frameon=False,
            size=POINT_SIZE,
            colorbar_loc='right',
        )
        ax.set_title(title, fontsize=16, fontweight='bold')
        ax.set_xlabel('')
        ax.set_ylabel('')

    plt.tight_layout()

    # Save
    pdf_path = os.path.join(OUTPUT_DIR, "Panel_3b_Normalized_Mutations_SBS2_A3A_A3B.pdf")
    png_path = os.path.join(OUTPUT_DIR, "Panel_3b_Normalized_Mutations_SBS2_A3A_A3B.png")
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.savefig(png_path, dpi=DPI, bbox_inches='tight')
    plt.close()

    log(f"\n  Saved: {pdf_path}")
    log(f"  Saved: {png_path}")

    # Print summary statistics for the figure legend
    log("\n  Panel 3b summary statistics:")
    basal_mask = adata.obs['final_annotation'] == 'basal cell'
    for col, label, _ in panels:
        basal_vals = adata.obs.loc[basal_mask, col]
        other_vals = adata.obs.loc[~basal_mask, col]
        log(f"    {label}:")
        log(f"      Basal cells — median: {basal_vals.median():.4f}, "
            f"mean: {basal_vals.mean():.4f}, max: {basal_vals.max():.4f}")
        log(f"      Other cells — median: {other_vals.median():.4f}, "
            f"mean: {other_vals.mean():.4f}, max: {other_vals.max():.4f}")


# =============================================================================
# STEP 3: FIGURE 3c — AGGREGATED SIGNATURE COMPARISON
# =============================================================================

def plot_figure_3c(adata):
    """
    Generate Figure 3c: Aggregated 96-trinucleotide-context mutation profile
    from highest-SBS2 basal cells compared to the COSMIC SBS2 reference.
    """
    log("\n" + "=" * 80)
    log("STEP 3: GENERATING FIGURE 3c — AGGREGATED SIGNATURE vs COSMIC SBS2")
    log("=" * 80)

    # Load COSMIC reference
    log(f"Loading COSMIC signatures: {COSMIC_FILE}")
    cosmic = pd.read_csv(COSMIC_FILE, sep='\t', index_col=0)
    cosmic_sbs2 = cosmic['SBS2']
    contexts = cosmic_sbs2.index.tolist()
    log(f"  Loaded {len(contexts)} trinucleotide contexts")

    # Load mutation matrix
    if not os.path.exists(MUTATION_MATRIX_PATH):
        # Try alternative path from the SigProfiler results
        alt_path = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/ALL_Basal_Cell_SNP_matrix_for_SigProfiler.txt"
        if os.path.exists(alt_path):
            mutation_matrix_path = alt_path
        else:
            log("  WARNING: Mutation matrix not found. Skipping Figure 3c.")
            return
    else:
        mutation_matrix_path = MUTATION_MATRIX_PATH

    log(f"Loading mutation matrix: {mutation_matrix_path}")
    mut_matrix = pd.read_csv(mutation_matrix_path, sep='\t', index_col=0)
    log(f"  Shape: {mut_matrix.shape[0]} contexts × {mut_matrix.shape[1]} cells")

    # Select highest-SBS2 basal cells
    basal_mask = adata.obs['final_annotation'] == 'basal cell'
    basal_cells = adata.obs.loc[basal_mask]
    sbs2_threshold = basal_cells['SBS2'].quantile(0.80)
    high_sbs2_cells = basal_cells[basal_cells['SBS2'] >= sbs2_threshold].index.tolist()

    # Filter to cells present in mutation matrix
    available_cells = [c for c in high_sbs2_cells if c in mut_matrix.columns]
    log(f"  High-SBS2 basal cells (top 20%): {len(high_sbs2_cells)}")
    log(f"  Present in mutation matrix: {len(available_cells)}")

    if len(available_cells) == 0:
        log("  WARNING: No high-SBS2 cells found in mutation matrix. Skipping.")
        return

    # Aggregate mutations
    aggregated = mut_matrix[available_cells].sum(axis=1)
    total_mutations = int(aggregated.sum())
    log(f"  Total mutations in aggregated profile: {total_mutations:,}")

    # Normalize to proportions
    aggregated_norm = aggregated / aggregated.sum()

    # Align contexts
    if not all(aggregated_norm.index == cosmic_sbs2.index):
        log("  Reordering contexts to match COSMIC...")
        common_contexts = [c for c in cosmic_sbs2.index if c in aggregated_norm.index]
        aggregated_norm = aggregated_norm.loc[common_contexts]
        cosmic_sbs2_aligned = cosmic_sbs2.loc[common_contexts]
    else:
        cosmic_sbs2_aligned = cosmic_sbs2

    # Calculate cosine similarity
    from scipy.spatial.distance import cosine
    cos_sim = 1 - cosine(aggregated_norm.values, cosmic_sbs2_aligned.values)
    log(f"  Cosine similarity to COSMIC SBS2: {cos_sim:.4f}")

    # Extract mutation types for coloring
    mutation_types = [ctx.split('[')[1].split(']')[0] for ctx in aggregated_norm.index]
    colors = [COSMIC_COLORS[mt] for mt in mutation_types]

    # Create figure: two panels stacked
    fig, axes = plt.subplots(2, 1, figsize=(16, 8))

    x = np.arange(len(aggregated_norm))

    # Panel 1: Aggregated profile from high-SBS2 basal cells
    ax1 = axes[0]
    ax1.bar(x, aggregated_norm.values * 100, color=colors, width=0.8, edgecolor='none')
    for i in range(1, 6):
        ax1.axvline(i * 16 - 0.5, color='gray', linestyle='--', alpha=0.4, linewidth=0.5)
    ax1.set_ylabel('% of mutations', fontsize=12)
    ax1.set_title(
        f'Aggregated High-SBS2 Basal Cells ({total_mutations:,} mutations, '
        f'cosine similarity = {cos_sim:.3f})',
        fontsize=14, fontweight='bold'
    )
    ax1.set_xlim(-0.5, 95.5)
    ax1.set_xticks([])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Add mutation type labels with colored rectangles
    y_max_1 = ax1.get_ylim()[1]
    mutation_categories = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
    for i, mut in enumerate(mutation_categories):
        rect_x = i * 16
        rect = Rectangle((rect_x - 0.5, y_max_1), 16, 0.04 * y_max_1,
                          facecolor=COSMIC_COLORS[mut], clip_on=False)
        ax1.add_patch(rect)
        ax1.text(rect_x + 7.5, y_max_1 * 1.06, mut,
                 ha='center', va='bottom', fontsize=11, fontweight='bold')

    # Panel 2: COSMIC SBS2 reference
    ax2 = axes[1]
    ax2.bar(x, cosmic_sbs2_aligned.values * 100, color=colors, width=0.8, edgecolor='none')
    for i in range(1, 6):
        ax2.axvline(i * 16 - 0.5, color='gray', linestyle='--', alpha=0.4, linewidth=0.5)
    ax2.set_ylabel('% of mutations', fontsize=12)
    ax2.set_title('COSMIC SBS2 Reference (APOBEC)', fontsize=14, fontweight='bold')
    ax2.set_xlim(-0.5, 95.5)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # Add mutation type x-axis labels
    tick_positions = [8, 24, 40, 56, 72, 88]
    tick_labels = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
    ax2.set_xticks(tick_positions)
    ax2.set_xticklabels(tick_labels, fontsize=11)

    plt.tight_layout()

    # Save
    pdf_path = os.path.join(OUTPUT_DIR, "Panel_3c_Aggregated_SBS2_vs_COSMIC.pdf")
    png_path = os.path.join(OUTPUT_DIR, "Panel_3c_Aggregated_SBS2_vs_COSMIC.png")
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.savefig(png_path, dpi=DPI, bbox_inches='tight')
    plt.close()

    log(f"\n  Saved: {pdf_path}")
    log(f"  Saved: {png_path}")


# =============================================================================
# STEP 4: SUPPLEMENTAL FIGURE Y — A3 FAMILY EXPRESSION ACROSS CELL TYPES
# =============================================================================

def plot_supplemental_figure_y(adata):
    """
    Generate Supplemental Figure 2: UMAP for each A3 family member showing
    expression patterns across all cell types.
    """
    log("\n" + "=" * 80)
    log("STEP 4: GENERATING SUPPLEMENTAL FIGURE Y — A3 FAMILY EXPRESSION")
    log("=" * 80)

    # Determine grid layout for 7 A3 members (2 rows × 4 cols, last spot empty)
    n_genes = len(A3_FAMILY)
    ncols = 4
    nrows = 2

    fig, axes = plt.subplots(nrows, ncols, figsize=(24, 12))
    axes_flat = axes.flatten()

    for idx, (gene_full, gene_short) in enumerate(zip(A3_FAMILY, A3_SHORT)):
        ax = axes_flat[idx]

        # Extract expression
        expr = get_gene_expression(adata, gene_full)
        if expr is not None:
            adata.obs[f'_plot_{gene_short}'] = expr

            sc.pl.umap(
                adata,
                color=f'_plot_{gene_short}',
                cmap=CMAP_EXPRESSION,
                ax=ax,
                show=False,
                title=gene_short,
                frameon=False,
                size=8,
                colorbar_loc='right',
            )
            ax.set_title(gene_short, fontsize=22, fontweight='bold')

            # Log statistics
            nonzero = (expr > 0).sum()
            log(f"  {gene_short}: {nonzero:,} cells with expression "
                f"({100 * nonzero / len(expr):.1f}%)")
        else:
            ax.text(0.5, 0.5, f'{gene_short}\nnot found',
                    ha='center', va='center', fontsize=14, transform=ax.transAxes)
            ax.set_title(gene_short, fontsize=22, fontweight='bold')

        ax.set_xlabel('')
        ax.set_ylabel('')

    # Turn off unused axes
    for idx in range(n_genes, nrows * ncols):
        axes_flat[idx].axis('off')

    fig.suptitle('APOBEC3 Family Expression Across Single-Cell Populations',
                 fontsize=20, fontweight='bold', y=1.02)

    plt.tight_layout()

    # Save
    pdf_path = os.path.join(OUTPUT_DIR, "Supplemental_Figure_2_A3_Family_Expression.pdf")
    png_path = os.path.join(OUTPUT_DIR, "Supplemental_Figure_2_A3_Family_Expression.png")
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.savefig(png_path, dpi=DPI, bbox_inches='tight')
    plt.close()

    log(f"\n  Saved: {pdf_path}")
    log(f"  Saved: {png_path}")

    # Clean up temporary columns
    for gene_short in A3_SHORT:
        col = f'_plot_{gene_short}'
        if col in adata.obs.columns:
            del adata.obs[col]


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def main():
    """Run the complete Figure 3 panel generation pipeline."""

    log("=" * 80)
    log("FIGURE 3 PANEL GENERATION PIPELINE")
    log("=" * 80)
    log(f"Output directory: {OUTPUT_DIR}")
    log("")

    # Step 1: Load data
    adata = load_data()

    # Step 2: Figure 3b — Mutation and expression UMAPs
    plot_figure_3b(adata)

    # Step 3: Figure 3c — Aggregated signature comparison
    plot_figure_3c(adata)

    # Step 4: Supplemental Figure 2 — A3 family expression
    plot_supplemental_figure_y(adata)

    # Summary
    log("\n" + "=" * 80)
    log("FIGURE 3 PANEL GENERATION COMPLETE")
    log("=" * 80)
    log(f"\nOutput files in {OUTPUT_DIR}:")
    for f in sorted(os.listdir(OUTPUT_DIR)):
        fpath = os.path.join(OUTPUT_DIR, f)
        size = os.path.getsize(fpath)
        log(f"  {f} ({size / 1024:.1f} KB)")

    log("\nFigure 3 panel mapping:")
    log("  Panel a  — From ClusterCatcher: UMAP_clusters + stacked_bar + UMAP_final_annotation")
    log("  Panel b  — Panel_3b_Normalized_Mutations_SBS2_A3A_A3B.pdf")
    log("  Panel c  — Panel_3c_Aggregated_SBS2_vs_COSMIC.pdf")
    log("  Supp. 1  — From ClusterCatcher: UMAP_popv_prediction.pdf")
    log("  Supp. 2  — Supplemental_Figure_2_A3_Family_Expression.pdf")
    log("")


if __name__ == "__main__":
    main()
