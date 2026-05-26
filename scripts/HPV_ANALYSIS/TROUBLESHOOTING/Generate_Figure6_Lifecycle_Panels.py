#!/usr/bin/env python3
"""
Generate_Figure6_Lifecycle_Panels.py
======================================
Figure 6: HPV16 Lifecycle States Drive Divergent Mutagenic Programs

Panels:
  A — UMAP: epithelial cells colored by total HPV16 reads (log1p)
  B — Normalized HPV16 gene fractions by lifecycle phase (grouped bar)
  C — Dot plot: host marker genes (ATM/productive vs immune/differentiation)
  D — Split violins: A3A vs A3B expression by population

INPUTS:
  - data/FIG_4/01_group_selection/three_group_assignments.tsv
  - data/FIG_4/00_input/adata_final.h5ad
  - data/FIG_6/03_hpv16_genome/per_cell_hpv16_gene_counts.tsv

OUTPUTS (to data/FIG_6/FIGURE_6_PANELS/):
  - Panel_6A_UMAP_HPV16_viral_load.pdf/.png
  - Panel_6B_HPV_lifecycle_fractions.pdf/.png
  - Panel_6C_host_marker_dotplot.pdf/.png
  - Panel_6D_A3A_vs_A3B_violins.pdf/.png
  - Figure6_composite.pdf/.png

Env: NETWORK
Usage: conda run -n NETWORK python Generate_Figure6_Lifecycle_Panels.py

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import seaborn as sns
from scipy.stats import mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"

# Input paths
THREE_GROUP_PATH = os.path.join(PROJECT_ROOT,
    "data/FIG_4/01_group_selection/three_group_assignments.tsv")
ADATA_PATH = os.path.join(PROJECT_ROOT,
    "data/FIG_4/00_input/adata_final.h5ad")
HPV_GENE_PATH = os.path.join(PROJECT_ROOT,
    "data/FIG_6/03_hpv16_genome/per_cell_hpv16_gene_counts.tsv")

# Output
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_6/FIGURE_6_PANELS")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================================================
# STYLE — consistent with Figure 4
# =============================================================================
DPI = 300

# Population colors (from Generate_Figure4_Panels.py)
COLOR_SBS2_HIGH = "#ed6a5a"   # coral
COLOR_CNV_HIGH  = "#F6D155"   # mustard yellow
COLOR_NORMAL    = "#4682b4"   # steelblue
COLOR_OTHER     = "#e0e0e0"   # light gray

POP_COLORS = {
    'SBS2_HIGH': COLOR_SBS2_HIGH,
    'CNV_HIGH':  COLOR_CNV_HIGH,
    'NORMAL':    COLOR_NORMAL,
}
POP_ORDER = ['SBS2_HIGH', 'CNV_HIGH', 'NORMAL']
POP_LABELS = {
    'SBS2_HIGH': 'SBS2-HIGH',
    'CNV_HIGH':  'CNV-HIGH',
    'NORMAL':    'Normal',
}

# Lifecycle phase colors
PHASE_COLORS = {
    'Maintenance':   '#2166ac',   # blue
    'Amplification': '#92c5de',   # light blue
    'Oncogene':      '#f4a582',   # salmon
    'Capsid':        '#b2182b',   # dark red
}

# Font sizes (28-32 range, no panel titles)
FONT_AXIS   = 30
FONT_TICK   = 26
FONT_LEGEND = 22
FONT_ANNOT  = 24
FONT_LABEL  = 32   # panel labels (A, B, C, D)
FONT_GENE   = 20   # gene names in dot plot

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

# HPV16 gene groupings
HPV16_PHASES = {
    'Maintenance':   ['E1', 'E2'],
    'Amplification': ['E4', 'E5'],
    'Oncogene':      ['E6', 'E7'],
    'Capsid':        ['L1', 'L2'],
}

# Host marker genes for dot plot (Panel C)
# Left block: ATM/productive replication (enriched in CNV-HIGH)
# Right block: Immune/differentiation (enriched in SBS2-HIGH)
DOTPLOT_GENES_LEFT = {
    'label': 'Productive Replication / ATM',
    'genes': ['CHEK2', 'BRCA1', 'MRE11', 'H2AX', 'CDC25C', 'CDK1', 'CCNB1', 'MKI67', 'CGAS'],
}
DOTPLOT_GENES_RIGHT = {
    'label': 'Immune Visibility / Differentiation',
    'genes': ['HLA-A', 'HLA-B', 'HLA-C', 'B2M', 'TAP1', 'IVL', 'CDH1', 'CDKN1A', 'MDM2'],
}

# =============================================================================
# LOGGING
# =============================================================================
report_lines = []
def log(msg=""):
    print(msg, flush=True)
    report_lines.append(str(msg))

def banner(title):
    log("")
    log("=" * 80)
    log(f"  {title}")
    log("=" * 80)

def save_fig(fig, name):
    for ext in ['pdf', 'png']:
        path = os.path.join(OUTPUT_DIR, f"{name}.{ext}")
        fig.savefig(path, dpi=DPI, bbox_inches='tight')
    log(f"  [SAVE] {name}.pdf/.png")
    plt.close(fig)

# =============================================================================
# HELPER: Extract gene expression
# =============================================================================
def get_expression(adata, gene_symbol):
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

def p_to_stars(p):
    if p < 1e-4: return '****'
    if p < 1e-3: return '***'
    if p < 0.01: return '**'
    if p < 0.05: return '*'
    return 'ns'

# =============================================================================
# STEP 0: LOAD DATA
# =============================================================================
banner("STEP 0: Load data")

# Population assignments
groups = pd.read_csv(THREE_GROUP_PATH, sep='\t')
sbs2_cells = set(groups.loc[groups['group'] == 'SBS2_HIGH', 'cell_barcode'])
cnv_cells  = set(groups.loc[groups['group'] == 'CNV_HIGH',  'cell_barcode'])
normal_cells = set(groups.loc[groups['group'] == 'NORMAL',  'cell_barcode'])
log(f"  Populations: {len(sbs2_cells)} SBS2-HIGH, {len(cnv_cells)} CNV-HIGH, {len(normal_cells)} Normal")

# AnnData
log("  Loading adata_final.h5ad...")
adata = sc.read_h5ad(ADATA_PATH)
log(f"  adata: {adata.shape}")

# Tag populations
adata.obs['population'] = 'other'
adata.obs.loc[adata.obs_names.isin(sbs2_cells), 'population'] = 'SBS2_HIGH'
adata.obs.loc[adata.obs_names.isin(cnv_cells), 'population'] = 'CNV_HIGH'
adata.obs.loc[adata.obs_names.isin(normal_cells), 'population'] = 'NORMAL'

# Subset to basal/epithelial cells for UMAP
basal_mask = adata.obs['final_annotation'] == 'basal cell'
adata_basal = adata[basal_mask].copy()
log(f"  Epithelial cells (basal annotation): {adata_basal.shape[0]}")

# Population subset for panels C and D
pop_mask = adata.obs['population'].isin(POP_ORDER)
adata_pop = adata[pop_mask].copy()
log(f"  Cells in three populations: {adata_pop.shape[0]}")

# HPV16 gene counts
hpv_genes = pd.read_csv(HPV_GENE_PATH, sep='\t', index_col=0)
log(f"  HPV16 gene counts: {hpv_genes.shape}")

# Tag HPV populations
hpv_genes['population'] = 'other'
hpv_genes.loc[hpv_genes.index.isin(sbs2_cells), 'population'] = 'SBS2_HIGH'
hpv_genes.loc[hpv_genes.index.isin(cnv_cells), 'population'] = 'CNV_HIGH'
hpv_genes.loc[hpv_genes.index.isin(normal_cells), 'population'] = 'NORMAL'

# Merge HPV total reads into adata_basal for UMAP
total_col = 'total_hpv16_genome_reads'
if total_col in hpv_genes.columns:
    hpv_total = hpv_genes[total_col].reindex(adata_basal.obs_names).fillna(0)
    adata_basal.obs['total_hpv16_reads'] = hpv_total.values
    log(f"  HPV16 reads merged into adata_basal: {(hpv_total > 0).sum()} cells with reads")


# =============================================================================
# PANEL A: UMAP — HPV16 viral load
# =============================================================================
banner("PANEL A: UMAP — HPV16 viral load")

fig, ax = plt.subplots(figsize=(14, 12))

umap_coords = adata_basal.obsm['X_umap']
hpv_vals = adata_basal.obs['total_hpv16_reads'].values.astype(float)
hpv_log = np.log1p(hpv_vals)

# Background: HPV-negative cells in gray
neg_mask = hpv_vals == 0
ax.scatter(umap_coords[neg_mask, 0], umap_coords[neg_mask, 1],
           c=COLOR_OTHER, s=3, alpha=0.2, rasterized=True)

# HPV-positive cells colored by viral load
pos_mask = hpv_vals > 0
pos_order = np.argsort(hpv_log[pos_mask])
sc_plot = ax.scatter(
    umap_coords[pos_mask, 0][pos_order],
    umap_coords[pos_mask, 1][pos_order],
    c=hpv_log[pos_mask][pos_order],
    cmap='magma',
    s=6, alpha=0.7, rasterized=True,
    vmin=0, vmax=np.percentile(hpv_log[pos_mask], 99),
)

cbar = plt.colorbar(sc_plot, ax=ax, shrink=0.6, pad=0.02)
cbar.set_label('HPV16 reads (log1p)', fontsize=FONT_AXIS, labelpad=10)
cbar.ax.tick_params(labelsize=FONT_TICK)

# Population boundary indicators (outline key regions)
for pop, color in [('SBS2_HIGH', COLOR_SBS2_HIGH), ('CNV_HIGH', COLOR_CNV_HIGH)]:
    pop_cells = adata_basal.obs['population'] == pop
    if pop_cells.any():
        pop_umap = umap_coords[pop_cells.values]
        # Plot sparse boundary markers
        ax.scatter(pop_umap[:, 0], pop_umap[:, 1],
                   facecolors='none', edgecolors=color, s=30,
                   alpha=0.4, linewidths=0.8, rasterized=True)

# Counts annotation
n_pos = pos_mask.sum()
n_neg = neg_mask.sum()
ax.text(0.02, 0.98,
        f'HPV16+: {n_pos:,}\nHPV16\u2212: {n_neg:,}',
        transform=ax.transAxes, fontsize=FONT_ANNOT,
        va='top', ha='left',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

ax.set_xlabel('UMAP 1', fontsize=FONT_AXIS)
ax.set_ylabel('UMAP 2', fontsize=FONT_AXIS)
ax.tick_params(labelsize=FONT_TICK)
ax.set_frame_on(False)

# Legend for population rings
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='none',
           markeredgecolor=COLOR_SBS2_HIGH, markersize=12, markeredgewidth=2,
           label='SBS2-HIGH'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='none',
           markeredgecolor=COLOR_CNV_HIGH, markersize=12, markeredgewidth=2,
           label='CNV-HIGH'),
]
ax.legend(handles=legend_elements, loc='lower right', fontsize=FONT_LEGEND,
          framealpha=0.9)

plt.tight_layout()
save_fig(fig, "Panel_6A_UMAP_HPV16_viral_load")


# =============================================================================
# PANEL B: Normalized HPV gene fractions by lifecycle phase
# =============================================================================
banner("PANEL B: Normalized HPV gene fractions by lifecycle phase")

fig, ax = plt.subplots(figsize=(14, 10))

# Compute fractions for each population (only cells in populations)
bar_data = {}
for pop in POP_ORDER:
    mask = hpv_genes['population'] == pop
    sub = hpv_genes[mask]
    total = sub[total_col].values + 0.5  # pseudocount

    fracs = {}
    for phase, genes in HPV16_PHASES.items():
        available = [g for g in genes if g in sub.columns]
        if available:
            fracs[phase] = (sub[available].sum(axis=1).values / total).mean()
        else:
            fracs[phase] = 0.0
    bar_data[pop] = fracs

# Grouped bar chart
phases = list(HPV16_PHASES.keys())
x = np.arange(len(phases))
width = 0.25
offsets = [-width, 0, width]

for i, pop in enumerate(POP_ORDER):
    vals = [bar_data[pop][phase] for phase in phases]
    bars = ax.bar(x + offsets[i], vals, width,
                  color=POP_COLORS[pop], edgecolor='#333333', linewidth=0.8,
                  label=POP_LABELS[pop], zorder=3)

    # Value labels on bars
    for bar_rect, val in zip(bars, vals):
        if val > 0.005:
            ax.text(bar_rect.get_x() + bar_rect.get_width()/2, bar_rect.get_height() + 0.003,
                    f'{val:.3f}', ha='center', va='bottom', fontsize=FONT_GENE, fontweight='bold')

# Significance brackets for key comparisons (SBS2 vs CNV)
# Maintenance: SBS2 > CNV (p=6.4e-14)
# Capsid: CNV > SBS2 (p<1e-60)
bracket_data = [
    (0, 'Maintenance', 6.42e-14),
    (3, 'Capsid', 1.05e-60),
]
for idx, phase_name, p_val in bracket_data:
    sbs2_val = bar_data['SBS2_HIGH'][phase_name]
    cnv_val = bar_data['CNV_HIGH'][phase_name]
    max_val = max(sbs2_val, cnv_val)
    y_bracket = max_val + 0.02

    # Draw bracket between SBS2 and CNV bars
    x1 = idx + offsets[0]
    x2 = idx + offsets[1]
    ax.plot([x1, x1, x2, x2], [y_bracket, y_bracket + 0.005, y_bracket + 0.005, y_bracket],
            color='black', linewidth=1.5)
    ax.text((x1 + x2)/2, y_bracket + 0.007, p_to_stars(p_val),
            ha='center', va='bottom', fontsize=FONT_ANNOT, fontweight='bold')

ax.set_xticks(x)
ax.set_xticklabels(phases, fontsize=FONT_AXIS)
ax.set_ylabel('Fraction of total HPV16 reads', fontsize=FONT_AXIS)
ax.tick_params(axis='y', labelsize=FONT_TICK)
ax.legend(fontsize=FONT_LEGEND, loc='upper right', framealpha=0.9)
ax.set_ylim(bottom=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='y', alpha=0.3, zorder=0)

plt.tight_layout()
save_fig(fig, "Panel_6B_HPV_lifecycle_fractions")


# =============================================================================
# PANEL C: Dot plot — host marker genes
# =============================================================================
banner("PANEL C: Host marker gene dot plot")

all_dotplot_genes = DOTPLOT_GENES_LEFT['genes'] + DOTPLOT_GENES_RIGHT['genes']
n_left = len(DOTPLOT_GENES_LEFT['genes'])
n_right = len(DOTPLOT_GENES_RIGHT['genes'])
n_genes = len(all_dotplot_genes)

# Compute mean expression and percent expressing for each gene x population
dot_data = []
for gene in all_dotplot_genes:
    expr = get_expression(adata_pop, gene)
    if expr is None:
        log(f"  WARNING: {gene} not found, skipping")
        continue

    for pop in POP_ORDER:
        mask = adata_pop.obs['population'] == pop
        vals = expr[mask.values]
        mean_expr = np.mean(vals)
        pct_expr = 100 * np.sum(vals > 0) / max(len(vals), 1)
        dot_data.append({
            'gene': gene, 'population': pop,
            'mean_expression': mean_expr,
            'pct_expressing': pct_expr,
        })

dot_df = pd.DataFrame(dot_data)
if len(dot_df) == 0:
    log("  ERROR: No genes found for dot plot")
else:
    # Build dot plot manually
    genes_found = [g for g in all_dotplot_genes if g in dot_df['gene'].values]
    n_found = len(genes_found)

    fig, ax = plt.subplots(figsize=(max(18, n_found * 1.2), 8))

    # Normalize mean expression to 0-1 per gene for color mapping
    for gene in genes_found:
        gmask = dot_df['gene'] == gene
        vals = dot_df.loc[gmask, 'mean_expression']
        vmin, vmax = vals.min(), vals.max()
        if vmax > vmin:
            dot_df.loc[gmask, 'norm_expr'] = (vals - vmin) / (vmax - vmin)
        else:
            dot_df.loc[gmask, 'norm_expr'] = 0.5

    # Plot
    for j, pop in enumerate(POP_ORDER):
        pop_data = dot_df[dot_df['population'] == pop]
        for i, gene in enumerate(genes_found):
            row = pop_data[pop_data['gene'] == gene]
            if len(row) == 0:
                continue
            pct = row['pct_expressing'].values[0]
            mean_e = row['mean_expression'].values[0]

            # Size proportional to percent expressing
            size = max(pct * 4, 15)  # scale factor

            # Color: use population color with alpha based on expression
            ax.scatter(i, j, s=size, c=POP_COLORS[pop],
                       edgecolors='#333333', linewidths=0.5, alpha=0.9, zorder=3)

            # Add percent text below each dot
            if pct >= 5:
                ax.text(i, j - 0.35, f'{pct:.0f}%', ha='center', va='top',
                        fontsize=FONT_GENE - 4, color='#555555')

    # Axis formatting
    ax.set_xticks(range(len(genes_found)))
    ax.set_xticklabels(genes_found, fontsize=FONT_GENE, rotation=45, ha='right')
    ax.set_yticks(range(len(POP_ORDER)))
    ax.set_yticklabels([POP_LABELS[p] for p in POP_ORDER], fontsize=FONT_AXIS)

    # Separator between left and right gene blocks
    n_left_found = len([g for g in DOTPLOT_GENES_LEFT['genes'] if g in genes_found])
    if n_left_found > 0 and n_left_found < len(genes_found):
        ax.axvline(n_left_found - 0.5, color='#999999', linestyle='--', linewidth=1.5, alpha=0.7)
        # Block labels
        ax.text(n_left_found/2 - 0.5, len(POP_ORDER) + 0.1,
                DOTPLOT_GENES_LEFT['label'],
                ha='center', va='bottom', fontsize=FONT_ANNOT - 2, fontstyle='italic')
        ax.text(n_left_found + (len(genes_found) - n_left_found)/2 - 0.5, len(POP_ORDER) + 0.1,
                DOTPLOT_GENES_RIGHT['label'],
                ha='center', va='bottom', fontsize=FONT_ANNOT - 2, fontstyle='italic')

    # Size legend
    for pct_val, label in [(10, '10%'), (50, '50%'), (90, '90%')]:
        ax.scatter([], [], s=pct_val * 4, c='#888888', edgecolors='#333333',
                   linewidths=0.5, label=f'{label} expressing')
    ax.legend(fontsize=FONT_LEGEND, loc='upper right', framealpha=0.9,
              title='% cells expressing', title_fontsize=FONT_LEGEND)

    ax.set_ylim(-0.7, len(POP_ORDER) + 0.6)
    ax.set_xlim(-0.7, len(genes_found) - 0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='x', alpha=0.2)

    plt.tight_layout()
    save_fig(fig, "Panel_6C_host_marker_dotplot")


# =============================================================================
# PANEL D: A3A vs A3B split violins
# =============================================================================
banner("PANEL D: A3A vs A3B split violins")

fig, axes = plt.subplots(1, 2, figsize=(16, 10), sharey=True)

for ax_idx, (gene, ax) in enumerate(zip(['APOBEC3A', 'APOBEC3B'], axes)):
    expr = get_expression(adata_pop, gene)
    if expr is None:
        log(f"  WARNING: {gene} not found")
        continue

    # Build data for violin
    violin_data = []
    for pop in POP_ORDER:
        mask = adata_pop.obs['population'] == pop
        vals = expr[mask.values]
        for v in vals:
            violin_data.append({'population': POP_LABELS[pop], 'expression': v})

    vdf = pd.DataFrame(violin_data)

    # Violin plot
    parts = ax.violinplot(
        [vdf.loc[vdf['population'] == POP_LABELS[pop], 'expression'].values for pop in POP_ORDER],
        positions=range(len(POP_ORDER)),
        showmeans=False, showmedians=False, showextrema=False,
    )

    # Color each violin
    for pc, pop in zip(parts['bodies'], POP_ORDER):
        pc.set_facecolor(POP_COLORS[pop])
        pc.set_edgecolor('#333333')
        pc.set_linewidth(1.2)
        pc.set_alpha(0.8)

    # Add box plot overlay for median/IQR
    bp = ax.boxplot(
        [vdf.loc[vdf['population'] == POP_LABELS[pop], 'expression'].values for pop in POP_ORDER],
        positions=range(len(POP_ORDER)),
        widths=0.15, patch_artist=True, showfliers=False,
        medianprops=dict(color='black', linewidth=2),
        whiskerprops=dict(color='black', linewidth=1),
        capprops=dict(color='black', linewidth=1),
    )
    for patch, pop in zip(bp['boxes'], POP_ORDER):
        patch.set_facecolor('white')
        patch.set_edgecolor('#333333')
        patch.set_linewidth(1)
        patch.set_alpha(0.8)

    # Mean annotations
    for i, pop in enumerate(POP_ORDER):
        vals = vdf.loc[vdf['population'] == POP_LABELS[pop], 'expression'].values
        mean_val = np.mean(vals)
        ax.scatter(i, mean_val, color='black', s=50, zorder=5, marker='D')

    # Significance: SBS2 vs CNV
    v_sbs2 = vdf.loc[vdf['population'] == POP_LABELS['SBS2_HIGH'], 'expression'].values
    v_cnv = vdf.loc[vdf['population'] == POP_LABELS['CNV_HIGH'], 'expression'].values
    if len(v_sbs2) > 5 and len(v_cnv) > 5:
        _, p = mannwhitneyu(v_sbs2, v_cnv, alternative='two-sided')
        y_max = max(np.percentile(v_sbs2, 99), np.percentile(v_cnv, 99))
        y_bracket = y_max * 1.05
        ax.plot([0, 0, 1, 1], [y_bracket, y_bracket * 1.02, y_bracket * 1.02, y_bracket],
                color='black', linewidth=1.5)
        ax.text(0.5, y_bracket * 1.03, p_to_stars(p),
                ha='center', va='bottom', fontsize=FONT_ANNOT, fontweight='bold')

    # Gene name as label (not title)
    alias = 'A3A' if 'A' in gene[-1:] else 'A3B'
    ax.text(0.5, 0.97, alias, transform=ax.transAxes,
            fontsize=FONT_AXIS + 2, fontweight='bold', ha='center', va='top')

    ax.set_xticks(range(len(POP_ORDER)))
    ax.set_xticklabels([POP_LABELS[p] for p in POP_ORDER], fontsize=FONT_TICK, rotation=30, ha='right')
    ax.tick_params(axis='y', labelsize=FONT_TICK)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if ax_idx == 0:
        ax.set_ylabel('Expression', fontsize=FONT_AXIS)

plt.tight_layout()
save_fig(fig, "Panel_6D_A3A_vs_A3B_violins")


# =============================================================================
# COMPOSITE FIGURE
# =============================================================================
banner("COMPOSITE: Assembling Figure 6")

fig = plt.figure(figsize=(32, 28))
gs = gridspec.GridSpec(2, 2, hspace=0.3, wspace=0.25)

# Panel labels
panel_labels = ['A', 'B', 'C', 'D']
panel_files = [
    "Panel_6A_UMAP_HPV16_viral_load.png",
    "Panel_6B_HPV_lifecycle_fractions.png",
    "Panel_6C_host_marker_dotplot.png",
    "Panel_6D_A3A_vs_A3B_violins.png",
]

for idx, (label, fname) in enumerate(zip(panel_labels, panel_files)):
    ax = fig.add_subplot(gs[idx // 2, idx % 2])
    img_path = os.path.join(OUTPUT_DIR, fname)
    if os.path.exists(img_path):
        img = plt.imread(img_path)
        ax.imshow(img)
    ax.axis('off')
    ax.text(-0.02, 1.02, label, transform=ax.transAxes,
            fontsize=FONT_LABEL + 4, fontweight='bold', va='top', ha='right')

save_fig(fig, "Figure6_composite")


# =============================================================================
# SAVE REPORT
# =============================================================================
report_path = os.path.join(OUTPUT_DIR, "figure6_lifecycle_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"\n  Report: {report_path}")
log(f"  Output: {OUTPUT_DIR}")
banner("FIGURE 6 GENERATION COMPLETE")
