#!/usr/bin/env python3
"""
Generate_Figure6_Panels.py
============================
Figure 6: HPV16 Lifecycle Genes Associate with Divergent Basal Cell
DNA Damage Profiles

Panels:
  6a — UMAP: basal cells colored by HPV16 raw UMI counts
  6b — UMAP: basal cells colored by early/late gene ratio (magma, gray for neg)
  6c — Violin plots: A3A, A3B, CNV, CytoTRACE2, SBS2, early/late ratio
        across SBS2-HIGH, Stealth CNV, Normal Control
  6d — Volcano: SBS2-HIGH vs Normal Control
  6e — Volcano: Stealth CNV vs Normal Control

Supplemental Figure 6:
  S6a — Volcano: SBS2-HIGH vs Stealth CNV
  S6b — GSEA dotplot: SBS2-HIGH vs Stealth CNV (spliceosome, ribosome, immune)

Inputs:
  - data/FIG_4/00_input/adata_final.h5ad
  - data/FIG_6/04_population_profiles_v2/revised_two_population_assignments.tsv
  - data/FIG_6/04_population_profiles_v2/DEG/normal_control_selection.tsv
  - data/FIG_6/04_population_profiles_v2/DEG/DEG_SBS2HIGH_vs_Normal.tsv
  - data/FIG_6/04_population_profiles_v2/DEG/DEG_StealthCNV_vs_Normal.tsv
  - data/FIG_6/04_population_profiles_v2/DEG/DEG_SBS2HIGH_vs_StealthCNV.tsv
  - data/FIG_6/04_population_profiles_v2/DEG/GSEA_SBS2HIGH_vs_StealthCNV.tsv
  - data/FIG_6/03_hpv16_genome/per_cell_hpv16_gene_counts.tsv
  - data/FIG_6/01_raw_hpv16_counts/basal_cell_master_table_with_raw_HPV16.tsv

Outputs (to data/FIG_6/FIGURE_6_PANELS/):
  - Panel_6a_UMAP_HPV16_UMI.pdf/.png
  - Panel_6b_UMAP_early_late_ratio.pdf/.png
  - Panel_6c_Violin_Population_Comparison.pdf/.png
  - Panel_6d_Volcano_SBS2HIGH_vs_Normal.pdf/.png
  - Panel_6e_Volcano_StealthCNV_vs_Normal.pdf/.png
  - Supplemental_6a_Volcano_SBS2HIGH_vs_StealthCNV.pdf/.png
  - Supplemental_6b_GSEA_Dotplot.pdf/.png

Env: NETWORK
Usage: conda run -n NETWORK python Generate_Figure6_Panels.py
Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
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
ADATA_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/adata_final.h5ad")
POP_PATH = os.path.join(PROJECT_ROOT, "data/FIG_6/04_population_profiles_v2/revised_two_population_assignments.tsv")
NORMAL_PATH = os.path.join(PROJECT_ROOT, "data/FIG_6/04_population_profiles_v2/DEG/normal_control_selection.tsv")
MASTER_PATH = os.path.join(PROJECT_ROOT, "data/FIG_6/01_raw_hpv16_counts/basal_cell_master_table_with_raw_HPV16.tsv")
HPV_GENE_PATH = os.path.join(PROJECT_ROOT, "data/FIG_6/03_hpv16_genome/per_cell_hpv16_gene_counts.tsv")
SIG_WEIGHTS_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/signature_weights_per_cell.txt")

# DEG/GSEA paths
DEG_DIR = os.path.join(PROJECT_ROOT, "data/FIG_6/04_population_profiles_v2/DEG")
DEG_HIGH_VS_NORMAL = os.path.join(DEG_DIR, "DEG_SBS2HIGH_vs_Normal.tsv")
DEG_STEALTH_VS_NORMAL = os.path.join(DEG_DIR, "DEG_StealthCNV_vs_Normal.tsv")
DEG_HIGH_VS_STEALTH = os.path.join(DEG_DIR, "DEG_SBS2HIGH_vs_StealthCNV.tsv")
GSEA_HIGH_VS_STEALTH = os.path.join(DEG_DIR, "GSEA_SBS2HIGH_vs_StealthCNV.tsv")

# Output
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_6/FIGURE_6_PANELS")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Style
DPI = 300
FONT_TITLE = 30
FONT_AXIS = 28
FONT_TICK = 22
FONT_LEGEND = 20
FONT_ANNOT = 18

# Population colors — consistent palette
COLOR_SBS2_HIGH = "#ed6a5a"       # coral
COLOR_STEALTH   = "#5b8e7d"       # muted teal
COLOR_NORMAL    = "#9bc1bc"        # light teal
COLOR_OTHER     = "#e0e0e0"        # light gray
POP_COLORS = {
    'SBS2_HIGH':      COLOR_SBS2_HIGH,
    'Stealth_CNV':    COLOR_STEALTH,
    'Normal_Control': COLOR_NORMAL,
}
POP_ORDER = ['SBS2_HIGH', 'Stealth_CNV', 'Normal_Control']
POP_LABELS = {
    'SBS2_HIGH':      'SBS2-HIGH',
    'Stealth_CNV':    'Stealth CNV',
    'Normal_Control': 'Normal Control',
}

# =============================================================================
# LOGGING
# =============================================================================
report_lines = []

def log(msg=""):
    print(msg, flush=True)
    report_lines.append(str(msg))

def banner(title, char="="):
    log("")
    log(char * 80)
    log(f"  {title}")
    log(char * 80)

def save_fig(fig, name):
    """Save figure as PDF and PNG."""
    for ext in ['pdf', 'png']:
        path = os.path.join(OUTPUT_DIR, f"{name}.{ext}")
        fig.savefig(path, dpi=DPI, bbox_inches='tight')
    log(f"  [SAVE] {name}.pdf/.png")
    plt.close(fig)

# =============================================================================
# STEP 0: LOAD ALL DATA
# =============================================================================
banner("STEP 0: Load all data")

# AnnData
log("  Loading adata_final.h5ad...")
adata = sc.read_h5ad(ADATA_PATH)
log(f"  adata: {adata.shape}")

# Subset to basal cells
basal_mask = adata.obs['final_annotation'].str.lower().str.contains('basal')
adata_basal = adata[basal_mask].copy()
log(f"  Basal cells: {adata_basal.shape[0]}")
if adata_basal.shape[0] == 0:
    log("  ERROR: No basal cells found. Unique annotations:")
    for ct, n in adata.obs['final_annotation'].value_counts().items():
        log(f"    {ct}: {n}")
    sys.exit(1)
del adata  # free memory

# Population assignments
pop = pd.read_csv(POP_PATH, sep='\t', index_col=0)
log(f"  Population assignments: {pop.shape}")

# Detect population column
pop_col = None
for candidate in ['population', 'group', 'Population', 'Group']:
    if candidate in pop.columns:
        pop_col = candidate
        break
if pop_col is None:
    for col in pop.columns:
        vals = pop[col].dropna().unique()
        if any('SBS2' in str(v) for v in vals):
            pop_col = col
            break
log(f"  Population column: '{pop_col}'")

# Normal control barcodes
normal_ctrl = pd.read_csv(NORMAL_PATH, sep='\t', index_col=0)
log(f"  Normal control cells: {len(normal_ctrl)}")

# Build group mapping for all basal cells
group_map = pd.Series('Other', index=adata_basal.obs_names)
if pop_col:
    for bc in pop.index:
        if bc in group_map.index:
            val = pop.loc[bc, pop_col]
            if pd.notna(val):
                group_map[bc] = str(val).replace(' ', '_').replace('-', '_')
for bc in normal_ctrl.index:
    if bc in group_map.index:
        group_map[bc] = 'Normal_Control'

adata_basal.obs['fig6_group'] = group_map
log(f"\n  Group counts:")
for g, c in adata_basal.obs['fig6_group'].value_counts().items():
    log(f"    {g}: {c}")

# Raw HPV16 counts
master = pd.read_csv(MASTER_PATH, sep='\t', index_col=0)
log(f"\n  Master table: {master.shape}")
adata_basal.obs['raw_HPV16'] = master['raw_HPV16'].reindex(adata_basal.obs_names).fillna(0).values

# HPV16 gene counts (early/late ratio)
hpv_genes = pd.read_csv(HPV_GENE_PATH, sep='\t', index_col=0)
log(f"  HPV16 gene table: {hpv_genes.shape}")

# Map early_late_ratio
if 'early_late_ratio' in hpv_genes.columns:
    adata_basal.obs['early_late_ratio'] = hpv_genes['early_late_ratio'].reindex(
        adata_basal.obs_names).fillna(0).values
else:
    # Compute from components
    early_cols = [c for c in ['E6', 'E7', 'E1', 'E2', 'E4', 'E5'] if c in hpv_genes.columns]
    late_cols = [c for c in ['L1', 'L2'] if c in hpv_genes.columns]
    early = hpv_genes[early_cols].sum(axis=1).reindex(adata_basal.obs_names).fillna(0)
    late = hpv_genes[late_cols].sum(axis=1).reindex(adata_basal.obs_names).fillna(0)
    adata_basal.obs['early_late_ratio'] = (early / (late + 0.5)).values

# SBS2 weight
sig_weights = pd.read_csv(SIG_WEIGHTS_PATH, sep='\t', index_col=0).T  # .T transposition!
log(f"  Signature weights (transposed): {sig_weights.shape}")
if 'SBS2' in sig_weights.columns:
    adata_basal.obs['SBS2_weight'] = sig_weights['SBS2'].reindex(adata_basal.obs_names).fillna(0).values
elif 'SBS2' in adata_basal.obs.columns:
    adata_basal.obs['SBS2_weight'] = adata_basal.obs['SBS2'].values
else:
    log("  WARNING: SBS2 weight not found, checking obs columns...")
    sbs2_candidates = [c for c in adata_basal.obs.columns if 'SBS2' in c.upper()]
    log(f"    Candidates: {sbs2_candidates}")
    if sbs2_candidates:
        adata_basal.obs['SBS2_weight'] = adata_basal.obs[sbs2_candidates[0]].values

# CNV and CytoTRACE2 should already be in obs
for col in ['cnv_score', 'CytoTRACE2_Score']:
    if col in adata_basal.obs.columns:
        log(f"  {col}: present, mean={adata_basal.obs[col].mean():.4f}")
    else:
        log(f"  WARNING: {col} not found in adata.obs")

# A3A and A3B expression
for gene in ['APOBEC3A', 'APOBEC3B']:
    if gene in adata_basal.var_names:
        gene_idx = list(adata_basal.var_names).index(gene)
        vals = adata_basal.X[:, gene_idx]
        if hasattr(vals, 'toarray'):
            vals = vals.toarray().flatten()
        else:
            vals = np.array(vals).flatten()
        adata_basal.obs[gene] = vals
        log(f"  {gene}: mean={vals.mean():.3f}")

log("\n  Data loading complete.")

# =============================================================================
# PANEL 6a: UMAP colored by HPV16 raw UMI counts
# =============================================================================
banner("PANEL 6a: UMAP — HPV16 raw UMI counts")

fig, ax = plt.subplots(figsize=(12, 10))

# Plot all basal cells, color by raw HPV16
hpv_vals = adata_basal.obs['raw_HPV16'].values.astype(float)
# Log-transform for visibility (many zeros)
hpv_log = np.log1p(hpv_vals)

# Sort so high values plot on top
order = np.argsort(hpv_log)
umap_coords = adata_basal.obsm['X_umap']

# Background (HPV-negative) in gray
neg_mask = hpv_vals == 0
ax.scatter(umap_coords[neg_mask, 0], umap_coords[neg_mask, 1],
           c='#e0e0e0', s=3, alpha=0.3, rasterized=True)

# HPV-positive colored
pos_mask = hpv_vals > 0
pos_order = np.argsort(hpv_log[pos_mask])
sc_plot = ax.scatter(
    umap_coords[pos_mask, 0][pos_order],
    umap_coords[pos_mask, 1][pos_order],
    c=hpv_log[pos_mask][pos_order],
    cmap='magma',
    s=5, alpha=0.7, rasterized=True,
    vmin=0, vmax=np.percentile(hpv_log[pos_mask], 99),
)

cbar = plt.colorbar(sc_plot, ax=ax, shrink=0.7, pad=0.02)
cbar.set_label('HPV16 UMI (log1p)', fontsize=FONT_AXIS, labelpad=10)
cbar.ax.tick_params(labelsize=FONT_TICK)

ax.set_title('HPV16 Viral Load in Basal Cells', fontsize=FONT_TITLE, pad=15)
ax.set_xlabel('UMAP 1', fontsize=FONT_AXIS)
ax.set_ylabel('UMAP 2', fontsize=FONT_AXIS)
ax.tick_params(labelsize=FONT_TICK)

# Count annotation
n_pos = pos_mask.sum()
n_neg = neg_mask.sum()
ax.text(0.02, 0.98, f'HPV16+ : {n_pos:,}\nHPV16− : {n_neg:,}',
        transform=ax.transAxes, fontsize=FONT_ANNOT,
        va='top', ha='left',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

plt.tight_layout()
save_fig(fig, "Panel_6a_UMAP_HPV16_UMI")

# =============================================================================
# PANEL 6b: UMAP colored by early/late gene ratio
# =============================================================================
banner("PANEL 6b: UMAP — Early/Late HPV16 gene ratio")

fig, ax = plt.subplots(figsize=(12, 10))

el_ratio = adata_basal.obs['early_late_ratio'].values.astype(float)
hpv_raw = adata_basal.obs['raw_HPV16'].values.astype(float)

# HPV-negative cells in gray
neg_mask = hpv_raw == 0
ax.scatter(umap_coords[neg_mask, 0], umap_coords[neg_mask, 1],
           c='#d0d0d0', s=3, alpha=0.3, rasterized=True, label='HPV16−')

# HPV-positive: color by log1p(early/late ratio)
pos_mask = hpv_raw > 0
el_log = np.log1p(el_ratio[pos_mask])
pos_order = np.argsort(el_log)

sc_plot = ax.scatter(
    umap_coords[pos_mask, 0][pos_order],
    umap_coords[pos_mask, 1][pos_order],
    c=el_log[pos_order],
    cmap='magma',
    s=5, alpha=0.7, rasterized=True,
    vmin=0, vmax=np.percentile(el_log, 99),
)

cbar = plt.colorbar(sc_plot, ax=ax, shrink=0.7, pad=0.02)
cbar.set_label('Early/Late Ratio (log1p)', fontsize=FONT_AXIS, labelpad=10)
cbar.ax.tick_params(labelsize=FONT_TICK)

ax.set_title('HPV16 Lifecycle Phase in Basal Cells', fontsize=FONT_TITLE, pad=15)
ax.set_xlabel('UMAP 1', fontsize=FONT_AXIS)
ax.set_ylabel('UMAP 2', fontsize=FONT_AXIS)
ax.tick_params(labelsize=FONT_TICK)

# Legend for gray
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#d0d0d0',
           markersize=10, label='HPV16−'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#fcfdbf',
           markersize=10, label='Late-dominant (low ratio)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#6a0d25',
           markersize=10, label='Early-dominant (high ratio)'),
]
ax.legend(handles=legend_elements, loc='lower right', fontsize=FONT_LEGEND,
          framealpha=0.9, markerscale=1.5)

plt.tight_layout()
save_fig(fig, "Panel_6b_UMAP_early_late_ratio")

# =============================================================================
# PANEL 6c: Violin plots across populations
# =============================================================================
banner("PANEL 6c: Violin plots — population comparison")

# Subset to the three defined groups only
mask_3groups = adata_basal.obs['fig6_group'].isin(POP_ORDER)
df_plot = adata_basal.obs.loc[mask_3groups].copy()
df_plot['Population'] = df_plot['fig6_group'].map(POP_LABELS)

# Define metrics to plot
metrics = [
    ('APOBEC3A',         'A3A Expression'),
    ('APOBEC3B',         'A3B Expression'),
    ('cnv_score',        'CNV Score'),
    ('CytoTRACE2_Score', 'CytoTRACE2\nStemness'),
    ('SBS2_weight',      'SBS2 Weight'),
    ('early_late_ratio', 'Early/Late\nRatio'),
]

# Filter to metrics that exist
metrics = [(col, label) for col, label in metrics if col in df_plot.columns]
n_metrics = len(metrics)

fig, axes = plt.subplots(1, n_metrics, figsize=(5 * n_metrics, 8))
if n_metrics == 1:
    axes = [axes]

pop_label_order = [POP_LABELS[p] for p in POP_ORDER]
palette = {POP_LABELS[k]: v for k, v in POP_COLORS.items()}

for i, (col, label) in enumerate(metrics):
    ax = axes[i]

    # Violin
    parts = sns.violinplot(
        data=df_plot, x='Population', y=col,
        order=pop_label_order, palette=palette,
        inner=None, linewidth=0.5, alpha=0.7, ax=ax,
        cut=0,
    )
    # Overlay boxplot
    sns.boxplot(
        data=df_plot, x='Population', y=col,
        order=pop_label_order, palette=palette,
        width=0.15, linewidth=1.5, fliersize=0, ax=ax,
        boxprops=dict(zorder=2),
        medianprops=dict(color='black', linewidth=2),
    )

    ax.set_ylabel(label, fontsize=FONT_AXIS)
    ax.set_xlabel('', fontsize=FONT_AXIS)
    ax.tick_params(axis='y', labelsize=FONT_TICK)
    ax.tick_params(axis='x', labelsize=FONT_TICK - 4, rotation=30)

    # Mann-Whitney p-values: SBS2-HIGH vs Stealth CNV
    g1 = df_plot.loc[df_plot['Population'] == POP_LABELS['SBS2_HIGH'], col].dropna()
    g2 = df_plot.loc[df_plot['Population'] == POP_LABELS['Stealth_CNV'], col].dropna()
    if len(g1) > 0 and len(g2) > 0:
        stat, pval = mannwhitneyu(g1, g2, alternative='two-sided')
        if pval < 1e-10:
            p_str = f"p<1e-10"
        elif pval < 0.001:
            p_str = f"p={pval:.1e}"
        else:
            p_str = f"p={pval:.3f}"

        y_max = df_plot[col].quantile(0.98)
        ax.text(0.5, 0.97, p_str, transform=ax.transAxes,
                fontsize=FONT_ANNOT - 2, ha='center', va='top',
                fontweight='bold')

    if i == 0:
        ax.set_title('Population Comparison', fontsize=FONT_TITLE, pad=15)

plt.tight_layout(w_pad=2)
save_fig(fig, "Panel_6c_Violin_Population_Comparison")

# =============================================================================
# VOLCANO PLOT HELPER
# =============================================================================
def plot_volcano(deg_path, title, save_name, n_label=15,
                 logfc_thresh=1.0, padj_thresh=0.05):
    """Generic volcano plot from DEG table."""
    banner(f"VOLCANO: {title}")

    deg = pd.read_csv(deg_path, sep='\t')
    log(f"  DEG table: {deg.shape}")
    log(f"  Columns: {list(deg.columns)}")

    # Detect column names
    logfc_col = None
    for c in ['logFC', 'log2FoldChange', 'logfoldchanges', 'log2_fold_change', 'avg_log2FC']:
        if c in deg.columns:
            logfc_col = c
            break
    pval_col = None
    for c in ['padj', 'pvals_adj', 'p_val_adj', 'FDR', 'adj_pval', 'Adjusted P-value']:
        if c in deg.columns:
            pval_col = c
            break
    gene_col = None
    for c in ['gene', 'names', 'Gene', 'gene_name', 'index']:
        if c in deg.columns:
            gene_col = c
            break
    if gene_col is None and deg.index.dtype == object:
        deg['gene'] = deg.index
        gene_col = 'gene'

    log(f"  logFC col: {logfc_col}, padj col: {pval_col}, gene col: {gene_col}")

    if logfc_col is None or pval_col is None:
        log("  ERROR: Cannot identify required columns. Skipping.")
        return

    deg = deg.dropna(subset=[logfc_col, pval_col])
    deg[pval_col] = deg[pval_col].clip(lower=1e-300)
    deg['_neglog10p'] = -np.log10(deg[pval_col])
    deg['_logFC'] = deg[logfc_col].astype(float)

    fig, ax = plt.subplots(figsize=(12, 10))

    # Color: significant up (red), significant down (blue), ns (gray)
    sig_up = (deg['_logFC'] > logfc_thresh) & (deg[pval_col] < padj_thresh)
    sig_down = (deg['_logFC'] < -logfc_thresh) & (deg[pval_col] < padj_thresh)
    ns = ~(sig_up | sig_down)

    ax.scatter(deg.loc[ns, '_logFC'], deg.loc[ns, '_neglog10p'],
               c='#d0d0d0', s=5, alpha=0.4, rasterized=True)
    ax.scatter(deg.loc[sig_up, '_logFC'], deg.loc[sig_up, '_neglog10p'],
               c=COLOR_SBS2_HIGH, s=8, alpha=0.6, rasterized=True, label=f'Up ({sig_up.sum()})')
    ax.scatter(deg.loc[sig_down, '_logFC'], deg.loc[sig_down, '_neglog10p'],
               c=COLOR_STEALTH, s=8, alpha=0.6, rasterized=True, label=f'Down ({sig_down.sum()})')

    # Label top genes
    sig_genes = deg[sig_up | sig_down].nlargest(n_label, '_neglog10p')
    texts = []
    for _, row in sig_genes.iterrows():
        ax.annotate(
            row[gene_col],
            xy=(row['_logFC'], row['_neglog10p']),
            fontsize=FONT_ANNOT - 4,
            fontweight='bold',
            ha='center', va='bottom',
            arrowprops=dict(arrowstyle='-', color='black', lw=0.5),
        )

    # Threshold lines
    ax.axhline(-np.log10(padj_thresh), color='black', ls='--', lw=0.8, alpha=0.5)
    ax.axvline(logfc_thresh, color='black', ls='--', lw=0.8, alpha=0.5)
    ax.axvline(-logfc_thresh, color='black', ls='--', lw=0.8, alpha=0.5)

    ax.set_xlabel('log₂ Fold Change', fontsize=FONT_AXIS)
    ax.set_ylabel('−log₁₀ Adjusted P-value', fontsize=FONT_AXIS)
    ax.set_title(title, fontsize=FONT_TITLE, pad=15)
    ax.tick_params(labelsize=FONT_TICK)
    ax.legend(fontsize=FONT_LEGEND, loc='upper right', framealpha=0.9, markerscale=2)

    plt.tight_layout()
    save_fig(fig, save_name)

# =============================================================================
# PANEL 6d: Volcano — SBS2-HIGH vs Normal
# =============================================================================
plot_volcano(
    DEG_HIGH_VS_NORMAL,
    'SBS2-HIGH vs Normal Control',
    'Panel_6d_Volcano_SBS2HIGH_vs_Normal'
)

# =============================================================================
# PANEL 6e: Volcano — Stealth CNV vs Normal
# =============================================================================
plot_volcano(
    DEG_STEALTH_VS_NORMAL,
    'Stealth CNV vs Normal Control',
    'Panel_6e_Volcano_StealthCNV_vs_Normal'
)

# =============================================================================
# SUPPLEMENTAL 6a: Volcano — SBS2-HIGH vs Stealth CNV
# =============================================================================
plot_volcano(
    DEG_HIGH_VS_STEALTH,
    'SBS2-HIGH vs Stealth CNV',
    'Supplemental_6a_Volcano_SBS2HIGH_vs_StealthCNV'
)

# =============================================================================
# SUPPLEMENTAL 6b: GSEA Dotplot — SBS2-HIGH vs Stealth CNV
# =============================================================================
banner("SUPPLEMENTAL 6b: GSEA Dotplot")

gsea = pd.read_csv(GSEA_HIGH_VS_STEALTH, sep='\t')
log(f"  GSEA table: {gsea.shape}")
log(f"  Columns: {list(gsea.columns)}")

# Detect columns
term_col = None
for c in ['Term', 'pathway', 'Name', 'Description', 'gene_set']:
    if c in gsea.columns:
        term_col = c
        break
nes_col = None
for c in ['NES', 'nes', 'normalized_enrichment_score']:
    if c in gsea.columns:
        nes_col = c
        break
fdr_col = None
for c in ['FDR q-val', 'fdr', 'FDR', 'padj', 'NOM p-val']:
    if c in gsea.columns:
        fdr_col = c
        break

log(f"  Term: {term_col}, NES: {nes_col}, FDR: {fdr_col}")

if term_col and nes_col and fdr_col:
    gsea = gsea.dropna(subset=[nes_col, fdr_col])
    gsea[fdr_col] = gsea[fdr_col].astype(float)
    gsea[nes_col] = gsea[nes_col].astype(float)

    # Select significant pathways, sorted by |NES|
    sig_gsea = gsea[gsea[fdr_col] < 0.05].copy()
    sig_gsea['_abs_nes'] = sig_gsea[nes_col].abs()
    sig_gsea = sig_gsea.sort_values('_abs_nes', ascending=True)

    # Highlight key pathways
    highlight_terms = ['spliceosome', 'ribosome', 'rna transport',
                       'allograft', 'autoimmune', 'antigen',
                       'proteasome', 'dna replication']

    # Take top 20 + any highlighted that are significant
    top20 = sig_gsea.tail(20)
    highlight_rows = sig_gsea[sig_gsea[term_col].str.lower().apply(
        lambda x: any(h in x for h in highlight_terms))]
    plot_gsea = pd.concat([top20, highlight_rows]).drop_duplicates(subset=[term_col])
    plot_gsea = plot_gsea.sort_values(nes_col)

    log(f"  Plotting {len(plot_gsea)} pathways")

    fig, ax = plt.subplots(figsize=(14, max(8, len(plot_gsea) * 0.45)))

    # Color by NES direction
    colors = [COLOR_STEALTH if x < 0 else COLOR_SBS2_HIGH for x in plot_gsea[nes_col]]
    sizes = 50 + 200 * (-np.log10(plot_gsea[fdr_col].clip(lower=1e-20)) / 20)

    # Clean term names (remove KEGG_ prefix, underscores)
    clean_terms = plot_gsea[term_col].str.replace('KEGG_', '').str.replace('_', ' ').str.title()

    ax.scatter(plot_gsea[nes_col], range(len(plot_gsea)),
               c=colors, s=sizes, alpha=0.8, edgecolors='black', linewidths=0.5)

    ax.set_yticks(range(len(plot_gsea)))
    ax.set_yticklabels(clean_terms, fontsize=FONT_TICK - 2)
    ax.set_xlabel('Normalized Enrichment Score (NES)', fontsize=FONT_AXIS)
    ax.set_title('KEGG GSEA: SBS2-HIGH vs Stealth CNV', fontsize=FONT_TITLE, pad=15)
    ax.tick_params(axis='x', labelsize=FONT_TICK)
    ax.axvline(0, color='black', ls='-', lw=1)

    # Direction labels
    ax.text(0.02, 0.02, '← Enriched in\nStealth CNV', transform=ax.transAxes,
            fontsize=FONT_ANNOT, color=COLOR_STEALTH, fontweight='bold', va='bottom')
    ax.text(0.98, 0.02, 'Enriched in\nSBS2-HIGH →', transform=ax.transAxes,
            fontsize=FONT_ANNOT, color=COLOR_SBS2_HIGH, fontweight='bold', va='bottom', ha='right')

    plt.tight_layout()
    save_fig(fig, "Supplemental_6b_GSEA_Dotplot")
else:
    log("  ERROR: Could not identify GSEA columns. Skipping dotplot.")

# =============================================================================
# SUMMARY
# =============================================================================
banner("FIGURE 6 GENERATION COMPLETE")

log(f"\nOutput directory: {OUTPUT_DIR}")
for f in sorted(os.listdir(OUTPUT_DIR)):
    fpath = os.path.join(OUTPUT_DIR, f)
    size = os.path.getsize(fpath)
    log(f"  {f} ({size / 1024:.1f} KB)")

log(f"\nFigure 6 panel mapping:")
log(f"  (a) Panel_6a_UMAP_HPV16_UMI — HPV16 viral load")
log(f"  (b) Panel_6b_UMAP_early_late_ratio — HPV16 lifecycle phase")
log(f"  (c) Panel_6c_Violin_Population_Comparison — A3A/A3B/CNV/stemness/SBS2/lifecycle")
log(f"  (d) Panel_6d_Volcano_SBS2HIGH_vs_Normal — DEG")
log(f"  (e) Panel_6e_Volcano_StealthCNV_vs_Normal — DEG")
log(f"  Supplemental 6a: Volcano SBS2-HIGH vs Stealth CNV")
log(f"  Supplemental 6b: GSEA dotplot")

# Save report
report_path = os.path.join(OUTPUT_DIR, "figure6_generation_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"\n  Report: {report_path}")
