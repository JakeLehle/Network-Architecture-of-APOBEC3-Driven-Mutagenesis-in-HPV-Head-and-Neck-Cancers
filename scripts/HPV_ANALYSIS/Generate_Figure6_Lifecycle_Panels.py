#!/usr/bin/env python3
"""
Generate_Figure6_Lifecycle_Panels.py  (v5)
===========================================
Figure 6: HPV16 Lifecycle States Drive Divergent Mutagenic Programs

Layout:
  ROW 1:
    Panel A - UMAP: epithelial cells colored by total HPV16 reads (log1p)
    Panel B - Baseline violins: SBS2 weight, CNV score, CytoTRACE2,
              APOBEC3A, APOBEC3B across three populations

  ROW 2:
    Panel C - Host marker dot plot (51 genes, 7 biological categories, horizontal)

  ROW 3:
    Panel D - HPV16 read distribution violin (all cells, threshold line at 8 UMI)
    Panel E - Summary boxes: HPV16+ cell counts per group (>= 8 UMI)
    Panel F - HPV16 viral gene fraction violins by lifecycle phase

INPUTS:
  - data/FIG_4/01_group_selection/three_group_assignments.tsv
  - data/FIG_4/00_input/adata_final.h5ad
  - data/FIG_4/00_input/signature_weights_per_cell.txt
  - data/FIG_6/01_raw_hpv16_counts/basal_cell_master_table_with_raw_HPV16.tsv
  - data/FIG_6/03_hpv16_genome/per_cell_hpv16_gene_counts.tsv

OUTPUTS (to data/FIG_6/FIGURE_6_PANELS/):
  - Panel_6A_UMAP_HPV16_viral_load.pdf/.png
  - Panel_6B_baseline_violins.pdf/.png
  - Panel_6C_host_marker_dotplot.pdf/.png
  - Panel_6DEF_HPV_overview.pdf/.png
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
from matplotlib.patches import Patch, FancyBboxPatch
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
import seaborn as sns
from scipy.stats import mannwhitneyu
from collections import OrderedDict
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
SIG_WEIGHTS_PATH = os.path.join(PROJECT_ROOT,
    "data/FIG_4/00_input/signature_weights_per_cell.txt")
MASTER_TABLE_PATH = os.path.join(PROJECT_ROOT,
    "data/FIG_6/01_raw_hpv16_counts/basal_cell_master_table_with_raw_HPV16.tsv")
HPV_GENE_PATH = os.path.join(PROJECT_ROOT,
    "data/FIG_6/03_hpv16_genome/per_cell_hpv16_gene_counts.tsv")

# HPV16 positivity threshold (L-method, Phase3)
HPV16_THRESHOLD = 8

# Output
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_6/FIGURE_6_PANELS")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================================================
# STYLE
# =============================================================================
DPI = 300

# Population colors (consistent with Figure 4)
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
    'Maintenance':   '#2166ac',
    'Amplification': '#92c5de',
    'Oncogene':      '#f4a582',
    'Capsid':        '#b2182b',
}

# Font sizes (28-34 range)
FONT_AXIS   = 30
FONT_TICK   = 26
FONT_LEGEND = 22
FONT_ANNOT  = 20
FONT_LABEL  = 34
FONT_GENE   = 18
FONT_CAT    = 20
FONT_PHASE  = 28
FONT_PVAL   = 14
FONT_STARS  = 20
FONT_BOX    = 24    # summary box text

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

# HPV16 gene groupings by lifecycle phase
HPV16_PHASES = OrderedDict([
    ('Maintenance',   ['E1', 'E2']),
    ('Amplification', ['E4', 'E5']),
    ('Oncogene',      ['E6', 'E7']),
    ('Capsid',        ['L1', 'L2']),
])

# Host marker genes for dot plot (Panel C) -- 7 categories, 51 genes
# Left-to-right: A3 enzymes -> immune -> differentiation -> damage/proliferation
DOTPLOT_CATEGORIES = OrderedDict([
    ('APOBEC/\nInnate Immune', [
        'APOBEC3A', 'APOBEC3B', 'APOBEC3C',
        'APOBEC3D', 'APOBEC3F', 'APOBEC3G', 'APOBEC3H',
        'CGAS', 'STING1',
    ]),
    ('Immune\nSignaling', [
        'STAT1', 'HLA-A', 'HLA-B', 'HLA-C',
        'IRF1', 'TAP1', 'B2M',
    ]),
    ('Differentiation', [
        'KRT5', 'KRT14', 'KRT1', 'KRT10', 'CDH1', 'IVL',
    ]),
    ('ATM/DNA\nDamage', [
        'ATM', 'CHEK2', 'CHEK1', 'BRCA1',
        'MRE11', 'RAD50', 'NBN', 'H2AX',
        'STAT5A', 'STAT5B',
    ]),
    ('Transformation/\nProliferation', [
        'CDKN2A', 'MCM7', 'CCNE1', 'MKI67', 'TOP2A',
        'PCNA', 'BRD4', 'MED1', 'E2F1', 'E2F2',
    ]),
    ('p53/Rb\nPathway', [
        'CDKN1A', 'BAX', 'MDM2', 'RB1', 'TP53',
    ]),
    ('G2/M\nArrest', [
        'CDC25A', 'CDC25C', 'CDK1', 'CCNB1',
    ]),
])

# Minimum cells required to run Mann-Whitney for a group
MIN_CELLS_FOR_STATS = 10

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
# HELPER: Extract gene expression from adata
# =============================================================================
def get_expression(adata, gene_symbol):
    """Extract expression vector for a gene."""
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

def format_p(p):
    """Format p-value for display on plot."""
    if p < 1e-99:
        return 'p<1e-99'
    elif p < 0.001:
        return f'p={p:.0e}'
    elif p < 0.05:
        return f'p={p:.3f}'
    else:
        return f'p={p:.2f}'


def make_violin_panel(ax, data_dict, pop_order, pop_colors, pop_labels,
                      ylabel='', show_stats=True):
    """
    Draw violin + box overlay for one metric across populations.
    All three pairwise comparisons shown with stars and p-values.
    Groups with < MIN_CELLS_FOR_STATS cells are marked N.D.
    """
    positions = list(range(len(pop_order)))
    plot_data = [data_dict[p] for p in pop_order]

    # Check which groups have enough data for violins
    valid_groups = [len(data_dict[p]) >= MIN_CELLS_FOR_STATS for p in pop_order]

    # Violin bodies (only for groups with enough data)
    if any(valid_groups):
        parts = ax.violinplot(
            [d if v else [0] for d, v in zip(plot_data, valid_groups)],
            positions=positions,
            showmeans=False, showmedians=False, showextrema=False)
        for idx_pc, (pc, pop) in enumerate(zip(parts['bodies'], pop_order)):
            if valid_groups[idx_pc]:
                pc.set_facecolor(pop_colors[pop])
                pc.set_edgecolor('#333333')
                pc.set_linewidth(1.2)
                pc.set_alpha(0.8)
            else:
                pc.set_visible(False)

    # Box overlay (only for valid groups)
    valid_data = [d if v else [0] for d, v in zip(plot_data, valid_groups)]
    bp = ax.boxplot(valid_data, positions=positions,
                    widths=0.15, patch_artist=True, showfliers=False,
                    medianprops=dict(color='black', linewidth=2),
                    whiskerprops=dict(color='black', linewidth=1),
                    capprops=dict(color='black', linewidth=1))
    for idx_bp, patch in enumerate(bp['boxes']):
        if valid_groups[idx_bp]:
            patch.set_facecolor('white')
            patch.set_edgecolor('#333333')
            patch.set_linewidth(1)
            patch.set_alpha(0.8)
        else:
            patch.set_visible(False)
            # Hide whiskers/caps/medians for invalid groups
            for element in ['whiskers', 'caps', 'medians']:
                # Each group has 2 whiskers, 2 caps, 1 median
                n_per = 2 if element != 'medians' else 1
                start = idx_bp * n_per
                end = start + n_per
                for line in bp[element][start:end]:
                    line.set_visible(False)

    # Mean diamonds (valid groups only)
    for i, pop in enumerate(pop_order):
        if valid_groups[i] and len(data_dict[pop]) > 0:
            mean_val = np.mean(data_dict[pop])
            ax.scatter(i, mean_val, color='black', s=40, zorder=5, marker='D')
        elif not valid_groups[i]:
            # Mark as N.D.
            ax.text(i, 0, 'N.D.', ha='center', va='center',
                    fontsize=FONT_PVAL + 2, fontstyle='italic', color='#888888')

    # Pairwise stat brackets
    if show_stats and len(pop_order) >= 3:
        pairs = [(0, 1), (1, 2), (0, 2)]

        # Compute y range from valid data only
        valid_maxes = [np.percentile(data_dict[p], 99.5)
                       for p in pop_order if len(data_dict[p]) > 0]
        valid_mins = [np.min(data_dict[p])
                      for p in pop_order if len(data_dict[p]) > 0]
        if valid_maxes:
            y_data_max = max(valid_maxes)
            y_data_min = min(valid_mins)
        else:
            y_data_max = 1
            y_data_min = 0
        y_range = max(y_data_max - y_data_min, 0.001)
        bracket_gap = 0.12 * y_range

        for bracket_idx, (i, j) in enumerate(pairs):
            v1 = data_dict[pop_order[i]]
            v2 = data_dict[pop_order[j]]

            # Both groups must have enough data
            if not (valid_groups[i] and valid_groups[j]):
                # N.D. bracket
                y_bar = y_data_max + (bracket_idx + 0.5) * bracket_gap
                y_tip = y_bar + 0.015 * y_range
                ax.plot([i, i, j, j],
                        [y_bar, y_tip, y_tip, y_bar],
                        color='#cccccc', linewidth=1.0, clip_on=False)
                ax.text((i + j) / 2.0, y_tip + 0.005 * y_range,
                        'N.D.', ha='center', va='bottom',
                        fontsize=FONT_PVAL, fontstyle='italic', color='#888888')
                continue

            if len(v1) > 5 and len(v2) > 5:
                _, p = mannwhitneyu(v1, v2, alternative='two-sided')
                y_bar = y_data_max + (bracket_idx + 0.5) * bracket_gap
                y_tip = y_bar + 0.015 * y_range

                ax.plot([i, i, j, j],
                        [y_bar, y_tip, y_tip, y_bar],
                        color='black', linewidth=1.2, clip_on=False)

                stars = p_to_stars(p)
                p_text = format_p(p)
                ax.text((i + j) / 2.0, y_tip + 0.005 * y_range,
                        f'{stars}  {p_text}',
                        ha='center', va='bottom',
                        fontsize=FONT_PVAL, fontweight='bold')

        # Extend y-axis
        ax.set_ylim(bottom=y_data_min - 0.05 * y_range,
                    top=y_data_max + (len(pairs) + 1.2) * bracket_gap)

    ax.set_xticks(positions)
    ax.set_xticklabels([pop_labels[p] for p in pop_order],
                       fontsize=FONT_TICK - 6, rotation=35, ha='right')
    ax.tick_params(axis='y', labelsize=FONT_TICK - 6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=FONT_AXIS - 6)


# =============================================================================
# STEP 0: LOAD ALL DATA
# =============================================================================
banner("STEP 0: Load all data")

# --- Population assignments ---
groups = pd.read_csv(THREE_GROUP_PATH, sep='\t')
sbs2_cells = set(groups.loc[groups['group'] == 'SBS2_HIGH', 'cell_barcode'])
cnv_cells  = set(groups.loc[groups['group'] == 'CNV_HIGH',  'cell_barcode'])
normal_cells = set(groups.loc[groups['group'] == 'NORMAL',  'cell_barcode'])
log(f"  Populations: {len(sbs2_cells)} SBS2-HIGH, {len(cnv_cells)} CNV-HIGH, "
    f"{len(normal_cells)} Normal")

# Cell-to-group lookup
cell_to_group = {}
for _, row in groups.iterrows():
    cell_to_group[row['cell_barcode']] = row['group']

# --- AnnData ---
log("  Loading adata_final.h5ad...")
adata = sc.read_h5ad(ADATA_PATH)
log(f"  adata: {adata.shape}")

adata.obs['population'] = 'other'
adata.obs.loc[adata.obs_names.isin(sbs2_cells), 'population'] = 'SBS2_HIGH'
adata.obs.loc[adata.obs_names.isin(cnv_cells), 'population'] = 'CNV_HIGH'
adata.obs.loc[adata.obs_names.isin(normal_cells), 'population'] = 'NORMAL'

basal_mask = adata.obs['final_annotation'] == 'basal cell'
adata_basal = adata[basal_mask].copy()
log(f"  Epithelial cells (basal annotation): {adata_basal.shape[0]}")

pop_mask = adata.obs['population'].isin(POP_ORDER)
adata_pop = adata[pop_mask].copy()
log(f"  Cells in three populations: {adata_pop.shape[0]}")

# --- Signature weights (for Panel B) ---
log(f"  Loading signature weights: {SIG_WEIGHTS_PATH}")
sig_weights = pd.read_csv(SIG_WEIGHTS_PATH, sep='\t', index_col=0)
log(f"  Signature weights shape: {sig_weights.shape}")

sbs2_col = None
for candidate in sig_weights.columns:
    if str(candidate).strip() == 'SBS2':
        sbs2_col = candidate
        break
if sbs2_col is None:
    for candidate in sig_weights.columns:
        if 'SBS2' in str(candidate):
            sbs2_col = candidate
            break
if sbs2_col is None and 'SBS2' in sig_weights.index:
    log("  Signature weights appear transposed, transposing...")
    sig_weights = sig_weights.T
    sbs2_col = 'SBS2'
log(f"  SBS2 column: '{sbs2_col}'")

# --- Master table (for raw HPV16 counts, Panel D/E) ---
log(f"  Loading master table: {MASTER_TABLE_PATH}")
master = pd.read_csv(MASTER_TABLE_PATH, sep='\t', index_col=0)
log(f"  Master table: {master.shape}")

# Tag master table with group labels
master['group'] = master.index.map(lambda x: cell_to_group.get(x, 'other'))
master_pop = master[master['group'].isin(POP_ORDER)].copy()
log(f"  Master table cells in three populations: {len(master_pop)}")
for pop in POP_ORDER:
    n = (master_pop['group'] == pop).sum()
    n_pos = ((master_pop['group'] == pop) &
             (master_pop['raw_HPV16'] >= HPV16_THRESHOLD)).sum()
    log(f"    {pop}: {n} total, {n_pos} HPV16+ (>= {HPV16_THRESHOLD} UMI)")

# --- HPV16 gene counts (Phase4 alignment, for Panel F) ---
hpv_genes = pd.read_csv(HPV_GENE_PATH, sep='\t', index_col=0)
log(f"  HPV16 gene counts: {hpv_genes.shape}")

hpv_genes['population'] = 'other'
hpv_genes.loc[hpv_genes.index.isin(sbs2_cells), 'population'] = 'SBS2_HIGH'
hpv_genes.loc[hpv_genes.index.isin(cnv_cells), 'population'] = 'CNV_HIGH'
hpv_genes.loc[hpv_genes.index.isin(normal_cells), 'population'] = 'NORMAL'

total_col = 'total_hpv16_genome_reads'
if total_col in hpv_genes.columns:
    hpv_total = hpv_genes[total_col].reindex(adata_basal.obs_names).fillna(0)
    adata_basal.obs['total_hpv16_reads'] = hpv_total.values
    log(f"  HPV16 reads merged into adata_basal: "
        f"{(hpv_total > 0).sum()} cells with reads")


# #############################################################################
#
#   PANEL A: UMAP -- HPV16 viral load
#
# #############################################################################
banner("PANEL A: UMAP -- HPV16 viral load")

fig, ax = plt.subplots(figsize=(14, 12))

umap_coords = adata_basal.obsm['X_umap']
hpv_vals = adata_basal.obs['total_hpv16_reads'].values.astype(float)
hpv_log = np.log1p(hpv_vals)

neg_mask = hpv_vals == 0
ax.scatter(umap_coords[neg_mask, 0], umap_coords[neg_mask, 1],
           c=COLOR_OTHER, s=5, alpha=0.15, rasterized=True)

pos_mask = hpv_vals > 0
pos_order = np.argsort(hpv_log[pos_mask])
sc_plot = ax.scatter(
    umap_coords[pos_mask, 0][pos_order],
    umap_coords[pos_mask, 1][pos_order],
    c=hpv_log[pos_mask][pos_order],
    cmap='magma', s=10, alpha=0.8, rasterized=True,
    vmin=0, vmax=np.percentile(hpv_log[pos_mask], 99),
)

cbar = plt.colorbar(sc_plot, ax=ax, shrink=0.6, pad=0.02)
cbar.set_label('HPV16 reads (log1p)', fontsize=FONT_AXIS, labelpad=10)
cbar.ax.tick_params(labelsize=FONT_TICK)

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
plt.tight_layout()
save_fig(fig, "Panel_6A_UMAP_HPV16_viral_load")


# #############################################################################
#
#   PANEL B: Baseline violins (SBS2, CNV, CytoTRACE, A3A, A3B)
#
# #############################################################################
banner("PANEL B: Baseline violins")

fig, axes = plt.subplots(1, 5, figsize=(40, 12))
plt.subplots_adjust(wspace=0.45)

# --- Metric 1: SBS2 signature weight ---
log("  Metric 1: SBS2 signature weight")
sbs2_data = {}
for pop in POP_ORDER:
    cells = sbs2_cells if pop == 'SBS2_HIGH' else (cnv_cells if pop == 'CNV_HIGH' else normal_cells)
    overlap = cells & set(sig_weights.index)
    if sbs2_col and len(overlap) > 0:
        sbs2_data[pop] = sig_weights.loc[list(overlap), sbs2_col].values.astype(float)
    else:
        sbs2_data[pop] = np.array([0.0])
    log(f"    {pop}: n={len(sbs2_data[pop])}, mean={np.mean(sbs2_data[pop]):.4f}")

make_violin_panel(axes[0], sbs2_data, POP_ORDER, POP_COLORS, POP_LABELS,
                  ylabel='SBS2 Weight')
axes[0].set_title('SBS2', fontsize=FONT_AXIS, fontweight='bold', pad=15)

# --- Metric 2: inferCNV score ---
log("  Metric 2: inferCNV score")
cnv_data = {}
for pop in POP_ORDER:
    mask = adata_pop.obs['population'] == pop
    cnv_data[pop] = adata_pop.obs.loc[mask, 'cnv_score'].values.astype(float)
    log(f"    {pop}: n={len(cnv_data[pop])}, mean={np.mean(cnv_data[pop]):.4f}")

make_violin_panel(axes[1], cnv_data, POP_ORDER, POP_COLORS, POP_LABELS,
                  ylabel='CNV Score')
axes[1].set_title('inferCNV', fontsize=FONT_AXIS, fontweight='bold', pad=15)

# --- Metric 3: CytoTRACE2 ---
log("  Metric 3: CytoTRACE2 score")
cyto_data = {}
for pop in POP_ORDER:
    mask = adata_pop.obs['population'] == pop
    cyto_data[pop] = adata_pop.obs.loc[mask, 'CytoTRACE2_Score'].values.astype(float)
    log(f"    {pop}: n={len(cyto_data[pop])}, mean={np.mean(cyto_data[pop]):.4f}")

make_violin_panel(axes[2], cyto_data, POP_ORDER, POP_COLORS, POP_LABELS,
                  ylabel='CytoTRACE2 Score')
axes[2].set_title('Stemness', fontsize=FONT_AXIS, fontweight='bold', pad=15)

# --- Metric 4: APOBEC3A ---
log("  Metric 4: APOBEC3A")
a3a_data = {}
a3a_expr = get_expression(adata_pop, 'APOBEC3A')
if a3a_expr is not None:
    for pop in POP_ORDER:
        mask = (adata_pop.obs['population'] == pop).values
        a3a_data[pop] = a3a_expr[mask]
        log(f"    {pop}: n={len(a3a_data[pop])}, mean={np.mean(a3a_data[pop]):.4f}")
else:
    log("  WARNING: APOBEC3A not found")
    for pop in POP_ORDER:
        a3a_data[pop] = np.array([0.0])

make_violin_panel(axes[3], a3a_data, POP_ORDER, POP_COLORS, POP_LABELS,
                  ylabel='Expression')
axes[3].set_title('A3A', fontsize=FONT_AXIS, fontweight='bold', pad=15)

# --- Metric 5: APOBEC3B ---
log("  Metric 5: APOBEC3B")
a3b_data = {}
a3b_expr = get_expression(adata_pop, 'APOBEC3B')
if a3b_expr is not None:
    for pop in POP_ORDER:
        mask = (adata_pop.obs['population'] == pop).values
        a3b_data[pop] = a3b_expr[mask]
        log(f"    {pop}: n={len(a3b_data[pop])}, mean={np.mean(a3b_data[pop]):.4f}")
else:
    log("  WARNING: APOBEC3B not found")
    for pop in POP_ORDER:
        a3b_data[pop] = np.array([0.0])

make_violin_panel(axes[4], a3b_data, POP_ORDER, POP_COLORS, POP_LABELS,
                  ylabel='Expression')
axes[4].set_title('A3B', fontsize=FONT_AXIS, fontweight='bold', pad=15)

legend_elements = [
    Patch(facecolor=POP_COLORS[p], edgecolor='#333333', label=POP_LABELS[p])
    for p in POP_ORDER
]
fig.legend(handles=legend_elements, loc='lower center',
           ncol=3, fontsize=FONT_LEGEND + 2, framealpha=0.9,
           bbox_to_anchor=(0.5, -0.01))
plt.subplots_adjust(bottom=0.15)
save_fig(fig, "Panel_6B_baseline_violins")


# #############################################################################
#
#   PANEL C: Host marker dot plot (51 genes, 7 categories, HORIZONTAL)
#
# #############################################################################
banner("PANEL C: Host marker gene dot plot")

all_genes = []
cat_boundaries = []
for cat_label, genes in DOTPLOT_CATEGORIES.items():
    start = len(all_genes)
    all_genes.extend(genes)
    end = len(all_genes)
    cat_boundaries.append((start, end, cat_label))

log(f"  Total genes in dot plot: {len(all_genes)}")

dot_records = []
genes_found = []
genes_missing = []

for gene in all_genes:
    expr = get_expression(adata_pop, gene)
    if expr is None:
        log(f"  WARNING: {gene} not found, skipping")
        genes_missing.append(gene)
        continue
    genes_found.append(gene)
    for pop in POP_ORDER:
        mask = (adata_pop.obs['population'] == pop).values
        vals = expr[mask]
        mean_expr = np.mean(vals)
        pct_expr = 100.0 * np.sum(vals > 0) / max(len(vals), 1)
        dot_records.append({
            'gene': gene, 'population': pop,
            'mean_expression': mean_expr, 'pct_expressing': pct_expr,
        })

dot_df = pd.DataFrame(dot_records)
log(f"  Genes found: {len(genes_found)}, missing: {len(genes_missing)}")
if genes_missing:
    log(f"  Missing: {genes_missing}")

n_found = len(genes_found)
n_pops = len(POP_ORDER)

fig_width = max(36, n_found * 0.72 + 6)
fig, ax = plt.subplots(figsize=(fig_width, 7))

cmap = plt.cm.YlOrRd

gene_min = {}
gene_max = {}
for gene in genes_found:
    gmask = dot_df['gene'] == gene
    vals = dot_df.loc[gmask, 'mean_expression']
    gene_min[gene] = vals.min()
    gene_max[gene] = vals.max()

SIZE_MIN = 40
SIZE_MAX = 400
PCT_MAX = 100.0

for j, pop in enumerate(POP_ORDER):
    pop_data = dot_df[dot_df['population'] == pop]
    for i, gene in enumerate(genes_found):
        row = pop_data[pop_data['gene'] == gene]
        if len(row) == 0:
            continue
        pct = row['pct_expressing'].values[0]
        mean_e = row['mean_expression'].values[0]
        size = SIZE_MIN + (pct / PCT_MAX) * (SIZE_MAX - SIZE_MIN)
        gmin = gene_min[gene]
        gmax = gene_max[gene]
        norm_val = (mean_e - gmin) / (gmax - gmin) if gmax > gmin else 0.5
        color = cmap(norm_val)
        ax.scatter(i, j, s=size, c=[color],
                   edgecolors='#555555', linewidths=0.6, zorder=3)

for start_orig, end_orig, cat_label in cat_boundaries:
    cat_genes_in_order = [g for g in all_genes[start_orig:end_orig]
                          if g in genes_found]
    if not cat_genes_in_order:
        continue
    first_idx = genes_found.index(cat_genes_in_order[0])
    last_idx = genes_found.index(cat_genes_in_order[-1])
    if first_idx > 0:
        ax.axvline(first_idx - 0.5, color='#bbbbbb', linestyle='-',
                   linewidth=1.2, alpha=0.7)
    mid_x = (first_idx + last_idx) / 2.0
    ax.text(mid_x, n_pops - 0.4, cat_label,
            fontsize=FONT_CAT, fontstyle='italic', va='bottom', ha='center',
            color='#333333')

ax.set_xticks(range(n_found))
ax.set_xticklabels(genes_found, fontsize=FONT_GENE, rotation=55, ha='right')
ax.set_yticks(range(n_pops))
ax.set_yticklabels([POP_LABELS[p] for p in POP_ORDER],
                   fontsize=FONT_AXIS - 2, fontweight='bold')
ax.set_xlim(-0.7, n_found - 0.3)
ax.set_ylim(-0.6, n_pops + 0.8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

legend_pcts = [10, 25, 50, 75, 100]
size_legend_elements = []
for pct in legend_pcts:
    s = SIZE_MIN + (pct / PCT_MAX) * (SIZE_MAX - SIZE_MIN)
    size_legend_elements.append(
        ax.scatter([], [], s=s, c='#888888', edgecolors='#555555',
                   linewidths=0.6, label=f'{pct}%'))
leg1 = ax.legend(handles=size_legend_elements,
                 title='% expressing', title_fontsize=FONT_LEGEND - 2,
                 fontsize=FONT_LEGEND - 4, loc='lower right',
                 framealpha=0.9, borderpad=1.0)
ax.add_artist(leg1)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(0, 1))
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.01, aspect=15,
                    orientation='vertical')
cbar.set_label('Scaled expression', fontsize=FONT_LEGEND)
cbar.set_ticks([0, 0.5, 1])
cbar.set_ticklabels(['Low', 'Mid', 'High'])
cbar.ax.tick_params(labelsize=FONT_LEGEND - 4)

plt.tight_layout()
save_fig(fig, "Panel_6C_host_marker_dotplot")


# #############################################################################
#
#   PANELS D + E + F: HPV overview row (combined figure)
#   D = read distribution violin with threshold
#   E = summary boxes (HPV16+ counts)
#   F = lifecycle fraction violins
#
# #############################################################################
banner("PANELS D/E/F: HPV16 overview row")

# Width ratios: D (1 violin set) : E (summary boxes) : F (8 violins)
fig = plt.figure(figsize=(48, 12))
gs = gridspec.GridSpec(1, 3, width_ratios=[1.2, 0.8, 5.0],
                       wspace=0.25)

# -------------------------------------------------------------------------
# PANEL D: HPV16 read distribution violin with threshold line
# -------------------------------------------------------------------------
log("  Panel D: HPV16 read distribution")
ax_d = fig.add_subplot(gs[0, 0])

# Get raw HPV16 counts for each group (all cells, including zeros)
hpv_dist_data = {}
for pop in POP_ORDER:
    mask = master_pop['group'] == pop
    vals = master_pop.loc[mask, 'raw_HPV16'].values.astype(float)
    hpv_dist_data[pop] = np.log1p(vals)
    n_pos = (vals >= HPV16_THRESHOLD).sum()
    log(f"    {pop}: n={len(vals)}, HPV16+ (>={HPV16_THRESHOLD}): {n_pos}, "
        f"mean raw={np.mean(vals):.1f}")

make_violin_panel(ax_d, hpv_dist_data, POP_ORDER, POP_COLORS, POP_LABELS,
                  ylabel='HPV16 UMI (log1p)', show_stats=True)

# Threshold line
thresh_log = np.log1p(HPV16_THRESHOLD)
ax_d.axhline(thresh_log, color='#333333', linestyle='--', linewidth=2.0,
             zorder=4, alpha=0.8)
# Label the threshold on the right side of the plot
ax_d.text(2.6, thresh_log, f'L-method\nthreshold\n({HPV16_THRESHOLD} UMI)',
          fontsize=FONT_PVAL, va='center', ha='left',
          color='#333333', fontstyle='italic')

ax_d.set_title('HPV16 Read\nDistribution', fontsize=FONT_AXIS - 2,
               fontweight='bold', pad=15)

# -------------------------------------------------------------------------
# PANEL E: Summary boxes (HPV16+ cell counts per group)
# -------------------------------------------------------------------------
log("  Panel E: HPV16+ summary boxes")
ax_e = fig.add_subplot(gs[0, 1])
ax_e.set_xlim(0, 1)
ax_e.set_ylim(0, 1)
ax_e.axis('off')

ax_e.text(0.5, 0.97, f'HPV16+ Cells\n(\u2265 {HPV16_THRESHOLD} UMI)',
          fontsize=FONT_BOX, fontweight='bold', ha='center', va='top',
          transform=ax_e.transAxes)

box_height = 0.22
box_width = 0.85
y_positions = [0.65, 0.40, 0.15]

for idx, pop in enumerate(POP_ORDER):
    mask = master_pop['group'] == pop
    n_total = mask.sum()
    n_pos = ((master_pop['group'] == pop) &
             (master_pop['raw_HPV16'] >= HPV16_THRESHOLD)).sum()
    pct = 100.0 * n_pos / n_total if n_total > 0 else 0.0

    y = y_positions[idx]
    x = (1 - box_width) / 2.0

    # Colored box
    rect = FancyBboxPatch((x, y), box_width, box_height,
                          boxstyle="round,pad=0.02",
                          facecolor=POP_COLORS[pop], edgecolor='#333333',
                          linewidth=1.5, alpha=0.85,
                          transform=ax_e.transAxes)
    ax_e.add_patch(rect)

    # Text inside
    ax_e.text(0.5, y + box_height * 0.6,
              POP_LABELS[pop], fontsize=FONT_BOX - 2,
              fontweight='bold', ha='center', va='center',
              transform=ax_e.transAxes, color='black')
    ax_e.text(0.5, y + box_height * 0.25,
              f'{n_pos} / {n_total}  ({pct:.1f}%)',
              fontsize=FONT_BOX - 4, ha='center', va='center',
              transform=ax_e.transAxes, color='#333333')

# -------------------------------------------------------------------------
# PANEL F: HPV16 lifecycle fraction violins
# -------------------------------------------------------------------------
log("  Panel F: HPV16 lifecycle fraction violins")

# Subset to HPV+ cells in the three populations
hpv_pop = hpv_genes[hpv_genes['population'].isin(POP_ORDER)].copy()
hpv_pop = hpv_pop[hpv_pop[total_col] > 0].copy()
log(f"  HPV+ cells (alignment-based) in three populations: {len(hpv_pop)}")
for pop in POP_ORDER:
    n = (hpv_pop['population'] == pop).sum()
    log(f"    {pop}: {n}")

# Compute per-cell fractions
all_hpv_genes = []
for phase_genes in HPV16_PHASES.values():
    all_hpv_genes.extend(phase_genes)

for gene in all_hpv_genes:
    if gene in hpv_pop.columns:
        hpv_pop[f'{gene}_frac'] = hpv_pop[gene] / hpv_pop[total_col]
    else:
        log(f"  WARNING: {gene} not in HPV count columns")
        hpv_pop[f'{gene}_frac'] = 0.0

# Nested gridspec for 8 violins within Panel F region
gs_f = gridspec.GridSpecFromSubplotSpec(1, 8, subplot_spec=gs[0, 2],
                                        wspace=0.40)

f_axes = []
gene_idx = 0
for phase_name, phase_genes in HPV16_PHASES.items():
    for gene in phase_genes:
        ax = fig.add_subplot(gs_f[0, gene_idx])
        f_axes.append(ax)
        frac_col = f'{gene}_frac'

        vdata = {}
        for pop in POP_ORDER:
            mask = hpv_pop['population'] == pop
            vdata[pop] = hpv_pop.loc[mask, frac_col].values.astype(float)

        make_violin_panel(ax, vdata, POP_ORDER, POP_COLORS, POP_LABELS,
                          ylabel='Fraction of HPV16 reads' if gene_idx == 0 else '',
                          show_stats=True)

        ax.set_title(gene, fontsize=FONT_AXIS - 2, fontweight='bold', pad=15)
        gene_idx += 1

# Phase headers above pairs
# Need positions after layout is computed
fig.canvas.draw()
phase_pairs = [(0, 1), (2, 3), (4, 5), (6, 7)]
for (left, right), (phase_name, _) in zip(phase_pairs, HPV16_PHASES.items()):
    pos_l = f_axes[left].get_position()
    pos_r = f_axes[right].get_position()
    mid_x = (pos_l.x0 + pos_r.x1) / 2.0
    top_y = pos_l.y1 + 0.04
    fig.text(mid_x, top_y, phase_name,
             fontsize=FONT_PHASE, fontweight='bold', ha='center', va='bottom',
             color=PHASE_COLORS[phase_name])

# Panel labels
ax_d_pos = ax_d.get_position()
fig.text(ax_d_pos.x0 - 0.02, ax_d_pos.y1 + 0.02, 'D',
         fontsize=FONT_LABEL, fontweight='bold', va='bottom', ha='right')
ax_e_pos = ax_e.get_position()
fig.text(ax_e_pos.x0 - 0.01, ax_e_pos.y1 + 0.02, 'E',
         fontsize=FONT_LABEL, fontweight='bold', va='bottom', ha='right')
pos_f0 = f_axes[0].get_position()
fig.text(pos_f0.x0 - 0.02, pos_f0.y1 + 0.02, 'F',
         fontsize=FONT_LABEL, fontweight='bold', va='bottom', ha='right')

# Population legend
legend_elements = [
    Patch(facecolor=POP_COLORS[p], edgecolor='#333333', label=POP_LABELS[p])
    for p in POP_ORDER
]
fig.legend(handles=legend_elements, loc='lower center',
           ncol=3, fontsize=FONT_LEGEND + 2, framealpha=0.9,
           bbox_to_anchor=(0.5, -0.01))

save_fig(fig, "Panel_6DEF_HPV_overview")


# #############################################################################
#
#   COMPOSITE FIGURE
#
# #############################################################################
banner("COMPOSITE: Assembling Figure 6")

panel_specs = [
    ("Panel_6A_UMAP_HPV16_viral_load.png", 'A'),
    ("Panel_6B_baseline_violins.png", 'B'),
    ("Panel_6C_host_marker_dotplot.png", 'C'),
    ("Panel_6DEF_HPV_overview.png", 'DEF'),
]

panel_imgs = {}
for fname, label in panel_specs:
    img_path = os.path.join(OUTPUT_DIR, fname)
    if os.path.exists(img_path):
        panel_imgs[label] = plt.imread(img_path)
        log(f"  Loaded {fname}: {panel_imgs[label].shape}")
    else:
        log(f"  WARNING: {fname} not found")

if len(panel_imgs) == 4:
    h_A = panel_imgs['A'].shape[0] / panel_imgs['A'].shape[1]
    h_B = panel_imgs['B'].shape[0] / panel_imgs['B'].shape[1]
    h_C = panel_imgs['C'].shape[0] / panel_imgs['C'].shape[1]
    h_DEF = panel_imgs['DEF'].shape[0] / panel_imgs['DEF'].shape[1]

    fig_w = 36
    row1_h = max(h_A, h_B) * fig_w * 0.45
    row2_h = h_C * fig_w
    row3_h = h_DEF * fig_w
    fig_h = min(row1_h + row2_h + row3_h, 65)

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs_comp = gridspec.GridSpec(3, 2,
                                height_ratios=[row1_h, row2_h, row3_h],
                                hspace=0.06, wspace=0.04)

    ax_a = fig.add_subplot(gs_comp[0, 0])
    ax_a.imshow(panel_imgs['A'])
    ax_a.axis('off')
    ax_a.text(-0.02, 1.02, 'A', transform=ax_a.transAxes,
              fontsize=FONT_LABEL + 4, fontweight='bold', va='top', ha='right')

    ax_b = fig.add_subplot(gs_comp[0, 1])
    ax_b.imshow(panel_imgs['B'])
    ax_b.axis('off')
    ax_b.text(-0.02, 1.02, 'B', transform=ax_b.transAxes,
              fontsize=FONT_LABEL + 4, fontweight='bold', va='top', ha='right')

    ax_c = fig.add_subplot(gs_comp[1, :])
    ax_c.imshow(panel_imgs['C'])
    ax_c.axis('off')
    ax_c.text(-0.01, 1.01, 'C', transform=ax_c.transAxes,
              fontsize=FONT_LABEL + 4, fontweight='bold', va='top', ha='right')

    ax_def = fig.add_subplot(gs_comp[2, :])
    ax_def.imshow(panel_imgs['DEF'])
    ax_def.axis('off')
    # D/E/F labels already baked into the panel image

    save_fig(fig, "Figure6_composite")
else:
    log("  ERROR: Not all panel images available, skipping composite")


# =============================================================================
# SAVE REPORT
# =============================================================================
report_path = os.path.join(OUTPUT_DIR, "figure6_lifecycle_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"\n  Report: {report_path}")
log(f"  Output: {OUTPUT_DIR}")
banner("FIGURE 6 GENERATION COMPLETE")
