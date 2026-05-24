#!/usr/bin/env python3
"""
Generate_Figure7_Panels.py
============================
Figure 7: Neoantigen Landscape of HPV16 Lifecycle Populations

REDESIGNED LAYOUT:
  Row 1: Panel A (burden + fusions) | Panel B (Venn diagram) | Panel C (ANXA1 peptides by region)
  Row 2: Panel D (ANXA1 gene track, full width, with hotspot region highlights)

Key design changes from v1:
  - Panel A: Added RNA fusion rates alongside neoantigens for SBS2/CNV
  - Panel B: Venn diagram replacing tier infographic blocks
  - Panel C: ANXA1 peptides grouped by two hotspot regions (pos 289, pos 321-327)
  - Panel D: Full-width gene track, no APOBEC coloring (all non-APOBEC confirmed),
             colored by binding strength, hotspot regions highlighted

Run in NETWORK conda env.
Save as PDF and PNG at 300 DPI.

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Circle
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import re
import os
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
MHC_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/03_mhc_binding")
FUSION_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/04_fusion_analysis")
SUMMARY_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/05_summary")
TIER_DIR = os.path.join(SUMMARY_DIR, "THERAPEUTIC_TIERS")
ANNOTATION_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/02_snpeff_annotation")

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/FIGURE_7_PANELS")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Font sizes (28-34 range)
FONT_SIZE = 30
FONT_SIZE_SMALL = 26
FONT_SIZE_TICK = 24
FONT_SIZE_LABEL = 22
FONT_SIZE_ANNOT = 20

# Colors (hex only) - three-group scheme matching other figures
COL_SBS2 = '#ed6a5a'    # coral
COL_CNV = '#F6D155'      # mustard yellow
COL_NORMAL = '#4682b4'   # steelblue
COL_SHARED = '#9B59B6'
COL_STRONG = '#C0392B'
COL_BINDER = '#F39C12'
COL_WEAK = '#BDC3C7'
COL_WT = '#95A5A6'
COL_REGION1 = '#E74C3C'   # pos 289 region
COL_REGION2 = '#2980B9'   # pos 321-327 region
COL_OTHER_MUT = '#7F8C8D'

# ANXA1 domain boundaries (UniProt P04083)
ANXA1_LENGTH = 346
ANXA1_DOMAINS = [
    (1, 41, 'N-term', '#F5B7B1'),
    (42, 109, 'Repeat 1', '#AED6F1'),
    (110, 122, 'Linker 1-2', '#D5DBDB'),
    (123, 196, 'Repeat 2', '#A9DFBF'),
    (197, 220, 'Linker 2-3', '#D5DBDB'),
    (221, 287, 'Repeat 3', '#F9E79F'),
    (288, 290, 'Linker 3-4', '#D5DBDB'),
    (291, 346, 'Repeat 4', '#D7BDE2'),
]

plt.rcParams.update({
    'font.size': FONT_SIZE_TICK,
    'axes.titlesize': FONT_SIZE,
    'axes.labelsize': FONT_SIZE_SMALL,
    'xtick.labelsize': FONT_SIZE_TICK,
    'ytick.labelsize': FONT_SIZE_TICK,
    'legend.fontsize': FONT_SIZE_LABEL,
    'font.family': 'sans-serif',
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# =============================================================================
# LOAD DATA
# =============================================================================
print("Loading data...", flush=True)

N_CELLS = 546

# Neoantigens
sbs2_neo = pd.read_csv(os.path.join(MHC_DIR, "SBS2_HIGH_neoantigens.tsv"), sep='\t')
cnv_neo = pd.read_csv(os.path.join(MHC_DIR, "CNV_HIGH_neoantigens.tsv"), sep='\t')

sbs2_total_neo = len(sbs2_neo)
cnv_total_neo = len(cnv_neo)

# Fusions
junction_summary = pd.read_csv(os.path.join(FUSION_DIR, "per_group_junction_summary.tsv"), sep='\t')

sbs2_fusions = junction_summary[junction_summary['group'] == 'SBS2_HIGH']['total_junctions'].iloc[0]
cnv_fusions = junction_summary[junction_summary['group'] == 'CNV_HIGH']['total_junctions'].iloc[0]

# Venn diagram gene sets
sbs2_genes = set(sbs2_neo['gene'].dropna().unique())
cnv_genes = set(cnv_neo['gene'].dropna().unique())
shared_genes = sbs2_genes & cnv_genes
sbs2_only_genes = sbs2_genes - cnv_genes
cnv_only_genes = cnv_genes - sbs2_genes

# Tier data
tier_df = pd.read_csv(os.path.join(TIER_DIR, "all_genes_tiered.tsv"), sep='\t')
n_1a = len(tier_df[tier_df['tier'] == '1A_hot_shared_escaped'])
n_1b = len(tier_df[tier_df['tier'] == '1B_hot_sbs2_specific'])
n_2 = len(tier_df[tier_df['tier'] == '2_cold_cnv_specific'])
n_3 = len(tier_df[tier_df['tier'] == '3_broad_coverage'])

# ANXA1 peptide data
sbs2_all_pep = pd.read_csv(os.path.join(MHC_DIR, "SBS2_HIGH_all_peptide_results.tsv"), sep='\t')
anxa1_all_pep = sbs2_all_pep[sbs2_all_pep['gene'] == 'ANXA1'].copy()
anxa1_binders = anxa1_all_pep[anxa1_all_pep['is_binder'] == True].copy()

# ANXA1 somatic variants
sbs2_spa = pd.read_csv(os.path.join(ANNOTATION_DIR, "SBS2_HIGH.somatic_protein_altering.tsv"), sep='\t')
anxa1_variants = sbs2_spa[sbs2_spa['gene'] == 'ANXA1'].copy()

print(f"  SBS2 neoantigens: {sbs2_total_neo}, CNV neoantigens: {cnv_total_neo}", flush=True)
print(f"  SBS2 fusions: {sbs2_fusions}, CNV fusions: {cnv_fusions}", flush=True)
print(f"  Venn: SBS2-only={len(sbs2_only_genes)}, shared={len(shared_genes)}, CNV-only={len(cnv_only_genes)}", flush=True)
print(f"  ANXA1 binder peptides: {len(anxa1_binders)}", flush=True)

# =============================================================================
# PARSE ANXA1 MUTATIONS
# =============================================================================
print("Parsing ANXA1 mutations...", flush=True)

anxa1_muts = []
for _, var in anxa1_variants.iterrows():
    hgvs_p = str(var.get('hgvs_p', ''))
    match = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', hgvs_p)
    if not match:
        continue
    pos = int(match.group(2))

    # Cross-reference with binding data
    if 'mut_pos_protein' in anxa1_all_pep.columns:
        pos_pep = anxa1_all_pep[anxa1_all_pep['mut_pos_protein'] == pos]
    else:
        pos_pep = pd.DataFrame()

    best_ic50 = pos_pep['mut_ic50'].min() if len(pos_pep) > 0 else 999
    n_strong = (pos_pep['mut_ic50'] < 50).sum() if len(pos_pep) > 0 else 0

    # Assign region
    if pos == 289:
        region = 'region1'
    elif 321 <= pos <= 327:
        region = 'region2'
    else:
        region = 'other'

    anxa1_muts.append({
        'position': pos,
        'hgvs_p': hgvs_p,
        'best_ic50': best_ic50,
        'n_strong': n_strong,
        'is_strong': best_ic50 < 50,
        'is_binder': best_ic50 < 500,
        'region': region,
    })

# Deduplicate to unique positions (keep best IC50)
pos_to_mut = {}
for m in anxa1_muts:
    p = m['position']
    if p not in pos_to_mut or m['best_ic50'] < pos_to_mut[p]['best_ic50']:
        pos_to_mut[p] = m
anxa1_muts_unique = sorted(pos_to_mut.values(), key=lambda x: x['position'])

print(f"  Unique mutation positions: {len(anxa1_muts_unique)}", flush=True)
for m in anxa1_muts_unique:
    print(f"    pos={m['position']}  IC50={m['best_ic50']:.1f}  region={m['region']}  strong={m['is_strong']}", flush=True)


# =============================================================================
# PANEL A: NEOANTIGEN + FUSION BURDEN (individual panel)
# =============================================================================
print("Panel A: Neoantigen + Fusion burden...", flush=True)

fig_a, ax_a = plt.subplots(figsize=(14, 10))

categories = ['Predicted\nneoantigens', 'RNA\nfusions']
sbs2_vals = [sbs2_total_neo, sbs2_fusions]
cnv_vals = [cnv_total_neo, cnv_fusions]

x = np.arange(len(categories))
width = 0.35

bars1 = ax_a.bar(x - width/2, sbs2_vals, width, color=COL_SBS2, label='SBS2-HIGH',
                 edgecolor='white', linewidth=1.5)
bars2 = ax_a.bar(x + width/2, cnv_vals, width, color=COL_CNV, label='CNV-HIGH',
                 edgecolor='white', linewidth=1.5)

for bar in bars1:
    h = bar.get_height()
    ax_a.text(bar.get_x() + bar.get_width()/2., h + max(sbs2_vals) * 0.01,
              f'{int(h):,}', ha='center', va='bottom', fontsize=FONT_SIZE_LABEL, fontweight='bold')
for bar in bars2:
    h = bar.get_height()
    ax_a.text(bar.get_x() + bar.get_width()/2., h + max(sbs2_vals) * 0.01,
              f'{int(h):,}', ha='center', va='bottom', fontsize=FONT_SIZE_LABEL, fontweight='bold')

# Fold change labels
for i in range(len(categories)):
    if cnv_vals[i] > 0:
        fold = sbs2_vals[i] / cnv_vals[i]
        mid_x = x[i]
        max_y = max(sbs2_vals[i], cnv_vals[i])
        ax_a.text(mid_x, max_y + max(max(sbs2_vals), max(cnv_vals)) * 0.08,
                  f'{fold:.1f}x', ha='center', va='bottom',
                  fontsize=FONT_SIZE_ANNOT, fontstyle='italic', color='#555555')

ax_a.set_ylabel('Total observed count', fontsize=FONT_SIZE_SMALL)
ax_a.set_xticks(x)
ax_a.set_xticklabels(categories, fontsize=FONT_SIZE_TICK)
ax_a.legend(fontsize=FONT_SIZE_LABEL, frameon=False, loc='upper right')
ax_a.set_title('Neoantigen and fusion burden', fontsize=FONT_SIZE, fontweight='bold', pad=20)
ax_a.set_ylim(0, max(max(sbs2_vals), max(cnv_vals)) * 1.3)

fig_a.tight_layout()
fig_a.savefig(os.path.join(OUTPUT_DIR, "Figure7_PanelA_Burden.pdf"), dpi=300, bbox_inches='tight')
fig_a.savefig(os.path.join(OUTPUT_DIR, "Figure7_PanelA_Burden.png"), dpi=300, bbox_inches='tight')
plt.close(fig_a)
print("  Panel A saved.", flush=True)


# =============================================================================
# PANEL B: VENN DIAGRAM (individual panel)
# =============================================================================
print("Panel B: Venn diagram...", flush=True)

fig_b, ax_b = plt.subplots(figsize=(14, 10))
ax_b.set_xlim(-6, 6)
ax_b.set_ylim(-5, 6)
ax_b.set_aspect('equal')
ax_b.axis('off')

# Draw two overlapping circles
r = 3.2
offset = 1.8  # controls overlap amount

circle_sbs2 = Circle((-offset, 0), r, fill=True, facecolor=COL_SBS2, alpha=0.15,
                      edgecolor=COL_SBS2, linewidth=3)
circle_cnv = Circle((offset, 0), r, fill=True, facecolor=COL_CNV, alpha=0.15,
                     edgecolor=COL_CNV, linewidth=3)
ax_b.add_patch(circle_sbs2)
ax_b.add_patch(circle_cnv)

# Region counts
ax_b.text(-offset - 1.2, 0.3, str(len(sbs2_only_genes)), fontsize=FONT_SIZE + 6,
          ha='center', va='center', fontweight='bold', color=COL_SBS2)
ax_b.text(-offset - 1.2, -0.6, 'SBS2-only', fontsize=FONT_SIZE_LABEL,
          ha='center', va='center', color=COL_SBS2)

ax_b.text(offset + 1.2, 0.3, str(len(cnv_only_genes)), fontsize=FONT_SIZE + 6,
          ha='center', va='center', fontweight='bold', color=COL_CNV)
ax_b.text(offset + 1.2, -0.6, 'CNV-only', fontsize=FONT_SIZE_LABEL,
          ha='center', va='center', color=COL_CNV)

ax_b.text(0, 0.3, str(len(shared_genes)), fontsize=FONT_SIZE + 6,
          ha='center', va='center', fontweight='bold', color=COL_SHARED)
ax_b.text(0, -0.6, 'shared', fontsize=FONT_SIZE_LABEL,
          ha='center', va='center', color=COL_SHARED)

# Circle labels
ax_b.text(-offset, r + 0.6, f'SBS2-HIGH\n({len(sbs2_genes)} genes)',
          fontsize=FONT_SIZE_SMALL, ha='center', va='bottom', fontweight='bold', color=COL_SBS2)
ax_b.text(offset, r + 0.6, f'CNV-HIGH\n({len(cnv_genes)} genes)',
          fontsize=FONT_SIZE_SMALL, ha='center', va='bottom', fontweight='bold', color=COL_CNV)

# Tier breakout annotations with arrows
# SBS2-only = Tier 1B (maintenance-phase targets)
ax_b.annotate(f'Tier 1B: {n_1b} maintenance\ntargets (ANXA1, S100A8)',
              xy=(-offset - 1.2, -1.2), xytext=(-offset - 1.5, -3.8),
              fontsize=FONT_SIZE_ANNOT, ha='center', va='top', color=COL_SBS2,
              arrowprops=dict(arrowstyle='->', color=COL_SBS2, lw=2))

# CNV-only = Tier 2 (productive-phase targets)
ax_b.annotate(f'Tier 2: {n_2} productive\ntargets (COX4I1, ENO1)',
              xy=(offset + 1.2, -1.2), xytext=(offset + 1.5, -3.8),
              fontsize=FONT_SIZE_ANNOT, ha='center', va='top', color=COL_CNV,
              arrowprops=dict(arrowstyle='->', color=COL_CNV, lw=2))

# Shared = Tier 1A (escaped) + Tier 3 (broad)
ax_b.annotate(f'Tier 1A: {n_1a} immune-escaped\nTier 3: {n_3} broad coverage',
              xy=(0, -1.2), xytext=(0, -3.8),
              fontsize=FONT_SIZE_ANNOT, ha='center', va='top', color=COL_SHARED,
              arrowprops=dict(arrowstyle='->', color=COL_SHARED, lw=2))

ax_b.set_title('Neoantigen gene overlap', fontsize=FONT_SIZE, fontweight='bold', pad=20)

fig_b.tight_layout()
fig_b.savefig(os.path.join(OUTPUT_DIR, "Figure7_PanelB_Venn.pdf"), dpi=300, bbox_inches='tight')
fig_b.savefig(os.path.join(OUTPUT_DIR, "Figure7_PanelB_Venn.png"), dpi=300, bbox_inches='tight')
plt.close(fig_b)
print("  Panel B saved.", flush=True)


# =============================================================================
# PANEL C: ANXA1 PEPTIDES BY REGION (individual panel)
# =============================================================================
print("Panel C: ANXA1 peptides by region...", flush=True)

if len(anxa1_binders) > 0:
    # Add IC50 drop column (positive = mutation improves binding)
    anxa1_binders['ic50_drop'] = anxa1_binders['wt_ic50'] - anxa1_binders['mut_ic50']

    # Split by region
    if 'mut_pos_protein' in anxa1_binders.columns:
        region1 = anxa1_binders[anxa1_binders['mut_pos_protein'] == 289].copy()
        region2 = anxa1_binders[anxa1_binders['mut_pos_protein'].isin([321, 322, 326, 327])].copy()
    else:
        region1 = pd.DataFrame()
        region2 = pd.DataFrame()

    # Deduplicate by peptide, keep largest IC50 drop
    if 'mut_peptide' in region1.columns and len(region1) > 0:
        region1 = region1.sort_values('ic50_drop', ascending=False).drop_duplicates(subset=['mut_peptide'], keep='first')
    if 'mut_peptide' in region2.columns and len(region2) > 0:
        region2 = region2.sort_values('ic50_drop', ascending=False).drop_duplicates(subset=['mut_peptide'], keep='first')

    # Take top peptides by IC50 drop, then reverse for display (biggest drop at top)
    n_r1 = min(7, len(region1))
    n_r2 = min(7, len(region2))
    top_r1 = region1.head(n_r1).sort_values('ic50_drop', ascending=True)   # smallest at bottom
    top_r2 = region2.head(n_r2).sort_values('ic50_drop', ascending=True)

    # Combine with a gap
    fig_c, ax_c = plt.subplots(figsize=(16, 12))

    # Build combined data with spacing
    all_peptides = []
    all_mut_ic50 = []
    all_wt_ic50 = []
    all_colors = []
    all_labels = []        # suffix text (allele, position)
    all_pep_seqs = []      # peptide sequence only
    all_mut_pos_pep = []   # mutation position within peptide (0-indexed)
    all_is_separator = []

    # Region 2 first (bottom)
    for _, row in top_r2.iterrows():
        pep = row.get('mut_peptide', '')
        allele = str(row.get('best_allele', '')).replace('HLA-', '')
        hgvs = str(row.get('hgvs_p', ''))
        match_h = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', hgvs)
        pos_label = match_h.group(2) if match_h else ''
        mut_in_pep = int(row.get('mut_position_in_peptide', -1)) if 'mut_position_in_peptide' in row.index else -1
        all_pep_seqs.append(pep)
        all_labels.append(f"  {allele}  (pos {pos_label})")
        all_mut_pos_pep.append(mut_in_pep)
        all_mut_ic50.append(row['mut_ic50'])
        all_wt_ic50.append(row.get('wt_ic50', 999))
        all_colors.append(COL_REGION2)
        all_is_separator.append(False)

    # Separator
    all_pep_seqs.append('')
    all_labels.append('')
    all_mut_pos_pep.append(-1)
    all_mut_ic50.append(0)
    all_wt_ic50.append(0)
    all_colors.append('none')
    all_is_separator.append(True)

    # Region 1 (top)
    for _, row in top_r1.iterrows():
        pep = row.get('mut_peptide', '')
        allele = str(row.get('best_allele', '')).replace('HLA-', '')
        hgvs = str(row.get('hgvs_p', ''))
        match_h = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', hgvs)
        pos_label = match_h.group(2) if match_h else ''
        mut_in_pep = int(row.get('mut_position_in_peptide', -1)) if 'mut_position_in_peptide' in row.index else -1
        all_pep_seqs.append(pep)
        all_labels.append(f"  {allele}  (pos {pos_label})")
        all_mut_pos_pep.append(mut_in_pep)
        all_mut_ic50.append(row['mut_ic50'])
        all_wt_ic50.append(row.get('wt_ic50', 999))
        all_colors.append(COL_REGION1)
        all_is_separator.append(False)

    y_pos = np.arange(len(all_labels))

    # WT bars (background)
    for i in range(len(all_labels)):
        if all_is_separator[i]:
            continue
        ax_c.barh(y_pos[i], min(all_wt_ic50[i], 2000),
                  height=0.6, color=COL_WT, alpha=0.4, edgecolor='white', linewidth=0.5)

    # Mutant bars
    for i in range(len(all_labels)):
        if all_is_separator[i]:
            continue
        color = all_colors[i]
        alpha = 0.9
        ax_c.barh(y_pos[i], all_mut_ic50[i],
                  height=0.6, color=color, alpha=alpha, edgecolor='white', linewidth=0.5)

    # Threshold lines
    ax_c.axvline(x=500, color='#333333', linestyle='--', linewidth=2, alpha=0.4, zorder=0)
    ax_c.axvline(x=50, color=COL_STRONG, linestyle=':', linewidth=2, alpha=0.4, zorder=0)

    # Region labels
    sep_idx = all_is_separator.index(True)
    r2_mid = sep_idx / 2
    r1_mid = sep_idx + 1 + (len(all_labels) - sep_idx - 1) / 2

    ax_c.text(max(all_wt_ic50) * 0.85, r1_mid, 'Pos 289\nregion',
              fontsize=FONT_SIZE_LABEL, ha='center', va='center',
              fontweight='bold', color=COL_REGION1, alpha=0.7,
              bbox=dict(boxstyle='round,pad=0.3', facecolor=COL_REGION1, alpha=0.1, edgecolor='none'))
    ax_c.text(max(all_wt_ic50) * 0.85, r2_mid, 'Pos 321-327\nregion',
              fontsize=FONT_SIZE_LABEL, ha='center', va='center',
              fontweight='bold', color=COL_REGION2, alpha=0.7,
              bbox=dict(boxstyle='round,pad=0.3', facecolor=COL_REGION2, alpha=0.1, edgecolor='none'))

    # Custom y-axis labels with highlighted mutated AA
    # Use monospace text drawn manually so we can color the mutated position
    ax_c.set_yticks(y_pos)
    ax_c.set_yticklabels(['' for _ in all_labels])  # blank out default labels

    label_fontsize = FONT_SIZE_ANNOT - 2
    highlight_color = '#D4380D'  # bold orange-red for mutated AA

    for i in range(len(all_labels)):
        if all_is_separator[i]:
            continue
        pep = all_pep_seqs[i]
        suffix = all_labels[i]
        mut_idx = all_mut_pos_pep[i]

        if mut_idx >= 0 and mut_idx < len(pep):
            # Three segments: before, mutated, after + suffix
            before = pep[:mut_idx]
            mutated = pep[mut_idx]
            after = pep[mut_idx + 1:]

            # Draw each segment at calculated x positions
            # Use transData for y, figure out x positioning relative to axis
            # We draw at a fixed x in axes coordinates using blended transform
            from matplotlib.transforms import blended_transform_factory
            trans = blended_transform_factory(ax_c.transAxes, ax_c.transData)

            # Approximate character width in axes fraction for monospace at this fontsize
            # Start position (right-aligned to axis edge)
            max_pep_len = max(len(s) for s in all_pep_seqs if s)
            max_suffix_len = max(len(s) for s in all_labels if s)
            total_chars = max_pep_len + max_suffix_len
            char_width = 0.018  # approximate monospace char width in axes fraction

            # Starting x position (from right edge of label area)
            x_start = -char_width * total_chars

            # Before mutation
            ax_c.text(x_start, y_pos[i], before,
                      fontsize=label_fontsize, fontfamily='monospace',
                      ha='left', va='center', color='#333333',
                      transform=trans, clip_on=False)
            # Mutated AA (highlighted)
            ax_c.text(x_start + char_width * len(before), y_pos[i], mutated,
                      fontsize=label_fontsize, fontfamily='monospace',
                      ha='left', va='center', color=highlight_color, fontweight='bold',
                      transform=trans, clip_on=False,
                      bbox=dict(boxstyle='round,pad=0.1', facecolor=highlight_color,
                                alpha=0.12, edgecolor='none'))
            # After mutation
            ax_c.text(x_start + char_width * (len(before) + 1), y_pos[i], after,
                      fontsize=label_fontsize, fontfamily='monospace',
                      ha='left', va='center', color='#333333',
                      transform=trans, clip_on=False)
            # Suffix (allele + position)
            ax_c.text(x_start + char_width * len(pep), y_pos[i], suffix,
                      fontsize=label_fontsize, fontfamily='monospace',
                      ha='left', va='center', color='#888888',
                      transform=trans, clip_on=False)
        else:
            # No mutation position info, draw plain
            from matplotlib.transforms import blended_transform_factory
            trans = blended_transform_factory(ax_c.transAxes, ax_c.transData)
            max_pep_len = max(len(s) for s in all_pep_seqs if s)
            max_suffix_len = max(len(s) for s in all_labels if s)
            total_chars = max_pep_len + max_suffix_len
            char_width = 0.018
            x_start = -char_width * total_chars

            ax_c.text(x_start, y_pos[i], pep + suffix,
                      fontsize=label_fontsize, fontfamily='monospace',
                      ha='left', va='center', color='#333333',
                      transform=trans, clip_on=False)

    ax_c.set_xlabel('Predicted IC50 (nM)', fontsize=FONT_SIZE_SMALL)
    ax_c.set_title('ANXA1 neoantigen peptides by hotspot region',
                   fontsize=FONT_SIZE, fontweight='bold', pad=20)
    ax_c.set_xlim(0, min(2000, max(all_wt_ic50) * 1.1))

    legend_c = [
        mpatches.Patch(facecolor=COL_REGION1, alpha=0.9, label='Pos 289 (Phe289Leu/Ile)'),
        mpatches.Patch(facecolor=COL_REGION2, alpha=0.9, label='Pos 321-327 cluster'),
        mpatches.Patch(facecolor=COL_WT, alpha=0.4, label='Wild-type IC50'),
    ]
    ax_c.legend(handles=legend_c, loc='lower right', fontsize=FONT_SIZE_LABEL,
                frameon=True, fancybox=True, framealpha=0.9)

    fig_c.tight_layout()
    fig_c.savefig(os.path.join(OUTPUT_DIR, "Figure7_PanelC_ANXA1_Peptides.pdf"), dpi=300, bbox_inches='tight')
    fig_c.savefig(os.path.join(OUTPUT_DIR, "Figure7_PanelC_ANXA1_Peptides.png"), dpi=300, bbox_inches='tight')
    plt.close(fig_c)
    print(f"  Panel C saved (R1: {n_r1}, R2: {n_r2} peptides).", flush=True)
else:
    print("  WARNING: No ANXA1 binders for Panel C", flush=True)


# =============================================================================
# PANEL D: ANXA1 GENE TRACK (individual panel, wide)
# =============================================================================
print("Panel D: ANXA1 gene track...", flush=True)

fig_d, ax_d = plt.subplots(figsize=(24, 6))

backbone_y = 0
backbone_height = 0.6
max_lollipop_height = 3.5

# Draw domain backbone
for start, end, name, color in ANXA1_DOMAINS:
    ax_d.add_patch(plt.Rectangle((start, backbone_y - backbone_height/2),
                                  end - start, backbone_height,
                                  facecolor=color, edgecolor='#333333',
                                  linewidth=1.5, zorder=2))
    mid = (start + end) / 2
    if (end - start) > 15:  # only label domains wide enough
        ax_d.text(mid, backbone_y, name, ha='center', va='center',
                  fontsize=FONT_SIZE_ANNOT, fontweight='bold', zorder=3)

# Highlight hotspot regions with shaded boxes
region1_box = plt.Rectangle((285, -backbone_height - 0.3), 10, backbone_height + max_lollipop_height + 1.0,
                              facecolor=COL_REGION1, alpha=0.08, edgecolor=COL_REGION1,
                              linewidth=2, linestyle='--', zorder=1)
ax_d.add_patch(region1_box)
ax_d.text(290, max_lollipop_height + 1.0, 'Pos 289\nregion',
          ha='center', va='bottom', fontsize=FONT_SIZE_ANNOT,
          fontweight='bold', color=COL_REGION1)

region2_box = plt.Rectangle((318, -backbone_height - 0.3), 14, backbone_height + max_lollipop_height + 1.0,
                              facecolor=COL_REGION2, alpha=0.08, edgecolor=COL_REGION2,
                              linewidth=2, linestyle='--', zorder=1)
ax_d.add_patch(region2_box)
ax_d.text(325, max_lollipop_height + 1.0, 'Pos 321-327\nregion',
          ha='center', va='bottom', fontsize=FONT_SIZE_ANNOT,
          fontweight='bold', color=COL_REGION2)

# Draw lollipop mutations
for m in anxa1_muts_unique:
    pos = m['position']

    # Height by binding strength (inverted IC50)
    if m['best_ic50'] < 500:
        height = max_lollipop_height * (500 - m['best_ic50']) / 500
        height = max(height, 1.0)
    else:
        height = 0.7

    # Color by region
    if m['region'] == 'region1':
        color = COL_REGION1
        edge = '#8B0000'
    elif m['region'] == 'region2':
        color = COL_REGION2
        edge = '#1A4F7A'
    else:
        color = COL_OTHER_MUT
        edge = '#4A4A4A'

    marker_size = 14 if m['is_strong'] else 10

    # Stem
    ax_d.plot([pos, pos], [backbone_height/2, backbone_height/2 + height],
              color=color, linewidth=2.5, alpha=0.7, zorder=4)
    # Head
    ax_d.scatter(pos, backbone_height/2 + height, s=marker_size**2,
                 color=color, edgecolors=edge, linewidth=2,
                 zorder=5, alpha=0.9)

    # Label strong binders with IC50
    if m['is_strong']:
        ax_d.text(pos, backbone_height/2 + height + 0.25,
                  f'{m["best_ic50"]:.0f}nM',
                  ha='center', va='bottom', fontsize=FONT_SIZE_ANNOT - 2,
                  fontweight='bold', color=color)

# Mutation labels below backbone
for m in anxa1_muts_unique:
    pos = m['position']
    # Short HGVS label
    match_h = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', m['hgvs_p'])
    if match_h:
        label = f"{match_h.group(1)}{match_h.group(2)}{match_h.group(3)}"
    else:
        label = m['hgvs_p']
    ax_d.text(pos, -backbone_height/2 - 0.3, label,
              ha='center', va='top', fontsize=FONT_SIZE_ANNOT - 4,
              rotation=45, color='#333333')

ax_d.set_xlim(-10, ANXA1_LENGTH + 15)
ax_d.set_ylim(-2.5, max_lollipop_height + 2.5)
ax_d.set_xlabel('ANXA1 protein position (amino acids)', fontsize=FONT_SIZE_SMALL)
ax_d.set_yticks([])
ax_d.spines['left'].set_visible(False)
ax_d.set_title(f'ANXA1 somatic mutations in SBS2-HIGH cells ({len(anxa1_muts_unique)} positions)',
               fontsize=FONT_SIZE, fontweight='bold', pad=20)

# Legend
legend_d = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=COL_REGION1,
               markersize=14, markeredgecolor='#8B0000', markeredgewidth=2,
               label='Pos 289 region'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=COL_REGION2,
               markersize=14, markeredgecolor='#1A4F7A', markeredgewidth=2,
               label='Pos 321-327 region'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=COL_OTHER_MUT,
               markersize=10, markeredgecolor='#4A4A4A', markeredgewidth=2,
               label='Other positions'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#555555',
               markersize=14, label='Strong binder (IC50 < 50nM)'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#555555',
               markersize=9, label='Binder (IC50 < 500nM)'),
]
ax_d.legend(handles=legend_d, loc='upper left', fontsize=FONT_SIZE_ANNOT,
            frameon=True, fancybox=True, framealpha=0.9, ncol=2)

fig_d.tight_layout()
fig_d.savefig(os.path.join(OUTPUT_DIR, "Figure7_PanelD_ANXA1_GeneTrack.pdf"), dpi=300, bbox_inches='tight')
fig_d.savefig(os.path.join(OUTPUT_DIR, "Figure7_PanelD_ANXA1_GeneTrack.png"), dpi=300, bbox_inches='tight')
plt.close(fig_d)
print("  Panel D saved.", flush=True)


# =============================================================================
# COMBINED FIGURE
# =============================================================================
print("Combined figure...", flush=True)

fig = plt.figure(figsize=(48, 24))
gs = gridspec.GridSpec(2, 3, height_ratios=[1.2, 0.8], hspace=0.35, wspace=0.3)

# --- Panel A (top left) ---
ax1 = fig.add_subplot(gs[0, 0])

categories_short = ['Predicted\nneoantigens', 'RNA\nfusions']
x = np.arange(len(categories_short))
width = 0.35

b1 = ax1.bar(x - width/2, sbs2_vals, width, color=COL_SBS2, label='SBS2-HIGH',
             edgecolor='white', linewidth=1.5)
b2 = ax1.bar(x + width/2, cnv_vals, width, color=COL_CNV, label='CNV-HIGH',
             edgecolor='white', linewidth=1.5)

for bar in b1:
    h = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., h + max(sbs2_vals) * 0.01,
             f'{int(h):,}', ha='center', va='bottom', fontsize=FONT_SIZE_LABEL, fontweight='bold')
for bar in b2:
    h = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., h + max(sbs2_vals) * 0.01,
             f'{int(h):,}', ha='center', va='bottom', fontsize=FONT_SIZE_LABEL, fontweight='bold')

for i in range(len(categories_short)):
    if cnv_vals[i] > 0:
        fold = sbs2_vals[i] / cnv_vals[i]
        mid_x = x[i]
        max_y = max(sbs2_vals[i], cnv_vals[i])
        ax1.text(mid_x, max_y + max(max(sbs2_vals), max(cnv_vals)) * 0.08,
                 f'{fold:.1f}x', ha='center', va='bottom',
                 fontsize=FONT_SIZE_ANNOT, fontstyle='italic', color='#555555')

ax1.set_ylabel('Total observed count', fontsize=FONT_SIZE_SMALL)
ax1.set_xticks(x)
ax1.set_xticklabels(categories_short, fontsize=FONT_SIZE_TICK)
ax1.legend(fontsize=FONT_SIZE_LABEL, frameon=False)
ax1.set_title('A   Neoantigen and fusion burden', fontsize=FONT_SIZE, fontweight='bold', loc='left')
ax1.set_ylim(0, max(max(sbs2_vals), max(cnv_vals)) * 1.3)

# --- Panel B (top center): Venn ---
ax2 = fig.add_subplot(gs[0, 1])
ax2.set_xlim(-6, 6)
ax2.set_ylim(-5.5, 5.5)
ax2.set_aspect('equal')
ax2.axis('off')

r = 3.0
offset = 1.7

circle1 = Circle((-offset, 0), r, fill=True, facecolor=COL_SBS2, alpha=0.15,
                  edgecolor=COL_SBS2, linewidth=3)
circle2 = Circle((offset, 0), r, fill=True, facecolor=COL_CNV, alpha=0.15,
                  edgecolor=COL_CNV, linewidth=3)
ax2.add_patch(circle1)
ax2.add_patch(circle2)

ax2.text(-offset - 1.1, 0.3, str(len(sbs2_only_genes)), fontsize=FONT_SIZE + 4,
         ha='center', va='center', fontweight='bold', color=COL_SBS2)
ax2.text(-offset - 1.1, -0.5, 'SBS2-only', fontsize=FONT_SIZE_ANNOT,
         ha='center', va='center', color=COL_SBS2)

ax2.text(offset + 1.1, 0.3, str(len(cnv_only_genes)), fontsize=FONT_SIZE + 4,
         ha='center', va='center', fontweight='bold', color=COL_CNV)
ax2.text(offset + 1.1, -0.5, 'CNV-only', fontsize=FONT_SIZE_ANNOT,
         ha='center', va='center', color=COL_CNV)

ax2.text(0, 0.3, str(len(shared_genes)), fontsize=FONT_SIZE + 4,
         ha='center', va='center', fontweight='bold', color=COL_SHARED)
ax2.text(0, -0.5, 'shared', fontsize=FONT_SIZE_ANNOT,
         ha='center', va='center', color=COL_SHARED)

ax2.text(-offset, r + 0.4, f'SBS2-HIGH ({len(sbs2_genes)})',
         fontsize=FONT_SIZE_SMALL, ha='center', va='bottom', fontweight='bold', color=COL_SBS2)
ax2.text(offset, r + 0.4, f'CNV-HIGH ({len(cnv_genes)})',
         fontsize=FONT_SIZE_SMALL, ha='center', va='bottom', fontweight='bold', color=COL_CNV)

# Tier breakouts
ax2.annotate(f'Tier 1B: {n_1b}\n(ANXA1, S100A8)',
             xy=(-offset - 1.1, -1.0), xytext=(-offset - 1.3, -3.8),
             fontsize=FONT_SIZE_ANNOT, ha='center', va='top', color=COL_SBS2,
             arrowprops=dict(arrowstyle='->', color=COL_SBS2, lw=2))
ax2.annotate(f'Tier 2: {n_2}\n(COX4I1, ENO1)',
             xy=(offset + 1.1, -1.0), xytext=(offset + 1.3, -3.8),
             fontsize=FONT_SIZE_ANNOT, ha='center', va='top', color=COL_CNV,
             arrowprops=dict(arrowstyle='->', color=COL_CNV, lw=2))
ax2.annotate(f'1A: {n_1a} escaped\n3: {n_3} broad',
             xy=(0, -1.0), xytext=(0, -3.8),
             fontsize=FONT_SIZE_ANNOT, ha='center', va='top', color=COL_SHARED,
             arrowprops=dict(arrowstyle='->', color=COL_SHARED, lw=2))

ax2.set_title('B   Neoantigen gene overlap', fontsize=FONT_SIZE, fontweight='bold', loc='left')

# --- Panel C (top right): ANXA1 peptides ---
ax3 = fig.add_subplot(gs[0, 2])

if len(anxa1_binders) > 0:
    y_pos_c = np.arange(len(all_labels))

    for i in range(len(all_labels)):
        if all_is_separator[i]:
            continue
        ax3.barh(y_pos_c[i], min(all_wt_ic50[i], 2000),
                 height=0.6, color=COL_WT, alpha=0.4, edgecolor='white', linewidth=0.5)

    for i in range(len(all_labels)):
        if all_is_separator[i]:
            continue
        ax3.barh(y_pos_c[i], all_mut_ic50[i],
                 height=0.6, color=all_colors[i], alpha=0.9, edgecolor='white', linewidth=0.5)

    ax3.axvline(x=500, color='#333333', linestyle='--', linewidth=2, alpha=0.4, zorder=0)
    ax3.axvline(x=50, color=COL_STRONG, linestyle=':', linewidth=2, alpha=0.4, zorder=0)

    # Region labels
    ax3.text(min(2000, max(all_wt_ic50)) * 0.82, r1_mid, 'Pos 289',
             fontsize=FONT_SIZE_ANNOT, ha='center', va='center',
             fontweight='bold', color=COL_REGION1, alpha=0.7)
    ax3.text(min(2000, max(all_wt_ic50)) * 0.82, r2_mid, 'Pos 321-327',
             fontsize=FONT_SIZE_ANNOT, ha='center', va='center',
             fontweight='bold', color=COL_REGION2, alpha=0.7)

    ax3.set_yticks(y_pos_c)
    ax3.set_yticklabels(['' for _ in all_labels])  # blank out default labels

    # Draw peptide labels with highlighted mutated AA
    from matplotlib.transforms import blended_transform_factory
    trans3 = blended_transform_factory(ax3.transAxes, ax3.transData)
    label_fs = FONT_SIZE_ANNOT - 4
    hl_color = '#D4380D'

    for i in range(len(all_labels)):
        if all_is_separator[i]:
            continue
        pep = all_pep_seqs[i]
        suffix = all_labels[i]
        mut_idx = all_mut_pos_pep[i]

        max_pl = max(len(s) for s in all_pep_seqs if s)
        max_sl = max(len(s) for s in all_labels if s)
        cw = 0.016
        xs = -cw * (max_pl + max_sl)

        if mut_idx >= 0 and mut_idx < len(pep):
            before = pep[:mut_idx]
            mutated = pep[mut_idx]
            after = pep[mut_idx + 1:]

            ax3.text(xs, y_pos_c[i], before,
                     fontsize=label_fs, fontfamily='monospace',
                     ha='left', va='center', color='#333333',
                     transform=trans3, clip_on=False)
            ax3.text(xs + cw * len(before), y_pos_c[i], mutated,
                     fontsize=label_fs, fontfamily='monospace',
                     ha='left', va='center', color=hl_color, fontweight='bold',
                     transform=trans3, clip_on=False,
                     bbox=dict(boxstyle='round,pad=0.1', facecolor=hl_color,
                               alpha=0.12, edgecolor='none'))
            ax3.text(xs + cw * (len(before) + 1), y_pos_c[i], after,
                     fontsize=label_fs, fontfamily='monospace',
                     ha='left', va='center', color='#333333',
                     transform=trans3, clip_on=False)
            ax3.text(xs + cw * len(pep), y_pos_c[i], suffix,
                     fontsize=label_fs, fontfamily='monospace',
                     ha='left', va='center', color='#888888',
                     transform=trans3, clip_on=False)
        else:
            ax3.text(xs, y_pos_c[i], pep + suffix,
                     fontsize=label_fs, fontfamily='monospace',
                     ha='left', va='center', color='#333333',
                     transform=trans3, clip_on=False)

    ax3.set_xlabel('Predicted IC50 (nM)', fontsize=FONT_SIZE_SMALL)
    ax3.set_xlim(0, min(2000, max(all_wt_ic50) * 1.1))

    legend_c_combined = [
        mpatches.Patch(facecolor=COL_REGION1, alpha=0.9, label='Pos 289'),
        mpatches.Patch(facecolor=COL_REGION2, alpha=0.9, label='Pos 321-327'),
        mpatches.Patch(facecolor=COL_WT, alpha=0.4, label='Wild-type'),
    ]
    ax3.legend(handles=legend_c_combined, loc='lower right', fontsize=FONT_SIZE_ANNOT,
               frameon=True, fancybox=True, framealpha=0.9)

ax3.set_title('C   ANXA1 peptides by hotspot region', fontsize=FONT_SIZE, fontweight='bold', loc='left')

# --- Panel D (bottom, full width): Gene track ---
ax4 = fig.add_subplot(gs[1, :])

for start, end, name, color in ANXA1_DOMAINS:
    ax4.add_patch(plt.Rectangle((start, -backbone_height/2), end - start, backbone_height,
                                 facecolor=color, edgecolor='#333333', linewidth=1.5, zorder=2))
    mid = (start + end) / 2
    if (end - start) > 15:
        ax4.text(mid, 0, name, ha='center', va='center',
                 fontsize=FONT_SIZE_ANNOT, fontweight='bold', zorder=3)

# Hotspot region highlights
r1_box = plt.Rectangle((285, -backbone_height - 0.3), 10, backbone_height + max_lollipop_height + 0.8,
                         facecolor=COL_REGION1, alpha=0.08, edgecolor=COL_REGION1,
                         linewidth=2, linestyle='--', zorder=1)
ax4.add_patch(r1_box)
ax4.text(290, max_lollipop_height + 0.7, 'Pos 289',
         ha='center', va='bottom', fontsize=FONT_SIZE_ANNOT,
         fontweight='bold', color=COL_REGION1)

r2_box = plt.Rectangle((318, -backbone_height - 0.3), 14, backbone_height + max_lollipop_height + 0.8,
                         facecolor=COL_REGION2, alpha=0.08, edgecolor=COL_REGION2,
                         linewidth=2, linestyle='--', zorder=1)
ax4.add_patch(r2_box)
ax4.text(325, max_lollipop_height + 0.7, 'Pos 321-327',
         ha='center', va='bottom', fontsize=FONT_SIZE_ANNOT,
         fontweight='bold', color=COL_REGION2)

# Mutations
for m in anxa1_muts_unique:
    pos = m['position']
    if m['best_ic50'] < 500:
        height = max_lollipop_height * (500 - m['best_ic50']) / 500
        height = max(height, 1.0)
    else:
        height = 0.7

    if m['region'] == 'region1':
        color = COL_REGION1
        edge = '#8B0000'
    elif m['region'] == 'region2':
        color = COL_REGION2
        edge = '#1A4F7A'
    else:
        color = COL_OTHER_MUT
        edge = '#4A4A4A'

    ms = 14 if m['is_strong'] else 10
    ax4.plot([pos, pos], [backbone_height/2, backbone_height/2 + height],
             color=color, linewidth=2.5, alpha=0.7, zorder=4)
    ax4.scatter(pos, backbone_height/2 + height, s=ms**2,
                color=color, edgecolors=edge, linewidth=2, zorder=5, alpha=0.9)

    if m['is_strong']:
        ax4.text(pos, backbone_height/2 + height + 0.2,
                 f'{m["best_ic50"]:.0f}nM', ha='center', va='bottom',
                 fontsize=FONT_SIZE_ANNOT - 2, fontweight='bold', color=color)

# Labels below
for m in anxa1_muts_unique:
    pos = m['position']
    match_h = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', m['hgvs_p'])
    if match_h:
        label = f"{match_h.group(1)}{match_h.group(2)}{match_h.group(3)}"
    else:
        label = m['hgvs_p']
    ax4.text(pos, -backbone_height/2 - 0.2, label,
             ha='center', va='top', fontsize=FONT_SIZE_ANNOT - 2,
             rotation=45, color='#333333')

ax4.set_xlim(-10, ANXA1_LENGTH + 15)
ax4.set_ylim(-2.2, max_lollipop_height + 2.0)
ax4.set_xlabel('ANXA1 protein position (amino acids)', fontsize=FONT_SIZE_SMALL)
ax4.set_yticks([])
ax4.spines['left'].set_visible(False)
ax4.set_title('D   ANXA1 somatic mutations in SBS2-HIGH cells',
              fontsize=FONT_SIZE, fontweight='bold', loc='left')

legend_d_combined = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=COL_REGION1,
               markersize=14, markeredgecolor='#8B0000', markeredgewidth=2,
               label='Pos 289 region'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=COL_REGION2,
               markersize=14, markeredgecolor='#1A4F7A', markeredgewidth=2,
               label='Pos 321-327 region'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=COL_OTHER_MUT,
               markersize=10, markeredgecolor='#4A4A4A', markeredgewidth=2,
               label='Other positions'),
]
ax4.legend(handles=legend_d_combined, loc='upper left', fontsize=FONT_SIZE_ANNOT,
           frameon=True, fancybox=True, framealpha=0.9, ncol=3)

# Save combined
fig.savefig(os.path.join(OUTPUT_DIR, "Figure7_Combined.pdf"), dpi=300, bbox_inches='tight')
fig.savefig(os.path.join(OUTPUT_DIR, "Figure7_Combined.png"), dpi=300, bbox_inches='tight')
plt.close(fig)

print(f"\nAll panels saved to: {OUTPUT_DIR}", flush=True)
print("Done.", flush=True)
