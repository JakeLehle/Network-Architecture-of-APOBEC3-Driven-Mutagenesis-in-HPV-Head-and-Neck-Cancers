#!/usr/bin/env python3
"""
Generate_Figure7_Panels.py
============================
Figure 7: Neoantigen Landscape of HPV16 Lifecycle Populations

Panel A: Neoantigen burden comparison (SBS2-HIGH vs CNV-HIGH)
Panel B: Therapeutic tier classification with key targets
Panel C: ANXA1 protein lollipop plot (mutation positions, APOBEC context, binding)
Panel D: ANXA1 top neoantigen peptides WT vs mutant IC50 comparison

Run in NETWORK conda env.
Save as PDF and PNG at 300 DPI.

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
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

# Style
FONT_SIZE = 28
FONT_SIZE_SMALL = 24
FONT_SIZE_TICK = 22
FONT_SIZE_LABEL = 20

# Colors (hex only)
COL_SBS2 = '#E74C3C'
COL_CNV = '#3498DB'
COL_NORMAL = '#95A5A6'
COL_SHARED = '#9B59B6'
COL_APOBEC = '#E74C3C'
COL_NON_APOBEC = '#2C3E50'
COL_STRONG = '#E74C3C'
COL_WEAK = '#F39C12'
COL_NONBIND = '#BDC3C7'
COL_WT = '#95A5A6'

# ANXA1 domain boundaries (UniProt P04083)
ANXA1_LENGTH = 346
ANXA1_DOMAINS = [
    (1, 41, 'N-term', '#F5B7B1'),
    (42, 109, 'Repeat 1', '#AED6F1'),
    (123, 196, 'Repeat 2', '#A9DFBF'),
    (221, 287, 'Repeat 3', '#F9E79F'),
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

sbs2_rank = pd.read_csv(os.path.join(MHC_DIR, "SBS2_HIGH_expression_weighted_ranking.tsv"), sep='\t')
cnv_rank = pd.read_csv(os.path.join(MHC_DIR, "CNV_HIGH_expression_weighted_ranking.tsv"), sep='\t')
tier_df = pd.read_csv(os.path.join(TIER_DIR, "all_genes_tiered.tsv"), sep='\t')

sbs2_neo = pd.read_csv(os.path.join(MHC_DIR, "SBS2_HIGH_neoantigens.tsv"), sep='\t')
anxa1_neo = sbs2_neo[sbs2_neo['gene'] == 'ANXA1'].copy()

sbs2_all_pep = pd.read_csv(os.path.join(MHC_DIR, "SBS2_HIGH_all_peptide_results.tsv"), sep='\t')
anxa1_all_pep = sbs2_all_pep[sbs2_all_pep['gene'] == 'ANXA1'].copy()

sbs2_spa = pd.read_csv(os.path.join(ANNOTATION_DIR, "SBS2_HIGH.somatic_protein_altering.tsv"), sep='\t')
anxa1_variants = sbs2_spa[sbs2_spa['gene'] == 'ANXA1'].copy()

cnv_neo = pd.read_csv(os.path.join(MHC_DIR, "CNV_HIGH_neoantigens.tsv"), sep='\t')

n_cells = 546
sbs2_total_neo = len(sbs2_neo)
cnv_total_neo = len(cnv_neo)
sbs2_strong = sbs2_neo['is_strong_binder'].sum() if 'is_strong_binder' in sbs2_neo.columns else 0
cnv_strong = cnv_neo['is_strong_binder'].sum() if 'is_strong_binder' in cnv_neo.columns else 0
sbs2_diff = sbs2_neo['is_differential'].sum() if 'is_differential' in sbs2_neo.columns else 0
cnv_diff = cnv_neo['is_differential'].sum() if 'is_differential' in cnv_neo.columns else 0

print(f"  SBS2 neoantigens: {sbs2_total_neo}, CNV neoantigens: {cnv_total_neo}", flush=True)
print(f"  ANXA1 neoantigen peptides: {len(anxa1_neo)}", flush=True)
print(f"  ANXA1 all peptides: {len(anxa1_all_pep)}", flush=True)
print(f"  ANXA1 somatic variants: {len(anxa1_variants)}", flush=True)

# =============================================================================
# PANEL A: NEOANTIGEN BURDEN COMPARISON
# =============================================================================
print("Panel A: Neoantigen burden...", flush=True)

fig_a, ax_a = plt.subplots(figsize=(14, 10))

categories = ['Neoantigens\n(IC50 < 500nM)', 'Strong binders\n(IC50 < 50nM)', 'Differential\n(new epitope)']
sbs2_vals = [sbs2_total_neo / n_cells, sbs2_strong / n_cells, sbs2_diff / n_cells]
cnv_vals = [cnv_total_neo / n_cells, cnv_strong / n_cells, cnv_diff / n_cells]

x = np.arange(len(categories))
width = 0.35

bars1 = ax_a.bar(x - width/2, sbs2_vals, width, color=COL_SBS2, label='SBS2-HIGH', edgecolor='white', linewidth=1.5)
bars2 = ax_a.bar(x + width/2, cnv_vals, width, color=COL_CNV, label='CNV-HIGH', edgecolor='white', linewidth=1.5)

for bar in bars1:
    height = bar.get_height()
    ax_a.text(bar.get_x() + bar.get_width()/2., height + 0.05,
              f'{height:.2f}', ha='center', va='bottom', fontsize=FONT_SIZE_LABEL, fontweight='bold')
for bar in bars2:
    height = bar.get_height()
    ax_a.text(bar.get_x() + bar.get_width()/2., height + 0.05,
              f'{height:.2f}', ha='center', va='bottom', fontsize=FONT_SIZE_LABEL, fontweight='bold')

for i in range(len(categories)):
    if cnv_vals[i] > 0:
        fold = sbs2_vals[i] / cnv_vals[i]
        mid_x = x[i]
        max_y = max(sbs2_vals[i], cnv_vals[i])
        ax_a.text(mid_x, max_y + 0.25, f'{fold:.1f}x', ha='center', va='bottom',
                  fontsize=FONT_SIZE_LABEL, fontstyle='italic', color='#555555')

ax_a.set_ylabel('Per cell', fontsize=FONT_SIZE_SMALL)
ax_a.set_xticks(x)
ax_a.set_xticklabels(categories, fontsize=FONT_SIZE_TICK)
ax_a.legend(fontsize=FONT_SIZE_LABEL, frameon=False, loc='upper right')
ax_a.set_title('Predicted neoantigen burden', fontsize=FONT_SIZE, fontweight='bold', pad=20)
ax_a.set_ylim(0, max(sbs2_vals) * 1.35)

fig_a.tight_layout()
fig_a.savefig(os.path.join(OUTPUT_DIR, "Figure7_PanelA_Burden.pdf"), dpi=300, bbox_inches='tight')
fig_a.savefig(os.path.join(OUTPUT_DIR, "Figure7_PanelA_Burden.png"), dpi=300, bbox_inches='tight')
plt.close(fig_a)
print("  Panel A saved.", flush=True)

# =============================================================================
# PANEL B: THERAPEUTIC TIER CLASSIFICATION
# =============================================================================
print("Panel B: Therapeutic tiers...", flush=True)

fig_b, ax_b = plt.subplots(figsize=(18, 12))
ax_b.set_xlim(0, 100)
ax_b.set_ylim(0, 100)
ax_b.axis('off')

tier_info = [
    {
        'name': 'Tier 1A: Immune-recognized shared targets',
        'subtitle': 'Shared neoantigens with CNV escape evidence',
        'count': len(tier_df[tier_df['tier'] == '1A_hot_shared_escaped']),
        'genes': 'KRT6A, KRT5, PI3, SERPINB2, HLA-A, HLA-C, SPRR3',
        'color': COL_SBS2,
        'y': 82, 'strategy': 'Hot tumors (high SBS2:CNV)',
    },
    {
        'name': 'Tier 1B: Maintenance-phase specific',
        'subtitle': 'SBS2-HIGH only neoantigens',
        'count': len(tier_df[tier_df['tier'] == '1B_hot_sbs2_specific']),
        'genes': 'ANXA1, MALL, EIF4A1, S100A8, LMO7, HSP90AB1',
        'color': '#C0392B',
        'y': 60, 'strategy': 'Hot tumors (high SBS2:CNV)',
    },
    {
        'name': 'Tier 2: Productive-phase specific',
        'subtitle': 'CNV-HIGH only neoantigens',
        'count': len(tier_df[tier_df['tier'] == '2_cold_cnv_specific']),
        'genes': 'COX4I1, S100A11, ENO1, HEXB, PCBP1, CLIC1',
        'color': COL_CNV,
        'y': 38, 'strategy': 'Cold tumors + checkpoint blockade',
    },
    {
        'name': 'Tier 3: Broad coverage backbone',
        'subtitle': 'Shared, no immune escape evidence',
        'count': len(tier_df[tier_df['tier'] == '3_broad_coverage']),
        'genes': 'MDK, SERPINB5, BTF3, ACTG1, AQP3, RPL13A',
        'color': COL_SHARED,
        'y': 16, 'strategy': 'All tumors (vaccine backbone)',
    },
]

for t in tier_info:
    y = t['y']
    ax_b.add_patch(FancyBboxPatch((2, y - 5), 3, 10,
                                   boxstyle="round,pad=0.3",
                                   facecolor=t['color'], edgecolor='none'))
    ax_b.text(8, y + 2, f"{t['name']} (n={t['count']})",
              fontsize=FONT_SIZE_SMALL, fontweight='bold', va='center')
    ax_b.text(8, y - 2.5, t['genes'],
              fontsize=FONT_SIZE_LABEL, va='center', color='#555555', style='italic')
    ax_b.text(95, y, t['strategy'],
              fontsize=FONT_SIZE_LABEL, va='center', ha='right',
              color=t['color'], fontweight='bold',
              bbox=dict(boxstyle='round,pad=0.3', facecolor=t['color'],
                        alpha=0.1, edgecolor=t['color'], linewidth=1.5))

ax_b.set_title('Therapeutic tier classification', fontsize=FONT_SIZE, fontweight='bold', pad=20)

fig_b.tight_layout()
fig_b.savefig(os.path.join(OUTPUT_DIR, "Figure7_PanelB_Tiers.pdf"), dpi=300, bbox_inches='tight')
fig_b.savefig(os.path.join(OUTPUT_DIR, "Figure7_PanelB_Tiers.png"), dpi=300, bbox_inches='tight')
plt.close(fig_b)
print("  Panel B saved.", flush=True)

# =============================================================================
# PANEL C: ANXA1 LOLLIPOP PLOT
# =============================================================================
print("Panel C: ANXA1 lollipop...", flush=True)

anxa1_muts = []
for _, var in anxa1_variants.iterrows():
    hgvs_p = str(var.get('hgvs_p', ''))
    match = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', hgvs_p)
    if not match:
        continue

    pos = int(match.group(2))
    ref_dna = str(var.get('ref', ''))
    alt_dna = str(var.get('alt', ''))

    is_apobec_type = False
    if (ref_dna == 'C' and alt_dna in ['T', 'G']) or \
       (ref_dna == 'G' and alt_dna in ['A', 'C']):
        is_apobec_type = True

    pos_neo = anxa1_neo[anxa1_neo['hgvs_p'] == hgvs_p] if 'hgvs_p' in anxa1_neo.columns else pd.DataFrame()
    if len(pos_neo) == 0 and 'mut_pos_protein' in anxa1_neo.columns:
        pos_neo = anxa1_neo[anxa1_neo['mut_pos_protein'] == pos]

    best_ic50 = pos_neo['mut_ic50'].min() if len(pos_neo) > 0 else 999

    anxa1_muts.append({
        'position': pos,
        'hgvs_p': hgvs_p,
        'is_apobec_type': is_apobec_type,
        'best_ic50': best_ic50,
        'is_strong': best_ic50 < 50,
        'is_binder': best_ic50 < 500,
    })

pos_to_mut = {}
for m in anxa1_muts:
    p = m['position']
    if p not in pos_to_mut or m['best_ic50'] < pos_to_mut[p]['best_ic50']:
        pos_to_mut[p] = m
anxa1_muts_unique = sorted(pos_to_mut.values(), key=lambda x: x['position'])

n_apobec = sum(1 for m in anxa1_muts_unique if m['is_apobec_type'])
print(f"  ANXA1 unique positions: {len(anxa1_muts_unique)}, APOBEC-type: {n_apobec}", flush=True)

fig_c, ax_c = plt.subplots(figsize=(20, 8))

backbone_y = 0
backbone_height = 0.8
max_lollipop_height = 6

for start, end, name, color in ANXA1_DOMAINS:
    ax_c.add_patch(plt.Rectangle((start, backbone_y - backbone_height/2),
                                  end - start, backbone_height,
                                  facecolor=color, edgecolor='#333333',
                                  linewidth=1.5, zorder=2))
    mid = (start + end) / 2
    ax_c.text(mid, backbone_y, name, ha='center', va='center',
              fontsize=FONT_SIZE_LABEL - 2, fontweight='bold', zorder=3)

for m in anxa1_muts_unique:
    pos = m['position']
    if m['best_ic50'] < 500:
        height = max_lollipop_height * (500 - m['best_ic50']) / 500
        height = max(height, 1.5)
    else:
        height = 1.0

    color = COL_APOBEC if m['is_apobec_type'] else COL_NON_APOBEC
    marker_edge = '#8B0000' if m['is_apobec_type'] else '#1A1A2E'
    marker_size = 14 if m['is_strong'] else 10

    ax_c.plot([pos, pos], [backbone_height/2, backbone_height/2 + height],
              color=color, linewidth=2, alpha=0.7, zorder=4)
    ax_c.scatter(pos, backbone_height/2 + height, s=marker_size**2,
                 color=color, edgecolors=marker_edge, linewidth=1.5,
                 zorder=5, alpha=0.9)

ax_c.set_xlim(-5, ANXA1_LENGTH + 10)
ax_c.set_ylim(-2, max_lollipop_height + 2.5)
ax_c.set_xlabel('ANXA1 protein position (amino acids)', fontsize=FONT_SIZE_SMALL)
ax_c.set_ylabel('MHC-I binding\nstrength', fontsize=FONT_SIZE_SMALL)
ax_c.set_title(f'ANXA1 somatic mutations in SBS2-HIGH cells ({len(anxa1_muts_unique)} unique positions)',
               fontsize=FONT_SIZE, fontweight='bold', pad=20)
ax_c.set_yticks([])
ax_c.text(-3, backbone_height/2 + max_lollipop_height, 'Strong', va='center', ha='right',
          fontsize=FONT_SIZE_LABEL, color='#555555')
ax_c.text(-3, backbone_height/2 + 1.0, 'Weak', va='center', ha='right',
          fontsize=FONT_SIZE_LABEL, color='#555555')

legend_elements = [
    mpatches.Patch(facecolor=COL_APOBEC, edgecolor='#8B0000', label=f'APOBEC-type C>T/C>G (n={n_apobec})'),
    mpatches.Patch(facecolor=COL_NON_APOBEC, edgecolor='#1A1A2E',
                   label=f'Other mutations (n={len(anxa1_muts_unique)-n_apobec})'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#555555',
               markersize=14, label='Strong binder (IC50 < 50nM)'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#555555',
               markersize=9, label='Binder (IC50 < 500nM)'),
]
ax_c.legend(handles=legend_elements, loc='upper right', fontsize=FONT_SIZE_LABEL - 2,
            frameon=True, fancybox=True, framealpha=0.9)
ax_c.spines['left'].set_visible(False)

fig_c.tight_layout()
fig_c.savefig(os.path.join(OUTPUT_DIR, "Figure7_PanelC_ANXA1_Lollipop.pdf"), dpi=300, bbox_inches='tight')
fig_c.savefig(os.path.join(OUTPUT_DIR, "Figure7_PanelC_ANXA1_Lollipop.png"), dpi=300, bbox_inches='tight')
plt.close(fig_c)
print("  Panel C saved.", flush=True)

# =============================================================================
# PANEL D: ANXA1 TOP PEPTIDES WT vs MUTANT IC50
# =============================================================================
print("Panel D: ANXA1 peptide comparison...", flush=True)

anxa1_binders = anxa1_all_pep[anxa1_all_pep['is_binder'] == True].copy()

if len(anxa1_binders) > 0:
    anxa1_binders['ic50_shift'] = anxa1_binders['wt_ic50'] - anxa1_binders['mut_ic50']
    anxa1_binders = anxa1_binders.sort_values('mut_ic50')
    anxa1_uniq = anxa1_binders.drop_duplicates(subset=['mut_peptide'], keep='first')
    top_peptides = anxa1_uniq.head(15).copy()
    top_peptides = top_peptides.sort_values('mut_ic50', ascending=False)

    fig_d, ax_d = plt.subplots(figsize=(16, 12))
    y_pos = np.arange(len(top_peptides))

    ax_d.barh(y_pos, top_peptides['wt_ic50'].clip(upper=2000),
              height=0.6, color=COL_WT, alpha=0.5, label='Wild-type IC50',
              edgecolor='white', linewidth=0.5)

    mut_colors = [COL_STRONG if ic50 < 50 else COL_WEAK for ic50 in top_peptides['mut_ic50']]
    ax_d.barh(y_pos, top_peptides['mut_ic50'],
              height=0.6, color=mut_colors, alpha=0.9, label='Mutant IC50',
              edgecolor='white', linewidth=0.5)

    ax_d.axvline(x=500, color='#333333', linestyle='--', linewidth=2, alpha=0.5, zorder=0)
    ax_d.text(510, len(top_peptides) - 0.5, 'Binding\nthreshold', fontsize=FONT_SIZE_LABEL - 2,
              va='top', color='#555555')

    ax_d.axvline(x=50, color=COL_STRONG, linestyle=':', linewidth=2, alpha=0.4, zorder=0)
    ax_d.text(55, len(top_peptides) - 0.5, 'Strong\nbinder', fontsize=FONT_SIZE_LABEL - 2,
              va='top', color=COL_STRONG)

    y_labels = []
    for _, row in top_peptides.iterrows():
        pep = row.get('mut_peptide', '')
        allele = row.get('best_allele', '')
        allele_short = allele.replace('HLA-', '') if isinstance(allele, str) else ''
        y_labels.append(f"{pep}  ({allele_short})")

    ax_d.set_yticks(y_pos)
    ax_d.set_yticklabels(y_labels, fontsize=FONT_SIZE_LABEL - 2, fontfamily='monospace')
    ax_d.set_xlabel('Predicted IC50 (nM)', fontsize=FONT_SIZE_SMALL)
    ax_d.set_title('ANXA1 top neoantigen peptides: mutant vs wild-type MHC-I binding',
                   fontsize=FONT_SIZE, fontweight='bold', pad=20)
    ax_d.set_xlim(0, min(2000, top_peptides['wt_ic50'].max() * 1.1))

    legend_elements = [
        mpatches.Patch(facecolor=COL_STRONG, alpha=0.9, label='Mutant (strong, IC50 < 50nM)'),
        mpatches.Patch(facecolor=COL_WEAK, alpha=0.9, label='Mutant (binder, IC50 < 500nM)'),
        mpatches.Patch(facecolor=COL_WT, alpha=0.5, label='Wild-type IC50'),
    ]
    ax_d.legend(handles=legend_elements, loc='lower right', fontsize=FONT_SIZE_LABEL,
                frameon=True, fancybox=True, framealpha=0.9)

    fig_d.tight_layout()
    fig_d.savefig(os.path.join(OUTPUT_DIR, "Figure7_PanelD_ANXA1_Peptides.pdf"), dpi=300, bbox_inches='tight')
    fig_d.savefig(os.path.join(OUTPUT_DIR, "Figure7_PanelD_ANXA1_Peptides.png"), dpi=300, bbox_inches='tight')
    plt.close(fig_d)
    print(f"  Panel D saved ({len(top_peptides)} peptides).", flush=True)
else:
    print("  WARNING: No ANXA1 binder peptides found for Panel D", flush=True)

# =============================================================================
# COMBINED FIGURE
# =============================================================================
print("Combined figure...", flush=True)

fig = plt.figure(figsize=(36, 32))
gs = gridspec.GridSpec(2, 2, hspace=0.3, wspace=0.25)

# Panel A (top left)
ax1 = fig.add_subplot(gs[0, 0])
categories_short = ['Neoantigens', 'Strong\nbinders', 'Differential']
x = np.arange(len(categories_short))
width = 0.35
b1 = ax1.bar(x - width/2, sbs2_vals, width, color=COL_SBS2, label='SBS2-HIGH', edgecolor='white', linewidth=1.5)
b2 = ax1.bar(x + width/2, cnv_vals, width, color=COL_CNV, label='CNV-HIGH', edgecolor='white', linewidth=1.5)
for bar in b1:
    h = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., h + 0.03, f'{h:.2f}', ha='center', va='bottom',
             fontsize=FONT_SIZE_LABEL, fontweight='bold')
for bar in b2:
    h = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., h + 0.03, f'{h:.2f}', ha='center', va='bottom',
             fontsize=FONT_SIZE_LABEL, fontweight='bold')
ax1.set_ylabel('Per cell', fontsize=FONT_SIZE_SMALL)
ax1.set_xticks(x)
ax1.set_xticklabels(categories_short, fontsize=FONT_SIZE_TICK)
ax1.legend(fontsize=FONT_SIZE_LABEL, frameon=False)
ax1.set_title('A   Predicted neoantigen burden', fontsize=FONT_SIZE, fontweight='bold', loc='left')
ax1.set_ylim(0, max(sbs2_vals) * 1.3)

# Panel B (top right)
ax2 = fig.add_subplot(gs[0, 1])
ax2.set_xlim(0, 100)
ax2.set_ylim(0, 100)
ax2.axis('off')
for t in tier_info:
    y = t['y']
    ax2.add_patch(FancyBboxPatch((1, y - 5), 3, 10,
                                  boxstyle="round,pad=0.3",
                                  facecolor=t['color'], edgecolor='none'))
    ax2.text(7, y + 2.5, f"{t['name']} (n={t['count']})",
             fontsize=FONT_SIZE_LABEL, fontweight='bold', va='center')
    ax2.text(7, y - 2, t['genes'],
             fontsize=FONT_SIZE_LABEL - 4, va='center', color='#555555', style='italic')
    ax2.text(97, y, t['strategy'],
             fontsize=FONT_SIZE_LABEL - 4, va='center', ha='right',
             color=t['color'], fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.3', facecolor=t['color'],
                       alpha=0.1, edgecolor=t['color'], linewidth=1.5))
ax2.set_title('B   Therapeutic tier classification', fontsize=FONT_SIZE, fontweight='bold', loc='left')

# Panel C (bottom left)
ax3 = fig.add_subplot(gs[1, 0])
for start, end, name, color in ANXA1_DOMAINS:
    ax3.add_patch(plt.Rectangle((start, -backbone_height/2), end - start, backbone_height,
                                 facecolor=color, edgecolor='#333333', linewidth=1.5, zorder=2))
    ax3.text((start + end)/2, 0, name, ha='center', va='center',
             fontsize=FONT_SIZE_LABEL - 4, fontweight='bold', zorder=3)
for m in anxa1_muts_unique:
    pos = m['position']
    if m['best_ic50'] < 500:
        height = max_lollipop_height * (500 - m['best_ic50']) / 500
        height = max(height, 1.5)
    else:
        height = 1.0
    color = COL_APOBEC if m['is_apobec_type'] else COL_NON_APOBEC
    marker_edge = '#8B0000' if m['is_apobec_type'] else '#1A1A2E'
    ms = 12 if m['is_strong'] else 8
    ax3.plot([pos, pos], [backbone_height/2, backbone_height/2 + height],
             color=color, linewidth=2, alpha=0.7, zorder=4)
    ax3.scatter(pos, backbone_height/2 + height, s=ms**2,
                color=color, edgecolors=marker_edge, linewidth=1.5, zorder=5, alpha=0.9)
ax3.set_xlim(-5, ANXA1_LENGTH + 10)
ax3.set_ylim(-2, max_lollipop_height + 2)
ax3.set_xlabel('ANXA1 protein position (aa)', fontsize=FONT_SIZE_SMALL)
ax3.set_yticks([])
ax3.spines['left'].set_visible(False)
ax3.set_title(f'C   ANXA1 mutations in SBS2-HIGH ({len(anxa1_muts_unique)} positions)',
              fontsize=FONT_SIZE, fontweight='bold', loc='left')
legend_c = [
    mpatches.Patch(facecolor=COL_APOBEC, edgecolor='#8B0000', label=f'APOBEC-type (n={n_apobec})'),
    mpatches.Patch(facecolor=COL_NON_APOBEC, edgecolor='#1A1A2E',
                   label=f'Other (n={len(anxa1_muts_unique)-n_apobec})'),
]
ax3.legend(handles=legend_c, loc='upper right', fontsize=FONT_SIZE_LABEL - 4, frameon=True)

# Panel D (bottom right)
if len(anxa1_binders) > 0:
    ax4 = fig.add_subplot(gs[1, 1])
    top_pep = anxa1_uniq.head(12).sort_values('mut_ic50', ascending=False)
    y_pos = np.arange(len(top_pep))
    ax4.barh(y_pos, top_pep['wt_ic50'].clip(upper=2000),
             height=0.6, color=COL_WT, alpha=0.5, edgecolor='white', linewidth=0.5)
    mut_c = [COL_STRONG if ic < 50 else COL_WEAK for ic in top_pep['mut_ic50']]
    ax4.barh(y_pos, top_pep['mut_ic50'],
             height=0.6, color=mut_c, alpha=0.9, edgecolor='white', linewidth=0.5)
    ax4.axvline(x=500, color='#333333', linestyle='--', linewidth=2, alpha=0.5, zorder=0)
    ax4.axvline(x=50, color=COL_STRONG, linestyle=':', linewidth=2, alpha=0.4, zorder=0)
    y_labs = []
    for _, row in top_pep.iterrows():
        pep = row.get('mut_peptide', '')
        allele_short = str(row.get('best_allele', '')).replace('HLA-', '')
        y_labs.append(f"{pep} ({allele_short})")
    ax4.set_yticks(y_pos)
    ax4.set_yticklabels(y_labs, fontsize=FONT_SIZE_LABEL - 4, fontfamily='monospace')
    ax4.set_xlabel('Predicted IC50 (nM)', fontsize=FONT_SIZE_SMALL)
    ax4.set_title('D   ANXA1 peptides: WT vs mutant binding',
                  fontsize=FONT_SIZE, fontweight='bold', loc='left')
    ax4.set_xlim(0, min(2000, top_pep['wt_ic50'].max() * 1.1))
    legend_d = [
        mpatches.Patch(facecolor=COL_STRONG, alpha=0.9, label='Mutant (strong)'),
        mpatches.Patch(facecolor=COL_WEAK, alpha=0.9, label='Mutant (binder)'),
        mpatches.Patch(facecolor=COL_WT, alpha=0.5, label='Wild-type'),
    ]
    ax4.legend(handles=legend_d, loc='lower right', fontsize=FONT_SIZE_LABEL - 4, frameon=True)

fig.savefig(os.path.join(OUTPUT_DIR, "Figure7_Combined.pdf"), dpi=300, bbox_inches='tight')
fig.savefig(os.path.join(OUTPUT_DIR, "Figure7_Combined.png"), dpi=300, bbox_inches='tight')
plt.close(fig)

print(f"\nAll panels saved to: {OUTPUT_DIR}", flush=True)
print("Done.", flush=True)
