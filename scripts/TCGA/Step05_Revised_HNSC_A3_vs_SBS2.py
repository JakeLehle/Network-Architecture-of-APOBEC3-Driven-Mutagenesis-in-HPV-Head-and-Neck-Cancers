#!/usr/bin/env python3
"""
Step05_Revised_HNSC_A3_vs_SBS2.py
====================================
Revised Figure 1 pipeline using new SigProfiler v3.4 SBS counts.

Loads the full TCGA expression matrix, extracts A3 family expression,
matches to SBS2 counts from our own SigProfiler pipeline, filters to
HNSC primary tumors, and generates the simplified two-panel Figure 1.

Barcode matching strategy (layered):
  1. DIRECT: WES barcode from original crosswalk file → new SigProfiler
  2. CASE_ID: 12-char patient ID match (TCGA-XX-XXXX), tumor samples only
  3. TRUNC16: 16-char truncation fallback

Panel A: A3A+A3B vs SBS2 scatter (necessary but not sufficient)
Panel B: A3A vs A3B colored by SBS2 (plateau / amplification)

Inputs:
  data/FIG_1/TCGA_master_FPKM_UQ.tsv           (full expression, 3-row header)
  data/FIG_1/TCGA_sample_metadata_final.tsv     (Project_ID for cancer type)
  data/FIG_1/Mutation_Table_Tumors_TCGA.tsv     (for WES barcode crosswalk col)
  SHARED/TCGA/VCF/SigProfiler_output/TCGA_SBS_signature_counts.tsv

Outputs:
  data/FIG_1/HNSC_A3_SBS2_matched_v2.tsv       (matched HNSC data table)
  data/FIG_1/FIGURE_1_PANELS/
    Panel_1a_A3sum_vs_SBS2.pdf/.png
    Panel_1b_A3A_vs_A3B_SBS2.pdf/.png
  data/FIG_1/matching_report.txt

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import spearmanr

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
SHARED_VCF = "/master/jlehle/SHARED/TCGA/VCF"

# Input paths
EXPRESSION_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_1", "TCGA_master_FPKM_UQ.tsv")
METADATA_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_1", "TCGA_sample_metadata_final.tsv")
ORIGINAL_MUT_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_1", "Mutation_Table_Tumors_TCGA.tsv")
NEW_COUNTS_PATH = os.path.join(SHARED_VCF, "SigProfiler_output", "TCGA_SBS_signature_counts.tsv")

# Output paths
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "FIG_1")
PANEL_DIR = os.path.join(OUTPUT_DIR, "FIGURE_1_PANELS")
os.makedirs(PANEL_DIR, exist_ok=True)

# A3 genes
A3_NUCLEAR = ['APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3H']
A3_ALL = ['APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D', 'APOBEC3F', 'APOBEC3G', 'APOBEC3H']

# Figure colors (from project conventions)
COLOR_SBS2_HIGH = "#ed6a5a"   # coral
COLOR_CREAM     = "#f4f1bb"   # cream (LOW SBS2)
COLOR_NORMAL    = "#9bc1bc"   # light teal
FONT_SIZE = 28

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
    """Save figure as both PDF and PNG at 300 DPI."""
    for ext in ['pdf', 'png']:
        fig.savefig(os.path.join(PANEL_DIR, f"{name}.{ext}"),
                    dpi=300, bbox_inches='tight')
    plt.close(fig)
    log(f"  Saved: {name}.pdf/.png")

# =============================================================================
# STEP 1: LOAD EXPRESSION MATRIX
# =============================================================================
banner("STEP 1: Load TCGA expression matrix")

log(f"  Loading: {EXPRESSION_PATH}")
log(f"  (This file has a 3-row header: row0=ENSG IDs, row1=gene symbols, row2=biotypes)")

# Read raw — first 3 rows are header info, data starts at row 4
raw = pd.read_csv(EXPRESSION_PATH, sep='\t', header=None, low_memory=False)
log(f"  Raw shape: {raw.shape}")

# Parse header structure
col_headers = raw.iloc[0].values   # Row 0: column names / ENSG IDs
gene_symbols = raw.iloc[1].values  # Row 1: gene symbols
gene_biotypes = raw.iloc[2].values # Row 2: gene biotypes

# Metadata columns (first 5)
meta_col_names = col_headers[:5]
log(f"  Metadata columns: {list(meta_col_names)}")

# Find Entity_ID column index
entity_id_idx = list(meta_col_names).index('Entity_ID') if 'Entity_ID' in meta_col_names else 4
log(f"  Entity_ID column index: {entity_id_idx}")

# Data rows (row 3 onward)
data = raw.iloc[3:].copy()
data.columns = range(len(data.columns))

# Extract Entity_IDs
entity_ids = data[entity_id_idx].astype(str).values
log(f"  Total samples: {len(entity_ids)}")

# Find A3 gene columns
gene_start_col = 5  # genes start after 5 metadata columns
a3_col_indices = {}
for gene in A3_ALL:
    matches = np.where(gene_symbols[gene_start_col:] == gene)[0]
    if len(matches) > 0:
        col_idx = matches[0] + gene_start_col
        a3_col_indices[gene] = col_idx
        log(f"  Found {gene} at column {col_idx}")
    else:
        log(f"  WARNING: {gene} not found in expression matrix")

# Extract A3 expression for all samples
a3_expr = pd.DataFrame({'Entity_ID': entity_ids})
for gene, col_idx in a3_col_indices.items():
    a3_expr[gene] = pd.to_numeric(data[col_idx].values, errors='coerce')

log(f"  A3 expression table: {a3_expr.shape}")

# =============================================================================
# STEP 2: ADD CANCER TYPE FROM METADATA
# =============================================================================
banner("STEP 2: Add cancer type from metadata")

metadata = pd.read_csv(METADATA_PATH, sep='\t')
log(f"  Metadata: {metadata.shape}")

# Map Entity_ID → Project_ID
meta_map = metadata.set_index('Entity_ID')['Project_ID'].to_dict()
a3_expr['Project_ID'] = a3_expr['Entity_ID'].map(meta_map)
log(f"  Mapped cancer type: {a3_expr['Project_ID'].notna().sum()} / {len(a3_expr)}")

# Filter to HNSC tumors
# HNSC tumors have barcode positions 14-15 = '01' (Primary Tumor)
a3_hnsc = a3_expr[a3_expr['Project_ID'] == 'TCGA-HNSC'].copy()
is_tumor = a3_hnsc['Entity_ID'].str[13:15].isin(['01', '02', '03', '04', '05', '06', '07', '08', '09'])
a3_hnsc = a3_hnsc[is_tumor].copy()
log(f"  HNSC primary tumors with A3 expression: {len(a3_hnsc)}")

# Compute A3A + A3B sum
if 'APOBEC3A' in a3_hnsc.columns and 'APOBEC3B' in a3_hnsc.columns:
    a3_hnsc['A3A_plus_A3B'] = a3_hnsc['APOBEC3A'] + a3_hnsc['APOBEC3B']
    log(f"  A3A+A3B: min={a3_hnsc['A3A_plus_A3B'].min():.1f}  max={a3_hnsc['A3A_plus_A3B'].max():.1f}")

# =============================================================================
# STEP 3: LOAD NEW SBS2 COUNTS
# =============================================================================
banner("STEP 3: Load new SigProfiler SBS2 counts")

log(f"  Loading: {NEW_COUNTS_PATH}")
raw_ct = pd.read_csv(NEW_COUNTS_PATH, sep='\t')

# Handle transposed format
first_col = raw_ct.columns[0]
first_vals = raw_ct[first_col].astype(str).tolist()
if any(v.startswith('SBS') for v in first_vals):
    log(f"  Transposed format — transposing...")
    raw_ct = raw_ct.set_index(first_col)
    new_ct = raw_ct.T.copy()
    new_ct.index.name = 'WES_Barcode'
    new_ct = new_ct.reset_index()
    # Drop non-sample rows
    bad = new_ct['WES_Barcode'].str.contains('Cancer_Type|^$', na=True, regex=True)
    new_ct = new_ct[~bad].copy()
    sbs_cols = [c for c in new_ct.columns if c.startswith('SBS')]
    for c in sbs_cols:
        new_ct[c] = pd.to_numeric(new_ct[c], errors='coerce')
else:
    new_ct = raw_ct
    new_ct = new_ct.rename(columns={first_col: 'WES_Barcode'})
    sbs_cols = [c for c in new_ct.columns if c.startswith('SBS')]

log(f"  New counts: {new_ct.shape}, {len(sbs_cols)} SBS signatures")
log(f"  SBS2: min={new_ct['SBS2'].min():.0f}  max={new_ct['SBS2'].max():.0f}")

# Extract Case_ID (first 12 chars) for matching
new_ct['Case_ID'] = new_ct['WES_Barcode'].str[:12]

# Filter new counts to HNSC using the MuTect2 manifest or Cancer_Type column
if 'Cancer_Type' in new_ct.columns:
    new_ct_hnsc = new_ct[new_ct['Cancer_Type'] == 'HNSC'].copy()
    log(f"  HNSC in new counts (Cancer_Type col): {len(new_ct_hnsc)}")
else:
    # Use manifest
    manifest_path = os.path.join(SHARED_VCF, "manifests", "TCGA_MuTect2_master_manifest.tsv")
    if os.path.exists(manifest_path):
        manifest = pd.read_csv(manifest_path, sep='\t')
        hnsc_wes_ids = set(manifest[manifest['Cancer_Type'] == 'HNSC']['Entity_ID'].astype(str))
        new_ct_hnsc = new_ct[new_ct['WES_Barcode'].isin(hnsc_wes_ids)].copy()
        log(f"  HNSC in new counts (manifest): {len(new_ct_hnsc)}")
    else:
        # Last resort: use Case_IDs from the expression metadata
        hnsc_case_ids = set(a3_hnsc['Entity_ID'].str[:12].values)
        new_ct_hnsc = new_ct[new_ct['Case_ID'].isin(hnsc_case_ids)].copy()
        log(f"  HNSC in new counts (Case_ID match): {len(new_ct_hnsc)}")

# Filter to tumor samples only in new counts (barcode pos 14-15)
new_ct_hnsc_tumor = new_ct_hnsc[new_ct_hnsc['WES_Barcode'].str[13:15].isin(
    ['01', '02', '03', '04', '05', '06', '07', '08', '09'])].copy()
log(f"  HNSC tumor samples in new counts: {len(new_ct_hnsc_tumor)}")

# =============================================================================
# STEP 4: LAYERED BARCODE MATCHING
# =============================================================================
banner("STEP 4: Layered barcode matching")

# Strategy 1: Direct WES barcode from original crosswalk
log(f"\n  --- Strategy 1: Direct WES barcode crosswalk ---")
orig_mut = pd.read_csv(ORIGINAL_MUT_PATH, sep='\t', usecols=[
    'TCGA_Gene_Expression_Entity_ID', 'Mutation_Signature__File_Orginal_Entity_ID'])
crosswalk = orig_mut.rename(columns={
    'TCGA_Gene_Expression_Entity_ID': 'RNA_Barcode',
    'Mutation_Signature__File_Orginal_Entity_ID': 'WES_Barcode'
})
# Filter crosswalk to HNSC RNA-seq samples
hnsc_rna_ids = set(a3_hnsc['Entity_ID'].values)
crosswalk_hnsc = crosswalk[crosswalk['RNA_Barcode'].isin(hnsc_rna_ids)].copy()
log(f"  HNSC crosswalk entries: {len(crosswalk_hnsc)}")

# Match crosswalk WES barcodes to new SigProfiler
new_ct_wes_ids = set(new_ct_hnsc_tumor['WES_Barcode'].values)
crosswalk_hnsc['direct_match'] = crosswalk_hnsc['WES_Barcode'].isin(new_ct_wes_ids)
n_direct = crosswalk_hnsc['direct_match'].sum()
log(f"  Direct WES matches: {n_direct}")

# Build matched table from direct matches
matched_direct = crosswalk_hnsc[crosswalk_hnsc['direct_match']].copy()
matched_direct = matched_direct.merge(
    new_ct_hnsc_tumor[['WES_Barcode', 'SBS2']],
    on='WES_Barcode', how='inner')
matched_direct = matched_direct.rename(columns={'SBS2': 'SBS2_new'})
log(f"  Direct matched with SBS2: {len(matched_direct)}")

# Strategy 2: Case_ID (12-char) matching for unmatched samples
log(f"\n  --- Strategy 2: Case_ID (12-char) matching ---")
already_matched_rna = set(matched_direct['RNA_Barcode'].values)
unmatched_rna = a3_hnsc[~a3_hnsc['Entity_ID'].isin(already_matched_rna)].copy()
log(f"  Unmatched HNSC RNA-seq samples: {len(unmatched_rna)}")

if len(unmatched_rna) > 0:
    unmatched_rna['Case_ID'] = unmatched_rna['Entity_ID'].str[:12]

    # For new counts, deduplicate on Case_ID (keep first tumor sample)
    new_ct_hnsc_dedup = new_ct_hnsc_tumor.drop_duplicates(subset='Case_ID', keep='first')

    case_matched = unmatched_rna[['Entity_ID', 'Case_ID']].merge(
        new_ct_hnsc_dedup[['WES_Barcode', 'Case_ID', 'SBS2']],
        on='Case_ID', how='inner')
    case_matched = case_matched.rename(columns={'Entity_ID': 'RNA_Barcode', 'SBS2': 'SBS2_new'})
    log(f"  Case_ID matches: {len(case_matched)}")
else:
    case_matched = pd.DataFrame(columns=['RNA_Barcode', 'WES_Barcode', 'SBS2_new'])

# Strategy 3: 16-char truncation fallback
log(f"\n  --- Strategy 3: 16-char truncation fallback ---")
all_matched_rna = set(matched_direct['RNA_Barcode'].values) | set(case_matched['RNA_Barcode'].values)
still_unmatched = a3_hnsc[~a3_hnsc['Entity_ID'].isin(all_matched_rna)].copy()
log(f"  Still unmatched: {len(still_unmatched)}")

if len(still_unmatched) > 0:
    still_unmatched['trunc16'] = still_unmatched['Entity_ID'].str[:16]
    new_ct_hnsc_tumor['trunc16'] = new_ct_hnsc_tumor['WES_Barcode'].str[:16]
    new_ct_dedup16 = new_ct_hnsc_tumor.drop_duplicates(subset='trunc16', keep='first')

    trunc_matched = still_unmatched[['Entity_ID', 'trunc16']].merge(
        new_ct_dedup16[['WES_Barcode', 'trunc16', 'SBS2']],
        on='trunc16', how='inner')
    trunc_matched = trunc_matched.rename(columns={'Entity_ID': 'RNA_Barcode', 'SBS2': 'SBS2_new'})
    trunc_matched = trunc_matched.drop(columns=['trunc16'])
    log(f"  16-char truncation matches: {len(trunc_matched)}")
else:
    trunc_matched = pd.DataFrame(columns=['RNA_Barcode', 'WES_Barcode', 'SBS2_new'])

# Combine all matches
log(f"\n  --- COMBINED MATCHING RESULTS ---")
all_matches = pd.concat([
    matched_direct[['RNA_Barcode', 'WES_Barcode', 'SBS2_new']],
    case_matched[['RNA_Barcode', 'WES_Barcode', 'SBS2_new']],
    trunc_matched[['RNA_Barcode', 'WES_Barcode', 'SBS2_new']]
], ignore_index=True)

# Deduplicate (prefer direct match)
all_matches = all_matches.drop_duplicates(subset='RNA_Barcode', keep='first')
log(f"  Total matched HNSC samples: {len(all_matches)}")
log(f"    Direct WES:      {len(matched_direct)}")
log(f"    Case_ID:         {len(case_matched)}")
log(f"    16-char trunc:   {len(trunc_matched)}")

# =============================================================================
# STEP 5: BUILD FINAL MATCHED TABLE
# =============================================================================
banner("STEP 5: Build final matched HNSC table")

# Merge A3 expression with SBS2
final = a3_hnsc.merge(all_matches, left_on='Entity_ID', right_on='RNA_Barcode', how='inner')
log(f"  Final matched table: {len(final)} HNSC tumors with both A3 expression and SBS2")

# Compute derived columns
final['A3A'] = final['APOBEC3A']
final['A3B'] = final['APOBEC3B']
final['A3A_plus_A3B'] = final['A3A'] + final['A3B']
final['SBS2'] = final['SBS2_new']

log(f"\n  A3A:  min={final['A3A'].min():.1f}  max={final['A3A'].max():.1f}  mean={final['A3A'].mean():.1f}")
log(f"  A3B:  min={final['A3B'].min():.1f}  max={final['A3B'].max():.1f}  mean={final['A3B'].mean():.1f}")
log(f"  SBS2: min={final['SBS2'].min():.0f}  max={final['SBS2'].max():.0f}  mean={final['SBS2'].mean():.1f}")
log(f"  Samples with SBS2 > 0: {(final['SBS2'] > 0).sum()}")

# Save matched table
final_path = os.path.join(OUTPUT_DIR, "HNSC_A3_SBS2_matched_v2.tsv")
final.to_csv(final_path, sep='\t', index=False)
log(f"\n  Saved: {final_path}")

# =============================================================================
# STEP 6: GENERATE FIGURE 1 PANELS
# =============================================================================
banner("STEP 6: Generate Figure 1 panels")

plt.rcParams.update({
    'font.size': FONT_SIZE,
    'axes.titlesize': FONT_SIZE + 2,
    'axes.labelsize': FONT_SIZE,
    'xtick.labelsize': FONT_SIZE - 4,
    'ytick.labelsize': FONT_SIZE - 4,
    'legend.fontsize': FONT_SIZE - 6,
})

# --- Panel 1a: A3A+A3B vs SBS2 ---
banner("  Panel 1a: A3A+A3B vs SBS2", char="-")

fig, ax = plt.subplots(figsize=(12, 10))

x = final['A3A_plus_A3B'].values
y = final['SBS2'].values

# Data-driven boundary: steepest SBS2/A3 ratio among low-A3 tumors
low_a3_mask = x <= 10
if low_a3_mask.sum() > 0:
    ratios = y[low_a3_mask] / np.maximum(x[low_a3_mask], 0.01)
    max_ratio = np.max(ratios) * 1.05  # 5% buffer
else:
    max_ratio = 10

median_sbs2 = np.median(y)

# Color assignment: coral=high A3 + high SBS2, cream=high A3 + low SBS2, teal=no high SBS2 without A3
colors = []
for xi, yi in zip(x, y):
    if yi > median_sbs2 and xi > 0:
        colors.append(COLOR_SBS2_HIGH)
    elif xi > 0:
        colors.append(COLOR_CREAM)
    else:
        colors.append(COLOR_NORMAL)

ax.scatter(x, y, c=colors, s=40, alpha=0.7, edgecolors='black', linewidths=0.3, rasterized=True)

# Reference lines
ax.axhline(median_sbs2, color='gray', linestyle='--', linewidth=1, alpha=0.6)
ax.text(ax.get_xlim()[1] * 0.02, median_sbs2 * 1.1, f'Median SBS2 = {median_sbs2:.0f}',
        fontsize=FONT_SIZE - 8, color='gray')

# Count per region
n_high = sum(1 for c in colors if c == COLOR_SBS2_HIGH)
n_low = sum(1 for c in colors if c == COLOR_CREAM)
n_none = sum(1 for c in colors if c == COLOR_NORMAL)

legend_elements = [
    Patch(facecolor=COLOR_SBS2_HIGH, edgecolor='black', label=f'A3 + SBS2 HIGH (n={n_high})'),
    Patch(facecolor=COLOR_CREAM, edgecolor='black', label=f'A3 + SBS2 LOW (n={n_low})'),
    Patch(facecolor=COLOR_NORMAL, edgecolor='black', label=f'Low A3 (n={n_none})'),
]
ax.legend(handles=legend_elements, loc='upper left', framealpha=0.9)

ax.set_xlabel('A3A + A3B Expression (FPKM-UQ)')
ax.set_ylabel('SBS2 Weight (mutation count)')
ax.set_title('A3 Expression vs SBS2 Mutagenesis\nin HNSCC (TCGA)')

rho, p = spearmanr(x, y)
ax.text(0.98, 0.02, f'Spearman rho = {rho:.3f}\np = {p:.2e}\nn = {len(final)}',
        transform=ax.transAxes, ha='right', va='bottom',
        fontsize=FONT_SIZE - 8, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()
save_fig(fig, "Panel_1a_A3sum_vs_SBS2")

# --- Panel 1b: A3A vs A3B colored by SBS2 ---
banner("  Panel 1b: A3A vs A3B colored by SBS2", char="-")

fig, ax = plt.subplots(figsize=(12, 10))

x2 = final['A3A'].values
y2 = final['A3B'].values
c2 = final['SBS2'].values

# Use log1p for color to handle zeros and wide range
c2_log = np.log1p(c2)

scatter = ax.scatter(x2, y2, c=c2_log, cmap='magma', s=40, alpha=0.7,
                     edgecolors='black', linewidths=0.3, rasterized=True)

cbar = plt.colorbar(scatter, ax=ax, shrink=0.8, pad=0.02)
# Set colorbar ticks to original scale
c_ticks = [0, 1, 5, 20, 50, 100, 200, 500]
c_ticks = [t for t in c_ticks if t <= c2.max()]
cbar.set_ticks([np.log1p(t) for t in c_ticks])
cbar.set_ticklabels([str(t) for t in c_ticks])
cbar.set_label('SBS2 Weight', fontsize=FONT_SIZE - 4)

ax.set_xlabel('APOBEC3A Expression (FPKM-UQ)')
ax.set_ylabel('APOBEC3B Expression (FPKM-UQ)')
ax.set_title('A3A vs A3B Expression\nColored by SBS2 Weight')

plt.tight_layout()
save_fig(fig, "Panel_1b_A3A_vs_A3B_SBS2")

# =============================================================================
# STEP 7: SUMMARY
# =============================================================================
banner("SUMMARY")

log(f"""
  INPUT:
    Expression matrix: {EXPRESSION_PATH}
    New SBS2 counts:   {NEW_COUNTS_PATH}
    Metadata:          {METADATA_PATH}
    Crosswalk:         {ORIGINAL_MUT_PATH}

  MATCHING:
    HNSC RNA-seq tumors:    {len(a3_hnsc)}
    HNSC WES tumors:        {len(new_ct_hnsc_tumor)}
    Matched (final):        {len(final)}
      Direct WES barcode:   {len(matched_direct)}
      Case_ID 12-char:      {len(case_matched)}
      16-char truncation:   {len(trunc_matched)}
    Unmatched RNA-seq:      {len(a3_hnsc) - len(final)}

  DATA:
    SBS2 range:  {final['SBS2'].min():.0f} - {final['SBS2'].max():.0f}
    SBS2 > 0:    {(final['SBS2'] > 0).sum()} / {len(final)}
    A3A range:   {final['A3A'].min():.1f} - {final['A3A'].max():.1f}
    A3B range:   {final['A3B'].min():.1f} - {final['A3B'].max():.1f}

  OUTPUT:
    Data:   {final_path}
    Panels: {PANEL_DIR}/Panel_1a_A3sum_vs_SBS2.pdf/.png
            {PANEL_DIR}/Panel_1b_A3A_vs_A3B_SBS2.pdf/.png
""")

# Confirm necessary-but-not-sufficient pattern
high_sbs2 = final[final['SBS2'] > median_sbs2]
low_a3_high_sbs2 = high_sbs2[high_sbs2['A3A_plus_A3B'] < 1.0]
log(f"  KEY CHECK: Tumors with SBS2 > median AND A3A+A3B < 1.0: {len(low_a3_high_sbs2)}")
if len(low_a3_high_sbs2) == 0:
    log(f"  → CONFIRMED: No tumors have high SBS2 without A3 expression")
    log(f"  → A3 expression is NECESSARY but NOT SUFFICIENT for SBS2")
else:
    log(f"  → WARNING: {len(low_a3_high_sbs2)} tumors have high SBS2 without A3")
    log(f"  → Inspect these samples manually")

# Save report
report_path = os.path.join(OUTPUT_DIR, "TROUBLESHOOTING", "figure1_v2_matching_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"\n  Report: {report_path}")

banner("FIGURE 1 PIPELINE COMPLETE")
