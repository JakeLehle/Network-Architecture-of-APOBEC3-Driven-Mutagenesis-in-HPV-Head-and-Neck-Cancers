#!/usr/bin/env python3
"""
Diagnostic_Figure7_Redesign.py
================================
Pull all numbers needed for Figure 7 redesign before writing figure code.

Sections:
  1. Panel A: Neoantigen burden (per cell) + RNA fusion burden (per cell)
  2. Panel B: Venn diagram gene sets (SBS2 vs CNV neoantigen genes)
  3. Panel C: ANXA1 top peptides (WT vs mutant IC50, identify key regions)
  4. Panel D: ANXA1 mutation positions + domain mapping for gene track

Run in NETWORK conda env.
"""

import os
import pandas as pd
import numpy as np

PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
MHC_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/03_mhc_binding")
FUSION_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/04_fusion_analysis")
SUMMARY_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/05_summary")
TIER_DIR = os.path.join(SUMMARY_DIR, "THERAPEUTIC_TIERS")
ANNOTATION_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/02_snpeff_annotation")

N_CELLS = 546  # per group

sep = "=" * 80

# =============================================================================
# 1. PANEL A: NEOANTIGEN + FUSION BURDEN
# =============================================================================
print(f"\n{sep}")
print("PANEL A: NEOANTIGEN + FUSION BURDEN (per cell)")
print(sep)

# --- Neoantigens ---
print("\n--- Neoantigens (germline-subtracted, IC50 < 500nM) ---")
for group in ['SBS2_HIGH', 'CNV_HIGH']:
    neo_file = os.path.join(MHC_DIR, f"{group}_neoantigens.tsv")
    if os.path.exists(neo_file):
        neo = pd.read_csv(neo_file, sep='\t')
        total = len(neo)
        strong = neo['is_strong_binder'].sum() if 'is_strong_binder' in neo.columns else 'N/A'
        diff = neo['is_differential'].sum() if 'is_differential' in neo.columns else 'N/A'
        unique_genes = neo['gene'].nunique() if 'gene' in neo.columns else 'N/A'
        print(f"\n  {group}:")
        print(f"    Total neoantigens (IC50 < 500nM): {total}")
        print(f"    Per cell: {total / N_CELLS:.2f}")
        print(f"    Strong binders (IC50 < 50nM): {strong}")
        print(f"    Strong per cell: {strong / N_CELLS:.2f}" if isinstance(strong, (int, float)) else "")
        print(f"    Differential (new epitope): {diff}")
        print(f"    Diff per cell: {diff / N_CELLS:.2f}" if isinstance(diff, (int, float)) else "")
        print(f"    Unique genes: {unique_genes}")
    else:
        print(f"  {group}: FILE NOT FOUND: {neo_file}")

# --- RNA Fusions ---
print("\n--- RNA Fusions ---")

# Check what fusion files exist
print("\n  Available fusion files:")
if os.path.exists(FUSION_DIR):
    for f in sorted(os.listdir(FUSION_DIR)):
        fpath = os.path.join(FUSION_DIR, f)
        size = os.path.getsize(fpath) / 1024
        print(f"    {f} ({size:.1f} KB)")
else:
    print(f"  FUSION_DIR NOT FOUND: {FUSION_DIR}")

# Try per_group_junction_summary first
junction_summary = os.path.join(FUSION_DIR, "per_group_junction_summary.tsv")
if os.path.exists(junction_summary):
    print(f"\n  per_group_junction_summary.tsv:")
    js = pd.read_csv(junction_summary, sep='\t')
    print(f"    Columns: {list(js.columns)}")
    print(f"    Shape: {js.shape}")
    print(js.to_string(index=False))

# Try all_filtered_junctions for per-cell calculation
all_junctions = os.path.join(FUSION_DIR, "all_filtered_junctions.tsv")
if os.path.exists(all_junctions):
    print(f"\n  all_filtered_junctions.tsv:")
    aj = pd.read_csv(all_junctions, sep='\t')
    print(f"    Columns: {list(aj.columns)}")
    print(f"    Shape: {aj.shape}")
    print(f"    First 3 rows:")
    print(aj.head(3).to_string(index=False))
    
    # Try to compute per-group totals
    if 'group' in aj.columns:
        print(f"\n    Per-group junction counts:")
        for grp in aj['group'].unique():
            n = len(aj[aj['group'] == grp])
            print(f"      {grp}: {n} junctions, {n / N_CELLS:.2f} per cell")
    elif 'population' in aj.columns:
        print(f"\n    Per-population junction counts:")
        for grp in aj['population'].unique():
            n = len(aj[aj['population'] == grp])
            print(f"      {grp}: {n} junctions, {n / N_CELLS:.2f} per cell")


# =============================================================================
# 2. PANEL B: VENN DIAGRAM GENE SETS
# =============================================================================
print(f"\n{sep}")
print("PANEL B: VENN DIAGRAM - NEOANTIGEN GENE OVERLAP")
print(sep)

sbs2_genes = set()
cnv_genes = set()

for group, gene_set in [('SBS2_HIGH', sbs2_genes), ('CNV_HIGH', cnv_genes)]:
    neo_file = os.path.join(MHC_DIR, f"{group}_neoantigens.tsv")
    if os.path.exists(neo_file):
        neo = pd.read_csv(neo_file, sep='\t')
        if 'gene' in neo.columns:
            gene_set.update(neo['gene'].dropna().unique())

shared = sbs2_genes & cnv_genes
sbs2_only = sbs2_genes - cnv_genes
cnv_only = cnv_genes - sbs2_genes

print(f"\n  SBS2-HIGH neoantigen genes: {len(sbs2_genes)}")
print(f"  CNV-HIGH neoantigen genes:  {len(cnv_genes)}")
print(f"  Shared genes:               {len(shared)}")
print(f"  SBS2-only genes:            {len(sbs2_only)}")
print(f"  CNV-only genes:             {len(cnv_only)}")

# Cross-reference with tiers
if os.path.exists(os.path.join(TIER_DIR, "all_genes_tiered.tsv")):
    tiers = pd.read_csv(os.path.join(TIER_DIR, "all_genes_tiered.tsv"), sep='\t')
    print(f"\n  Tier file columns: {list(tiers.columns)}")
    print(f"\n  Tier breakdown:")
    for tier_name in sorted(tiers['tier'].unique()):
        n = len(tiers[tiers['tier'] == tier_name])
        print(f"    {tier_name}: {n} genes")
    
    # Check: do tier assignments match our Venn counts?
    print(f"\n  Consistency check (Venn vs Tiers):")
    t1a = set(tiers[tiers['tier'] == '1A_hot_shared_escaped']['gene']) if 'gene' in tiers.columns else set()
    t1b = set(tiers[tiers['tier'] == '1B_hot_sbs2_specific']['gene']) if 'gene' in tiers.columns else set()
    t2 = set(tiers[tiers['tier'] == '2_cold_cnv_specific']['gene']) if 'gene' in tiers.columns else set()
    t3 = set(tiers[tiers['tier'] == '3_broad_coverage']['gene']) if 'gene' in tiers.columns else set()
    
    # gene column might be named differently
    if 'gene' not in tiers.columns:
        print(f"    NOTE: 'gene' column not found. Available: {list(tiers.columns)}")
        # try gene_symbol
        gcol = 'gene_symbol' if 'gene_symbol' in tiers.columns else tiers.columns[0]
        print(f"    Using column: {gcol}")
        t1a = set(tiers[tiers['tier'] == '1A_hot_shared_escaped'][gcol])
        t1b = set(tiers[tiers['tier'] == '1B_hot_sbs2_specific'][gcol])
        t2 = set(tiers[tiers['tier'] == '2_cold_cnv_specific'][gcol])
        t3 = set(tiers[tiers['tier'] == '3_broad_coverage'][gcol])
    
    print(f"    Tier 1A (shared+escaped):  {len(t1a)}  — in shared: {len(t1a & shared)}, in sbs2_only: {len(t1a & sbs2_only)}, in cnv_only: {len(t1a & cnv_only)}")
    print(f"    Tier 1B (SBS2-specific):   {len(t1b)}  — in shared: {len(t1b & shared)}, in sbs2_only: {len(t1b & sbs2_only)}, in cnv_only: {len(t1b & cnv_only)}")
    print(f"    Tier 2  (CNV-specific):    {len(t2)}   — in shared: {len(t2 & shared)}, in sbs2_only: {len(t2 & sbs2_only)}, in cnv_only: {len(t2 & cnv_only)}")
    print(f"    Tier 3  (broad, no escape):{len(t3)}   — in shared: {len(t3 & shared)}, in sbs2_only: {len(t3 & sbs2_only)}, in cnv_only: {len(t3 & cnv_only)}")

print(f"\n  Top 20 shared genes (alphabetical): {sorted(list(shared))[:20]}")
print(f"  Top 20 SBS2-only genes: {sorted(list(sbs2_only))[:20]}")
print(f"  Top 20 CNV-only genes: {sorted(list(cnv_only))[:20]}")


# =============================================================================
# 3. PANEL C: ANXA1 TOP PEPTIDES
# =============================================================================
print(f"\n{sep}")
print("PANEL C: ANXA1 TOP PEPTIDES (WT vs MUTANT IC50)")
print(sep)

all_pep_file = os.path.join(MHC_DIR, "SBS2_HIGH_all_peptide_results.tsv")
if os.path.exists(all_pep_file):
    all_pep = pd.read_csv(all_pep_file, sep='\t')
    anxa1_pep = all_pep[all_pep['gene'] == 'ANXA1'].copy()
    
    print(f"\n  ANXA1 peptide columns: {list(anxa1_pep.columns)}")
    print(f"  Total ANXA1 peptide predictions: {len(anxa1_pep)}")
    
    # Filter to binders
    if 'is_binder' in anxa1_pep.columns:
        binders = anxa1_pep[anxa1_pep['is_binder'] == True].copy()
    elif 'mut_ic50' in anxa1_pep.columns:
        binders = anxa1_pep[anxa1_pep['mut_ic50'] < 500].copy()
    else:
        binders = anxa1_pep.copy()
    
    print(f"  ANXA1 binders (IC50 < 500nM): {len(binders)}")
    
    if len(binders) > 0 and 'mut_ic50' in binders.columns:
        binders = binders.sort_values('mut_ic50')
        
        # Show top 20 unique peptides
        if 'mut_peptide' in binders.columns:
            uniq = binders.drop_duplicates(subset=['mut_peptide'], keep='first')
        else:
            uniq = binders
        
        print(f"  Unique binder peptides: {len(uniq)}")
        print(f"\n  Top 20 ANXA1 peptides by mutant IC50:")
        
        show_cols = [c for c in ['mut_peptide', 'wt_peptide', 'mut_ic50', 'wt_ic50', 
                                  'best_allele', 'mut_pos_protein', 'hgvs_p',
                                  'is_strong_binder', 'is_differential'] if c in uniq.columns]
        print(uniq.head(20)[show_cols].to_string(index=False))
        
        # Identify the two strongest peptide REGIONS (by protein position)
        if 'mut_pos_protein' in uniq.columns:
            print(f"\n  Strong binders (IC50 < 50nM) by protein position:")
            strong = uniq[uniq['mut_ic50'] < 50]
            if len(strong) > 0:
                for _, row in strong.iterrows():
                    pos = row.get('mut_pos_protein', '?')
                    pep = row.get('mut_peptide', '?')
                    ic50 = row.get('mut_ic50', '?')
                    wt = row.get('wt_ic50', '?')
                    allele = row.get('best_allele', '?')
                    hgvs = row.get('hgvs_p', '?')
                    print(f"      pos={pos}  IC50={ic50:.1f}nM  WT={wt:.1f}nM  {hgvs}  {pep}  {allele}")
            else:
                print("      No strong binders found")
            
            # Group by position to find hotspot regions
            print(f"\n  Binder count by protein position (top 10 most targeted positions):")
            pos_counts = uniq.groupby('mut_pos_protein').agg(
                n_peptides=('mut_peptide', 'count'),
                best_ic50=('mut_ic50', 'min'),
                hgvs=('hgvs_p', 'first')
            ).sort_values('n_peptides', ascending=False)
            print(pos_counts.head(10).to_string())
        
        # IC50 shift analysis
        if 'wt_ic50' in uniq.columns:
            uniq_top = uniq.head(20).copy()
            uniq_top['ic50_shift'] = uniq_top['wt_ic50'] - uniq_top['mut_ic50']
            uniq_top['fold_change'] = uniq_top['wt_ic50'] / uniq_top['mut_ic50']
            print(f"\n  IC50 shift for top 20 (WT - Mut, positive = mutation improves binding):")
            shift_cols = [c for c in ['mut_peptide', 'mut_ic50', 'wt_ic50', 'ic50_shift', 
                                       'fold_change', 'mut_pos_protein', 'hgvs_p'] if c in uniq_top.columns]
            print(uniq_top[shift_cols].to_string(index=False))


# =============================================================================
# 4. PANEL D: ANXA1 MUTATION POSITIONS FOR GENE TRACK
# =============================================================================
print(f"\n{sep}")
print("PANEL D: ANXA1 MUTATION POSITIONS + DOMAINS")
print(sep)

import re

spa_file = os.path.join(ANNOTATION_DIR, "SBS2_HIGH.somatic_protein_altering.tsv")
if os.path.exists(spa_file):
    spa = pd.read_csv(spa_file, sep='\t')
    anxa1_var = spa[spa['gene'] == 'ANXA1'].copy()
    
    print(f"\n  ANXA1 somatic protein-altering variants: {len(anxa1_var)}")
    print(f"  Columns: {list(anxa1_var.columns)}")
    
    # Parse positions from hgvs_p
    positions = []
    for _, row in anxa1_var.iterrows():
        hgvs = str(row.get('hgvs_p', ''))
        match = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', hgvs)
        if match:
            pos = int(match.group(2))
            ref = row.get('ref', '')
            alt = row.get('alt', '')
            is_apobec = (ref == 'C' and alt in ['T', 'G']) or (ref == 'G' and alt in ['A', 'C'])
            positions.append({
                'position': pos,
                'hgvs_p': hgvs,
                'ref': ref,
                'alt': alt,
                'is_apobec_type': is_apobec,
            })
    
    pos_df = pd.DataFrame(positions)
    if len(pos_df) > 0:
        print(f"\n  Parsed positions: {len(pos_df)}")
        print(f"  Unique positions: {pos_df['position'].nunique()}")
        print(f"  APOBEC-type (C>T/C>G or rev comp): {pos_df['is_apobec_type'].sum()}")
        print(f"  Non-APOBEC: {(~pos_df['is_apobec_type']).sum()}")
        
        # Domain mapping
        ANXA1_DOMAINS = [
            (1, 41, 'N-term'),
            (42, 109, 'Repeat 1'),
            (110, 122, 'Linker 1-2'),
            (123, 196, 'Repeat 2'),
            (197, 220, 'Linker 2-3'),
            (221, 287, 'Repeat 3'),
            (288, 290, 'Linker 3-4'),
            (291, 346, 'Repeat 4'),
        ]
        
        def get_domain(pos):
            for start, end, name in ANXA1_DOMAINS:
                if start <= pos <= end:
                    return name
            return 'Unknown'
        
        pos_df['domain'] = pos_df['position'].apply(get_domain)
        
        print(f"\n  Mutations by domain:")
        for dom in pos_df['domain'].unique():
            sub = pos_df[pos_df['domain'] == dom]
            n_apobec = sub['is_apobec_type'].sum()
            print(f"    {dom}: {len(sub)} mutations ({n_apobec} APOBEC-type)")
        
        print(f"\n  All ANXA1 mutations (sorted by position):")
        print(pos_df.sort_values('position').to_string(index=False))
        
        # Cross-reference with binding data
        if os.path.exists(all_pep_file):
            all_pep = pd.read_csv(all_pep_file, sep='\t')
            anxa1_pep = all_pep[all_pep['gene'] == 'ANXA1'].copy()
            
            print(f"\n  Mutation positions with neoantigen binding data:")
            for pos in sorted(pos_df['position'].unique()):
                if 'mut_pos_protein' in anxa1_pep.columns:
                    pos_pep = anxa1_pep[anxa1_pep['mut_pos_protein'] == pos]
                else:
                    pos_pep = pd.DataFrame()
                
                n_pep = len(pos_pep)
                best_ic50 = pos_pep['mut_ic50'].min() if n_pep > 0 else float('inf')
                n_strong = (pos_pep['mut_ic50'] < 50).sum() if n_pep > 0 else 0
                hgvs = pos_df[pos_df['position'] == pos]['hgvs_p'].iloc[0]
                apobec = pos_df[pos_df['position'] == pos]['is_apobec_type'].iloc[0]
                domain = pos_df[pos_df['position'] == pos]['domain'].iloc[0]
                
                flag = " *** STRONG" if best_ic50 < 50 else ""
                print(f"    pos={pos:3d}  {domain:10s}  {hgvs:20s}  APOBEC={apobec}  "
                      f"peptides={n_pep:2d}  best_IC50={best_ic50:7.1f}  strong={n_strong}{flag}")

else:
    print(f"  FILE NOT FOUND: {spa_file}")

print(f"\n{sep}")
print("DIAGNOSTIC COMPLETE")
print(sep)
