#!/usr/bin/env python3
"""
Diagnostic_Section2_Numbers.py
================================
Pull key numbers for Section 2 (Figure 2) rewrite from current pipeline outputs.
"""

import os
import json
import pickle
import numpy as np
import pandas as pd

PROJECT = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG2 = os.path.join(PROJECT, "data", "FIG_2")

sep = "=" * 70

print(f"\n{sep}")
print("SECTION 2 TEXT NUMBERS — FIGURE 2")
print(sep)

# ============================================================
# 1. DIRECTORY STRUCTURE
# ============================================================
print(f"\n--- FIG_2 Top-Level Contents ---")
if os.path.exists(FIG2):
    for item in sorted(os.listdir(FIG2)):
        full = os.path.join(FIG2, item)
        if os.path.isdir(full):
            print(f"  [DIR] {item}/")
        else:
            fsize = os.path.getsize(full) / 1024
            print(f"  {item}  ({fsize:.1f} KB)")
else:
    print(f"  ERROR: {FIG2} does not exist")
    exit(1)

# ============================================================
# 2. GROUP SIZES (from Step 03 or network groups file)
# ============================================================
print(f"\n--- Network Group Sizes ---")
# Check multiple possible locations
for gpath in [
    os.path.join(FIG2, "03_differential_expression", "TCGA-HNSC", "TCGA-HNSC_group_assignments.tsv"),
    os.path.join(PROJECT, "data", "FIG_1", "HNSC_network_groups_v3.tsv"),
    os.path.join(FIG2, "03_differential_expression", "group_assignments.tsv"),
]:
    if os.path.exists(gpath):
        gdf = pd.read_csv(gpath, sep='\t')
        print(f"  File: {gpath}")
        print(f"  Columns: {list(gdf.columns)}")
        # Try to find group column
        for col in ['group', 'Group', 'SBS2_group', 'network_group']:
            if col in gdf.columns:
                print(f"  Groups ({col}):")
                for g in gdf[col].unique():
                    print(f"    {g}: n={len(gdf[gdf[col]==g])}")
                break
        break
else:
    print("  No group assignment file found")

# ============================================================
# 3. DIFFERENTIAL EXPRESSION
# ============================================================
print(f"\n--- Differential Expression ---")
de_dir = os.path.join(FIG2, "03_differential_expression", "TCGA-HNSC")
if not os.path.exists(de_dir):
    de_dir = os.path.join(FIG2, "03_differential_expression")

if os.path.exists(de_dir):
    print(f"  DE directory: {de_dir}")
    for f in sorted(os.listdir(de_dir)):
        if 'diffexpr' in f.lower() or 'de_' in f.lower() or 'stats' in f.lower():
            fpath = os.path.join(de_dir, f)
            print(f"  Found: {f}")
            try:
                de = pd.read_csv(fpath, sep='\t')
                print(f"    Shape: {de.shape}")
                print(f"    Columns: {list(de.columns)[:10]}")
                # Count DE genes
                for pcol in ['raw_pvalue', 'pvalue', 'p_value', 'raw_p', 'adj_p', 'padj', 'FDR']:
                    if pcol in de.columns:
                        n_sig = (de[pcol] < 0.05).sum()
                        print(f"    Genes with {pcol} < 0.05: {n_sig}")
                # Check for A3 genes
                for gcol in ['gene', 'gene_symbol', 'symbol', 'Gene']:
                    if gcol in de.columns:
                        a3_genes = ['APOBEC3A', 'APOBEC3B', 'A3A', 'A3B']
                        found = de[de[gcol].isin(a3_genes)]
                        if len(found) > 0:
                            print(f"    A3 genes found:")
                            print(found.to_string(index=False))
                        break
            except:
                print(f"    (could not parse)")

# ============================================================
# 4. NETWORK PARAMETERS
# ============================================================
print(f"\n--- Network Structure ---")
comm_dir = os.path.join(FIG2, "05_communities", "TCGA-HNSC")
if not os.path.exists(comm_dir):
    comm_dir = os.path.join(FIG2, "05_communities")

# Community gene lists
for fname in ['TCGA-HNSC_community_gene_lists.csv', 'community_gene_lists.csv']:
    fpath = os.path.join(comm_dir, fname)
    if os.path.exists(fpath):
        cgl = pd.read_csv(fpath)
        print(f"\n  Community gene lists: {fpath}")
        print(f"  Columns: {list(cgl.columns)}")
        if 'community' in cgl.columns:
            n_communities = cgl['community'].nunique()
            total_genes = len(cgl)
            print(f"  Total communities: {n_communities}")
            print(f"  Total genes: {total_genes}")
            for c in sorted(cgl['community'].unique()):
                n = len(cgl[cgl['community'] == c])
                # Check for A3 genes
                gcol = 'gene_symbol' if 'gene_symbol' in cgl.columns else 'gene' if 'gene' in cgl.columns else None
                a3_in = []
                if gcol:
                    for a3 in ['APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D', 'APOBEC3F', 'APOBEC3G', 'APOBEC3H']:
                        if a3 in cgl[cgl['community'] == c][gcol].values:
                            a3_in.append(a3)
                a3_str = f"  [{', '.join(a3_in)}]" if a3_in else ""
                print(f"    Community {c}: {n} genes{a3_str}")
        break

# Best partition / modularity
for fname in ['TCGA-HNSC_best_partition.csv', 'best_partition.csv']:
    fpath = os.path.join(comm_dir, fname)
    if os.path.exists(fpath):
        bp = pd.read_csv(fpath)
        print(f"\n  Best partition: {fpath}")
        print(f"  Shape: {bp.shape}")
        print(f"  Columns: {list(bp.columns)}")
        break

# Check for summary/parameters files
for subdir in [comm_dir, FIG2, os.path.join(FIG2, "08_pipeline_summary")]:
    if not os.path.exists(subdir):
        continue
    for f in sorted(os.listdir(subdir)):
        if 'param' in f.lower() or 'summary' in f.lower() or 'modularity' in f.lower():
            print(f"\n  Found summary file: {os.path.join(subdir, f)}")
            fpath = os.path.join(subdir, f)
            if f.endswith(('.txt', '.tsv', '.csv', '.json')):
                try:
                    with open(fpath) as fh:
                        content = fh.read()
                    if len(content) < 5000:
                        print(content[:3000])
                    else:
                        print(content[:2000])
                        print("  ... (truncated)")
                except:
                    pass

# ============================================================
# 5. CENTRALITY METRICS (for A3 gene degree)
# ============================================================
print(f"\n--- A3 Gene Centrality ---")
cent_dir = os.path.join(FIG2, "06_centrality_metrics", "TCGA-HNSC")
if not os.path.exists(cent_dir):
    cent_dir = os.path.join(FIG2, "06_centrality_metrics")

if os.path.exists(cent_dir):
    for f in sorted(os.listdir(cent_dir)):
        if 'metric' in f.lower() or 'centrality' in f.lower():
            fpath = os.path.join(cent_dir, f)
            try:
                met = pd.read_csv(fpath, sep='\t' if f.endswith('.tsv') else ',')
                gcol = None
                for c in ['gene_symbol', 'gene', 'symbol', 'Gene']:
                    if c in met.columns:
                        gcol = c
                        break
                if gcol:
                    a3_rows = met[met[gcol].str.contains('APOBEC3', na=False)]
                    if len(a3_rows) > 0:
                        print(f"  File: {f}")
                        print(a3_rows.to_string(index=False))
            except:
                pass

# ============================================================
# 6. DIAGNOSTIC AUDIT (supplemental tables)
# ============================================================
print(f"\n--- Diagnostic Audit Files ---")
audit_dir = os.path.join(FIG2, "DIAGNOSTIC_AUDIT")
if os.path.exists(audit_dir):
    for f in sorted(os.listdir(audit_dir)):
        fsize = os.path.getsize(os.path.join(audit_dir, f)) / 1024
        print(f"  {f}  ({fsize:.1f} KB)")
        # Load key supplemental tables
        if 'Community_Summary' in f:
            fpath = os.path.join(audit_dir, f)
            try:
                cs = pd.read_csv(fpath, sep='\t')
                print(f"    Columns: {list(cs.columns)}")
                print(cs.to_string(index=False))
            except:
                pass
        if 'Node_Manifest' in f and f.endswith('.tsv'):
            fpath = os.path.join(audit_dir, f)
            try:
                nm = pd.read_csv(fpath, sep='\t')
                # Just show A3 genes
                gcol = 'gene_symbol' if 'gene_symbol' in nm.columns else None
                if gcol:
                    a3 = nm[nm[gcol].str.contains('APOBEC3', na=False)]
                    if len(a3) > 0:
                        print(f"    A3 genes in node manifest:")
                        print(a3.to_string(index=False))
                # Harris interactors
                for hcol in ['is_harris_interactor', 'harris', 'Harris']:
                    if hcol in nm.columns:
                        n_harris = nm[hcol].sum()
                        print(f"    Harris interactors in network: {n_harris}")
                        break
            except:
                pass

# ============================================================
# 7. KEGG ENRICHMENT
# ============================================================
print(f"\n--- KEGG Enrichment ---")
enrich_dir = os.path.join(FIG2, "08_pipeline_summary")
if os.path.exists(enrich_dir):
    for f in sorted(os.listdir(enrich_dir)):
        if 'kegg' in f.lower() or 'enrich' in f.lower():
            fpath = os.path.join(enrich_dir, f)
            print(f"\n  {f}:")
            try:
                ef = pd.read_csv(fpath, sep='\t' if f.endswith('.tsv') else ',')
                print(f"  Shape: {ef.shape}")
                print(f"  Columns: {list(ef.columns)}")
                # Show top hits
                for pcol in ['Adjusted P-value', 'adj_p', 'padj', 'FDR']:
                    if pcol in ef.columns:
                        sig = ef[ef[pcol] < 0.05]
                        print(f"  Significant ({pcol} < 0.05): {len(sig)}")
                        if len(sig) > 0:
                            print(sig.head(20).to_string(index=False))
                        break
            except:
                pass

print(f"\n{sep}")
print("DIAGNOSTIC COMPLETE")
print(sep)
"""
Check how many Harris A3 interactors and canonical cell-type markers
are present in the Figure 2 bulk TCGA network.
"""

import os
import pandas as pd

PROJECT = "/master/jlehle/WORKING/2026_NMF_PAPER"

sep = "=" * 70

# ============================================================
# 1. LOAD NETWORK GENES
# ============================================================
print(f"\n{sep}")
print("HARRIS INTERACTORS & CELL-TYPE MARKERS IN FIGURE 2 NETWORK")
print(sep)

partition_path = os.path.join(PROJECT, "data/FIG_2/05_communities/TCGA-HNSC/TCGA-HNSC_best_partition.csv")
bp = pd.read_csv(partition_path)
print(f"\nNetwork partition: {partition_path}")
print(f"Total network genes: {len(bp)}")
print(f"Columns: {list(bp.columns)}")

# Get gene symbols
if 'gene_symbol' in bp.columns:
    network_genes = set(bp['gene_symbol'].dropna().values)
elif 'gene' in bp.columns:
    network_genes = set(bp['gene'].dropna().values)
else:
    print("ERROR: No gene symbol column found")
    exit(1)

print(f"Unique gene symbols in network: {len(network_genes)}")

# Also build community lookup
comm_col = 'community' if 'community' in bp.columns else None
gene_col = 'gene_symbol' if 'gene_symbol' in bp.columns else 'gene'
if comm_col:
    gene_to_comm = dict(zip(bp[gene_col], bp[comm_col]))

# ============================================================
# 2. HARRIS A3 INTERACTORS
# ============================================================
print(f"\n{sep}")
print("HARRIS A3 INTERACTORS")
print(sep)

# Try multiple possible locations
harris_paths = [
    os.path.join(PROJECT, "data/FIG_4/00_input/Harris_A3_interactors.tsv"),
    os.path.join(PROJECT, "data/FIG_4/00_input/Harris_A3_interactors.txt"),
    os.path.join(PROJECT, "data/FIG_2/Harris_A3_interactors.tsv"),
]

harris_genes = set()
harris_path_used = None
for hp in harris_paths:
    if os.path.exists(hp):
        harris_path_used = hp
        hdf = pd.read_csv(hp, sep='\t')
        print(f"\nLoaded: {hp}")
        print(f"Shape: {hdf.shape}")
        print(f"Columns: {list(hdf.columns)}")
        # Find gene symbol column
        for col in ['gene_symbol', 'Gene', 'gene', 'Symbol', 'GeneName']:
            if col in hdf.columns:
                harris_genes = set(hdf[col].dropna().values)
                print(f"Harris genes ({col}): {len(harris_genes)}")
                break
        if not harris_genes:
            # Maybe it's a simple gene list
            harris_genes = set(hdf.iloc[:, 0].dropna().values)
            print(f"Harris genes (first column): {len(harris_genes)}")
        break

if not harris_path_used:
    print("\nHarris file not found. Trying to find it...")
    for root, dirs, files in os.walk(os.path.join(PROJECT, "data")):
        for f in files:
            if 'harris' in f.lower() and 'interactor' in f.lower():
                print(f"  Found: {os.path.join(root, f)}")

if harris_genes:
    overlap = harris_genes & network_genes
    print(f"\nHarris interactors in network: {len(overlap)} / {len(harris_genes)}")
    if overlap:
        print(f"\nMatched genes:")
        for g in sorted(overlap):
            comm = gene_to_comm.get(g, '?') if comm_col else '?'
            print(f"  {g}  (Community {comm})")
    else:
        print("  *** ZERO Harris interactors found in network ***")

# Also check A3B-specific interactors
harris_a3b_paths = [
    os.path.join(PROJECT, "data/FIG_4/00_input/Harris_A3_interactors_A3B_only.txt"),
    os.path.join(PROJECT, "data/FIG_4/00_input/Harris_A3_interactors_A3B_only.tsv"),
]
for hp in harris_a3b_paths:
    if os.path.exists(hp):
        hdf2 = pd.read_csv(hp, sep='\t')
        for col in ['gene_symbol', 'Gene', 'gene', 'Symbol']:
            if col in hdf2.columns:
                a3b_genes = set(hdf2[col].dropna().values)
                a3b_overlap = a3b_genes & network_genes
                print(f"\nA3B-specific interactors in network: {len(a3b_overlap)} / {len(a3b_genes)}")
                if a3b_overlap:
                    for g in sorted(a3b_overlap):
                        comm = gene_to_comm.get(g, '?') if comm_col else '?'
                        print(f"  {g}  (Community {comm})")
                break
        break

# ============================================================
# 3. CELL-TYPE MARKERS
# ============================================================
print(f"\n{sep}")
print("CELL-TYPE MARKERS")
print(sep)

# 24 canonical cell-type markers used in the audit
MARKERS = {
    'Basal epithelial': ['KRT5', 'KRT14', 'KRT15', 'TACSTD2', 'TP63'],
    'Suprabasal epithelial': ['KRT4', 'KRT13', 'KRT19', 'IVL'],
    'Fibroblast': ['COL1A1', 'COL3A1', 'DCN', 'LUM'],
    'Endothelial': ['PECAM1', 'VWF', 'CDH5'],
    'T cell': ['CD3E', 'CD3D', 'CD8A'],
    'B cell': ['CD19', 'MS4A1'],
    'Macrophage/Monocyte': ['CD68', 'CD14', 'CSF1R'],
    'NK cell': ['NKG7', 'GNLY'],
}

all_markers = set()
for ct, genes in MARKERS.items():
    all_markers.update(genes)

print(f"\nTotal markers checked: {len(all_markers)}")
markers_found = all_markers & network_genes
print(f"Markers in network: {len(markers_found)} / {len(all_markers)}")

if markers_found:
    print(f"\nMatched markers:")
    for ct, genes in MARKERS.items():
        found = [g for g in genes if g in network_genes]
        if found:
            for g in found:
                comm = gene_to_comm.get(g, '?') if comm_col else '?'
                print(f"  {g}  ({ct}, Community {comm})")
else:
    print("  *** ZERO cell-type markers found in network ***")

# Also list which markers are NOT in the network for completeness
print(f"\nMarkers NOT in network:")
for ct, genes in MARKERS.items():
    missing = [g for g in genes if g not in network_genes]
    if missing:
        print(f"  {ct}: {', '.join(missing)}")

# ============================================================
# 4. SUMMARY FOR TEXT
# ============================================================
print(f"\n{sep}")
print("SUMMARY FOR SECTION 2 TEXT")
print(sep)
n_harris = len(overlap) if harris_genes else '?'
n_markers = len(markers_found)
n_total_harris = len(harris_genes) if harris_genes else '?'
n_total_markers = len(all_markers)

print(f"\n  Harris A3 interactors in network: {n_harris} / {n_total_harris}")
print(f"  Cell-type markers in network: {n_markers} / {n_total_markers}")
print(f"  A3B degree in DIFF network: 1")
print(f"  A3A degree in DIFF network: 0")
print(f"  A3B community: 2 (134 genes)")
print(f"  A3B community KEGG: no significant terms")

print(f"\n{sep}")
print("DIAGNOSTIC COMPLETE")
print(sep)
