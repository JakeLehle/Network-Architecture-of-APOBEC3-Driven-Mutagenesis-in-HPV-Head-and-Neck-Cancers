#!/usr/bin/env python3
"""
Diagnostic_Section4_Numbers.py
================================
Pull every number needed for the ~550-word Section 4 rewrite.
Covers: three pairwise network stats, Harris interactor recovery,
activating/inhibiting chains, A3 wall, directional coherence.
"""

import os
import sys
import numpy as np
import pandas as pd

PROJECT = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG_4 = os.path.join(PROJECT, "data", "FIG_4")

sep = "=" * 70

print(f"\n{sep}")
print("SECTION 4 TEXT NUMBERS — FIGURE 4")
print(sep)

# ============================================================
# 1. DIRECTORY STRUCTURE
# ============================================================
print(f"\n--- FIG_4 Directory Structure ---")
if os.path.exists(FIG_4):
    for item in sorted(os.listdir(FIG_4)):
        full = os.path.join(FIG_4, item)
        if os.path.isdir(full):
            print(f"  [DIR] {item}/")
            # Show subdirectories one level deep
            for sub in sorted(os.listdir(full)):
                subfull = os.path.join(full, sub)
                if os.path.isdir(subfull):
                    print(f"    [DIR] {sub}/")
                elif os.path.getsize(subfull) > 1024:
                    print(f"    {sub}  ({os.path.getsize(subfull)/1024:.1f} KB)")

# ============================================================
# 2. LOAD HARRIS INTERACTORS
# ============================================================
print(f"\n{sep}")
print("HARRIS A3 INTERACTORS")
print(sep)

harris_path = os.path.join(FIG_4, "00_input", "Harris_A3_interactors.txt")
if os.path.exists(harris_path):
    harris_df = pd.read_csv(harris_path, sep='\t')
    harris_genes = set(harris_df['gene_symbol'].dropna().values)
    print(f"  Total Harris interactors: {len(harris_genes)}")
    a3b_specific = set(harris_df[harris_df['A3B_interactor'] == True]['gene_symbol'].values) \
        if 'A3B_interactor' in harris_df.columns else set()
    print(f"  A3B-specific: {len(a3b_specific)}")
else:
    print(f"  Harris file not found at {harris_path}")
    harris_genes = set()

# ============================================================
# 3. FIND AND ANALYZE EACH NETWORK
# ============================================================
NETWORK_NAMES = [
    ("SBS2_VS_CNV", "SBS2-HIGH vs CNV-HIGH"),
    ("SBS2_VS_NORMAL", "SBS2-HIGH vs NORMAL"),
    ("CNV_VS_NORMAL", "CNV-HIGH vs NORMAL"),
    ("NETWORK_A", "Network A"),
    ("NETWORK_B", "Network B"),
    ("NETWORK_C", "Network C"),
    ("NETWORK_SBS2_VS_NORMAL", "SBS2 vs NORMAL"),
    ("NETWORK_SBS2_VS_CNV", "SBS2 vs CNV"),
    ("NETWORK_CNV_VS_NORMAL", "CNV vs NORMAL"),
]

def find_network_dirs():
    """Search for network output directories across FIG_4."""
    found = {}
    for root, dirs, files in os.walk(FIG_4):
        for d in dirs:
            for name, label in NETWORK_NAMES:
                if name.lower() == d.lower():
                    full = os.path.join(root, d)
                    # Check if this looks like a network output dir
                    has_partition = any('partition' in f.lower() or 'community' in f.lower()
                                       for f in os.listdir(full) if os.path.isfile(os.path.join(full, f)))
                    has_subdir = any(os.path.isdir(os.path.join(full, s))
                                    for s in os.listdir(full)
                                    if 'communit' in s.lower() or 'centrality' in s.lower() or 'pipeline' in s.lower())
                    if has_partition or has_subdir:
                        found[name] = full
                        break
        # Also check for pipeline step directories
        for f in files:
            if 'selected_parameters' in f.lower():
                # This directory contains network output
                parent_name = os.path.basename(root)
                for name, label in NETWORK_NAMES:
                    if name.lower() in parent_name.lower() or name.lower() in root.lower():
                        found[name] = root
                        break
    return found

net_dirs = find_network_dirs()
print(f"\n{sep}")
print("NETWORK DIRECTORIES FOUND")
print(sep)
for name, path in sorted(net_dirs.items()):
    print(f"  {name}: {path}")

# Also search more broadly
print(f"\n--- Searching for parameter files ---")
for root, dirs, files in os.walk(FIG_4):
    for f in files:
        if 'selected_parameters' in f.lower() or 'pipeline_summary' in f.lower():
            print(f"  {os.path.join(root, f)}")

print(f"\n--- Searching for partition files ---")
for root, dirs, files in os.walk(FIG_4):
    for f in files:
        if 'best_partition' in f.lower() or 'community_gene_list' in f.lower():
            print(f"  {os.path.join(root, f)}")

# ============================================================
# 4. ANALYZE EACH NETWORK
# ============================================================

def analyze_network(net_dir, name, harris_genes):
    """Extract key metrics from a network output directory."""
    print(f"\n{'='*60}")
    print(f"  NETWORK: {name}")
    print(f"  Directory: {net_dir}")
    print(f"{'='*60}")

    # Find parameter file
    for root, dirs, files in os.walk(net_dir):
        for f in files:
            if 'selected_parameters' in f.lower():
                fpath = os.path.join(root, f)
                print(f"\n  Parameters ({f}):")
                with open(fpath) as fh:
                    for line in fh:
                        print(f"    {line.strip()}")

    # Find partition file
    partition = None
    for root, dirs, files in os.walk(net_dir):
        for f in files:
            if 'best_partition' in f.lower() and f.endswith('.csv'):
                fpath = os.path.join(root, f)
                partition = pd.read_csv(fpath)
                print(f"\n  Partition: {fpath}")
                print(f"    Total genes: {len(partition)}")
                print(f"    Columns: {list(partition.columns)}")

                # Community counts
                if 'community' in partition.columns:
                    n_comm = partition['community'].nunique()
                    sizes = partition.groupby('community').size().sort_values(ascending=False)
                    print(f"    Communities: {n_comm}")
                    print(f"    Largest 5: {dict(sizes.head())}")

                    # Gene symbol column
                    gcol = 'gene_symbol' if 'gene_symbol' in partition.columns else 'gene'

                    # Harris interactors
                    net_genes = set(partition[gcol].dropna().values)
                    harris_in = harris_genes & net_genes
                    print(f"\n    Harris interactors in network: {len(harris_in)} / {len(harris_genes)}")

                    # Harris by community
                    if harris_in:
                        gene_comm = dict(zip(partition[gcol], partition['community']))
                        harris_comm = {}
                        for g in harris_in:
                            c = gene_comm.get(g, '?')
                            harris_comm.setdefault(c, []).append(g)
                        for c in sorted(harris_comm.keys()):
                            genes = harris_comm[c]
                            print(f"      Community {c}: {len(genes)} - {', '.join(sorted(genes)[:10])}")

                    # A3 genes
                    a3_in = [g for g in net_genes if 'APOBEC3' in g]
                    print(f"\n    A3 genes in network: {a3_in}")
                    for a3g in a3_in:
                        c = gene_comm.get(a3g, '?')
                        print(f"      {a3g}: Community {c}")
                break
        if partition is not None:
            break

    # Find centrality/degree for A3 genes
    for root, dirs, files in os.walk(net_dir):
        for f in files:
            if 'diff_metrics' in f.lower() or 'centrality' in f.lower():
                if f.endswith(('.csv', '.tsv')):
                    fpath = os.path.join(root, f)
                    try:
                        sep_char = '\t' if f.endswith('.tsv') else ','
                        met = pd.read_csv(fpath, sep=sep_char)
                        gcol = None
                        for c in ['gene_symbol', 'gene', 'symbol']:
                            if c in met.columns:
                                gcol = c
                                break
                        if gcol and 'degree' in met.columns:
                            a3_rows = met[met[gcol].str.contains('APOBEC3', na=False)]
                            if len(a3_rows) > 0:
                                print(f"\n    A3 degree ({f}):")
                                for _, row in a3_rows.iterrows():
                                    print(f"      {row[gcol]}: degree={row['degree']}")
                    except:
                        pass

    # Find KEGG enrichment
    for root, dirs, files in os.walk(net_dir):
        for f in files:
            if 'kegg' in f.lower() and 'summary' in f.lower():
                fpath = os.path.join(root, f)
                try:
                    kdf = pd.read_csv(fpath, sep=',' if f.endswith('.csv') else '\t')
                    if 'n_kegg_significant' in kdf.columns:
                        sig = kdf[kdf['n_kegg_significant'] > 0]
                        print(f"\n    KEGG significant communities: {len(sig)}")
                        for _, row in sig.iterrows():
                            print(f"      Community {row.get('community','?')} ({row.get('n_genes','?')} genes): "
                                  f"{row.get('top_kegg_term','?')} (p={row.get('top_kegg_pvalue','?'):.2e})")
                except:
                    pass

    return partition

for name, net_dir in sorted(net_dirs.items()):
    analyze_network(net_dir, name, harris_genes)

# ============================================================
# 5. CHAIN GENES AND A3 WALL
# ============================================================
print(f"\n{sep}")
print("CHAIN GENES AND A3 WALL")
print(sep)

# Search for chain-related output files
for root, dirs, files in os.walk(FIG_4):
    for f in files:
        if 'chain' in f.lower() or 'wall' in f.lower() or 'concordan' in f.lower():
            fpath = os.path.join(root, f)
            fsize = os.path.getsize(fpath) / 1024
            print(f"  {os.path.relpath(fpath, FIG_4)}  ({fsize:.1f} KB)")
            if f.endswith(('.tsv', '.csv', '.txt')) and fsize < 50:
                try:
                    with open(fpath) as fh:
                        content = fh.read()
                    if len(content) < 3000:
                        print(content)
                    else:
                        print(content[:2000])
                        print("  ... (truncated)")
                except:
                    pass

# ============================================================
# 6. SESSION HANDOFF FILES
# ============================================================
print(f"\n{sep}")
print("SESSION HANDOFF FILES")
print(sep)

for root, dirs, files in os.walk(PROJECT):
    for f in files:
        if 'handoff' in f.lower() and ('figure4' in f.lower() or 'fig4' in f.lower() or 'concordance' in f.lower()):
            fpath = os.path.join(root, f)
            print(f"\n  {fpath}")
            try:
                with open(fpath) as fh:
                    content = fh.read()
                # Print first 3000 chars
                print(content[:3000])
                if len(content) > 3000:
                    print("  ... (truncated)")
            except:
                pass

print(f"\n{sep}")
print("DIAGNOSTIC COMPLETE")
print(sep)
