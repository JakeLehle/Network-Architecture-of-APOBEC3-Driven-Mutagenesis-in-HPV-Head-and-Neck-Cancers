#!/usr/bin/env python3
"""
Reparse_Chimeric_With_GenePairs.py
=====================================
Re-parse STAR Chimeric.out.junction files to save junction-level gene pair
table (geneA, geneB, barcode, group, chrom_a, chrom_b) that was not saved
by the original parse_chimeric_junctions_filtered.py.

Uses identical logic and filters as the original parser:
  1. Unique mapped (junction_type 1 or 2)
  2. Alignment score ratio >= 0.80
  3. Both breakpoints in annotated genes
  4. Both genes protein-coding
  5. No artifact genes
  6. Inter-chromosomal or same-chrom >1Mb

Then runs enrichment analysis on group-exclusive recurrent gene pairs.

Inputs:
  - STAR Chimeric.out.junction files in star_chimeric/SRR*/
  - GTF annotation (genes_unzipped.gtf)
  - Barcode lists per group (barcodes_*.tsv)

Outputs (to fusion_detection/gsea_diagnostic/):
  - all_filtered_junctions.tsv (full junction-level table)
  - *_exclusive_pairs.tsv, *_exclusive_genes.txt
  - *_enrichment_*.tsv
  - reparse_chimeric_report.txt

Env: NETWORK (needs gseapy; does NOT need STAR/pysam)
Usage: conda run -n NETWORK python Reparse_Chimeric_With_GenePairs.py
Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"

STAR_DIR = os.path.join(PROJECT_ROOT,
    "data/FIG_6/05_neoantigen/fusion_detection/star_chimeric")
INPUT_DIR = os.path.join(PROJECT_ROOT,
    "data/FIG_6/05_neoantigen/inputs")
FUSION_DIR = os.path.join(PROJECT_ROOT,
    "data/FIG_6/05_neoantigen/fusion_detection")
GTF_PATH = "/master/jlehle/WORKING/SC/ref/GRCh38/genes/genes_unzipped.gtf"

OUTPUT_DIR = os.path.join(FUSION_DIR, "gsea_diagnostic")
os.makedirs(OUTPUT_DIR, exist_ok=True)

GROUPS = ['SBS2_HIGH', 'Stealth_CNV', 'Normal_Control']
MIN_ALN_SCORE_RATIO = 0.80
MIN_DISTANCE_SAME_CHROM = 1_000_000
MIN_CELLS_RECURRENT = 2

# Artifact gene detection (matching original parser exactly)
ARTIFACT_PREFIXES = [
    'RPL', 'RPS', 'MRPL', 'MRPS', 'MT-',
    'IGH', 'IGK', 'IGL',
    'TRA', 'TRB', 'TRD', 'TRG',
    'LINC', 'LOC', 'SNOR', 'SNAR', 'MIR', 'ENSG',
]
ARTIFACT_GENES = {'MALAT1', 'NEAT1', 'XIST', 'ACTB', 'GAPDH'}

def is_artifact_gene(g):
    if g in ('intergenic', 'unknown') or g in ARTIFACT_GENES:
        return True
    for p in ARTIFACT_PREFIXES:
        if g.startswith(p):
            return g not in ('TRA2A', 'TRA2B')
    return False

# =============================================================================
# LOGGING
# =============================================================================
report_lines = []

def log(msg=""):
    print(msg, flush=True)
    report_lines.append(str(msg))

def banner(title, char="="):
    log(""); log(char * 80); log(f"  {title}"); log(char * 80)

# =============================================================================
# STEP 1: LOAD BARCODE -> GROUP MAPPING
# =============================================================================
banner("STEP 1: Load barcode -> group mapping")

bc_to_group = {}
for group in GROUPS:
    bc_path = os.path.join(INPUT_DIR, f"barcodes_{group}.tsv")
    if os.path.exists(bc_path):
        barcodes = open(bc_path).read().strip().split("\n")
        for bc in barcodes:
            bc_to_group[bc] = group
        log(f"  {group}: {len(barcodes)} barcodes from {bc_path}")
    else:
        log(f"  WARNING: {bc_path} not found")

log(f"  Total target barcodes: {len(bc_to_group)}")

# Build SRR-specific lookup: srr -> {16bp_raw -> full_barcode}
srr_bc_map = defaultdict(dict)
for bc in bc_to_group:
    parts = bc.split("-")
    if len(parts) >= 3:
        raw_16 = parts[0]
        srr = parts[2]
        srr_bc_map[srr][raw_16] = bc

log(f"  SRR IDs with target barcodes: {len(srr_bc_map)}")

# =============================================================================
# STEP 2: LOAD GTF GENE ANNOTATIONS
# =============================================================================
banner("STEP 2: Load GTF gene annotations")

gene_intervals = defaultdict(list)
gene_biotype = {}
n_genes = 0

with open(GTF_PATH) as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 9 or fields[2] != "gene":
            continue
        chrom = fields[0]
        start, end = int(fields[3]), int(fields[4])
        gene_name, biotype = "unknown", "unknown"
        for attr in fields[8].split(";"):
            attr = attr.strip()
            if attr.startswith("gene_name"):
                gene_name = attr.split('"')[1] if '"' in attr else attr.split(" ")[-1]
            elif attr.startswith("gene_biotype") or attr.startswith("gene_type"):
                biotype = attr.split('"')[1] if '"' in attr else attr.split(" ")[-1]
        gene_intervals[chrom].append((start, end, gene_name))
        gene_biotype[gene_name] = biotype
        n_genes += 1

log(f"  Loaded {n_genes} gene annotations across {len(gene_intervals)} chromosomes")
n_pc = sum(1 for v in gene_biotype.values() if v == "protein_coding")
log(f"  Protein-coding genes: {n_pc}")

def pos_to_gene(chrom, pos):
    """Map genomic position to gene name (linear scan, same as original)."""
    for s, e, name in gene_intervals.get(chrom, []):
        if s <= pos <= e:
            return name
    return "intergenic"

# =============================================================================
# STEP 3: PARSE CHIMERIC JUNCTION FILES
# =============================================================================
banner("STEP 3: Parse STAR Chimeric.out.junction files")

junction_files = sorted(glob.glob(os.path.join(STAR_DIR, "SRR*", "Chimeric.out.junction")))
log(f"  Found {len(junction_files)} junction files")

cascade = Counter()
all_junctions = []  # This is the key addition: save every passing junction
sample_summaries = []

for jfile in junction_files:
    srr = os.path.basename(os.path.dirname(jfile))
    target_bcs = srr_bc_map.get(srr, {})
    if not target_bcs:
        continue

    n_total = n_target = n_passed = 0

    with open(jfile) as f:
        for line in f:
            if line.startswith("#") or line.startswith("chr_donorA"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 22:
                continue
            try:
                pos_a = int(fields[1])
                pos_b = int(fields[4])
            except ValueError:
                continue

            n_total += 1
            cascade['total'] += 1

            # Cell barcode from field[20]
            raw_cb = fields[20]
            full_bc = target_bcs.get(raw_cb)
            if full_bc is None:
                continue
            n_target += 1
            cascade['in_target'] += 1

            chrom_a, chrom_b = fields[0], fields[3]

            # FILTER 1: Unique mapped (type 1 or 2)
            try:
                jtype = int(fields[6])
            except ValueError:
                continue
            if jtype == 0:
                continue
            cascade['unique_mapped'] += 1

            # FILTER 2: Alignment score ratio >= 0.80
            try:
                max_s = int(fields[15])
                chim_s = int(fields[17])
                ratio = chim_s / max_s if max_s > 0 else 0
            except (ValueError, IndexError):
                ratio = 0
            if ratio < MIN_ALN_SCORE_RATIO:
                continue
            cascade['good_aln_score'] += 1

            # FILTER 3: Both genic
            gene_a = pos_to_gene(chrom_a, pos_a)
            gene_b = pos_to_gene(chrom_b, pos_b)
            if gene_a == "intergenic" or gene_b == "intergenic":
                continue
            cascade['both_genic'] += 1

            # FILTER 4: Both protein-coding
            if gene_biotype.get(gene_a, "") != "protein_coding" or \
               gene_biotype.get(gene_b, "") != "protein_coding":
                continue
            cascade['both_protein_coding'] += 1

            # FILTER 5: No artifacts
            if is_artifact_gene(gene_a) or is_artifact_gene(gene_b):
                continue
            cascade['no_artifacts'] += 1

            # FILTER 6: Inter-chrom or >1Mb same-chrom
            is_inter = chrom_a != chrom_b
            is_long = (not is_inter) and abs(pos_a - pos_b) > MIN_DISTANCE_SAME_CHROM
            if not (is_inter or is_long):
                continue
            cascade['inter_or_longrange'] += 1
            n_passed += 1

            # SAVE the junction with gene pair info
            group = bc_to_group.get(full_bc, "Unknown")
            all_junctions.append({
                'barcode': full_bc,
                'group': group,
                'geneA': gene_a,
                'geneB': gene_b,
                'chrom_a': chrom_a,
                'chrom_b': chrom_b,
                'pos_a': pos_a,
                'pos_b': pos_b,
                'junction_type': jtype,
                'aln_score_ratio': round(ratio, 3),
                'srr': srr,
            })

    sample_summaries.append({'srr': srr, 'total': n_total,
                             'target': n_target, 'passed': n_passed})
    log(f"  {srr}: {n_total:,} total -> {n_target:,} target -> {n_passed:,} passed")

# Save cascade
log(f"\n  Filter cascade:")
for step in ['total', 'in_target', 'unique_mapped', 'good_aln_score',
             'both_genic', 'both_protein_coding', 'no_artifacts', 'inter_or_longrange']:
    pct = 100 * cascade[step] / cascade['total'] if cascade['total'] > 0 else 0
    log(f"    {step:25s}: {cascade[step]:>12,} ({pct:.3f}%)")

# Save junction-level table
jdf = pd.DataFrame(all_junctions)
junction_path = os.path.join(OUTPUT_DIR, "all_filtered_junctions.tsv")
jdf.to_csv(junction_path, sep='\t', index=False)
log(f"\n  Saved {len(jdf)} junctions to {junction_path}")

# Per-group summary
log(f"\n  Per-group junction counts:")
for g in GROUPS:
    n = (jdf['group'] == g).sum() if len(jdf) > 0 else 0
    log(f"    {g}: {n}")

# =============================================================================
# STEP 4: FIND GROUP-EXCLUSIVE RECURRENT GENE PAIRS
# =============================================================================
banner("STEP 4: Group-exclusive recurrent gene pairs")

# Build canonical pairs
if len(jdf) > 0:
    jdf['_pair'] = jdf.apply(lambda r: tuple(sorted([r['geneA'], r['geneB']])), axis=1)

pair_cells = {g: defaultdict(set) for g in GROUPS}
pair_counts = {g: Counter() for g in GROUPS}

if len(jdf) > 0:
    for _, row in jdf.iterrows():
        grp = row['group']
        if grp not in GROUPS:
            continue
        pair_cells[grp][row['_pair']].add(row['barcode'])
        pair_counts[grp][row['_pair']] += 1

recurrent = {}
for grp in GROUPS:
    recurrent[grp] = {p for p, cells in pair_cells[grp].items()
                       if len(cells) >= MIN_CELLS_RECURRENT}
    log(f"  {grp}: {len(pair_counts[grp])} unique pairs, "
        f"{len(recurrent[grp])} recurrent (>={MIN_CELLS_RECURRENT} cells)")

exclusive = {}
for grp in GROUPS:
    others = set().union(*(recurrent[g] for g in GROUPS if g != grp))
    exclusive[grp] = recurrent[grp] - others
    log(f"  {grp} exclusive: {len(exclusive[grp])} pairs")

# Save exclusive pairs
for grp in GROUPS:
    rows = [{'geneA': p[0], 'geneB': p[1],
             'n_cells': len(pair_cells[grp][p]),
             'n_events': pair_counts[grp][p]} for p in sorted(exclusive[grp])]
    df = pd.DataFrame(rows).sort_values('n_cells', ascending=False) if rows else pd.DataFrame()
    df.to_csv(os.path.join(OUTPUT_DIR, f"{grp}_exclusive_pairs.tsv"), sep='\t', index=False)
    log(f"\n  {grp}: {len(df)} exclusive pairs saved")
    for _, r in df.head(10).iterrows():
        log(f"    {r['geneA']} — {r['geneB']} ({r['n_cells']} cells, {r['n_events']} events)")

# =============================================================================
# STEP 5: EXTRACT UNIQUE GENES
# =============================================================================
banner("STEP 5: Unique genes from exclusive pairs")

exclusive_genes = {}
for grp in GROUPS:
    genes = set()
    for p in exclusive[grp]:
        genes.update(p)
    exclusive_genes[grp] = sorted(genes)
    log(f"  {grp}: {len(exclusive_genes[grp])} unique genes")
    with open(os.path.join(OUTPUT_DIR, f"{grp}_exclusive_genes.txt"), 'w') as f:
        f.write('\n'.join(exclusive_genes[grp]))

# =============================================================================
# STEP 6: ENRICHMENT ANALYSIS
# =============================================================================
banner("STEP 6: Enrichment analysis (gseapy)")

try:
    import gseapy as gp
    HAS_GSEAPY = True
    log("  gseapy loaded")
except ImportError:
    HAS_GSEAPY = False
    log("  gseapy not available, skipping enrichment")

if HAS_GSEAPY:
    LIBS = ['KEGG_2021_Human', 'Reactome_2022', 'GO_Biological_Process_2023']
    for grp in ['SBS2_HIGH', 'Stealth_CNV']:
        gl = exclusive_genes.get(grp, [])
        if len(gl) < 5:
            log(f"\n  {grp}: {len(gl)} genes, too few for enrichment")
            continue
        log(f"\n  --- {grp} ({len(gl)} genes) ---")
        for lib in LIBS:
            try:
                enr = gp.enrichr(gene_list=gl, gene_sets=lib, organism='human',
                                 outdir=None, no_plot=True)
                res = enr.results
                sig = res[res['Adjusted P-value'] < 0.05]
                log(f"    {lib}: {len(sig)} significant / {len(res)} total")
                res.to_csv(os.path.join(OUTPUT_DIR,
                    f"{grp}_enrichment_{lib.split('_')[0]}.tsv"), sep='\t', index=False)
                show = sig if len(sig) > 0 else res
                for _, r in show.head(5).iterrows():
                    pre = "" if r['Adjusted P-value'] < 0.05 else "(ns) "
                    log(f"      {pre}{r['Term']}: p_adj={r['Adjusted P-value']:.2e}, "
                        f"genes={r['Genes']}")
            except Exception as e:
                log(f"    {lib}: ERROR - {e}")

# =============================================================================
# STEP 7: NOTABLE GENES
# =============================================================================
banner("STEP 7: Notable genes in exclusive fusions")

spliceosome_genes = {'HNRNPA1', 'HNRNPA2B1', 'HNRNPC', 'SRSF1', 'SRSF2', 'SRSF3',
                     'SFPQ', 'SNRPD2', 'SNRPC', 'U2AF1', 'U2AF2', 'SF3B1',
                     'PRPF8', 'DDX5', 'DDX17', 'NONO'}
immune_genes = {'HLA-A', 'HLA-B', 'HLA-C', 'B2M', 'TAP1', 'TAP2', 'CD274',
                'IDO1', 'PARP1', 'TP53'}

for grp in ['SBS2_HIGH', 'Stealth_CNV']:
    gs = set(exclusive_genes.get(grp, []))
    sh = gs & spliceosome_genes
    ih = gs & immune_genes
    if sh:
        log(f"  {grp} spliceosome genes in fusions: {sorted(sh)}")
    if ih:
        log(f"  {grp} immune genes in fusions: {sorted(ih)}")
    if not sh and not ih:
        log(f"  {grp}: no spliceosome/immune gene hits in exclusive fusions")

# =============================================================================
# STEP 8: SUPPLEMENTAL TABLE — ALL EXCLUSIVE PAIRS COMBINED
# =============================================================================
banner("STEP 8: Supplemental table")

all_exclusive = []
for grp in GROUPS:
    for p in exclusive[grp]:
        all_exclusive.append({
            'group': grp,
            'geneA': p[0],
            'geneB': p[1],
            'n_cells': len(pair_cells[grp][p]),
            'n_events': pair_counts[grp][p],
        })
supp_df = pd.DataFrame(all_exclusive).sort_values(['group', 'n_cells'], ascending=[True, False])
supp_path = os.path.join(OUTPUT_DIR, "Supplemental_Table_Exclusive_Fusion_Pairs.tsv")
supp_df.to_csv(supp_path, sep='\t', index=False)
log(f"  Saved supplemental table: {supp_path} ({len(supp_df)} total pairs)")

# Save report
report_path = os.path.join(OUTPUT_DIR, "reparse_chimeric_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"\n  Report: {report_path}")

banner("COMPLETE")
