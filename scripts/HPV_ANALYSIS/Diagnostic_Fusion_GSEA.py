#!/usr/bin/env python3
"""
Diagnostic_Fusion_GSEA.py
==========================
Run enrichment analysis on group-exclusive chimeric fusion gene partners.
Reads all_filtered_junctions.tsv produced by Reparse_Chimeric_With_GenePairs.py.

Env: NETWORK
Usage: conda run -n NETWORK python Diagnostic_Fusion_GSEA.py
"""

import os, sys
import pandas as pd
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
JUNCTIONS_PATH = os.path.join(PROJECT_ROOT,
    "data/FIG_6/05_neoantigen/fusion_detection/gsea_diagnostic/all_filtered_junctions.tsv")
OUTPUT_DIR = os.path.join(PROJECT_ROOT,
    "data/FIG_6/05_neoantigen/fusion_detection/gsea_diagnostic")
os.makedirs(OUTPUT_DIR, exist_ok=True)

GROUPS = ['SBS2_HIGH', 'Stealth_CNV', 'Normal_Control']
MIN_CELLS_RECURRENT = 2

report_lines = []
def log(msg=""):
    print(msg, flush=True)
    report_lines.append(str(msg))
def banner(title, char="="):
    log(""); log(char*80); log(f"  {title}"); log(char*80)

# ── STEP 1: LOAD JUNCTIONS ──────────────────────────────────────────────────
banner("STEP 1: Load filtered junctions")

if not os.path.exists(JUNCTIONS_PATH):
    log(f"  ERROR: {JUNCTIONS_PATH} not found.")
    log("  Run Reparse_Chimeric_With_GenePairs.py first.")
    sys.exit(1)

jdf = pd.read_csv(JUNCTIONS_PATH, sep='\t')
log(f"  Junctions: {jdf.shape}")
log(f"  Columns: {list(jdf.columns)}")
log(f"  Groups: {jdf['group'].value_counts().to_dict()}")

# ── STEP 2: BUILD GENE PAIR COUNTS ──────────────────────────────────────────
banner("STEP 2: Build gene pair counts per group")

jdf['_pair'] = jdf.apply(lambda r: tuple(sorted([str(r['geneA']), str(r['geneB'])])), axis=1)

pair_cells = {g: defaultdict(set) for g in GROUPS}
pair_counts = {g: Counter() for g in GROUPS}

for _, row in jdf.iterrows():
    grp = row['group']
    if grp not in GROUPS:
        continue
    pair_cells[grp][row['_pair']].add(row['barcode'])
    pair_counts[grp][row['_pair']] += 1

for g in GROUPS:
    log(f"  {g}: {len(pair_counts[g])} unique pairs, "
        f"{sum(pair_counts[g].values())} total events")

# ── STEP 3: EXCLUSIVE RECURRENT PAIRS ───────────────────────────────────────
banner("STEP 3: Group-exclusive recurrent pairs")

recurrent = {}
for grp in GROUPS:
    recurrent[grp] = {p for p, cells in pair_cells[grp].items()
                       if len(cells) >= MIN_CELLS_RECURRENT}
    log(f"  {grp} recurrent (>={MIN_CELLS_RECURRENT} cells): {len(recurrent[grp])}")

exclusive = {}
for grp in GROUPS:
    others = set().union(*(recurrent[g] for g in GROUPS if g != grp))
    exclusive[grp] = recurrent[grp] - others
    log(f"  {grp} exclusive: {len(exclusive[grp])}")

for grp in GROUPS:
    rows = [{'geneA': p[0], 'geneB': p[1],
             'n_cells': len(pair_cells[grp][p]),
             'n_events': pair_counts[grp][p]} for p in sorted(exclusive[grp])]
    df = pd.DataFrame(rows).sort_values('n_cells', ascending=False) if rows else pd.DataFrame()
    df.to_csv(os.path.join(OUTPUT_DIR, f"{grp}_exclusive_pairs.tsv"), sep='\t', index=False)
    log(f"\n  {grp}: {len(df)} exclusive pairs")
    for _, r in df.head(10).iterrows():
        log(f"    {r['geneA']} — {r['geneB']} ({r['n_cells']} cells, {r['n_events']} events)")

# ── STEP 4: UNIQUE GENES ────────────────────────────────────────────────────
banner("STEP 4: Unique genes from exclusive pairs")

exclusive_genes = {}
for grp in GROUPS:
    genes = set()
    for p in exclusive[grp]:
        genes.update(p)
    exclusive_genes[grp] = sorted(genes)
    log(f"  {grp}: {len(exclusive_genes[grp])} unique genes")
    with open(os.path.join(OUTPUT_DIR, f"{grp}_exclusive_genes.txt"), 'w') as f:
        f.write('\n'.join(exclusive_genes[grp]))

# ── STEP 5: ENRICHMENT ──────────────────────────────────────────────────────
banner("STEP 5: Enrichment analysis (gseapy)")

try:
    import gseapy as gp
    HAS_GSEAPY = True
    log("  gseapy loaded")
except ImportError:
    HAS_GSEAPY = False
    log("  gseapy not available, skipping")

if HAS_GSEAPY:
    LIBS = ['KEGG_2021_Human', 'Reactome_2022', 'GO_Biological_Process_2023']
    for grp in ['SBS2_HIGH', 'Stealth_CNV']:
        gl = exclusive_genes.get(grp, [])
        if len(gl) < 5:
            log(f"\n  {grp}: {len(gl)} genes, too few")
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

# ── STEP 6: NOTABLE GENES ───────────────────────────────────────────────────
banner("STEP 6: Notable genes in exclusive fusions")

spliceosome = {'HNRNPA1','HNRNPA2B1','HNRNPC','SRSF1','SRSF2','SRSF3',
               'SFPQ','SNRPD2','SNRPC','U2AF1','U2AF2','SF3B1',
               'PRPF8','DDX5','DDX17','NONO'}
immune = {'HLA-A','HLA-B','HLA-C','B2M','TAP1','TAP2','CD274','PARP1','TP53'}

for grp in ['SBS2_HIGH', 'Stealth_CNV']:
    gs = set(exclusive_genes.get(grp, []))
    sh, ih = gs & spliceosome, gs & immune
    if sh: log(f"  {grp} spliceosome: {sorted(sh)}")
    if ih: log(f"  {grp} immune: {sorted(ih)}")
    if not sh and not ih: log(f"  {grp}: no spliceosome/immune hits")

# ── STEP 7: SHARED PAIRS SUMMARY ────────────────────────────────────────────
banner("STEP 7: Shared recurrent pairs")

shared_disease = recurrent.get('SBS2_HIGH', set()) & recurrent.get('Stealth_CNV', set())
log(f"  SBS2_HIGH ∩ Stealth_CNV: {len(shared_disease)} shared recurrent pairs")
shared_all = shared_disease & recurrent.get('Normal_Control', set())
log(f"  All three groups: {len(shared_all)} shared recurrent pairs")

if len(shared_disease) > 0:
    log(f"\n  Top shared SBS2_HIGH ∩ Stealth_CNV pairs:")
    shared_rows = []
    for p in shared_disease:
        h_cells = len(pair_cells['SBS2_HIGH'][p])
        s_cells = len(pair_cells['Stealth_CNV'][p])
        shared_rows.append({'geneA': p[0], 'geneB': p[1],
                            'SBS2_HIGH_cells': h_cells, 'Stealth_CNV_cells': s_cells})
    sdf = pd.DataFrame(shared_rows).sort_values('SBS2_HIGH_cells', ascending=False)
    sdf.to_csv(os.path.join(OUTPUT_DIR, "shared_disease_recurrent_pairs.tsv"),
               sep='\t', index=False)
    for _, r in sdf.head(10).iterrows():
        log(f"    {r['geneA']} — {r['geneB']} "
            f"(HIGH:{r['SBS2_HIGH_cells']}, Stealth:{r['Stealth_CNV_cells']})")

# Save report
report_path = os.path.join(OUTPUT_DIR, "fusion_gsea_diagnostic_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"\n  Report: {report_path}")
banner("DIAGNOSTIC COMPLETE")
