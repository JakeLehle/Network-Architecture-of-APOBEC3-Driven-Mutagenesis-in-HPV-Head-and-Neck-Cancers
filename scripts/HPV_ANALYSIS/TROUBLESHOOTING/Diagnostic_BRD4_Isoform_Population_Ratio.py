#!/usr/bin/env python3
"""
Diagnostic_BRD4_Isoform_Population_Ratio.py   (Module 2, stage 2)
=================================================================
READ-ONLY. Population-stratified pseudobulk BRD4 L : S(a) : S(b) ratio.

Stage 1 established that per-cell isoform resolution is not supported by depth,
but that S(a) is cleanly recoverable pseudobulk (corrected run: L=125, S(a)=42,
S(b)=3 on one sample). This asks the question that matters for Chiang's model:
does the L : S(a) ratio shift across the lifecycle axis,
    SBS2-HIGH (maintenance)  vs  CNV-HIGH (productive)  vs  NORMAL?

Approach
--------
Reuses the CCDS-anchored grouping and isoform-unique intervals from stage 1
(L=CCDS12328, S(a)=CCDS46004, S(b)=CCDS82307; both L transcript models grouped).
Molecules are classified by GENOMIC overlap with the isoform-unique intervals,
deduplicated by (CB, UB) per sample, restricted to the curated group cells, and
tallied by population, pooled across all contributing samples.

Barcode join (confirmed from the repo): three_group_assignments.tsv barcodes are
'SEQUENCE-1-SRR#####'. SRR = last hyphen field (which BAM to open); raw CB =
'SEQUENCE-1' (matches the BAM CB tag). Group membership is the cell filter.

Test
----
2x3 chi-square on the (L, S(a)) x (SBS2, CNV, NORMAL) count table, plus pairwise
Fisher exact (BH-corrected). This is a PSEUDOBULK test on pooled molecule counts;
pooling across cells ignores cell-level clustering, so it is reported as
suggestive at the population level, not as a per-cell-grade result. S(b) is
reported descriptively (rare form), not tested.

INPUTS (read-only)
------------------
  GTF    : SC/ref/GRCh38/genes/genes_unzipped.gtf
  groups : 2026_NMF_PAPER/data/FIG_4/01_group_selection/three_group_assignments.tsv
  BAMs   : .../GSE173468/<SRR>/<SRR>_S1_L001_/outs/possorted_genome_bam.bam

OUTPUT (to 2026_NMF_PAPER/data/FIG_6/DIAGNOSTIC_BRD4_ISOFORM/)
--------------------------------------------------------------
  brd4_population_ratio.txt              full report
  brd4_population_ratio_by_srr.tsv       per-SRR x population contribution table
  brd4_population_ratio.pdf/.png         composition + L/S(a) ratio bars

Env: needs pysam (NEOANTIGEN).
  conda run -n NEOANTIGEN python Diagnostic_BRD4_Isoform_Population_Ratio.py

Author: Jake Lehle
"""

import os
import re
from collections import defaultdict, OrderedDict

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    import pysam
except ImportError:
    import sys
    sys.exit("ERROR: pysam not found. Run in the NEOANTIGEN env.")

# =============================================================================
# CONFIG
# =============================================================================
REF_DIR    = "/master/jlehle/WORKING/SC/ref/GRCh38"
GTF_PATH   = os.path.join(REF_DIR, "genes/genes_unzipped.gtf")

GSE_ROOT   = ("/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/"
              "results_NMF_v0.1.1/fastq/GSE173468")
THREE_GROUP_PATH = ("/master/jlehle/WORKING/2026_NMF_PAPER/"
                    "data/FIG_4/01_group_selection/three_group_assignments.tsv")
OUTPUT_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_6/DIAGNOSTIC_BRD4_ISOFORM"
os.makedirs(OUTPUT_DIR, exist_ok=True)

GENE_NAME = "BRD4"
LOCUS_PAD = 2000
POP_ORDER = ['SBS2_HIGH', 'CNV_HIGH', 'NORMAL']
POP_LABELS = {'SBS2_HIGH': 'SBS2-HIGH', 'CNV_HIGH': 'CNV-HIGH', 'NORMAL': 'Normal'}
POP_COLORS = {'SBS2_HIGH': '#ed6a5a', 'CNV_HIGH': '#F6D155', 'NORMAL': '#5B7C99'}

# authoritative CCDS anchors (BRD4 review / UniProt)
CCDS_L, CCDS_SA, CCDS_SB = "CCDS12328", "CCDS46004", "CCDS82307"

MIN_SA_FOR_TEST = 10   # per-population S(a) floor below which we call it underpowered
DPI = 300

# =============================================================================
# LOGGING
# =============================================================================
report = []
def log(msg=""):
    print(msg, flush=True)
    report.append(str(msg))
def banner(t):
    log(""); log("=" * 78); log(f"  {t}"); log("=" * 78)

# =============================================================================
# INTERVAL + GTF HELPERS  (mirror stage 1 exactly)
# =============================================================================
def merge(ivs):
    ivs = sorted([iv for iv in ivs if iv[1] > iv[0]])
    if not ivs:
        return []
    out = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]

def subtract(A, B):
    A = merge(A); B = merge(B); out = []
    for s, e in A:
        cur = s
        for bs, be in B:
            if be <= cur or bs >= e:
                continue
            if bs > cur:
                out.append((cur, min(bs, e)))
            cur = max(cur, be)
            if cur >= e:
                break
        if cur < e:
            out.append((cur, e))
    return [iv for iv in out if iv[1] > iv[0]]

def total_bp(ivs):
    return sum(e - s for s, e in ivs)

def overlaps_any(b0, b1, ivs):
    for s, e in ivs:
        if b0 < e and s < b1:
            return True
    return False

def attr(field, key):
    m = re.search(key + r' "([^"]+)"', field)
    return m.group(1) if m else None

def build_isoform_intervals(gtf_path, gene_name):
    """Parse BRD4, group by CCDS, return isoform-unique intervals (0-based)."""
    txs = OrderedDict()
    tx_meta = {}
    seqname = strand = None
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#') or '\t' not in line:
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or ('gene_name "%s"' % gene_name) not in f[8]:
                continue
            if seqname is None:
                seqname, strand = f[0], f[6]
            tid = attr(f[8], 'transcript_id')
            if f[2] == 'transcript' and tid:
                tx_meta[tid] = attr(f[8], 'ccds_id') or attr(f[8], 'ccdsid')
            elif f[2] == 'exon' and tid:
                txs.setdefault(tid, []).append((int(f[3]) - 1, int(f[4])))
    groups = {'L': [], 'Sa': [], 'Sb': []}
    background = []
    for tid in txs:
        c = (tx_meta.get(tid) or '').split('.')[0]
        if   c == CCDS_L:  groups['L'].append(tid)
        elif c == CCDS_SA: groups['Sa'].append(tid)
        elif c == CCDS_SB: groups['Sb'].append(tid)
        else:              background.append(tid)

    def gex(ids):
        ivs = []
        for tid in ids:
            ivs += txs[tid]
        return merge(ivs)

    L_ex, Sa_ex, Sb_ex, bg_ex = (gex(groups['L']), gex(groups['Sa']),
                                 gex(groups['Sb']), gex(background))
    uniq = {
        'L':  subtract(L_ex,  merge(Sa_ex + Sb_ex + bg_ex)),
        'Sa': subtract(Sa_ex, merge(L_ex + Sb_ex + bg_ex)),
        'Sb': subtract(Sb_ex, merge(L_ex + Sa_ex + bg_ex)),
    }
    all_ex = merge(L_ex + Sa_ex + Sb_ex + bg_ex)
    locus_lo = min(s for s, _ in all_ex) - LOCUS_PAD
    locus_hi = max(e for _, e in all_ex) + LOCUS_PAD
    return seqname, strand, groups, uniq, locus_lo, locus_hi

def bh_adjust(pvals):
    p = np.asarray(pvals, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = ~np.isnan(p); n = int(ok.sum())
    if n == 0:
        return q
    idx = np.where(ok)[0]; pv = p[idx]
    order = np.argsort(pv); ranked = pv[order]
    raw = ranked * n / np.arange(1, n + 1)
    q_sorted = np.clip(np.minimum.accumulate(raw[::-1])[::-1], 0, 1)
    qb = np.empty(n); qb[order] = q_sorted; q[idx] = qb
    return q

# =============================================================================
# STEP 0: barcode -> (SRR, raw CB, population)
# =============================================================================
banner("STEP 0: Load groups + build SRR/CB/population map")
groups_df = pd.read_csv(THREE_GROUP_PATH, sep='\t')
log(f"  three_group_assignments: {len(groups_df)} cells")
log(f"  columns: {list(groups_df.columns)}")

# --- resolve a patient / donor label per cell ------------------------------
# The CMH in the Figure 6 harness must stratify by PATIENT, not by SRR: samples
# nest within patients, and only one patient (SC001) contributes to both tumour
# populations, so an SRR-stratified test can silently collapse onto one donor.
MASTER_TABLE_PATH = ("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_6/"
                     "01_raw_hpv16_counts/basal_cell_master_table_with_raw_HPV16.tsv")
# The master table calls it 'subject id', with a space. Ordered so a true
# subject-level column always wins over anything sample-level.
PATIENT_CANDIDATES = ['subject id', 'subject_id', 'subject',
                      'patient', 'patient_id', 'donor', 'donor_id']
RUN_CANDIDATES = ['run_accession', 'srr', 'run']

def _find_col(cols, cands):
    low = {str(c).strip().lower(): c for c in cols}
    for c in cands:
        if c in low:
            return low[c]
    return None

cell_to_patient = {}
pcol = _find_col(groups_df.columns, PATIENT_CANDIDATES)
if pcol:
    cell_to_patient = dict(zip(groups_df['cell_barcode'], groups_df[pcol].astype(str)))
    log(f"  patient column found in group table: '{pcol}'")
else:
    log("  no patient column in the group table; trying the master table ...")
    try:
        _m = pd.read_csv(MASTER_TABLE_PATH, sep='\t', index_col=0, nrows=5)
        mcol = _find_col(_m.columns, PATIENT_CANDIDATES)
        rcol = _find_col(_m.columns, RUN_CANDIDATES)
        if mcol:
            # Read by NAME, not by position. usecols=[0, idx+1] assumes the
            # index sits at column 0 and nothing has shifted; a column name
            # containing a space makes that assumption harder to verify.
            # The table is 52k rows, so reading it whole costs nothing.
            _mfull = pd.read_csv(MASTER_TABLE_PATH, sep='\t', index_col=0)
            cell_to_patient = {str(k): str(v)
                               for k, v in _mfull[mcol].to_dict().items()}
            log(f"  patient column found in the master table: '{mcol}' "
                f"({len(set(cell_to_patient.values()))} distinct subjects)")

            # Cross-check: the SRR parsed out of the barcode should equal the
            # run_accession recorded for that cell. This is the one assumption
            # the whole barcode join rests on, and it has never been tested.
            if rcol:
                _run = {str(k): str(v) for k, v in _mfull[rcol].to_dict().items()}
                checked = mismatched = 0
                for _srr, _cells in srr_cell_pop.items():
                    for _raw_cb in _cells:
                        _bc = f"{_raw_cb}-{_srr}"
                        if _bc in _run:
                            checked += 1
                            if _run[_bc] != _srr:
                                mismatched += 1
                if checked == 0:
                    log("  WARNING: no barcode matched the master table index; "
                        "the barcode join may be wrong. Do not trust the "
                        "patient labels below.")
                elif mismatched:
                    log(f"  ERROR: {mismatched}/{checked} barcodes disagree with "
                        f"'{rcol}'. The SEQ-1-SRR parse is not reliable.")
                else:
                    log(f"  barcode-derived SRR matches '{rcol}' for all "
                        f"{checked} checked cells")
        else:        
            log(f"  master table columns: {list(_m.columns)}")
            log("  NO PATIENT COLUMN ANYWHERE. The CMH will fall back to SRR")
            log("  stratification, which is NOT equivalent. Resolve this before")
            log("  any isoform number goes into the manuscript.")
    except Exception as e:
        log(f"  master table lookup failed ({e}); falling back to SRR")

srr_cell_pop = defaultdict(dict)   # SRR -> {raw_CB: population}
n_per_pop = defaultdict(int)
bad = 0
for bc, grp in zip(groups_df['cell_barcode'], groups_df['group']):
    parts = bc.split('-')
    if len(parts) < 3:
        bad += 1
        continue
    srr = parts[-1]                       # e.g. SRR14340900
    raw_cb = "-".join(parts[:-1])         # e.g. AAACCTGAGAAACCAT-1
    srr_cell_pop[srr][raw_cb] = grp
    n_per_pop[grp] += 1
if bad:
    log(f"  WARNING: {bad} barcodes did not match SEQ-1-SRR format")
log(f"  SRRs contributing group cells: {len(srr_cell_pop)}")
for p in POP_ORDER:
    log(f"    {POP_LABELS[p]:<10s}: {n_per_pop.get(p, 0)} cells")

# =============================================================================
# STEP 1: isoform-unique intervals
# =============================================================================
banner("STEP 1: BRD4 isoform-unique intervals (CCDS-anchored)")
seqname, strand, iso_groups, uniq, locus_lo, locus_hi = build_isoform_intervals(GTF_PATH, GENE_NAME)
log(f"  BRD4 seqname '{seqname}', strand '{strand}'")
for g, anchor in (('L', CCDS_L), ('Sa', CCDS_SA), ('Sb', CCDS_SB)):
    log(f"    {g:<3s} [{anchor}]: {iso_groups[g] or 'NONE'}  "
        f"unique={total_bp(uniq[g])} bp")
if total_bp(uniq['L']) == 0 or total_bp(uniq['Sa']) == 0:
    import sys
    sys.exit("ERROR: L or S(a) unique interval empty; cannot compute ratio.")

# =============================================================================
# STEP 2: walk each SRR, classify molecules restricted to group cells
# =============================================================================
banner("STEP 2: Per-sample discriminating-molecule counting (group cells only)")
# pop -> class -> count ; and per-SRR breakdown
pop_counts = {p: {'L': 0, 'Sa': 0, 'Sb': 0} for p in POP_ORDER}
pop_cells_hit = {p: set() for p in POP_ORDER}
srr_rows = []
n_srr_ok = 0
missing_bam = []
no_index = []
for srr in sorted(srr_cell_pop):
    bam_path = os.path.join(GSE_ROOT, srr, f"{srr}_S1_L001_", "outs",
                            "possorted_genome_bam.bam")
    cellmap = srr_cell_pop[srr]           # raw_CB -> population
    if not os.path.exists(bam_path):
        missing_bam.append(srr)
        log(f"  MISSING BAM: {srr}  ({len(cellmap)} group cells skipped)")
        continue
    bam = pysam.AlignmentFile(bam_path, "rb")
    if not bam.has_index():               # fetch() requires a .bai/.csi
        no_index.append(srr)
        bam.close()
        sys.exit(f"ERROR: no .bai for {srr}. These counts feed a manuscript "
                 f"number, so a silently skipped sample is not acceptable. "
                 f"Index it with the conda samtools and re-run.")    
    # resolve contig for this BAM
    refs = set(bam.references)
    contig = next((c for c in (seqname, f"chr{seqname}", seqname.replace("chr", ""))
                   if c in refs), None)
    if contig is None:
        log(f"  {srr}: contig {seqname} not in BAM; skipped")
        bam.close()
        continue
    n_srr_ok += 1
    mol_class = defaultdict(set)          # (CB,UB) -> classes, this SRR
    for r in bam.fetch(contig, locus_lo, locus_hi):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        if not (r.has_tag('CB') and r.has_tag('UB')):
            continue
        cb = r.get_tag('CB')
        if cb not in cellmap:             # group cells only
            continue
        cls = set()
        for b0, b1 in r.get_blocks():
            for k in ('L', 'Sa', 'Sb'):
                if uniq[k] and overlaps_any(b0, b1, uniq[k]):
                    cls.add(k)
        if cls:
            mol_class[(cb, r.get_tag('UB'))] |= cls
    # resolve molecules for this SRR
    srr_tally = {p: {'L': 0, 'Sa': 0, 'Sb': 0} for p in POP_ORDER}
    for (cb, ub), cls in mol_class.items():
        if len(cls) != 1:
            continue
        c = next(iter(cls))
        pop = cellmap[cb]
        pop_counts[pop][c] += 1
        srr_tally[pop][c] += 1
        pop_cells_hit[pop].add((srr, cb))
    for p in POP_ORDER:
        t = srr_tally[p]
        if t['L'] + t['Sa'] + t['Sb'] > 0:
            # patient label for this SRR. Take the modal patient across its
            # group cells; warn loudly if one SRR maps to more than one.
            pats = set()
            for raw_cb, gp in cellmap.items():
                full_bc = f"{raw_cb}-{srr}"
                if full_bc in cell_to_patient:
                    pats.add(cell_to_patient[full_bc])
                elif raw_cb in cell_to_patient:
                    pats.add(cell_to_patient[raw_cb])
            if len(pats) > 1:
                log(f"  WARNING: {srr} maps to multiple patients {sorted(pats)}; "
                    f"using the first. Check the barcode join.")
            patient = sorted(pats)[0] if pats else 'UNKNOWN'
            srr_rows.append({'SRR': srr, 'patient': patient, 'population': p,
                             'L': t['L'], 'Sa': t['Sa'], 'Sb': t['Sb'],
                             'n_group_cells': sum(1 for v in cellmap.values() if v == p)})            
    bam.close()
log(f"  processed {n_srr_ok}/{len(srr_cell_pop)} sample BAMs")
if no_index:
    log(f"  {len(no_index)} sample(s) SKIPPED for a missing .bai index: {no_index}")
    log(f"    to include them, index with the conda samtools (system one is broken), "
        f"then re-run, e.g.:")
    log(f"    for s in {' '.join(no_index)}; do samtools index "
        f"{GSE_ROOT}/$s/${{s}}_S1_L001_/outs/possorted_genome_bam.bam; done")
if missing_bam:
    log(f"  {len(missing_bam)} sample(s) had no BAM at all: {missing_bam}")
if n_srr_ok < len(srr_cell_pop):
    log(f"  NOTE: pooled counts below reflect only the {n_srr_ok} indexed samples; "
        f"index the skipped ones and re-run for full coverage before trusting power.")

# =============================================================================
# STEP 3: pooled per-population ratio + tests
# =============================================================================
banner("STEP 3: Pooled per-population pseudobulk ratio")
log(f"  {'population':<12s} {'L':>6s} {'S(a)':>6s} {'S(b)':>6s} "
    f"{'L/S(a)':>8s} {'S(a)frac':>9s} {'cells hit':>10s}")
for p in POP_ORDER:
    c = pop_counts[p]
    tot = c['L'] + c['Sa']
    ratio = (c['L'] / c['Sa']) if c['Sa'] > 0 else float('inf')
    safrac = (c['Sa'] / tot) if tot > 0 else float('nan')
    log(f"  {POP_LABELS[p]:<12s} {c['L']:>6d} {c['Sa']:>6d} {c['Sb']:>6d} "
        f"{ratio:>8.2f} {safrac:>9.3f} {len(pop_cells_hit[p]):>10d}")

# 2x3 chi-square on (L, Sa) x populations
table = np.array([[pop_counts[p]['L'] for p in POP_ORDER],
                  [pop_counts[p]['Sa'] for p in POP_ORDER]], dtype=float)
log("")
thin = [p for p in POP_ORDER if pop_counts[p]['Sa'] < MIN_SA_FOR_TEST]
if thin:
    log(f"  NOTE: S(a) below {MIN_SA_FOR_TEST} in {', '.join(POP_LABELS[p] for p in thin)}; "
        f"treat the test as underpowered / suggestive only.")
if (table.sum(axis=0) > 0).all() and table.sum() > 0:
    try:
        chi2, pchi, dof, _ = chi2_contingency(table)
        log(f"  2x3 chi-square (L vs S(a) across populations): "
            f"chi2={chi2:.2f}, dof={dof}, p={pchi:.3g}")
    except ValueError as e:
        log(f"  chi-square not computable: {e}")
    # pairwise Fisher on (L, Sa)
    pairs = [('SBS2_HIGH', 'CNV_HIGH'), ('CNV_HIGH', 'NORMAL'), ('SBS2_HIGH', 'NORMAL')]
    praw = []
    for a, b in pairs:
        t2 = [[pop_counts[a]['L'], pop_counts[b]['L']],
              [pop_counts[a]['Sa'], pop_counts[b]['Sa']]]
        try:
            _, pf = fisher_exact(t2)
        except ValueError:
            pf = np.nan
        praw.append(pf)
    q = bh_adjust(praw)
    log("  pairwise Fisher (L vs S(a)), BH-corrected:")
    for (a, b), pr, qq in zip(pairs, praw, q):
        log(f"    {POP_LABELS[a]} vs {POP_LABELS[b]}: p={pr:.3g}, q={qq:.3g}")
else:
    log("  insufficient counts for a contingency test.")

# =============================================================================
# STEP 4: outputs
# =============================================================================
banner("STEP 4: Write outputs")
srr_df = pd.DataFrame(srr_rows)
srr_path = os.path.join(OUTPUT_DIR, "brd4_population_ratio_by_srr.tsv")
srr_df.to_csv(srr_path, sep='\t', index=False)
log(f"  [SAVE] {srr_path}")

# plot: composition (left) + L/S(a) ratio (right)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
x = np.arange(len(POP_ORDER))
Lv = np.array([pop_counts[p]['L'] for p in POP_ORDER], float)
Sav = np.array([pop_counts[p]['Sa'] for p in POP_ORDER], float)
Sbv = np.array([pop_counts[p]['Sb'] for p in POP_ORDER], float)
tot = np.maximum(Lv + Sav + Sbv, 1)
ax1.bar(x, Lv / tot, color='#4c72b0', label='L')
ax1.bar(x, Sav / tot, bottom=Lv / tot, color='#dd8452', label='S(a)')
ax1.bar(x, Sbv / tot, bottom=(Lv + Sav) / tot, color='#b0b0b0', label='S(b)')
ax1.set_xticks(x); ax1.set_xticklabels([POP_LABELS[p] for p in POP_ORDER], fontsize=10)
ax1.set_ylabel('isoform fraction (pseudobulk)', fontsize=11)
ax1.set_title('BRD4 isoform composition', fontsize=11)
ax1.legend(fontsize=9, frameon=False)
for i, p in enumerate(POP_ORDER):
    ax1.text(i, 1.02, f"n={int(Lv[i]+Sav[i]+Sbv[i])}", ha='center', fontsize=8)

ratio = [ (pop_counts[p]['L'] / pop_counts[p]['Sa']) if pop_counts[p]['Sa'] > 0 else np.nan
          for p in POP_ORDER ]
ax2.bar(x, ratio, color=[POP_COLORS[p] for p in POP_ORDER])
ax2.set_xticks(x); ax2.set_xticklabels([POP_LABELS[p] for p in POP_ORDER], fontsize=10)
ax2.set_ylabel('L / S(a)', fontsize=11)
ax2.set_title('BRD4 L : S(a) ratio', fontsize=11)
for i, rv in enumerate(ratio):
    if not np.isnan(rv):
        ax2.text(i, rv, f"{rv:.2f}", ha='center', va='bottom', fontsize=9)
plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(os.path.join(OUTPUT_DIR, f"brd4_population_ratio.{ext}"),
                dpi=DPI, bbox_inches='tight')
plt.close(fig)
log(f"  [SAVE] brd4_population_ratio.pdf/.png")

rp = os.path.join(OUTPUT_DIR, "brd4_population_ratio.txt")
with open(rp, 'w') as fh:
    fh.write("\n".join(report))
print(f"\n[SAVE] {rp}")
print("BRD4 POPULATION RATIO COMPLETE")
