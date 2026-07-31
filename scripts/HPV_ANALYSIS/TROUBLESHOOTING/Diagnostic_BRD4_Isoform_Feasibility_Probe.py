#!/usr/bin/env python3
"""
Diagnostic_BRD4_Isoform_Feasibility_Probe.py   (Module 2, stage 1)
==================================================================
READ-ONLY. Answers the two questions that gate the BRD4 L/S isoform module
before any per-cell machinery is built:

  Q1. Are BRD4-L, BRD4-S(a), BRD4-S(b) resolvable as distinct transcripts in the
      transcriptome GTF Cell Ranger actually used, so we can define isoform-unique
      genomic intervals?
  Q2. On a real 3' possorted BAM, how many discriminating UMIs land in those
      unique intervals, pseudobulk and per cell? i.e. is per-cell L/S resolution
      feasible, or do we report pseudobulk only?

Nothing is trusted from memory. The contig name (19 vs chr19), the BRD4 transcript
identities, and the unique intervals are all resolved at runtime from the GTF and
BAM header. Molecules are classified by GENOMIC overlap with the isoform-unique
intervals, NOT by Cell Ranger's gene tags, because the discriminating reads sit in
3' UTR extensions that CR's xf/RE tags may not flag as transcriptomic.

Chemistry is confirmed 3' (Jake), so each isoform's poly-A-proximal reads pile in
its own unique 3' region: BRD4 is minus-strand, so the 3' end (poly-A) is at the
LOW coordinate; BRD4-L terminates lowest, BRD4-S(a)'s unique terminal exon / 3' UTR
sits at a HIGHER coordinate. That is exactly where the reads concentrate, which is
why this is worth probing.

INPUTS (read-only)
------------------
  GTF        : SC/ref/GRCh38/genes/genes_unzipped.gtf
  reference  : SC/ref/GRCh38/reference.json   (build cross-check vs Ensembl 115)
  one sample : .../GSE173468/<SRR>/<SRR>_S1_L001_/outs/
                 possorted_genome_bam.bam (+ .bai)
                 filtered_feature_bc_matrix/barcodes.tsv.gz

OUTPUT (to 2026_NMF_PAPER/data/FIG_6/DIAGNOSTIC_BRD4_ISOFORM/)
--------------------------------------------------------------
  brd4_isoform_feasibility_<SRR>.txt   full report
  brd4_percell_umi_hist_<SRR>.pdf/.png per-cell discriminating-UMI histogram

Env: needs pysam (NEOANTIGEN or ClusterCatcher). Not scanpy.
Usage:
  conda run -n NEOANTIGEN python Diagnostic_BRD4_Isoform_Feasibility_Probe.py [SRR_ID]

Author: Jake Lehle
"""

import os
import sys
import gzip
import json
import re
from collections import defaultdict, OrderedDict

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    import pysam
except ImportError:
    sys.exit("ERROR: pysam not found. Run in the NEOANTIGEN (or ClusterCatcher) env.")

# =============================================================================
# CONFIG
# =============================================================================
REF_DIR    = "/master/jlehle/WORKING/SC/ref/GRCh38"
GTF_PATH   = os.path.join(REF_DIR, "genes/genes_unzipped.gtf")
REF_JSON   = os.path.join(REF_DIR, "reference.json")

GSE_ROOT   = ("/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/"
              "results_NMF_v0.1.1/fastq/GSE173468")
SAMPLE     = sys.argv[1] if len(sys.argv) > 1 else "SRR14340883"

OUTPUT_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_6/DIAGNOSTIC_BRD4_ISOFORM"
os.makedirs(OUTPUT_DIR, exist_ok=True)

GENE_NAME  = "BRD4"
LOCUS_PAD  = 2000          # bp padding around the BRD4 exon span for BAM fetch
MAX_TAG_SCAN = 5000        # reads to scan for the CB/UB tag-presence report
EXPECTED_ENSEMBL = "115"   # neoantigen proteome build, for the cross-check note
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
# INTERVAL HELPERS  (all half-open, 0-based)
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
    """Return intervals in A not covered by B."""
    A = merge(A)
    B = merge(B)
    out = []
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

# =============================================================================
# GTF PARSING
# =============================================================================
def attr(field, key):
    m = re.search(key + r' "([^"]+)"', field)
    return m.group(1) if m else None

def parse_brd4(gtf_path, gene_name):
    """Return (seqname, strand, {tx_id: {'exons':[(s,e)0-based], 'name','biotype','ccds'}})."""
    txs = OrderedDict()
    seqname = strand = None
    tx_meta = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#') or '\t' not in line:
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            if ('gene_name "%s"' % gene_name) not in f[8]:
                continue
            feat = f[2]
            s0 = int(f[3]) - 1          # GTF 1-based inclusive -> 0-based half-open
            e0 = int(f[4])
            if seqname is None:
                seqname, strand = f[0], f[6]
            tid = attr(f[8], 'transcript_id')
            if feat == 'transcript' and tid:
                tx_meta[tid] = {
                    'name':    attr(f[8], 'transcript_name'),
                    'biotype': attr(f[8], 'transcript_biotype') or attr(f[8], 'transcript_type'),
                    'ccds':    attr(f[8], 'ccds_id') or attr(f[8], 'ccdsid'),
                }
            elif feat == 'exon' and tid:
                txs.setdefault(tid, []).append((s0, e0))
    # assemble
    out = OrderedDict()
    for tid, exons in txs.items():
        exons = merge(exons)
        meta = tx_meta.get(tid, {})
        out[tid] = {
            'exons': exons,
            'name': meta.get('name'),
            'biotype': meta.get('biotype'),
            'ccds': meta.get('ccds'),
            'lo': min(s for s, _ in exons),
            'hi': max(e for _, e in exons),
            'n_exons': len(exons),
            'span': max(e for _, e in exons) - min(s for s, _ in exons),
        }
    return seqname, strand, out

# =============================================================================
# STEP 0: reference build
# =============================================================================
banner(f"STEP 0: Reference build + sample ({SAMPLE})")
if os.path.exists(REF_JSON):
    with open(REF_JSON) as fh:
        rj = json.load(fh)
    for k in ('genomes', 'version', 'mkref_version', 'input_gtf_files',
              'input_fasta_files', 'assembly', 'annotation'):
        if k in rj:
            log(f"  reference.json[{k}] = {rj[k]}")
    log(f"  NOTE: neoantigen proteome build is Ensembl {EXPECTED_ENSEMBL}; if the "
        f"transcriptome annotation above differs, BRD4 transcript IDs may not match "
        f"the proteome and we anchor by CCDS/structure instead.")
else:
    log(f"  WARNING: {REF_JSON} not found; cannot report build.")

SAMPLE_OUTS = os.path.join(GSE_ROOT, SAMPLE, f"{SAMPLE}_S1_L001_", "outs")
BAM_PATH    = os.path.join(SAMPLE_OUTS, "possorted_genome_bam.bam")
WL_PATH     = os.path.join(SAMPLE_OUTS, "filtered_feature_bc_matrix", "barcodes.tsv.gz")
for p in (GTF_PATH, BAM_PATH, WL_PATH):
    log(f"  {'OK ' if os.path.exists(p) else 'MISSING'}  {p}")
if not (os.path.exists(BAM_PATH) and os.path.exists(WL_PATH)):
    sys.exit("ERROR: BAM or whitelist missing for this sample.")

# =============================================================================
# STEP 1: BRD4 transcripts + unique intervals
# =============================================================================
banner("STEP 1: BRD4 transcripts and isoform-unique intervals")
seqname, strand, txs = parse_brd4(GTF_PATH, GENE_NAME)
if not txs:
    sys.exit(f"ERROR: no BRD4 transcripts found in {GTF_PATH}. The filtered GTF may "
             f"not carry BRD4 isoforms; fall back to the full Ensembl GTF.")
log(f"  BRD4 on GTF seqname '{seqname}', strand '{strand}', {len(txs)} transcript(s)")
log(f"  minus-strand => 3' (poly-A) end is at the LOW coordinate\n")
log(f"  {'transcript_id':<20s} {'name':<14s} {'biotype':<16s} {'ccds':<12s} "
    f"{'exons':>5s} {'span':>7s} {'3p_end':>10s}")
for tid, t in txs.items():
    three_p = t['lo'] if strand == '-' else t['hi']
    log(f"  {tid:<20s} {str(t['name']):<14s} {str(t['biotype']):<16s} "
        f"{str(t['ccds']):<12s} {t['n_exons']:>5d} {t['span']:>7d} {three_p:>10d}")

# ---- CCDS-anchored isoform grouping (authoritative; replaces the structural guess) ----
# From the BRD4 review / UniProt:
#   BRD4-L    = CCDS12328 (O60885-1)  canonical long form
#   BRD4-S(a) = CCDS46004 (O60885-2)  the "BRD4-S" of interest
#   BRD4-S(b) = CCDS82307 (O60885-3)  rare / stress-induced
# Multiple transcript models can share one CCDS (same CDS, different UTRs); they MUST
# be grouped, or sibling transcripts subtract each other's unique intervals to zero.
CCDS_L, CCDS_SA, CCDS_SB = "CCDS12328", "CCDS46004", "CCDS82307"
groups = {'L': [], 'Sa': [], 'Sb': []}
background = []
for tid, t in txs.items():
    c = (t['ccds'] or '').split('.')[0]
    if   c == CCDS_L:  groups['L'].append(tid)
    elif c == CCDS_SA: groups['Sa'].append(tid)
    elif c == CCDS_SB: groups['Sb'].append(tid)
    else:              background.append(tid)

log("\n  CCDS-anchored grouping:")
for g, anchor in (('L', CCDS_L), ('Sa', CCDS_SA), ('Sb', CCDS_SB)):
    names = [f"{tid}({txs[tid]['name']})" for tid in groups[g]]
    log(f"    {g:<3s} [{anchor}]: {names if names else 'NONE FOUND'}")
log(f"    background (subtracted, not scored): {len(background)} transcript(s)")
for g, anchor in (('L', CCDS_L), ('Sa', CCDS_SA), ('Sb', CCDS_SB)):
    if not groups[g]:
        log(f"  WARNING: no transcript carries {anchor} for {g}; check GTF annotation.")

def group_exons(ids):
    ivs = []
    for tid in ids:
        ivs += txs[tid]['exons']
    return merge(ivs)

L_ex  = group_exons(groups['L'])
Sa_ex = group_exons(groups['Sa'])
Sb_ex = group_exons(groups['Sb'])
bg_ex = group_exons(background)

# each group's unique intervals = its exons minus the union of the other two CCDS
# groups PLUS all background transcripts, so a read there is isoform-specific
L_unique  = subtract(L_ex,  merge(Sa_ex + Sb_ex + bg_ex))
Sa_unique = subtract(Sa_ex, merge(L_ex + Sb_ex + bg_ex))
Sb_unique = subtract(Sb_ex, merge(L_ex + Sa_ex + bg_ex))

log(f"\n  unique interval sizes (bp):  L={total_bp(L_unique)}  "
    f"S(a)={total_bp(Sa_unique)}  S(b)={total_bp(Sb_unique)}")
if total_bp(L_unique) == 0:
    log("  WARNING: L-unique is empty; confirm both CCDS12328 models grouped together.")
if total_bp(Sa_unique) == 0:
    log("  WARNING: S(a)-unique is empty; L and S(a) not separable from this grouping. "
        "Inspect the transcript table above; may need the full Ensembl GTF.")

# locus fetch window
locus_lo = min(t['lo'] for t in txs.values()) - LOCUS_PAD
locus_hi = max(t['hi'] for t in txs.values()) + LOCUS_PAD

# =============================================================================
# STEP 2: BAM contig + tag presence
# =============================================================================
banner("STEP 2: BAM contig resolution + CB/UB tag check")
bam = pysam.AlignmentFile(BAM_PATH, "rb")
refs = set(bam.references)
contig = None
for cand in (seqname, f"chr{seqname}", seqname.replace("chr", "")):
    if cand in refs:
        contig = cand
        break
if contig is None:
    sys.exit(f"ERROR: could not match GTF seqname '{seqname}' to any BAM contig.")
log(f"  GTF seqname '{seqname}' -> BAM contig '{contig}'")
log(f"  fetch window: {contig}:{locus_lo}-{locus_hi}")

n_scan = n_cb = n_ub = 0
for r in bam.fetch(contig, locus_lo, locus_hi):
    if r.is_unmapped or r.is_secondary or r.is_supplementary:
        continue
    n_scan += 1
    if r.has_tag('CB'):
        n_cb += 1
    if r.has_tag('UB'):
        n_ub += 1
    if n_scan >= MAX_TAG_SCAN:
        break
if n_scan == 0:
    sys.exit("ERROR: no primary reads in the BRD4 locus for this sample.")
log(f"  scanned {n_scan} primary reads: CB present {100*n_cb/n_scan:.1f}%, "
    f"UB present {100*n_ub/n_scan:.1f}%")
if n_cb == 0 or n_ub == 0:
    sys.exit("ERROR: CB or UB tags absent; cannot do UMI-level counting on this BAM.")

# =============================================================================
# STEP 3: classify molecules by unique-interval overlap; dedup by (CB, UB)
# =============================================================================
banner("STEP 3: Discriminating-molecule counting")
whitelist = set()
with gzip.open(WL_PATH, 'rt') as fh:
    for line in fh:
        whitelist.add(line.strip())
log(f"  whitelist barcodes: {len(whitelist)}")

# per (CB,UB) -> set of classes its reads overlap
mol_class = defaultdict(set)
n_reads_used = 0            # discriminating reads passing the whitelist
n_disc_total = 0           # discriminating reads with CB+UB, pre-whitelist
n_cb_not_in_wl = 0         # discriminating reads whose CB is absent from the whitelist
for r in bam.fetch(contig, locus_lo, locus_hi):
    if r.is_unmapped or r.is_secondary or r.is_supplementary:
        continue
    if not (r.has_tag('CB') and r.has_tag('UB')):
        continue
    blocks = r.get_blocks()          # 0-based half-open reference blocks
    cls = set()
    for b0, b1 in blocks:
        if overlaps_any(b0, b1, L_unique):
            cls.add('L')
        if Sa_unique and overlaps_any(b0, b1, Sa_unique):
            cls.add('Sa')
        if Sb_unique and overlaps_any(b0, b1, Sb_unique):
            cls.add('Sb')
    if not cls:
        continue
    n_disc_total += 1
    cb = r.get_tag('CB')
    if cb not in whitelist:
        n_cb_not_in_wl += 1
        continue
    ub = r.get_tag('UB')
    mol_class[(cb, ub)] |= cls
    n_reads_used += 1

if n_disc_total > 0 and n_reads_used == 0:
    log(f"  WARNING: {n_disc_total} discriminating reads had CB+UB but 0 passed the "
        f"whitelist ({n_cb_not_in_wl} CBs absent). Likely a barcode-format mismatch "
        f"between the BAM CB tag and barcodes.tsv.gz, NOT a depth problem. Check formats.")
elif n_disc_total > 0:
    log(f"  discriminating reads: {n_disc_total} total, {n_reads_used} in whitelist "
        f"({100*n_cb_not_in_wl/n_disc_total:.1f}% dropped as non-whitelist barcodes)")

# resolve each molecule
pb = {'L': 0, 'Sa': 0, 'Sb': 0, 'ambiguous': 0}
percell = defaultdict(int)
for (cb, ub), cls in mol_class.items():
    if len(cls) == 1:
        c = next(iter(cls))
        pb[c] += 1
        percell[cb] += 1
    else:
        pb['ambiguous'] += 1

log(f"  reads overlapping a unique interval: {n_reads_used}")
log(f"  discriminating molecules (unique-class UMIs):")
log(f"    L    = {pb['L']}")
log(f"    S(a) = {pb['Sa']}")
log(f"    S(b) = {pb['Sb']}")
log(f"    ambiguous (multi-class, excluded) = {pb['ambiguous']}")
denom = pb['L'] + pb['Sa'] + pb['Sb']
if denom > 0:
    log(f"  pseudobulk ratio  L : S(a) : S(b)  =  "
        f"{pb['L']/denom:.3f} : {pb['Sa']/denom:.3f} : {pb['Sb']/denom:.3f}")
    if pb['Sa'] + pb['Sb'] > 0:
        log(f"  L / S total = {pb['L'] / (pb['Sa'] + pb['Sb']):.2f}")
else:
    log("  no discriminating molecules recovered.")

# =============================================================================
# STEP 4: per-cell depth + feasibility verdict
# =============================================================================
banner("STEP 4: Per-cell depth + feasibility verdict")
counts = np.array(list(percell.values())) if percell else np.array([])
n_wl = len(whitelist)
def frac_ge(n):
    return (np.sum(counts >= n) / n_wl) if n_wl else 0.0
if counts.size:
    log(f"  cells with >=1 discriminating UMI: {counts.size} "
        f"({100*counts.size/n_wl:.1f}% of {n_wl} whitelist cells)")
    log(f"  among those cells: median={np.median(counts):.0f}, "
        f"90th pct={np.percentile(counts,90):.0f}, max={counts.max():.0f}")
    log(f"  fraction of ALL whitelist cells with >=1 / >=3 / >=5 disc. UMIs: "
        f"{100*frac_ge(1):.1f}% / {100*frac_ge(3):.1f}% / {100*frac_ge(5):.1f}%")
else:
    log("  no cells with discriminating UMIs.")

log("")
log("  FEASIBILITY VERDICT (stage 1, single sample):")
if denom == 0 or total_bp(Sa_unique) == 0:
    log("    - Isoforms NOT separable here (empty S(a)-unique or no molecules). "
        "Revisit the transcript list / GTF before proceeding.")
elif frac_ge(3) >= 0.20:
    log("    - Per-cell L/S looks PLAUSIBLE (>=20% of cells carry >=3 disc. UMIs). "
        "Proceed to stage 2: population-stratified per-cell ratio across samples.")
elif denom >= 50:
    log("    - Per-cell too sparse, but pseudobulk is solid. Report population-level "
        "L : S(a) only; per-cell resolution is not supported by depth.")
else:
    log("    - Even pseudobulk is thin on this sample. Re-run on a higher-BRD4 sample "
        "before deciding; long-read (MinION) is the clean route if scRNA can't resolve it.")

# histogram
if counts.size:
    fig, ax = plt.subplots(figsize=(7, 5))
    mx = int(counts.max())
    ax.hist(counts, bins=range(1, max(mx + 2, 3)), color='#7a4fa3',
            edgecolor='white', align='left')
    ax.set_xlabel('discriminating UMIs per cell', fontsize=12)
    ax.set_ylabel('cells', fontsize=12)
    ax.set_title(f'BRD4 isoform-discriminating UMIs per cell\n{SAMPLE} '
                 f'(n={counts.size} cells with >=1)', fontsize=11)
    ax.tick_params(labelsize=10)
    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUTPUT_DIR, f"brd4_percell_umi_hist_{SAMPLE}.{ext}"),
                    dpi=DPI, bbox_inches='tight')
    plt.close(fig)
    log(f"\n  [SAVE] brd4_percell_umi_hist_{SAMPLE}.pdf/.png")

bam.close()

# =============================================================================
# WRITE REPORT
# =============================================================================
rp = os.path.join(OUTPUT_DIR, f"brd4_isoform_feasibility_{SAMPLE}.txt")
with open(rp, 'w') as fh:
    fh.write("\n".join(report))
print(f"\n[SAVE] {rp}")
print("BRD4 ISOFORM FEASIBILITY PROBE COMPLETE")
