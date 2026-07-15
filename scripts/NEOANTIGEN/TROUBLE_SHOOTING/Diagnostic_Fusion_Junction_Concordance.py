#!/usr/bin/env python3
"""
Diagnostic_Fusion_Junction_Concordance.py
============================================
Prove the from-scratch STAR chimeric regeneration reproduces the original run.

Compares the freshly generated Chimeric.out.junction files (built by Step05a
into data/FIG_7/04_fusion_analysis/star_chimeric/) against the pre-split backup
set, per sample, at the JUNCTION-SET level.

Why set-level, not byte-level: a from-scratch index is never byte-identical to
the backed-up one, and multithreaded STAR emits junction lines in
nondeterministic order. So byte diff is meaningless. What must be reproducible
is the SET of chimeric junctions called. For each sample we build the set of
junction keys (donor chrom/pos/strand, acceptor chrom/pos/strand, junction type)
and report Jaccard overlap plus the fraction of the backup's junctions that the
new run recovers.

Interpretation:
  - High overlap  -> our from-scratch process reproduces the original result;
                     we own the process, not just a backed-up file.
  - Low overlap   -> something drifted (STAR version, index, params); catch it
                     here, before it reaches a reviewer.

This is a diagnostic only. It reads both junction trees and writes a report;
it changes nothing in the pipeline.

Env: NETWORK or NEOANTIGEN (pandas only; no STAR/pysam needed).
Usage: python Diagnostic_Fusion_Junction_Concordance.py
Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import glob
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"

NEW_DIR = os.path.join(PROJECT_ROOT,
    "data/FIG_7/04_fusion_analysis/star_chimeric")
BACKUP_DIR = os.path.join(PROJECT_ROOT,
    "data/_BACKUP_FIG6_full_20260625_105330/FIG_6/05_neoantigen/"
    "fusion_detection/star_chimeric")

OUTPUT_DIR = os.path.join(PROJECT_ROOT,
    "data/FIG_7/TROUBLESHOOTING/fusion_junction_concordance")

# A sample is flagged if the fraction of backup junctions recovered falls below
# this. Set intentionally high; this is a reproducibility check, not a fuzzy one.
RECOVERY_FLAG_THRESHOLD = 0.95

report_lines = []

def log(msg=""):
    print(msg, flush=True)
    report_lines.append(str(msg))

def banner(title, char="="):
    log("")
    log(char * 80)
    log(f"  {title}")
    log(char * 80)


def load_junction_keys(path):
    """Return a set of junction keys from a Chimeric.out.junction file.

    Key = (chr_donorA, brkpt_donorA, strand_donorA,
           chr_acceptorB, brkpt_acceptorB, strand_acceptorB, junction_type)
    Skips the header line and the chimOutJunctionFormat=1 trailing '#' summary.
    """
    keys = set()
    with open(path) as f:
        for line in f:
            if line.startswith('#') or line.startswith('chr_donorA'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 7:
                continue
            keys.add(tuple(fields[0:7]))
    return keys


def sample_ids(root):
    ids = set()
    for d in glob.glob(os.path.join(root, "SRR*")):
        if os.path.isfile(os.path.join(d, "Chimeric.out.junction")):
            ids.add(os.path.basename(d))
    return ids


def main():
    banner("Diagnostic: STAR chimeric junction concordance (new vs backup)")
    log(f"  {datetime.now().isoformat(timespec='seconds')}")
    log(f"  New:    {NEW_DIR}")
    log(f"  Backup: {BACKUP_DIR}")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    if not os.path.isdir(NEW_DIR):
        log(f"  [FATAL] New junction dir not found: {NEW_DIR}")
        log( "          Run Step05a first.")
        return
    if not os.path.isdir(BACKUP_DIR):
        log(f"  [FATAL] Backup junction dir not found: {BACKUP_DIR}")
        return

    new_ids = sample_ids(NEW_DIR)
    bak_ids = sample_ids(BACKUP_DIR)
    shared = sorted(new_ids & bak_ids)

    log(f"\n  Samples with junctions: new={len(new_ids)}, backup={len(bak_ids)}, "
        f"comparable={len(shared)}")
    only_new = sorted(new_ids - bak_ids)
    only_bak = sorted(bak_ids - new_ids)
    if only_new:
        log(f"  Only in NEW (no backup to compare): {', '.join(only_new)}")
    if only_bak:
        log(f"  Only in BACKUP (not regenerated):   {', '.join(only_bak)}")

    if not shared:
        log("\n  [FATAL] No overlapping samples to compare.")
        return

    banner("Per-sample concordance", char="-")
    log(f"  {'sample':16s} {'new':>8s} {'backup':>8s} {'shared':>8s} "
        f"{'Jaccard':>8s} {'bkp_recov':>9s}  flag")
    log(f"  {'-'*16} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*9}  ----")

    rows = []
    flagged = []
    for srr in shared:
        new_keys = load_junction_keys(os.path.join(NEW_DIR, srr, "Chimeric.out.junction"))
        bak_keys = load_junction_keys(os.path.join(BACKUP_DIR, srr, "Chimeric.out.junction"))
        inter = new_keys & bak_keys
        union = new_keys | bak_keys
        jacc = len(inter) / len(union) if union else 1.0
        recov = len(inter) / len(bak_keys) if bak_keys else 1.0
        flag = "" if recov >= RECOVERY_FLAG_THRESHOLD else "  <-- LOW"
        if recov < RECOVERY_FLAG_THRESHOLD:
            flagged.append(srr)
        log(f"  {srr:16s} {len(new_keys):8d} {len(bak_keys):8d} {len(inter):8d} "
            f"{jacc:8.4f} {recov:9.4f}{flag}")
        rows.append({'sample': srr, 'n_new': len(new_keys), 'n_backup': len(bak_keys),
                     'n_shared': len(inter), 'jaccard': round(jacc, 4),
                     'backup_recovery': round(recov, 4)})

    # -- Summary -------------------------------------------------------------
    banner("Summary", char="-")
    mean_jacc = sum(r['jaccard'] for r in rows) / len(rows)
    mean_recov = sum(r['backup_recovery'] for r in rows) / len(rows)
    tot_new = sum(r['n_new'] for r in rows)
    tot_bak = sum(r['n_backup'] for r in rows)
    tot_shared = sum(r['n_shared'] for r in rows)
    log(f"  Samples compared:            {len(rows)}")
    log(f"  Mean Jaccard:                {mean_jacc:.4f}")
    log(f"  Mean backup recovery:        {mean_recov:.4f}")
    log(f"  Pooled backup recovery:      {tot_shared/tot_bak if tot_bak else 1.0:.4f} "
        f"({tot_shared:,}/{tot_bak:,})")
    log(f"  Pooled new total junctions:  {tot_new:,}")
    if flagged:
        log(f"\n  [WARNING] {len(flagged)} sample(s) below {RECOVERY_FLAG_THRESHOLD:.2f} "
            f"recovery: {', '.join(flagged)}")
        log( "            Investigate before trusting regenerated fusion numbers.")
    else:
        log(f"\n  All samples >= {RECOVERY_FLAG_THRESHOLD:.2f} backup recovery. "
            f"Regeneration reproduces the original junction set.")

    # -- Write report + table ------------------------------------------------
    try:
        import pandas as pd
        pd.DataFrame(rows).to_csv(
            os.path.join(OUTPUT_DIR, "junction_concordance_by_sample.tsv"),
            sep='\t', index=False)
    except Exception:
        pass
    with open(os.path.join(OUTPUT_DIR, "fusion_junction_concordance_report.txt"), 'w') as fh:
        fh.write('\n'.join(report_lines))
    log(f"\n  Wrote: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
