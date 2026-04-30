#!/usr/bin/env python3
"""
TROUBLESHOOT_KEGG_Enrichment.py
================================
Standalone KEGG enrichment for all communities with aggressive rate
limiting and verbose diagnostics. Does NOT depend on network_config.

Reads community gene lists, runs Enrichr one community at a time with
long pauses between requests, and saves results to the pipeline summary
directory.

Diagnostic goals:
  1. Confirm gseapy exception type and message format on 429
  2. Test whether longer inter-request delays prevent 429s
  3. Fill in all missing KEGG results in one clean run

Usage:
    conda run -n NETWORK python TROUBLESHOOT_KEGG_Enrichment.py

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os, sys, time, traceback
import pandas as pd

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
CANCER_TYPE = "TCGA-HNSC"

COMMUNITY_GENE_LISTS = os.path.join(
    PROJECT_ROOT, "data/FIG_2/05_communities", CANCER_TYPE,
    f"{CANCER_TYPE}_community_gene_lists.csv"
)
REPORT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_2/08_pipeline_summary")
KEGG_SUMMARY_PATH = os.path.join(REPORT_DIR, f"{CANCER_TYPE}_KEGG_enrichment_summary.csv")

# Aggressive rate limiting
BETWEEN_DELAY = 15      # seconds between successful requests
MAX_RETRIES = 5
RETRY_DELAYS = [30, 60, 120, 240, 300]  # fixed schedule, not exponential

# =============================================================================
# HELPERS
# =============================================================================
def log(msg=""):
    print(msg, flush=True)


def run_enrichr_verbose(gene_list, community_id, max_retries=MAX_RETRIES):
    """
    Run Enrichr with verbose diagnostics on every exception.
    Returns (results_df_or_None, error_msg_or_None)
    """
    import gseapy as gp

    for attempt in range(max_retries):
        log(f"    Attempt {attempt+1}/{max_retries}...")

        try:
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=["KEGG_2021_Human"],
                organism="human",
                outdir=None,
                no_plot=True,
                verbose=False,
            )
            log(f"    SUCCESS on attempt {attempt+1}")
            return enr.results, None

        except Exception as e:
            err_type = type(e).__name__
            err_str = str(e)
            log(f"    EXCEPTION on attempt {attempt+1}:")
            log(f"      Type: {err_type}")
            log(f"      Message: {err_str}")
            log(f"      Has '429' in message: {'429' in err_str}")
            log(f"      Has '429' in repr: {'429' in repr(e)}")

            # Print full traceback for first failure of each community
            if attempt == 0:
                log(f"      Full traceback:")
                tb_lines = traceback.format_exception(type(e), e, e.__traceback__)
                for tb_line in tb_lines:
                    for sub_line in tb_line.rstrip().split('\n'):
                        log(f"        {sub_line}")

            # Check if it looks like a rate limit (broad matching)
            is_rate_limit = any(marker in err_str.lower() for marker in [
                '429', 'rate limit', 'too many requests', 'retry'
            ])

            if is_rate_limit and attempt < max_retries - 1:
                delay = RETRY_DELAYS[min(attempt, len(RETRY_DELAYS) - 1)]
                log(f"    Rate limit detected, waiting {delay}s before retry...")
                time.sleep(delay)
            elif not is_rate_limit:
                log(f"    Non-rate-limit error, not retrying")
                return None, f"{err_type}: {err_str}"
            else:
                log(f"    Final retry exhausted")
                return None, f"Rate limited after {max_retries} retries: {err_str}"

    return None, f"Exhausted {max_retries} retries"


# =============================================================================
# MAIN
# =============================================================================
log("=" * 70)
log("  TROUBLESHOOT KEGG Enrichment")
log("=" * 70)

# Load community gene lists
comm_df = pd.read_csv(COMMUNITY_GENE_LISTS)
community_genes = {}
for _, row in comm_df.iterrows():
    cid = int(row['community'])
    genes = [g.strip() for g in str(row['genes']).split(';') if g.strip()]
    community_genes[cid] = genes

log(f"\n  Communities: {len(community_genes)}")
for cid in sorted(community_genes.keys()):
    log(f"    C{cid}: {len(community_genes[cid])} genes")

# Check gseapy
log(f"\n  Checking gseapy...")
try:
    import gseapy as gp
    log(f"  gseapy version: {gp.__version__}")
except ImportError:
    log("  FATAL: gseapy not installed")
    sys.exit(1)

# Quick connectivity test with a small known gene list
log(f"\n  Testing Enrichr connectivity with a small gene list...")
test_genes = ["TP53", "BRCA1", "EGFR"]
test_result, test_err = run_enrichr_verbose(test_genes, community_id=-1, max_retries=2)
if test_result is not None:
    log(f"  Connectivity test PASSED ({len(test_result)} results)")
else:
    log(f"  Connectivity test FAILED: {test_err}")
    log(f"  Enrichr may be down or blocked from this network")
    log(f"  Continuing anyway in case it's transient...")

log(f"\n  Waiting {BETWEEN_DELAY}s after connectivity test...")
time.sleep(BETWEEN_DELAY)

# Run enrichment for each community
log(f"\n{'='*70}")
log(f"  Running KEGG enrichment for all {len(community_genes)} communities")
log(f"  Inter-request delay: {BETWEEN_DELAY}s")
log(f"  Max retries per community: {MAX_RETRIES}")
log(f"  Retry delays: {RETRY_DELAYS}")
log(f"{'='*70}")

all_results = []
sorted_comms = sorted(community_genes.keys())

for i, cid in enumerate(sorted_comms):
    genes = community_genes[cid]
    n_genes = len(genes)

    log(f"\n  [{i+1}/{len(sorted_comms)}] Community {cid} ({n_genes} genes)")
    log(f"    First 5 genes: {', '.join(genes[:5])}")

    if n_genes < 3:
        log(f"    SKIP: too few genes")
        all_results.append({
            'community': cid,
            'n_genes': n_genes,
            'n_mapped': n_genes,
            'top_kegg_term': 'N/A (too few genes)',
            'top_kegg_pvalue': None,
            'top_kegg_genes': None,
            'n_kegg_significant': 0,
            'status': 'skipped',
        })
        continue

    results_df, err = run_enrichr_verbose(genes, cid)

    if results_df is not None and len(results_df) > 0:
        sig = results_df[results_df['Adjusted P-value'] < 0.05]
        n_sig = len(sig)
        top = results_df.iloc[0]

        log(f"    Top KEGG: {top['Term']} (adj p={top['Adjusted P-value']:.2e})")
        log(f"    Significant terms (adj p < 0.05): {n_sig}")

        if n_sig > 0:
            for _, r in sig.head(5).iterrows():
                log(f"      {r['Term']} (adj p={r['Adjusted P-value']:.2e}, "
                    f"genes: {r['Genes']})")

        all_results.append({
            'community': cid,
            'n_genes': n_genes,
            'n_mapped': n_genes,
            'top_kegg_term': top['Term'],
            'top_kegg_pvalue': top['Adjusted P-value'],
            'top_kegg_genes': top['Genes'],
            'n_kegg_significant': n_sig,
            'status': 'success',
        })

        # Save per-community CSV
        csv_path = os.path.join(REPORT_DIR, f"{CANCER_TYPE}_community_{cid:02d}_KEGG.csv")
        results_df.to_csv(csv_path, index=False)
        log(f"    Saved: {os.path.basename(csv_path)}")

    elif results_df is not None and len(results_df) == 0:
        log(f"    No KEGG results returned (empty DataFrame)")
        all_results.append({
            'community': cid,
            'n_genes': n_genes,
            'n_mapped': n_genes,
            'top_kegg_term': 'No results',
            'top_kegg_pvalue': None,
            'top_kegg_genes': None,
            'n_kegg_significant': 0,
            'status': 'empty',
        })

    else:
        log(f"    FAILED: {err}")
        all_results.append({
            'community': cid,
            'n_genes': n_genes,
            'n_mapped': n_genes,
            'top_kegg_term': f'ERROR: {str(err)[:80]}',
            'top_kegg_pvalue': None,
            'top_kegg_genes': None,
            'n_kegg_significant': 0,
            'status': 'error',
        })

    # Wait between communities (skip after last)
    if i < len(sorted_comms) - 1:
        log(f"    Waiting {BETWEEN_DELAY}s before next community...")
        time.sleep(BETWEEN_DELAY)


# =============================================================================
# SUMMARY AND SAVE
# =============================================================================
log(f"\n{'='*70}")
log(f"  RESULTS SUMMARY")
log(f"{'='*70}")

results_df_all = pd.DataFrame(all_results)

n_success = (results_df_all['status'] == 'success').sum()
n_empty = (results_df_all['status'] == 'empty').sum()
n_error = (results_df_all['status'] == 'error').sum()
n_skip = (results_df_all['status'] == 'skipped').sum()

log(f"  Success: {n_success}")
log(f"  Empty (no results): {n_empty}")
log(f"  Error (failed): {n_error}")
log(f"  Skipped (too few genes): {n_skip}")

if n_error > 0:
    log(f"\n  Failed communities:")
    for _, row in results_df_all[results_df_all['status'] == 'error'].iterrows():
        log(f"    C{int(row['community'])}: {row['top_kegg_term']}")

# Save as the canonical KEGG summary (overwrite)
save_df = results_df_all.drop(columns=['status'])
save_df.to_csv(KEGG_SUMMARY_PATH, index=False)
log(f"\n  Saved KEGG summary: {KEGG_SUMMARY_PATH}")

# Also save a diagnostic copy with status column
diag_path = os.path.join(
    PROJECT_ROOT, "data/FIG_2/DIAGNOSTIC_AUDIT",
    f"{CANCER_TYPE}_KEGG_troubleshoot_results.tsv"
)
os.makedirs(os.path.dirname(diag_path), exist_ok=True)
results_df_all.to_csv(diag_path, sep='\t', index=False)
log(f"  Saved diagnostic: {diag_path}")

if n_error == 0:
    log(f"\n  ALL COMMUNITIES COMPLETE")
else:
    log(f"\n  {n_error} communities still failed.")
    log(f"  Check the traceback output above for exception details.")
    log(f"  If Enrichr is consistently rate-limiting, consider:")
    log(f"    1. Increasing BETWEEN_DELAY (currently {BETWEEN_DELAY}s)")
    log(f"    2. Running at a different time of day (lower Enrichr traffic)")
    log(f"    3. Using a local KEGG database instead of the Enrichr API")

log(f"\n{'='*70}")
log(f"  DONE")
log(f"{'='*70}")
