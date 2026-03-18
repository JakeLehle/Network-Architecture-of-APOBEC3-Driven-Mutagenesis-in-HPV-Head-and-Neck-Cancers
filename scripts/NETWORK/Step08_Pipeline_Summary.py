#!/usr/bin/env python3
"""
Step08_Pipeline_Summary_and_Enrichment.py

Gather all pipeline parameters, per-step statistics, and community
gene lists. Run KEGG pathway enrichment on each community using
gseapy's enrichr interface. Output a comprehensive summary report
for writing the results and methods sections.

Run AFTER the full pipeline (Steps 01–07) has completed.

Dependencies:
    pip install gseapy    (for KEGG enrichment)

Usage:
    conda run -n NETWORK python Step08_Pipeline_Summary_and_Enrichment.py
"""

import os
import json
import pickle
import numpy as np
import pandas as pd
import networkx as nx
from datetime import datetime

from network_config import (
    # Directories
    FIG2_ROOT,
    DIR_01_CLEANED, DIR_02_MERGED, DIR_03_DIFFEXPR,
    DIR_04_NETWORKS, DIR_05_COMMUNITIES, DIR_06_CENTRALITY,
    # Cancer types
    CANCER_TYPES,
    # Clinical
    CLINICAL_COLS, N_CLINICAL,
    # A3 genes
    A3_GENES, A3_ID_TO_ALIAS, A3_ALIAS_TO_ID,
    # Biomarkers
    BIOMARKERS,
    # Step 03 params
    MIN_SAMPLES_DETECTED,
    A3_SUM_PERCENTILE, SBS2_HIGH_PERCENTILE, SBS2_LOW_PERCENTILE,
    RAW_P_THRESHOLD, LOGFC_THRESHOLD, FORCE_KEEP_A3,
    # Step 04 params
    CORRELATION_METHOD, CORR_THRESHOLD, DIFF_THRESHOLD,
    MIN_GROUP_SIZE,
    # Step 05 params
    COMMUNITY_METHOD, COMMUNITY_RESOLUTIONS, RUNS_PER_RESOLUTION,
    COMMUNITY_BASE_SEED, USE_LARGEST_COMPONENT,
    TARGET_BIG_COMMUNITIES, MIN_COMMUNITY_SIZE,
    # Utilities
    banner, log, ensure_dir
)


# =============================================================================
# OUTPUT DIRECTORY
# =============================================================================
REPORT_DIR = ensure_dir(os.path.join(FIG2_ROOT, "08_pipeline_summary"))
TIMESTAMP = datetime.now().strftime("%Y-%m-%d_%H%M%S")

# Symbol mapping
with open(os.path.join(DIR_01_CLEANED, "ensg_to_symbol.json")) as f:
    ensg_to_symbol = json.load(f)

def sym(ensg):
    return ensg_to_symbol.get(str(ensg), str(ensg))


# =============================================================================
# SECTION 1 — Pipeline Configuration Parameters
# =============================================================================
banner("[SECTION 1] Pipeline Configuration Parameters")

config_lines = []

def add(label, value, section=None):
    """Append a key-value pair to the config report."""
    line = f"  {label}: {value}"
    config_lines.append(line)
    log(line)

log("--- Step 03: Group Definition & Differential Expression ---")
add("A3 sum percentile (A3A+A3B)", f"{A3_SUM_PERCENTILE} (median = top {int((1-A3_SUM_PERCENTILE)*100)}%)")
add("SBS2 HIGH percentile", f"top {int((1-SBS2_HIGH_PERCENTILE)*100)}% within high-A3 group")
add("SBS2 LOW percentile", f"bottom {int(SBS2_LOW_PERCENTILE*100)}% within high-A3 group")
add("Min samples detected per gene", MIN_SAMPLES_DETECTED)
add("Differential test", "Wilcoxon rank-sum on log1p(FPKM-UQ)")
add("Raw p-value threshold", RAW_P_THRESHOLD)
add("Log fold-change threshold", LOGFC_THRESHOLD)
add("Force-keep A3 genes", FORCE_KEEP_A3)

log("\n--- Step 04: Correlation Network Construction ---")
add("Correlation method", CORRELATION_METHOD)
add("TOP/BOTTOM |rho| threshold", CORR_THRESHOLD)
add("DIFF |delta-rho| threshold", DIFF_THRESHOLD)

log("\n--- Step 05: Community Detection ---")
add("Community method", COMMUNITY_METHOD)
add("Resolutions tested", COMMUNITY_RESOLUTIONS)
add("Runs per resolution", RUNS_PER_RESOLUTION)
add("Base seed", COMMUNITY_BASE_SEED)
add("Use largest component", USE_LARGEST_COMPONENT)
add("Target big communities", TARGET_BIG_COMMUNITIES)
add("Min community size", MIN_COMMUNITY_SIZE)


# =============================================================================
# SECTION 2 — Per-Cancer-Type Pipeline Statistics
# =============================================================================

for cancer_type in CANCER_TYPES:

    banner(f"[SECTION 2] Pipeline Statistics — {cancer_type}")

    report = {}
    report["cancer_type"] = cancer_type
    report["timestamp"] = TIMESTAMP

    # ---- Step 03: Group sizes and A3/SBS2 thresholds ----
    log("\n--- Step 03: Group Definitions ---")

    # Load HIGH and LOW group data
    high_pkl = os.path.join(DIR_03_DIFFEXPR, cancer_type,
                            f"{cancer_type}_SBS2_HIGH_group.pkl")
    low_pkl = os.path.join(DIR_03_DIFFEXPR, cancer_type,
                           f"{cancer_type}_SBS2_LOW_group.pkl")

    if os.path.exists(high_pkl) and os.path.exists(low_pkl):
        with open(high_pkl, "rb") as f:
            df_high = pickle.load(f)
        with open(low_pkl, "rb") as f:
            df_low = pickle.load(f)

        n_high = df_high.shape[0]
        n_low = df_low.shape[0]

        report["n_SBS2_HIGH"] = n_high
        report["n_SBS2_LOW"] = n_low

        log(f"  SBS2 HIGH group: {n_high} tumors")
        log(f"  SBS2 LOW group:  {n_low} tumors")

        # Extract A3 and SBS2 stats from groups
        for alias, ensg in A3_ALIAS_TO_ID.items():
            if alias in ["A3A", "A3B"]:
                if ensg in df_high.columns:
                    high_vals = pd.to_numeric(df_high[ensg], errors="coerce")
                    low_vals = pd.to_numeric(df_low[ensg], errors="coerce")
                    report[f"{alias}_median_HIGH"] = round(float(high_vals.median()), 4)
                    report[f"{alias}_median_LOW"] = round(float(low_vals.median()), 4)
                    log(f"  {alias} median — HIGH: {high_vals.median():.4f}, LOW: {low_vals.median():.4f}")
                elif alias in df_high.columns:
                    high_vals = pd.to_numeric(df_high[alias], errors="coerce")
                    low_vals = pd.to_numeric(df_low[alias], errors="coerce")
                    report[f"{alias}_median_HIGH"] = round(float(high_vals.median()), 4)
                    report[f"{alias}_median_LOW"] = round(float(low_vals.median()), 4)
                    log(f"  {alias} median — HIGH: {high_vals.median():.4f}, LOW: {low_vals.median():.4f}")

        # SBS2 stats
        if "SBS2" in df_high.columns:
            report["SBS2_median_HIGH"] = round(float(df_high["SBS2"].median()), 6)
            report["SBS2_median_LOW"] = round(float(df_low["SBS2"].median()), 6)
            report["SBS2_range_HIGH"] = f"{df_high['SBS2'].min():.6f} – {df_high['SBS2'].max():.6f}"
            report["SBS2_range_LOW"] = f"{df_low['SBS2'].min():.6f} – {df_low['SBS2'].max():.6f}"
            log(f"  SBS2 median — HIGH: {df_high['SBS2'].median():.6f}, LOW: {df_low['SBS2'].median():.6f}")
    else:
        log("  [WARNING] Group pickle files not found")

    # ---- Step 03: Differential expression gene counts ----
    log("\n--- Step 03: Differential Expression ---")

    stats_csv = os.path.join(DIR_03_DIFFEXPR, cancer_type,
                             f"{cancer_type}_diffexpr_stats.csv")
    sel_csv = os.path.join(DIR_03_DIFFEXPR, cancer_type,
                           f"{cancer_type}_selected_genes.csv")
    sel_filtered_csv = os.path.join(DIR_03_DIFFEXPR, cancer_type,
                                    f"{cancer_type}_selected_genes_filtered.csv")

    if os.path.exists(stats_csv):
        stats_df = pd.read_csv(stats_csv)
        report["n_genes_tested"] = len(stats_df)
        report["n_genes_significant_raw"] = int((stats_df["p_value"] < RAW_P_THRESHOLD).sum())
        log(f"  Genes tested: {len(stats_df)}")
        log(f"  Significant (raw p < {RAW_P_THRESHOLD}): {report['n_genes_significant_raw']}")

    if os.path.exists(sel_csv):
        sel_df = pd.read_csv(sel_csv)
        report["n_genes_selected"] = len(sel_df)
        log(f"  Genes selected (before filtering): {len(sel_df)}")

    if os.path.exists(sel_filtered_csv):
        sel_filt_df = pd.read_csv(sel_filtered_csv)
        report["n_genes_selected_filtered"] = len(sel_filt_df)
        log(f"  Genes selected (after network-readiness filter): {len(sel_filt_df)}")

    # ---- Step 04: Correlation network stats ----
    log("\n--- Step 04: Correlation Networks ---")

    for graph_label in ["G_top_noiso", "G_bot_noiso", "G_diff_noiso"]:
        gpickle_path = os.path.join(DIR_04_NETWORKS, cancer_type,
                                    "graph_objects",
                                    f"{cancer_type}_{graph_label}.gpickle")
        if os.path.exists(gpickle_path):
            with open(gpickle_path, "rb") as f:
                G = pickle.load(f)
            n_nodes = G.number_of_nodes()
            n_edges = G.number_of_edges()
            avg_deg = 2 * n_edges / n_nodes if n_nodes > 0 else 0
            report[f"{graph_label}_nodes"] = n_nodes
            report[f"{graph_label}_edges"] = n_edges
            report[f"{graph_label}_avg_degree"] = round(avg_deg, 2)
            log(f"  {graph_label}: {n_nodes} nodes, {n_edges} edges, avg degree {avg_deg:.1f}")

    # ---- Step 05: Community detection ----
    log("\n--- Step 05: Community Detection ---")

    part_csv = os.path.join(DIR_05_COMMUNITIES, cancer_type,
                            f"{cancer_type}_best_partition.csv")
    sweep_csv = os.path.join(DIR_05_COMMUNITIES, cancer_type,
                             f"{cancer_type}_resolution_sweep.csv")

    partition_df = None
    if os.path.exists(part_csv):
        partition_df = pd.read_csv(part_csv)
        n_communities = partition_df["community"].nunique()
        report["n_communities"] = n_communities
        report["n_genes_in_communities"] = len(partition_df)
        log(f"  Communities: {n_communities}")
        log(f"  Genes in communities: {len(partition_df)}")

        # Per-community sizes
        comm_sizes = partition_df.groupby("community").size().sort_values(ascending=False)
        for c, size in comm_sizes.items():
            report[f"community_{c}_size"] = int(size)
            log(f"    Community {c}: {size} genes")

    if os.path.exists(sweep_csv):
        sweep_df = pd.read_csv(sweep_csv)
        # Find best resolution (highest modularity or user-selected)
        if "modularity_mean" in sweep_df.columns:
            best_row = sweep_df.loc[sweep_df["modularity_mean"].idxmax()]
            report["best_resolution"] = float(best_row["resolution"])
            report["best_modularity"] = round(float(best_row["modularity_mean"]), 4)
            report["best_ARI"] = round(float(best_row.get("ARI_mean", 0)), 4)
            log(f"  Best resolution: {report['best_resolution']} "
                f"(modularity={report['best_modularity']}, ARI={report['best_ARI']})")

    # ---- Check biomarker presence in communities ----
    log("\n--- Biomarker Presence in Communities ---")

    biomarker_presence = {}
    if partition_df is not None:
        comm_genes_set = set(partition_df["gene"])
        for ensg in BIOMARKERS:
            symbol = sym(ensg)
            present = ensg in comm_genes_set
            if present:
                comm = int(partition_df[partition_df["gene"] == ensg]["community"].iloc[0])
                biomarker_presence[symbol] = f"Community {comm}"
            else:
                biomarker_presence[symbol] = "NOT IN ANY COMMUNITY"
            log(f"  {symbol} ({ensg}): {biomarker_presence[symbol]}")

        n_present = sum(1 for v in biomarker_presence.values() if "Community" in v)
        report["biomarkers_in_communities"] = n_present
        report["biomarkers_total"] = len(BIOMARKERS)
        log(f"\n  Biomarkers in communities: {n_present}/{len(BIOMARKERS)}")

    # ---- Check A3 gene presence in communities ----
    log("\n--- A3 Gene Presence in Communities ---")

    a3_presence = {}
    if partition_df is not None:
        for ensg in A3_GENES:
            alias = A3_ID_TO_ALIAS.get(ensg, ensg)
            present = ensg in comm_genes_set
            if present:
                comm = int(partition_df[partition_df["gene"] == ensg]["community"].iloc[0])
                a3_presence[alias] = f"Community {comm}"
            else:
                a3_presence[alias] = "NOT IN ANY COMMUNITY"
            log(f"  {alias} ({ensg}): {a3_presence[alias]}")


    # ================================================================
    # SECTION 3 — KEGG Pathway Enrichment Per Community
    # ================================================================
    banner(f"[SECTION 3] KEGG Enrichment — {cancer_type}")

    try:
        import gseapy as gp
        HAS_GSEAPY = True
        log("gseapy loaded successfully")
    except ImportError:
        HAS_GSEAPY = False
        log("[WARNING] gseapy not installed — skipping KEGG enrichment")
        log("  Install with: pip install gseapy")

    enrichment_results = []

    if HAS_GSEAPY and partition_df is not None:

        # Group genes by community
        comm_to_genes = {}
        for _, row in partition_df.iterrows():
            c = int(row["community"])
            g = row["gene"]
            comm_to_genes.setdefault(c, []).append(g)

        # Convert ENSG IDs to symbols for enrichment
        for c in sorted(comm_to_genes.keys()):
            genes_ensg = comm_to_genes[c]
            genes_symbols = [sym(g) for g in genes_ensg]
            # Remove any that are still ENSG (unmapped)
            genes_symbols_clean = [g for g in genes_symbols
                                   if not g.startswith("ENSG")]

            n_total = len(genes_ensg)
            n_mapped = len(genes_symbols_clean)

            log(f"\n  Community {c}: {n_total} genes ({n_mapped} mapped to symbols)")

            if n_mapped < 3:
                log(f"    [SKIP] Too few mapped genes for enrichment")
                enrichment_results.append({
                    "community": c,
                    "n_genes": n_total,
                    "n_mapped": n_mapped,
                    "top_kegg_term": "N/A (too few genes)",
                    "top_kegg_pvalue": None,
                    "top_kegg_genes": None,
                    "n_kegg_significant": 0
                })
                continue

            try:
                enr = gp.enrichr(
                    gene_list=genes_symbols_clean,
                    gene_sets=["KEGG_2021_Human"],
                    organism="human",
                    outdir=None,
                    no_plot=True,
                    verbose=False
                )

                results_df = enr.results

                if len(results_df) > 0:
                    # Filter significant results
                    sig_results = results_df[results_df["Adjusted P-value"] < 0.05]
                    n_sig = len(sig_results)

                    # Top result
                    top = results_df.iloc[0]
                    top_term = top["Term"]
                    top_pval = top["Adjusted P-value"]
                    top_genes = top["Genes"]

                    log(f"    Top KEGG: {top_term} (adj.p = {top_pval:.2e})")
                    log(f"    Significant KEGG terms (adj.p < 0.05): {n_sig}")

                    if n_sig > 0:
                        log(f"    Top 5 significant terms:")
                        for _, r in sig_results.head(5).iterrows():
                            log(f"      {r['Term']} (adj.p = {r['Adjusted P-value']:.2e}, "
                                f"genes: {r['Genes']})")

                    enrichment_results.append({
                        "community": c,
                        "n_genes": n_total,
                        "n_mapped": n_mapped,
                        "top_kegg_term": top_term,
                        "top_kegg_pvalue": top_pval,
                        "top_kegg_genes": top_genes,
                        "n_kegg_significant": n_sig
                    })

                    # Save full enrichment results for this community
                    comm_enr_path = os.path.join(REPORT_DIR,
                                                  f"{cancer_type}_community_{c:02d}_KEGG.csv")
                    results_df.to_csv(comm_enr_path, index=False)

                else:
                    log(f"    No KEGG results returned")
                    enrichment_results.append({
                        "community": c,
                        "n_genes": n_total,
                        "n_mapped": n_mapped,
                        "top_kegg_term": "No results",
                        "top_kegg_pvalue": None,
                        "top_kegg_genes": None,
                        "n_kegg_significant": 0
                    })

            except Exception as e:
                log(f"    [ERROR] Enrichment failed: {e}")
                enrichment_results.append({
                    "community": c,
                    "n_genes": n_total,
                    "n_mapped": n_mapped,
                    "top_kegg_term": f"ERROR: {str(e)[:50]}",
                    "top_kegg_pvalue": None,
                    "top_kegg_genes": None,
                    "n_kegg_significant": 0
                })

        # Save enrichment summary
        enr_summary_df = pd.DataFrame(enrichment_results)
        enr_summary_path = os.path.join(REPORT_DIR,
                                         f"{cancer_type}_KEGG_enrichment_summary.csv")
        enr_summary_df.to_csv(enr_summary_path, index=False)
        log(f"\n[SAVE] KEGG enrichment summary -> {enr_summary_path}")


    # ================================================================
    # SECTION 4 — Compile Full Report
    # ================================================================
    banner(f"[SECTION 4] Full Report — {cancer_type}")

    report_path = os.path.join(REPORT_DIR, f"{cancer_type}_pipeline_report.txt")

    with open(report_path, "w") as f:

        f.write("=" * 80 + "\n")
        f.write(f"PIPELINE SUMMARY REPORT — {cancer_type}\n")
        f.write(f"Generated: {TIMESTAMP}\n")
        f.write("=" * 80 + "\n\n")

        # --- Configuration ---
        f.write("PIPELINE CONFIGURATION\n")
        f.write("-" * 40 + "\n")
        for line in config_lines:
            f.write(line + "\n")
        f.write("\n")

        # --- Group Definitions ---
        f.write("GROUP DEFINITIONS\n")
        f.write("-" * 40 + "\n")
        for key in ["n_SBS2_HIGH", "n_SBS2_LOW",
                     "A3A_median_HIGH", "A3A_median_LOW",
                     "A3B_median_HIGH", "A3B_median_LOW",
                     "SBS2_median_HIGH", "SBS2_median_LOW",
                     "SBS2_range_HIGH", "SBS2_range_LOW"]:
            if key in report:
                f.write(f"  {key}: {report[key]}\n")
        f.write("\n")

        # --- Gene Selection ---
        f.write("GENE SELECTION\n")
        f.write("-" * 40 + "\n")
        for key in ["n_genes_tested", "n_genes_significant_raw",
                     "n_genes_selected", "n_genes_selected_filtered"]:
            if key in report:
                f.write(f"  {key}: {report[key]}\n")
        f.write("\n")

        # --- Network Statistics ---
        f.write("NETWORK STATISTICS\n")
        f.write("-" * 40 + "\n")
        for graph_label in ["G_top_noiso", "G_bot_noiso", "G_diff_noiso"]:
            if f"{graph_label}_nodes" in report:
                f.write(f"  {graph_label}:\n")
                f.write(f"    Nodes: {report[f'{graph_label}_nodes']}\n")
                f.write(f"    Edges: {report[f'{graph_label}_edges']}\n")
                f.write(f"    Avg degree: {report[f'{graph_label}_avg_degree']}\n")
        f.write("\n")

        # --- Community Detection ---
        f.write("COMMUNITY DETECTION\n")
        f.write("-" * 40 + "\n")
        for key in ["n_communities", "n_genes_in_communities",
                     "best_resolution", "best_modularity", "best_ARI"]:
            if key in report:
                f.write(f"  {key}: {report[key]}\n")
        f.write("\n  Community sizes:\n")
        for key, val in sorted(report.items()):
            if key.startswith("community_") and key.endswith("_size"):
                f.write(f"    {key}: {val}\n")
        f.write("\n")

        # --- Biomarker Presence ---
        f.write("BIOMARKER PRESENCE IN COMMUNITIES\n")
        f.write("-" * 40 + "\n")
        for symbol, status in biomarker_presence.items():
            f.write(f"  {symbol}: {status}\n")
        f.write(f"\n  Summary: {report.get('biomarkers_in_communities', 'N/A')}"
                f"/{report.get('biomarkers_total', 'N/A')} biomarkers found in communities\n\n")

        # --- A3 Gene Presence ---
        f.write("A3 GENE PRESENCE IN COMMUNITIES\n")
        f.write("-" * 40 + "\n")
        for alias, status in a3_presence.items():
            f.write(f"  {alias}: {status}\n")
        f.write("\n")

        # --- KEGG Enrichment ---
        f.write("KEGG PATHWAY ENRICHMENT PER COMMUNITY\n")
        f.write("-" * 40 + "\n")
        if enrichment_results:
            for er in enrichment_results:
                f.write(f"\n  Community {er['community']} "
                        f"({er['n_genes']} genes, {er['n_mapped']} mapped):\n")
                f.write(f"    Top KEGG term: {er['top_kegg_term']}\n")
                if er['top_kegg_pvalue'] is not None:
                    f.write(f"    Adjusted p-value: {er['top_kegg_pvalue']:.2e}\n")
                f.write(f"    Significant terms (adj.p < 0.05): {er['n_kegg_significant']}\n")
                if er['top_kegg_genes']:
                    f.write(f"    Top term genes: {er['top_kegg_genes']}\n")
        else:
            f.write("  No enrichment results (gseapy not available or no communities)\n")

        f.write("\n" + "=" * 80 + "\n")
        f.write("END OF REPORT\n")
        f.write("=" * 80 + "\n")

    log(f"\n[SAVE] Full report -> {report_path}")

    # Also save the report dict as JSON for programmatic access
    json_path = os.path.join(REPORT_DIR, f"{cancer_type}_pipeline_report.json")
    # Convert any non-serializable types
    report_serializable = {}
    for k, v in report.items():
        if isinstance(v, (np.integer,)):
            report_serializable[k] = int(v)
        elif isinstance(v, (np.floating,)):
            report_serializable[k] = float(v)
        else:
            report_serializable[k] = v

    with open(json_path, "w") as f:
        json.dump(report_serializable, f, indent=2)
    log(f"[SAVE] Report JSON -> {json_path}")


# =============================================================================
# FINAL SUMMARY
# =============================================================================
banner("STEP 08 COMPLETE")

log(f"\nOutput directory: {REPORT_DIR}")
log(f"\nFiles generated per cancer type:")
log(f"  {{cancer_type}}_pipeline_report.txt    — Human-readable full report")
log(f"  {{cancer_type}}_pipeline_report.json   — Machine-readable report")
log(f"  {{cancer_type}}_KEGG_enrichment_summary.csv — One row per community")
log(f"  {{cancer_type}}_community_XX_KEGG.csv  — Full KEGG results per community")

log(f"\nKey questions this report addresses:")
log(f"  1. What parameters defined the HIGH and LOW SBS2 groups?")
log(f"  2. How many genes survived differential expression filtering?")
log(f"  3. What were the network sizes at each threshold?")
log(f"  4. How many communities were detected and at what resolution?")
log(f"  5. Are cell-type marker genes present in any community?")
log(f"  6. What KEGG pathways are enriched in each community?")
log(f"  7. Are these known pathways or potentially novel A3-associated networks?")

banner("ALL DONE")
