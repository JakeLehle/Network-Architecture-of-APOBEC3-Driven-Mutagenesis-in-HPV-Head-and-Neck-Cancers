#!/usr/bin/env python3
"""
Generate_High_SBS2_VCF_List.py
===============================
Generates a collaborator-ready TSV listing all VCF files for HIGH SBS2 group
patients across the 7 testable TCGA cancer types.

Replicates the identical group selection logic from
Pan_Cancer_Germline_Enrichment.py, then maps HIGH group patients to their
MuTect2 annotated VCF files via the GDC download manifest.

Output columns:
  cancer_type, case_id, sbs2_count, a3_expression, vcf_filename, full_path

Usage:
  python Generate_High_SBS2_VCF_List.py
"""

import os
import numpy as np
import pandas as pd

# =============================================================================
# CONFIGURATION — matches Pan_Cancer_Germline_Enrichment.py exactly
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
SHARED_VCF = "/master/jlehle/SHARED/TCGA/VCF"

EXPRESSION_PATH   = os.path.join(PROJECT_ROOT, "data", "FIG_1", "TCGA_master_FPKM_UQ.tsv")
METADATA_PATH     = os.path.join(PROJECT_ROOT, "data", "FIG_1", "TCGA_sample_metadata_final.tsv")
CROSSWALK_PATH    = os.path.join(PROJECT_ROOT, "data", "FIG_1", "Mutation_Table_Tumors_TCGA.tsv")
NEW_COUNTS_PATH   = os.path.join(SHARED_VCF, "SigProfiler_output", "TCGA_SBS_signature_counts.tsv")
FEASIBILITY_PATH  = os.path.join(PROJECT_ROOT, "data", "FIG_GERMLINE", "pan_cancer_feasibility.tsv")
MANIFEST_PATH     = os.path.join(SHARED_VCF, "manifests", "TCGA_MuTect2_VCF_only_manifest.tsv")
VCF_BASE_DIR      = os.path.join(SHARED_VCF, "MuTect2_Annotated")

OUTPUT_PATH       = os.path.join(PROJECT_ROOT, "data", "FIG_GERMLINE",
                                 "high_sbs2_vcf_file_list.tsv")

# Group selection parameters — MUST MATCH network_config.py / Step03
A3_SUM_PERCENTILE  = 0.50
SBS2_GROUP_FRACTION = 0.25
MIN_GROUP_SIZE     = 8

# A3 gene list for expression sum
A3_GENES = ["APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
            "APOBEC3F", "APOBEC3G", "APOBEC3H"]

# =============================================================================
# MAIN
# =============================================================================
def main():
    print("=" * 72)
    print("  Generate HIGH SBS2 VCF File List for Collaborator")
    print("=" * 72)

    # -------------------------------------------------------------------------
    # 1. Load expression data
    # -------------------------------------------------------------------------
    print("\nLoading expression data...")
    expr = pd.read_csv(EXPRESSION_PATH, sep="\t", index_col=0)
    print(f"  Expression matrix: {expr.shape[0]} genes x {expr.shape[1]} samples")

    # Compute A3 sum (A3A + A3B used for filtering, matching enrichment script)
    a3_present = [g for g in ["APOBEC3A", "APOBEC3B"] if g in expr.index]
    if len(a3_present) < 2:
        print(f"  WARNING: Only found {a3_present} in expression matrix")
    a3_sum = expr.loc[a3_present].sum(axis=0)

    # Also compute full A3 sum for reporting
    a3_all_present = [g for g in A3_GENES if g in expr.index]
    a3_full_sum = expr.loc[a3_all_present].sum(axis=0)

    # -------------------------------------------------------------------------
    # 2. Load SBS2 counts
    # -------------------------------------------------------------------------
    print("Loading SBS2 signature counts...")
    sbs = pd.read_csv(NEW_COUNTS_PATH, sep="\t", index_col=0)
    if "SBS2" not in sbs.columns:
        raise ValueError("SBS2 column not found in signature counts file")
    sbs2 = sbs["SBS2"]
    print(f"  SBS2 data: {len(sbs2)} samples")

    # -------------------------------------------------------------------------
    # 3. Load crosswalk (maps expression barcodes to mutation barcodes)
    # -------------------------------------------------------------------------
    print("Loading barcode crosswalk...")
    crosswalk = pd.read_csv(CROSSWALK_PATH, sep="\t")
    print(f"  Crosswalk: {len(crosswalk)} entries")

    # -------------------------------------------------------------------------
    # 4. Load metadata for cancer type assignment
    # -------------------------------------------------------------------------
    print("Loading sample metadata...")
    meta = pd.read_csv(METADATA_PATH, sep="\t")
    print(f"  Metadata: {len(meta)} entries")

    # -------------------------------------------------------------------------
    # 5. Load testable cancer types from Phase 1
    # -------------------------------------------------------------------------
    print("Loading testable cancer types...")
    with open(os.path.join(PROJECT_ROOT, "data", "FIG_GERMLINE",
                           "testable_cancer_types.txt")) as f:
        testable = [line.strip() for line in f if line.strip()]
    print(f"  Testable cancers: {len(testable)} -> {testable}")

    # -------------------------------------------------------------------------
    # 6. Load VCF manifest
    # -------------------------------------------------------------------------
    print("Loading VCF manifest...")
    manifest = pd.read_csv(MANIFEST_PATH, sep="\t")
    manifest_vcf = manifest[manifest["Is_VCF"] == True].copy()
    print(f"  Manifest total: {len(manifest)} files")
    print(f"  VCF files only: {len(manifest_vcf)} files")

    # -------------------------------------------------------------------------
    # 7. For each testable cancer, identify HIGH group and map to VCFs
    # -------------------------------------------------------------------------
    print("\n" + "=" * 72)
    print("  Identifying HIGH SBS2 groups per cancer")
    print("=" * 72)

    all_results = []

    for cancer in sorted(testable):
        print(f"\n  --- {cancer} ---")

        # Get samples for this cancer type from metadata
        cancer_meta = meta[meta["cancer_type"] == cancer]
        cancer_barcodes = set(cancer_meta.iloc[:, 0])  # first column = barcode

        # Find samples present in both expression and crosswalk
        # Build barcode mapping: expression barcode -> mutation barcode
        xwalk_cancer = crosswalk[crosswalk.iloc[:, 0].isin(cancer_barcodes)]

        # Get expression barcodes that have SBS2 data via crosswalk
        expr_barcodes = []
        expr_to_mut = {}
        for _, row in xwalk_cancer.iterrows():
            expr_bc = row.iloc[0]
            mut_bc = row.iloc[1] if len(row) > 1 else expr_bc
            if expr_bc in a3_sum.index and mut_bc in sbs2.index:
                expr_barcodes.append(expr_bc)
                expr_to_mut[expr_bc] = mut_bc

        if len(expr_barcodes) == 0:
            print(f"    No overlapping samples, skipping")
            continue

        # Step 1: Filter to tumors above median A3A+A3B expression
        a3_vals = a3_sum[expr_barcodes]
        median_a3 = a3_vals.median()
        high_a3_barcodes = a3_vals[a3_vals >= median_a3].index.tolist()
        n_high_a3 = len(high_a3_barcodes)

        # Step 2: Rank by SBS2
        sbs2_vals = pd.Series(
            {bc: sbs2[expr_to_mut[bc]] for bc in high_a3_barcodes}
        ).sort_values()

        # Step 3: Group size
        n_per_group = int(np.floor(n_high_a3 * SBS2_GROUP_FRACTION))

        if n_per_group < MIN_GROUP_SIZE:
            print(f"    n_per_group={n_per_group} < {MIN_GROUP_SIZE}, skipping")
            continue

        # Step 4: LOW = bottom n, HIGH = top n
        low_group = sbs2_vals.head(n_per_group).index.tolist()
        high_group = sbs2_vals.tail(n_per_group).index.tolist()

        print(f"    Samples with expression+SBS2: {len(expr_barcodes)}")
        print(f"    Above median A3: {n_high_a3}")
        print(f"    n_per_group: {n_per_group}")
        print(f"    HIGH group SBS2 range: {sbs2_vals[high_group].min():.0f} - {sbs2_vals[high_group].max():.0f}")
        print(f"    HIGH group SBS2 median: {sbs2_vals[high_group].median():.0f}")

        # Extract Case_IDs for HIGH group (first 12 characters of barcode)
        for bc in high_group:
            case_id = bc[:12]  # TCGA-XX-XXXX
            all_results.append({
                "cancer_type": cancer,
                "case_id": case_id,
                "expression_barcode": bc,
                "sbs2_count": int(sbs2[expr_to_mut[bc]]),
                "a3ab_expression": round(a3_sum[bc], 2),
            })

    # -------------------------------------------------------------------------
    # 8. Map to VCF files via manifest
    # -------------------------------------------------------------------------
    print("\n" + "=" * 72)
    print("  Mapping HIGH group patients to VCF files")
    print("=" * 72)

    results_df = pd.DataFrame(all_results)
    print(f"  Total HIGH group patients across all cancers: {len(results_df)}")

    # Build Case_ID lookup from manifest
    manifest_vcf["Case_ID_clean"] = manifest_vcf["Case_ID"].str.strip()

    # Merge
    merged = results_df.merge(
        manifest_vcf[["Cancer_Type", "Case_ID", "File_Name"]],
        left_on="case_id",
        right_on="Case_ID",
        how="left"
    )

    # Build full path
    merged["full_path"] = merged.apply(
        lambda row: os.path.join(VCF_BASE_DIR, f"TCGA-{row['cancer_type']}",
                                  row["File_Name"])
        if pd.notna(row["File_Name"]) else "NOT_FOUND",
        axis=1
    )

    # Report mapping success
    n_mapped = merged["File_Name"].notna().sum()
    n_unmapped = merged["File_Name"].isna().sum()
    print(f"  Successfully mapped to VCF: {n_mapped}")
    if n_unmapped > 0:
        print(f"  WARNING: {n_unmapped} patients could not be mapped to a VCF file")
        unmapped = merged[merged["File_Name"].isna()]["case_id"].unique()
        print(f"    Unmapped Case_IDs: {list(unmapped)[:10]}...")

    # Clean up output columns
    output = merged[[
        "cancer_type", "case_id", "expression_barcode",
        "sbs2_count", "a3ab_expression", "File_Name", "full_path"
    ]].rename(columns={"File_Name": "vcf_filename"})

    output = output.sort_values(["cancer_type", "sbs2_count"],
                                 ascending=[True, False])

    # -------------------------------------------------------------------------
    # 9. Write output
    # -------------------------------------------------------------------------
    output.to_csv(OUTPUT_PATH, sep="\t", index=False)
    print(f"\n  Output written to: {OUTPUT_PATH}")
    print(f"  Total rows: {len(output)}")

    # Summary per cancer
    print("\n  Per-cancer summary:")
    print(f"  {'Cancer':<8} {'Patients':>10} {'VCF files':>10} {'SBS2 median':>12}")
    print(f"  {'-'*8} {'-'*10} {'-'*10} {'-'*12}")
    for cancer in sorted(testable):
        cancer_rows = output[output["cancer_type"] == cancer]
        n_patients = cancer_rows["case_id"].nunique()
        n_vcfs = cancer_rows["vcf_filename"].notna().sum()
        med_sbs2 = cancer_rows["sbs2_count"].median()
        print(f"  {cancer:<8} {n_patients:>10} {n_vcfs:>10} {med_sbs2:>12.0f}")

    print(f"\n  Total unique patients: {output['case_id'].nunique()}")
    print(f"  Total VCF files: {output['vcf_filename'].notna().sum()}")
    print("\n  Done.")


if __name__ == "__main__":
    main()
