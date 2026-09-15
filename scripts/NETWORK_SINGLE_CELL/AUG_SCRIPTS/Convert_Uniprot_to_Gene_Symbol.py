#!/usr/bin/env python3
"""
Convert_Uniprot_to_Gene_Symbol.py
==================================

Converts UniProt accessions to gene symbols using the UniProt REST API.
Reads ALL unique UniProt accessions from the Jang et al. (2024) mmc15.xlsx
supplementary table, checks against any existing mapping file, queries
UniProt for missing entries, and outputs a complete deduplicated mapping.

The output mapping file is used by Extract_Harris_A3_Interactors.py to
convert high-confidence AP-MS interactors to gene symbols without any
manual/hardcoded name guessing.

Input:
  - mmc15.xlsx (Jang et al. 2024 supplementary data)
  - uniprot_to_gene_symbol_mapping.tsv (existing mapping, if present)

Output:
  - uniprot_to_gene_symbol_mapping.tsv (complete, deduplicated)

Usage:
  python Convert_Uniprot_to_Gene_Symbol.py

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
import sys
import time
import requests
import pandas as pd
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"

# Input: mmc15.xlsx for extracting all unique accessions
MMC15_PATH = os.path.join(BASE_DIR, "data", "FIG_4", "00_input", "mmc15.xlsx")

# Output: mapping file (same location as before)
MAPPING_FILE = os.path.join(BASE_DIR, "scripts", "uniprot_to_gene_symbol_mapping.tsv")

# UniProt REST API
UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb/search"
BATCH_SIZE = 100  # Query up to 100 accessions per API call
SLEEP_BETWEEN_BATCHES = 1.0  # seconds


# =============================================================================
# HELPERS
# =============================================================================

def log(msg):
    timestamp = datetime.now().strftime('%H:%M:%S')
    print(f"[{timestamp}] {msg}", flush=True)


def query_uniprot_batch(accessions):
    """
    Query UniProt REST API for a batch of accessions.
    Returns list of dicts with UniProt_Accession, UniProt_Entry_Name, Gene_Symbol.
    """
    results = []

    # Build OR query: (accession:P12345) OR (accession:Q67890) ...
    query_parts = [f"(accession:{acc})" for acc in accessions]
    query_str = " OR ".join(query_parts)

    params = {
        'query': query_str,
        'format': 'tsv',
        'fields': 'accession,id,gene_primary',
        'size': 500,
    }

    try:
        response = requests.get(UNIPROT_API_URL, params=params, timeout=30)
        response.raise_for_status()

        lines = response.text.strip().split('\n')
        if len(lines) <= 1:
            return results

        # Parse TSV response (header + data lines)
        header = lines[0].split('\t')
        for line in lines[1:]:
            fields = line.split('\t')
            if len(fields) >= 3:
                acc = fields[0].strip()
                entry_name = fields[1].strip()
                gene_symbol = fields[2].strip()

                if gene_symbol:  # Only include if gene symbol is present
                    results.append({
                        'UniProt_Accession': acc,
                        'UniProt_Entry_Name': entry_name,
                        'Gene_Symbol': gene_symbol,
                    })

    except requests.exceptions.RequestException as e:
        log(f"    API error: {e}")

    return results


# =============================================================================
# MAIN
# =============================================================================

def main():
    log("=" * 80)
    log("UNIPROT ACCESSION → GENE SYMBOL CONVERSION")
    log("=" * 80)

    # -----------------------------------------------------------------
    # STEP 1: Extract ALL unique accessions from mmc15.xlsx
    # -----------------------------------------------------------------
    log(f"\nSTEP 1: Extracting unique UniProt accessions from mmc15.xlsx")

    if not os.path.exists(MMC15_PATH):
        log(f"  ERROR: {MMC15_PATH} not found")
        sys.exit(1)

    df = pd.read_excel(MMC15_PATH, sheet_name='Saint and CompPASS Scores')
    df = df.dropna(subset=['Prey.x'])

    # Get all unique accessions and their entry names
    prey_info = df[['Prey.x', 'PreyGene']].drop_duplicates()
    all_accessions = set(prey_info['Prey.x'].astype(str).str.strip().unique())

    # Also build accession → entry name mapping from the file
    acc_to_entry_from_file = {}
    for _, row in prey_info.iterrows():
        acc = str(row['Prey.x']).strip()
        entry = str(row['PreyGene']).strip()
        acc_to_entry_from_file[acc] = entry

    log(f"  Unique accessions in mmc15.xlsx: {len(all_accessions)}")

    # -----------------------------------------------------------------
    # STEP 2: Load existing mapping file (if present)
    # -----------------------------------------------------------------
    log(f"\nSTEP 2: Loading existing mapping file")

    existing_mappings = {}
    if os.path.exists(MAPPING_FILE):
        existing_df = pd.read_csv(MAPPING_FILE, sep='\t')
        log(f"  Existing mappings: {len(existing_df)}")

        for _, row in existing_df.iterrows():
            acc = str(row['UniProt_Accession']).strip()
            existing_mappings[acc] = {
                'UniProt_Accession': acc,
                'UniProt_Entry_Name': str(row['UniProt_Entry_Name']).strip(),
                'Gene_Symbol': str(row['Gene_Symbol']).strip(),
            }
    else:
        log(f"  No existing mapping file found — will create from scratch")

    # -----------------------------------------------------------------
    # STEP 3: Identify missing accessions
    # -----------------------------------------------------------------
    log(f"\nSTEP 3: Identifying accessions needing conversion")

    already_mapped = set(existing_mappings.keys())
    missing = all_accessions - already_mapped
    log(f"  Already mapped: {len(already_mapped & all_accessions)}")
    log(f"  Need conversion: {len(missing)}")

    if len(missing) == 0:
        log("  All accessions already mapped — nothing to do")
    else:
        log(f"\n  Missing accessions:")
        for acc in sorted(missing):
            entry = acc_to_entry_from_file.get(acc, '?')
            log(f"    {acc} ({entry})")

    # -----------------------------------------------------------------
    # STEP 4: Query UniProt API for missing accessions
    # -----------------------------------------------------------------
    if len(missing) > 0:
        log(f"\nSTEP 4: Querying UniProt REST API for {len(missing)} accessions")

        missing_list = sorted(missing)
        new_mappings = []

        for i in range(0, len(missing_list), BATCH_SIZE):
            batch = missing_list[i:i + BATCH_SIZE]
            log(f"  Batch {i // BATCH_SIZE + 1}: querying {len(batch)} accessions...")

            batch_results = query_uniprot_batch(batch)
            new_mappings.extend(batch_results)

            log(f"    Retrieved: {len(batch_results)} mappings")

            if i + BATCH_SIZE < len(missing_list):
                time.sleep(SLEEP_BETWEEN_BATCHES)

        log(f"\n  Total new mappings from API: {len(new_mappings)}")

        # Add new mappings to existing
        for mapping in new_mappings:
            acc = mapping['UniProt_Accession']
            existing_mappings[acc] = mapping

        # Check for still-unmapped
        mapped_accs = set(existing_mappings.keys())
        still_missing = missing - mapped_accs
        if still_missing:
            log(f"\n  WARNING: {len(still_missing)} accessions could not be resolved by UniProt API:")
            for acc in sorted(still_missing):
                entry = acc_to_entry_from_file.get(acc, '?')
                log(f"    {acc} ({entry})")
    else:
        log(f"\nSTEP 4: Skipped — no missing accessions")

    # -----------------------------------------------------------------
    # STEP 5: Build complete mapping and deduplicate
    # -----------------------------------------------------------------
    log(f"\nSTEP 5: Building final mapping file")

    # Convert to DataFrame
    all_rows = list(existing_mappings.values())
    final_df = pd.DataFrame(all_rows)

    # Deduplicate by accession (keep first)
    n_before = len(final_df)
    final_df = final_df.drop_duplicates(subset='UniProt_Accession', keep='first')
    n_after = len(final_df)
    if n_before != n_after:
        log(f"  Removed {n_before - n_after} duplicate accessions")

    # Sort by accession
    final_df = final_df.sort_values('UniProt_Accession').reset_index(drop=True)

    log(f"  Final mapping entries: {len(final_df)}")

    # Verify coverage of mmc15.xlsx accessions
    final_accs = set(final_df['UniProt_Accession'])
    coverage = all_accessions & final_accs
    uncovered = all_accessions - final_accs
    log(f"  Coverage of mmc15.xlsx: {len(coverage)} / {len(all_accessions)} "
        f"({100 * len(coverage) / len(all_accessions):.1f}%)")
    if uncovered:
        log(f"  Still unmapped: {len(uncovered)}")

    # -----------------------------------------------------------------
    # STEP 6: Save
    # -----------------------------------------------------------------
    log(f"\nSTEP 6: Saving mapping file")

    final_df.to_csv(MAPPING_FILE, sep='\t', index=False)
    log(f"  Saved: {MAPPING_FILE}")
    log(f"  Entries: {len(final_df)}")

    # Summary
    log(f"\n{'='*80}")
    log("CONVERSION COMPLETE")
    log(f"{'='*80}")
    log(f"  Mapping file: {MAPPING_FILE}")
    log(f"  Total entries: {len(final_df)}")
    log(f"  mmc15.xlsx coverage: {len(coverage)} / {len(all_accessions)}")
    if uncovered:
        log(f"  Unmapped accessions: {len(uncovered)}")
        log(f"  (These may be obsolete/deprecated UniProt entries)")


if __name__ == "__main__":
    main()
