#!/usr/bin/env python3
"""
Diagnostic_HLA_panel.py

Confirm the HLA class I allele panel actually used in the MHCflurry runs by
pulling the unique alleles out of the per-group peptide-result files. Verifies
the count (expected 10) and lists them so the Methods/Results text can name the
real set rather than asserting "10 common alleles" without specifics.

Read-only. Run in NETWORK env on titan.
"""

import os
import pandas as pd

MHC_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_7/03_mhc_binding"
ALLELE_CANDIDATES = ["allele", "hla", "hla_allele", "mhc_allele",
                     "best_allele", "Allele", "HLA"]

panel = {}
for group in ["SBS2_HIGH", "CNV_HIGH"]:
    for suffix in ["_all_peptide_results.tsv", "_neoantigens.tsv"]:
        path = os.path.join(MHC_DIR, f"{group}{suffix}")
        if not os.path.exists(path):
            continue
        df = pd.read_csv(path, sep="\t")
        col = next((c for c in ALLELE_CANDIDATES if c in df.columns), None)
        if col is None:
            print(f"  [no allele column in {os.path.basename(path)}; "
                  f"columns: {list(df.columns)}]")
            continue
        alleles = sorted(df[col].dropna().astype(str).unique())
        panel.setdefault(group, set()).update(alleles)
        print(f"{group}{suffix}: column '{col}', {len(alleles)} unique alleles")

print("\n" + "=" * 60)
all_alleles = sorted(set().union(*panel.values())) if panel else []
print(f"Union of alleles across files: {len(all_alleles)}")
for a in all_alleles:
    print(f"   {a}")
print("=" * 60)
print(f"COUNT CHECK: {len(all_alleles)} alleles "
      f"({'matches text (10)' if len(all_alleles) == 10 else 'DOES NOT match text claim of 10'})")
