#!/usr/bin/env python3
"""
HNSC_Somatic_Enrichment_Analysis.py
====================================
Somatic mutation enrichment between SBS2-HIGH and SBS2-LOW HNSC tumors.

Tests whether tumors with high APOBEC mutational burden acquire different
somatic mutations (gene-level) compared to tumors with low burden, when
both groups have comparable A3 expression.

Phases:
  0.  Load HIGH/LOW groups from network pipeline Step03 output
  0B. Load reference gene sets (Harris, communities, DDR, chromatin)
  1.  Load per-cancer MAF, filter to PASS somatic, compute TMB
  2.  Track A: gene-level burden (non-silent, damaging, silent control)
      + TMB-adjusted logistic regression
  2C. Frequency tier binning + per-tier gene set cross-reference + KEGG
  3.  Track B: A3-specific somatic catalog
  4.  APOBEC trinucleotide context from per-sample MAF files
  5.  Track C: gene set enrichment (hypergeometric + KEGG via Enrichr)

Groups use the SAME Step03 selection as the network pipeline to ensure
identical patient sets across all analyses.

Usage:
    conda run -n NETWORK python HNSC_Somatic_Enrichment_Analysis.py

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os, sys, time, json, glob, pickle, warnings
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm

warnings.filterwarnings("ignore", category=FutureWarning)

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
SHARED_VCF = "/master/jlehle/SHARED/TCGA/VCF"

# Step03 group outputs (from updated network pipeline)
CANCER_TYPE = "TCGA-HNSC"
DIR_03 = os.path.join(PROJECT_ROOT, "data", "FIG_2", "03_differential_expression", CANCER_TYPE)
HIGH_GROUP_PKL = os.path.join(DIR_03, f"{CANCER_TYPE}_SBS2_HIGH_group.pkl")
LOW_GROUP_PKL = os.path.join(DIR_03, f"{CANCER_TYPE}_SBS2_LOW_group.pkl")

# Per-cancer MAF (full annotations, all variant types)
MAF_PATH = os.path.join(SHARED_VCF, "consolidated", "per_cancer", "TCGA-HNSC_mutations.maf.tsv")

# Per-sample MAF directory (for CONTEXT column / trinucleotide analysis)
PER_SAMPLE_MAF_DIR = os.path.join(SHARED_VCF, "MuTect2_Annotated", "TCGA-HNSC")

# Gene symbol mapping
ENSG_TO_SYMBOL_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_2",
                                     "01_cleaned_expression", "ensg_to_symbol.json")

# Network community gene lists
COMMUNITY_GENE_LISTS_PATH = os.path.join(
    PROJECT_ROOT, "data", "FIG_2", "05_communities", CANCER_TYPE,
    f"{CANCER_TYPE}_community_gene_lists.csv"
)

# Harris A3 interactor lists
HARRIS_ALL_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_4", "00_input",
                                "Harris_A3_interactors.txt")
HARRIS_A3B_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_4", "00_input",
                                "Harris_A3_interactors_A3B_only.txt")

# Output
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "FIG_1", "SOMATIC_ENRICHMENT")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Track A parameters
TRACK_A_MIN_CARRIERS = 3    # gene must be mutated in >= 3 patients total
BH_THRESHOLD = 0.10

# Variant classification groups
NONSILENT_CLASSES = [
    "Missense_Mutation", "Nonsense_Mutation", "Splice_Site",
    "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
    "Nonstop_Mutation", "Translation_Start_Site"
]
DAMAGING_CLASSES = [
    "Nonsense_Mutation", "Splice_Site", "Frame_Shift_Del", "Frame_Shift_Ins",
    "Nonstop_Mutation"
]
SILENT_CLASSES = ["Silent"]

# A3 gene symbols
A3_GENE_SYMBOLS = [
    'APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D',
    'APOBEC3F', 'APOBEC3G', 'APOBEC3H',
]

# MAF columns to read
MAF_COLS = [
    'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position',
    'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Type',
    'Variant_Classification', 'FILTER', 'Tumor_Sample_Barcode',
    'dbSNP_RS', 'IMPACT', 'Consequence', 'HGVSp_Short',
    'SIFT', 'PolyPhen', 'EXON', 'INTRON'
]

# Per-sample MAF columns (for CONTEXT analysis)
CONTEXT_MAF_COLS = [
    'Tumor_Sample_Barcode', 'Hugo_Symbol', 'Reference_Allele',
    'Tumor_Seq_Allele2', 'Variant_Type', 'FILTER', 'CONTEXT'
]

# SBS96 classification (from Run_SigProfiler.py)
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

# DDR gene list (curated from Knijnenburg et al. 2018 TCGA DDR paper)
DDR_GENES = [
    # Base excision repair (BER)
    "OGG1", "MUTYH", "NEIL1", "NEIL2", "NEIL3", "MPG", "NTHL1", "SMUG1",
    "TDG", "UNG", "APEX1", "APEX2", "PNKP", "XRCC1", "PARP1", "PARP2",
    "LIG3", "POLB", "FEN1",
    # Nucleotide excision repair (NER)
    "XPA", "XPC", "RAD23B", "CETN2", "DDB1", "DDB2", "ERCC1", "ERCC2",
    "ERCC3", "ERCC4", "ERCC5", "ERCC6", "ERCC8", "RPA1", "RPA2", "RPA3",
    "LIG1",
    # Mismatch repair (MMR)
    "MLH1", "MLH3", "MSH2", "MSH3", "MSH6", "PMS1", "PMS2", "EXO1",
    "PCNA",
    # Homologous recombination (HR)
    "BRCA1", "BRCA2", "RAD51", "RAD51B", "RAD51C", "RAD51D", "XRCC2",
    "XRCC3", "RAD54L", "RAD54B", "NBN", "MRE11", "RAD50", "PALB2",
    "BRIP1", "BLM", "WRN", "RECQL", "RECQL4", "RECQL5", "ATM", "ATR",
    "ATRIP", "CHEK1", "CHEK2", "MDC1", "RNF8", "RNF168", "TP53BP1",
    "RBBP8",
    # Non-homologous end joining (NHEJ)
    "XRCC4", "XRCC5", "XRCC6", "LIG4", "DCLRE1C", "PRKDC", "NHEJ1",
    # Fanconi anemia
    "FANCA", "FANCB", "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG",
    "FANCI", "FANCL", "FANCM", "UBE2T", "SLX4",
    # Translesion synthesis
    "REV1", "REV3L", "POLH", "POLI", "POLK", "POLQ",
]

# Chromatin remodelers / SWI-SNF + broader epigenetic modifiers
CHROMATIN_GENES = [
    # SWI/SNF complex
    "SMARCA2", "SMARCA4", "SMARCB1", "SMARCC1", "SMARCC2", "SMARCD1",
    "SMARCD2", "SMARCD3", "SMARCE1", "ARID1A", "ARID1B", "ARID2",
    "PBRM1", "BRD7", "BRD9", "DPF1", "DPF2", "DPF3", "PHF10",
    # Other chromatin modifiers
    "KDM6A", "KMT2A", "KMT2C", "KMT2D", "SETD2", "NSD1", "NSD2",
    "EZH2", "SUZ12", "EED", "EP300", "CREBBP", "KAT6A", "KAT6B",
    "DNMT1", "DNMT3A", "DNMT3B", "TET1", "TET2", "TET3",
    "IDH1", "IDH2", "ATRX", "DAXX", "H3-3A", "H3-3B",
    "CHD1", "CHD2", "CHD3", "CHD4", "CHD5", "CHD6", "CHD7", "CHD8",
]

# =============================================================================
# LOGGING
# =============================================================================
report_lines = []
def log(msg=""):
    print(msg, flush=True); report_lines.append(str(msg))
def banner(title, char="="):
    log(""); log(char * 80); log(f"  {title}"); log(char * 80)


# =============================================================================
# PHASE 0: Load groups from Step03
# =============================================================================
banner("PHASE 0: Load HIGH/LOW groups from network pipeline Step03")

# Verify required files exist
for path, desc in [
    (HIGH_GROUP_PKL, "HIGH group pkl"),
    (LOW_GROUP_PKL, "LOW group pkl"),
    (MAF_PATH, "Per-cancer MAF"),
    (ENSG_TO_SYMBOL_PATH, "ENSG-to-symbol mapping"),
    (COMMUNITY_GENE_LISTS_PATH, "Community gene lists"),
    (HARRIS_ALL_PATH, "Harris A3 interactors (all)"),
    (HARRIS_A3B_PATH, "Harris A3B interactors"),
]:
    if not os.path.exists(path):
        raise FileNotFoundError(f"REQUIRED FILE MISSING: {desc}\n  {path}")
    log(f"  Found: {desc}")

# Load groups
high_df = pd.read_pickle(HIGH_GROUP_PKL)
low_df = pd.read_pickle(LOW_GROUP_PKL)
log(f"\n  HIGH group: {len(high_df)} tumors")
log(f"  LOW group:  {len(low_df)} tumors")

# Extract patient IDs (first 12 chars of Entity_ID)
high_patients = set(high_df["Entity_ID"].str[:12].values)
low_patients = set(low_df["Entity_ID"].str[:12].values)
all_patients = high_patients | low_patients
n_high = len(high_patients)
n_low = len(low_patients)

log(f"  HIGH patients: {n_high}")
log(f"  LOW patients:  {n_low}")

if n_high != n_low:
    log(f"  WARNING: unequal group sizes ({n_high} vs {n_low})")

# SBS2 distribution sanity check
if "SBS2" in high_df.columns:
    log(f"\n  SBS2 distribution:")
    log(f"    HIGH: median={high_df['SBS2'].median():.0f}, "
        f"min={high_df['SBS2'].min():.0f}, max={high_df['SBS2'].max():.0f}")
    log(f"    LOW:  median={low_df['SBS2'].median():.0f}, "
        f"min={low_df['SBS2'].min():.0f}, max={low_df['SBS2'].max():.0f}")

# Check A3 expression
for a3_col in ["A3A", "A3B", "APOBEC3A", "APOBEC3B"]:
    if a3_col in high_df.columns:
        log(f"    {a3_col} HIGH mean={high_df[a3_col].mean():.1f}, "
            f"LOW mean={low_df[a3_col].mean():.1f}")


# =============================================================================
# PHASE 0B: Load Reference Gene Sets
# =============================================================================
# Loaded early so all downstream phases (2C tier cross-referencing, 5A/5B
# gene set enrichment) can use the same definitions without duplication.
# =============================================================================
banner("PHASE 0B: Load Reference Gene Sets")

# 1. Harris A3 interactors (all)
with open(HARRIS_ALL_PATH) as f:
    harris_all = set(line.strip() for line in f if line.strip())
log(f"  Harris A3 interactors (all): {len(harris_all)} genes")

# 2. Harris A3B-specific
with open(HARRIS_A3B_PATH) as f:
    harris_a3b = set(line.strip() for line in f if line.strip())
log(f"  Harris A3B interactors: {len(harris_a3b)} genes")

# 3. Network community genes
with open(ENSG_TO_SYMBOL_PATH) as f:
    ensg_to_symbol = json.load(f)

comm_df = pd.read_csv(COMMUNITY_GENE_LISTS_PATH)
community_gene_sets = {}
all_network_genes = set()

for _, row in comm_df.iterrows():
    comm_id = row['community']
    gene_list = str(row['genes']).split(';')

    # Detect format: if genes start with ENSG, convert to symbols;
    # if they're already symbols (e.g., CASP8, APOBEC3B), use directly.
    sample_genes = [g.strip() for g in gene_list[:5] if g.strip()]
    is_ensg = any(g.startswith('ENSG') for g in sample_genes)

    symbols = set()
    if is_ensg:
        for ensg in gene_list:
            sym = ensg_to_symbol.get(ensg.strip(), '')
            if sym:
                symbols.add(sym)
    else:
        # Already gene symbols
        symbols = set(g.strip() for g in gene_list if g.strip())

    if symbols:
        community_gene_sets[f"Community_{comm_id}"] = symbols
        all_network_genes |= symbols

gene_format = "ENSG (converted)" if any(
    g.startswith('ENSG') for g in str(comm_df.iloc[0]['genes']).split(';')[:3]
) else "gene symbols (direct)"
log(f"  Network community gene format: {gene_format}")
log(f"  Network communities loaded: {len(community_gene_sets)}")
for comm_name, genes in sorted(community_gene_sets.items()):
    log(f"    {comm_name}: {len(genes)} genes")
log(f"  Total network genes (all communities): {len(all_network_genes)}")

if len(community_gene_sets) == 0:
    raise ValueError(
        "Community gene loading produced 0 communities!\n"
        f"  File: {COMMUNITY_GENE_LISTS_PATH}\n"
        f"  Columns: {comm_df.columns.tolist()}\n"
        f"  First gene entry: {comm_df.iloc[0]['genes'][:100] if len(comm_df) > 0 else 'N/A'}"
    )

# 4. DDR genes
ddr_set = set(DDR_GENES)
log(f"  DDR genes: {len(ddr_set)}")

# 5. Chromatin remodelers
chromatin_set = set(CHROMATIN_GENES)
log(f"  Chromatin remodelers: {len(chromatin_set)}")

# Consolidated dictionary for cross-referencing (used by Phase 2C and Phase 5A)
REFERENCE_GENE_SETS = {
    'Harris_A3_all': harris_all,
    'Harris_A3B_specific': harris_a3b,
    'DDR_genes': ddr_set,
    'Chromatin_remodelers': chromatin_set,
    'Network_all_communities': all_network_genes,
}
# Add individual communities
for comm_name, comm_genes in community_gene_sets.items():
    REFERENCE_GENE_SETS[comm_name] = comm_genes

log(f"\n  Reference gene sets consolidated: {len(REFERENCE_GENE_SETS)} sets")


# =============================================================================
# PHASE 1: Load MAF, filter to PASS somatic, compute TMB
# =============================================================================
banner("PHASE 1: Load MAF and compute TMB")

log(f"  Reading: {MAF_PATH}")
maf = pd.read_csv(MAF_PATH, sep='\t', usecols=MAF_COLS, low_memory=False)
log(f"  Total rows: {len(maf):,}")

# Patient ID from WES barcode
maf['patient_id'] = maf['Tumor_Sample_Barcode'].str[:12]

# Filter to our patients
maf_ours = maf[maf['patient_id'].isin(all_patients)].copy()
log(f"  Rows for our patients: {len(maf_ours):,}")
log(f"  Patients found in MAF: {maf_ours['patient_id'].nunique()}")

# Check for patients missing from MAF
missing_high = high_patients - set(maf_ours['patient_id'].unique())
missing_low = low_patients - set(maf_ours['patient_id'].unique())
if missing_high:
    log(f"  WARNING: {len(missing_high)} HIGH patients missing from MAF: {missing_high}")
if missing_low:
    log(f"  WARNING: {len(missing_low)} LOW patients missing from MAF: {missing_low}")

# Filter to PASS somatic variants
somatic_mask = (maf_ours['FILTER'] == 'PASS')
somatic = maf_ours[somatic_mask].copy()
log(f"\n  PASS somatic variants: {len(somatic):,}")
log(f"  Non-PASS (filtered out): {len(maf_ours) - len(somatic):,}")

# Variant type breakdown
log(f"\n  Variant types in PASS somatic:")
for vt, count in somatic['Variant_Type'].value_counts().items():
    log(f"    {vt}: {count:,}")

# Variant classification breakdown
log(f"\n  Top 15 Variant_Classification in PASS somatic:")
for vc, count in somatic['Variant_Classification'].value_counts().head(15).items():
    log(f"    {vc}: {count:,}")

# Assign group labels
somatic['group'] = somatic['patient_id'].map(
    lambda p: 'HIGH' if p in high_patients else ('LOW' if p in low_patients else 'OTHER')
)

# --- TMB per patient ---
banner("PHASE 1B: Tumor Mutation Burden (TMB)")

tmb_all = somatic.groupby('patient_id').size().reset_index(name='TMB_all')
tmb_snp = somatic[somatic['Variant_Type'] == 'SNP'].groupby('patient_id').size().reset_index(name='TMB_snp')
tmb_nonsilent = somatic[somatic['Variant_Classification'].isin(NONSILENT_CLASSES)].groupby(
    'patient_id').size().reset_index(name='TMB_nonsilent')

tmb = tmb_all.merge(tmb_snp, on='patient_id', how='left').merge(
    tmb_nonsilent, on='patient_id', how='left').fillna(0)
tmb['group'] = tmb['patient_id'].map(
    lambda p: 'HIGH' if p in high_patients else 'LOW'
)

tmb_high = tmb[tmb['group'] == 'HIGH']
tmb_low = tmb[tmb['group'] == 'LOW']

for metric in ['TMB_all', 'TMB_snp', 'TMB_nonsilent']:
    h_vals = tmb_high[metric].values
    l_vals = tmb_low[metric].values
    _, p = mannwhitneyu(h_vals, l_vals, alternative='two-sided')
    log(f"\n  {metric}:")
    log(f"    HIGH: mean={h_vals.mean():.0f}, median={np.median(h_vals):.0f}, "
        f"range={h_vals.min():.0f}-{h_vals.max():.0f}")
    log(f"    LOW:  mean={l_vals.mean():.0f}, median={np.median(l_vals):.0f}, "
        f"range={l_vals.min():.0f}-{l_vals.max():.0f}")
    log(f"    HIGH/LOW ratio: {h_vals.mean()/l_vals.mean():.2f}x")
    log(f"    Mann-Whitney p = {p:.2e}")

tmb.to_csv(os.path.join(OUTPUT_DIR, "HNSC_somatic_tmb_comparison.tsv"),
           sep='\t', index=False)
log(f"\n  Saved: HNSC_somatic_tmb_comparison.tsv")


# =============================================================================
# PHASE 2: Track A -- Gene-level burden tests
# =============================================================================
def run_gene_burden_test(somatic_df, variant_classes, label, high_patients, low_patients,
                         n_high, n_low, min_carriers=TRACK_A_MIN_CARRIERS):
    """Run Fisher's exact gene-level carrier enrichment test."""
    log(f"\n  --- Track A: {label} ---")

    # Filter to specified variant classes
    filtered = somatic_df[somatic_df['Variant_Classification'].isin(variant_classes)].copy()
    log(f"  Variants in class: {len(filtered):,}")

    # Gene-level carrier status per patient
    gene_patient_carrier = filtered.groupby(['Hugo_Symbol', 'patient_id']).size().reset_index(name='n_muts')
    gene_carriers = gene_patient_carrier.groupby('Hugo_Symbol')['patient_id'].apply(set).to_dict()

    results = []
    for gene, carriers in gene_carriers.items():
        n_total = len(carriers)
        if n_total < min_carriers:
            continue

        a = len(carriers & high_patients)   # carriers in HIGH
        b = n_high - a                       # non-carriers in HIGH
        c = len(carriers & low_patients)    # carriers in LOW
        d = n_low - c                        # non-carriers in LOW

        odds_ratio, pval = fisher_exact([[a, b], [c, d]], alternative='two-sided')

        results.append({
            'gene': gene,
            'n_carriers_high': a,
            'n_carriers_low': c,
            'n_carriers_total': n_total,
            'freq_high': a / n_high if n_high > 0 else 0,
            'freq_low': c / n_low if n_low > 0 else 0,
            'odds_ratio': odds_ratio,
            'pval': pval,
            'direction': 'HIGH' if a > c else ('LOW' if c > a else 'EQUAL'),
        })

    df = pd.DataFrame(results)
    if len(df) > 0:
        _, pval_bh, _, _ = multipletests(df['pval'], method='fdr_bh')
        df['pval_bh'] = pval_bh
        df = df.sort_values('pval')

    log(f"  Genes tested: {len(df)}")
    n_nom = (df['pval'] < 0.05).sum() if len(df) > 0 else 0
    n_bh = (df['pval_bh'] < BH_THRESHOLD).sum() if len(df) > 0 else 0
    log(f"  Nominal p < 0.05: {n_nom}")
    log(f"  BH < {BH_THRESHOLD}: {n_bh}")

    if len(df) > 0:
        log(f"\n  Top 20 genes by p-value ({label}):")
        log(f"  {'Gene':15s} {'HIGH':>5s} {'LOW':>5s} {'Total':>6s} "
            f"{'OR':>8s} {'p':>10s} {'BH':>10s} {'Dir':>6s}")
        log(f"  {'-'*15} {'-'*5} {'-'*5} {'-'*6} {'-'*8} {'-'*10} {'-'*10} {'-'*6}")
        for _, row in df.head(20).iterrows():
            log(f"  {row['gene']:15s} {row['n_carriers_high']:>5d} "
                f"{row['n_carriers_low']:>5d} {row['n_carriers_total']:>6d} "
                f"{row['odds_ratio']:>8.2f} {row['pval']:>10.4f} "
                f"{row['pval_bh']:>10.4f} {row['direction']:>6s}")

    return df


banner("PHASE 2: Track A -- Gene-Level Burden Tests")

# A1: Non-silent
track_a1 = run_gene_burden_test(somatic, NONSILENT_CLASSES, "Non-silent",
                                 high_patients, low_patients, n_high, n_low)
track_a1.to_csv(os.path.join(OUTPUT_DIR, "HNSC_somatic_track_A_nonsilent.tsv"),
                sep='\t', index=False)

# A2: Damaging / LOF
track_a2 = run_gene_burden_test(somatic, DAMAGING_CLASSES, "Damaging (LOF)",
                                 high_patients, low_patients, n_high, n_low)
track_a2.to_csv(os.path.join(OUTPUT_DIR, "HNSC_somatic_track_A_damaging.tsv"),
                sep='\t', index=False)

# A3: Silent (control)
track_a3 = run_gene_burden_test(somatic, SILENT_CLASSES, "Silent (control)",
                                 high_patients, low_patients, n_high, n_low)
track_a3.to_csv(os.path.join(OUTPUT_DIR, "HNSC_somatic_track_A_silent.tsv"),
                sep='\t', index=False)

# --- TMB-adjusted logistic regression for non-silent hits ---
banner("PHASE 2B: TMB-Adjusted Logistic Regression (Non-Silent)")

tmb_dict = tmb.set_index('patient_id')['TMB_all'].to_dict()

if len(track_a1) > 0:
    # Build patient-level data
    nonsilent = somatic[somatic['Variant_Classification'].isin(NONSILENT_CLASSES)].copy()
    gene_patient_carrier_ns = nonsilent.groupby(['Hugo_Symbol', 'patient_id']).size().reset_index(name='n_muts')
    gene_carriers_ns = gene_patient_carrier_ns.groupby('Hugo_Symbol')['patient_id'].apply(set).to_dict()

    lr_results = []
    for _, row in track_a1.iterrows():
        gene = row['gene']
        carriers = gene_carriers_ns.get(gene, set())

        # Build per-patient DataFrame for regression
        reg_rows = []
        for p in all_patients:
            reg_rows.append({
                'patient_id': p,
                'group_high': 1 if p in high_patients else 0,
                'carrier': 1 if p in carriers else 0,
                'log_tmb': np.log1p(tmb_dict.get(p, 0)),
            })
        reg_df = pd.DataFrame(reg_rows)

        # Logistic regression: carrier ~ group_high + log_tmb
        try:
            X = sm.add_constant(reg_df[['group_high', 'log_tmb']])
            y = reg_df['carrier']
            # Skip if no variation in carrier status
            if y.sum() == 0 or y.sum() == len(y):
                lr_results.append({
                    'gene': gene,
                    'fisher_p': row['pval'],
                    'fisher_or': row['odds_ratio'],
                    'lr_group_coef': np.nan,
                    'lr_group_p': np.nan,
                    'lr_tmb_coef': np.nan,
                    'lr_tmb_p': np.nan,
                    'lr_converged': False,
                })
                continue

            model = sm.Logit(y, X)
            result = model.fit(disp=0, maxiter=100, method='bfgs')

            lr_results.append({
                'gene': gene,
                'fisher_p': row['pval'],
                'fisher_or': row['odds_ratio'],
                'lr_group_coef': result.params.get('group_high', np.nan),
                'lr_group_p': result.pvalues.get('group_high', np.nan),
                'lr_tmb_coef': result.params.get('log_tmb', np.nan),
                'lr_tmb_p': result.pvalues.get('log_tmb', np.nan),
                'lr_converged': result.mle_retvals.get('converged', False) if hasattr(result, 'mle_retvals') else True,
            })
        except Exception:
            lr_results.append({
                'gene': gene,
                'fisher_p': row['pval'],
                'fisher_or': row['odds_ratio'],
                'lr_group_coef': np.nan,
                'lr_group_p': np.nan,
                'lr_tmb_coef': np.nan,
                'lr_tmb_p': np.nan,
                'lr_converged': False,
            })

    lr_df = pd.DataFrame(lr_results)

    # BH correction on LR group p-values
    valid_lr_p = lr_df['lr_group_p'].dropna()
    if len(valid_lr_p) > 1:
        _, lr_bh, _, _ = multipletests(valid_lr_p, method='fdr_bh')
        lr_df.loc[valid_lr_p.index, 'lr_group_p_bh'] = lr_bh
    else:
        lr_df['lr_group_p_bh'] = np.nan

    lr_df = lr_df.sort_values('lr_group_p')
    lr_df.to_csv(os.path.join(OUTPUT_DIR, "HNSC_somatic_track_A_tmb_adjusted.tsv"),
                 sep='\t', index=False)

    # Compare Fisher vs LR
    if len(lr_df) > 0:
        n_fisher_sig = (lr_df['fisher_p'] < 0.05).sum()
        n_lr_sig = (lr_df['lr_group_p'] < 0.05).sum()
        n_lr_bh = (lr_df['lr_group_p_bh'] < 0.10).sum()
        n_both = ((lr_df['fisher_p'] < 0.05) & (lr_df['lr_group_p'] < 0.05)).sum()
        n_fisher_only = ((lr_df['fisher_p'] < 0.05) & (lr_df['lr_group_p'] >= 0.05)).sum()

        log(f"\n  TMB adjustment comparison:")
        log(f"    Fisher p < 0.05: {n_fisher_sig}")
        log(f"    LR group p < 0.05: {n_lr_sig}")
        log(f"    LR group BH < 0.10: {n_lr_bh}")
        log(f"    Both Fisher + LR significant: {n_both}")
        log(f"    Fisher-only (TMB-driven): {n_fisher_only}")

        log(f"\n  Top 15 by TMB-adjusted p:")
        log(f"  {'Gene':15s} {'FisherP':>10s} {'FisherOR':>10s} "
            f"{'LR_grp_p':>10s} {'LR_BH':>10s} {'LR_grp_coef':>12s} {'LR_tmb_p':>10s}")
        for _, row in lr_df.head(15).iterrows():
            log(f"  {row['gene']:15s} {row['fisher_p']:>10.4f} "
                f"{row['fisher_or']:>10.2f} {row['lr_group_p']:>10.4f} "
                f"{row.get('lr_group_p_bh', np.nan):>10.4f} "
                f"{row['lr_group_coef']:>12.3f} {row['lr_tmb_p']:>10.4f}")

        # --- Focused depletion report (negative LR coef = depleted in HIGH) ---
        lr_sig = lr_df[lr_df['lr_group_p'] < 0.05].copy()
        depleted = lr_sig[lr_sig['lr_group_coef'] < 0].copy()
        enriched = lr_sig[lr_sig['lr_group_coef'] > 0].copy()

        log(f"\n  === ENRICHED in HIGH (TMB-adjusted, LR p < 0.05) ===")
        log(f"  Genes: {len(enriched)}")
        if len(enriched) > 0:
            log(f"  {'Gene':15s} {'HIGH':>5s} {'LOW':>5s} {'LR_p':>10s} {'LR_BH':>10s} {'Coef':>8s}")
            for _, row in enriched.iterrows():
                gene = row['gene']
                # Get carrier counts from Track A1
                a1_row = track_a1[track_a1['gene'] == gene]
                n_h = a1_row['n_carriers_high'].values[0] if len(a1_row) > 0 else '?'
                n_l = a1_row['n_carriers_low'].values[0] if len(a1_row) > 0 else '?'
                log(f"  {gene:15s} {str(n_h):>5s} {str(n_l):>5s} "
                    f"{row['lr_group_p']:>10.4f} {row.get('lr_group_p_bh', np.nan):>10.4f} "
                    f"{row['lr_group_coef']:>8.3f}")

        log(f"\n  === DEPLETED in HIGH (TMB-adjusted, LR p < 0.05) ===")
        log(f"  Genes: {len(depleted)}")
        log(f"  (These genes are mutated LESS in HIGH after controlling for TMB,")
        log(f"   suggesting their loss may be incompatible with high APOBEC activity)")
        if len(depleted) > 0:
            log(f"  {'Gene':15s} {'HIGH':>5s} {'LOW':>5s} {'LR_p':>10s} {'LR_BH':>10s} {'Coef':>8s}")
            for _, row in depleted.iterrows():
                gene = row['gene']
                a1_row = track_a1[track_a1['gene'] == gene]
                n_h = a1_row['n_carriers_high'].values[0] if len(a1_row) > 0 else '?'
                n_l = a1_row['n_carriers_low'].values[0] if len(a1_row) > 0 else '?'
                log(f"  {gene:15s} {str(n_h):>5s} {str(n_l):>5s} "
                    f"{row['lr_group_p']:>10.4f} {row.get('lr_group_p_bh', np.nan):>10.4f} "
                    f"{row['lr_group_coef']:>8.3f}")


# =============================================================================
# PHASE 2C: Frequency Tier + Group Bias Binning
# =============================================================================
# Parallels the single-cell tier analysis (Tier2B_SNP_Pattern_Analysis.py)
# so that bulk and SC results use comparable vocabulary.
#
# Frequency tiers (based on HIGH group carrier rate):
#   Highly retained: >=75% of HIGH (>=40/53)
#   Common:          >=50% of HIGH (>=27/53)
#   Moderate:        >=25% of HIGH (>=14/53)
#   Uncommon:        >=10% of HIGH (>=6/53)
#
# Group bias (fraction of total carriers in HIGH):
#   HIGH-exclusive:    n_low == 0
#   HIGH-predominant:  >=75% carriers in HIGH
#   Balanced:          25-75% carriers in HIGH
#   LOW-predominant:   <25% carriers in HIGH
#   LOW-exclusive:     n_high == 0
# =============================================================================
banner("PHASE 2C: Frequency Tier + Group Bias Binning")

if len(track_a1) > 0:
    # Frequency tier thresholds (based on n_high)
    TIER_THRESHOLDS = [
        ('Highly_retained', int(np.ceil(n_high * 0.75))),
        ('Common',          int(np.ceil(n_high * 0.50))),
        ('Moderate',        int(np.ceil(n_high * 0.25))),
        ('Uncommon',        int(np.ceil(n_high * 0.10))),
    ]

    log(f"  Frequency tier thresholds (n_high={n_high}):")
    for tier_name, thresh in TIER_THRESHOLDS:
        log(f"    {tier_name}: >= {thresh} carriers in HIGH (>= {100*thresh/n_high:.0f}%)")

    def assign_frequency_tier(n_carriers_high):
        for tier_name, thresh in TIER_THRESHOLDS:
            if n_carriers_high >= thresh:
                return tier_name
        return 'Below_threshold'

    def assign_group_bias(n_h, n_l):
        total = n_h + n_l
        if total == 0:
            return 'No_carriers'
        if n_l == 0:
            return 'HIGH-exclusive'
        if n_h == 0:
            return 'LOW-exclusive'
        frac_high = n_h / total
        if frac_high >= 0.75:
            return 'HIGH-predominant'
        elif frac_high <= 0.25:
            return 'LOW-predominant'
        else:
            return 'Balanced'

    track_a1['frequency_tier'] = track_a1['n_carriers_high'].apply(assign_frequency_tier)
    track_a1['group_bias'] = track_a1.apply(
        lambda r: assign_group_bias(r['n_carriers_high'], r['n_carriers_low']), axis=1
    )

    # Report tier distribution
    log(f"\n  Frequency tier distribution (non-silent genes):")
    for tier_name, _ in TIER_THRESHOLDS:
        tier_genes = track_a1[track_a1['frequency_tier'] == tier_name]
        if len(tier_genes) > 0:
            log(f"\n  --- {tier_name} (n={len(tier_genes)} genes) ---")

            # Group bias breakdown within this tier
            bias_counts = tier_genes['group_bias'].value_counts()
            for bias, count in bias_counts.items():
                log(f"    {bias}: {count} genes")

            # Show HIGH-exclusive and HIGH-predominant genes in this tier
            interesting = tier_genes[tier_genes['group_bias'].isin(
                ['HIGH-exclusive', 'HIGH-predominant']
            )].sort_values('n_carriers_high', ascending=False)
            if len(interesting) > 0:
                log(f"\n    HIGH-biased genes in {tier_name} tier:")
                log(f"    {'Gene':15s} {'HIGH':>5s} {'LOW':>5s} {'Bias':>18s} "
                    f"{'FisherP':>10s} {'LR_P':>10s}")
                for _, row in interesting.head(15).iterrows():
                    lr_p_val = ''
                    if len(lr_df) > 0:
                        lr_match = lr_df[lr_df['gene'] == row['gene']]
                        if len(lr_match) > 0:
                            lr_p_val = f"{lr_match.iloc[0]['lr_group_p']:.4f}"
                    log(f"    {row['gene']:15s} {row['n_carriers_high']:>5d} "
                        f"{row['n_carriers_low']:>5d} {row['group_bias']:>18s} "
                        f"{row['pval']:>10.4f} {lr_p_val:>10s}")

            # Show LOW-exclusive and LOW-predominant genes in this tier
            low_interesting = tier_genes[tier_genes['group_bias'].isin(
                ['LOW-exclusive', 'LOW-predominant']
            )].sort_values('n_carriers_low', ascending=False)
            if len(low_interesting) > 0:
                log(f"\n    LOW-biased genes in {tier_name} tier:")
                log(f"    {'Gene':15s} {'HIGH':>5s} {'LOW':>5s} {'Bias':>18s} "
                    f"{'FisherP':>10s} {'LR_P':>10s}")
                for _, row in low_interesting.head(15).iterrows():
                    lr_p_val = ''
                    if len(lr_df) > 0:
                        lr_match = lr_df[lr_df['gene'] == row['gene']]
                        if len(lr_match) > 0:
                            lr_p_val = f"{lr_match.iloc[0]['lr_group_p']:.4f}"
                    log(f"    {row['gene']:15s} {row['n_carriers_high']:>5d} "
                        f"{row['n_carriers_low']:>5d} {row['group_bias']:>18s} "
                        f"{row['pval']:>10.4f} {lr_p_val:>10s}")

    # Re-save Track A1 with tier and bias columns
    track_a1.to_csv(os.path.join(OUTPUT_DIR, "HNSC_somatic_track_A_nonsilent.tsv"),
                    sep='\t', index=False)
    log(f"\n  Updated Track A1 with frequency_tier and group_bias columns")

    # Summary cross-tab
    log(f"\n  Cross-tab: Frequency Tier x Group Bias")
    tier_order = ['Highly_retained', 'Common', 'Moderate', 'Uncommon']
    bias_order = ['HIGH-exclusive', 'HIGH-predominant', 'Balanced', 'LOW-predominant', 'LOW-exclusive']
    ct = pd.crosstab(track_a1['frequency_tier'], track_a1['group_bias'])
    ct = ct.reindex(index=[t for t in tier_order if t in ct.index],
                    columns=[b for b in bias_order if b in ct.columns], fill_value=0)
    log(f"\n  {'Tier':20s} " + " ".join(f"{b:>18s}" for b in ct.columns))
    for tier in ct.index:
        vals = " ".join(f"{ct.loc[tier, b]:>18d}" for b in ct.columns)
        log(f"  {tier:20s} {vals}")

    # -----------------------------------------------------------------
    # PHASE 2C (cont.): Per-tier gene set cross-referencing + KEGG
    # -----------------------------------------------------------------
    # For each frequency tier, check overlap with reference gene sets
    # and run KEGG enrichment. Parallels SC Analyze_SNP_Tier_Genes.py.
    # -----------------------------------------------------------------
    banner("PHASE 2C+: Per-Tier Gene Set Cross-Reference + KEGG")

    # Reference sets to cross-reference (top-level only, skip individual communities)
    TIER_XREF_SETS = {
        'Harris_A3_all': harris_all,
        'Harris_A3B_specific': harris_a3b,
        'DDR_genes': ddr_set,
        'Chromatin_remodelers': chromatin_set,
        'Network_all_communities': all_network_genes,
    }
    # Also include the A3B community (Community_4) individually since it is
    # the focal community for the paper
    if 'Community_4' in community_gene_sets:
        TIER_XREF_SETS['Community_4_A3B'] = community_gene_sets['Community_4']

    tier_xref_rows = []

    for tier_name, _ in TIER_THRESHOLDS:
        tier_genes_set = set(
            track_a1[track_a1['frequency_tier'] == tier_name]['gene'].values
        )
        n_tier = len(tier_genes_set)
        if n_tier == 0:
            continue

        log(f"\n  --- {tier_name} tier ({n_tier} genes) ---")

        for gs_name, gs_genes in TIER_XREF_SETS.items():
            overlap = tier_genes_set & gs_genes
            n_overlap = len(overlap)
            tier_xref_rows.append({
                'frequency_tier': tier_name,
                'gene_set': gs_name,
                'n_tier_genes': n_tier,
                'n_set_genes': len(gs_genes),
                'n_overlap': n_overlap,
                'overlap_genes': ';'.join(sorted(overlap)) if overlap else '',
            })
            if n_overlap > 0:
                log(f"    {gs_name}: {n_overlap} overlap "
                    f"({', '.join(sorted(overlap))})")

        # KEGG enrichment for this tier (if enough genes)
        TIER_KEGG_MIN_GENES = 5
        if n_tier >= TIER_KEGG_MIN_GENES:
            try:
                import gseapy as gp
                tier_gene_list = sorted(list(tier_genes_set))
                enrichr_res = gp.enrichr(
                    gene_list=tier_gene_list,
                    gene_sets=['KEGG_2021_Human'],
                    organism='human',
                    outdir=os.path.join(OUTPUT_DIR, f"enrichr_kegg_{tier_name}"),
                    no_plot=True,
                )
                kegg_tier = enrichr_res.results.sort_values('Adjusted P-value')
                n_sig = (kegg_tier['Adjusted P-value'] < 0.05).sum()
                log(f"    KEGG ({tier_name}): {len(kegg_tier)} pathways tested, "
                    f"{n_sig} adj p < 0.05")
                if n_sig > 0:
                    for _, krow in kegg_tier[kegg_tier['Adjusted P-value'] < 0.05].head(5).iterrows():
                        log(f"      {krow['Term'][:55]:55s} "
                            f"{krow['Overlap']:>8s}  adj p={krow['Adjusted P-value']:.4f}")
                kegg_tier['frequency_tier'] = tier_name
                kegg_tier.to_csv(
                    os.path.join(OUTPUT_DIR, f"HNSC_somatic_tier_kegg_{tier_name}.tsv"),
                    sep='\t', index=False,
                )
            except ImportError:
                log(f"    KEGG ({tier_name}): skipped (gseapy not installed)")
            except Exception as e:
                log(f"    KEGG ({tier_name}): failed ({e})")
        else:
            log(f"    KEGG ({tier_name}): skipped (only {n_tier} genes, "
                f"need >= {TIER_KEGG_MIN_GENES})")

    # Save cross-reference table
    if tier_xref_rows:
        xref_df = pd.DataFrame(tier_xref_rows)
        xref_path = os.path.join(OUTPUT_DIR, "HNSC_somatic_tier_geneset_xref.tsv")
        xref_df.to_csv(xref_path, sep='\t', index=False)
        log(f"\n  Saved: {os.path.basename(xref_path)}")

else:
    log("  No Track A1 results for binning")


# =============================================================================
# PHASE 3: Track B -- A3-specific somatic catalog
# =============================================================================
banner("PHASE 3: Track B -- A3-Specific Somatic Catalog")

a3_somatic = somatic[somatic['Hugo_Symbol'].isin(A3_GENE_SYMBOLS)].copy()
log(f"  PASS somatic variants in A3 genes: {len(a3_somatic)}")

if len(a3_somatic) > 0:
    log(f"\n  Per-gene breakdown:")
    for gene in A3_GENE_SYMBOLS:
        gdf = a3_somatic[a3_somatic['Hugo_Symbol'] == gene]
        if len(gdf) > 0:
            n_h = gdf[gdf['group'] == 'HIGH']['patient_id'].nunique()
            n_l = gdf[gdf['group'] == 'LOW']['patient_id'].nunique()
            log(f"    {gene}: {len(gdf)} variants ({n_h} HIGH patients, {n_l} LOW patients)")
            for _, row in gdf.iterrows():
                log(f"      {row['patient_id']} | {row['Variant_Classification']} | "
                    f"{row.get('HGVSp_Short', '')} | {row['group']}")

    a3_somatic.to_csv(os.path.join(OUTPUT_DIR, "HNSC_somatic_track_B_a3_variants.tsv"),
                      sep='\t', index=False)
else:
    log("  No somatic variants found in A3 genes")


# =============================================================================
# PHASE 4: APOBEC trinucleotide context from per-sample MAFs
# =============================================================================
banner("PHASE 4: APOBEC Trinucleotide Context Analysis")

log(f"  Reading per-sample MAF files from: {PER_SAMPLE_MAF_DIR}")

# Find MAF files for our patients
maf_gz_files = sorted(glob.glob(os.path.join(PER_SAMPLE_MAF_DIR, "*.maf.gz")))
log(f"  Total MAF.gz files in directory: {len(maf_gz_files)}")

# Map patient ID -> per-sample MAF file
patient_to_maf = {}
for maf_file in maf_gz_files:
    try:
        # Read just first few rows to get barcode
        peek = pd.read_csv(maf_file, sep='\t', comment='#', nrows=2,
                           usecols=['Tumor_Sample_Barcode'], low_memory=False)
        if len(peek) > 0:
            bc = peek['Tumor_Sample_Barcode'].iloc[0]
            pid = bc[:12]
            if pid in all_patients:
                patient_to_maf[pid] = maf_file
    except Exception:
        continue

log(f"  Matched MAF files for our patients: {len(patient_to_maf)}")
missing_maf = all_patients - set(patient_to_maf.keys())
if missing_maf:
    log(f"  WARNING: {len(missing_maf)} patients without per-sample MAF files")

# Read CONTEXT for each patient, classify APOBEC vs non-APOBEC
gene_context_counts = {}  # gene -> {'apobec': count, 'non_apobec': count, 'total': count}

n_processed = 0
n_total_pass_snps = 0

for pid, maf_file in patient_to_maf.items():
    try:
        df = pd.read_csv(maf_file, sep='\t', comment='#',
                         usecols=CONTEXT_MAF_COLS, low_memory=False)

        # PASS SNPs only
        mask = (
            (df['Variant_Type'] == 'SNP') &
            (df['FILTER'] == 'PASS') &
            (df['CONTEXT'].str.len() == 11)
        )
        pass_snps = df[mask]
        n_total_pass_snps += len(pass_snps)

        for _, row in pass_snps.iterrows():
            gene = row['Hugo_Symbol']
            ctx = row['CONTEXT']
            ref = row['Reference_Allele']
            alt = row['Tumor_Seq_Allele2']

            if gene not in gene_context_counts:
                gene_context_counts[gene] = {'apobec': 0, 'non_apobec': 0, 'total': 0}

            gene_context_counts[gene]['total'] += 1

            # Extract trinucleotide: positions [4],[5],[6] (0-indexed)
            five_prime = ctx[4]
            context_ref = ctx[5]
            three_prime = ctx[6]

            # APOBEC context: TCW motif (5' base is T, ref is C)
            # After pyrimidine convention: C>T or C>G at TC context
            is_apobec = False
            if context_ref == ref:
                # Check both orientations
                if ref == 'C' and five_prime == 'T' and alt in ('T', 'G'):
                    is_apobec = True
                elif ref == 'G' and three_prime == 'A' and alt in ('A', 'C'):
                    # Reverse complement: G at position with 3'=A means
                    # on opposite strand: 5'=T, ref=C
                    is_apobec = True

            if is_apobec:
                gene_context_counts[gene]['apobec'] += 1
            else:
                gene_context_counts[gene]['non_apobec'] += 1

        n_processed += 1

    except Exception as e:
        log(f"  WARNING: Failed to read {os.path.basename(maf_file)}: {e}")

log(f"\n  Patients processed: {n_processed}")
log(f"  Total PASS SNPs across all patients: {n_total_pass_snps:,}")
log(f"  Genes with context data: {len(gene_context_counts)}")

# Build context summary table
context_rows = []
for gene, counts in gene_context_counts.items():
    total = counts['total']
    apobec = counts['apobec']
    frac = apobec / total if total > 0 else 0
    context_rows.append({
        'gene': gene,
        'n_mutations_total': total,
        'n_apobec_context': apobec,
        'n_non_apobec_context': counts['non_apobec'],
        'frac_apobec_context': frac,
    })

context_df = pd.DataFrame(context_rows).sort_values('n_mutations_total', ascending=False)
context_df.to_csv(os.path.join(OUTPUT_DIR, "HNSC_somatic_apobec_context_by_gene.tsv"),
                  sep='\t', index=False)

# Merge context info into Track A1 results
if len(track_a1) > 0 and len(context_df) > 0:
    ctx_dict = context_df.set_index('gene')['frac_apobec_context'].to_dict()
    track_a1['frac_apobec_context'] = track_a1['gene'].map(ctx_dict)
    # Re-save with context column
    track_a1.to_csv(os.path.join(OUTPUT_DIR, "HNSC_somatic_track_A_nonsilent.tsv"),
                    sep='\t', index=False)
    log(f"  Updated Track A1 with frac_apobec_context column")

# Report overall APOBEC context rate
total_apobec = sum(c['apobec'] for c in gene_context_counts.values())
total_all = sum(c['total'] for c in gene_context_counts.values())
log(f"\n  Overall APOBEC context rate: {total_apobec:,}/{total_all:,} "
    f"({100*total_apobec/total_all:.1f}%)" if total_all > 0 else "")


# =============================================================================
# PHASE 5: Track C -- Gene set enrichment + KEGG
# =============================================================================
# NOTE: Reference gene sets (harris_all, harris_a3b, community_gene_sets,
# ddr_set, chromatin_set, REFERENCE_GENE_SETS) were loaded in Phase 0B.
# =============================================================================
banner("PHASE 5: Track C -- Gene Set and Pathway Enrichment")

# --- Hypergeometric gene set tests ---
banner("PHASE 5A: Hypergeometric Gene Set Tests")

# Use Track A1 (non-silent) nominal hits
if len(track_a1) > 0:
    nominal_hits = set(track_a1[track_a1['pval'] < 0.05]['gene'].values)
    all_tested = set(track_a1['gene'].values)
    n_hits = len(nominal_hits)
    n_tested = len(all_tested)
    log(f"  Nominal hits (p < 0.05): {n_hits} / {n_tested} genes tested")

    from scipy.stats import hypergeom

    enrichment_results = []
    for gs_name, gs_genes in REFERENCE_GENE_SETS.items():
        # Genes in gene set that were actually tested
        gs_tested = gs_genes & all_tested
        gs_hits = gs_genes & nominal_hits
        n_gs_tested = len(gs_tested)
        n_gs_hits = len(gs_hits)

        if n_gs_tested == 0:
            continue

        # Hypergeometric test
        # P(X >= n_gs_hits) where X ~ Hypergeom(N=n_tested, K=n_hits, n=n_gs_tested)
        pval = hypergeom.sf(n_gs_hits - 1, n_tested, n_hits, n_gs_tested)

        enrichment_results.append({
            'gene_set': gs_name,
            'n_in_set': len(gs_genes),
            'n_tested': n_gs_tested,
            'n_hits': n_gs_hits,
            'expected': n_hits * n_gs_tested / n_tested if n_tested > 0 else 0,
            'fold_enrichment': (n_gs_hits / n_gs_tested) / (n_hits / n_tested)
                                if n_hits > 0 and n_gs_tested > 0 else 0,
            'pval': pval,
            'hit_genes': ';'.join(sorted(gs_hits)) if gs_hits else '',
        })

    enrich_df = pd.DataFrame(enrichment_results).sort_values('pval')
    if len(enrich_df) > 0:
        _, pval_bh, _, _ = multipletests(enrich_df['pval'], method='fdr_bh')
        enrich_df['pval_bh'] = pval_bh

    enrich_df.to_csv(os.path.join(OUTPUT_DIR, "HNSC_somatic_track_C_geneset_enrichment.tsv"),
                     sep='\t', index=False)

    log(f"\n  Gene set enrichment results:")
    log(f"  {'Gene Set':30s} {'Tested':>7s} {'Hits':>5s} {'Exp':>6s} "
        f"{'Fold':>6s} {'p':>10s} {'BH':>10s}")
    log(f"  {'-'*30} {'-'*7} {'-'*5} {'-'*6} {'-'*6} {'-'*10} {'-'*10}")
    for _, row in enrich_df.iterrows():
        log(f"  {row['gene_set']:30s} {row['n_tested']:>7d} {row['n_hits']:>5d} "
            f"{row['expected']:>6.1f} {row['fold_enrichment']:>6.2f} "
            f"{row['pval']:>10.4f} {row['pval_bh']:>10.4f}")
else:
    log("  No Track A1 results to run enrichment on")

# --- KEGG enrichment via gseapy/Enrichr ---
banner("PHASE 5B: KEGG Pathway Enrichment (Enrichr)")

if len(track_a1) > 0 and n_hits > 0:
    try:
        import gseapy as gp

        gene_list = sorted(list(nominal_hits))
        log(f"  Submitting {len(gene_list)} genes to Enrichr...")

        enrichr_result = gp.enrichr(
            gene_list=gene_list,
            gene_sets=['KEGG_2021_Human'],
            organism='human',
            outdir=os.path.join(OUTPUT_DIR, "enrichr_kegg"),
            no_plot=True,
        )

        kegg_df = enrichr_result.results
        kegg_df = kegg_df.sort_values('Adjusted P-value')
        kegg_df.to_csv(os.path.join(OUTPUT_DIR, "HNSC_somatic_track_C_kegg_enrichr.tsv"),
                       sep='\t', index=False)

        n_kegg_sig = (kegg_df['Adjusted P-value'] < 0.05).sum()
        log(f"  KEGG pathways tested: {len(kegg_df)}")
        log(f"  KEGG pathways (adj p < 0.05): {n_kegg_sig}")

        if len(kegg_df) > 0:
            log(f"\n  Top 15 KEGG pathways:")
            log(f"  {'Pathway':50s} {'Overlap':>10s} {'AdjP':>10s}")
            log(f"  {'-'*50} {'-'*10} {'-'*10}")
            for _, row in kegg_df.head(15).iterrows():
                log(f"  {row['Term'][:50]:50s} {row['Overlap']:>10s} "
                    f"{row['Adjusted P-value']:>10.4f}")

    except ImportError:
        log("  WARNING: gseapy not installed, skipping KEGG enrichment")
        log("  Install with: pip install gseapy")
    except Exception as e:
        log(f"  WARNING: Enrichr query failed: {e}")
        log("  This may be due to network restrictions on the HPC")
        log("  You can run the Enrichr step separately if needed")
else:
    log("  No nominal hits to submit to Enrichr")


# =============================================================================
# PHASE 6: Summary report
# =============================================================================
banner("PHASE 6: Summary Report")

log(f"""
  HNSC SOMATIC ENRICHMENT ANALYSIS COMPLETE
  ==========================================

  Groups: {n_high} HIGH + {n_low} LOW = {n_high + n_low} patients
  Source: Step03 output from updated network pipeline
  MAF: {MAF_PATH}

  PASS somatic variants: {len(somatic):,}
  Genes tested (non-silent): {len(track_a1)}
  Genes tested (damaging):   {len(track_a2)}
  Genes tested (silent):     {len(track_a3)}

  Track A1 (non-silent) nominal p < 0.05: {(track_a1['pval'] < 0.05).sum() if len(track_a1) > 0 else 0}
  Track A1 BH < {BH_THRESHOLD}: {(track_a1['pval_bh'] < BH_THRESHOLD).sum() if len(track_a1) > 0 else 0}
  Track A2 (damaging) nominal p < 0.05: {(track_a2['pval'] < 0.05).sum() if len(track_a2) > 0 else 0}
  Track A3 (silent) nominal p < 0.05: {(track_a3['pval'] < 0.05).sum() if len(track_a3) > 0 else 0}

  A3 somatic variants found: {len(a3_somatic)}
  APOBEC context rate: {100*total_apobec/total_all:.1f}% (of PASS SNPs)

  Output directory: {OUTPUT_DIR}
""")

# Save report
report_path = os.path.join(OUTPUT_DIR, "HNSC_somatic_enrichment_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Report: {report_path}")

banner("PIPELINE COMPLETE")
