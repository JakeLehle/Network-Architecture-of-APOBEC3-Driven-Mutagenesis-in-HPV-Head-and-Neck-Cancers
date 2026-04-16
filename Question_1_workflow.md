# Question 1 Workflow

## Question 1: What is the fundamental relationship between A3 expression and observed A3-induced mutations in HNSCC?

### Rationale

APOBEC3 (A3) cytidine deaminases with access to the nuclear compartment (A3A, A3B, A3C, A3H) are capable of editing genomic DNA, producing characteristic C>T and C>G mutations captured by the COSMIC SBS2 mutational signature. However, the presence of A3 enzyme expression alone does not guarantee that A3-driven mutagenesis will occur in a given tumor. This question establishes the foundational observation of the paper: **A3A and A3B expression is necessary but not sufficient for SBS2 mutagenesis in HNSCC**, implying the existence of unknown cofactors that regulate A3 enzymatic activity on genomic DNA. Additionally, this question dissects the individual and combinatorial contributions of A3A and A3B, establishing that A3B provides a baseline mutagenic state that A3A amplifies.

### Key Results (Revised, April 2026)

| Metric | Value |
|--------|-------|
| HNSC tumors matched (RNA-seq + WES) | 502 |
| SBS2 range | 0 – 728 mutations |
| SBS2 median | 16 |
| Tumors with SBS2 > 0 | 326 (64.9%) |
| Coral region (A3 + SBS2 HIGH) | 251 |
| Cream region (A3 + SBS2 LOW) | 251 |
| Teal region (no A3 + high SBS2) | 0 |
| Low-A3 outliers (A3A+A3B < 1.0, SBS2 > median) | 3 (0.6%) |
| Median-split quadrant medians | 10 → 15.5 → 18.5 → 22.0 |
| Both-high vs neither p-value | 1.54 × 10⁻⁴ (Wilcoxon) |
| SBS signature source | SigProfilerAssignment, COSMIC v3.4, exome mode |

### Data Sources

| Source | Description |
|--------|-------------|
| [GDC Data Portal](https://portal.gdc.cancer.gov/) | TCGA RNA-seq gene expression (STAR - Counts, FPKM-UQ) across 33 cancer types |
| [GDC Data Portal](https://portal.gdc.cancer.gov/) | TCGA MuTect2 Annotated Somatic Mutation MAFs across 33 cancer types (controlled access) |
| SigProfilerAssignment | COSMIC v3.4 SBS signature weights computed from PASS-filtered somatic SNVs |

### Directory Structure

```
scripts/TCGA/
├── DOWNLOAD_TCGA_COUNT_DATA.sh              # SLURM Job 1: RNA-seq data acquisition
├── Download_RNA-seq_Counts_TCGA.R           #   Step 1 — Download from GDC
├── Organize_RNA-seq_Counts_TCGA.R           #   Step 2 — Organize & build manifest
├── Master_TCGA_RNA-seq_Counts_Table.R       #   Step 3 — Build master expression tables
│
├── ANALYZE_TCGA_COUNT_DATA.sh               # SLURM Job 2: Original analysis (Steps 4-8)
├── Prep_mutation_analysis_files.R           #   Step 4 — Match A3 expression ↔ SBS signatures
├── Patient_Level_HNSCC_TCGA_A3s_vs_SBS2.R  #   Step 5 (ORIGINAL) — 6-panel Fig 1a
├── Diagnostic_A3C_A3H_bystander.R          #   Step 6 (ORIGINAL) — ROC/AUC (Fig 1b)
├── Patient_Level_HNSCC_TCGA_3D_A3s.R       #   Step 7 (ORIGINAL) — 3D scatter (Fig 1c)
├── A3A_A3B_additive_SBS2.R                 #   Step 8 (ORIGINAL) — Additive model (Fig 1d-f)
│
├── Step05_Revised_HNSC_A3_vs_SBS2.py       #   Step 5 (REVISED) — 3-panel Fig 1 + supplemental
│                                            #     Uses new SigProfiler v3.4 SBS2 counts
│                                            #     Replaces Steps 5-8 R scripts for figure gen
│
├── VCF/                                     # Branch A: Somatic mutation & SBS signature pipeline
│   ├── Diagnostic_TCGA_SNV_Availability.R   #   Step 0 — Check GDC data availability
│   ├── Download_MuTect2_VCFs_TCGA.R         #   Step 1 — Download MuTect2 MAFs (10,939 samples)
│   ├── Download_Pindel_VCFs_TCGA.R          #   Step 2 — Download Pindel MAFs
│   ├── Organize_SNV_VCFs_TCGA.R             #   Step 3 — Build unified manifest
│   ├── Consolidate_MAFs_TCGA.R              #   Step 4 — Consolidate to pan-cancer master MAF
│   ├── Run_SigProfiler.py                   #   Step 5 — Build SBS96 matrix + COSMIC assignment
│   └── SIG_PROFILER_SCRIPT.sh               #   SLURM batch script for Steps 4-5
│
└── TROUBLESHOOTING/
    ├── Diagnostic_Check_SBS2_Weight_Normalization.py   # Verified original SBS2 weights are
    │                                                    # absolute counts (CV=4.22), not ratios
    ├── Diagnostic_Compare_SBS_Weight_Sources.py        # Compared original vs new SigProfiler
    │                                                    # outputs (Spearman rho=0.90 for HNSC)
    ├── Diagnostic_Verify_Barcode_Matching.py           # Audited all 504 RNA-seq↔WES barcode
    │                                                    # pairs; verified 3 low-A3 outliers
    ├── Diagnostic_Figure1_Text_Numbers.py              # Extracted all text-ready numbers
    ├── TCGA_troubleshooting.R
    ├── Compare_A3_Expression_Sources.R
    ├── FILTER_MASTER_TCGA_TSV.R
    ├── Global_TCGA_A3s_vs_SBS2.R
    ├── Patient_Level_HNSCC_TCGA_Network_Analysis.R
    ├── A3_Contribution_to_SBS2.R
    └── Diagnose_nonlinear_A3_Relationships.R
```

### Pipeline Overview

The analysis has two branches that converge in the revised Step 5:

```
 BRANCH A: RNA-seq Expression (SLURM Jobs 1-2)
 ──────────────────────────────────────────────
   Step 1-3: Download & build TCGA_master_FPKM_UQ.tsv
   Step 4: Extract A3 expression, match to mutation data
       │
       └──→ TCGA_master_FPKM_UQ.tsv (11,505 samples)
            TCGA_sample_metadata_final.tsv (Project_ID)
            Mutation_Table_Tumors_TCGA.tsv (barcode crosswalk)


 BRANCH B: Somatic Mutation SBS Signatures (scripts/TCGA/VCF/)
 ──────────────────────────────────────────────────────────────
   Step 0: Diagnostic — check GDC availability
   Step 1: Download MuTect2 MAFs (10,939 samples, 33 cancer types)
   Step 2: Download Pindel MAFs (for XRCC4 analysis, not used here)
   Step 3: Organize & verify files, build unified manifest
   Step 4: Consolidate per-sample MAFs → pan-cancer master MAF
   Step 5: Build SBS96 matrix from CONTEXT column → SigProfilerAssignment
       │
       └──→ TCGA_SBS_signature_counts.tsv (absolute mutation counts per sig)
            TCGA_SBS_signature_weights.tsv (normalized fractions)
            TCGA_MuTect2_master_manifest.tsv (Entity_ID → Cancer_Type)


 CONVERGENCE: Revised Step 5 (Step05_Revised_HNSC_A3_vs_SBS2.py)
 ────────────────────────────────────────────────────────────────
   Inputs from both branches:
     TCGA_master_FPKM_UQ.tsv ← Branch A
     TCGA_SBS_signature_counts.tsv ← Branch B
     TCGA_sample_metadata_final.tsv ← Branch A (cancer type)
     Mutation_Table_Tumors_TCGA.tsv ← Branch A (barcode crosswalk)
     TCGA_MuTect2_master_manifest.tsv ← Branch B (HNSC WES IDs)
       │
       ▼
   1. Load full expression matrix, extract A3 genes
   2. Filter to HNSC primary tumors (522 RNA-seq samples)
   3. Load & transpose SigProfiler counts (512 HNSC WES samples)
   4. Layered barcode matching:
       a. Direct WES crosswalk (426 matches)
       b. Case_ID with sample-type constraint (76 matches)
       c. Patient-level deduplication → 502 final matches
   5. Build matched table with A3A, A3B, SBS2
   6. Generate Figure 1 panels:
       │
       ├── Panel 1a: A3A+A3B vs SBS2 scatter
       │     Polygon background regions (teal/coral/cream)
       │     Mint y-axis line, broken x-axis
       │
       ├── Panel 1b: A3A vs A3B colored by SBS2
       │     Depth-sorted, broken x-axis, magma colormap
       │
       ├── Panel 1c: Boxplot + 2×2 heatmap
       │     4 quadrants by median A3A/A3B
       │     Medians: 10 → 15.5 → 18.5 → 22.0
       │
       └── Supplemental: Low-A3 zoom
             4 lowest-A3 tumors in high-SBS2 group labeled
       │
       ▼
   FIGURE 1
```

### Environments

| Environment | Scripts | Purpose |
|-------------|---------|---------|
| `RNA-seq_NovoGene` | Steps 1-4 (R), VCF Steps 0-4 (R) | TCGA data download, processing |
| `SComatic` | VCF Step 5 (Python) | SigProfilerAssignment |
| `NETWORK` | Revised Step 5 (Python) | Figure generation, diagnostics |

---

### Revised Step 5: Figure 1 Generation (CURRENT)

**Script:** `Step05_Revised_HNSC_A3_vs_SBS2.py`

The primary Figure 1 generation script. Replaces the original R scripts (Steps 5-8) with a single Python script that loads new SigProfiler v3.4 data and generates the simplified 3-panel figure.

**Dependencies:** `pandas`, `numpy`, `matplotlib`, `scipy`

**Input:**
- `TCGA_master_FPKM_UQ.tsv` — full pan-cancer expression matrix (from Step 3)
- `TCGA_sample_metadata_final.tsv` — Project_ID for cancer type filtering
- `Mutation_Table_Tumors_TCGA.tsv` — barcode crosswalk (RNA-seq ↔ WES)
- `TCGA_SBS_signature_counts.tsv` — SigProfiler v3.4 absolute SBS counts (from VCF Step 5)
- `TCGA_MuTect2_master_manifest.tsv` — HNSC WES Entity_IDs

**Output (→ `data/FIG_1/`):**
- `HNSC_A3_SBS2_matched_v2.tsv` — 502 matched tumors with all A3 expression + SBS2 + region + quadrant
- `FIGURE_1_PANELS/Panel_1a_A3sum_vs_SBS2.pdf/.png`
- `FIGURE_1_PANELS/Panel_1b_A3A_vs_A3B_SBS2.pdf/.png`
- `FIGURE_1_PANELS/Panel_1c_Boxplot_Heatmap.pdf/.png`
- `FIGURE_1_PANELS/Supplemental_Low_A3_High_SBS2_Zoom.pdf/.png`
- `TROUBLESHOOTING/figure1_final_pipeline_report.txt`

**Barcode matching strategy:**
1. **Direct crosswalk** (primary): Uses the `Mutation_Signature__File_Orginal_Entity_ID` column in the original mutation file to map RNA-seq barcodes to WES barcodes
2. **Case_ID matching** (fallback): 12-character patient ID match, constrained to same sample type code (both primary tumor = 01, both metastatic = 06) to prevent cross-matching
3. **Patient deduplication**: If a patient has both a direct and Case_ID match, only the direct match is kept

**Barcode verification (Diagnostic_Verify_Barcode_Matching.py):**
- 502/502 patient IDs match
- 502/502 both tumor samples
- 502/502 correct analyte types (RNA vs DNA)
- 3 sample-type mismatches (Primary ↔ Metastatic) excluded by the type constraint
- 3 low-A3/high-SBS2 outliers verified as correctly matched

**Usage:**
```bash
conda run -n NETWORK python scripts/TCGA/Step05_Revised_HNSC_A3_vs_SBS2.py
```

---

### VCF Pipeline (scripts/TCGA/VCF/)

Independent pipeline that downloads all TCGA somatic mutation data and recomputes COSMIC SBS signature weights. This replaced a pre-processed signature file of unknown provenance.

**Scale:** 10,939 samples across 33 cancer types, 87.6 million total mutations, 6.58 million PASS-filtered somatic SNVs.

**Key output:** `TCGA_SBS_signature_counts.tsv` — absolute mutation count attributions per COSMIC v3.4 signature per sample. Used as input to the revised Step 5.

**Comparison with original file (Diagnostic_Compare_SBS_Weight_Sources.py):**
- Original: 8,465 samples, 65 signatures (older COSMIC version), integer counts
- New: 10,936 samples, 86 signatures (COSMIC v3.4), integer counts
- HNSC-specific Spearman rho = 0.90 (strong rank agreement despite different decompositions)
- Decision: Use new counts for documented provenance and current COSMIC version

For full VCF pipeline documentation, see the `TCGA_Analysis` reference document.

---

### Original R Scripts (Steps 5-8, SUPERSEDED)

These R scripts produced the original 6-panel Figure 1. They are retained in the repository for reference but are superseded by `Step05_Revised_HNSC_A3_vs_SBS2.py` for the current manuscript.

| Script | Original Panel | Status |
|--------|---------------|--------|
| `Patient_Level_HNSCC_TCGA_A3s_vs_SBS2.R` | Fig 1a (scatter) | Superseded by Step05 Panel 1a |
| `Diagnostic_A3C_A3H_bystander.R` | Fig 1b (ROC) | Dropped per PI directive (4/2) |
| `Patient_Level_HNSCC_TCGA_3D_A3s.R` | Fig 1c (3D) | Dropped per PI directive (4/2) |
| `A3A_A3B_additive_SBS2.R` | Fig 1d-f (additive) | Superseded by Step05 Panels 1b-c |

---

### Troubleshooting Scripts

| Script | Purpose | Key Finding |
|--------|---------|-------------|
| `Diagnostic_Check_SBS2_Weight_Normalization.py` | Check if original SBS2 weights are ratio-based | CV=4.22, absolute counts confirmed |
| `Diagnostic_Compare_SBS_Weight_Sources.py` | Compare original vs new SigProfiler outputs | HNSC rho=0.90, different COSMIC versions |
| `Diagnostic_Verify_Barcode_Matching.py` | Audit all RNA-seq↔WES barcode pairs | 502/502 pass, 3 outliers verified |
| `Diagnostic_Figure1_Text_Numbers.py` | Extract text-ready numbers for results section | All quadrant medians, p-values, region counts |

---

### Figure 1 Summary (Revised)

**Title:** A3A and A3B expression is necessary but not sufficient for SBS2 mutagenesis in HNSCC.

| Panel | Content | Script |
|-------|---------|--------|
| a | HNSC summed A3A+A3B vs SBS2 scatter with colored regions | `Step05_Revised_HNSC_A3_vs_SBS2.py` |
| b | A3A vs A3B per tumor, colored by SBS2 (depth-sorted) | `Step05_Revised_HNSC_A3_vs_SBS2.py` |
| c | Box-and-whisker + 2×2 heatmap of median-split quadrants | `Step05_Revised_HNSC_A3_vs_SBS2.py` |
| Supp | Zoomed view of low-A3/high-SBS2 boundary (4 tumors labeled) | `Step05_Revised_HNSC_A3_vs_SBS2.py` |

**Narrative arc:**
1. A3A+A3B expression is necessary but not sufficient for SBS2 (Panel a)
2. Both A3A and A3B contribute, with SBS2 signal concentrated where both are high (Panel b)
3. A3B provides a baseline that A3A amplifies, with the highest burden requiring both (Panel c)
4. Unknown cofactors must regulate A3 activity → **motivates Question 2 (network analysis)**
