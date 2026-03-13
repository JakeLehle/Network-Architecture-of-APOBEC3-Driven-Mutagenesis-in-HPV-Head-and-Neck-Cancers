## Question 1: What is the fundamental relationship between A3 expression and observed A3-induced mutations in HNSCC?

### Rationale

APOBEC3 (A3) cytidine deaminases with access to the nuclear compartment (A3A, A3B, A3C, A3H) are capable of editing genomic DNA, producing characteristic C>T and C>G mutations captured by the COSMIC SBS2 mutational signature. However, the presence of A3 enzyme expression alone does not guarantee that A3-driven mutagenesis will occur in a given tumor. This question establishes the foundational observation of the paper: **A3 expression is necessary but not sufficient for SBS2 mutagenesis in HNSCC**, implying the existence of unknown cofactors that regulate A3 enzymatic activity on genomic DNA.

### Data Sources

| Source | Description |
|--------|-------------|
| [GDC Data Portal](https://portal.gdc.cancer.gov/) | TCGA RNA-seq gene expression (STAR - Counts workflow) across 33 cancer types |
| [COSMIC SBS Signatures](https://cancer.sanger.ac.uk/signatures/) | SBS mutational signature weights derived from TCGA whole-exome sequencing, provided as `Mutation_Table_Tumors_TCGA.tsv` |

### Directory Structure

```
scripts/TCGA/
├── DOWNLOAD_TCGA_COUNT_DATA.sh          # SLURM Job 1: Data acquisition
├── Download_RNA-seq_Counts_TCGA.R       #   Step 1 — Download from GDC
├── Organize_RNA-seq_Counts_TCGA.R       #   Step 2 — Organize & build manifest
├── Master_TCGA_RNA-seq_Counts_Table.R   #   Step 3 — Build master expression tables
│
├── ANALYZE_TCGA_COUNT_DATA.sh           # SLURM Job 2: Analysis & visualization
├── Prep_mutation_analysis_files.R       #   Step 4 — Match A3 expression ↔ SBS signatures
├── Patient_Level_HNSCC_TCGA_A3s_vs_SBS2.R   #   Step 5 — HNSC patient-level scatter
├── Patient_Level_HNSCC_TCGA_3D_A3s.R        #   Step 6 — HNSC 3D individual A3s
│
└── TROUBLESHOOTING/                     # Supplementary/diagnostic scripts
    ├── Compare_A3_Expression_Sources.R
    ├── FILTER_MASTER_TCGA_TSV.R
    ├── Global_TCGA_A3s_vs_SBS2.R
    ├── Patient_Level_HNSCC_TCGA_Network_Analysis.R
    └── TCGA_troubleshooting.R
```

### Pipeline Overview

The analysis is executed as two sequential SLURM batch jobs on an HPC cluster:

```
 SLURM Job 1: DOWNLOAD_TCGA_COUNT_DATA.sh
 ─────────────────────────────────────────
   Step 1: Download from GDC
       │
       ▼
   Step 2: Organize files & build manifest
       │
       ▼
   Step 3: Build pan-cancer master expression tables
                    │
                    ▼
          TCGA_master_FPKM_UQ.tsv


 SLURM Job 2: ANALYZE_TCGA_COUNT_DATA.sh
 ─────────────────────────────────────────
          TCGA_master_FPKM_UQ.tsv
          Mutation_Table_Tumors_TCGA.tsv  ← external input
                    │
                    ▼
   Step 4: Match A3 expression ↔ SBS mutation signatures
                    │
               ┌────┴────┐
               ▼         ▼
   Step 5: Patient-   Step 6: 3D
   level A3 vs SBS2   individual
   scatter            A3s vs SBS2
               │         │
               └────┬────┘
                    ▼
                FIGURE 1
```

---

### SLURM Job 1: Data Acquisition

**Batch script:** `DOWNLOAD_TCGA_COUNT_DATA.sh`

Submits Steps 1–3 as a serial pipeline. Downloads all TCGA RNA-seq data, organizes the files, and builds the pan-cancer master expression tables. This is a long-running job (up to 7 days walltime) due to GDC download speeds.

```bash
Rscript Download_RNA-seq_Counts_TCGA.R
Rscript Organize_RNA-seq_Counts_TCGA.R
Rscript Master_TCGA_RNA-seq_Counts_Table.R
```

**Environment:** `RNA-seq_NovoGene` conda environment, `ulimit -s unlimited`

---

#### Step 1: Download TCGA RNA-seq Data

**Script:** `Download_RNA-seq_Counts_TCGA.R`

Downloads STAR-Counts RNA-seq gene expression quantification files from the GDC for all 33 TCGA cancer types using the `TCGAbiolinks` R package. The script includes retry logic for failed downloads and skips cancer types that have already been successfully downloaded.

**Dependencies:** `TCGAbiolinks`, `SummarizedExperiment`

**Input:**
- GDC API (internet connection required)

**Output:**
- `GDCdata/TCGA-{CANCER_TYPE}/Transcriptome_Profiling/Gene_Expression_Quantification/` — raw per-sample `.tsv` files organized by cancer type

**Key details:**
- Downloads all 33 TCGA projects (BRCA, LUAD, LUSC, ..., HNSC, ..., UCEC)
- Each sample file contains raw counts, TPM, FPKM, and FPKM-UQ quantifications
- Files are identified by GDC UUID

---

#### Step 2: Organize Files and Build Master Manifest

**Script:** `Organize_RNA-seq_Counts_TCGA.R`

Queries GDC metadata to build a master manifest that maps each file's UUID to the full 28-character TCGA aliquot barcode (Entity_ID), Case_ID, cancer type, and tissue status (Tumor vs. Normal). Reorganizes the raw downloaded files into a structured directory.

**Dependencies:** `TCGAbiolinks`, `dplyr`, `tidyr`

**Input:**
- `GDCdata/` directory from Step 1

**Output:**
- `TCGA_RNAseq_master_manifest.tsv` — manifest mapping UUID → Entity_ID → Case_ID → cancer type → tissue status
- `organized_counts/TCGA-{CANCER_TYPE}/{Tumor|Normal}/` — reorganized expression files with informative filenames (`{Case_ID}_{UUID}.rna_seq.augmented_star_gene_counts.tsv`)

**Key details:**
- Entity_ID is the full 28-character TCGA aliquot barcode (e.g., `TCGA-BH-A18H-01A-11R-A12D-07`), which serves as the unique sample identifier throughout the pipeline
- The `cases` column from TCGAbiolinks metadata is used as the source for Entity_IDs
- Tissue status is determined from the sample type code embedded in the barcode (positions 14–15): `01`–`09` = Tumor, `10`–`19` = Normal

---

#### Step 3: Build Pan-Cancer Master Expression Tables

**Script:** `Master_TCGA_RNA-seq_Counts_Table.R`

Reads all organized per-sample expression files and assembles three pan-cancer master expression matrices (TPM, FPKM, FPKM-UQ), each containing all genes across all samples from all 33 cancer types.

**Dependencies:** `data.table`

**Input:**
- `organized_counts/` directory from Step 2
- `TCGA_RNAseq_master_manifest.tsv` from Step 2

**Output:**
- `TCGA_master_TPM.tsv`
- `TCGA_master_FPKM.tsv`
- `TCGA_master_FPKM_UQ.tsv` ← **primary expression matrix used in all downstream analyses**

**File structure (all three files share this format):**
```
Row 1:  Column headers — Project_ID, Tissue_Type, Case_ID, File_ID, Entity_ID, ENSG00000000003.15, ...
Row 2:  Gene symbols   — NA, NA, NA, NA, NA, TSPAN6, ...
Row 3:  Gene biotypes  — NA, NA, NA, NA, NA, protein_coding, ...
Row 4+: Sample data    — TCGA-BRCA, Tumor, TCGA-BH-A18H, UUID, Entity_ID, 1234.56, ...
```

**Key details:**
- FPKM-UQ (upper quartile normalized FPKM) is used as the primary normalization for all downstream analyses because it is robust to highly expressed outlier genes
- Metadata columns (1–5) are embedded alongside the expression matrix for self-contained sample tracking
- Gene identifiers are versioned ENSG IDs (e.g., `ENSG00000000003.15`)

---

### SLURM Job 2: Analysis & Visualization

**Batch script:** `ANALYZE_TCGA_COUNT_DATA.sh`

Submits Steps 4–6 as a serial pipeline. Matches A3 expression to mutation signatures and generates the Figure 1 panels. Requires 48 GB memory.

```bash
Rscript Prep_mutation_analysis_files.R
Rscript Patient_Level_HNSCC_TCGA_A3s_vs_SBS2.R
Rscript Patient_Level_HNSCC_TCGA_3D_A3s.R
```

**Environment:** `RNA-seq_NovoGene` conda environment, `ulimit -s unlimited`, `--mem=48G`

---

#### Step 4: Match A3 Expression with SBS Mutation Signatures

**Script:** `Prep_mutation_analysis_files.R`

Extracts expression values for the four nuclear A3 family members (A3A, A3B, A3C, A3H) from the pan-cancer FPKM-UQ master table and matches them to COSMIC SBS mutational signature weights by Entity_ID. Produces two paired output files containing only samples present in both datasets.

**Dependencies:** `data.table`, `dplyr`

**Input:**
- `TCGA_master_FPKM_UQ.tsv` from Step 3
- `Mutation_Table_Tumors_TCGA.tsv` — COSMIC SBS signature weights per TCGA sample (external input)

**Output:**
- `A3_Expression_FPKM_UQ_Matched.tsv` — columns: `Project_ID`, `Tissue_Type`, `Entity_ID`, `APOBEC3A`, `APOBEC3B`, `APOBEC3C`, `APOBEC3H`
- `Mutation_Signatures_Matched.tsv` — all COSMIC SBS signature columns, matched to the same sample set

**Key details:**
- A3 genes are identified by gene symbol from Row 2 of the master file
- Matching is performed on Entity_ID (expression file) ↔ `TCGA_Gene_Expression_Entity_ID` (mutation file)
- Only samples present in **both** files are retained, ensuring a clean one-to-one mapping
- These two matched files serve as the shared input for Steps 5 and 6, as well as for subsequent questions in the paper (Questions 2+)

---

#### Step 5: HNSC Patient-Level Nuclear A3s vs. SBS2

**Script:** `Patient_Level_HNSCC_TCGA_A3s_vs_SBS2.R`

Generates a patient-level scatter plot for HNSC tumors only, showing the summed nuclear A3 expression (A3A + A3B + A3C + A3H) vs. SBS2 weight for each individual tumor. This is the key plot that reveals the three populations:

1. **High A3 / High SBS2** — tumors with both A3 expression and active A3-driven mutagenesis
2. **High A3 / Low SBS2** — tumors with A3 expression but no evidence of A3-driven mutagenesis
3. **No Low-A3 / High-SBS2 tumors** — SBS2 mutations do not occur without A3 expression

**Dependencies:** `data.table`, `dplyr`, `ggplot2`

**Input:**
- `A3_Expression_FPKM_UQ_Matched.tsv` from Step 4
- `Mutation_Signatures_Matched.tsv` from Step 4

**Output:**
- HNSC patient-level scatter plot (summed Nuclear A3s vs. SBS2 per tumor)

**Key details:**
- Filters for `Project_ID == "TCGA-HNSC"` and `Tissue_Type == "Tumor"` only
- Nuclear A3 sum = A3A + A3B + A3C + A3H
- The absence of tumors in the low-A3/high-SBS2 quadrant is the central observation: **A3 expression is necessary but not sufficient for SBS2 mutagenesis**

---

#### Step 6: HNSC 3D Individual A3 Expression vs. SBS2

**Script:** `Patient_Level_HNSCC_TCGA_3D_A3s.R`

Generates a 3D scatter plot of individual A3 family member expression (A3A × A3B × A3C) for each HNSC tumor, with points colored by SBS2 mutational signature weight. This decomposes the summed A3 signal from Step 5 to show the relative contribution of each family member.

**Dependencies:** `data.table`, `dplyr`, `plotly`, `htmlwidgets`, `scatterplot3d`, `RColorBrewer`

**Input:**
- `A3_Expression_FPKM_UQ_Matched.tsv` from Step 4
- `Mutation_Signatures_Matched.tsv` from Step 4

**Output:**
- `Figure_HNSC_3D_Individual_A3_vs_SBS2.html` — interactive 3D plot (plotly)
- `Figure_HNSC_3D_Individual_A3_vs_SBS2.pdf` — static 3D plot
- `Figure_HNSC_3D_Individual_A3_vs_SBS2.png` — static 3D plot
- `HNSC_Patient_Level_Individual_A3_SBS2_Data.tsv` — patient-level data table
- Per-A3 Spearman correlation results with SBS2

**Key details:**
- X = A3A, Y = A3B, Z = A3C (all FPKM-UQ)
- Color = SBS2 weight (viridis colormap, low → high)
- A3H is excluded from the 3D visualization due to consistently low expression in HNSC
- Calculates Spearman correlations for each individual A3 vs. SBS2

---

### Figure 1 Summary

Together, Steps 5 and 6 produce the panels of **Figure 1**:

| Panel | Content | Script |
|-------|---------|--------|
| A | HNSC patient-level summed nuclear A3s vs. SBS2 scatter | Step 5: `Patient_Level_HNSCC_TCGA_A3s_vs_SBS2.R` |
| B | HNSC 3D individual A3A/A3B/A3C colored by SBS2 | Step 6: `Patient_Level_HNSCC_TCGA_3D_A3s.R` |

**Conclusion from Figure 1:** A3 expression is required but not sufficient to drive SBS2 mutagenesis in HNSCC, suggesting that additional cofactors modulate A3 enzymatic activity on genomic DNA. This motivates **Question 2**: identifying those cofactors through differential network analysis.

---

### Troubleshooting Scripts

The following scripts in `scripts/TCGA/TROUBLESHOOTING/` are not part of the main pipeline but were used during development for diagnostics and validation:

| Script | Purpose |
|--------|---------|
| `TCGA_troubleshooting.R` | Verifies TCGAbiolinks metadata column structure and confirms Entity_ID barcode format |
| `Compare_A3_Expression_Sources.R` | Compares A3 expression values across TPM, FPKM, and FPKM-UQ normalizations to confirm data provenance |
| `FILTER_MASTER_TCGA_TSV.R` | Filters pan-cancer master table to HNSC only; superseded by direct filtering in Step 5 |
| `Global_TCGA_A3s_vs_SBS2.R` | Pan-cancer scatter of average A3 expression vs. average SBS2 per cancer type; exploratory analysis |
| `Patient_Level_HNSCC_TCGA_Network_Analysis.R` | Early prototype of differential network analysis in R; superseded by the dedicated network pipeline (Question 2) |
