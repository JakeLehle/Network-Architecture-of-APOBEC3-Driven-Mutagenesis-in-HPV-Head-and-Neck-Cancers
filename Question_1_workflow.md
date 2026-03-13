## Question 1: What is the fundamental relationship between A3 expression and observed A3-induced mutations in HNSCC?

### Rationale

APOBEC3 (A3) cytidine deaminases with access to the nuclear compartment (A3A, A3B, A3C, A3H) are capable of editing genomic DNA, producing characteristic C>T and C>G mutations captured by the COSMIC SBS2 mutational signature. However, the presence of A3 enzyme expression alone does not guarantee that A3-driven mutagenesis will occur in a given tumor. This question establishes the foundational observation of the paper: **A3 expression is necessary but not sufficient for SBS2 mutagenesis in HNSCC**, implying the existence of unknown cofactors that regulate A3 enzymatic activity on genomic DNA. Additionally, this question dissects the individual contributions of each A3 family member, establishing that A3B provides a baseline mutagenic state that A3A amplifies, while A3C and A3H are bystanders reflecting immune cell infiltration rather than independent mutagenic activity.

### Data Sources

| Source | Description |
|--------|-------------|
| [GDC Data Portal](https://portal.gdc.cancer.gov/) | TCGA RNA-seq gene expression (STAR - Counts workflow) across 33 cancer types |
| [COSMIC SBS Signatures](https://cancer.sanger.ac.uk/signatures/) | SBS mutational signature weights derived from TCGA whole-exome sequencing, provided as `Mutation_Table_Tumors_TCGA.tsv` |

### Directory Structure

```
scripts/TCGA/
├── DOWNLOAD_TCGA_COUNT_DATA.sh              # SLURM Job 1: Data acquisition
├── Download_RNA-seq_Counts_TCGA.R           #   Step 1 — Download from GDC
├── Organize_RNA-seq_Counts_TCGA.R           #   Step 2 — Organize & build manifest
├── Master_TCGA_RNA-seq_Counts_Table.R       #   Step 3 — Build master expression tables
│
├── ANALYZE_TCGA_COUNT_DATA.sh               # SLURM Job 2: Analysis & visualization
├── Prep_mutation_analysis_files.R           #   Step 4 — Match A3 expression ↔ SBS signatures
├── Patient_Level_HNSCC_TCGA_A3s_vs_SBS2.R  #   Step 5 — HNSC patient-level scatter (Fig 1a)
├── Diagnostic_A3C_A3H_bystander.R           #   Step 6 — ROC/AUC & bystander analysis (Fig 1b)
├── Patient_Level_HNSCC_TCGA_3D_A3s.R       #   Step 7 — 3D A3A×A3B×A3C scatter (Fig 1c)
├── A3A_A3B_additive_SBS2.R                  #   Step 8 — A3B saturation & additive model (Fig 1d,e,f)
│
└── TROUBLESHOOTING/                         # Supplementary/diagnostic scripts
    ├── A3_Contribution_to_SBS2.R
    ├── Compare_A3_Expression_Sources.R
    ├── Diagnose_nonlinear_A3_Relationships.R
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
                    ▼
          A3_Expression_FPKM_UQ_Matched.tsv
          Mutation_Signatures_Matched.tsv
                    │
       ┌────────────┼────────────┬────────────┐
       ▼            ▼            ▼            ▼
   Step 5:      Step 6:      Step 7:      Step 8:
   A3 vs SBS2   ROC/AUC      3D A3s       A3B/A3A
   scatter       bystander    scatter      additive
   (Fig 1a)     (Fig 1b)     (Fig 1c)     (Fig 1d,e,f)
       │            │            │            │
       └────────────┴────────────┴────────────┘
                         │
                      FIGURE 1
```

### Environment

All scripts are executed using the `RNA-seq_NovoGene` conda environment (see `environments/RNA-seq_NovoGene.yml`).

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

Submits Steps 4–8 as a serial pipeline. Matches A3 expression to mutation signatures and generates all Figure 1 panels. Requires 48 GB memory.

```bash
Rscript Prep_mutation_analysis_files.R
Rscript Patient_Level_HNSCC_TCGA_A3s_vs_SBS2.R
Rscript Diagnostic_A3C_A3H_bystander.R
Rscript Patient_Level_HNSCC_TCGA_3D_A3s.R
Rscript A3A_A3B_additive_SBS2.R
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
- These two matched files serve as the shared input for all subsequent steps (5–8), as well as for the network analysis in Question 2

---

#### Step 5: HNSC Patient-Level Nuclear A3s vs. SBS2 → Figure 1a

**Script:** `Patient_Level_HNSCC_TCGA_A3s_vs_SBS2.R`

Generates a patient-level scatter plot for HNSC tumors showing summed nuclear A3 expression (A3A + A3B + A3C + A3H) vs. SBS2 weight, with three data-driven colored background regions highlighting the key biological populations.

**Dependencies:** `data.table`, `dplyr`, `ggplot2`, `ggbreak`

**Input:**
- `A3_Expression_FPKM_UQ_Matched.tsv` from Step 4
- `Mutation_Signatures_Matched.tsv` from Step 4

**Output:**
- `Figure_HNSC_Patient_Level_A3_vs_SBS2.pdf` / `.png`
- `HNSC_Patient_Level_A3_SBS2_Data.tsv` — patient-level data with region assignments

**Key details:**
- Filters for TCGA-HNSC primary tumors only
- Three colored regions delineated by a data-driven diagonal (slope computed from the steepest SBS2/A3 ratio in the left zone with 5% buffer) and a horizontal line at the median SBS2:
  - **Coral (#ed6a5a):** A3 present + SBS2 above median — tumors with active A3-driven mutagenesis
  - **Cream (#f4f1bb):** A3 present + SBS2 below median — tumors where A3s are expressed but mutagenesis is suppressed
  - **Teal (#9bc1bc):** No A3 expression — contains zero tumors with elevated SBS2 (A3 is necessary)
- Points colored to match their region; legend above the plot with per-region sample counts
- X-axis broken at 350–600 FPKM-UQ to accommodate outlier without compressing the main data cloud
- **Conclusion:** A3 expression is necessary but not sufficient for SBS2 mutagenesis

---

#### Step 6: A3C/A3H Bystander Analysis — ROC/AUC → Figure 1b

**Script:** `Diagnostic_A3C_A3H_bystander.R`

Evaluates the independent predictive capacity of each A3 family member for SBS2 status using three complementary approaches: (A) fraction of high-SBS2 tumors with elevated expression, (B) ROC/AUC classification, and (C) expression shift between SBS2 groups. Also includes a co-expression heatmap. The ROC curve panel is used as Figure 1b.

**Dependencies:** `data.table`, `dplyr`, `ggplot2`, `patchwork`, `viridis`, `pROC`

**Input:**
- `A3_Expression_FPKM_UQ_Matched.tsv` from Step 4
- `Mutation_Signatures_Matched.tsv` from Step 4

**Output:**
- `Diag_Bridge_A3C_A3H_Bystander_Main.pdf` / `.png` — composite of all three approaches
- `Diag_Bridge_A3C_A3H_Bystander_Supp.pdf` / `.png` — supplementary panels (AUC bar, effect size, co-expression)
- Individual panel PDFs (Approaches A, B, C, bonus heatmap)
- `Diag_ApproachA_Scene_Presence.tsv`, `Diag_ApproachB_AUC_Results.tsv`, `Diag_ApproachC_Expression_Shift.tsv`

**Key details:**
- Tumors split at median SBS2 into SBS2-High and SBS2-Low groups
- ROC curves (Figure 1b) show A3A (AUC = 0.599) and A3B (AUC = 0.579) modestly above chance, A3C (AUC = 0.578) marginal, and A3H (AUC = 0.519) at chance level
- DeLong tests compare each A3's AUC to A3B
- Co-expression heatmap reveals A3C/A3H clustering (shared immune cell source) distinct from A3A/A3B
- **Conclusion:** A3H can be excluded; A3A, A3B, and A3C require further investigation via 3D visualization. The overall low AUC values for all family members reinforce that non-A3 cofactors drive SBS2 status.

---

#### Step 7: HNSC 3D Individual A3 Expression vs. SBS2 → Figure 1c

**Script:** `Patient_Level_HNSCC_TCGA_3D_A3s.R`

Generates a 3D scatter plot of individual A3 family member expression (A3A × A3B × A3C) for each HNSC tumor, with points colored by SBS2 weight. A3H is excluded based on Step 6. A3C is retained because its AUC could not exclude it, though its expression is expected to reflect immune infiltration rather than intrinsic mutagenic activity.

**Dependencies:** `data.table`, `dplyr`, `plotly`, `htmlwidgets`, `viridis`, `scatterplot3d`

**Input:**
- `A3_Expression_FPKM_UQ_Matched.tsv` from Step 4
- `Mutation_Signatures_Matched.tsv` from Step 4

**Output:**
- `Figure_HNSC_3D_Individual_A3_vs_SBS2.html` — interactive 3D plot (plotly)
- `Figure_HNSC_3D_Individual_A3_vs_SBS2_static.pdf` / `.png` — static 3D plots
- `HNSC_Patient_Level_Individual_A3_SBS2_Data.tsv` — patient-level data table
- `HNSC_Individual_A3_SBS2_Correlations.tsv` — per-A3 Spearman correlations with SBS2

**Key details:**
- X = A3A, Y = A3B, Z = A3C (all FPKM-UQ); Color = SBS2 weight (viridis)
- The plot is rotated to reveal that along the A3B axis, tumors show a consistent moderate SBS2 burden (blue to dark green), while tumors with concurrent high A3A expression exhibit extreme SBS2 accumulation (bright yellow)
- **Conclusion:** A3B establishes a baseline mutagenic state; A3A amplifies it. A3C does not independently influence SBS2 and likely reflects tumor-infiltrating lymphocytes. This motivates the focused A3B/A3A analysis in Step 8.

---

#### Step 8: A3B Saturation and A3A/A3B Additive Model → Figure 1d, 1e, 1f

**Script:** `A3A_A3B_additive_SBS2.R`

Formally characterizes the non-linear A3B–SBS2 relationship and quantifies the additive contributions of A3A and A3B through a threshold sweep and 2×2 stratification model. Produces the final three panels of Figure 1.

**Dependencies:** `data.table`, `dplyr`, `ggplot2`, `patchwork`, `viridis`

**Input:**
- `A3_Expression_FPKM_UQ_Matched.tsv` from Step 4
- `Mutation_Signatures_Matched.tsv` from Step 4

**Output:**
- `Figure_A3B_A3A_Additive_SBS2.pdf` / `.png` — composite of all three panels
- `Panel_Threshold_Sweep.pdf` — Figure 1d
- `Panel_Ordered_Boxplot.pdf` — Figure 1e
- `Panel_2x2_Heatmap.pdf` — Figure 1f
- `A3B_Threshold_Sweep_Correlations.tsv`, `A3B_A3A_Quadrant_Stats.tsv`

**Key details:**

- **Panel 1d — Threshold sweep:** Spearman correlation between A3B and SBS2 computed separately for tumors below and above successive A3B quantile thresholds (20th–80th percentile in 5% steps). Both lines hold steady at rho = 0.1–0.2 until the 80th percentile, where the above-threshold correlation collapses to negative values, demonstrating A3B saturation.

- **Panel 1e — Ordered boxplot:** Tumors stratified into four groups by median A3B and A3A expression, ordered left-to-right by ascending median SBS2. Diamond markers connected by a trend line trace the additive progression: A3B-/A3A- (median SBS2 ≈ 4.5) → A3B+/A3A- (≈ 11) → A3B-/A3A+ (≈ 15) → A3B+/A3A+ (≈ 17). Pairwise Wilcoxon tests with BH correction.

- **Panel 1f — 2×2 heatmap:** Median SBS2 weight per A3B×A3A quadrant with viridis color scale and sample sizes.

- **Conclusion:** A3B drives a baseline accumulation of SBS2 mutations that saturates at high expression. A3A amplifies SBS2 beyond this ceiling. The highest mutagenic burden requires both enzymes — alongside as-yet-unidentified cofactors that license A3 access to genomic DNA.

---

### Figure 1 Summary

**Title:** Nuclear APOBEC3 expression defines the mutational landscape of HNSCC but does not fully explain SBS2 accumulation.

| Panel | Content | Script | Step |
|-------|---------|--------|------|
| a | HNSC patient-level summed A3s vs. SBS2 scatter with colored regions | `Patient_Level_HNSCC_TCGA_A3s_vs_SBS2.R` | 5 |
| b | ROC curves: per-A3 predictive capacity for SBS2 status | `Diagnostic_A3C_A3H_bystander.R` | 6 |
| c | 3D scatter: A3A × A3B × A3C colored by SBS2 | `Patient_Level_HNSCC_TCGA_3D_A3s.R` | 7 |
| d | A3B threshold sweep: saturation of A3B–SBS2 correlation | `A3A_A3B_additive_SBS2.R` | 8 |
| e | Ordered boxplot: SBS2 across A3B×A3A quadrants | `A3A_A3B_additive_SBS2.R` | 8 |
| f | 2×2 heatmap: median SBS2 per A3B×A3A quadrant | `A3A_A3B_additive_SBS2.R` | 8 |

**Narrative arc:**
1. A3 expression is necessary but not sufficient for SBS2 (panel a)
2. No single A3 strongly predicts SBS2; A3H excluded, A3C marginal (panel b)
3. A3B provides steady-state SBS2; A3A drives extreme values; A3C is a bystander (panel c)
4. A3B's contribution saturates at high expression (panel d)
5. A3A and A3B work additively — the highest SBS2 requires both (panels e, f)
6. Unknown cofactors must regulate A3 activity → **motivates Question 2 (network analysis)**

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
| `A3_Contribution_to_SBS2.R` | Initial exploration of individual A3 family member contributions to SBS2 (partial correlations, regression, conditional analysis); findings refined into Steps 6 and 8 |
| `Diagnose_nonlinear_A3_Relationships.R` | Diagnostic testing four approaches (threshold, log-transform, sequential ANOVA, gate/driver model) to characterize non-linear A3B–SBS2 relationship; best approaches incorporated into Step 8 |
