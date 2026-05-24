# Question 6 Workflow

## Question 6: How does HPV16 viral lifecycle stage determine the divergent mutagenic programs observed in APOBEC3-active epithelial cells?

### Rationale

Figures 3-4 established that epithelial cells bifurcate into two populations with distinct genomic damage profiles: SBS2-HIGH cells accumulate A3A-driven point mutations while CNV-HIGH cells accumulate chromosomal instability with A3B-dominant expression. The Figure 4 network analysis identified a co-expression chain linking A3A to immune and inflammatory genes (CCL20, LCN2, SMOX via the Harris interactors RALY and HNRNPA2B1) that was completely disrupted in CNV-HIGH cells (100% negative DIFF edges on A3A). Question 5 traced somatic variants in high-contributor patients to HPV infection and antigen presentation pathways, motivating a direct investigation of how HPV16 infection context shapes these divergent phenotypes.

This question addresses a gap in the HPV field: while HPV lifecycle stages are well characterized in stratified epithelium (Doorbar et al., 2012; Moody and Laimins, 2010; McBride and Warburton, 2017), the relationship between viral lifecycle phase and APOBEC mutational outcomes has not been resolved at single-cell resolution. The central finding is that **SBS2-HIGH cells carry HPV16 in E2-regulated maintenance while CNV-HIGH cells carry HPV16 in ATM-dependent productive infection**, and this lifecycle distinction organizes the entire host response including which APOBEC enzyme is active (A3A vs A3B), which class of genomic damage accumulates (SBS2 vs CNV), and whether the cell is immunologically visible or engaged in viral replication. HPV gene expression was normalized as fractions of total viral reads per cell to control for the approximately two-fold higher viral load in CNV-HIGH cells (reflecting episomal amplification during productive infection). A panel of 53 host marker genes drawn from the Doorbar, Laimins, and McBride laboratories validated predictions from established HPV biology, including ATM/CHK2-dependent G2/M arrest for productive replication (Moody and Laimins, 2009), elevated antigen presentation in maintenance-phase cells, and functional p53 signaling consistent with intact E2 suppression of E6/E7.

**Note:** The neoantigen analysis (formerly Phase 5B) has been refactored into a standalone six-step pipeline in `scripts/NEOANTIGEN/`. See `scripts/NEOANTIGEN/NEOANTIGEN_WALKTHROUGH.md` for that pipeline's documentation.

---

## Directory Structure

```
scripts/HPV_ANALYSIS/
├── Phase1_HPV16_Data_Inventory.py           # Phase 1: per-cell HPV16 data consolidation
├── Phase1_HPV16_Data_Inventory_v2.py        # Phase 1 (v2): updated inventory
├── Phase1_5_Index_Diagnostic.py             # Phase 1.5: index alignment diagnostic
├── Phase2_Raw_HPV16_Counts.py               # Phase 2: raw UMI count recovery from Kraken2
├── Phase3_HPV16_Populations_and_Genome.py   # Phase 3: HPV16 threshold, population definition, genome probing
├── Phase4_HPV16_Genome_Alignment.py         # Phase 4: minimap2 alignment to HPV16, per-cell gene counts
├── Phase5A_Population_Consolidation.py      # Phase 5A: multi-axis population profiling (v1, superseded)
├── Phase5A_Revised_Population_Discovery.py  # Phase 5A: data-driven population selection (v2, superseded)
├── Phase5A_v2_DEG_Analysis.py               # Phase 5A: DE and GSEA between populations (superseded)
├── Generate_Figure6_Lifecycle_Panels.py     # Figure 6 panel generation (current)
├── Generate_Figure6_Panels.py               # Figure 6 panel generation (pre-lifecycle reframe, archived)
├── parse_chimeric_junctions.py              # Chimeric junction parsing utility
├── sample_list.txt                          # Sample manifest (31 GSE173468 samples)
├── RUN_HPV_ANALYSIS.sh                      # SLURM master runner (Phases 1-4)
│
└── TROUBLESHOOTING/
    └── Diagnostic_HPV_Lifecycle_Markers.py  # Lifecycle marker diagnostic (3 analyses)
```

---

## Input Data

| Input | Path | Description |
|-------|------|-------------|
| AnnData (final) | `data/FIG_4/00_input/adata_final.h5ad` | 155,650 cells, 27,736 genes (ClusterCatcher output) |
| AnnData (viral) | `data/FIG_6/00_input/adata_v_pp.h5ad` | Viral detection object with Kraken2 annotations |
| Signature weights | `data/FIG_4/00_input/signature_weights_per_cell.txt` | Per-cell SBS signature weights (SigProfiler) |
| Three-group assignments | `data/FIG_4/01_group_selection/three_group_assignments.tsv` | 546 SBS2-HIGH, 546 CNV-HIGH, 546 NORMAL |
| Kraken2 matrices | `SC/fastq/.../kraken2_filtered_feature_bc_matrix/` | Per-sample sparse matrices with raw viral UMI counts |
| Unmapped BAMs | `SC/fastq/.../possorted_genome_bam_unmapped.bam` | Unmapped reads with cell barcode (CB) tags |
| Unmapped FASTQs | `SC/fastq/.../possorted_genome_bam_unmapped.fq` | Unmapped reads in FASTQ format for alignment |
| HPV16 reference | `data/FIG_6/03_hpv16_genome/HPV16_NC_001526.4.fa` | NC_001526.4, 7,906 bp (downloaded by Phase 4) |

---

## Pipeline Overview

The pipeline has five phases. Phases 1-4 extract and process per-cell viral data (population-independent). Phase 5A defined the original two-population model and is now superseded by the Figure 4 three-group selection. The lifecycle diagnostic and figure generation scripts use Phase 4 output combined with the Figure 4 population definitions.

```
Phase 1 ─► Phase 1.5 ─► Phase 2 ─► Phase 3 ─► Phase 4 ──────────────┐
(inventory)  (index dx)   (raw UMI)  (threshold)  (genome alignment)   │
                                                                        ▼
Figure 4 populations ──► Diagnostic_HPV_Lifecycle_Markers.py ──► Generate_Figure6_Lifecycle_Panels.py
(three_group_assignments.tsv)     (marker profiling)                    (figure panels)
```

### Execution

Phases 1-4 run sequentially via SLURM:

```bash
cd scripts/HPV_ANALYSIS/
sbatch RUN_HPV_ANALYSIS.sh
```

The lifecycle diagnostic and figure generation run interactively:

```bash
conda run -n NETWORK python TROUBLESHOOTING/Diagnostic_HPV_Lifecycle_Markers.py
conda run -n NETWORK python Generate_Figure6_Lifecycle_Panels.py
```

---

## Script Descriptions

### Phase 1: Data Inventory

**`Phase1_HPV16_Data_Inventory.py`** consolidates all per-cell measurements into a master table for downstream analysis. Loads `adata_final.h5ad` and merges HPV16 status from the viral detection object, SBS2 weights from signature refitting, CNV scores, CytoTRACE2 stemness, and cell annotations. Outputs cross-tabulations of HPV16 status against SBS2 group, cancer cell status, and cnv_leiden clusters.

**`Phase1_HPV16_Data_Inventory_v2.py`** is the updated version with revised column mappings and additional diagnostic outputs.

| Output | Description |
|--------|-------------|
| `basal_cell_master_table.tsv` | One row per epithelial cell, all measurements |
| `diagnostic_report.txt` | Comprehensive profiling report |
| `HPV16_distribution_by_patient.tsv` | Per-patient HPV16 counts |

### Phase 1.5: Index Alignment Diagnostic

**`Phase1_5_Index_Diagnostic.py`** runs alignment diagnostics between the Kraken2 viral detection indices and the main single-cell object to verify barcode mapping integrity before proceeding to raw count extraction in Phase 2.

### Phase 2: Raw HPV16 Count Recovery

**`Phase2_Raw_HPV16_Counts.py`** recovers actual UMI counts from Kraken2 per-sample sparse matrices. The normalized counts in `adata_v_pp.h5ad` (CPM + log1p) saturate at 13.82, making HPV16 expression effectively binary. This script reads the raw 10x-format matrices to get true integer counts per cell.

For each of 31 GSE173468 samples, the script reads `matrix.mtx.gz`, `barcodes.tsv.gz`, and `features.tsv.gz` from the Kraken2 filtered feature barcode matrix directory. It identifies the HPV16 feature row (taxonomy ID 333760 or name "Human papillomavirus 16"), extracts per-cell counts, and maps barcodes to the integrated single-cell object.

| Output | Description |
|--------|-------------|
| `raw_HPV16_counts_all_cells.tsv` | Per-cell raw counts (HPV16+ cells only) |
| `basal_cell_master_table_with_raw_HPV16.tsv` | Master table with raw counts merged |
| `per_sample_hpv16_summary.tsv` | Per-sample HPV16 detection statistics |

### Phase 3: Population Definition and Genome Probing

**`Phase3_HPV16_Populations_and_Genome.py`** applies L-method piecewise linear regression to the sorted raw HPV16 count distribution to identify a data-driven threshold separating background noise from genuine viral signal. The optimal breakpoint was identified at 8 UMI reads per cell, yielding three classification tiers: HPV16-negative (0), HPV16-ambiguous (1-7), and HPV16-positive (8+).

This script also probes for Kraken2 intermediate files (unmapped BAMs, FASTQ files) needed for Phase 4 genome alignment and determines the alignment approach.

**Note:** The population definitions from this script (HPV-neg/low/high crossed with cnv_leiden) are superseded by the three-group selection in Figure 4 Step00B (`three_group_assignments.tsv`). The HPV16 threshold and genome probing outputs remain valid.

| Output | Description |
|--------|-------------|
| `population_assignments.tsv` | Per-cell HPV status and population labels |
| `elbow_detection_residuals.tsv` | L-method threshold statistics |
| `approach_determination.txt` | Genome alignment approach decision |

### Phase 4: HPV16 Genome Alignment

**`Phase4_HPV16_Genome_Alignment.py`** is the core alignment script. It downloads the HPV16 reference genome (NC_001526.4, 7,906 bp), builds a minimap2 index, aligns unmapped reads from all samples to the viral genome, extracts cell barcodes from the unmapped BAM files, and maps alignment positions to HPV16 gene annotations.

HPV16 gene regions follow standard annotations: E6 (83-559), E7 (562-858), E1 (865-2813), E2 (2756-3852), E4 (3332-3619), E5 (3849-4100), L2 (4236-5657), L1 (5560-7155), URR (7156-7906 and 1-82). The primary output is a per-cell gene count matrix (76,978 cells) that serves as input for all downstream lifecycle analysis.

| Output | Description |
|--------|-------------|
| `per_cell_hpv16_gene_counts.tsv` | Cells x HPV16 genes matrix (76,978 cells) |
| `per_sample_alignment_summary.tsv` | Per-sample alignment statistics |
| `HPV16_NC_001526.4.fa` | Downloaded reference genome |
| `HPV16_NC_001526.4.mmi` | Minimap2 index |

### Phase 5A: Population Profiling (Superseded)

**`Phase5A_Population_Consolidation.py`** and **`Phase5A_Revised_Population_Discovery.py`** defined the original two-population model ("Pop1_Mutagenic" / "Pop2_Stealth") using multi-axis concordance scoring across A3A dominance, HPV lifecycle stage, CNV status, and differentiation. These scripts are superseded by the three-group selection in Figure 4 Step00B, which defines SBS2-HIGH, CNV-HIGH, and NORMAL populations using composite scores optimized for the network analysis. Scripts are retained for provenance.

**`Phase5A_v2_DEG_Analysis.py`** performed differential expression and KEGG GSEA between the old two-population model. Key findings (spliceosome enrichment in CNV-HIGH, immune/allograft pathways in SBS2-HIGH) are consistent with the updated lifecycle analysis.

### Lifecycle Diagnostic

**`TROUBLESHOOTING/Diagnostic_HPV_Lifecycle_Markers.py`** is the primary diagnostic that informed the Figure 6 redesign. It runs three analyses using the Figure 4 three-group populations:

**Diagnostic 1** compares HPV16 gene expression across the four lifecycle phases (maintenance, amplification, oncogene, capsid) in both raw counts and normalized fractions. The normalization (gene reads / total HPV reads per cell) controls for the two-fold viral load difference between populations and revealed that SBS2-HIGH cells dedicate more reads to maintenance genes while CNV-HIGH cells dedicate more to capsid genes, with identical oncogene (E6/E7) fractions.

**Diagnostic 2** computes integration proxy metrics: E2:(E6+E7) ratio, E2 fraction of total reads, E6E7 fraction of total reads, L1L2 fraction, total HPV16 reads (copy number proxy), and the early/late ratio. The E2 fraction being higher in SBS2-HIGH (7.0% vs 4.2%, p=0.052) argues against integration with E2 disruption and supports intact E2 regulation of the maintenance lifecycle.

**Diagnostic 3** profiles 53 host marker genes organized by biological process: transformation indicators (CDKN2A, MCM7, CCNE1, MKI67), p53/Rb readouts (CDKN1A, MDM2), ATM cascade (ATM, CHEK2, BRCA1, MRE11, H2AX), G2/M arrest (CDC25C, CDK1, CCNB1), immune signaling (HLA-A/B/C, B2M, TAP1), differentiation (IVL, KRT5, KRT14), and APOBEC family (A3A-H, CGAS, STING1). Gene selection was informed by Doorbar et al. (2012), Moody and Laimins (2009), and McBride and Warburton (2017).

### Figure Generation

**`Generate_Figure6_Lifecycle_Panels.py`** produces the four main Figure 6 panels:

| Panel | Content |
|-------|---------|
| A | UMAP of epithelial cells colored by total HPV16 reads (log1p, magma) with population boundary indicators |
| B | Grouped bar chart of normalized HPV16 gene fractions by lifecycle phase across three populations |
| C | Dot plot of host marker genes: ATM/productive replication (left) vs immune/differentiation (right) |
| D | Violin + box plots of A3A and A3B expression by population |

**`Generate_Figure6_Panels.py`** is the pre-lifecycle reframe version of the figure generation script. Retained for reference but superseded by the lifecycle panels script above.

### Utility Scripts

**`parse_chimeric_junctions.py`** parses STAR chimeric junction output files. Used by the fusion analysis pipeline in `scripts/NEOANTIGEN/Step05_Fusion_Analysis.py`.

**`sample_list.txt`** contains the 31 GSE173468 sample accessions used across all phases.

**`RUN_HPV_ANALYSIS.sh`** is the SLURM master runner for Phases 1-4.

---

## Output Data Structure

```
data/FIG_6/
├── 00_diagnostic_inventory/          # Phase 1: per-cell profiling
│   ├── basal_cell_master_table.tsv
│   └── diagnostic_report.txt
├── 00_input/
│   └── adata_v_pp.h5ad               # Viral detection object
├── 01_raw_hpv16_counts/              # Phase 2: raw UMI counts
│   ├── basal_cell_master_table_with_raw_HPV16.tsv
│   ├── raw_HPV16_counts_all_cells.tsv
│   └── per_sample_hpv16_summary.tsv
├── 03_hpv16_genome/                  # Phase 4: genome alignment
│   ├── per_cell_hpv16_gene_counts.tsv   ← PRIMARY OUTPUT (76,978 cells)
│   ├── HPV16_NC_001526.4.fa
│   └── per_sample_alignment_summary.tsv
├── DIAGNOSTIC_LIFECYCLE_MARKERS/     # Lifecycle diagnostic outputs
│   ├── diagnostic_lifecycle_report.txt
│   ├── hpv16_gene_by_phase_and_population.tsv
│   ├── hpv16_gene_fractions_normalized.tsv
│   ├── integration_proxy_metrics.tsv
│   ├── host_marker_expression_summary.tsv
│   └── host_marker_per_cell_values.tsv
└── FIGURE_6_PANELS/                  # Figure 6 panels
    ├── Panel_6A_UMAP_HPV16_viral_load.pdf/.png
    ├── Panel_6B_HPV_lifecycle_fractions.pdf/.png
    ├── Panel_6C_host_marker_dotplot.pdf/.png
    ├── Panel_6D_A3A_vs_A3B_violins.pdf/.png
    └── Figure6_composite.pdf/.png
```

---

## Key Results

| Finding | Evidence |
|---------|----------|
| CNV-HIGH has 2x viral load | 111.8 vs 57.0 total HPV16 reads per cell (p=8.45e-11) |
| SBS2-HIGH is in maintenance phase | 25.7% of reads to E1/E2 vs 14.1% in CNV-HIGH (p=6.42e-14) |
| CNV-HIGH is in productive phase | 23.4% of reads to L1/L2 vs 9.1% in SBS2-HIGH (p=1.05e-60) |
| E6/E7 dosage is identical | 0.77% vs 0.95% of total reads (not the driver) |
| E2 is intact in SBS2-HIGH | E2 fraction 7.0% vs 4.2% (p=0.052, argues against integration) |
| ATM cascade in CNV-HIGH | CHEK2, BRCA1, MRE11, H2AX, CDC25C, CDK1, CCNB1 all enriched |
| Immune visibility in SBS2-HIGH | HLA-A (p=5e-39), HLA-B (p=1e-26), B2M (p=1e-48) |
| Differentiation in SBS2-HIGH | IVL (p=6.75e-66), CDH1 (p=6.58e-12) |
| p53 active in SBS2-HIGH | CDKN1A/p21 (p=7.57e-8), MDM2 (p=3.12e-7) |
| A3A dominates SBS2-HIGH | 6.46 vs 3.40 (p=1.44e-90) |
| A3B dominates CNV-HIGH | 4.23 vs 2.21 (p=8.27e-36) |

### Biological Model

SBS2 and CNV represent sequential consequences of HPV lifecycle progression. During E2-regulated maintenance, immune surveillance activates interferon-driven A3A mutagenesis through the co-expression network identified in Figure 4, and the resulting point mutations accumulate as SBS2. As cells transition to productive infection, ATM-dependent DNA damage repair corrects point mutations while replication stress generates chromosomal instability. A cell's position in the HPV lifecycle determines both the type of APOBEC enzyme active (A3A vs A3B) and the class of genomic damage that accumulates (SBS2 vs CNV).

---

## Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| minimap2 | 2.28 | HPV16 genome alignment |
| pysam | 0.22+ | BAM file manipulation and CB tag extraction |
| scanpy | 1.9.6 | Differential expression, gene extraction |
| gseapy | 1.1.1 | KEGG gene set enrichment |

Conda environment: `NETWORK` (shared with Figure 4 network analysis)

---

## Related Pipelines

The neoantigen analysis (formerly Phase 5B in this directory) has been refactored into a standalone six-step pipeline:

- Location: `scripts/NEOANTIGEN/`
- Documentation: `scripts/NEOANTIGEN/NEOANTIGEN_WALKTHROUGH.md`
- Pipeline: Step01 (prep) through Step06 (integration) plus figure generation
- Produces: Figure 7 (neoantigen landscape)
