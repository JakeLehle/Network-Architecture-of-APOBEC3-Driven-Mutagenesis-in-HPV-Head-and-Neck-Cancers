# Question 3 Workflow

## Question 3: At single-cell resolution, are SBS2 mutational signatures localized to a specific cell type in HNSCC?

### Rationale

Figure 2 established that differential co-expression network analysis on A3-matched bulk HNSCC tumors identifies 13 gene communities associated with SBS2 mutagenic activity. Two communities were of particular biological interest: Community 2 (160 genes) containing A3B alongside KRT5, a canonical basal epithelial cell marker, and Community 0 (288 genes) containing A3G with an immune-associated anti-correlation pattern. However, only 2 of 24 canonical cell-type markers appeared in any community, meaning bulk data fundamentally could not resolve whether these communities reflect cell-type-specific programs. This question uses single-cell RNA-seq to test two specific predictions from the bulk analysis: (1) that SBS2 mutagenesis originates specifically from basal epithelial cells (predicted by Community 2: A3B + KRT5), and (2) that A3C and A3G expression originate from immune cell populations rather than tumor cells (predicted by the bystander analysis in Figure 1 and the immune-associated Community 0 in Figure 2).

All single-cell data were processed using two custom pipelines: **SRAscraper** (https://github.com/JakeLehle/SRAscraper) for automated data acquisition from the Sequence Read Archive, and **ClusterCatcher** (v1.3.0; https://github.com/JakeLehle/ClusterCatcher) for end-to-end single-cell analysis including Cell Ranger alignment, adaptive QC filtering, popV cell type annotation, cancer cell detection (CytoTRACE2 + inferCNV), SComatic somatic mutation calling, and semi-supervised mutational signature deconvolution. For a complete description of these pipelines, see their respective GitHub repositories. The scripts in this repository focus on the figure generation steps specific to the paper.

### Data Sources

| Source | Description |
|--------|-------------|
| [SRA / GEO](https://www.ncbi.nlm.nih.gov/geo/) | HNSCC scRNA-seq from tumor and adjacent normal tissue (downloaded via SRAscraper) |
| [COSMIC SBS Signatures](https://cancer.sanger.ac.uk/signatures/) | COSMIC v3.4 SBS reference signatures (`COSMIC_v3.4_SBS_GRCh38.txt`) |
| [10x Genomics](https://support.10xgenomics.com/) | GRCh38 Cell Ranger reference genome (`refdata-gex-GRCh38-2020-A`) |
| [Tabula Sapiens](https://huggingface.co/popV/tabula_sapiens_All_Cells) | Pre-trained popV HubModel for cell type annotation |

### Directory Structure

```
scripts/SINGLE_CELL/
├── RUN_CLUSTERCATCHER_PIPELINE.sh           # SLURM batch script — install, configure, run ClusterCatcher
├── Step01_Generate_Figure3_Panels.py        # Step 1 — Generate Figure 3b, 3c, Supplemental Figure 2
│
├── ClusterCatcher/                          # Full ClusterCatcher pipeline (submodule/copy)
│   ├── README.md                            #   Comprehensive pipeline documentation
│   ├── environment.yml                      #   Conda environment specification
│   ├── pyproject.toml                       #   Package metadata
│   ├── setup.py                             #   Installation script
│   ├── cli/                                 #   CLI commands (sample-information, create-config, run-config)
│   └── snakemake_wrapper/                   #   Snakemake pipeline
│       ├── Snakefile                        #     Master workflow definition
│       ├── config.yaml                      #     Configuration template
│       ├── scripts/                         #     Per-module Python scripts
│       │   ├── scanpy_qc_annotation.py      #       QC filtering + popV annotation
│       │   ├── cancer_cell_detection.py     #       CytoTRACE2 + inferCNV consensus
│       │   ├── kraken2_viral_detection.py   #       Viral pathogen detection (optional)
│       │   ├── scomatic_mutation_calling.py #       SComatic somatic variant calling
│       │   └── signature_analysis.py        #       NNLS signature deconvolution
│       ├── envs/                            #     Per-rule conda environments
│       └── rules/                           #     Modular Snakemake rules
│
└── TROUBLESHOOTING/                         # Original development scripts (not part of pipeline)
    ├── 2025-05-29_SC_Cluster_Annotation.py  #   Original popV annotation + cluster refinement
    ├── 2025-06-26_HNC_SC_Cancer_Detection.py#   Original CytoTRACE2 + inferCNV script
    ├── 2025-12-07_SigProfiler.py            #   Original signature deconvolution + SBS2 analysis
    ├── 2025-12-07_HNC_Virus_FOR_SABCS.py    #   Original HPV analysis + UMAP plotting
    ├── 2026-02-20_HNSCC_Get_Marker_Genes.py #   Marker gene validation (supplemental)
    ├── Cell_ranger_following_SRAscraper_download.py
    ├── mg2sc_all_virus.py
    └── SComatic_script_tmp_v3.py
```

### External Dependencies

| Software | Version | Purpose | Installation |
|----------|---------|---------|--------------|
| **ClusterCatcher** | >=1.3.0 | End-to-end scRNA-seq pipeline | `git clone https://github.com/JakeLehle/ClusterCatcher.git && pip install -e .` |
| **SRAscraper** | latest | Automated SRA data download | `git clone https://github.com/JakeLehle/SRAscraper.git && pip install -e .` |
| **Cell Ranger** | >=7.0 | FASTQ alignment & counting | [10x Genomics](https://support.10xgenomics.com/) |
| **SComatic** | latest | Somatic mutation calling | `git clone https://github.com/cortes-ciriano-lab/SComatic.git` |
| **CytoTRACE2** | latest | Stemness/potency scoring | `git clone https://github.com/digitalcytometry/cytotrace2.git` |
| **Snakemake** | >=7.0 | Pipeline orchestration | `conda install -c bioconda snakemake` |

---

### Pipeline Overview

Unlike Questions 1 and 2, which use step-by-step R/Python scripts, Question 3 relies on the ClusterCatcher Snakemake pipeline for all data processing. The pipeline produces a fully annotated AnnData object (`adata_final.h5ad`) containing cell type annotations, cancer/normal cell status, somatic mutations, and per-cell signature weights. A single post-pipeline Python script then generates the remaining figure panels not produced by ClusterCatcher.

**SLURM Job: `RUN_CLUSTERCATCHER_PIPELINE.sh`**
- Verifies all dependencies
- Generates pipeline configuration with all analysis parameters
- Executes the 7-module Snakemake pipeline
- Verifies output completeness

**Post-pipeline: `Step01_Generate_Figure3_Panels.py`**
- Loads `adata_final.h5ad`
- Generates Figure 3b, Figure 3c, and Supplemental Figure 2

---

### ClusterCatcher Pipeline Modules (executed via SLURM script)

The following modules are executed sequentially by Snakemake. Full documentation is in the [ClusterCatcher README](https://github.com/JakeLehle/ClusterCatcher).

#### Module 1: Cell Ranger Alignment

**Script:** `snakemake_wrapper/scripts/cellranger_count.py`

Aligns FASTQ files to GRCh38 using Cell Ranger, producing filtered feature-barcode matrices and position-sorted BAM files with cell barcode (CB) tags.

**Key parameters:**
- Reference: `refdata-gex-GRCh38-2020-A`
- Chemistry: auto-detect
- Expected cells: 10,000
- BAM generation: enabled (required for SComatic)

**Output:** `cellranger/{sample}/outs/` — filtered matrices, BAM files, web summaries

---

#### Module 2–3: QC, Filtering, and Cell Type Annotation

**Script:** `snakemake_wrapper/scripts/scanpy_qc_annotation.py`

Performs adaptive QC filtering, doublet removal, and automated cell type annotation with cluster-level refinement.

**Key parameters:**
- QC mode: adaptive (MAD-based outlier detection)
- MAD thresholds: 5 MADs for total counts and genes detected
- Maximum mitochondrial: 20%
- Minimum genes per cell: 200
- Minimum cells per gene: 20
- Doublet detection: Scrublet (rate = 0.08)
- Annotation model: popV with Tabula Sapiens All Cells (inference mode)
- Normalization: CPM + log1p (post-annotation)
- Clustering: Leiden (resolution = 1.0, seed = 42)
- Final annotation: cluster-level assignment using highest summed normalized popV confidence score

**Output:**
- `annotation/adata_annotated.h5ad` — annotated AnnData
- `annotation/figures/UMAP_popv_prediction.pdf` → **Supplemental Figure 1**
- `annotation/figures/UMAP_clusters.pdf` → **Figure 3a** (cluster numbers)
- `annotation/figures/stacked_bar_cluster_composition.pdf` → **Figure 3a** (stacked bar)
- `annotation/figures/UMAP_final_annotation.pdf` → **Figure 3a** (final cell types)

---

#### Module 4: Cancer Cell Detection

**Script:** `snakemake_wrapper/scripts/cancer_cell_detection.py`

Identifies cancer vs. normal cells using dual-model consensus of CytoTRACE2 stemness scoring and inferCNV copy number variation analysis.

**Key parameters:**
- CytoTRACE2: species = human, max 200,000 cells per chunk
- inferCNV: window size = 250 genes
- Agreement: alpha = 0.5, minimum Spearman rho = 0.5

**Output:**
- `dysregulation/adata_cancer_detected.h5ad` — AnnData with `Final_cancer_cell_status`
- `dysregulation/figures/` — CytoTRACE2, inferCNV, and agreement visualizations

---

#### Module 5: Viral Detection (Optional)

**Script:** `snakemake_wrapper/scripts/kraken2_viral_detection.py`

Detects viral pathogens (including HPV16) in unmapped reads using Kraken2. Produces per-cell organism count matrices for downstream integration.

**Output:** `viral/viral_detection_summary.tsv`, per-cell viral counts

---

#### Module 6: Somatic Mutation Calling

**Script:** `snakemake_wrapper/scripts/scomatic_mutation_calling.py`

Calls somatic mutations at single-cell resolution using SComatic. Reads are pooled across cells of the same annotated cell type to achieve sufficient depth for variant calling against GRCh38.

**Key parameters:**
- Minimum coverage: 5 reads
- Minimum cells: 5
- Minimum base quality: 30
- Minimum mapping quality: 30
- Germline filtering: variants in >1 cell type removed
- Additional filtering: Panel of Normals, RNA editing sites, mappable regions BED

**Output:**
- `mutations/all_samples.single_cell_genotype.filtered.tsv` — filtered somatic mutations
- `mutations/CombinedCallableSites/complete_callable_sites.tsv` — callable sites per cell
- `mutations/cell_annotations.tsv` — cell type assignments for SComatic

---

#### Module 7: Mutational Signature Deconvolution

**Script:** `snakemake_wrapper/scripts/signature_analysis.py`

Converts somatic mutations to a 96-trinucleotide-context matrix and deconvolves into COSMIC reference signatures using semi-supervised NNLS refitting.

**Key parameters:**
- COSMIC version: v3.4 (GRCh38)
- HNSCC signature set: SBS1, SBS2, SBS4, SBS5, SBS7a, SBS7b, SBS13, SBS16, SBS17a, SBS17b, SBS18, SBS29, SBS39, SBS40, SBS44
- Core signatures (always retained): SBS2, SBS13, SBS5
- Signature selection: scree plot elbow detection on reconstruction error (Frobenius norm)
- Normalization: mutations per cell normalized by callable sites (depth ≥ 5)

**Output:**
- `signatures/signature_weights_per_cell.txt` — per-cell signature weights (rows: signatures, columns: cells)
- `signatures/adata_final.h5ad` — master AnnData with all annotations, mutations, and signature weights
- `signatures/mutation_matrix_96context.txt` — 96-context mutation count matrix

---

### Step 1: Generate Figure 3 Panels (Post-Pipeline)

**Script:** `Step01_Generate_Figure3_Panels.py`

Loads the final AnnData object from ClusterCatcher and generates the figure panels that require custom visualization beyond what the pipeline produces. Run interactively or via SLURM after the ClusterCatcher pipeline completes.

**Dependencies:** `scanpy`, `matplotlib`, `numpy`, `pandas`, `scipy`

**Input:**
- `adata_final.h5ad` from ClusterCatcher Module 7
- `COSMIC_v3.4_SBS_GRCh38.txt`
- `signature_weights_per_cell.txt`
- Mutation matrix (96-context)

**Output (→ `data/FIG_3/`):**
- `Panel_3b_Normalized_Mutations_SBS2_A3A_A3B.pdf/.png` — Figure 3b
- `Panel_3c_Aggregated_SBS2_vs_COSMIC.pdf/.png` — Figure 3c
- `Supplemental_Figure_2_A3_Family_Expression.pdf/.png` — Supplemental Figure 2

**Figure 3b — Mutation burden and APOBEC3 expression (4-panel UMAP):**
- Panel 1: Normalized somatic mutation count per cell (total mutations ÷ callable sites with depth ≥ 5). Basal epithelial cells show dramatically elevated mutation burden relative to all other cell types.
- Panel 2: SBS2 signature weight per cell from NNLS refitting. High SBS2 co-localizes with high-mutation basal cells, confirming APOBEC-driven mutagenesis is basal-cell-specific.
- Panel 3: APOBEC3A expression. Enriched in the basal cell cluster, concentrated in a specific subregion.
- Panel 4: APOBEC3B expression. Also enriched in basal cells but partially distinct subregion. The overlap zone between A3A and A3B corresponds to the cells with the highest mutation burden and SBS2 weights.
- Prints per-panel summary statistics (basal vs. other cells) to console for reporting.

**Figure 3c — Aggregated signature validation (96-context bar plot):**
- Selects basal cells in the top 20th percentile of SBS2 weight.
- Aggregates their 96-trinucleotide-context mutation profiles by summing across all selected cells.
- Normalizes to proportions and compares to the pure COSMIC SBS2 reference signature.
- Two stacked bar plots: aggregated profile (top) and COSMIC SBS2 (bottom).
- Computes and displays cosine similarity in the panel title.
- Standard COSMIC 6-class coloring (C>A, C>G, C>T, T>A, T>C, T>G).

**Supplemental Figure 2 — A3 family expression (7-panel UMAP):**
- 2×4 grid of UMAP embeddings colored by expression of each A3 family member: A3A, A3B, A3C, A3D, A3F, A3G, A3H.
- Gene expression extracted by symbol lookup (handles gene_symbol, feature_name, or direct var_names indexing).
- Confirms A3A/A3B in basal cells, A3C/A3G in immune populations, A3H low/diffuse.
- Prints per-gene non-zero cell counts to console.

---

### Figure 3 Summary

**Title:** Single-cell resolution localizes SBS2 mutagenesis to basal epithelial cells co-expressing A3A and A3B.

| Panel | Content | Source | Step |
|-------|---------|--------|------|
| a | Stacked bar (cluster composition) + cluster UMAP + final annotation UMAP | ClusterCatcher Module 2–3 | Pipeline |
| b | 4-panel UMAP: normalized mutations, SBS2, A3A, A3B | `Step01_Generate_Figure3_Panels.py` | 1 |
| c | Aggregated 96-context signature vs. COSMIC SBS2 | `Step01_Generate_Figure3_Panels.py` | 1 |

| Supplemental | Content | Source | Step |
|--------------|---------|--------|------|
| Supp. Fig. 1 | Per-cell popV prediction UMAP | ClusterCatcher Module 2–3 | Pipeline |
| Supp. Fig. 2 | A3 family expression across cell types (7 UMAPs) | `Step01_Generate_Figure3_Panels.py` | 1 |

**Narrative arc:**
1. 12 cell populations resolved from HNSCC scRNA-seq; cluster-level annotation validated by popV confidence scores (panel a, Supplemental Figure 1)
2. Basal epithelial cells carry dramatically elevated normalized somatic mutation burden (panel b)
3. SBS2 signature weight co-localizes with high-mutation basal cells — APOBEC mutagenesis is basal-cell-specific (panel b)
4. A3A and A3B are both enriched in basal cells with partially distinct subregions; the overlap zone corresponds to the highest SBS2 (panel b)
5. This validates the bulk TCGA additive model (Fig. 1e,f) at single-cell resolution and localizes it to a specific subpopulation
6. A3C and A3G localize to immune populations, confirming bystander/immune-compartment signals from Figures 1–2 (Supplemental Figure 2)
7. Aggregated basal cell mutation profile matches COSMIC SBS2 reference, confirming genuine APOBEC activity (panel c)
8. **Motivates Question 4:** apply differential network analysis directly to single-cell basal epithelial expression data to recapitulate and refine cofactor communities at single-cell resolution

### Key Results

| Metric | Value |
|--------|-------|
| Cell populations resolved | 12 |
| Annotation method | popV (Tabula Sapiens) + cluster-level refinement |
| Cancer detection | CytoTRACE2 + inferCNV dual-model consensus |
| Somatic mutation method | SComatic (cell-type-pooled, germline-filtered) |
| Signature method | Semi-supervised NNLS against COSMIC v3.4 |
| Core signatures | SBS2, SBS13, SBS5 |
| Cell type with highest mutation burden | Basal epithelial cells |
| SBS2 enrichment | Basal epithelial cells (co-localized with A3A+A3B) |
| A3A expression | Basal cell cluster (subregion) |
| A3B expression | Basal cell cluster (partially distinct subregion) |
| A3C expression | Immune cell populations (TILs) |
| A3G expression | Immune cell populations |
| A3H expression | Low/diffuse across all cell types |

### Troubleshooting Scripts

The following scripts in `scripts/SINGLE_CELL/TROUBLESHOOTING/` are not part of the main pipeline but were used during development. They contain the original analysis code from which ClusterCatcher modules were derived:

| Script | Purpose |
|--------|---------|
| `2025-05-29_SC_Cluster_Annotation.py` | Original popV annotation + Leiden clustering + weighted cluster-based refinement |
| `2025-06-26_HNC_SC_Cancer_Detection.py` | Original CytoTRACE2 + inferCNV cancer cell detection |
| `2025-12-07_SigProfiler.py` | Original semi-supervised signature deconvolution + SBS2 post-processing + control cell selection |
| `2025-12-07_HNC_Virus_FOR_SABCS.py` | Original HPV16 analysis + APOBEC expression UMAPs + differential expression |
| `2026-02-20_HNSCC_Get_Marker_Genes.py` | Marker gene validation (Wilcoxon rank-sum DE per cell type) |
| `SComatic_script_tmp_v3.py` | Original multi-phase SComatic mutation calling pipeline |
