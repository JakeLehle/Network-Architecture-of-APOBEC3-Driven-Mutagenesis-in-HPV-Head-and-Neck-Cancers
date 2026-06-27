# Question 1 Workflow

## Question 1: At single-cell resolution, where does A3-driven mutagenesis localize, and do A3A and A3B drive divergent mutational programs?

### Rationale

A3A and A3B are necessary but not sufficient for SBS2, and bulk tissue cannot say why because it averages over two layers of heterogeneity at once. First, A3A and A3B are not confined to tumor cells. A3A is expressed in macrophages and dendritic cells and A3B in plasma and CD10+ B cells, so a signature measured in bulk cannot be cleanly attributed to the cells that produced it. Second, even a pure population of tumor epithelial cells is not uniform, because the infecting virus sits at a different phase of its lifecycle from one cell to the next. As a result, no tumor presents an uncontaminated SBS2 signal in bulk, and the comparison against normal tissue cannot be made there.

This question moves to single-cell resolution to resolve both layers. Using the HPV16-positive HNSCC dataset GSE173468, with tumor and matched normal-adjacent tissue, we call somatic mutations and assign SBS2 to each cell, then ask two things: whether A3-driven mutagenesis localizes to a specific epithelial compartment, and whether A3A and A3B drive divergent classes of genomic damage within that compartment, point mutations (SBS2) versus copy-number change (CNV). This section establishes the localization and the divergence that the rest of the paper builds on.

### Data Sources

| Source | Description |
|--------|-------------|
| [SRA / GEO](https://www.ncbi.nlm.nih.gov/geo/) | GSE173468, HPV16-positive HNSCC scRNA-seq from tumor and matched normal-adjacent tissue (downloaded via SRAscraper) |
| [COSMIC SBS Signatures](https://cancer.sanger.ac.uk/signatures/) | COSMIC v3.4 SBS reference signatures (`COSMIC_v3.4_SBS_GRCh38.txt`) |
| [10x Genomics](https://support.10xgenomics.com/) | GRCh38 Cell Ranger reference genome (`refdata-gex-GRCh38-2020-A`) |
| [Tabula Sapiens](https://huggingface.co/popV/tabula_sapiens_All_Cells) | Pre-trained popV HubModel for cell type annotation |

### Directory Structure

```
scripts/SINGLE_CELL/
├── ClusterCatcher/                          # Full pipeline (engine for Figs 1, 2, and early supplementals)
│   ├── README.md                            #   Comprehensive pipeline documentation (parameters live here)
│   ├── environment.yml
│   ├── cli/                                 #   sample-information, create-config, run-config
│   └── snakemake_wrapper/                   #   Snakefile, per-module scripts, per-rule envs, rules
├── Run_Cluster_Catcher_Pipeline.sh          # SLURM runner: install, configure, and execute ClusterCatcher
├── Generate_Figure3_A3_Dominance_Slider.py  # Figure 3: A3A/A3B dominance divergence slider
├── Step02_Supplemental_Marker_Validation.py # Supplemental: classical-marker validation of popV annotation
├── Step01_Generate_Figure3_Panels.py        # Panel helper; its UMAP/COSMIC panels feed the ClusterCatcher-derived Fig 1/2 assembly (not a defined figure step)
└── TROUBLESHOOTING/                          # Diagnostics and development scripts (not documented here)
```

### External Dependencies

| Software | Version | Purpose | Installation |
|----------|---------|---------|--------------|
| **ClusterCatcher** | >=1.3.0 | End-to-end scRNA-seq pipeline | `git clone https://github.com/Diako-Lab/ClusterCatcher.git && pip install -e .` |
| **SRAscraper** | latest | Automated SRA data download | `git clone https://github.com/Diako-Lab/SRAscraper.git && pip install -e .` |
| **Cell Ranger** | >=7.0 | FASTQ alignment and counting | [10x Genomics](https://support.10xgenomics.com/) |
| **SComatic** | latest | Somatic mutation calling | `git clone https://github.com/cortes-ciriano-lab/SComatic.git` |
| **CytoTRACE2** | latest | Stemness/potency scoring | `git clone https://github.com/digitalcytometry/cytotrace2.git` |
| **Snakemake** | >=7.0 | Pipeline orchestration | `conda install -c bioconda snakemake` |

---

### Pipeline Overview

This question is engine-plus-two-analyses. ClusterCatcher does all of the heavy data processing and produces a single fully annotated AnnData object, `adata_final.h5ad`, that carries cell-type annotations, cancer/normal status, per-cell somatic mutations, and per-cell signature weights. Figures 1 and 2 and the early supplementals (cell-type annotation, A3 family localization, the mutation, SBS2, A3A, A3B and CNV UMAPs, and the aggregated COSMIC comparison) are assembled from these ClusterCatcher outputs separately, so there is no single repo script that "makes Figure 1" or "makes Figure 2."

Two bespoke analyses live in this directory as standalone main-dir scripts:

- **`Generate_Figure3_A3_Dominance_Slider.py`** produces Figure 3, the within-tissue SBS2-versus-CNV divergence axis split by A3A versus A3B dominance.
- **`Step02_Supplemental_Marker_Validation.py`** produces the supplemental classical-marker validation of the popV annotations.

`Step01_Generate_Figure3_Panels.py` remains in the directory as a panel helper. Its UMAP and COSMIC panels feed the ClusterCatcher-derived Figure 1 and 2 assembly rather than standing as a defined figure step, so it is not written up as one below.

---

### ClusterCatcher: the engine

ClusterCatcher (v1.3.0) runs as a seven-module Snakemake pipeline launched by `Run_Cluster_Catcher_Pipeline.sh`. Full parameter documentation is in the [ClusterCatcher README](https://github.com/Diako-Lab/ClusterCatcher); the condensed view of what each module contributes to this question is below.

| Module | Script | Produces |
|--------|--------|----------|
| 1. Alignment | `cellranger_count.py` | GRCh38 alignment, filtered feature-barcode matrices, CB-tagged BAMs |
| 2-3. QC + annotation | `scanpy_qc_annotation.py` | Adaptive MAD QC, Scrublet doublet removal, popV annotation with cluster-level refinement; annotation and composition UMAPs that feed Figure 1 and Supplemental Figure 1 |
| 4. Cancer detection | `cancer_cell_detection.py` | CytoTRACE2 + inferCNV dual-model consensus; `cnv_score` and stemness per cell, feeding Supplemental Figure 5 |
| 5. Viral detection | `kraken2_viral_detection.py` | Per-cell organism counts including HPV16 (used downstream in Question 4) |
| 6. Mutation calling | `scomatic_mutation_calling.py` | Cell-type-pooled, germline-filtered somatic variants, per-cell callable sites |
| 7. Signature deconvolution | `signature_analysis.py` | 96-context matrix, semi-supervised NNLS refit against COSMIC v3.4, per-cell signature weights; writes `adata_final.h5ad` |

Key processing choices that matter for interpretation: somatic calls are pooled within each annotated cell type and normalized to each cell's callable sites (depth >= 5); signature refitting always retains SBS2, SBS13, and SBS5 and adds other HNSCC-relevant COSMIC signatures by scree-plot elbow on reconstruction error; the basal annotation in `adata.obs` is lowercase `basal cell`.

**Engine outputs used in this question:** `signatures/adata_final.h5ad` and `signatures/signature_weights_per_cell.txt`.

---

### Figure 3: A3 dominance divergence slider

**Script:** `Generate_Figure3_A3_Dominance_Slider.py`

This is the one bespoke analysis behind a main figure in this section. It tests whether the A3A-to-SBS2 and A3B-to-CNV link, which is visible in the Figure 2 UMAPs but weak at the level of global expression, sharpens when cells are split by which enzyme dominates.

**Approach:**
- Restrict to basal cells expressing either enzyme (A3A + A3B > 0).
- Define an A3A dominance fraction, A3A / (A3A + A3B); cells above 0.5 are A3A-dominant, below 0.5 A3B-dominant.
- Place each cell on a within-tissue axis, z(SBS2) minus z(CNV), where SBS2 and the inferCNV `cnv_score` are each standardized within that tissue's expressing-either population. Standardizing within tissue keeps "center" at each tissue's own median so the tumor/normal marginal differences do not push normal off-center artefactually.
- Plot four rows: A3A-dominant tumor, A3B-dominant tumor, A3A-dominant normal-adjacent, A3B-dominant normal-adjacent. The right pole is SBS2-high, the left pole is CNV-high (productive).

**Design notes carried in the script header (worth keeping in mind for the text):**
- Stemness (CytoTRACE2) is deliberately not blended into the productive pole, because the A3A fraction anti-correlates with stemness in both tumor and normal, so stemness is not tumor-specific. The productive pole is CNV alone.
- The per-enzyme A3 > 0 conditioning is gone; conditioning on expressing cells deleted the co-occurrence signal. Dominance is defined among cells expressing either enzyme.
- Absolute A3 level is not the axis driver. SBS2 does not separate the two enzymes at the level of absolute expression; dominance does.

**Input:** `data/FIG_4/00_input/adata_final.h5ad`, `data/FIG_4/00_input/signature_weights_per_cell.txt`. `TARGET_CELL_TYPE = "basal cell"`.

**Output (-> `data/FIG_3/figures/`):**
- `Figure3_A3_dominance_slider.pdf` / `.png` (300 DPI). SBS2 pole `#ed6a5a`, CNV pole `#F6D155`.

**Honest caveat for the text:** the separation is carried mainly by CNV, low CNV on the A3A-dominant side, rather than by high SBS2, since SBS2 is shared between the two enzymes. This fits the maintenance (low-CNV) versus productive (high-CNV) framing that Question 4 develops.

---

### Supplemental: classical-marker annotation validation

**Script:** `Step02_Supplemental_Marker_Validation.py`

Validates the popV cell-type annotations with classical markers, providing an audit trail for how the 12 populations were assigned.

**Approach:**
1. Load `adata_final.h5ad` and confirm `final_annotation` in `.obs`.
2. Run Wilcoxon rank-sum differential expression to find the top 20 markers per cell type.
3. Save the full marker table; log where each curated marker ranks within the top 20.
4. Generate a UMAP grid of two curated markers per cell type.

**Input:** `adata_final.h5ad`.

**Output (-> `data/FIG_3/`):**
- `Supplemental_Marker_Genes_Top20_Per_CellType.tsv`, `selected_marker_genes.tsv`
- `Supplemental_Figure_Marker_Validation.pdf` / `.png` (Supplemental Figure 3)

**Curated markers (12 cell types x 2):**

| Cell Type | Marker 1 | Marker 2 |
|-----------|----------|----------|
| CD4+ T cell | *IL7R* | *CD3E* |
| B cell | *CD79A* | *MS4A1* |
| CD8+ T cell | *CD8A* | *CCL5* |
| Regulatory T cell | *TIGIT* | *CTLA4* |
| Macrophage | *CD68* | *SPI1* |
| Myeloid DC | *LAMP3* | *CCR7* |
| Plasmacytoid DC | *IL3RA* | *LILRA4* |
| Mast cell | *TPSAB1* | *CPA3* |
| Fibroblast | *DCN* | *COL1A1* |
| Smooth muscle | *TAGLN* | *RGS5* |
| Basal cell | *KRT5* | *TACSTD2* |
| Endothelial cell | *PECAM1* | *VWF* |

---

### Figure and Supplemental Summary

| Figure | Content | Source |
|--------|---------|--------|
| Fig 1 | Cell-type annotation and A3A/A3B localization to the basal compartment | Assembled from ClusterCatcher annotation outputs |
| Fig 2 | A3-association multipanel: normalized mutations, SBS2, A3A, A3B UMAPs; SBS2 and CNV associations; aggregated 96-context versus COSMIC SBS2 | Assembled from ClusterCatcher outputs (mutations, weights, expression, inferCNV) |
| Fig 3 | A3A/A3B dominance divergence slider (within-tissue SBS2-vs-CNV axis) | `Generate_Figure3_A3_Dominance_Slider.py` |
| Supp Fig 1 | Per-cell popV prediction UMAP | ClusterCatcher annotation |
| Supp Fig 2 | A3 family expression across cell types | Assembled from ClusterCatcher outputs |
| Supp Fig 3 | Classical-marker validation UMAPs (12 cell types x 2) | `Step02_Supplemental_Marker_Validation.py` |
| Supp Fig 5 | CytoTRACE2 and inferCNV per cell | ClusterCatcher cancer detection |

**Narrative arc:**
1. Twelve cell populations resolve from the HPV16-positive HNSCC data; cluster-level popV annotation is validated by classical markers (Fig 1, Supplemental Figures 1 and 3).
2. A3A and A3B expression concentrates almost exclusively in basal epithelial cells, with only low-level A3A in macrophages and myeloid dendritic cells (Fig 1, Fig 2a-c).
3. Among basal cells, SBS2 burden associates with A3A more than A3B, but the global association is weak, and the aggregated trinucleotide profile of the highest-SBS2 cells matches the COSMIC SBS2 reference (Fig 2d-f).
4. A subset of SBS2-positive basal cells carries no detectable A3A or A3B at capture, consistent with pulsatile A3 induction.
5. Elevated CNV and stemness map instead to an A3B-associated basal subpopulation, so A3A-linked point mutations and A3B-linked chromosomal instability occupy separable subregions of the same compartment (Fig 2g-h).
6. Splitting basal cells by A3 dominance and placing them on a within-tissue SBS2-vs-CNV axis resolves the divergence that global expression obscures: A3A-dominant cells shift toward SBS2, A3B-dominant cells toward CNV, and only in tumor (Fig 3).
7. **Motivates Question 2:** apply differential co-expression network analysis to these populations to resolve the cofactor biology behind the divergence.

### Key Results

| Metric | Value |
|--------|-------|
| Total cells | 155,650 (129,828 tumor; 25,822 normal-adjacent) |
| Samples / patients | 44 / 14 |
| Cell populations resolved | 12 |
| Basal epithelial cells | 52,126 (33.5%) |
| A3A in basal | mean 1.62, 25.7% positive |
| A3B in basal | mean 1.69, 33.2% positive |
| A3A (low-level, off-target) | macrophages (mean 0.62, 10.7%), myeloid DC (mean 0.40, 7.3%) |
| Basal cells with somatic calls + signature refitting | 31,912 |
| Basal cells with SBS2 > 0 | 5,911 (18.5%) |
| Spearman A3A vs SBS2 (all basal) | rho = 0.149, p = 3.06e-158 |
| Spearman A3B vs SBS2 (all basal) | rho = 0.052, p = 1.86e-20 |
| Spearman A3B vs CNV (all basal) | rho = 0.147, p = 9.85e-154 |
| Spearman A3A vs CNV (all basal) | rho = -0.199, p = 3.35e-283 |
| SBS2+ basal cells with no detectable A3A/A3B | 1,899 (32.1%) |
| Divergence axis, A3A-dominant tumor | median +0.60 (n = 11,252) |
| Divergence axis, A3B-dominant tumor | median -0.77 (n = 12,438) |
| Divergence axis, normal-adjacent | overlapping near center (A3A-dom n = 119, A3B-dom n = 57) |
| A3A fraction vs CNV (tumor vs normal) | rho = -0.44 vs +0.07 |
| A3A fraction vs SBS2 (tumor vs normal) | rho = +0.05 vs +0.04 (weak in both) |

| Method | Detail |
|--------|--------|
| Annotation | popV (Tabula Sapiens) + cluster-level refinement |
| Cancer detection | CytoTRACE2 + inferCNV dual-model consensus |
| Somatic mutation calling | SComatic (cell-type-pooled, germline-filtered) |
| Signature deconvolution | Semi-supervised NNLS against COSMIC v3.4 |
| Core signatures (always retained) | SBS2, SBS13, SBS5 |
