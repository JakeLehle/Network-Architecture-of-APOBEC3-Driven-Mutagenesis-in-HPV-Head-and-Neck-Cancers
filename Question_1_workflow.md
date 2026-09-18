# Question 1 Workflow

**Maps to manuscript Results 4.1** — *Single-cell analysis localizes A3-driven mutagenesis to basal epithelial cells and resolves divergent A3A:SBS2 and A3B:CNA linked programs.*

> **Terminology note.** The manuscript writes copy-number alteration as **CNA**. The code predates that decision and uses `CNV` throughout (`cnv_score`, `CNV_HIGH`, `NETWORK_CNV_VS_NORMAL`). The two refer to the same quantity. Column and directory names below are given as they appear on disk.

## Question 1: At single-cell resolution, where does A3-driven mutagenesis localize, and do A3A and A3B drive divergent mutational programs?

### Rationale

A3A and A3B are necessary but not sufficient for SBS2, and bulk tissue cannot say why because it averages over two layers of heterogeneity at once. First, A3A and A3B are not confined to tumor cells. A3A is expressed in macrophages and dendritic cells and A3B in plasma cells and some B cells, so a signature measured in bulk cannot be cleanly attributed to the cells that produced it. Second, even a pure population of tumor epithelial cells is not uniform, because the infecting virus sits at a different phase of its lifecycle from one cell to the next. As a result, no tumor presents an uncontaminated SBS2 signal in bulk, and the comparison against normal tissue cannot be made there.

This question moves to single-cell resolution to resolve both layers. Using the HPV16-positive HNSCC dataset GSE173468, with tumor and matched normal-adjacent tissue, we call somatic mutations and assign SBS2 to each cell, then ask two things: whether A3-driven mutagenesis localizes to a specific epithelial compartment, and whether A3A and A3B drive divergent classes of genomic damage within that compartment, point mutations (SBS2) versus copy-number change (CNA). This section establishes the localization and the divergence that the rest of the paper builds on.

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
│   ├── setup.py, pyproject.toml, MANIFEST.in, LICENSE
│   ├── cli/                                 #   cli, sample_information, create_config, run_config
│   └── snakemake_wrapper/
│       ├── Snakefile
│       ├── config.yaml
│       ├── envs/                            #   sc_pre, sc_post, scomatic, kraken2, signatures
│       ├── rules/                           #   scomatic_signatures.smk
│       └── scripts/                         #   cellranger_count, scanpy_qc_annotation,
│                                            #   cancer_cell_detection, kraken2_viral_detection,
│                                            #   scomatic_mutation_calling, SingleCellGenotype,
│                                            #   signature_analysis, viral_integration,
│                                            #   summarize_viral_detection, generate_summary
├── Run_Cluster_Catcher_Pipeline.sh          # SLURM runner: install, configure, and execute ClusterCatcher
├── Generate_Figure3_A3_Dominance_Slider.py  # Figure 3: A3A/A3B dominance divergence slider
├── Diagnostic_SBS13_Presence_and_REV1.py    # Results 4.1: SBS13 quantification and the REV1/UNG bypass test
├── Step02_Supplemental_Marker_Validation.py # Supplemental Figure 3: classical-marker validation of popV annotation
├── Step01_Generate_Figure3_Panels.py        # Panel helper; feeds the Fig 1/2 assembly (see naming note below)
└── TROUBLESHOOTING/                          # Diagnostics and development scripts (not documented here)
```

Each Snakemake rule runs in its own conda environment (`envs/`), so the five environment files are the authoritative record of what each module was executed with. `SingleCellGenotype.py` is the vendored SComatic helper that writes the per-cell genotype table used later by the neoantigen analysis.

> **Naming note.** `Step01_Generate_Figure3_Panels.py` is named for an earlier figure order. Its UMAP and COSMIC panels now feed **Figures 1 and 2**, not Figure 3. The file name was left unchanged so that existing output paths and logs stay valid.

> **Directory note.** `data/FIG_4/` holds the inputs and outputs for manuscript **Figure 5** (the network figure), and `data/FIG_3/` holds Figure 3. The `FIG_N` directory numbers were fixed before the figure order settled and do not track the manuscript numbering. `adata_final.h5ad` and `signature_weights_per_cell.txt` live under `data/FIG_4/00_input/` and are read by scripts across every question.

### External Dependencies

| Software | Version | Purpose | Installation |
|----------|---------|---------|--------------|
| **ClusterCatcher** | >=1.3.0 | End-to-end scRNA-seq pipeline | `git clone https://github.com/Diako-Lab/ClusterCatcher.git && pip install -e .` |
| **SRAscraper** | latest | Automated SRA data download | `git clone https://github.com/Diako-Lab/SRAscraper.git && pip install -e .` |
| **Cell Ranger** | >=7.0 | FASTQ alignment and counting | [10x Genomics](https://support.10xgenomics.com/) |
| **SComatic** | latest | Somatic mutation calling | `git clone https://github.com/cortes-ciriano-lab/SComatic.git` |
| **CytoTRACE2** | latest | Stemness/potency scoring | `git clone https://github.com/digitalcytometry/cytotrace2.git` |
| **Snakemake** | >=7.0 | Pipeline orchestration | `conda install -c bioconda snakemake` |

Conda environments: ClusterCatcher runs under its own `environment.yml`. The bespoke scripts in this directory run under `sc_pre` (scanpy, anndata, Python 3.11) or `NETWORK`; either works, since they only read the AnnData object and the weights table.

---

### Pipeline Overview

This question is engine-plus-three-analyses. ClusterCatcher does all of the heavy data processing and produces a single fully annotated AnnData object, `adata_final.h5ad`, that carries cell-type annotations, cancer/normal status, per-cell somatic mutations, and per-cell signature weights. Figures 1 and 2 and the early supplementals (cell-type annotation, A3 family localization, the mutation, SBS2, A3A, A3B and CNV UMAPs, and the aggregated COSMIC comparison) are assembled from these ClusterCatcher outputs separately, so there is no single repo script that "makes Figure 1" or "makes Figure 2."

Three bespoke analyses live in this directory as standalone main-dir scripts:

- **`Generate_Figure3_A3_Dominance_Slider.py`** produces Figure 3, the within-tissue SBS2-versus-CNA divergence axis split by A3A versus A3B dominance.
- **`Diagnostic_SBS13_Presence_and_REV1.py`** supports the SBS13 sentences in Results 4.1, quantifying SBS13 against SBS2 and testing whether the deficit is biological, technical, or a limit of the deconvolution.
- **`Step02_Supplemental_Marker_Validation.py`** produces Supplemental Figure 3, the classical-marker validation of the popV annotations.

`Step01_Generate_Figure3_Panels.py` remains in the directory as a panel helper. Its UMAP and COSMIC panels feed the ClusterCatcher-derived Figure 1 and 2 assembly rather than standing as a defined figure step, so it is not written up as one below.

---

### ClusterCatcher: the engine

ClusterCatcher (v1.3.0) runs as a seven-module Snakemake pipeline launched by `Run_Cluster_Catcher_Pipeline.sh`. Full parameter documentation is in the [ClusterCatcher README](https://github.com/Diako-Lab/ClusterCatcher); the condensed view of what each module contributes to this question is below.

| Module | Script | Produces |
|--------|--------|----------|
| 1. Alignment | `cellranger_count.py` | GRCh38 alignment, filtered feature-barcode matrices, CB-tagged BAMs |
| 2-3. QC + annotation | `scanpy_qc_annotation.py` | Adaptive MAD QC, Scrublet doublet removal, popV annotation with cluster-level refinement; annotation and composition UMAPs that feed Figure 1 and Supplemental Figures 1-2 |
| 4. Cancer detection | `cancer_cell_detection.py` | CytoTRACE2 + inferCNVpy dual-model consensus; `cnv_score` and stemness per cell, feeding **Supplemental Figure 4** |
| 5. Viral detection | `kraken2_viral_detection.py` | Per-cell organism counts including HPV16 (`raw_HPV16`; used downstream in Question 4) |
| 6. Mutation calling | `scomatic_mutation_calling.py` | Cell-type-pooled, germline-filtered somatic variants, per-cell callable sites |
| 7. Signature deconvolution | `signature_analysis.py` | 96-context matrix, semi-supervised NMF refit against COSMIC v3.4, per-cell signature weights; writes `adata_final.h5ad` |

Key processing choices that matter for interpretation:

- Somatic calls are pooled within each annotated cell type and normalized to each cell's callable sites (depth >= 5).
- Signature refitting always retains SBS2, SBS13, and SBS5, and adds other HNSCC-relevant COSMIC signatures by scree-plot elbow on reconstruction error. **Fifteen signatures were retained in the final fit.**
- `signature_weights_per_cell.txt` is written **signatures x cells** and requires `.T` on load.
- The basal annotation in `adata.obs` is lowercase `basal cell`.
- `inferCNVpy` is run with window 250 and normal cells as reference.

**Engine outputs used in this question:** `signatures/adata_final.h5ad` and `signatures/signature_weights_per_cell.txt`, staged to `data/FIG_4/00_input/`.

---

### Figure 3: A3 dominance divergence slider

**Script:** `Generate_Figure3_A3_Dominance_Slider.py`

This is the one bespoke analysis behind a main figure in this section. It tests whether the A3A-to-SBS2 and A3B-to-CNA link, which is visible in the Figure 2 UMAPs but weak at the level of global expression, sharpens when cells are split by which enzyme dominates.

**Approach:**
- Restrict to basal cells expressing either enzyme (A3A + A3B > 0).
- Define an A3A dominance fraction, A3A / (A3A + A3B); cells above 0.5 are A3A-dominant, below 0.5 A3B-dominant. **Cells sitting at exactly 0.5 fall into neither arm and are excluded** (1,014 of 24,880 expressing basal cells). The four reported group sizes therefore sum to 23,866 rather than to the expressing total.
- Place each cell on a within-tissue axis, z(SBS2) minus z(CNV), where SBS2 and the inferCNVpy `cnv_score` are each standardized within that tissue's expressing-either population. Standardizing within tissue keeps "center" at each tissue's own median so the tumor/normal marginal differences do not push normal off-center artefactually.
- Plot four rows: A3A-dominant tumor, A3B-dominant tumor, A3A-dominant normal-adjacent, A3B-dominant normal-adjacent. The right pole is SBS2-high, the left pole is CNA-high (productive).

**Design notes carried in the script header (worth keeping in mind for the text):**
- Stemness (CytoTRACE2) is deliberately not blended into the productive pole, because the A3A fraction anti-correlates with stemness in both tumor and normal, so stemness is not tumor-specific. The productive pole is CNA alone.
- The per-enzyme A3 > 0 conditioning is gone; conditioning on expressing cells deleted the co-occurrence signal. Dominance is defined among cells expressing either enzyme.
- Absolute A3 level is not the axis driver. SBS2 does not separate the two enzymes at the level of absolute expression; dominance does.

**Input:** `data/FIG_4/00_input/adata_final.h5ad`, `data/FIG_4/00_input/signature_weights_per_cell.txt`. `TARGET_CELL_TYPE = "basal cell"`.

**Output (-> `data/FIG_3/figures/`):**
- `Figure3_A3_dominance_slider.pdf` / `.png` (300 DPI). SBS2 pole `#ed6a5a`, CNA pole `#F6D155`.

**Honest caveat for the text:** the separation is carried mainly by CNA, low CNA on the A3A-dominant side, rather than by high SBS2, since SBS2 is shared between the two enzymes. This fits the maintenance (low-CNA) versus productive (high-CNA) framing that Question 4 develops.

---

### Results 4.1: SBS13 quantification and the REV1 test

**Script:** `Diagnostic_SBS13_Presence_and_REV1.py`

SBS2 and SBS13 are both APOBEC3-associated and arise from the same lesion: A3 deaminates C to U, and replication over the uracil yields C>T (SBS2), whereas UNG excision followed by REV1/POL-zeta bypass of the abasic site yields C>G (SBS13). This script establishes why the paper reports SBS2 and not SBS13, and rules out the obvious mechanistic explanation.

**Approach:** rank all fifteen fitted signatures by mean weight in basal cells; test SBS2/SBS13 co-occurrence per cell; compare SBS13 across the three populations; correlate both signatures against A3A and A3B; and profile REV1, REV3L, MAD2L2, UNG, SMUG1, TDG and APEX1 against POLH and POLK as similarly-expressed controls.

**Findings:**

| Observation | Value |
|---|---|
| SBS13 rank among 15 fitted signatures | **15 of 15** by mean weight |
| SBS13 in basal cells | mean weight 0.040, 7.7% of cells positive |
| SBS2 in basal cells, for comparison | mean weight 0.186, 18.5% positive |
| Pooled SBS13:SBS2 ratio | 0.215 |
| SBS13 vs A3A expression | rho = 0.006, **p = 0.32** (not significant) |
| SBS2 vs A3A expression, for comparison | rho = 0.149, p = 3.06e-158 |
| SBS13 across populations | SBS2-HIGH 0.073, CNV-HIGH 0.056, NORMAL 0.051 |
| REV1 in basal cells | mean 1.045, 22.0% positive, **highest of all 12 cell types** |
| REV1 across populations | SBS2-HIGH 1.226, CNV-HIGH 1.232, NORMAL 1.418 |

**The REV1 hypothesis is not supported.** REV1 is basal-enriched rather than depleted, exceeds both POLH (0.507) and POLK (0.503), and shows no gradient across the three populations. The whole bypass arm (UNG 1.43, TDG 1.65, APEX1 2.60, MAD2L2 1.41) is likewise present. Limited translesion capacity does not explain the low SBS13 weight.

**What the manuscript claims, and what it deliberately does not.** Results 4.1 states that SBS13 was detected, carried the lowest mean weight of the fifteen signatures retained, was not enriched in any basal population, and showed no association with A3A. It makes **no** claim that SBS13 is biologically absent or low in HPV-positive HNSCC, which would be wrong. Three observations argue that the per-cell SBS13 weight is not a reliable readout at this mutation depth, and they belong in a reviewer response rather than the paper:

1. SBS13 does not track A3A at all (p = 0.32), where SBS2 does strongly, on the same cells.
2. Among carriers, SBS2 and SBS13 are strongly anti-correlated (rho = -0.586), the fingerprint of NMF splitting collinear signatures that share the TCW context.
3. CNV-HIGH cells carry zero SBS2 by construction yet the second-highest SBS13 mean.

**Deciding test, not yet run.** STEP 5 of the script extracts the raw C>T and C>G call counts straight from the SComatic table, independent of the NMF. If the raw C>G:C>T ratio is normal while the fitted SBS13:SBS2 ratio is not, the deficit is a deconvolution limit. That step currently skips: `comment="#"` consumes the header row so the REF/ALT columns fail to resolve. Fix by counting `##` lines and passing `skiprows` instead.

**Input:** `adata_final.h5ad`, `signature_weights_per_cell.txt`, `three_group_assignments.tsv`, SComatic filtered genotype table.

**Output (-> `data/FIG_2/DIAGNOSTIC_SBS13/`):** `SBS13_vs_SBS2_UMAP`, `SBS13_vs_SBS2_weight_distribution`, `REV1_TLS_expression_basal` (PDF + PNG, 300 DPI), plus `signature_summary_basal.tsv`, `sbs13_by_celltype.tsv`, `tls_gene_expression.tsv`.

---

### Supplemental Figure 3: classical-marker annotation validation

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
| Fig 2 | A3-association multipanel: normalized mutations, SBS2, A3A, A3B UMAPs (a-c); SBS2 and A3 associations (d-e); aggregated 96-context versus COSMIC SBS2 (f); CNA and stemness (g-h) | Assembled from ClusterCatcher outputs (mutations, weights, expression, inferCNVpy) |
| Fig 3 | A3A/A3B dominance divergence slider (within-tissue SBS2-vs-CNA axis) | `Generate_Figure3_A3_Dominance_Slider.py` |
| Supp Fig 1 | Per-cell popV prediction UMAP | ClusterCatcher annotation |
| Supp Fig 2 | A3 family expression across cell types | Assembled from ClusterCatcher outputs |
| Supp Fig 3 | Classical-marker validation UMAPs (12 cell types x 2) | `Step02_Supplemental_Marker_Validation.py` |
| **Supp Fig 4** | CytoTRACE2 and inferCNVpy per cell | ClusterCatcher cancer detection |

> **Corrected from the previous version of this walkthrough**, which listed the CytoTRACE2/inferCNV panel as Supplemental Figure 5. Manuscript Supplemental Figure 5 belongs to Question 2 (the differential co-expression network supplement). Results 4.1 cites Supp. Figs. 1-3 for the annotation and Supp. Fig. 4 for CytoTRACE2/inferCNVpy.

**Narrative arc:**
1. Twelve cell populations resolve from the HPV16-positive HNSCC data; cluster-level popV annotation is validated by classical markers (Fig 1, Supplemental Figures 1-3).
2. A3A and A3B expression concentrates almost exclusively in basal epithelial cells, with only low-level A3A in macrophages and myeloid dendritic cells (Fig 1, Fig 2a-c).
3. Among basal cells, SBS2 burden associates with A3A more than A3B, but the global association is weak, and the aggregated trinucleotide profile of the highest-SBS2 cells matches the COSMIC SBS2 reference (Fig 2d-f).
4. Most SBS2-positive basal cells carry no detectable A3A (56.7%) or A3B (49.8%) at capture, and a third carry neither (32.1%), consistent with pulsatile A3 induction or transcript dropout.
5. SBS13, the second A3-associated signature, is detected but is the weakest of the fifteen fitted signatures, is not enriched in any basal population, and shows no association with A3A, so SBS2 is used as the readout throughout.
6. Elevated CNA and stemness map instead to an A3B-associated basal subpopulation, so A3A-linked point mutations and A3B-linked chromosomal instability occupy separable subregions of the same compartment (Fig 2g-h).
7. Splitting basal cells by A3 dominance and placing them on a within-tissue SBS2-vs-CNA axis resolves the divergence that global expression obscures: A3A-dominant cells shift toward SBS2, A3B-dominant cells toward CNA, and only in tumor (Fig 3).
8. **Motivates Question 2:** apply differential co-expression network analysis to these populations to resolve the cofactor biology behind the divergence.

### Key Results

All values below were re-verified against on-disk pipeline output in September 2026 and match the submitted manuscript.

| Metric | Value |
|--------|-------|
| Total cells | 155,650 (129,828 tumor; 25,822 normal-adjacent) |
| Samples / patients | 44 / 14 |
| Cell populations resolved | 12 |
| Basal epithelial cells | 52,126 (33.5%) |
| Basal cells by source | 51,572 tumor / 554 normal-adjacent |
| A3A in basal | mean 1.62, 25.7% positive |
| A3B in basal | mean 1.69, 33.2% positive |
| A3A (low-level, off-target) | macrophages (mean 0.62, 10.7%), myeloid DC (mean 0.40, 7.3%) |
| Basal cells with somatic calls + signature refitting | 31,912 |
| Basal cells with SBS2 > 0 | 5,911 (18.5%) |
| Signatures retained in the final fit | 15 |
| Spearman A3A vs SBS2 (all basal) | rho = 0.149, p = 3.06e-158 |
| Spearman A3B vs SBS2 (all basal) | rho = 0.052, p = 1.86e-20 |
| Spearman A3B vs CNA (all basal) | rho = 0.147, p = 9.85e-154 |
| Spearman A3A vs CNA (all basal) | rho = -0.199, p = 3.35e-283 |
| SBS2+ basal cells with no detectable A3A | 3,351 (56.7%) |
| SBS2+ basal cells with no detectable A3B | 2,943 (49.8%) |
| SBS2+ basal cells with neither enzyme detected | 1,899 (32.1%) |
| SBS13 mean weight / prevalence in basal | 0.040 / 7.7% (rank 15 of 15) |
| SBS13 vs A3A | rho = 0.006, p = 0.32 (ns) |
| A3-expressing basal cells (A3A + A3B > 0) | 24,880 |
| Excluded at exactly 0.5 dominance | 1,014 (leaving 23,866) |
| Divergence axis, A3A-dominant tumor | median +0.60 (n = 11,252) |
| Divergence axis, A3B-dominant tumor | median -0.77 (n = 12,438) |
| Divergence axis, normal-adjacent | overlapping near center (A3A-dom n = 119, A3B-dom n = 57) |
| A3A fraction vs CNA (tumor vs normal) | rho = -0.44 vs +0.07 |
| A3A fraction vs SBS2 (tumor vs normal) | rho = +0.05 vs +0.04 (weak in both) |

| Method | Detail |
|--------|--------|
| Annotation | popV (Tabula Sapiens) + cluster-level refinement |
| Cancer detection | CytoTRACE2 + inferCNVpy dual-model consensus (window 250, normal reference) |
| Somatic mutation calling | SComatic (cell-type-pooled, germline-filtered, depth >= 5) |
| Signature deconvolution | Semi-supervised NMF against COSMIC v3.4 |
| Core signatures (always retained) | SBS2, SBS13, SBS5 |
| Signature weights file | `signature_weights_per_cell.txt`, signatures x cells, requires `.T` on load |

---

### Retired Figure 3 variants

Three alternate versions of the Figure 3 analysis remain under `TROUBLESHOOTING/` and are **not** the published figure. They are listed here only so that nobody mistakes one for the canonical script:

| File | What it was |
|---|---|
| `Generate_Figure3_A3_Lean_Slider.py` | earlier slider using a lean score rather than the z(SBS2) − z(CNV) axis |
| `Generate_Figure3_A3_SBS2_CNV_Scatter.py` | scatter form of the same contrast, replaced by the slider |
| `Diagnostic_Fig3_Lean_Decompose.py` | decomposition of the retired lean score |

`Generate_Figure3_A3_Dominance_Slider.py` in the main directory is the only script that produces manuscript Figure 3.

Two further diagnostics in `TROUBLESHOOTING/` are **candidate replacements pending review**, not retired: `Diagnostic_Fig3_A3_Plane_SBS2_CNV.py` (2×2 A3A/A3B dot plane) and `Diagnostic_Fig3_Dominance_Boxplots_SBS2_CNV.py` (dominance-split boxplots). Either could be promoted to Figure 3; both reuse the slider's plumbing and colors.

### Numbers verified outside the main directory

The September 2026 audit reproduced every Results 4.1 number from two scripts that live in `TROUBLESHOOTING/` and are therefore not documented above. Recorded here so the provenance is not lost:

- `Diagnostic_section3_numbers.py` — cell counts, per-cell-type A3 means and prevalence, the four Spearman correlations, the pulsatile percentages, CytoTRACE2 and CNA summaries.
- `Diagnostic_Fig3_Dominance_Statistics.py` — the four group sizes, the four A3A-fraction correlations, plus patient-level paired tests not used in the manuscript.

