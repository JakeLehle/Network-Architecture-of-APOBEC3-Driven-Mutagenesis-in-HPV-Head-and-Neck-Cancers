## Question 2: Can differential co-expression network analysis identify gene communities associated with elevated SBS2 mutagenesis in A3-matched HNSCC tumors?

### Rationale

Figure 1 established that A3 expression is necessary but not sufficient for SBS2 mutagenesis in HNSCC, and that the highest mutational burden requires the coordinated activity of both A3B and A3A alongside as-yet-unidentified cofactors. This question exploits a key observation from Figure 1: tumors with similar A3A+A3B expression have vastly different SBS2 levels. By constructing Spearman co-expression networks separately for high-SBS2 and low-SBS2 tumors matched for A3 expression, then computing the differential network (HIGH − LOW), we can identify gene communities whose co-expression relationships are specifically rewired in the high-mutagenesis condition. These communities are candidate cofactor networks that regulate A3 enzymatic access to genomic DNA.

The pipeline was originally developed by Dr. Mohadeseh Soleimanpour (Texas Biomedical Research Institute) as a monolithic Python script. For this paper, the pipeline was refactored into modular steps (Step01–Step08), with adjustments to gene feature selection, community detection parameters, and figure generation to align with the manuscript narrative.

### Data Sources

| Source | Description |
|--------|-------------|
| `TCGA_master_FPKM_UQ.tsv` | Pan-cancer FPKM-UQ expression matrix generated in Question 1 (Step 3) |
| `Mutation_Table_Tumors_TCGA.tsv` | COSMIC SBS signature weights (same external input as Question 1) |

### Directory Structure

```
scripts/NETWORK/
├── RUN_NETWORK_PIPELINE.sh              # SLURM batch script (Steps 01–07)
├── network_config.py                    # Centralized pipeline parameters
│
├── Step01_Load_Clean_TCGA.py            # Step 1 — Load & clean expression data
├── Step02_Merge_SBS_Signatures.py       # Step 2 — Merge expression with SBS signatures
├── Step03_Differential_Expression.py    # Step 3 — Group definition & differential expression
├── Step04_Correlation_Networks.py       # Step 4 — Build TOP/BOTTOM/DIFF correlation networks
├── Step04.1_Sweep_DIFF_Threshold.py     # Step 4.1 — Diagnostic: sweep DIFF thresholds
├── Step05_Community_Detection.py        # Step 5 — Leiden community detection
├── Step06_Centrality_Metrics.py         # Step 6 — Per-gene centrality metrics
├── Step07_Generate_Figure2_Panels.py    # Step 7 — Generate Figure 2 panels (a,b,c)
├── Step08_Pipeline_Summary.py           # Step 8 — Pipeline summary & KEGG enrichment
│
└── TROUBLESHOOTING/
    ├── network_community_analysis_pipeline.py  # Original monolithic pipeline (Dr. Soleimanpour)
    └── network_utils.py                        # Original utility functions
```

### Configuration

All pipeline parameters are centralized in `network_config.py`:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `A3_SUM_PERCENTILE` | 0.50 | Keep tumors with A3A+A3B above median |
| `SBS2_HIGH_PERCENTILE` | 0.75 | Top 25% SBS2 within high-A3 = SBS2-HIGH group |
| `SBS2_LOW_PERCENTILE` | 0.25 | Bottom 25% SBS2 within high-A3 = SBS2-LOW group |
| `MIN_SAMPLES_DETECTED` | 20 | Gene must be detected in ≥20 samples |
| `RAW_P_THRESHOLD` | 0.05 | Wilcoxon p-value cutoff for gene selection |
| `LOGFC_THRESHOLD` | 0 | No fold-change filter (bulk effect sizes too small) |
| `FORCE_KEEP_A3` | True | Always retain A3 genes in the network |
| `CORRELATION_METHOD` | spearman | Rank-based correlation |
| `CORR_THRESHOLD` | 0.80 | \|rho\| threshold for TOP/BOTTOM network edges |
| `DIFF_THRESHOLD` | 0.70 | \|delta-rho\| threshold for DIFF network edges |
| `COMMUNITY_METHOD` | leiden | Leiden algorithm (RBConfigurationVertexPartition) |
| `COMMUNITY_RESOLUTIONS` | [0.2, 0.4, 0.6, 0.8, 1.0] | Resolution sweep range |
| `RUNS_PER_RESOLUTION` | 15 | Stability assessment runs per resolution |
| `TARGET_BIG_COMMUNITIES` | 14 | Max communities before merging |
| `MIN_COMMUNITY_SIZE` | 10 | Communities smaller than this are merged |

### Pipeline Overview

The analysis is executed as a single SLURM batch job:

```
 SLURM Job: RUN_NETWORK_PIPELINE.sh
 ───────────────────────────────────
   TCGA_master_FPKM_UQ.tsv  ← from Question 1
   Mutation_Table_Tumors_TCGA.tsv  ← external input
       │
       ▼
   Step 1: Load & clean TCGA expression
       │
       ▼
   Step 2: Merge with SBS mutation signatures
       │
       ▼
   Step 3: Define A3-matched groups + differential expression
       │
       ▼
   Step 4: Build TOP/BOTTOM/DIFF correlation networks
       │    └── Step 4.1: Diagnostic threshold sweep (optional)
       ▼
   Step 5: Leiden community detection on DIFF network
       │
       ▼
   Step 6: Centrality metrics (degree, betweenness, closeness, eigenvector)
       │
       ▼
   Step 7: Generate Figure 2 panels
       │
       ▼
   Step 8: Pipeline summary + KEGG enrichment
       │
       ▼
   FIGURE 2
```

### Environment

All scripts are executed using the `NETWORK` conda environment (see `environments/Network.yml` or `scripts/NETWORK/enviornment_NETWORK.yml`).

---

### SLURM Job: Network Analysis Pipeline

**Batch script:** `RUN_NETWORK_PIPELINE.sh`

Submits Steps 1–8 as a serial pipeline. Requires 500 GB memory and 32 CPUs for the correlation matrix computations.

```bash
conda run -n NETWORK python Step01_Load_Clean_TCGA.py
conda run -n NETWORK python Step02_Merge_SBS_Signatures.py
conda run -n NETWORK python Step03_Differential_Expression.py
conda run -n NETWORK python Step04_Correlation_Networks.py
conda run -n NETWORK python Step04.1_Sweep_DIFF_Threshold.py
conda run -n NETWORK python Step05_Community_Detection.py
conda run -n NETWORK python Step06_Centrality_Metrics.py
conda run -n NETWORK python Step07_Generate_Figure2_Panels.py
conda run -n NETWORK python Step08_Pipeline_Summary.py
```

**Environment:** `NETWORK` conda environment, `--mem=500G`, `--cpus-per-task=32`

---

#### Step 1: Load and Clean TCGA Expression Data

**Script:** `Step01_Load_Clean_TCGA.py`

Loads the pan-cancer FPKM-UQ master table generated in Question 1, parses the three-row header structure (column headers, gene symbols, gene biotypes), removes pseudoautosomal Y-linked (_PAR_Y) gene columns, cleans versioned ENSG identifiers (removing `.15` suffixes), resolves duplicated gene columns, and builds an ENSG-to-symbol mapping dictionary.

**Input:**
- `TCGA_master_FPKM_UQ.tsv` from Question 1 Step 3

**Output (→ `data/FIG_2/01_cleaned_expression/`):**
- Cleaned expression matrix (pickle)
- `ensg_to_symbol.json` — gene ID to symbol mapping used throughout the pipeline

**Key details:**
- Filters for tumor samples only (barcode positions 14–15: `01`)
- Clinical metadata columns (Project_ID, Tissue_Type, Case_ID, File_ID, Entity_ID) preserved separately
- All downstream analysis uses cleaned ENSG IDs (no version suffix)

---

#### Step 2: Merge with SBS Mutation Signatures

**Script:** `Step02_Merge_SBS_Signatures.py`

Merges the cleaned TCGA expression data with COSMIC SBS mutational signature weights using Entity_ID matching. Removes duplicated samples that arise from many-to-many barcode matches.

**Input:**
- Cleaned expression data from Step 1
- `Mutation_Table_Tumors_TCGA.tsv` — COSMIC SBS signature weights (external input)

**Output (→ `data/FIG_2/02_merged_data/`):**
- Merged expression + SBS signature dataframe (pickle)

**Key details:**
- Entity_IDs shortened to 4-part TCGA barcodes for matching (`TCGA-XX-YYYY-ZZZZ`)
- Duplicated barcode keys are removed entirely (both rows) to ensure one-to-one mapping
- Filters for specified cancer type(s) in `CANCER_TYPES` (default: `["TCGA-HNSC"]`)

---

#### Step 3: Group Definition and Differential Expression

**Script:** `Step03_Differential_Expression.py`

Defines A3-matched SBS2-HIGH and SBS2-LOW tumor groups, then performs differential expression analysis between them. This step is the single source of truth for group assignments used by all downstream steps.

**Input:**
- Merged data from Step 2

**Output (→ `data/FIG_2/03_differential_expression/{cancer_type}/`):**
- `{ct}_SBS2_HIGH_group.pkl` / `{ct}_SBS2_LOW_group.pkl` — group dataframes
- `{ct}_diffexpr_stats.csv` — per-gene Wilcoxon test results
- `{ct}_selected_genes.csv` — genes passing significance threshold
- `{ct}_selected_genes_filtered.csv` — genes passing network-readiness filter

**Key details:**
- A3 matching: sum A3A+A3B FPKM-UQ per tumor, keep top 50% (above median)
- Within high-A3 tumors: SBS2-HIGH = top quartile, SBS2-LOW = bottom quartile
- Result: 53 HIGH and 53 LOW tumors with comparable A3A/A3B expression (A3A: 30.5 vs 22.5; A3B: 9.5 vs 8.6 FPKM-UQ)
- Differential test: Wilcoxon rank-sum on log1p(FPKM-UQ)
- Gene selection: raw p < 0.05 → 2,810 genes selected (from 19,101 tested)
- A3 family genes force-retained regardless of significance
- Network-readiness filter ensures genes are present and non-constant in both groups

---

#### Step 4: Correlation Network Construction

**Script:** `Step04_Correlation_Networks.py`

Builds Spearman co-expression networks for the SBS2-HIGH (TOP) and SBS2-LOW (BOTTOM) groups independently, then computes the differential network.

**Input:**
- Group dataframes and selected gene list from Step 3

**Output (→ `data/FIG_2/04_correlation_networks/{cancer_type}/`):**
- `corr_matrices/` — TOP, BOTTOM, DIFF correlation matrices (pickle + CSV)
- `heatmaps/` — correlation distribution plots
- `edge_lists/` — network edge lists for Cytoscape/Gephi
- `graph_objects/` — NetworkX graph pickles (G_top, G_bot, G_diff, with `_noiso` variants)

**Key details:**
- Spearman correlation computed across the 2,810 selected genes for each group
- TOP and BOTTOM networks: edges retained where |rho| ≥ 0.80
- DIFF network: element-wise subtraction (TOP − BOTTOM), edges retained where |delta-rho| ≥ 0.70
- The threshold captures gene pairs whose co-expression changed substantially between conditions
- Edge attributes: `weight` (signed delta-rho) and `abs_weight` (unsigned magnitude)
- Isolated nodes (degree = 0) removed → DIFF network: 1,530 nodes, 3,682 edges, avg degree 4.8

**Step 4.1 — Diagnostic threshold sweep (optional):**

**Script:** `Step04.1_Sweep_DIFF_Threshold.py`

Sweeps DIFF thresholds from 0.3 to 0.9, reporting edge count, node count, average degree, density, and community structure at each. Used to select the 0.70 threshold as the optimal balance between network connectivity and community interpretability.

---

#### Step 5: Community Detection

**Script:** `Step05_Community_Detection.py`

Runs Leiden community detection on the largest connected component of the DIFF network across multiple resolutions, evaluates partition stability, and merges small communities.

**Input:**
- DIFF network graph from Step 4

**Output (→ `data/FIG_2/05_communities/{cancer_type}/`):**
- `{ct}_resolution_sweep.csv` — per-resolution modularity, ARI, NMI metrics
- `{ct}_best_partition.csv` — final gene-to-community assignments
- `{ct}_community_gene_lists.csv` — genes per community
- `{ct}_G_comm.gpickle` — community graph object

**Key details:**
- Algorithm: Leiden (RBConfigurationVertexPartition) with `abs_weight` as edge weight
- Leiden groups genes by the **magnitude** of co-expression change, regardless of direction
- Resolution sweep: 0.2–1.0, with 15 runs per resolution (seed = 42) for stability
- Best resolution: 1.0 (modularity = 0.53, ARI = 0.46)
- Communities with < 10 genes merged
- Final result: **13 communities, 1,437 genes**
- Community sizes range from 19 (Community 18) to 288 (Community 0)

---

#### Step 6: Centrality Metrics

**Script:** `Step06_Centrality_Metrics.py`

Computes per-gene network centrality metrics for the TOP, BOTTOM, and DIFF networks. Identifies hub genes whose network role changes between conditions.

**Input:**
- Graph objects from Step 4
- Community assignments from Step 5

**Output (→ `data/FIG_2/06_centrality_metrics/{cancer_type}/`):**
- `{ct}_{TOP|BOTTOM|DIFF}_metrics.csv` — per-gene centrality tables
- `{ct}_hub_genes_DIFF.csv` — top hub genes by multiple metrics

**Key details:**
- Metrics computed: degree, betweenness centrality, closeness centrality, eigenvector centrality, strength (sum of |edge weights|)
- Hub genes = genes ranking highly across multiple metrics in the DIFF network
- Used for labeling top genes in Figure 2c network plot

---

#### Step 7: Generate Figure 2 Panels

**Script:** `Step07_Generate_Figure2_Panels.py`

Generates publication-quality Figure 2 panels from pipeline outputs. Run after Steps 01–06 are complete.

**Input:**
- All outputs from Steps 01–06

**Output (→ `data/FIG_2/FIGURE_2_PANELS/{cancer_type}/`):**
- `{ct}_Panel_A_selection.pdf/.png` — Figure 2a
- `{ct}_Panel_B_dual_heatmap.pdf/.png` — Figure 2b
- `{ct}_Panel_C_network_full.pdf/.png` — Figure 2c
- `community_zooms/` — per-community zoomed network subplots

**Panel A — Tumor selection plot (Figure 2a):**
- Scatter of summed A3A+A3B expression vs SBS2 weight for all 425 HNSC tumors
- HIGH (n = 53, coral) and LOW (n = 53, cream) groups highlighted
- Colored background regions matching Figure 1a convention
- Data-driven broken x-axis for outlier accommodation

**Panel B — Dual correlation heatmap (Figure 2b):**
- Side-by-side heatmaps of the TOP (high-SBS2 group) and DIFF (TOP − BOTTOM) Spearman correlation matrices
- 1,437 genes ordered by hierarchical clustering within and between communities
- Within-community ordering: average-linkage clustering on 1 − |corr| distance
- Between-community ordering: average-linkage on mean inter-community |corr|
- Plasma colormap with white community boundary lines
- Community labels annotated on right margin

**Panel C — Exploded community network (Figure 2c):**
- Two-tier community-aware layout: super-graph positions community centers (`INTER_SCALE = 6.0`), local spring layout positions nodes within each community (`INTRA_SCALE = 1.0`)
- Intra-community edges colored by differential correlation sign (red = gained, blue = lost)
- Inter-community edges shown as dashed gray lines
- Top 3 hub genes per community labeled (font_size = 14)
- A3 genes highlighted in gold with enlarged labels (font_size = 20)
- Legend with community labels (font_size = 18)
- Per-community zoomed subplots saved separately

---

#### Step 8: Pipeline Summary and KEGG Enrichment

**Script:** `Step08_Pipeline_Summary.py`

Gathers all pipeline parameters, per-step statistics, and community gene lists. Runs KEGG pathway enrichment on each community. Outputs a comprehensive report for writing the results and methods sections. Run separately after reviewing pipeline output.

**Dependencies:** `gseapy` (KEGG enrichment via Enrichr API; requires internet access)

**Input:**
- All outputs from Steps 01–06
- `network_config.py` parameters

**Output (→ `data/FIG_2/08_pipeline_summary/`):**
- `{ct}_pipeline_report.txt` — human-readable full report
- `{ct}_pipeline_report.json` — machine-readable report
- `{ct}_KEGG_enrichment_summary.csv` — one row per community with top KEGG term
- `{ct}_community_XX_KEGG.csv` — full KEGG enrichment results per community

**Key details:**
- Reports group sizes, A3A/A3B/SBS2 medians, gene counts at each filtering stage
- Reports network node/edge counts for TOP/BOTTOM/DIFF graphs
- Reports community sizes, resolution, modularity
- Cross-references 24 cell-type marker genes (12 cell types × 2 markers) against community assignments
- Cross-references A3 family gene community assignments
- KEGG enrichment via gseapy.enrichr (KEGG_2021_Human, adj. p < 0.05 significance)

---

### Figure 2 Summary

**Title:** Differential co-expression networks in A3-matched HNSCC tumors reveal gene communities associated with SBS2 mutagenic activity.

| Panel | Content | Script | Step |
|-------|---------|--------|------|
| a | A3-matched tumor selection: SBS2-HIGH vs SBS2-LOW groups | `Step07_Generate_Figure2_Panels.py` | 7 |
| b | Dual heatmap (TOP + DIFF) with hierarchical clustering by community | `Step07_Generate_Figure2_Panels.py` | 7 |
| c | Exploded Leiden community network with hub gene labels | `Step07_Generate_Figure2_Panels.py` | 7 |

**Narrative arc:**
1. A3-matched groups isolate cofactor-driven differences from A3 expression variation (panel a)
2. 13 gene communities identified with robust modularity (panels b, c)
3. Community 2 contains A3B + KRT5 (basal cell marker) — cofactor network may be basal-cell-specific
4. Community 0 contains A3G with striking anti-correlation pattern and loss of correlation in high-SBS2 — potential immune evasion signature
5. Most communities don't map to canonical KEGG pathways → novel A3-associated regulatory modules
6. Only 2/24 cell-type markers found in any community → bulk data cannot resolve cell-type specificity
7. **Motivates Question 3:** single-cell resolution needed to assign communities to cell types and validate A3B's role in basal epithelial cells

### Key Results from Step 8 Summary

| Metric | Value |
|--------|-------|
| Total HNSC tumors | 425 |
| SBS2-HIGH group | 53 tumors |
| SBS2-LOW group | 53 tumors |
| A3A median (HIGH / LOW) | 30.5 / 22.5 FPKM-UQ |
| A3B median (HIGH / LOW) | 9.5 / 8.6 FPKM-UQ |
| SBS2 median (HIGH / LOW) | 74.0 / 0.0 |
| Genes tested | 19,101 |
| Genes selected (differential) | 2,810 |
| DIFF network nodes | 1,530 |
| DIFF network edges | 3,682 |
| Communities | 13 (1,437 genes) |
| Modularity | 0.53 |
| Cell-type markers in communities | 2 / 24 (KRT5 in C2, RGS5 in C4) |
| A3 genes in communities | A3B in C2, A3G in C0 |

### KEGG Enrichment Highlights

| Community | Size | Top KEGG Term | adj. p | Notable |
|-----------|------|---------------|--------|---------|
| 0 | 288 | Thyroid hormone synthesis | 0.35 | Not significant; contains A3G |
| 1 | 238 | Dilated cardiomyopathy | 0.019 | 7 significant terms; immune signaling (phagosome, TB) |
| 2 | 160 | Herpes simplex virus 1 | 0.13 | Not significant; contains A3B + KRT5 |
| 3 | 123 | Herpes simplex virus 1 | 0.14 | Not significant |
| 4 | 122 | Hypertrophic cardiomyopathy | 0.21 | Not significant; contains RGS5 |
| 5 | 118 | MAPK signaling pathway | 0.74 | Not significant |
| 6 | 113 | Herpes simplex virus 1 | 0.011 | ZNF protein cluster |
| 7 | 70 | Cardiac muscle contraction | 0.011 | 3 significant terms |
| 8 | 68 | Neurotrophin signaling | 0.52 | Not significant |
| 9 | 53 | Colorectal cancer | 0.11 | Not significant |
| 10 | 44 | Herpes simplex virus 1 | 7.2×10⁻⁵ | ZNF protein cluster (strongest enrichment) |
| 11 | 21 | Glycosaminoglycan biosynthesis | 0.027 | 1 significant term |
| 18 | 19 | Protein digestion/absorption | 0.14 | Not significant |

---

### Troubleshooting Scripts

The following scripts in `scripts/NETWORK/TROUBLESHOOTING/` are not part of the main pipeline but were used during development:

| Script | Purpose |
|--------|---------|
| `network_community_analysis_pipeline.py` | Original monolithic pipeline by Dr. Mohadeseh Soleimanpour; refactored into Steps 01–08 for modularity and reproducibility |
| `network_utils.py` | Original utility function library; functions integrated into individual step scripts and `network_config.py` |
