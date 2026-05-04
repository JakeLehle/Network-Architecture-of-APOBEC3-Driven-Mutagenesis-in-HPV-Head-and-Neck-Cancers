## Question 2: Can differential co-expression network analysis identify gene communities associated with elevated SBS2 mutagenesis in A3-matched HNSCC tumors?

### Rationale

Figure 1 established that A3 expression is necessary but not sufficient for SBS2 mutagenesis in HNSCC, and that the highest mutational burden requires the coordinated activity of both A3B and A3A alongside as-yet-unidentified cofactors. This question exploits a key observation from Figure 1: tumors with similar A3A+A3B expression have vastly different SBS2 levels. By constructing Spearman co-expression networks separately for high-SBS2 and low-SBS2 tumors matched for A3 expression, then computing the differential network (HIGH - LOW), we can identify gene communities whose co-expression relationships are specifically rewired in the high-mutagenesis condition. These communities are candidate cofactor networks that regulate A3 enzymatic access to genomic DNA.

The pipeline was originally developed by Dr. Mohadeseh Soleimanpour (Texas Biomedical Research Institute) as a monolithic Python script. For this paper, the pipeline was refactored into modular steps (Step01-Step08), with adjustments to gene feature selection, community detection parameters, and figure generation to align with the manuscript narrative.

### Data Sources

| Source | Description |
|--------|-------------|
| `TCGA_master_FPKM_UQ.tsv` | Pan-cancer FPKM-UQ expression matrix generated in Question 1 (Step 3) |
| `Mutation_Table_Tumors_TCGA.tsv` | COSMIC SBS signature weights (same external input as Question 1) |

### Directory Structure

```
scripts/NETWORK/
├── RUN_NETWORK_PIPELINE.sh              # SLURM batch script (Steps 01-08)
├── network_config.py                    # Centralized pipeline parameters
│
├── Step01_Load_Clean_TCGA.py            # Step 1 -- Load & clean expression data
├── Step02_Merge_SBS_Signatures.py       # Step 2 -- Merge expression with SBS signatures
├── Step03_Differential_Expression.py    # Step 3 -- Group definition & differential expression
├── Step04_Correlation_Networks.py       # Step 4 -- Build TOP/BOTTOM/DIFF correlation networks
├── Step04.1_Sweep_DIFF_Threshold.py     # Step 4.1 -- Diagnostic: sweep DIFF thresholds
├── Step05_Community_Detection.py        # Step 5 -- Leiden community detection
├── Step06_Centrality_Metrics.py         # Step 6 -- Per-gene centrality metrics
├── Step07_Generate_Figure2_Panels.py    # Step 7 -- Original figure generation (superseded)
├── Step08_Pipeline_Summary.py           # Step 8 -- Pipeline summary & KEGG enrichment
│
├── Compute_Node_Importance_Scores.py    # Post-pipeline: intra/inter community node scoring
├── Generate_Figure2_Panels.py           # Post-pipeline: publication Figure 2 panels
│
├── enviornment_NETWORK.yml              # Conda environment specification
│
└── TROUBLESHOOTING/
    ├── network_community_analysis_pipeline.py  # Original monolithic pipeline (Dr. Soleimanpour)
    ├── network_utils.py                        # Original utility functions
    ├── Diagnostic_Network_Community_Audit.py   # Cross-references Harris interactors, DDR, markers
    ├── TROUBLESHOOT_KEGG_Enrichment.py         # Standalone KEGG retry script (429 rate-limit fix)
    ├── network_config.py.bak                   # Config backup
    ├── Step02_Merge_SBS_Signatures.py.bak      # Step02 backup
    └── Step08_Pipeline_Summary.py.bak          # Step08 backup
```

### Configuration

All pipeline parameters are centralized in `network_config.py`:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `A3_SUM_PERCENTILE` | 0.50 | Keep tumors with A3A+A3B above median |
| `SBS2_HIGH_PERCENTILE` | 0.75 | Top 25% SBS2 within high-A3 = SBS2-HIGH group |
| `SBS2_LOW_PERCENTILE` | 0.25 | Bottom 25% SBS2 within high-A3 = SBS2-LOW group |
| `MIN_SAMPLES_DETECTED` | 20 | Gene must be detected in >=20 samples |
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

The analysis is executed as a single SLURM batch job (Steps 01-08), followed by post-pipeline node scoring and figure generation:

```
 SLURM Job: RUN_NETWORK_PIPELINE.sh
 ───────────────────────────────────
   TCGA_master_FPKM_UQ.tsv  <- from Question 1
   Mutation_Table_Tumors_TCGA.tsv  <- external input
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
   Step 7: Generate Figure 2 panels (original, superseded)
       │
       ▼
   Step 8: Pipeline summary + KEGG enrichment
       │
       ▼
   Post-pipeline (run manually):
       │
       ├── Compute_Node_Importance_Scores.py  (intra/inter scoring)
       │
       └── Generate_Figure2_Panels.py  (publication figures)
            │
            ▼
        FIGURE 2
```

### Environment

All scripts are executed using the `NETWORK` conda environment (see `environments/Network.yml` or `scripts/NETWORK/enviornment_NETWORK.yml`).

---

### SLURM Job: Network Analysis Pipeline

**Batch script:** `RUN_NETWORK_PIPELINE.sh`

Submits Steps 1-8 as a serial pipeline. Requires 500 GB memory and 32 CPUs for the correlation matrix computations.

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

**Output (-> `data/FIG_2/01_cleaned_expression/`):**
- Cleaned expression matrix (pickle)
- `ensg_to_symbol.json` -- gene ID to symbol mapping used throughout the pipeline

**Key details:**
- Filters for tumor samples only (barcode positions 14-15: `01`)
- Clinical metadata columns (Project_ID, Tissue_Type, Case_ID, File_ID, Entity_ID) preserved separately
- All downstream analysis uses cleaned ENSG IDs (no version suffix)

---

#### Step 2: Merge with SBS Mutation Signatures

**Script:** `Step02_Merge_SBS_Signatures.py`

Merges the cleaned TCGA expression data with COSMIC SBS mutational signature weights using Entity_ID matching. Removes duplicated samples that arise from many-to-many barcode matches.

**Input:**
- Cleaned expression data from Step 1
- `Mutation_Table_Tumors_TCGA.tsv` -- COSMIC SBS signature weights (external input)

**Output (-> `data/FIG_2/02_merged_with_SBS/`):**
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

**Output (-> `data/FIG_2/03_differential_expression/{cancer_type}/`):**
- `{ct}_SBS2_HIGH_group.pkl` / `{ct}_SBS2_LOW_group.pkl` -- group dataframes
- `{ct}_diffexpr_stats.csv` -- per-gene Wilcoxon test results
- `{ct}_selected_genes.csv` -- genes passing significance threshold
- `{ct}_selected_genes_filtered.csv` -- genes passing network-readiness filter

**Key details:**
- A3 matching: sum A3A+A3B FPKM-UQ per tumor, keep top 50% (above median)
- Within high-A3 tumors: SBS2-HIGH = top quartile, SBS2-LOW = bottom quartile
- Result: 53 HIGH and 53 LOW tumors with comparable A3A/A3B expression
- Differential test: Wilcoxon rank-sum on log1p(FPKM-UQ)
- Gene selection: raw p < 0.05 -> 1,643 DE genes selected
- A3 family genes force-retained regardless of significance
- Network-readiness filter ensures genes are present and non-constant in both groups

---

#### Step 4: Correlation Network Construction

**Script:** `Step04_Correlation_Networks.py`

Builds Spearman co-expression networks for the SBS2-HIGH (TOP) and SBS2-LOW (BOTTOM) groups independently, then computes the differential network.

**Input:**
- Group dataframes and selected gene list from Step 3

**Output (-> `data/FIG_2/04_correlation_networks/{cancer_type}/`):**
- `corr_matrices/` -- TOP, BOTTOM, DIFF correlation matrices (pickle)
- `heatmaps/` -- correlation distribution plots
- `edge_lists/` -- network edge lists for Cytoscape/Gephi
- `graph_objects/` -- NetworkX graph pickles (G_top, G_bot, G_diff, with `_noiso` variants)

**Key details:**
- Spearman correlation computed across the 1,643 selected genes for each group
- TOP and BOTTOM networks: edges retained where |rho| >= 0.80
- DIFF network: element-wise subtraction (TOP - BOTTOM), edges retained where |delta-rho| >= 0.70
- The threshold captures gene pairs whose co-expression changed substantially between conditions
- Edge attributes: `weight` (signed delta-rho) and `abs_weight` (unsigned magnitude)
- Isolated nodes (degree = 0) removed -> DIFF network: 734 nodes, 1,144 edges

**Step 4.1 -- Diagnostic threshold sweep (optional):**

**Script:** `Step04.1_Sweep_DIFF_Threshold.py`

Sweeps DIFF thresholds from 0.3 to 0.9, reporting edge count, node count, average degree, density, and community structure at each. Used to select the 0.70 threshold as the optimal balance between network connectivity and community interpretability.

---

#### Step 5: Community Detection

**Script:** `Step05_Community_Detection.py`

Runs Leiden community detection on the largest connected component (LCC) of the DIFF network across multiple resolutions, evaluates partition stability, and merges small communities.

**Input:**
- DIFF network graph from Step 4

**Output (-> `data/FIG_2/05_communities/{cancer_type}/`):**
- `{ct}_resolution_sweep.csv` -- per-resolution modularity, ARI, NMI metrics
- `{ct}_best_partition.csv` -- final gene-to-community assignments
- `{ct}_community_gene_lists.csv` -- genes per community
- `{ct}_G_comm.gpickle` -- community graph object (LCC only)

**Key details:**
- Algorithm: Leiden (RBConfigurationVertexPartition) with `abs_weight` as edge weight
- Leiden groups genes by the **magnitude** of co-expression change, regardless of direction
- Resolution sweep: 0.2-1.0, with 15 runs per resolution (seed = 42) for stability
- Best resolution: 0.80 (modularity = 0.61)
- Communities with < 10 genes merged
- Final result: **14 communities, 616 genes** (LCC of the DIFF network)
- Community sizes range from 6 (Community 14) to 88 (Community 0)

---

#### Step 6: Centrality Metrics

**Script:** `Step06_Centrality_Metrics.py`

Computes per-gene network centrality metrics for the TOP, BOTTOM, and DIFF networks. Identifies hub genes whose network role changes between conditions.

**Input:**
- Graph objects from Step 4
- Community assignments from Step 5

**Output (-> `data/FIG_2/06_centrality_metrics/{cancer_type}/`):**
- `{ct}_{TOP|BOTTOM|DIFF}_metrics.csv` -- per-gene centrality tables
- `{ct}_hub_genes_DIFF.csv` -- top hub genes by multiple metrics

**Key details:**
- Metrics computed: degree, betweenness centrality, closeness centrality, eigenvector centrality, strength (sum of |edge weights|)
- Hub genes = genes ranking highly across multiple metrics in the DIFF network
- Used for labeling top genes in Figure 2 network panels

---

#### Step 7: Generate Figure 2 Panels (Original, Superseded)

**Script:** `Step07_Generate_Figure2_Panels.py`

The original figure generation step included in the SLURM pipeline. This script has been **superseded** by `Generate_Figure2_Panels.py` (see below) which uses intra_score-based node importance sizing. Step07 remains in the pipeline for backward compatibility but its outputs are no longer used for the manuscript.

---

#### Step 8: Pipeline Summary and KEGG Enrichment

**Script:** `Step08_Pipeline_Summary.py`

Gathers all pipeline parameters, per-step statistics, and community gene lists. Runs KEGG pathway enrichment on each community via the Enrichr API. Outputs a comprehensive report for writing the results and methods sections. Run separately after reviewing pipeline output.

**Dependencies:** `gseapy` (KEGG enrichment via Enrichr API; requires internet access)

**Input:**
- All outputs from Steps 01-06
- `network_config.py` parameters

**Output (-> `data/FIG_2/08_pipeline_summary/`):**
- `{ct}_pipeline_report.txt` -- human-readable full report
- `{ct}_pipeline_report.json` -- machine-readable report
- `{ct}_KEGG_enrichment_summary.csv` -- one row per community with top KEGG term
- `{ct}_community_XX_KEGG.csv` -- full KEGG enrichment results per community

**Key details:**
- Reports group sizes, A3A/A3B/SBS2 medians, gene counts at each filtering stage
- Reports network node/edge counts for TOP/BOTTOM/DIFF graphs
- Reports community sizes, resolution, modularity
- Cross-references 24 cell-type marker genes (12 cell types x 2 markers) against community assignments
- Cross-references A3 family gene community assignments
- KEGG enrichment via gseapy.enrichr (KEGG_2021_Human, adj. p < 0.05 significance)
- Enrichr retry logic with 15-second inter-community delay to avoid 429 rate limiting

---

### Post-Pipeline: Node Importance Scoring

**Script:** `Compute_Node_Importance_Scores.py`

Computes dual intra-community and inter-community importance scores for every node in the LCC. These scores drive node sizing and label sizing in the publication Figure 2 panels.

**Input:**
- `{ct}_G_comm.gpickle` -- community graph from Step 5
- `{ct}_best_partition.csv` -- community assignments from Step 5

**Output (-> `data/FIG_2/DIAGNOSTIC_AUDIT/`):**
- `Figure2_Node_Sizing.tsv` -- per-gene sizing table with columns: `gene_symbol`, `ensg_id`, `community`, `intra_score`, `inter_score`, `intra_degree`, `inter_degree`, `total_degree`
- `Supp_Table_Node_Scores.tsv` -- full scoring detail for supplemental
- `Node_Score_Summary_Report.txt` -- human-readable summary

**Scoring framework:**
- **INTRA score** ("local hub"): intra-community degree + strength + within-community eigenvector centrality, normalized 0-1
- **INTER score** ("bridge"): inter-community degree + strength + n_communities_connected + global betweenness centrality, normalized 0-1
- Edge classification: 776 intra-community (72.3%), 297 inter-community (27.7%)
- Top local hubs: PDCD1LG2 (C0, intra=0.997), MLF1 (C2, intra=0.967), GBP5 (C1, intra=0.967)
- Top bridges: MLF1 (inter=0.987, 8 communities), PDCD1LG2 (inter=0.998, 7 communities)
- Community 4 top hub: FAM78B (intra=5, inter=9, 5 communities)
- A3B: degree=1, bottom of all rankings (passenger in bulk network)

---

### Post-Pipeline: Publication Figure 2 Panels

**Script:** `Generate_Figure2_Panels.py`

Generates publication-quality Figure 2 panels using intra_score-based node and label sizing from `Figure2_Node_Sizing.tsv`. This script supersedes `Step07_Generate_Figure2_Panels.py`.

**Input:**
- All outputs from Steps 01-06
- `data/FIG_2/DIAGNOSTIC_AUDIT/Figure2_Node_Sizing.tsv` from `Compute_Node_Importance_Scores.py`

**Output (-> `data/FIG_2/FIGURE_2_PANELS/{cancer_type}/`):**
- `{ct}_Panel_A_selection.pdf/.png` -- Figure 2a
- `{ct}_Panel_B_network_full.pdf/.png` -- Figure 2b
- `{ct}_Panel_C_C4_zoom.pdf/.png` -- Figure 2c
- `{ct}_Supplement_methods_overview.pdf/.png` -- Supplemental methods graphical abstract
- `community_zooms/` -- per-community zoomed network subplots (PNG + PDF)

**Panel A -- Tumor selection plot (Figure 2a):**
- Scatter of summed A3A+A3B expression vs SBS2 weight for all 426 HNSC tumors
- HIGH (n=53, coral) and LOW (n=53, teal) groups highlighted
- Point size s=100 for all points (consistent with GitHub walkthrough)
- Colored background regions matching Figure 1a convention
- Data-driven broken x-axis for outlier accommodation
- Font sizes: axis labels 30pt, tick labels 22pt, legend 20pt

**Panel B -- Full community network (Figure 2b):**
- Two-tier community-aware exploded layout: super-graph positions community centers (`INTER_SCALE = 6.0`), local spring layout positions nodes within each community (`INTRA_SCALE = 1.0`)
- Node size scaled by intra_score: `100 + 1160 * intra_score`
- Label font size scaled by intra_score: `15 + 28 * intra_score`
- Top 3 hub genes per community labeled (by intra_score), A3B and A3H always labeled
- Label repulsion algorithm prevents overlap (300 iterations, push_strength=0.08)
- Labels drawn with semi-transparent white background boxes (alpha=0.42) and shifted slightly above nodes
- Intra-community edges colored by differential correlation sign (red/firebrick = gained, blue/steelblue = lost), width `0.8 + 4.0 * |w|`
- Inter-community edges shown as dashed gray lines
- A3 gene nodes drawn with coral fill, enlarged 1.8x, with red ring
- Community center annotations (Cx) positioned above each cluster

**Panel C -- Community 4 (A3B) zoom (Figure 2c):**
- Expanded spring layout (k=4.0/sqrt(n), 500 iterations, scale=2.0)
- Node size: `750 + 5000 * intra_score`
- Label font size: `18 + 28 * intra_score` (A3B forced to 30pt)
- All nodes labeled (community is ~50 genes)
- Edge thickness doubled for zoom readability: `2.0 + 10.0 * |w|`
- Label repulsion (250 iterations, push_strength=0.07)
- Title includes top KEGG enrichment (Wnt signaling / Basal cell carcinoma)

**Supplement -- Methods graphical abstract:**
- Four panels in a row: HIGH correlation heatmap, LOW correlation heatmap, DIFF correlation heatmap, mini network
- Operator annotations between panels: " - ", " = ", " >> "
- Tells the visual story: "this minus this equals this, which becomes the network"
- All text 24-32pt range
- Genes ordered by hierarchical clustering within and between communities

**Per-community zooms (supplemental):**
- All 14 communities with >= 5 nodes plotted individually
- Node size: `500 + 3500 * intra_score`
- Label font: `15 + 21 * intra_score`
- Edge thickness: `1.6 + 8.0 * |w|`
- Saved as both PNG and PDF at 300 DPI

**Style conventions (all panels):**
- All hex color codes (no named colors)
- Font sizes 28-34pt for main text
- Saved as both PDF and PNG at 300 DPI
- Label background alpha = 0.42 for node visibility
- Initial label offset [0, 0.16] above nodes

---

### Diagnostic Audit

**Script:** `TROUBLESHOOTING/Diagnostic_Network_Community_Audit.py`

Cross-references the 14 network communities against external gene sets to characterize community composition.

**Output (-> `data/FIG_2/DIAGNOSTIC_AUDIT/`):**
- `Supp_Table_Network_Community_Summary.tsv`
- `Supp_Table_Node_Manifest.tsv`
- `Network_Community_Audit_Report.txt`

**Key findings:**
- Harris A3 interactor cross-reference: 0/175 total in network, 0/52 A3B-specific
- DDR genes in network: 3 (FANCA in C0, CHEK1 in C5, ATRIP in C7), scattered across communities
- Chromatin remodelers: 0 in network
- Cell-type markers: 0/24 in any community
- Interpretation: Harris interactors are absent from the bulk co-expression network but appear at single-cell resolution (HNRNPA2B1, HSPD1, RPL5 in SC C0), consistent with post-transcriptional cofactors that only emerge when cell-type heterogeneity is resolved

---

### Figure 2 Summary

**Title:** Differential co-expression networks in A3-matched HNSCC tumors reveal gene communities associated with SBS2 mutagenic activity.

| Panel | Content | Script |
|-------|---------|--------|
| a | A3-matched tumor selection: SBS2-HIGH vs SBS2-LOW groups | `Generate_Figure2_Panels.py` |
| b | Full exploded Leiden community network with intra_score-based hub labels | `Generate_Figure2_Panels.py` |
| c | Community 4 (A3B) zoomed network with all gene labels | `Generate_Figure2_Panels.py` |
| Supplement | Methods graphical abstract: HIGH -> LOW -> DIFF -> network | `Generate_Figure2_Panels.py` |

**Narrative arc:**
1. A3-matched groups isolate cofactor-driven differences from A3 expression variation (panel a)
2. 14 gene communities identified with robust modularity (0.61) in the differential network (panel b)
3. A3B lands in Community 4 (50 genes, degree=1) alongside Wnt signaling genes (WNT5B, WNT7B, WNT16) and hub FAM78B, suggesting a basal cell identity context rather than a direct regulatory mechanism (panel c)
4. A3H lands in Community 11 (25 genes) with glycosaminoglycan biosynthesis enrichment
5. A3B is a passenger in the bulk network (degree=1, bottom of all importance rankings), consistent with cofactors being post-transcriptional
6. 0/175 Harris A3 interactors found in any community; they only appear at single-cell resolution, suggesting bulk data masks cell-type-specific cofactor interactions
7. 0/24 cell-type markers in any community, confirming bulk data cannot resolve cell-type specificity
8. **Motivates Question 3:** single-cell resolution needed to assign communities to cell types, validate A3B's basal cell context, and recover the missing interactors

### Key Results

| Metric | Value |
|--------|-------|
| Total HNSC tumors | 426 |
| SBS2-HIGH group | 53 tumors |
| SBS2-LOW group | 53 tumors |
| Genes tested | 19,101 |
| DE genes selected | 1,643 |
| DIFF network nodes | 734 |
| DIFF network edges | 1,144 |
| LCC nodes | 616 |
| Communities | 14 |
| Modularity | 0.61 |
| Resolution | 0.80 |
| Intra-community edges | 776 (72.3%) |
| Inter-community edges | 297 (27.7%) |
| Cell-type markers in communities | 0 / 24 |
| Harris A3 interactors in communities | 0 / 175 |
| A3 genes in communities | A3B in C4 (degree=1), A3H in C11 (degree=4) |

### KEGG Enrichment Highlights

| Community | Size | A3 Gene | Top KEGG Term | Adj p | Sig Terms |
|-----------|------|---------|---------------|-------|-----------|
| C0 | 88 | -- | Cardiac muscle contraction | 0.003 | 5 |
| C1 | 79 | -- | Pertussis | 0.346 | 0 |
| C2 | 77 | -- | Arrhythmogenic RV cardiomyopathy | 0.363 | 0 |
| C3 | 57 | -- | Cortisol synthesis/secretion | 0.367 | 0 |
| C4 | 50 | A3B | Basal cell carcinoma | 0.017 | 6 |
| C5 | 47 | -- | p53 signaling pathway | 0.042 | 2 |
| C6 | 41 | -- | PPAR signaling | 0.385 | 0 |
| C7 | 36 | -- | Insulin resistance | 0.056 | 0 |
| C8 | 35 | -- | Mismatch repair | 0.354 | 0 |
| C9 | 33 | -- | Osteoclast differentiation | 0.005 | 1 |
| C10 | 26 | -- | Ribosome | 0.003 | 2 |
| C11 | 25 | A3H | Glycosaminoglycan biosynthesis | 0.026 | 1 |
| C12 | 16 | -- | Breast cancer | 0.142 | 0 |
| C14 | 6 | -- | Glioma | 0.067 | 0 |

### Community 4 Key Genes

| Gene | Role | Intra Degree | Inter Degree |
|------|------|-------------|-------------|
| FAM78B | Top hub + bridge | 5 | 9 (5 communities) |
| FBXO32 | Top hub (tied) | 5 | -- |
| WNT5B | Wnt signaling | -- | -- |
| WNT7B | Wnt signaling | -- | -- |
| WNT16 | Wnt signaling | -- | -- |
| APOBEC3B | Passenger | 1 | 0 |

---

### Troubleshooting Scripts

The following scripts in `scripts/NETWORK/TROUBLESHOOTING/` are not part of the main pipeline but were used during development and diagnostic analysis:

| Script | Purpose |
|--------|---------|
| `network_community_analysis_pipeline.py` | Original monolithic pipeline by Dr. Mohadeseh Soleimanpour; refactored into Steps 01-08 for modularity and reproducibility |
| `network_utils.py` | Original utility function library; functions integrated into individual step scripts and `network_config.py` |
| `Diagnostic_Network_Community_Audit.py` | Cross-references communities against Harris interactors (175 total, 52 A3B-specific), DDR genes, chromatin remodelers, and cell-type markers; all null in bulk network |
| `TROUBLESHOOT_KEGG_Enrichment.py` | Standalone script to diagnose and fix Enrichr 429 rate limiting; led to 15-second inter-community delay in Step08 |
| `network_config.py.bak` | Configuration backup from before pipeline rerun |
| `Step02_Merge_SBS_Signatures.py.bak` | Step02 backup |
| `Step08_Pipeline_Summary.py.bak` | Step08 backup (before Enrichr retry logic) |
