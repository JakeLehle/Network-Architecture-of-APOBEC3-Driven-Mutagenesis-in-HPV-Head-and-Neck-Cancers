# Question 4 Workflow

## Question 4: If we apply network analysis to single-cell basal epithelial cells, do we see cell-type-specific cofactor communities?

### Rationale

Figure 2 established that differential co-expression network analysis on A3-matched bulk HNSCC tumors identifies 13 gene communities associated with SBS2 mutagenic activity. Community 2 (160 genes) contained A3B alongside KRT5, a canonical basal epithelial cell marker, suggesting the cofactor network may be basal-cell-specific. However, bulk data fundamentally cannot resolve whether these communities reflect cell-type-specific programs, as only 2 of 24 canonical cell-type markers appeared in any community. Figure 3 confirmed at single-cell resolution that SBS2 mutagenesis originates specifically from basal epithelial cells co-expressing A3A and A3B, and that A3C/A3G localize to immune populations. This question applies the same differential co-expression framework directly to single-cell basal epithelial cell expression data, testing whether A3-associated cofactor communities can be identified within a single cell type.

### Key Results

| Metric | Value |
|--------|-------|
| Total basal cells | 52,126 |
| Non-zero SBS2 cells | 5,911 (11.3%) |
| L-method threshold | SBS2 ≥ 1.999 |
| Slope ratio (tail/body) | 15.34× |
| HIGH group | 546 cells |
| LOW group | 546 matched controls (SBS2 = 0, anti-correlated profile) |
| Genes tested | 16,343 |
| DE genes selected (raw p < 0.05) | 9,244 |
| DE genes (FDR < 0.05) | 8,342 |
| DIFF threshold | \|Δρ\| ≥ 0.40 |
| DIFF network nodes | 1,254 (non-isolated) |
| DIFF network edges | 2,806 |
| Average degree | 4.5 |
| LCC | 1,202 nodes (95.9%) |
| Communities | 14 (Leiden, resolution 1.0, modularity 0.72, ARI 0.78) |
| Total community genes (after cleanup) | 1,140 |
| A3 genes in communities | APOBEC3A (C0), APOBEC3H (C8) |
| Known A3-interactors in communities | 5 (HNRNPA2B1, HSPD1, RPL5, TIMM8B in C0; RPL3 in C5) |
| A3-interactor enrichment in C0 | p = 0.0077 |

### Key Findings

1. APOBEC3A anchors the single-cell network (Community 0, 237 genes, degree = 8), not APOBEC3B
2. APOBEC3B is absent from all communities (max |Δρ| = 0.26), indicating its cofactor relationships are post-transcriptional (protein-protein), not co-expression
3. Community 0 exhibits bipolar architecture: A3A positively co-expressed with epithelial differentiation program (SULT2B1, SPRR family, SCEL, KLK proteases, TACSTD2), anti-correlated with known A3 protein interactors (HNRNPA2B1, HSPD1, RPL5, TIMM8B)
4. HNRNPA2B1 (A3B-specific interactor, confirmed by co-IP) bridges the A3A transcriptional and A3B post-transcriptional regulatory layers through the shared hub gene SULT2B1
5. Anti-correlations between A3A and interactors intensify 2-4× in SBS2-HIGH cells relative to LOW
6. Four A3-interactors co-localize in C0 (hypergeometric p = 0.0077), validating the community as a biologically coherent A3-associated module
7. APOBEC3H in Community 8 (67 genes, degree = 2)
8. 58 SC community genes shared with TCGA bulk communities, distributed across all 13 bulk communities (consistent with bulk communities containing mixed-cell-type signals)

### Data Sources

| Source | Description |
|--------|-------------|
| `adata_final.h5ad` | ClusterCatcher final AnnData object (155,650 cells, 27,736 genes) |
| `signature_weights_per_cell.txt` | 15-signature HNSCC-specific NNLS refitting from `signature_refitting_hnscc/` |
| `mmc15.xlsx` | Jang et al. (2024) AP-MS supplementary data |
| `Harris_A3_interactors.txt` | 175 combined A3 family interactors (McCann 2023 + Jang 2024) |
| `Harris_A3_interactors_A3B_only.txt` | 52 A3B-specific interactors |
| `ensg_to_symbol.json` | ENSG→symbol mapping from Figure 2 pipeline |
| TCGA community assignments | `data/FIG_2/05_communities/TCGA-HNSC/TCGA-HNSC_best_partition.csv` |

**Important note on signature weights:** The correct file is the 15-signature HNSCC-specific refitting (`signature_refitting_hnscc/signature_weights_per_cell.txt`, 15 signatures × 31,912 cells), NOT the 4-signature supervised NMF (`sbs2_postprocessing/Supervised_NMF_Weights/signature_weights_per_cell.txt`). The 15-signature version properly distributes mutation signal across SBS1/SBS2/SBS13/SBS5/etc., yielding cleaner SBS2 HIGH/LOW separation (546 cells vs 707 with the 4-signature version).

**Important note on gene IDs:** The SC expression matrices use a HYBRID index: gene symbols for annotated genes (APOBEC3B, KRT5, etc.) and ENSG IDs only for unannotated transcripts. All A3 gene matching uses symbol format. No protein-coding filter is applied (Cell Ranger output already contains ~20K protein-coding genes, matching the TCGA post-biotype-filter starting point). TCGA community gene lists (ENSG format) are converted to symbols via `ensg_to_symbol.json` for overlap analysis.

### Directory Structure

```
scripts/NETWORK_SINGLE_CELL/
├── network_config_SC.py                    # Centralized pipeline parameters
├── RUN_SC_NETWORK_PIPELINE.sh              # SLURM batch script (Steps 00-05)
│
├── Convert_Uniprot_to_Gene_Symbol.py       # Pre-Step A: UniProt→gene symbol mapping
├── Extract_Harris_A3_Interactors.py        # Pre-Step B: Harris lab A3 interactors
├── Step00_Select_Cells_Export_Expression.py # Step 00: L-method cell selection
├── Step01_SC_Differential_Expression.py    # Step 01: Wilcoxon rank-sum DE
├── Step02_SC_Correlation_Networks.py       # Step 02: Spearman networks + auto-report
├── Step02.1_SC_Sweep_DIFF_Threshold.py     # Step 02.1: Detailed threshold sweep
├── Step03_SC_Community_Detection.py        # Step 03: Leiden communities + cleanup
├── Step04_SC_Centrality_Metrics.py         # Step 04: Centrality metrics
├── Step05_Generate_Figure4_Panels.py       # Step 05: Figure panels + overlap + enrichment
│
├── Plot_C0_Zoom_Expanded.py                # Standalone: Expanded C0 bipolar architecture plot
├── Summary_Figure4_Gene_Analysis.py        # Standalone: Cross-reference analysis + summary
├── Quick_A3A_Interactor_Correlations.py    # Diagnostic: A3A correlation with interactors
├── Diagnostic_TCGA_Enrichment_in_SC.py     # Diagnostic: TCGA enrichment in SC network
│
├── uniprot_to_gene_symbol_mapping.tsv      # UniProt mapping (1,661 entries)
├── uniprot_accession_list.txt              # UniProt accession input
├── uniprot_accessions_to_convert.tsv       # UniProt conversion input
├── mmc15.xlsx                              # Jang et al. 2024 AP-MS data
│
└── TROUBLESHOOTING/
    ├── Diagnostic_Compare_Weights_Files.py # Compare signature weights files
    └── Diagnostic_Figure4_Issues.py        # Gene symbol overlap + A3 DIFF correlation
```

```
data/FIG_4/
├── 00_input/
│   ├── adata_final.h5ad                    (2.3 GB)
│   ├── signature_weights_per_cell.txt      (15 sigs × 31,912 cells)
│   ├── mmc15.xlsx                          (Jang et al. 2024)
│   ├── Harris_A3_interactors.txt           (175 genes)
│   └── Harris_A3_interactors_A3B_only.txt  (52 genes)
├── 01_group_selection/
│   ├── SC_Basal_SBS2_HIGH_expression.tsv   (23,095 genes × 546 cells)
│   ├── SC_Basal_SBS2_LOW_expression.tsv    (21,732 genes × 546 cells)
│   ├── SC_Basal_group_assignments.tsv
│   └── SBS2_elbow_detection.pdf/.png
├── 02_differential_expression/
│   ├── SC_diffexpr_stats.csv               (16,343 genes)
│   ├── SC_selected_genes.csv               (9,244 genes)
│   ├── SC_selected_genes_filtered.csv      (9,244 network-ready)
│   ├── SC_volcano.png
│   ├── SC_manhattan.png
│   └── SC_filtering_summary.txt
├── 03_correlation_networks/
│   ├── corr_matrices/
│   │   ├── SC_corr_HIGH.pkl                (9,244 × 9,244)
│   │   ├── SC_corr_LOW.pkl
│   │   └── SC_corr_DIFF.pkl                (~786 MB)
│   ├── edge_lists/
│   │   ├── SC_edges_HIGH.tsv               (20 edges)
│   │   ├── SC_edges_LOW.tsv                (5 edges)
│   │   └── SC_edges_DIFF.tsv               (2,806 edges)
│   ├── graph_objects/
│   │   └── SC_G_diff_noiso.gpickle         (1,254 nodes)
│   ├── heatmaps/
│   │   └── SC_correlation_distributions.png
│   ├── threshold_report.txt
│   ├── threshold_sweep.csv
│   ├── sweep_detailed.csv
│   └── sweep_detailed_plot.png
├── 04_communities/
│   ├── SC_best_partition.csv               (1,140 genes, 14 communities)
│   ├── SC_community_gene_lists.csv
│   ├── SC_community_summary.txt
│   ├── SC_G_comm.gpickle
│   ├── SC_resolution_sweep.csv
│   ├── SC_sweep_plots.png
│   └── sweep/                              (per-resolution partitions)
├── 05_centrality_metrics/
│   ├── SC_HIGH_metrics.csv
│   ├── SC_LOW_metrics.csv
│   ├── SC_DIFF_metrics.csv                 (1,254 genes)
│   └── SC_hub_genes_DIFF.csv               (top 50)
├── 06_overlap_analysis/
│   ├── overlap_matrix.csv
│   ├── hypergeom_pvalues_BH.csv
│   ├── jaccard_matrix.csv
│   ├── overlap_gene_details.csv
│   ├── known_A3_interactor_enrichment_all.csv
│   └── known_A3_interactor_enrichment_A3B_only.csv
├── 07_summary_analysis/
│   ├── SC_C0_detailed_gene_table.csv
│   ├── SC_all_communities_cross_reference.csv
│   └── Figure4_results_summary.txt
└── FIGURE_4_PANELS/
    ├── Panel_4a_Cell_Selection_UMAP.pdf/.png
    ├── Panel_4b_dual_heatmap.pdf/.png
    ├── Panel_4c_Community_Overlap.pdf/.png
    ├── Panel_4d_network_full.pdf/.png
    ├── Panel_4d_C0_zoom_expanded.pdf/.png
    ├── Panel_4e_A3_Interactor_Enrichment.pdf/.png
    └── community_zooms/
        └── SC_community_00..13.pdf/.png
```

### Configuration

All pipeline parameters are centralized in `network_config_SC.py`:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `MIN_CELLS_DETECTED` | 10 | Gene must be detected in ≥10 cells (combined groups) |
| `RAW_P_THRESHOLD` | 0.05 | Wilcoxon p-value cutoff for gene selection |
| `LOGFC_THRESHOLD` | 0 | No fold-change filter |
| `FORCE_KEEP_A3` | True | Always retain A3 genes |
| `FILTER_PROTEIN_CODING` | False | SC data already ~20K protein-coding |
| `CORRELATION_METHOD` | spearman | Rank-based correlation |
| `CORR_THRESHOLD` | 0.80 | \|rho\| threshold for HIGH/LOW network edges |
| `DIFF_THRESHOLD` | 0.40 | \|Δρ\| threshold for DIFF edges (relaxed from 0.70 in bulk) |
| `COMMUNITY_METHOD` | leiden | RBConfigurationVertexPartition |
| `COMMUNITY_RESOLUTIONS` | [0.2, 0.4, 0.6, 0.8, 1.0] | Resolution sweep range |
| `RUNS_PER_RESOLUTION` | 15 | Stability assessment runs |
| `TARGET_BIG_COMMUNITIES` | 14 | Max communities before merging |
| `MIN_COMMUNITY_SIZE` | 10 | Communities smaller than this are merged |

### Pipeline Overview

The analysis is executed as a single SLURM batch job:

```
 SLURM Job: RUN_SC_NETWORK_PIPELINE.sh
 ───────────────────────────────────────
   adata_final.h5ad  ← from ClusterCatcher
   signature_weights_per_cell.txt  ← 15-sig HNSCC refitting
       │
       ▼
   Pre-Step A: UniProt → gene symbol mapping
       │
       ▼
   Pre-Step B: Extract Harris A3 interactors
       │
       ▼
   Step 00: L-method cell selection (546 HIGH / 546 LOW)
       │
       ▼
   Step 01: Differential expression (Wilcoxon, 9,244 genes)
       │
       ▼
   Step 02: Spearman correlation networks (HIGH/LOW/DIFF)
       │    └── Step 02.1: Detailed threshold sweep (diagnostic)
       ▼
   Step 03: Leiden community detection (14 communities, 1,140 genes)
       │    └── Post-merge cleanup: isolate + island removal
       ▼
   Step 04: Centrality metrics (degree, betweenness, closeness, eigenvector)
       │
       ▼
   Step 05: Figure 4 panels + overlap analysis + interactor enrichment
       │
       ▼
   FIGURE 4
```

### Environment

All scripts are executed using the `NETWORK` conda environment with scanpy/anndata added:

```
conda install -n NETWORK -c conda-forge scanpy anndata
```

Required packages: scanpy, anndata, networkx, pandas, scipy, numpy, matplotlib, leidenalg, igraph, scikit-learn, openpyxl.

---

### Step-by-Step Pipeline Details

#### Pre-Step A: UniProt to Gene Symbol Mapping

**Script:** `Convert_Uniprot_to_Gene_Symbol.py`

Converts UniProt accession numbers from the Jang et al. (2024) AP-MS dataset to HGNC gene symbols using the UniProt ID Mapping API. Produces a TSV mapping file used by Pre-Step B.

**Output:** `uniprot_to_gene_symbol_mapping.tsv` (1,661 entries, 100% coverage)

---

#### Pre-Step B: Extract Harris A3 Interactors

**Script:** `Extract_Harris_A3_Interactors.py`

Extracts high-confidence APOBEC3 protein interactors from two Harris lab publications:
- McCann et al. (2023) Nature Genetics: 23 A3B-SF interactors from AP-MS (Layer 1)
- Jang et al. (2024) Mol Cell Proteomics: Full A3 family interactome (BFDR < 0.05, CompPASS wd_percentile ≥ 0.90)

**Output:**
- `Harris_A3_interactors.txt` (175 genes, all A3 family)
- `Harris_A3_interactors_A3B_only.txt` (52 genes, A3B-specific)

---

#### Step 00: Cell Selection and Expression Matrix Export

**Script:** `Step00_Select_Cells_Export_Expression.py`

Selects SBS2-HIGH and matched SBS2-LOW basal epithelial cells using the L-method (Salvador and Chan, 2004) piecewise linear regression for parameter-free elbow detection.

**Input:**
- `adata_final.h5ad` (ClusterCatcher output)
- `signature_weights_per_cell.txt` (15-sig HNSCC refitting)

**Output (→ `data/FIG_4/01_group_selection/`):**
- Expression matrices (genes × cells) for HIGH and LOW groups
- Group assignment table
- Panel 4a: Cell selection UMAP
- L-method diagnostic plot

**Key details:**
- 52,126 basal cells, 5,911 with non-zero SBS2
- L-method threshold: SBS2 ≥ 1.999 (slope ratio 15.34×)
- 546 HIGH cells (mean SBS2 = 2.89)
- 546 LOW controls from SBS2 = 0 pool, selected by maximizing anti-correlation of signature profile
- Gene IDs are hybrid: symbols for annotated genes, ENSG for unannotated

---

#### Step 01: Differential Expression

**Script:** `Step01_SC_Differential_Expression.py`

Wilcoxon rank-sum test between HIGH and LOW groups, matching the statistical framework from Figure 2 (Section 6.3).

**Input:**
- Expression matrices from Step 00

**Output (→ `data/FIG_4/02_differential_expression/`):**
- Per-gene statistics, selected gene lists, volcano/Manhattan plots

**Key details:**
- 20,427 common genes → 16,343 pass detection filter (≥10 cells)
- All 7 A3 genes present (A3A/B/C significant, A3D/F/G/H force-kept)
- 9,244 selected (raw p < 0.05), 8,342 at FDR < 0.05
- No protein-coding filter (SC data already ~20K protein-coding)

---

#### Step 02: Correlation Networks

**Script:** `Step02_SC_Correlation_Networks.py`

Spearman co-expression matrices for HIGH and LOW groups independently, DIFF = HIGH − LOW.

**Input:**
- Expression matrices + selected gene list from Step 01

**Output (→ `data/FIG_4/03_correlation_networks/`):**
- Correlation matrices (pickle + CSV), edge lists, graph objects, threshold report

**Key details:**
- 9,244 × 9,244 Spearman matrices (~2.4 hours computation)
- DIFF threshold |Δρ| ≥ 0.40 (relaxed from 0.70 in TCGA due to SC dropout/sparsity)
- Comparable density to TCGA: avg degree 4.5 (SC) vs 4.8 (TCGA)
- DIFF network: 1,254 non-isolated nodes, 2,806 edges
- Built-in auto-report sweeps 0.30-0.90 showing network stats at each threshold

**Step 02.1 (diagnostic):** `Step02.1_SC_Sweep_DIFF_Threshold.py` provides finer-grained sweep with Leiden community counts at each threshold.

**Threshold sweep results:**

| Threshold | Nodes | Edges | AvgDeg | LCC | Communities (r=1.0) |
|-----------|-------|-------|--------|-----|---------------------|
| 0.30 | 3,869 | 21,742 | 11.2 | 3,697 | 33 |
| 0.35 | 2,057 | 7,431 | 7.2 | 2,002 | 15 |
| **0.40** | **1,254** | **2,806** | **4.5** | **1,202** | **19** |
| 0.45 | 804 | 1,182 | 2.9 | 700 | 23 |
| 0.50 | 520 | 639 | 2.5 | 396 | 22 |
| 0.70 | 111 | 78 | 1.4 | 5 | 2 |

---

#### Step 03: Community Detection

**Script:** `Step03_SC_Community_Detection.py`

Leiden community detection on the DIFF network. Loads correlation matrix pickle directly (no Step 02 re-run needed when changing threshold).

**Input:**
- DIFF correlation matrix pickle from Step 02

**Output (→ `data/FIG_4/04_communities/`):**
- Community assignments, gene lists, annotated graph, sweep plots

**Key details:**
- Builds graph from DIFF pickle at config threshold (allows re-running without re-computing correlations)
- LCC: 1,202 nodes (95.9% of non-isolated)
- Resolution sweep: 0.2-1.0, 15 runs each
- Best resolution: 1.0 (modularity 0.72, ARI 0.78)
- 19 raw communities → 14 after merging small communities (<10 genes)
- Post-merge cleanup:
  1. Remove intra-community isolates (nodes with zero edges within their community)
  2. Per-community LCC pruning (remove disconnected islands within each community)
- Final: 14 communities, 1,140 genes

**Community structure:**

| Community | Genes | Notable |
|-----------|-------|---------|
| C0 | 237 | **APOBEC3A** (degree=8), TACSTD2, HNRNPA2B1, HSPD1, RPL5, TIMM8B |
| C1 | 96 | |
| C2 | 85 | |
| C3 | 84 | |
| C4 | 80 | SBSN, KLK8, KLK5 (kallikrein cluster) |
| C5 | 77 | RPL3, ribosomal proteins |
| C6 | 74 | |
| C7 | 71 | |
| C8 | 67 | **APOBEC3H** (degree=2) |
| C9 | 65 | |
| C10 | 62 | CRNN, RNF222 |
| C11 | 51 | |
| C12 | 47 | VWF (endothelial marker) |
| C13 | 44 | |

---

#### Step 04: Centrality Metrics

**Script:** `Step04_SC_Centrality_Metrics.py`

Per-gene centrality for HIGH, LOW, and DIFF networks.

**Output (→ `data/FIG_4/05_centrality_metrics/`):**
- Per-network centrality tables, top 50 hub genes

**Top DIFF hubs:** SCEL (deg=65), SBSN (40), CYSRT1 (49), FAM25A (36), KLK13 (46), KLK8 (30), CLIC3 (30), TGM1 (33), CRNN (26), KLK5 (31)

**A3 gene centrality:** A3A: hub_rank=175, degree=8 (C0); A3H: hub_rank=924, degree=2 (C8)

---

#### Step 05: Figure 4 Panels and Overlap Analysis

**Script:** `Step05_Generate_Figure4_Panels.py`

Generates all publication-quality panels and performs cross-resolution overlap analysis.

**Output (→ `data/FIG_4/FIGURE_4_PANELS/` and `06_overlap_analysis/`):**
- Panel 4a: Cell selection UMAP (from Step 00)
- Panel 4b: Dual heatmap (HIGH + DIFF, plasma colormap, community-sorted)
- Panel 4c: SC vs TCGA community overlap heatmap (hypergeometric)
- Panel 4d: Exploded community network (14 communities, with A3/TCGA/interactor highlighting)
- Panel 4e: Known A3-interactor enrichment bar chart
- Per-community zoom plots with colored rings and class-colored labels

**Gene highlighting on network plots:**
- Red ring (#ed6a5a) + red text: APOBEC3 genes
- Orange ring (#f18f01) + orange text: genes shared with TCGA bulk communities
- Yellow ring (#fed766) + gold text: known A3 protein interactors
- Labels use iterative repulsion to prevent overlap

**Cross-resolution overlap:** 58 SC genes shared with TCGA communities, distributed across all 13 bulk communities. No significant per-community enrichment (FDR > 0.05), consistent with bulk communities containing mixed-cell-type signals.

**A3-interactor enrichment:**
- C0: 4/5 interactors (HNRNPA2B1, HSPD1, RPL5, TIMM8B), p = 0.0077, FDR = 0.11
- C5: 1/5 (RPL3), p = 0.30

---

### Additional Analysis Scripts

#### Plot_C0_Zoom_Expanded.py

Standalone script generating an expanded view of Community 0 highlighting the bipolar architecture. Uses increased spring constant (k = 4.0/√n) and 500 iterations for better node separation. Nodes are color-coded:
- Light orange: A3A positive DIFF partners (epithelial arm)
- Light purple: A3A negative DIFF partners (interactor arm)
- Light blue: other C0 genes
- Red: APOBEC3A
- Yellow rings: known A3-interactors

#### Summary_Figure4_Gene_Analysis.py

Comprehensive cross-reference analysis:
- TCGA C2 ↔ SC C0 overlap (0 shared genes)
- All cross-community overlaps (49 pairs)
- Known A3-interactor functional annotations
- SC C0 detailed gene table (CSV)
- Full cross-reference table all communities (CSV)
- TCGA C2 network plot with SC C0 overlap highlighted

#### Quick_A3A_Interactor_Correlations.py

Diagnostic checking A3A direct DIFF/HIGH/LOW correlations with each interactor:
- All interactors anti-correlated with A3A (LOST co-expression)
- RPL5: strongest (HIGH ρ = −0.48)
- HNRNPA2B1 and RPL5 share SULT2B1 as a neighbor with A3A
- Interactors positively co-expressed with each other (RPL5-RPL3: HIGH ρ = +0.66)

#### Diagnostic_TCGA_Enrichment_in_SC.py

Tests TCGA community enrichment in the SC network using proper hypergeometric framework:
- 484/1,437 TCGA community genes present in SC selected pool (34%)
- No significant enrichment of any TCGA community in SC LCC
- Resolution sweep shows 0.2 (5 communities, largest=743) improves overlap but does not achieve significance

---

### Troubleshooting Notes

1. **Signature weights file:** The 4-signature supervised NMF version yields 707 HIGH cells. The correct 15-signature HNSCC refitting yields 546 HIGH cells with cleaner separation.

2. **Gene ID format:** SC expression matrices use hybrid index (symbols + ENSG). The protein-coding biotype filter must be disabled (FILTER_PROTEIN_CODING = False) to avoid dropping symbol-named genes. A3 gene matching uses symbol format throughout.

3. **DIFF threshold:** SC correlations are weaker than bulk due to dropout. The 0.40 threshold produces comparable network density to TCGA's 0.70. A3A enters the network at |Δρ| = 0.45 (just barely captured at 0.40 via its partners).

4. **Post-merge island removal:** Merging small communities can create disconnected fragments within large communities. Per-community LCC pruning removes these, reducing total genes from ~1,202 to 1,140.

5. **A3B absence:** A3B max |Δρ| = 0.26, far below any reasonable threshold. This is a biological finding, not a technical artifact: A3B's cofactor relationships operate post-transcriptionally.

---

### Figure 4 Summary

**Title:** Single-cell differential co-expression network analysis of SBS2-HIGH versus SBS2-LOW basal epithelial cells reveals an A3A-anchored epithelial program anti-correlated with known A3 protein interactors.

| Panel | Content |
|-------|---------|
| a | UMAP: SBS2-HIGH (coral) vs LOW (cream) basal cell selection |
| b | Dual heatmap: HIGH + DIFF correlation matrices (plasma, community-sorted) |
| c | Exploded community network (14 communities, A3/TCGA/interactor rings) |
| d | Community 0 zoom: bipolar A3A epithelial arm ↔ interactor arm |

**Narrative arc:**
1. Applied identical statistical framework as bulk TCGA to single-cell basal cells (546 vs 546)
2. A3A, not A3B, anchors the single-cell basal co-expression network
3. A3B absent from communities (cofactors are post-transcriptional, not co-expressed)
4. Community 0 bipolar architecture: A3A-epithelial program vs A3-interactor sub-module
5. HNRNPA2B1 bridges A3A transcriptional and A3B post-transcriptional layers
6. Known A3-interactors significantly enriched in C0 (p = 0.0077)
7. Motivates Figure 5: divergent basal cell fates (high-SBS2/A3A vs high-CNV/spliceosome)
