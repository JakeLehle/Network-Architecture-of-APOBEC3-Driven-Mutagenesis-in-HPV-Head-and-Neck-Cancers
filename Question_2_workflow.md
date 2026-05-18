## Question 2: Can differential co-expression network analysis identify gene communities associated with elevated SBS2 mutagenesis in A3-matched HNSCC tumors?

### Rationale

Figure 1 established that A3 expression is necessary but not sufficient for SBS2 mutagenesis in HNSCC, and that the highest mutational burden requires the coordinated activity of both A3B and A3A alongside as-yet-unidentified cofactors. This question exploits a key observation from Figure 1: tumors with similar A3A+A3B expression have vastly different SBS2 levels. By constructing Spearman co-expression networks separately for high-SBS2 and low-SBS2 tumors matched for A3 expression, then computing the differential network (HIGH - LOW), we can identify gene communities whose co-expression relationships are specifically rewired in the high-mutagenesis condition.

The pipeline was originally developed by Dr. Mohadeseh Soleimanpour (Texas Biomedical Research Institute) as a monolithic Python script. For this paper, the pipeline was refactored into modular steps (Step01-Step08), with adjustments to gene feature selection, community detection parameters, and figure generation.

### Data Sources

| Source | Description |
|--------|-------------|
| `TCGA_master_FPKM_UQ.tsv` | Pan-cancer FPKM-UQ expression matrix from Question 1 |
| `TCGA_SBS_signature_counts.tsv` | SigProfiler v3.4 COSMIC SBS counts (86 signatures) |
| `Mutation_Table_Tumors_TCGA.tsv` | RNA-to-WES barcode crosswalk for SBS merge |

### Directory Structure

```
scripts/FIG_2/
├── RUN_NETWORK_PIPELINE.sh              # SLURM batch script (Steps 01-08)
├── network_config.py                    # Centralized parameters (V4)
│
├── Step01_Load_Clean_TCGA.py            # Load & clean expression
├── Step02_Merge_SBS_Signatures.py       # Merge expression with SBS signatures
├── Step03_Differential_Expression.py    # Group definition & DE (fixed FDR, A3A/B force-keep)
├── Step04_Correlation_Networks.py       # TOP/BOTTOM/DIFF correlation networks
├── Step04.1_Sweep_DIFF_Threshold.py     # Diagnostic threshold sweep
├── Step05_Community_Detection.py        # Max frag rate threshold + full-network Leiden
├── Step06_Centrality_Metrics.py         # Per-gene centrality metrics
├── Step07_Generate_Figure2_Panels.py    # Original figure gen (superseded)
├── Step08_Pipeline_Summary.py           # Summary & KEGG enrichment
│
├── Compute_Node_Importance_Scores.py    # Post-pipeline: intra/inter scoring
├── Generate_Figure2_Panels.py           # Post-pipeline: publication figures (V4 two-tier)
│
└── TROUBLESHOOTING/
    ├── network_community_analysis_pipeline.py  # Original monolithic (Dr. Soleimanpour)
    ├── Diagnostic_Network_Community_Audit.py   # Harris/DDR/marker cross-reference
    ├── KEGG_A3_Neighborhood_Diagnostic.py      # Community + A3 sub-network KEGG (NEW V4)
    ├── RUN_KEGG_DIAGNOSTIC.sh                  # SLURM wrapper for KEGG diagnostic
    └── Diagnostic_DE_Threshold_Comparison.py   # BH-FDR comparison (raw vs adj)
```

### Configuration (V4)

| Parameter | Value | Notes |
|-----------|-------|-------|
| `FORCE_KEEP_A3` | True | A3A and A3B only (not full family) |
| `RAW_P_THRESHOLD` | 0.05 | TCGA stays raw p (0 genes pass FDR) |
| `DIFF_THRESHOLD` | 0.70 | Fallback only; auto-selection preferred |
| `SWEEP_THRESHOLDS` | 0.30-0.90 step 0.05 | For max fragmentation rate |
| `COMMUNITY_RESOLUTIONS` | [0.1-0.8 step 0.1] | Composite score selection |
| `USE_LARGEST_COMPONENT` | False | Full-network Leiden |
| `MIN_COMMUNITY_SIZE` | 10 | Satellites preserved below this |
| `TARGET_BIG_COMMUNITIES` | 14 | Within-large-component merge |

### Pipeline Overview (V4)

```
Step 1: Load & clean TCGA expression (10,603 tumors, 60,616 genes)
    │
Step 2: Merge with SigProfiler v3.4 SBS counts (8,276 matched, 426 HNSC)
    │
Step 3: A3-matched groups (53/53) + DE (raw p<0.05, A3A/B force-kept -> 1,638 genes)
    │
Step 4: Spearman TOP/BOTTOM/DIFF matrices (1,638 x 1,638)
    │
Step 5: Max fragmentation rate -> 0.65 | Full-network Leiden -> res 0.50
    │    (1,131 nodes, 2,599 edges, 11 main + 33 satellite communities)
    │
Step 6: Centrality metrics
    │
Step 8: KEGG enrichment
    │
Post: Node scores + Figure panels + KEGG A3 neighborhood diagnostic
```

### Step 3: V4 Changes

- BH-FDR bug fixed (gene order scramble in monotonicity enforcement step)
- Force-keep scoped to A3A + A3B only (not all 7 A3 genes)
- A3A: p=0.193 (fails DE, force-kept); A3B: p=0.483 (fails DE, force-kept)
- This is expected: A3-controlled group design prevents A3 from being DE by construction
- 1,636 DE genes + 2 force-kept = 1,638 total

### Step 5: V4 Major Rewrite (Three Changes)

**1. Max Fragmentation Rate Threshold:**
Among thresholds where A3A+A3B have degree >= 1, compute forward delta-comp between consecutive steps. Select the upper threshold of the interval with the largest positive delta-comp. Tiebreaker: lower upper threshold.

TCGA result: 0.65 (delta-comp=+25 at 0.60->0.65 interval). At 0.65: 1,131 nodes, 2,599 edges, A3A degree=1, A3B degree=4.

**2. Full-Network Leiden:**
Leiden runs on all 34 connected components (not LCC only). Large components subdivided; small satellites become their own communities.

**3. Component-Aware Merge:**
Satellites (entire small components) preserved as-is. Only within-large-component communities go through top-k merge. Gene list CSV includes `is_satellite` column.

### Step 5: Threshold Sweep (TCGA)

| Thresh | Nodes | Edges | Comp | A3A deg | A3B deg | Selected? |
|--------|-------|-------|------|---------|---------|-----------|
| 0.55 | 1,605 | 12,415 | 1 | 7 | 11 | |
| 0.60 | 1,459 | 5,926 | 9 | 3 | 8 | |
| **0.65** | **1,131** | **2,599** | **34** | **1** | **4** | **YES (delta=+25)** |
| 0.70 | 733 | 1,140 | 48 | 0 | 1 | A3A lost |

### KEGG Results (V4)

| Community | Genes | A3 | Sig KEGG | Top Term |
|-----------|-------|----|----------|----------|
| C0 | 281 | -- | 1 | Cardiac muscle contraction (0.049) |
| C1 | 257 | -- | 2 | Dilated cardiomyopathy (0.047) |
| C2 | 134 | A3B | 0 | Chronic myeloid leukemia (0.79) |
| C3 | 130 | -- | 0 | Neutrophil extracellular trap (0.61) |
| C4 | 121 | -- | 0 | Gastric cancer (0.25) |
| C5 | 44 | -- | 0 | Glycosaminoglycan biosynthesis (0.24) |
| C6 | 33 | -- | 1 | HSV1 infection (3.36e-05) |
| C7 | 23 | -- | 1 | Basal transcription factors (0.018) |
| C8 | 16 | -- | 14 | Epithelial cell signaling H. pylori (9.17e-04) |
| C43 | 2 | A3A | -- | Too small for enrichment |

**A3B sub-network KEGG:** Positive (72 genes): 0 sig. Negative (61 genes): 0 sig (ribosome at 0.056).

### Figure 2 Panels (V4)

| Panel | Content | Script |
|-------|---------|--------|
| A | Methods overview: HIGH - LOW = DIFF -> network | `Generate_Figure2_Panels.py` |
| B | Full network, two-tier layout (main + satellite ring) | `Generate_Figure2_Panels.py` |
| B2 | LCC-only backup (drops satellites) | `Generate_Figure2_Panels.py` |
| C | Community zoom insets (A3B comm, epithelial signaling, ZNF) | `Generate_Figure2_Panels.py` |
| Supp | SBS2 vs A3 sum selection plot | `Generate_Figure2_Panels.py` |

### Key Results (V4)

| Metric | Value |
|--------|-------|
| HNSC tumors | 426 |
| Groups | 53 HIGH / 53 LOW |
| DE genes | 1,638 (1,636 + 2 force-kept) |
| DIFF threshold | 0.65 (max fragmentation rate) |
| Network | 1,131 nodes, 2,599 edges |
| Components | 34 (LCC = 1,054, 93.2%) |
| Resolution | 0.50 (composite: mod x ARI x evenness) |
| Modularity | 0.49 |
| Communities | 11 main + 33 satellite |
| A3B | C2 (134 genes, degree=4), 0 sig KEGG |
| A3A | C43 (2-gene satellite, degree=1) |
| Harris interactors | 0/175 |
| Cell-type markers | 1/24 |

### Narrative (V4)

Bulk TCGA network does not resolve A3-specific pathway biology. A3B's 134-gene community has no significant KEGG enrichment at any level (community, positive sub-network, negative sub-network). A3A is in a 2-gene satellite. 0/175 Harris interactors recovered. The null result motivates single-cell resolution in Figure 4.
