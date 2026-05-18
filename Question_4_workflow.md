## Question 4: Does single-cell differential co-expression network analysis resolve A3-specific cofactor biology that bulk data cannot?

### Rationale

Figure 2 established that bulk TCGA network analysis cannot resolve A3-specific pathway biology: A3B's community had no significant KEGG enrichment, 0/175 Harris interactors were recovered, and only 1/24 cell-type markers appeared in any community. Figure 3 confirmed that SBS2 mutagenesis localizes to basal epithelial cells co-expressing A3A and A3B, and revealed that high SBS2 weight associates with A3A while high CNV burden associates with A3B, suggesting divergent mutagenic fates within the basal compartment. This question applies the same differential co-expression framework to single-cell data, using a three-population design (SBS2-HIGH, CNV-HIGH, NORMAL) to construct three pairwise networks that capture distinct aspects of A3-driven mutagenesis.

### Key Results (V4: Three-Network Summary)

| Metric | SBS2 vs CNV | SBS2 vs NORMAL | CNV vs NORMAL |
|--------|-------------|----------------|---------------|
| DE genes (FDR<0.05) | 5,107 | 3,907 | 5,805 |
| DIFF threshold | 0.60 | 0.40 | 0.40 |
| Network genes | 240 | 2,948 | 5,174 |
| Edges | 747 | 23,754 | 124,279 |
| Communities (main+sat) | 5+3 | 6+17 | 5+11 |
| Resolution | 0.70 | 0.70 | 0.70 |
| Modularity | 0.37 | 0.34 | 0.25 |
| A3 in network | A3A, A3B | A3A, A3B | A3A, A3B, A3C, A3G |
| A3 community | Satellite (3 genes) | C0 (1,054 genes) | C0 (2,210 genes) |
| Harris in network | 3 | 54 | 116 |
| TCGA overlap | 15 | 149 | 232 |
| RALY rank (Harris) | Top hub C0 | Top Harris C0 | Top Harris C2 |

### Data Sources

| Source | Description |
|--------|-------------|
| `adata_final.h5ad` | ClusterCatcher AnnData (155,650 cells, 27,736 genes) |
| `three_group_assignments.tsv` | SBS2-HIGH/CNV-HIGH/NORMAL (546 each) |
| `Harris_A3_interactors.txt` | 174 A3 interactors (Harris lab) |
| `Harris_A3_interactors_A3B_only.txt` | 51 A3B-specific interactors |
| `ensg_to_symbol.json` | From Figure 2 pipeline |

### Directory Structure

```
scripts/
├── network_config_SC.py                       # Centralized SC parameters (V4)
├── RUN_THREE_NETWORK_PIPELINE.sh              # SLURM wrapper for all 3 networks (V4)
│
├── Step01_SC_Differential_Expression.py       # scanpy rank_genes_groups, FDR<0.05 (V4 rewrite)
├── Step02_SC_Correlation_Networks.py          # Spearman HIGH/LOW/DIFF matrices
├── Step03_SC_Community_Detection.py           # Max frag rate + full-network Leiden (V4 rewrite)
├── Step04_SC_Centrality_Metrics.py            # Centrality metrics
├── Compute_Node_Importance_Scores_SC.py       # Intra/inter scoring (V4 bug fixes)
│
├── Step00B_Three_Group_Selection_and_Export.py # Define 3 populations + export expression
│
├── NETWORK_SINGLE_CELL/                       # Legacy single-network scripts (superseded)
│
└── TROUBLESHOOTING/
    ├── KEGG_A3_Neighborhood_Diagnostic.py     # Community + A3 sub-network KEGG (all 4 networks)
    └── RUN_KEGG_DIAGNOSTIC.sh                 # SLURM wrapper
```

```
data/FIG_4/
├── 00_input/
│   ├── adata_final.h5ad                    (2.3 GB)
│   ├── Harris_A3_interactors.txt           (174 genes)
│   └── Harris_A3_interactors_A3B_only.txt  (51 genes)
├── 01_group_selection/
│   ├── three_group_assignments.tsv         (1,638 cells: 546 x 3)
│   ├── NETWORK_SBS2_VS_CNV/               (per-network expression matrices)
│   ├── NETWORK_SBS2_VS_NORMAL/
│   └── NETWORK_CNV_VS_NORMAL/
├── NETWORK_SBS2_VS_CNV/
│   ├── 02_differential_expression/
│   ├── 03_correlation_networks/
│   ├── 04_communities/
│   ├── 05_centrality_metrics/
│   └── DIAGNOSTIC_AUDIT/
├── NETWORK_SBS2_VS_NORMAL/
│   └── (same structure)
├── NETWORK_CNV_VS_NORMAL/
│   └── (same structure)
└── FIGURE_4_PANELS/
```

### Configuration (V4)

| Parameter | Value | Notes |
|-----------|-------|-------|
| `FORCE_KEEP_A3` | False | A3 passes FDR naturally in SC |
| `RAW_P_THRESHOLD` | 0.05 | Used for DE stats, but FDR<0.05 for selection |
| `DIFF_THRESHOLD` | 0.40 | Fallback; auto-selection preferred |
| `SWEEP_THRESHOLDS` | 0.30-0.90 step 0.05 | Max fragmentation rate |
| `COMMUNITY_RESOLUTIONS` | [0.1-0.8 step 0.1] | Composite score |
| `USE_LARGEST_COMPONENT` | False | Full-network Leiden |
| `MIN_COMMUNITY_SIZE` | 10 | Satellites preserved |
| `MIN_CELLS_DETECTED` | 10 | Gene detection filter |

### Pipeline Overview (V4: Three-Network)

```
 RUN_THREE_NETWORK_PIPELINE.sh
 ─────────────────────────────
 For each of 3 comparisons (SBS2_VS_CNV, SBS2_VS_NORMAL, CNV_VS_NORMAL):
     │
     ├── Copy expression matrices to standard locations
     ├── Clean previous outputs
     │
     ├── Step01: scanpy rank_genes_groups (FDR<0.05, no force-keep)
     │           Takes comparison name as command-line argument
     │
     ├── Step02: Spearman correlation matrices (HIGH/LOW/DIFF)
     │
     ├── Step03: Max fragmentation rate threshold + full-network Leiden
     │           Component-aware merge preserves satellites
     │
     ├── Step04: Centrality metrics
     │
     ├── Step04B: Node importance scores (Compute_Node_Importance_Scores_SC.py)
     │
     └── Move outputs to NETWORK_{comparison}/
```

### Step 01: V4 Rewrite (scanpy DEG)

- Replaced manual gene-by-gene mannwhitneyu loop with scanpy `rank_genes_groups`
- Loads adata_final.h5ad directly (not separate expression TSVs)
- Reads three_group_assignments.tsv for cell labels
- Takes comparison name as command-line argument (SBS2_VS_CNV, etc.)
- NO additional log normalization (data already processed)
- Uses scanpy's internal BH-FDR
- Selects at FDR < 0.05 (not raw p < 0.05)
- No force-keep: A3A/A3B pass FDR naturally in all 3 comparisons

**A3 DE results across comparisons:**

| Gene | SBS2 vs CNV | SBS2 vs NORMAL | CNV vs NORMAL |
|------|-------------|----------------|---------------|
| A3A | FDR=6.3e-87, FC=+4.47 | FDR=3.7e-109, FC=+7.73 | FDR=2.6e-24, FC=+3.25 |
| A3B | FDR=1.2e-32, FC=-3.06 | FDR=1.6e-15, FC=+2.58 | FDR=4.6e-92, FC=+5.65 |

Note: A3A dominates in SBS2_VS_NORMAL; A3B dominates in CNV_VS_NORMAL; opposite directions in SBS2_VS_CNV.

### Step 03: V4 Rewrite (Threshold + Clustering)

Same max fragmentation rate + full-network Leiden + component-aware merge as TCGA Step05.

**Threshold selections:**

| Network | Threshold | Reason | A3A deg | A3B deg |
|---------|-----------|--------|---------|---------|
| SBS2_VS_CNV | 0.60 | delta-comp=+2 at 0.55->0.60 | 1 | 2 |
| SBS2_VS_NORMAL | 0.40 | delta-comp=+17 at 0.35->0.40 | 55 | 16 |
| CNV_VS_NORMAL | 0.40 | delta-comp=+11 at 0.35->0.40 | 10 | 239 |

All three independently selected resolution 0.70.

### Compute_Node_Importance_Scores_SC.py: V4 Bug Fixes

- Fixed `load_harris_interactors()` unpacking (two separate calls, not tuple)
- Fixed TCGA cross-reference (convert ENSG to symbols via `load_ensg_to_symbol()`)
- Added `load_ensg_to_symbol` to imports

### Per-Network Results

#### NETWORK_SBS2_VS_CNV (Divergent Fates, 240 genes)

- A3A+A3B+PLAUR in 3-gene satellite (C5)
- A3A-A3B positively correlated within satellite (rho=+0.74)
- But opposite correlations OUTSIDE: A3A positive with PRRX1/UBD/KRT15, A3B negative with those same genes
- RALY: top hub of entire network (C0, degree=127, intra=0.926)
- Main communities: C0 mineral absorption (MT1M/MT1G), C1 ribosome (adj p=5.22e-11), C2 ribosome + oxidative phosphorylation
- Only 3 Harris interactors: RALY (C0), PTGES3 (C1), RPL5 (C2)

#### NETWORK_SBS2_VS_NORMAL (Mutagenic Entry, 2,948 genes)

- A3A+A3B in C0 (1,054 genes) with TACSTD2 and CD68
- A3A positive sub-network (54 genes): 7 sig KEGG including antigen processing (adj p=3.50e-03), spliceosome (0.022), protein processing in ER (0.028). HSP core: HSP90AB1, HSPA1A, HSPA1B
- A3B positive sub-network (62 genes): 9 sig KEGG including legionellosis (3.33e-03), protein processing in ER (7.54e-03). Same HSP core
- 54 Harris interactors across C0 (16), C1 (17), C2 (20)
- C1 (995 genes): endocytosis (2.54e-07), ubiquitin proteolysis (4.78e-05), H. pylori epithelial signaling (4.04e-04)
- C2 (778 genes): ribosome (6.21e-50), HIF-1 signaling (9.45e-04)
- C3 (32 genes): interferon response module (IFIH1, RSAD2, OAS1, MX2, ISG15)
- C4 (32 genes): DNA replication (1.72e-05), cell cycle

#### NETWORK_CNV_VS_NORMAL (Productive Infection, 5,174 genes)

- A3A+A3B+A3C+A3G all in C0 (2,210 genes)
- A3B degree=239, intra_score=0.969 (genuine hub)
- A3B anti-correlated with epithelial adhesion: DSG3 (-0.58), DSC2 (-0.56), KRT16 (-0.59), NOTCH3 (-0.64)
- A3B positive sub-network (41 genes): 0 sig KEGG
- 116 Harris interactors with 61 in A3 community (C0)
- A3B-specific interactors scoring high: PFDN2, MAP4, HNRNPUL1, SRSF10, RAD23B, SNRPB2, RNF126
- C1 (1,424 genes): HPV infection (adj p=7.41e-07), focal adhesion, pathways in cancer
- C2 (1,141 genes): ribosome (7.60e-33), apoptosis, protein processing in ER
- C3 (359 genes): cell cycle (1.22e-10), DNA replication (8.19e-07)
- RALY: top Harris interactor (C2, degree=375, intra=0.988)

### Cross-Network Patterns

1. **RALY** is the most consistent finding across all analyses (top Harris hub in all 3 SC networks + TCGA)
2. **Harris interactor recovery** scales with network complexity: 3 -> 54 -> 116
3. **A3A/A3B divergence** captured most specifically in SBS2_VS_CNV (opposite outside-satellite correlations)
4. **Proteotoxic stress** (HSP90AB1, HSPA1A/B) co-activated with both A3A and A3B in SBS2_VS_NORMAL
5. **Epithelial identity loss** (DSG3, DSC2, KRT16) anti-correlated with A3B in CNV_VS_NORMAL
6. **HPV infection pathway** appears only in CNV_VS_NORMAL (C1), the comparison capturing productive infection
7. **Interferon response** module appears only in SBS2_VS_NORMAL (C3)
8. **All three** independently selected resolution 0.70

### Figure 4 Panels (Planned)

| Panel | Content |
|-------|---------|
| Central | Three-population UMAP (SBS2-HIGH / CNV-HIGH / NORMAL) |
| Left | SBS2_VS_CNV: A3 satellite zoom |
| Top | SBS2_VS_NORMAL: A3 community zoom |
| Right | CNV_VS_NORMAL: A3 community zoom |
| Arrows | Double-sided between UMAP and each network |
| Supp | Full networks + non-A3 community zooms |

### Narrative Arc (V4)

1. Three-population framework captures divergent mutagenic fates (SBS2-HIGH/A3A vs CNV-HIGH/A3B)
2. SBS2_VS_CNV: A3A and A3B co-localize in a 3-gene satellite but have mirror-image external correlations, capturing the divergence at its most specific resolution
3. SBS2_VS_NORMAL: A3 activation coincides with proteotoxic stress response (antigen processing, spliceosome, ER protein processing), motivating the neoantigen analysis
4. CNV_VS_NORMAL: A3B is a genuine network hub (degree=239) anti-correlated with epithelial identity, with 116 Harris interactors recovered including A3B-specific ones involved in DNA damage response and RNA processing
5. RALY is the single most consistent finding across all analyses
6. Bulk TCGA recovered 0/175 Harris interactors; single-cell recovers up to 116, validating the resolution gain
