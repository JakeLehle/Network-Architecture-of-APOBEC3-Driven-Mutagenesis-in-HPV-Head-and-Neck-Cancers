# Question 2 Workflow

## Question 2: Does single-cell differential co-expression network analysis resolve A3 cofactor biology when each tumor population is compared to normal tissue?

### Rationale

Question 1 resolved two divergent programs in the basal compartment, A3A-linked SBS2 point mutations and A3B-linked copy-number change. Differential expression direction alone cannot take that further: a gene moving in the same direction as tumorigenesis could belong to a coordinated A3 cofactor program or could just be part of the broad transcriptional shift of becoming a tumor, and direction says nothing about how a gene tracks A3A or A3B. To assign functional roles to the known A3 interactors, we build differential co-expression networks comparing each tumor population to normal-adjacent basal cells and look for chains of genes that gain co-expression with the enzymes in tumor (candidate activators) or lose it (candidate brakes).

The lens here is tumor versus normal. An earlier SBS2-HIGH versus CNV-HIGH comparison was explored but dropped when the paper refocused on the tumor-to-normal contrast, so this question runs two networks: SBS2-HIGH vs NORMAL and CNV-HIGH vs NORMAL.

### Population design

Three basal populations of 546 cells each are defined in `Step00B`, and the two networks each compare a tumor population to NORMAL:

- **SBS2-HIGH (n=546):** from SBS2 > 0 cancer-tissue basal cells, by a composite of SBS2 weight (40%), low CNV (20%), low CytoTRACE2 stemness (20%), and high A3A/(A3A+A3B) fraction (20%).
- **CNV-HIGH (n=546):** from SBS2 == 0 cancer-tissue basal cells, by a mirror composite of high CNV (40%), total A3 (A3A+A3B) range-matched to the SBS2-HIGH mean (20%), high stemness (20%), and high A3B/(A3A+A3B) fraction (20%). After the CNV-HIGH reselection the HPV16 late-gene and anti-correlated-profile terms were removed, so selection carries no viral input and is not circular with the Question 4 virus analysis.
- **NORMAL (n=546):** random sample of 546 of the 554 normal-adjacent basal cells.

Total A3 expression is range-matched between the two tumor groups so they differ in A3A/A3B composition rather than in summed enzyme level.

### Data Sources

| Source | Description |
|--------|-------------|
| `adata_final.h5ad` | ClusterCatcher AnnData (155,650 cells, 27,736 genes) |
| `three_group_assignments.tsv` | SBS2-HIGH / CNV-HIGH / NORMAL, 546 each (1,638 cells) |
| `Harris_A3_interactors.txt` | 174 catalogued A3 interactors |
| `Harris_A3_interactors_A3B_only.txt` | A3B-specific interactor subset |
| `mmc15.xlsx` | AP-MS interactor source table (built into the interactor list) |

### Directory Structure

```
scripts/NETWORK_SINGLE_CELL/
├── network_config_SC.py                        # Centralized network parameters
├── RUN_THREE_NETWORK_PIPELINE.sh               # Compute + diagnostics (SLURM); plotting deferred
├── MAKE_FIGURES.sh                             # Plotting stage (separate from compute)
│
├── Step00B_Three_Group_Selection_and_Export.py # Define 3 populations (546 each), export per-network matrices
├── Step01_SC_Differential_Expression.py        # scanpy rank_genes_groups, FDR < 0.05
├── Step02_SC_Correlation_Networks.py           # Spearman HIGH / LOW / DIFF matrices
├── Step03_SC_Community_Detection.py            # Auto DIFF threshold + full-network Leiden
├── Step04_SC_Centrality_Metrics.py            # Centrality metrics
├── Compute_Node_Importance_Scores_SC.py        # Intra / inter community importance scoring
├── Diagnostic_A3_Interactor_Concordance.py     # Interactor concordance across the two networks (chain inputs)
├── Diagnostic_Chain_Validation_SBS2_VS_CNV.py  # Activating / inhibiting chain validation across the two tumor networks
│
├── Extract_Harris_A3_Interactors.py            # Build the 174 A3 interactor list
├── Convert_Uniprot_to_Gene_Symbol.py           # UniProt -> gene symbol mapping for the interactor list
├── mmc15.xlsx                                  # AP-MS interactor source table
├── uniprot_accession_list.txt                  # Interactor-list intermediates
├── uniprot_accessions_to_convert.tsv
├── uniprot_to_gene_symbol_mapping.tsv
│
├── Generate_Figure4_Panels.py                  # Figure 4 panels (plotting)
├── Generate_Figure4_Supplemental.py            # Supplemental network figure (plotting)
└── TROUBLESHOOTING/                            # Diagnostics and prior drafts (not documented here)
```

Note: `Step05_Generate_Figure4_Panels.py` is also present and overlaps `Generate_Figure4_Panels.py`; collapse to one canonical panel script on the per-script pass.

### Configuration (`network_config_SC.py`)

| Parameter | Value | Notes |
|-----------|-------|-------|
| DE selection | FDR < 0.05 | scanpy `rank_genes_groups`, internal BH-FDR |
| `FORCE_KEEP_A3` | False | A3A and A3B pass FDR naturally in both networks |
| DIFF threshold | auto | max fragmentation rate over a 0.30-0.90 sweep (fallback 0.40) |
| Community detection | Leiden | full-network, resolution by composite-score sweep, component-aware merge |
| `MIN_COMMUNITY_SIZE` | 10 | satellites preserved |
| `MIN_CELLS_DETECTED` | 10 | gene detection filter |

### Pipeline: two stages

The compute and the plotting are deliberately separated so the figure can be tweaked without rerunning the networks.

**Stage A, compute and diagnostics (`RUN_THREE_NETWORK_PIPELINE.sh`):**
```
Step00B  group selection + export (3 populations, per-network matrices)
   │
   ├── for SBS2_VS_NORMAL and CNV_VS_NORMAL:
   │      Step01  scanpy DE (FDR < 0.05)
   │      Step02  Spearman HIGH / LOW / DIFF correlation matrices
   │      Step03  auto DIFF threshold + full-network Leiden communities
   │      Step04  centrality metrics
   │      Compute_Node_Importance_Scores_SC  intra / inter scoring
   │
   └── Diagnostics: A3 interactor concordance + chain validation
          (these produce the activating / inhibiting chain composition that
           the figure design and the manuscript narrative are built on)
```

**Stage B, plotting (`MAKE_FIGURES.sh`):** `Generate_Figure4_Panels.py` and `Generate_Figure4_Supplemental.py`.

The interactor list itself is built once by `Extract_Harris_A3_Interactors.py` and `Convert_Uniprot_to_Gene_Symbol.py` from `mmc15.xlsx`, yielding the 174-gene catalogue used throughout.

### Key Results

A3A and A3B both pass DE naturally in each network (A3A log2 fold change 7.7 in SBS2-HIGH vs normal, falling to 1.2 in CNV-HIGH vs normal, reflecting the sharp drop in A3A induction once cells enter the productive state).

Of the 174 catalogued A3 interactors, relative to normal-adjacent tissue 52 were upregulated in SBS2-HIGH cells, 129 in CNV-HIGH, and 51 in both; only three (*HSPD1*, *RPL3*, *RPL5*) were lower in both tumor states, and none switched direction between them, so no single interactor stood out as a brake on mutagenesis from direction alone.

| Metric | SBS2-HIGH vs NORMAL | CNV-HIGH vs NORMAL |
|--------|---------------------|--------------------|
| Network genes | 2,948 | 4,886 |
| Gene groups | 23 | 38 |
| Interactors recovered (of 174) | 54 | 109 |
| A3 arrangement | A3A and A3B share a gene group | A3A and A3B separate |
| DIFF threshold / resolution | 0.40 / 0.70 | auto-selected |

**Activating chains (RALY-anchored, conserved in both):**
- SBS2-HIGH: *RALY* anchors a five-gene chain with *LCN2*, *KRT24*, *LINC00278*, and *UTY*; a second interactor, *HNRNPA2B1*, anchors a smaller chain with *CCL20* (Fig. 4b).
- CNV-HIGH: *RALY* anchors the largest chain, a thirteen-gene module of translation and metabolic genes (*CPNE1*, *EIF6*, *CA2*, *DYNLRB1*) within the A3A group (Fig. 4d).
- The partners turn over between states but *RALY* holds the same role even though CNV-HIGH cells carry no SBS2 and A3A induction has fallen sharply, marking *RALY* as a tumor-conserved candidate coactivator rather than a tracker of enzyme level.

**Inhibiting chain (no catalogued interactor):**
- A roughly seventy-gene epithelial differentiation and cornification program in the A3A group of the CNV-HIGH network (*PRSS3*, *CLIC3*, *MAB21L4*, *SBSN*, *CYSRT1*, and the *SPRR* family), tightly co-expressed in normal tissue and losing coherence in tumor, with A3A at its edge through *SPRR1A*, *SPRR2D*, and *RAB11A* (Fig. 4c). This suggests loss of a normal differentiation program as cells accumulate CNV.

**The A3 wall:**
- Nearly every edge directly linking A3A or A3B to its neighbors is negative, so each enzyme's co-expression with its partners is stronger in normal tissue even though both genes are induced in tumor. Not one enzyme edge is positive. The enzymes are induced but co-expression-decoupled from the programs they sit among, consistent with pulsatile, cell-to-cell-variable A3 bursts.

### Figure 4 Panels

| Panel | Content | Source |
|-------|---------|--------|
| a | Three-population UMAP (SBS2-HIGH / CNV-HIGH / NORMAL) | `Step00B` + plotting |
| b | SBS2-HIGH activating chain (*RALY* five-gene + *HNRNPA2B1* / *CCL20*) | SBS2_VS_NORMAL network |
| c | CNV-HIGH inhibiting cornification chain | CNV_VS_NORMAL network |
| d | CNV-HIGH *RALY* activating module | CNV_VS_NORMAL network |
| Supp Fig 6 | Full networks clustered into gene groups, non-A3 community zooms | both networks |

### Narrative Arc

1. DE direction alone cannot separate a coordinated cofactor program from the broad shift of tumorigenesis, so differential co-expression networks are built, each tumor population against normal-adjacent (Fig. 4a).
2. The two networks recover 54 and 109 of the 174 catalogued A3 interactors; A3A and A3B share a gene group in SBS2-HIGH but separate in CNV-HIGH.
3. A *RALY*-anchored activating program is conserved across both tumor states, with *HNRNPA2B1* anchoring *CCL20* in SBS2-HIGH (Fig. 4b, 4d), nominating *RALY* and *HNRNPA2B1* as functional coactivators.
4. A cornification differentiation program forms the clearest inhibiting chain, lost as cells accumulate CNV (Fig. 4c).
5. The A3 wall shows the enzymes induced but decoupled from stable co-expression, consistent with pulsatile A3 induction.
6. **Motivates Question 3:** test whether the SBS2-HIGH program is a conserved cross-patient signal or driven by a few individuals, before reading the networks as biology.
