# Question 2 Workflow

**Maps to manuscript Results 4.2** — *Single-cell network analysis identifies co-expression programs that may activate A3 mutagenic activity.*

> **Two naming traps in this section.** First, the manuscript writes copy-number alteration as **CNA**; the code predates that decision and uses `CNV` throughout (`cnv_score`, `CNV_HIGH`, `NETWORK_CNV_VS_NORMAL`). Second, and more confusing, **`scripts/NETWORK_SINGLE_CELL/` and `data/FIG_4/` produce manuscript Figure 5**, not Figure 4. Manuscript Figure 4 is the differential-expression and KEGG panel. The `FIG_4` and `Generate_Figure4_*` names were fixed before the figure order settled and were left alone so that existing output paths and logs stay valid. Panel and directory names below are given as they appear on disk; manuscript figure numbers are stated explicitly wherever they differ.

## Question 2: Does single-cell differential co-expression network analysis resolve A3 cofactor biology when each tumor population is compared to normal tissue?

### Rationale

Question 1 resolved two divergent programs in the basal compartment, A3A-linked SBS2 point mutations and A3B-linked copy-number change. Differential expression direction alone cannot take that further: a gene moving in the same direction as tumorigenesis could belong to a coordinated A3 cofactor program or could just be part of the broad transcriptional shift of becoming a tumor, and direction says nothing about how a gene tracks A3A or A3B. To assign functional roles to the known A3 interactors, we build differential co-expression networks comparing each tumor population to normal-adjacent basal cells and look for chains of genes that gain co-expression with the enzymes in tumor (candidate activators) or lose it (candidate brakes).

The lens here is tumor versus normal. An earlier SBS2-HIGH versus CNA-HIGH comparison was explored but dropped when the paper refocused on the tumor-to-normal contrast, so this question runs two networks: SBS2-HIGH vs NORMAL and CNV-HIGH vs NORMAL. The `NETWORK_SBS2_VS_CNV/` directory still exists under `data/FIG_4/` and is **not** used by the manuscript.

### Population design

Three basal populations of 546 cells each are defined in `Step00B`, and the two networks each compare a tumor population to NORMAL:

- **SBS2-HIGH (n=546):** from SBS2 > 0 cancer-tissue basal cells, by a composite of SBS2 weight (40%), low CNA (20%), low CytoTRACE2 stemness (20%), and high A3A/(A3A+A3B) fraction (20%).
- **CNV-HIGH (n=546):** from SBS2 == 0 cancer-tissue basal cells, by a mirror composite of high CNA (40%), total A3 (A3A+A3B) range-matched to the SBS2-HIGH mean by Gaussian proximity (20%), high stemness (20%), and high A3B/(A3A+A3B) fraction (20%). After the CNV-HIGH reselection the HPV16 late-gene and anti-correlated-profile terms were removed, so selection carries no viral input and is not circular with the Question 4 virus analysis.
- **NORMAL (n=546):** random sample of 546 of the 554 normal-adjacent basal cells.

Total A3 expression is range-matched between the two tumor groups so they differ in A3A/A3B composition rather than in summed enzyme level.

> **Two discrepancies between the code and Methods, both open.** The selection score uses `A3A/(A3A+A3B+0.01)`, with a pseudocount that Methods 6.3 does not mention. And the NORMAL draw is described as a random sample with no stated seed; `RANDOM_SEED = 42` exists in config but it is unconfirmed whether `Step00B` passes it. An unseeded draw would make the networks irreproducible from the public code.

### Data Sources

| Source | Description |
|--------|-------------|
| `adata_final.h5ad` | ClusterCatcher AnnData (155,650 cells, 27,736 genes) |
| `three_group_assignments.tsv` | SBS2-HIGH / CNV-HIGH / NORMAL, 546 each (1,638 cells) |
| `Harris_A3_interactors.txt` | 174 catalogued A3 interactors |
| `Harris_A3_interactors_A3B_only.txt` | A3B-specific interactor subset (51 genes) |
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
├── Step04_SC_Centrality_Metrics.py             # Centrality metrics
├── Compute_Node_Importance_Scores_SC.py        # Intra / inter community importance scoring
├── Diagnostic_A3_Interactor_Concordance.py     # Interactor concordance and chain enumeration
├── Diagnostic_Chain_Validation_SBS2_VS_CNV.py  # Chain behaviour across the two tumor populations
│
├── Extract_Harris_A3_Interactors.py            # Build the 174 A3 interactor list
├── Convert_Uniprot_to_Gene_Symbol.py           # UniProt -> gene symbol mapping for the interactor list
├── mmc15.xlsx                                  # AP-MS interactor source table
├── uniprot_accession_list.txt                  # Interactor-list intermediates
├── uniprot_accessions_to_convert.tsv
├── uniprot_to_gene_symbol_mapping.tsv
│
├── Generate_Figure4_Panels.py                  # Manuscript FIGURE 5 panels (plotting)
├── Generate_Figure4_Supplemental.py            # Manuscript Supplemental Figure 5 (plotting)
├── Step05_Generate_Figure4_Panels.py           # SUPERSEDED, see note
└── TROUBLESHOOTING/                            # Diagnostics and prior drafts (not documented here)
```

> **`Step05_Generate_Figure4_Panels.py` is not the current panel script.** It is a different file from `Generate_Figure4_Panels.py` (20,783 bytes versus 23,622) and was not run for the submitted figure. Only `Generate_Figure4_Panels.py` is called by `MAKE_FIGURES.sh`. Collapse to one on the next cleanup pass.

> **`TROUBLESHOOTING/` holds four complete copies of the pipeline** under `FIRST_DRAFT_SC_NETWORK/`, `SECOND_DRAFT_SC_NETWORK/`, `THIRD_DRAFT_SC_NETWORK/`, `FOURTH_DRAFT_SC_NETWORK/`, plus a `BACKUP/`. Every one contains a file named `Step03_SC_Community_Detection.py`. None is current. Run only from the main directory.

### Configuration (`network_config_SC.py`)

| Parameter | Value | Notes |
|-----------|-------|-------|
| DE selection | FDR < 0.05 | scanpy `rank_genes_groups`, internal BH-FDR |
| `FORCE_KEEP_A3` | False | A3A and A3B pass FDR naturally in both networks |
| DIFF threshold | auto | max fragmentation-rate criterion over a threshold sweep |
| `COMMUNITY_RESOLUTIONS` | `[0.2, 0.4, 0.6, 0.8, 1.0]` | **the resolution grid; see the note below** |
| Community detection | Leiden | full-network, resolution by modularity x ARI x evenness sweep, component-aware merge |
| `RUNS_PER_RESOLUTION` | 15 | with `COMMUNITY_BASE_SEED = 42` |
| `MIN_COMMUNITY_SIZE` | 10 | satellites preserved |
| `MIN_CELLS_DETECTED` | 10 | gene detection filter |

> **The module count depends on the resolution grid, so the grid belongs with the parameters.** An earlier build of this pipeline swept resolution in 0.1 steps and selected 0.70, giving 23 communities for the SBS2 network. The current grid has no 0.70, so 0.80 is the best available and gives **25**. The DIFF matrix is byte-identical between the two builds (`md5sum` on `SC_corr_DIFF.pkl` and `SC_diffexpr_stats.csv`), so this is a selection-grid change, not a data change and not run-to-run instability. The manuscript reports 25.

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
          (the concordance run enumerates the activating and inhibiting
           chains that the figure design and the manuscript narrative use)
```

**Stage B, plotting (`MAKE_FIGURES.sh`):** `Generate_Figure4_Panels.py` and `Generate_Figure4_Supplemental.py`.

The interactor list itself is built once by `Extract_Harris_A3_Interactors.py` and `Convert_Uniprot_to_Gene_Symbol.py` from `mmc15.xlsx`, yielding the 174-gene catalogue used throughout.

> **The reference matrix is still named `LOW`.** `Step02` writes `SC_corr_HIGH.pkl`, `SC_corr_LOW.pkl` and `SC_corr_DIFF.pkl`, where `LOW` is the NORMAL population, carried over from the retired two-group pipeline. `DIFF = HIGH − LOW` is therefore `HIGH − NORMAL`, as Methods 6.3 describes.

### Key Results

All values below were re-verified against on-disk output in September 2026 by two independent readers, `Diagnostic_A3_Interactor_Concordance.py` and `TROUBLESHOOTING/Harvest_Section4.2_numbers.py`, which agree on every field.

| Metric | SBS2-HIGH vs NORMAL | CNV-HIGH vs NORMAL |
|--------|---------------------|--------------------|
| DE genes at FDR < 0.05 | 3,877 | 6,802 |
| Network genes | 2,948 | 4,886 |
| Gene groups | **25** | 38 |
| Interactors recovered (of 174) | 54 | 109 |
| DIFF threshold | 0.40 | 0.45 |
| Leiden resolution | **0.80** | **0.80** |
| Modularity / ARI | 0.342 / 0.637 | 0.232 / 0.563 |
| Components / satellites | 18 / 17 | 34 / 33 |
| Fragmentation margin at the chosen threshold | +17 components | +26 components |
| A3A degree / community | 55 / C0 | 15 / C1 |
| A3B degree / community | 16 / C0 | 136 / C0 |
| A3 arrangement | A3A and A3B share a gene group | A3A and A3B separate |
| A3A DE (adj p, log2FC) | 3.66e-109, 7.73 | 3.00e-05, 1.22 |
| A3B DE (adj p, log2FC) | 1.64e-15, 2.58 | 2.80e-139, 6.69 |

The threshold selection reproduces independently. Component counts by threshold for SBS2 are 1, 1, **18**, 28, 29 across 0.30 to 0.50, so the largest jump is into 0.40; for CNV they are 2, 8, **34**, 39, 30 across 0.35 to 0.55, so the largest jump is into 0.45. Both match the max-fragmentation-rate criterion in Methods 6.3.

A3A and A3B both pass DE naturally in each network. A3A log2 fold change falls from 7.73 in SBS2-HIGH to 1.22 in CNV-HIGH, reflecting the sharp drop in A3A induction once cells enter the productive state, while A3B rises from 2.58 to 6.69.

**Activating chains (RALY-anchored, conserved in both):**

C0 of the SBS2 network contains **six** activating chains totalling 16 genes. The three that reach the manuscript are:

| Chain | Members | Interactor |
|---|---|---|
| act-0 | *RALY*, *LCN2*, *KRT24*, *LINC00278*, *CHMP4B*, *UTY* | *RALY* |
| act-3 | *HNRNPA2B1*, *CCL20* | *HNRNPA2B1* |
| act-1 | *RRAD*, *SMOX* | none |

The remaining three (act-2 *TMSB4X*/*EIF1AX*, act-4 *CDKN2A*/*TM4SF1*, act-5 *MRPL47*/*NTS*) are two-gene chains that the manuscript does not name.

- **CNV-HIGH:** *RALY* anchors the largest chain, a thirteen-gene module of translation and metabolic genes (*CPNE1*, *EIF6*, *CA2*, *DYNLRB1*, *ATOX1*, *C7orf50*, *DDAH2*, *ERGIC3*, *KLC3*, *MYL6*, *RHOD*, *ROMO1*) within the A3A group.
- The partners turn over between states but *RALY* holds the same role even though CNV-HIGH cells carry no SBS2 and A3A induction has fallen sharply, marking *RALY* as a tumor-conserved candidate coactivator rather than a tracker of enzyme level.

> **`CHMP4B` entered act-0 with the resolution change.** The chain was five genes under the earlier build and is six now. `patient_config.ACTIVATING_CHAIN_GENES` was updated from 9 to 10 genes in September 2026 and the LOPO scores rescored accordingly; see the Question 3 walkthrough.

> **A caveat worth having ready.** *RALY* sits at 20q11.2 and *CHMP4B* is also at 20q11, so their co-expression may reflect shared local regulation rather than a common functional program. 20q is a recurrent copy-gain region in HPV-associated cancers. The defence is specific: this pairing appears in the SBS2-HIGH network, and SBS2-HIGH cells are selected for **low** CNA, so a 20q amplification explanation is weak there. The same pairing in the CNV-HIGH network is more exposed. Separately, *LINC00278* and *UTY* are both Y-linked and will co-vary with patient sex in any mixed-sex cohort.

**Inhibiting chain (no catalogued interactor):**
- A **70-node** epithelial differentiation and cornification program in the A3A group of the CNV-HIGH network (*PRSS3*, *CLIC3*, *MAB21L4*, *SBSN*, *CYSRT1*, and the *SPRR* family), tightly co-expressed in normal tissue and losing coherence in tumor, with A3A at its edge through *SPRR1A*, *SPRR2D*, and *RAB11A*. This suggests loss of a normal differentiation program as cells accumulate CNA.
- The 70 counts A3A itself as a node; the chain-validation script reports the same object as 69 non-A3 genes plus A3A. Both are correct under their own convention.

**The A3 wall:**
- **Every** edge directly linking A3A or A3B to a first-degree neighbour is negative, in all four enzyme-by-network profiles. Zero positive edges out of 51, 10, 9 and 67 respectively. The enzymes are induced but co-expression-decoupled from the programs they sit among, consistent with pulsatile, cell-to-cell-variable A3 bursts.
- Group-level wall fractions: C0 of the SBS2 network is 85.7% wall edges; in the CNV network C0 (A3B) is 94.2% and C1 (A3A) is 84.5%.

> **The wall audit is module-scoped.** The concordance diagnostic profiles only the edges inside the A3-containing community, which covers 137 of the 222 total A3 edges across both networks. The largest gap is A3B in the CNV network, where 67 of 136 edges were profiled. Confirm against the full graph before quoting "every edge ... throughout both networks" in a response to a reviewer; if any unprofiled edge is positive, the claim needs scoping to the A3-containing modules.

**Chain behaviour across the two tumor populations (`Diagnostic_Chain_Validation_SBS2_VS_CNV.py`):**

This script asks whether chains defined against NORMAL also separate SBS2-HIGH from CNV-HIGH. **All three directional predictions fail.** Activating-chain internal DIFF comes back at mean −0.105 with 42.9% positive where positive was predicted; inhibiting-chain DIFF at +0.065 with 30.4% negative where negative was predicted; A3B edges 79 of 96 positive where negative was predicted.

This contradicts nothing in the manuscript, because every chain claim is scoped to a tumor-versus-NORMAL contrast and the SBS2_VS_CNV network was dropped. But it is a live reviewer question ("do your activating chains distinguish your two tumor populations?"), the honest answer is no, and it independently justifies dropping that third network. Worth having in the reply-letter folder rather than being surprised by it.

### Figure mapping

| Manuscript figure | Content | Source |
|---|---|---|
| **Fig. 4a** | KEGG pathway enrichment, up and down sets per comparison | see open item below |
| **Fig. 4b** | Differential expression, SBS2-HIGH and CNV-HIGH vs NORMAL | see open item below |
| **Fig. 5 (Center)** | Three-population UMAP (SBS2-HIGH / CNV-HIGH / NORMAL) | `Generate_Figure4_Panels.py`, `Panel_a_UMAP` |
| **Fig. 5 (Left)** | SBS2-HIGH concordant chain subnetwork | `Panel_b_*`, SBS2_VS_NORMAL, focus A3A, 2 hops |
| **Fig. 5 (Right)** | CNV-HIGH concordant chain subnetwork | `Panel_c_*`, CNV_VS_NORMAL, focus A3A, 2 hops |
| **Supp. Fig. 5** | Correlation and DIFF matrices, global networks, A3-community insets | `Generate_Figure4_Supplemental.py` |

The panel script also writes `Panel_d_CNV_A3B_chains` (CNV_VS_NORMAL focused on A3B, C0) and separate top/bottom activator and inhibitor variants of panels b and c. Not all are used in the composite.

> **Panels draw chains plus bridge edges, so node counts in the figure exceed the chain sizes in the text.** `Panel_b_top_SBS2_A3A_activator` reports 16 activating nodes and **22 bridge edges**, rendering as 40 nodes and 97 edges. That is why act-0 and act-5 appear joined in the left panel even though they are separate concordant components. The legend states this; a reader counting nodes off the figure will not otherwise match the six-gene chain in the text.

> **Supplemental Figure 5 shows a subset of communities.** The global-network panels drop communities of six or fewer genes, keeping 8 of 25 for SBS2 and 5 of 38 for CNV. The 54 and 109 interactor counts in the text are **network-wide** and include interactors in communities the panel never draws. Anyone counting interactors off the figure will come up short.

### Supplementary files

| File | Content |
|---|---|
| Supp. File 1 | SBS2-HIGH vs NORMAL differential expression table |
| Supp. File 2 | CNV-HIGH vs NORMAL differential expression table |
| Supp. Files 3 and 4 | DIFF correlation matrices, one per network |

> The DIFF matrices are dense gene-by-gene tables spanning the Step01 network-ready gene set, so they are large. `03_correlation_networks/edge_lists/SC_edges_DIFF.tsv` carries the same information above threshold in a fraction of the space and is the practical route if file size becomes a problem.

### Narrative Arc

1. DE direction alone cannot separate a coordinated cofactor program from the broad shift of tumorigenesis, so differential co-expression networks are built, each tumor population against normal-adjacent.
2. The two networks recover 54 and 109 of the 174 catalogued A3 interactors; A3A and A3B share a gene group in SBS2-HIGH but separate in CNV-HIGH.
3. Every edge from either enzyme to a first-degree neighbour is negative, so edge sign cannot be used to identify regulators, and the analysis extends to coherent chains beyond the first-degree neighbourhood.
4. A *RALY*-anchored activating program is conserved across both tumor states, with *HNRNPA2B1* anchoring *CCL20* in SBS2-HIGH, nominating *RALY* and *HNRNPA2B1* as candidate coactivators.
5. A cornification differentiation program forms the clearest inhibiting chain, lost as cells accumulate CNA.
6. The A3 wall shows the enzymes induced but decoupled from stable co-expression, consistent with pulsatile A3 induction.
7. **Motivates Question 3:** test whether the SBS2-HIGH program is a conserved cross-patient signal or driven by a few individuals, before reading the networks as biology.

