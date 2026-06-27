# Question 4 Workflow: HPV16 Lifecycle (Figure 6)

## Question 4: How does HPV16 viral lifecycle stage organize the divergent mutagenic programs in A3-active epithelial cells?

### Rationale

Question 1 (Figs 1-3) resolved two epithelial states with distinct genomic damage, A3A-linked SBS2 point mutations and A3B-linked copy-number change, and Question 2 (Fig 4) traced an activating co-expression chain linking A3A to immune and inflammatory genes (*CCL20*, *LCN2*, *SMOX* through the interactors *RALY* and *HNRNPA2B1*) that was decoupled in CNV-HIGH cells. Question 3 (Fig 5) then traced the variants distinguishing high-contributor patients to HPV16-associated genes rather than the A3 program. Those threads point at the virus, so this question asks directly how HPV16 infection context shapes the divergence.

The lifecycle stages of high-risk HPV are well characterized in stratified epithelium (Doorbar et al. 2012; Moody and Laimins 2010; McBride and Warburton 2017), but the link between lifecycle phase and APOBEC mutational outcome has not been resolved at single-cell resolution. The central finding is that SBS2-HIGH cells carry HPV16 in maintenance while CNV-HIGH cells carry it in productive infection, and that this lifecycle position organizes the host response: which enzyme is active (A3A vs A3B), which class of damage accumulates (SBS2 vs CNV), and whether the cell is immune-visible or replicating virus. Viral gene usage is expressed as a fraction of total HPV16 reads per cell to control for the roughly 2.6-fold higher viral load in CNV-HIGH cells.

**Note:** The neoantigen analysis (formerly Phase 5B in this directory) is now a standalone six-step pipeline in `scripts/NEOANTIGEN/`; see `Question_5_workflow.md` (Figure 7).

---

## Directory Structure

```
scripts/HPV_ANALYSIS/
├── RUN_HPV_ANALYSIS.sh                                # SLURM master runner (Phases 1-4)
├── sample_list.txt                                   # 31 GSE173468 sample accessions
│
├── Phase1_HPV16_Data_Inventory.py                    # Phase 1: per-cell HPV16 consolidation
├── Phase1_HPV16_Data_Inventory_v2.py                 # Phase 1 (v2): revised mappings
├── Phase1_5_Index_Diagnostic.py                      # Phase 1.5: barcode index alignment check
├── Phase2_Raw_HPV16_Counts.py                        # Phase 2: raw UMI recovery from Kraken2
├── Phase3_HPV16_Populations_and_Genome.py            # Phase 3: L-method threshold + genome probe
├── Phase4_HPV16_Genome_Alignment.py                  # Phase 4: minimap2 alignment, per-cell gene counts
│
├── Phase5A_Population_Consolidation.py               # Phase 5A: old two-population model (superseded)
├── Phase5A_Revised_Population_Discovery.py           # Phase 5A: data-driven discovery (superseded)
├── Phase5A_v2_DEG_Analysis.py                        # Phase 5A: DE/GSEA on old model (superseded)
│
├── Diagnostic_Figure6_HostMarkers_and_IntegrationProxy.py  # Section 4.4 diagnostic + text-number audit
├── Diagnose_NORMAL_HPV16_ORF_reads.py                # NORMAL off-ORF read check, raw vs genome agreement
├── Generate_Figure6_Lifecycle_Panels.py             # Figure 6 panel generation (current)
│
├── parse_chimeric_junctions.py                       # STAR chimeric junction parser (used by NEOANTIGEN)
└── TROUBLESHOOTING/                                  # Diagnostics and pre-rerun figure scripts (not documented here)
```

---

## Input Data

| Input | Path | Description |
|-------|------|-------------|
| AnnData (final) | `data/FIG_4/00_input/adata_final.h5ad` | 155,650 cells, 27,736 genes (ClusterCatcher output) |
| AnnData (viral) | `data/FIG_6/00_input/adata_v_pp.h5ad` | Kraken2 viral detection object |
| Signature weights | `data/FIG_4/00_input/signature_weights_per_cell.txt` | Per-cell SBS signature weights |
| Three-group assignments | `data/FIG_4/01_group_selection/three_group_assignments.tsv` | 546 SBS2-HIGH, 546 CNV-HIGH, 546 NORMAL |
| Kraken2 matrices | `SC/fastq/.../kraken2_filtered_feature_bc_matrix/` | Raw viral UMI counts per sample |
| Unmapped BAMs / FASTQs | `SC/fastq/.../possorted_genome_bam_unmapped.*` | Unmapped reads with CB tags, for viral alignment |
| HPV16 reference | `data/FIG_6/03_hpv16_genome/HPV16_NC_001526.4.fa` | NC_001526.4, 7,906 bp (downloaded by Phase 4) |

---

## Pipeline Overview

Phases 1-4 extract and process per-cell viral data and are population-independent. The Phase 5A population model is superseded by the Figure 4 three-group selection. The diagnostic and figure scripts combine Phase 4 viral gene counts with the Figure 4 populations.

```
Phase 1 -> Phase 1.5 -> Phase 2 -> Phase 3 -> Phase 4
(inventory) (index dx)  (raw UMI) (threshold) (genome alignment)
                                                     │
        three_group_assignments.tsv (Figure 4) ──────┤
                                                     ▼
   Diagnostic_Figure6_HostMarkers_and_IntegrationProxy.py
                                                     ▼
            Generate_Figure6_Lifecycle_Panels.py
```

Phases 1-4 run via SLURM (`sbatch RUN_HPV_ANALYSIS.sh`); the diagnostic and figure generation run interactively in the `NETWORK` environment.

---

## Script Descriptions

### Phases 1-2: Inventory and Raw Counts

**`Phase1_HPV16_Data_Inventory.py`** (and `_v2`) consolidate per-cell measurements into a master table: HPV16 status, SBS2 weights, CNV scores, stemness, and annotations, with cross-tabulations of HPV16 status against population and cluster. **`Phase1_5_Index_Diagnostic.py`** verifies barcode mapping integrity between the Kraken2 indices and the main object before raw extraction.

**`Phase2_Raw_HPV16_Counts.py`** recovers true integer UMI counts from the Kraken2 per-sample sparse matrices, since the normalized counts in `adata_v_pp.h5ad` saturate and read as effectively binary. It identifies the HPV16 feature (taxonomy 333760) per sample and maps counts back to the integrated object.

### Phase 3: Threshold

**`Phase3_HPV16_Populations_and_Genome.py`** applies L-method piecewise linear regression to the sorted raw HPV16 count distribution, placing the breakpoint at 8 UMIs and splitting cells into HPV16-negative (0), HPV16-ambiguous (1-7), and HPV16-positive (8+). It also probes for the unmapped-read files needed by Phase 4. The old HPV-tier-by-cluster population definitions here are superseded by the Figure 4 three-group selection; the threshold and genome-probe outputs remain in use.

### Phase 4: Viral Genome Alignment

**`Phase4_HPV16_Genome_Alignment.py`** is the core alignment step. It downloads the HPV16 reference (NC_001526.4, 7,906 bp), builds a minimap2 index, aligns unmapped reads from all samples to the viral genome, recovers cell barcodes, and maps positions to HPV16 gene regions: E6 (83-559), E7 (562-858), E1 (865-2813), E2 (2756-3852), E4 (3332-3619), E5 (3849-4100), L2 (4236-5657), L1 (5560-7155), and URR (7156-7906 and 1-82). Its `per_cell_hpv16_gene_counts.tsv` (76,978 cells) feeds all downstream lifecycle analysis.

### Phase 5A (Superseded)

**`Phase5A_Population_Consolidation.py`**, **`Phase5A_Revised_Population_Discovery.py`**, and **`Phase5A_v2_DEG_Analysis.py`** defined and profiled the original two-population model (Pop1_Mutagenic / Pop2_Stealth). They are superseded by the Figure 4 three-group selection and retained for provenance; their DE findings (spliceosome enrichment in the CNV side, immune pathways in the SBS2 side) remain consistent with the current analysis.

### Section 4.4 Diagnostic

**`Diagnostic_Figure6_HostMarkers_and_IntegrationProxy.py`** is the primary diagnostic and the text-number audit harness for Section 4.4. On the gated HPV16-positive set (raw HPV16 >= 8 and total > 0) it reproduces the figure's lifecycle-fraction computation exactly (per-cell gene fractions, 10,000-permutation test, BH-FDR within the 8-gene and 4-phase families), reports the URR / ORF / intergenic breakdown both pooled and per-cell-mean (confirming the prose should cite the pooled two-thirds), computes the integration-proxy metrics, profiles the host-marker panel on the ungated 546/546/546 groups, and then diffs every Section 4.4 number against freshly computed values and prints MATCH / DIFF / OUT-OF-SCOPE per claim. Para-1 numbers (basal enrichment, tier counts, Fisher) are marked out of scope with their correct source (the Phase 3 L-method).

**`Diagnose_NORMAL_HPV16_ORF_reads.py`** is a text-only check confirming how many NORMAL HPV16-positive cells have zero ORF reads, whether off-ORF reads are a NORMAL-only artifact or universal, and whether the raw Cell Ranger feature and the Phase 4 genome-aligned total agree.

### Figure Generation

**`Generate_Figure6_Lifecycle_Panels.py`** produces the Figure 6 panels from the Phase 4 gene counts crossed with the Figure 4 populations: a UMAP of epithelial cells colored by HPV16 load with tier indicators, the normalized lifecycle gene and phase fractions (Panel F, the source of the Section 4.4 fraction numbers), the host-marker dot plot contrasting the ATM/productive program against the immune/differentiation program, the per-population viral load, and the A3A versus A3B distributions. Panel F gates to HPV16-positive cells and runs the same permutation and BH-FDR scheme as the diagnostic.

### Utility

**`parse_chimeric_junctions.py`** parses STAR chimeric junction output and is used by `scripts/NEOANTIGEN/Step05_Fusion_Analysis.py`. **`sample_list.txt`** holds the 31 GSE173468 accessions. **`RUN_HPV_ANALYSIS.sh`** runs Phases 1-4.

---

## Key Results

| Finding | Evidence |
|---------|----------|
| HPV16 reads are basal-restricted | 94.6% of HPV16-positive cells are basal epithelial |
| HPV16 tiers (L-method, 8 UMIs) | negative 22,153 (42.5%), ambiguous 14,046 (26.9%), positive 15,927 (30.6%) |
| HPV16 presence alone does not predict SBS2-HIGH | Fisher OR = 1.01, p = 0.91 |
| URR dominates reads (pooled) | 63.5 / 63.5 / 64.9% in SBS2-HIGH / CNV-HIGH / NORMAL |
| CNV-HIGH carries higher viral load | 235.1 vs 90.1 total HPV16 reads per positive cell, 2.6-fold, q = 1.7e-14 |
| SBS2-HIGH is in maintenance | more maintenance reads, driven by E1 (q = 2.0e-4) |
| CNV-HIGH is in productive infection | more capsid L1 and L2 (both q = 2.0e-4); higher E5 (q = 2.0e-4) |
| E6/E7 dosage does not differ | under 1% of reads each, q = 0.10 |
| E2 intact, episomal in both | E2 comparable, q = 0.11, no E2 loss |
| Immune visibility, SBS2-HIGH | B2M q = 4.0e-100, HLA-A q = 5.8e-42, HLA-B q = 1.1e-15, TAP1 q = 1.2e-4 |
| Differentiation, SBS2-HIGH | *IVL* 2.68 vs 0.09, q = 3.1e-70 |
| ATM-dependent DDR, CNV-HIGH | *BRCA1* q = 1.2e-8, *H2AX* q = 4.6e-11, *CHEK2* q = 1.6e-5 |
| G2/M arrest, CNV-HIGH | *CDK1* q = 4.9e-18, *CCNB1* q = 9.9e-27, *CDC25C* q = 3.3e-12 |
| Proliferative basal identity, CNV-HIGH | *KRT14* q = 4.8e-57, *MKI67* q = 2.6e-19 |
| A3A dominates SBS2-HIGH | 6.46 vs 2.08, q = 5.5e-134 |
| A3B dominates CNV-HIGH | 4.95 vs 2.21, q = 8.3e-66 |

### Biological Model

A cell's position in the HPV16 lifecycle sets both which A3 enzyme is active and the class of damage that accumulates. In E2-regulated maintenance, SBS2-HIGH cells are immune-visible and differentiating, with A3A dominant and point mutations accumulating as SBS2. In productive infection, CNV-HIGH cells run an ATM-dependent damage response and G2/M arrest with A3B dominant and chromosomal instability accumulating. The productive program runs here in cells that keep a proliferative basal keratin identity (*KRT5*, *KRT14* high) without terminal differentiation (*IVL* near-absent), so it is uncoupled from keratinocyte differentiation. SBS2-HIGH and CNV-HIGH are therefore two states along a maintenance-to-productive continuum, and the SBS2:CNV ratio offers a molecular estimate of a tumor's age along that axis.

---

## Output Structure

```
data/FIG_6/
├── 00_input/                       adata_v_pp.h5ad
├── 00_diagnostic_inventory/        Phase 1 master table + report
├── 01_raw_hpv16_counts/            Phase 2 raw UMI counts
├── 03_hpv16_genome/                Phase 4: per_cell_hpv16_gene_counts.tsv (PRIMARY, 76,978 cells)
├── DIAGNOSTIC_LIFECYCLE_MARKERS/   diagnostic fractions, integration proxy, host markers
└── FIGURE_6_PANELS/                Figure 6 panels + composite
```

---

## Dependencies and Relationships

`NETWORK` conda environment (minimap2 2.28, pysam, scanpy, gseapy). Figure 4 supplies the three-group assignments and adata. Figure 6 establishes the lifecycle framework, and the immune-visible versus immune-evasive split it defines feeds directly into the neoantigen tiers of Question 5 (Figure 7).
