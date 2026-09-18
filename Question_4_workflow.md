# Question 4 Workflow: HPV16 Lifecycle (Figure 6)

**Maps to manuscript Results 4.3** — *HPV16 lifecycle stages are associated with distinct mutagenic programs in epithelial cells.*

> **Section numbering changed.** Earlier drafts carried this analysis as Section 4.4, and that label survives in several script docstrings and in the diagnostic's own banner text. In the submitted manuscript it is **Results 4.3**, and 4.4 is the patient-determinants section. Where a script prints "Section 4.4", read "Results 4.3".

> **Order relative to Question 3 is inverted.** The manuscript presents the lifecycle analysis (4.3) **before** the patient determinants (4.4), because the lifecycle stages are what the patient-level thresholds are built on. This walkthrough is numbered 4 and Question 3 is numbered 3, but Question 3 depends on the per-cell viral tables produced here and cannot run before it.

> **CNA versus CNV.** The manuscript writes copy-number alteration as CNA; the code uses `CNV` throughout.

## Question 4: How does HPV16 viral lifecycle stage organize the divergent mutagenic programs in A3-active epithelial cells?

### Rationale

Question 1 (Figs 1-3) resolved two epithelial states with distinct genomic damage, A3A-linked SBS2 point mutations and A3B-linked copy-number change. Question 2 (Fig 4, Fig 5) traced an activating co-expression chain linking A3A to immune and inflammatory genes (*CCL20*, *LCN2*, *SMOX* through the interactors *RALY* and *HNRNPA2B1*) that was decoupled in CNV-HIGH cells.

In the submitted manuscript, the motivation for turning to the virus comes from those two results directly: SBS2-HIGH cells up-regulate viral-recognition programs while CNA-HIGH cells down-regulate antigen processing and presentation (Fig. 4), and inflammatory genes including *LCN2* and *CCL20* sit inside the SBS2-HIGH co-expression network (Fig. 5).

> An earlier version of this walkthrough motivated the section through Question 3, on the grounds that the variants distinguishing high-contributor patients mapped to HPV16-associated genes. **That variant analysis is not in the submitted manuscript**, and the section order has since reversed, so that framing no longer applies.

The lifecycle stages of high-risk HPV are well characterized in stratified epithelium (Schichl and Doorbar 2025, ref 52; Doorbar et al. 2012, ref 53), but the link between lifecycle phase and APOBEC mutational outcome has not been resolved at single-cell resolution. The central finding is that SBS2-HIGH cells carry HPV16 in maintenance while CNV-HIGH cells carry it in productive infection, and that this lifecycle position organizes the host response: which enzyme is active (A3A vs A3B), which class of damage accumulates (SBS2 vs CNA), and whether the cell is immune-visible or replicating virus. Viral gene usage is expressed as a fraction of total HPV16 reads per cell to control for the roughly 2.6-fold higher viral load in CNV-HIGH cells.

**Note:** The neoantigen analysis (formerly Phase 5B in this directory) is now a standalone pipeline in `scripts/NEOANTIGEN/`; see `Question_5_workflow.md` (Figure 7).

---

## Directory Structure

```
scripts/HPV_ANALYSIS/
├── RUN_HPV_ANALYSIS.sh                                # SLURM master runner (Phases 1-4)
├── sample_list.txt                                    # 31 GSE173468 sample accessions
│
├── Phase1_HPV16_Data_Inventory.py                     # Phase 1: per-cell HPV16 consolidation
├── Phase1_HPV16_Data_Inventory_v2.py                  # Phase 1 (v2): revised mappings
├── Phase1_5_Index_Diagnostic.py                       # Phase 1.5: barcode index alignment check
├── Phase2_Raw_HPV16_Counts.py                         # Phase 2: raw UMI recovery from Kraken2
├── Phase3_HPV16_Populations_and_Genome.py             # Phase 3: L-method threshold + genome probe
├── Phase4_HPV16_Genome_Alignment.py                   # Phase 4: minimap2 alignment, per-cell gene counts
│
├── Phase5A_Population_Consolidation.py                # Phase 5A: old two-population model (superseded)
├── Phase5A_Revised_Population_Discovery.py            # Phase 5A: data-driven discovery (superseded)
├── Phase5A_v2_DEG_Analysis.py                         # Phase 5A: DE/GSEA on old model (superseded)
│
├── Diagnostic_Figure6_HostMarkers_and_IntegrationProxy.py  # Results 4.3 diagnostic + number audit
├── Diagnostic_Figure6_PanelC_Marker_Audit.py          # Panel C tier membership audit
├── Diagnose_NORMAL_HPV16_ORF_reads.py                 # NORMAL off-ORF read check
├── Generate_Figure6_Lifecycle_Panels.py               # Figure 6 panel generation (current)
│
├── parse_chimeric_junctions.py                        # STAR chimeric junction parser (used by NEOANTIGEN)
└── TROUBLESHOOTING/                                   # Diagnostics and pre-rerun figure scripts
```

> **`TROUBLESHOOTING/` here is load-bearing for a different section.** Five scripts in it produce most of the numbers in Results **4.4**: `Diagnostic_Patient_HPV16_Load.py`, `Diagnostic_Patient_Lifecycle_Phase.py`, `Diagnostic_Patient_Determinants_Table.py`, `Diagnostic_CNV_High_Patient_Contributions.py`, plus the shared `contribution.py` and a second, diverged `patient_config.py`. They write into `data/FIG_5/00_diagnostics/`. See the Question 3 walkthrough.
>
> `TROUBLESHOOTING/` also holds `Generate_Supp_CellCycle.py` (**Supplemental Figure 6**), `Generate_Supp_Contribution_Virus.py` (**Supp. Fig. 7a-b**) and `Generate_Supp_Enzyme_Conjunction.py` (**Supp. Fig. 9**). Three supplemental figures in the paper are produced from a directory the walkthrough convention excludes. Worth promoting these to the main directory.

---

## Input Data

| Input | Path | Description |
|-------|------|-------------|
| AnnData (final) | `data/FIG_4/00_input/adata_final.h5ad` | 155,650 cells, 27,736 genes (ClusterCatcher output) |
| AnnData (viral) | `data/FIG_6/00_input/adata_v_pp.h5ad` | Kraken2 viral detection object |
| Signature weights | `data/FIG_4/00_input/signature_weights_per_cell.txt` | Per-cell SBS weights, requires `.T` on load |
| Three-group assignments | `data/FIG_4/01_group_selection/three_group_assignments.tsv` | 546 SBS2-HIGH, 546 CNV-HIGH, 546 NORMAL |
| Kraken2 matrices | `SC/fastq/.../kraken2_filtered_feature_bc_matrix/` | Raw viral UMI counts per sample |
| Unmapped BAMs / FASTQs | `SC/fastq/.../possorted_genome_bam_unmapped.*` | Unmapped reads with CB tags, for viral alignment |
| HPV16 reference | `data/FIG_6/03_hpv16_genome/HPV16_NC_001526.4.fa` | NC_001526.4, 7,906 bp (downloaded by Phase 4) |

---

## Pipeline Overview

Phases 1-4 extract and process per-cell viral data and are population-independent. The Phase 5A population model is superseded by the Question 2 three-group selection. The diagnostic and figure scripts combine Phase 4 viral gene counts with those populations.

```
Phase 1 -> Phase 1.5 -> Phase 2 -> Phase 3 -> Phase 4
(inventory) (index dx)  (raw UMI) (threshold) (genome alignment)
                                                     │
        three_group_assignments.tsv (Question 2) ────┤
                                                     ▼
   Diagnostic_Figure6_HostMarkers_and_IntegrationProxy.py
                                                     ▼
            Generate_Figure6_Lifecycle_Panels.py
```

Phases 1-4 run via SLURM (`sbatch RUN_HPV_ANALYSIS.sh`); the diagnostic and figure generation run interactively in the `NETWORK` environment.

---

## Script Descriptions

### Phases 1-2: Inventory and Raw Counts

**`Phase1_HPV16_Data_Inventory.py`** (and `_v2`) consolidate per-cell measurements into a master table: HPV16 status, SBS2 weights, CNA scores, stemness, and annotations, with cross-tabulations of HPV16 status against population and cluster. **`Phase1_5_Index_Diagnostic.py`** verifies barcode mapping integrity between the Kraken2 indices and the main object before raw extraction.

**`Phase2_Raw_HPV16_Counts.py`** recovers true integer UMI counts from the Kraken2 per-sample sparse matrices, since the normalized counts in `adata_v_pp.h5ad` saturate and read as effectively binary. It identifies the HPV16 feature (taxonomy 333760) per sample and maps counts back to the integrated object, writing `raw_HPV16`.

### Phase 3: Threshold

**`Phase3_HPV16_Populations_and_Genome.py`** applies L-method piecewise linear regression to the sorted raw HPV16 count distribution, placing the breakpoint at 8 UMIs and splitting cells into HPV16-negative (0), HPV16-ambiguous (1-7), and HPV16-positive (8+). It also probes for the unmapped-read files needed by Phase 4. The old HPV-tier-by-cluster population definitions here are superseded; the threshold and genome-probe outputs remain in use.

This is also the source of the three numbers the Figure 6 diagnostic marks **out of scope**: the 94.6% basal restriction, the tier counts, and the Fisher OR of 1.01.

### Phase 4: Viral Genome Alignment

**`Phase4_HPV16_Genome_Alignment.py`** is the core alignment step. It downloads the HPV16 reference (NC_001526.4, 7,906 bp), builds a minimap2 index, aligns unmapped reads from all samples to the viral genome, recovers cell barcodes, and maps positions to HPV16 gene regions: E6 (83-559), E7 (562-858), E1 (865-2813), E2 (2756-3852), E4 (3332-3619), E5 (3849-4100), L2 (4236-5657), L1 (5560-7155), and URR (7156-7906 and 1-82). Its `per_cell_hpv16_gene_counts.tsv` (76,978 cells) feeds all downstream lifecycle analysis.

> **E4 reads are zero in every cell.** This is expected and documented in Methods 6.4: the spliced E1^E4 transcript does not map to the standalone E4 coordinates, so the amplification signal reflects E5 alone. Do not read the zero as a pipeline failure.

### Phase 5A (Superseded)

**`Phase5A_Population_Consolidation.py`**, **`Phase5A_Revised_Population_Discovery.py`**, and **`Phase5A_v2_DEG_Analysis.py`** defined and profiled the original two-population model (Pop1_Mutagenic / Pop2_Stealth). They are superseded by the three-group selection and retained for provenance. Their column names (`two_pop`, `lifecycle_stage`, `refined_group` with values like `Stealth_CNV` and `late_dominant`) still ship inside `population_assignments.tsv` and the FIG_6 profile tables, and use terminology the paper has retired ("productive", not "late"). Worth a cleanup pass before the repo gets traffic.

### Results 4.3 Diagnostic

**`Diagnostic_Figure6_HostMarkers_and_IntegrationProxy.py`** is the primary diagnostic and the text-number audit harness. On the gated HPV16-positive set (`raw_HPV16 >= 8` and total genome reads > 0, giving **197 / 446 / 8** cells) it reproduces the figure's lifecycle-fraction computation exactly (per-cell gene fractions, 10,000-permutation test, BH-FDR within the 8-gene and 4-phase families), reports the URR / ORF / intergenic breakdown both pooled and per-cell-mean, computes integration-proxy metrics, profiles the host-marker panel on the ungated 546/546/546 groups, and diffs every Results 4.3 number against freshly computed values, printing MATCH / DIFF / OUT-OF-SCOPE per claim.

**How to run it (two passes, not one).** The script is assertion-based against a hardcoded `CLAIMS` list and emits a corrected block at the end:

1. Run once. Expect DIFF rows for any claim whose stored value is stale.
2. Copy the emitted block from `DIAGNOSTIC_LIFECYCLE_MARKERS/emitted_claims_block.py` over the `CLAIMS` list in the script.
3. Re-run and confirm ALL MATCH. Only then are the values safe to quote.

**The BH family moved from 57 genes to 59.** v4 locked the family at 57 to match the panel; v5 uses 59. Every host-gene q shifted by roughly 3.5% as a result, and the September 2026 run confirmed the **manuscript carries the 59-gene values**. Any q-value from a pre-v5 run of this script is superseded. `Generate_Figure6_Lifecycle_Panels.py` must carry the same 8-tier / 59-gene `DOTPLOT_CATEGORIES`, or the panel's significance marks and the paragraph's q-values are computed against different denominators.

> The emitted `CLAIMS` block also contains five `brd4_*` entries. Those belong to the BRD4 isoform analysis, which is **reply-letter material for Cheng-Ming Chiang and is not in the manuscript**. Comment them so a future run does not send someone hunting for a manuscript sentence that does not exist.

**`Diagnostic_Figure6_PanelC_Marker_Audit.py`** audits which genes sit in which Panel C tier. A copy also exists in `TROUBLESHOOTING/`; the main-directory one is current.

**`Diagnose_NORMAL_HPV16_ORF_reads.py`** is a text-only check confirming how many NORMAL HPV16-positive cells have zero ORF reads, whether off-ORF reads are a NORMAL-only artifact or universal, and whether the raw Cell Ranger feature and the Phase 4 genome-aligned total agree.

### Figure Generation

**`Generate_Figure6_Lifecycle_Panels.py`** produces the Figure 6 panels from the Phase 4 gene counts crossed with the three-group populations.

**Panel names in the script do not match the manuscript panel letters.** The mapping:

| Script name | Manuscript panel | Content |
|---|---|---|
| — | Fig. 6a | UMAP of epithelial cells by HPV16 load |
| Panel B family (18 tests) | Fig. 6b | HPV16 UMI violins per population |
| — | Fig. 6c | HPV16-positive cell counts, Fisher OR = 1.01, p = 0.91 |
| **Panel F** | **Fig. 6d** | lifecycle gene and phase fractions; source of the fraction numbers |
| Panel B family | Fig. 6e | SBS2, CNA, stemness, A3A, A3B violins |
| **Panel C** | **Fig. 6f** | host-marker dot plot, 59 genes in 8 tiers |

Panel F gates to HPV16-positive cells and runs the same permutation and BH-FDR scheme as the diagnostic. The A3A and A3B q-values sit in the **Panel B family of 18**, a separate correction family from the 59-gene Panel C family, which is why the diagnostic lists them as outside the Panel C audit.

### Utility

**`parse_chimeric_junctions.py`** parses STAR chimeric junction output and is used by `scripts/NEOANTIGEN/Step04b_Fusion_Analysis.py`. (An earlier version of this document named `Step05_Fusion_Analysis.py`, which does not exist; the neoantigen pipeline was renumbered.) **`sample_list.txt`** holds the 31 GSE173468 accessions. **`RUN_HPV_ANALYSIS.sh`** runs Phases 1-4.

---

## Key Results

All values below are from the 59-gene BH family and were re-verified against on-disk output in September 2026. **Every host-gene q-value in the previous version of this document came from the retired 57-gene family and is superseded.**

### Viral

| Finding | Evidence |
|---------|----------|
| HPV16 reads are basal-restricted | 94.6% of HPV16-positive cells are basal epithelial |
| HPV16 tiers (L-method, 8 UMIs) | negative 22,153 (42.5%), ambiguous 14,046 (26.9%), positive 15,927 (30.6%) |
| HPV16 presence alone does not predict SBS2-HIGH | Fisher OR = 1.01, p = 0.91 |
| Gated HPV16-positive set | 197 SBS2-HIGH / 446 CNV-HIGH / 8 NORMAL |
| URR dominates reads (pooled) | 63.5 / 63.5 / 64.9% |
| CNV-HIGH carries higher viral load | 235.1 vs 90.1 reads per positive cell, 2.6-fold, q = 1.7e-14 |
| SBS2-HIGH is in maintenance | 25.9% vs 13.2%, stage q = 1.3e-4, driven by E1 (q = 2.0e-4) |
| CNV-HIGH is in productive infection | capsid 17.8% vs 10.3%, stage q = 1.3e-4; L1 and L2 each q = 2.0e-4; E5 q = 2.0e-4 |
| E6/E7 dosage does not differ | 0.53% and 0.79%, phase q = 0.10 |
| E2 intact, episomal in both | q = 0.11, no E2 loss |

> **The two tumor URR values are genuinely different numbers.** Pooled values are 63.535% and 63.4928%, which both round to 63.5. The per-cell-mean versions differ visibly (63.3 and 67.1), so the prose must cite **pooled**, which it does. Worth knowing, because the coincidence looks like a copy error.

### Host

Of 59 host genes in 8 tiers, **47 differ significantly** between SBS2-HIGH and CNV-HIGH.

| Finding | Evidence (59-gene family) |
|---------|---------|
| Immune visibility, SBS2-HIGH | B2M q = 9.0e-100, HLA-A q = 7.5e-42, HLA-B q = 1.2e-15, HLA-C q = 5.9e-8, TAP1 q = 1.0e-4 |
| Interferon effectors, SBS2-HIGH | 7 of 8 significant; ISG15 is the exception (q = 0.20) |
| Upstream sensing does not separate | STAT1 q = 0.31, IRF1 q = 0.62, STAT2 q = 0.13, DDX58 q = 0.09 |
| Differentiation, SBS2-HIGH | *IVL* 2.68 vs 0.09, q = 5.2e-70 |
| DDR, CNV-HIGH | *CHEK2* q = 1.4e-5, *BRCA1* q = 1.1e-8, *NBN* q = 1.3e-3, *H2AX* q = 4.3e-11 |
| ATR arm, CNV-HIGH | *TOPBP1* q = 8.4e-14, *CHEK1* q = 7.2e-11 |
| BET proteins, CNV-HIGH | *BRD2* q = 1.6e-9, *BRD3* q = 1.8e-25, *BRD4* q = 1.3e-11 |
| G2/M arrest, CNV-HIGH | *CDC25A* q = 2.4e-16, *CDC25C* q = 3.0e-12, *CDK1* q = 4.7e-18, *CCNB1* q = 9.8e-27 |
| Cell-cycle re-entry, CNV-HIGH | 8 of 9 significant; *MCM7* q = 1.7e-47, *TOP2A* q = 1.2e-25, *MKI67* q = 2.4e-19 |
| Proliferative basal identity, CNV-HIGH | *KRT14* q = 8.9e-57, *KRT5* q = 4.4e-7 |
| A3A dominates SBS2-HIGH | 6.46 vs 2.08 (Panel B family, q = 2.66e-135) |
| A3B dominates CNV-HIGH | 4.95 vs 2.21 (Panel B family, q = 8.06e-67) |
| Cell-cycle phase follows the split | SBS2-HIGH 63.9% G1; CNV-HIGH 81.7% S or G2/M, 41.8% G2/M (Supp. Fig. 6a) |

Three tier-membership choices worth knowing, since they are not obvious from the gene names. *STAT5A* and *STAT5B* sit in the DDR tier because STAT-5 drives TopBP1 transcription and also regulates the ATM response in HPV-positive keratinocytes (ref 49). *CASP7* is there because caspase cleavage of E1 is required for genome amplification (ref 50). *NSD2* is there as a BRD4-L-recruited repair factor (ref 51). All three are stated in the Figure 6 legend.

### Biological Model

A cell's position in the HPV16 lifecycle sets both which A3 enzyme is active and the class of damage that accumulates. In E2-regulated maintenance, SBS2-HIGH cells are immune-visible and differentiating, with A3A dominant and point mutations accumulating as SBS2. In productive infection, CNV-HIGH cells run a damage response spanning both the ATM and ATR arms with G2/M arrest, A3B dominant, and chromosomal instability accumulating. The productive program runs here in cells that keep a proliferative basal keratin identity (*KRT5*, *KRT14* high) without terminal differentiation (*IVL* near-absent), so it is uncoupled from keratinocyte differentiation. SBS2-HIGH and CNV-HIGH are therefore two states along a maintenance-to-productive continuum, and the SBS2:CNA ratio offers a molecular estimate of a tumor's position along that axis.

---

## Output Structure

```
data/FIG_6/
├── 00_input/                       adata_v_pp.h5ad
├── 00_diagnostic_inventory/        Phase 1 master table + report
├── 01_raw_hpv16_counts/            Phase 2 raw UMI counts
├── 02_populations/                 population_assignments.tsv (all 52,126 basal cells,
│                                   with raw_HPV16 + n_counts; used by the Q3 load analysis)
├── 03_hpv16_genome/                Phase 4: per_cell_hpv16_gene_counts.tsv (PRIMARY, 76,978 cells)
│                                   hpv16_gene_by_population.tsv (29,921 cells)
├── 04_population_profiles/         two-population profiles (superseded terminology)
├── DIAGNOSTIC_LIFECYCLE_MARKERS/   diagnostic fractions, integration proxy, host markers,
│                                   emitted_claims_block.py
├── DIAGNOSTIC_BRD4_ISOFORM/        BRD4-L vs BRD4-S ratio (reply-letter only)
└── FIGURE_6_PANELS/                Figure 6 panels + composite
```

---

## Dependencies and Relationships

`NETWORK` conda environment (minimap2 2.28, pysam, scanpy, gseapy). Question 2 supplies the three-group assignments and the AnnData. This question establishes the lifecycle framework; the immune-visible versus immune-evasive split it defines feeds the neoantigen tiers of Question 5, and its per-cell viral tables feed the patient determinants of Question 3.

---

## Open items

- **Confirm `DOTPLOT_CATEGORIES` in the figure script matches the 59-gene / 8-tier structure.** If it does not, Panel C's significance marks and the Results 4.3 q-values are computed on different denominators.
- **Comment or remove the five `brd4_*` entries** from the diagnostic's `CLAIMS` list; they are reply-letter values, not manuscript claims.
- **Promote the three supplemental-figure generators** out of `TROUBLESHOOTING/`: `Generate_Supp_CellCycle.py` (Supp. Fig. 6), `Generate_Supp_Contribution_Virus.py` (Supp. Fig. 7a-b), `Generate_Supp_Enzyme_Conjunction.py` (Supp. Fig. 9).
- **Retired column names still ship** in `population_assignments.tsv` and the FIG_6 profile tables: `two_pop` (`Pop1_Mutagenic` / `Pop2_Stealth`), `lifecycle_stage` (`early_dominant` / `late_dominant`), `refined_group` (`Stealth_CNV` / `OVERLAP`), `pop1_score`, and four `axis_*` flags. A reader will find "late_dominant" and "Stealth_CNV" and wonder why they do not match the paper.
- **`Diagnostic_Figure6_PanelC_Marker_Audit.py` is duplicated** in the main directory and in `TROUBLESHOOTING/`. Diff and keep one.
- **Two independent tests of the E6/E7 result disagree on significance.** The Panel F phase permutation gives q = 0.10 and the Mann-Whitney on `E6E7_frac_of_total` gives BH q = 5.55e-5. The diagnostic prints the conflict and recommends citing the effect size (under 1% of viral reads in both groups) rather than "no difference." The manuscript does this correctly; noted so the discrepancy is not rediscovered.
