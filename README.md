# Network Architecture of APOBEC3-Driven Mutagenesis in HPV Head and Neck Cancers

## Overview

APOBEC3 (A3) cytidine deaminases are major drivers of somatic mutagenesis in head and neck squamous cell carcinoma (HNSCC), producing the characteristic SBS2 mutational signature. However, A3 expression alone fails to predict mutational burden: approximately 51% of A3-high tumors accumulate little SBS2. This paper uses differential co-expression network analysis across bulk TCGA and single-cell RNA-seq data to identify the transcriptional cofactors, HPV lifecycle context, and neoantigen consequences of APOBEC3-driven mutagenesis.

The central finding is that the HPV16 viral lifecycle organizes the entire host mutagenic program. During E2-regulated maintenance, interferon-driven A3A mutagenesis generates SBS2 point mutations through a co-expression network anchored by the RNA-binding proteins RALY and HNRNPA2B1. As cells transition to productive infection, ATM-dependent replication stress shifts mutagenesis to A3B-driven chromosomal instability. This lifecycle framework explains why A3 expression is necessary but not sufficient for SBS2: the cofactor network must be co-activated, and that activation depends on the cell's position in the viral lifecycle.

## Manuscript Structure

The paper is organized as eight figures, each addressing a specific question. Every figure has a corresponding analysis directory with a detailed walkthrough file documenting the pipeline architecture, scripts, input/output files, key results, and troubleshooting.

| Figure | Question | Analysis Directory | Walkthrough |
|--------|----------|-------------------|-------------|
| 1 | [A3 expression vs SBS2 mutagenesis](#figure-1-a3-expression-vs-sbs2-mutagenesis) | [`scripts/TCGA/`](scripts/TCGA/) | [`Question_1_workflow.md`](scripts/TCGA/Question_1_workflow.md) |
| 2 | [Bulk co-expression network](#figure-2-bulk-tcga-co-expression-network) | [`scripts/NETWORK/`](scripts/NETWORK/) | [`Question_2_workflow.md`](scripts/NETWORK/Question_2_workflow.md) |
| 3 | [Single-cell resolution](#figure-3-single-cell-resolution) | [`scripts/SINGLE_CELL/`](scripts/SINGLE_CELL/) | [`Question_3_workflow.md`](scripts/SINGLE_CELL/Question_3_workflow.md) |
| 4 | [Single-cell network analysis](#figure-4-single-cell-network-analysis) | [`scripts/NETWORK_SINGLE_CELL/`](scripts/NETWORK_SINGLE_CELL/) | [`Question_4_workflow.md`](scripts/NETWORK_SINGLE_CELL/Question_4_workflow.md) |
| 5 | [Patient-specific effects](#figure-5-patient-specific-effects) | [`scripts/PATIENT_SPECIFIC_EFFECTS/`](scripts/PATIENT_SPECIFIC_EFFECTS/) | [`Question_5_workflow.md`](scripts/PATIENT_SPECIFIC_EFFECTS/Question_5_workflow.md) |
| 6 | [HPV lifecycle](#figure-6-hpv-lifecycle) | [`scripts/HPV_ANALYSIS/`](scripts/HPV_ANALYSIS/) | [`Question_6_workflow.md`](scripts/HPV_ANALYSIS/Question_6_workflow.md) |
| 7 | [Neoantigen landscape](#figure-7-neoantigen-landscape) | [`scripts/NEOANTIGEN/`](scripts/NEOANTIGEN/) | [`NEOANTIGEN_WALKTHROUGH.md`](scripts/NEOANTIGEN/NEOANTIGEN_WALKTHROUGH.md) |

## Figure Summaries

### Figure 1: A3 Expression vs SBS2 Mutagenesis

**Question:** What is the fundamental relationship between A3 expression and observed A3-induced mutations in HNSCC?

A3A and A3B expression is necessary but not sufficient for SBS2 mutagenesis. A3A drives the SBS2 association in bulk (rho=0.149); A3B contributes minimally (rho=0.052) and saturates at the 79th percentile. Germline SNP analysis across 6,938 variants found zero BH-significant hits. Somatic enrichment analysis found that mutations distinguishing HIGH from LOW tumors (CASP8, HLA-A) reflect consequences of mutagenesis, not causes. Together, these negative results narrow the cofactor search to transcriptional regulation.

**Pipeline:** Three-branch architecture (RNA-seq expression, somatic mutation signatures, network group selection) converging on Figure 1 generation and enrichment analyses. See [`Question_1_workflow.md`](scripts/TCGA/Question_1_workflow.md).

### Figure 2: Bulk TCGA Co-expression Network

**Question:** Can differential co-expression network analysis identify gene communities associated with elevated SBS2 mutagenesis in A3-matched HNSCC tumors?

Differential co-expression network (HIGH minus LOW, 53 vs 53 tumors matched for A3 sum) identifies 11 main communities plus 33 satellites (1,131 genes, 2,599 edges). A3B lands in a 134-gene community with no significant KEGG enrichment. Zero Harris A3 interactors and only 1/24 cell-type markers are recovered. Bulk data cannot resolve A3-specific pathway biology, motivating the single-cell approach.

**Pipeline:** Eight-step modular pipeline (clean, merge, DE, correlation, threshold sweep, Leiden, centrality, figure generation). See [`Question_2_workflow.md`](scripts/NETWORK/Question_2_workflow.md).

### Figure 3: Single-Cell Resolution

**Question:** At single-cell resolution, are SBS2 mutational signatures localized to a specific cell type in HNSCC?

SBS2 mutagenesis localizes to basal epithelial cells co-expressing A3A and A3B. High SBS2 weight associates with A3A while high CNV burden associates with A3B, revealing divergent mutagenic fates within the basal compartment. 32.1% of SBS2-positive basal cells show no detectable A3A/A3B expression, consistent with pulsatile APOBEC induction. All single-cell data processed with ClusterCatcher (v1.3.0) and SRAscraper.

**Pipeline:** ClusterCatcher Snakemake pipeline (Cell Ranger, popV annotation, CytoTRACE2 + inferCNV, SComatic, SigProfiler) plus two post-pipeline figure generation scripts. See [`Question_3_workflow.md`](scripts/SINGLE_CELL/Question_3_workflow.md).

### Figure 4: Single-Cell Network Analysis

**Question:** Does single-cell differential co-expression network analysis resolve A3-specific cofactor biology that bulk data cannot?

Three pairwise networks (SBS2 vs CNV, SBS2 vs NORMAL, CNV vs NORMAL) recover what bulk cannot: RALY is the top hub across all three networks, 54-116 Harris A3 interactors are recovered (vs 0 in bulk), and a directionally coherent activating chain links A3A to immune/inflammatory genes (CCL20, LCN2, SMOX) through RALY and HNRNPA2B1. A3A and A3B show 100% negative DIFF edges (the "A3 wall"), consistent with their roles in opposing mutagenic programs.

**Pipeline:** Unified three-network pipeline using the same threshold selection and Leiden community detection framework as Figure 2. See [`Question_4_workflow.md`](scripts/NETWORK_SINGLE_CELL/Question_4_workflow.md).

### Figure 5: Patient-Specific Effects

**Question:** Is the SBS2-HIGH population a generalizable transcriptional program or an artifact of individual patients?

Three patients contribute 67.6% of SBS2-HIGH cells but are not transcriptionally distinct (PCA silhouette = 0.125). The co-expression network survives leave-one-patient-out removal with the A3 wall intact at 100%. What distinguishes high-contributor patients is somatic mutations in HPV infection and antigen presentation genes (HLA-B, HLA-C, TAP1), motivating the HPV lifecycle analysis.

**Pipeline:** Three-tier analysis (expression profiling, somatic variant analysis, network sensitivity testing) with six main scripts plus figure generation. See [`Question_5_workflow.md`](scripts/PATIENT_SPECIFIC_EFFECTS/Question_5_workflow.md).

### Figure 6: HPV Lifecycle

**Question:** How does HPV16 viral lifecycle stage determine the divergent mutagenic programs observed in APOBEC3-active epithelial cells?

SBS2-HIGH cells carry HPV16 in E2-regulated maintenance (25.7% reads to E1/E2, immune-visible, A3A-dominant). CNV-HIGH cells carry HPV16 in ATM-dependent productive infection (23.4% reads to L1/L2, 2x viral load, A3B-dominant). E6/E7 dosage is identical between groups. E2 is intact (not disrupted by integration). A panel of 53 host marker genes validates predictions from established HPV biology (Doorbar, Laimins, McBride labs).

**Pipeline:** Four-phase HPV data extraction pipeline (inventory, raw UMI recovery, threshold detection, genome alignment) plus lifecycle diagnostic and figure generation. See [`Question_6_workflow.md`](scripts/HPV_ANALYSIS/Question_6_workflow.md).

### Figure 7: Neoantigen Landscape

**Question:** How does HPV lifecycle stage shape the neoantigen landscape and therapeutic vulnerability of each population?

SBS2-HIGH cells produce 1.77x more predicted neoantigens per cell (2,361 vs 1,331). ANXA1 is the top target (IC50 = 13.1 nM, 42 neoantigen peptides, 99.5% cell coverage). The neoantigen landscape partitions into four therapeutic tiers: Tier 1A (22 shared + escaped, hot tumor targets), Tier 1B (276 SBS2-specific, ANXA1 leads), Tier 2 (135 CNV-specific + checkpoint blockade), Tier 3 (83 broad coverage backbone). RNA fusion rates are similar across groups; specific fusion partners drive asymmetric immune escape.

**Pipeline:** Six-step pipeline (prep, SnpEff annotation, MHCflurry binding with proteome-based peptides, expression-weighted ranking, fusion analysis, immune escape integration) plus figure generation. See [`NEOANTIGEN_WALKTHROUGH.md`](scripts/NEOANTIGEN/NEOANTIGEN_WALKTHROUGH.md).

## Computational Environment

All analyses were performed on the titan/zeus HPC cluster at Texas Biomedical Research Institute using SLURM job scheduling.

| Environment | Primary Use |
|-------------|------------|
| `NETWORK` | Network analysis, figure generation, scanpy, matplotlib |
| `NEOANTIGEN` | SnpEff annotation, MHCflurry binding prediction |
| `sc_pre` | Single-cell preprocessing (ClusterCatcher) |
| `SComatic` | Somatic mutation calling, SigProfiler |
| `RNA-seq_NovoGene` | TCGA data download and processing (R) |

## Data

All TCGA data were accessed through the GDC API under dbGaP authorization. Single-cell RNA-seq data are from GSE173468 (31 HPV+ HNSCC samples). The Ensembl GRCh38 release 115 reference proteome is used for neoantigen peptide generation.

## Authors

Jake D. Lehle, Mohadeseh Soleimanpour, Niloofar Haghjoo, Colin Rorex, Feng Li, Armando Mendez, Diako Ebrahimi

Texas Biomedical Research Institute
