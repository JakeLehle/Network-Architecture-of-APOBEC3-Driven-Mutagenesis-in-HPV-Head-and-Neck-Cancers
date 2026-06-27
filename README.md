# Deciphering the Network Architecture of APOBEC3-Driven Mutagenesis in HPV+ Head and Neck Cancers

## Overview

APOBEC3A (A3A) and APOBEC3B (A3B) are controlled in normal epithelium, but human papillomavirus (HPV) infection dismantles that control and turns the enzymes against the host genome in head and neck squamous cell carcinoma (HNSCC). A3 expression is necessary but not sufficient for the SBS2 mutational signature, and bulk tissue cannot resolve why, because it averages across immune-compartment A3 expression and across tumor cells caught at different stages of the viral lifecycle. This study profiles single cells from HPV16-positive HNSCC tumors and matched normal-adjacent tissue (GSE173468) to resolve the cofactors, lifecycle context, and neoantigen consequences of A3-driven mutagenesis.

Within the basal epithelial compartment, A3A-linked SBS2 point mutations and A3B-linked copy-number change (CNV) diverge into two states. Differential co-expression networks built tumor against normal support a functional coactivator role for the A3 interactors *RALY* and *HNRNPA2B1*, and the enzymes themselves are co-expression-decoupled from their partners in tumor (the A3 wall), consistent with pulsatile A3 induction. The SBS2 and CNV programs track the maintenance and productive stages of the HPV16 lifecycle, so the SBS2:CNV ratio offers a molecular estimate of tumor age, and the neoantigens of the immune-visible SBS2-HIGH and immune-evasive CNV-HIGH states define tiered candidates for lifecycle-matched mRNA vaccines.

The cohort-level TCGA validation of these findings is tracked separately in a companion repository; this repository covers the single-cell study only.

## Manuscript Structure

The paper is organized as seven figures across five analysis questions. Each question has an analysis directory under `scripts/` and a walkthrough documenting the pipeline architecture, scripts, inputs and outputs, and key results.

| Question | Figures | Analysis Directory | Walkthrough |
|----------|---------|--------------------|-------------|
| 1. [Single-cell localization and divergence](#question-1-single-cell-localization-and-divergence) | 1-3 | [`scripts/SINGLE_CELL/`](scripts/SINGLE_CELL/) | [`Question_1_workflow.md`](Question_1_workflow.md) |
| 2. [Single-cell co-expression networks](#question-2-single-cell-co-expression-networks) | 4 | [`scripts/NETWORK_SINGLE_CELL/`](scripts/NETWORK_SINGLE_CELL/) | [`Question_2_workflow.md`](Question_2_workflow.md) |
| 3. [Patient-specific effects](#question-3-patient-specific-effects) | 5 | [`scripts/PATIENT_SPECIFIC_EFFECTS/`](scripts/PATIENT_SPECIFIC_EFFECTS/) | [`Question_3_workflow.md`](Question_3_workflow.md) |
| 4. [HPV16 lifecycle](#question-4-hpv16-lifecycle) | 6 | [`scripts/HPV_ANALYSIS/`](scripts/HPV_ANALYSIS/) | [`Question_4_workflow.md`](Question_4_workflow.md) |
| 5. [Neoantigen landscape](#question-5-neoantigen-landscape) | 7 | [`scripts/NEOANTIGEN/`](scripts/NEOANTIGEN/) | [`Question_5_workflow.md`](Question_5_workflow.md) |

## Figure Summaries

### Question 1: Single-Cell Localization and Divergence

**Figures 1-3.** From 155,650 cells across 14 patients, A3A and A3B expression concentrates almost exclusively in the 52,126 basal epithelial cells. SBS2 burden localizes to that compartment, associating with A3A more than A3B, while CNV and stemness map instead to an A3B-associated subpopulation. Splitting basal cells by A3 dominance and placing them on a within-tissue SBS2-versus-CNV axis resolves the divergence that global expression obscures, and only in tumor. A subset of SBS2-positive cells carry no detectable A3A or A3B at capture, consistent with pulsatile A3 induction. All single-cell processing runs through ClusterCatcher (v1.3.0) and SRAscraper. See [`Question_1_workflow.md`](Question_1_workflow.md).

### Question 2: Single-Cell Co-expression Networks

**Figure 4.** Differential co-expression networks compare each tumor population to normal-adjacent basal cells. The SBS2-HIGH network (2,948 genes, 23 gene groups) and the CNV-HIGH network (4,886 genes, 38 gene groups) recover 54 and 109 of the 174 catalogued A3 interactors. A *RALY*-anchored activating program is conserved across both states, with *HNRNPA2B1* anchoring *CCL20* in SBS2-HIGH, nominating both as functional coactivators. A cornification differentiation program forms the clearest inhibiting chain, and every edge linking A3A or A3B to its neighbors is negative (the A3 wall), consistent with pulsatile induction. See [`Question_2_workflow.md`](Question_2_workflow.md).

### Question 3: Patient-Specific Effects

**Figure 5.** Three patients contribute 74.1% of the SBS2-HIGH cells, but they are not transcriptionally distinct: their cells intermingle with everyone else's (silhouette = 0.021), and the network survives leave-one-patient-out removal with the A3 wall intact at 100%. What distinguishes the high contributors is somatic variants in HPV16-associated genes rather than the A3 program, several of which sit inside the network (*HLA-C*, *MX1*, *MDM2*), pointing toward the virus. See [`Question_3_workflow.md`](Question_3_workflow.md).

### Question 4: HPV16 Lifecycle

**Figure 6.** SBS2-HIGH cells carry HPV16 in maintenance (more E1-driven maintenance reads, immune-visible, A3A-dominant), while CNV-HIGH cells carry it in productive infection (more L1/L2 capsid reads, an ATM-dependent damage response, G2/M arrest, A3B-dominant), at roughly 2.6-fold higher viral load. E6/E7 dosage does not differ and E2 is intact, so both populations carry episomal virus. A cell's lifecycle position sets both the active enzyme and the class of damage, and the SBS2:CNV ratio estimates tumor age along the maintenance-to-productive axis. See [`Question_4_workflow.md`](Question_4_workflow.md).

### Question 5: Neoantigen Landscape

**Figure 7.** SBS2-HIGH cells produce 1.77-fold more predicted neoantigens than CNV-HIGH cells (2,361 vs 1,331) on an A3-skewed mutational background. The 516 neoantigen-producing genes partition into three priority tiers: Tier 1, 105 broadly conserved shared targets led by *MDK* (with a 22-gene escape subset under direct immune selection); Tier 2, 276 SBS2-specific hot-tumor targets led by *ANXA1*; and Tier 3, 135 CNV-specific cold-tumor targets. *ANXA1* is the top single target, expressed in 99.5% of SBS2-HIGH cells with a best IC50 of 13.1 nM. See [`Question_5_workflow.md`](Question_5_workflow.md).

## Computational Environment

All analyses ran on the titan/zeus HPC cluster at Texas Biomedical Research Institute under SLURM.

| Environment | Primary Use |
|-------------|-------------|
| `ClusterCatcher` | Single-cell preprocessing, annotation, mutation calling, signatures (bundles its own sub-environments) |
| `NETWORK` | Co-expression networks, patient and lifecycle analysis, figure generation |
| `NEOANTIGEN` | SnpEff annotation and MHCflurry binding prediction |

## Data

Single-cell RNA-seq data are from GSE173468: 14 patients, 44 samples, 155,650 cells, of which 52,126 are basal epithelial. Somatic signatures are refit against COSMIC v3.4, and neoantigen peptides are generated against the Ensembl GRCh38 release 115 reference proteome. Raw reads are retrieved with SRAscraper and processed through the ClusterCatcher pipeline.

## Authors

Jake D. Lehle, Mohadeseh Soleimanpour, Niloofar Haghjoo, Colin Rorex, Feng Li, Armando Mendez, Rachael Rodriguez, Elizabeth Sommer, Diako Ebrahimi

Texas Biomedical Research Institute
