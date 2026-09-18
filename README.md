# Deciphering the Network Architecture of APOBEC3-Driven Mutagenesis in HPV+ Head and Neck Cancers

Analysis code for the single-cell study of APOBEC3-driven mutagenesis in HPV16-positive head and neck squamous cell carcinoma.

## Overview

APOBEC3A (A3A) and APOBEC3B (A3B) are controlled in normal epithelium, but human papillomavirus (HPV) infection dismantles that control and turns the enzymes against the host genome in head and neck squamous cell carcinoma (HNSCC). A3 expression is necessary but not sufficient for the SBS2 mutational signature, and bulk tissue cannot resolve why, because it averages across immune-compartment A3 expression and across tumor cells caught at different stages of the viral lifecycle. This study profiles single cells from HPV16-positive HNSCC tumors and matched normal-adjacent tissue (GSE173468) to resolve the cofactors, lifecycle context, and neoantigen consequences of A3-driven mutagenesis.

Within the basal epithelial compartment, A3A-linked SBS2 point mutations and A3B-linked copy-number alteration (CNA) diverge into two states. Differential co-expression networks built tumor against normal nominate the A3 interactors *RALY* and *HNRNPA2B1* as candidate coactivators, and the enzymes themselves are co-expression-decoupled from their partners in tumor (the A3 wall), consistent with pulsatile A3 induction. The SBS2 and CNA programs track the maintenance and productive stages of the HPV16 lifecycle, so the SBS2:CNA ratio may offer a molecular estimate of tumor age, and the neoantigens of the immune-visible SBS2-HIGH and immune-evasive CNA-HIGH states define tiered candidates for lifecycle-matched mRNA vaccines.

The cohort-level TCGA validation of these findings is tracked separately in a companion repository; this repository covers the single-cell study only.

---

## Before you start: three naming traps

The repository predates several decisions made late in writing. Nothing below is a bug, but each will mislead a reader who assumes the names track the paper.

**1. `CNV` in the code is `CNA` in the paper.** The manuscript writes copy-number alteration as CNA. Every column, directory and group name in the code uses `CNV` (`cnv_score`, `CNV_HIGH`, `NETWORK_CNV_VS_NORMAL`). They are the same quantity.

**2. `data/FIG_N/` directories do not track manuscript figure numbers.**

| Data directory | Holds |
|---|---|
| `data/FIG_3/` | Figure 3 |
| `data/FIG_4/` | **Figure 5** (networks) and the Figure 4 inputs |
| `data/FIG_5/` | **Supplementary Figures 7-9** (patient effects) |
| `data/FIG_6/` | Figure 6 |
| `data/FIG_7/` | Figure 7 |

Consequently `Generate_Figure4_Panels.py` produces manuscript **Figure 5**, and `Generate_Supplemental_Patient_Effects.py` writes files named `Supp_Panel_*` into `FIGURE_5_PANELS/`.

**3. Question numbering and manuscript section order diverge in two places.** Question 3 maps to Results **4.4** and Question 4 maps to Results **4.3**. The manuscript presents the HPV16 lifecycle first, because the lifecycle stages are what the patient-level thresholds are built on. Question 3 also depends on per-cell viral tables produced by Question 4, so it cannot be run first.

---

## Repository structure

The analysis is organized as five questions. Each has a directory under `scripts/` and a walkthrough documenting architecture, scripts, inputs and outputs, and key results.

| Question | Manuscript section | Figures | Analysis directory | Walkthrough |
|----------|--------------------|---------|--------------------|-------------|
| 1. [Single-cell localization and divergence](#question-1-single-cell-localization-and-divergence) | Results 4.1 | 1-3, Supp. 1-4 | [`scripts/SINGLE_CELL/`](scripts/SINGLE_CELL/) | [`Question_1_workflow.md`](Question_1_workflow.md) |
| 2. [Single-cell co-expression networks](#question-2-single-cell-co-expression-networks) | Results 4.2 | 4, 5, Supp. 5 | [`scripts/NETWORK_SINGLE_CELL/`](scripts/NETWORK_SINGLE_CELL/) | [`Question_2_workflow.md`](Question_2_workflow.md) |
| 3. [Patient-specific effects](#question-3-patient-specific-effects) | Results 4.4 | Supp. 7-9 | [`scripts/PATIENT_SPECIFIC_EFFECTS/`](scripts/PATIENT_SPECIFIC_EFFECTS/) | [`Question_3_workflow.md`](Question_3_workflow.md) |
| 4. [HPV16 lifecycle](#question-4-hpv16-lifecycle) | Results 4.3 | 6, Supp. 6 | [`scripts/HPV_ANALYSIS/`](scripts/HPV_ANALYSIS/) | [`Question_4_workflow.md`](Question_4_workflow.md) |
| 5. [Neoantigen landscape](#question-5-neoantigen-landscape) | Results 4.5 | 7 | [`scripts/NEOANTIGEN/`](scripts/NEOANTIGEN/) | [`Question_5_workflow.md`](Question_5_workflow.md) |

### Arriving from the manuscript

| Manuscript figure | Produced by |
|---|---|
| Fig. 1, Fig. 2 | assembled from ClusterCatcher outputs (`scripts/SINGLE_CELL/`) |
| Fig. 3 | `SINGLE_CELL/Generate_Figure3_A3_Dominance_Slider.py` |
| Fig. 4 (DE and KEGG) | see Question 2 open items |
| Fig. 5 (networks) | `NETWORK_SINGLE_CELL/Generate_Figure4_Panels.py` |
| Fig. 6 (HPV16 lifecycle) | `HPV_ANALYSIS/Generate_Figure6_Lifecycle_Panels.py` |
| Fig. 7 (neoantigens) | `NEOANTIGEN/Generate_Figure7_Panels.py` |
| Supp. Fig. 3 | `SINGLE_CELL/Step02_Supplemental_Marker_Validation.py` |
| Supp. Fig. 5 | `NETWORK_SINGLE_CELL/Generate_Figure4_Supplemental.py` |
| Supp. Fig. 6 | `HPV_ANALYSIS/TROUBLESHOOTING/Generate_Supp_CellCycle.py` |
| Supp. Fig. 7 | `PATIENT_SPECIFIC_EFFECTS/Generate_Supplemental_Patient_Effects.py` + `HPV_ANALYSIS/TROUBLESHOOTING/Generate_Supp_Contribution_Virus.py` |
| Supp. Fig. 8 | `PATIENT_SPECIFIC_EFFECTS/Generate_Supplemental_Patient_Effects.py` (Panel D) |
| Supp. Fig. 9 | `HPV_ANALYSIS/TROUBLESHOOTING/Generate_Supp_Enzyme_Conjunction.py` |

Three supplemental figures are currently produced from `TROUBLESHOOTING/` directories, which the walkthroughs otherwise exclude. Promoting them is an open item.

---

## Figure summaries

### Question 1: Single-Cell Localization and Divergence

**Results 4.1, Figures 1-3.** From 155,650 cells across 14 patients and 44 samples, A3A and A3B expression concentrates almost exclusively in the 52,126 basal epithelial cells (A3A mean 1.62, 25.7% positive; A3B mean 1.69, 33.2% positive). Among 31,912 basal cells with somatic calls, 5,911 (18.5%) carry SBS2, which associates with A3A (rho = 0.149) more than A3B (rho = 0.052). CNA and stemness map instead to an A3B-associated subpopulation, with A3B correlating positively with CNA (rho = 0.147) and A3A negatively (rho = -0.199).

Most SBS2-positive cells carry no detectable A3A (56.7%) or A3B (49.8%) at capture, and 32.1% carry neither, consistent with pulsatile A3 induction or transcript dropout. SBS13, the second A3-associated signature, is detected but is the weakest of the fifteen fitted signatures and shows no association with A3A (rho = 0.006, p = 0.32), so SBS2 is the readout throughout.

Splitting basal cells by A3 dominance and placing them on a within-tissue SBS2-versus-CNA axis resolves the divergence that global expression obscures: A3A-dominant cells shift toward SBS2 (median +0.60) and A3B-dominant toward CNA (median -0.77), in tumor only. All single-cell processing runs through ClusterCatcher (v1.3.0) and SRAscraper.

### Question 2: Single-Cell Co-expression Networks

**Results 4.2, Figures 4 and 5.** Three size-matched basal populations of 546 cells (SBS2-HIGH, CNV-HIGH, NORMAL) are compared by differential expression and by differential co-expression. The SBS2-HIGH network (2,948 genes, **25** gene groups) and the CNV-HIGH network (4,886 genes, 38 gene groups) recover 54 and 109 of the 174 catalogued A3 interactors.

Every edge linking A3A or A3B to a first-degree neighbour is negative in both networks, zero positive out of 222 (the A3 wall), so edge sign cannot identify regulators and the analysis extends to coherent chains beyond the first-degree neighbourhood. A *RALY*-anchored activating chain is conserved across both states (six genes in SBS2-HIGH with *LCN2*, *KRT24*, *LINC00278*, *CHMP4B* and *UTY*; thirteen genes in CNV-HIGH), with *HNRNPA2B1* anchoring *CCL20* in SBS2-HIGH, nominating both interactors as candidate coactivators. A 70-node cornification and differentiation program forms the clearest inhibiting chain.

### Question 3: Patient-Specific Effects

**Results 4.4, Supplementary Figures 7-9.** Both tumor populations are concentrated in a few patients: three (SC013, SC029, SC001) contribute 74.0% of the SBS2-HIGH cells and two (SC027, SC001) contribute 82.1% of the CNV-HIGH cells, with only SC001 in both.

The co-expression program nevertheless survives leave-one-patient-out reconstruction with the **A3 wall 100% intact in every run**. Activating-chain genes are largely retained when SC029 or SC001 is removed (9 of 10 each) and less so when SC013 is removed (3 of 10), the patient supplying 40.5% of the population.

Expression level does not predict contribution; prevalence does. The highest A3A expressor in the cohort contributes to neither fate, mean expression among expressing cells varies only 1.4-fold across patients, and the fraction of cells expressing varies far more. Three patient-level conditions together separate the contributors: viral load, viral lifecycle direction, and prevalence of the matching A3 enzyme (Fisher exact p = 0.0027 and p = 0.011).

> An earlier version of this analysis reported a somatic-variant-sharing narrative (HC-exclusive variants, KEGG HPV enrichment, network overlap through *HLA-C*, *MX1* and *MDM2*). **That work is not in the submitted manuscript.** The scripts remain and are documented as analyses performed; see the Question 3 walkthrough.

### Question 4: HPV16 Lifecycle

**Results 4.3, Figure 6.** HPV16 reads localize to the basal compartment (94.6% of positive cells), with an L-method breakpoint at 8 UMIs giving 15,927 HPV16-positive basal cells. SBS2-HIGH cells carry HPV16 in maintenance (maintenance reads 25.9% vs 13.2%, driven by E1; immune-visible antigen presentation and interferon effectors; A3A-dominant), while CNV-HIGH cells carry it in productive infection (capsid reads 17.8% vs 10.3%; a damage response spanning both the ATM and ATR arms; G2/M arrest; A3B-dominant) at roughly 2.6-fold higher viral load.

Of 59 host genes across eight functional tiers, 47 differ significantly between the two populations. E6/E7 dosage does not differ, accounting for under 1% of viral reads in both, and E2 is intact with no evidence of the loss that accompanies integration, so both populations carry episomal virus. A cell's lifecycle position appears to set both the active enzyme and the class of damage, and the SBS2:CNA ratio may estimate tumor position along the maintenance-to-productive axis.

### Question 5: Neoantigen Landscape

**Results 4.5, Figure 7.** SBS2-HIGH cells produce 1.82-fold more predicted neoantigen-forming mutations than CNV-HIGH cells (560 vs 308, p = 9.50e-18), an excess that survives normalization for sequencing depth (0.343 vs 0.184 per 1,000 UMI, BH p = 2.15e-06). RNA fusion burden does not track the immune divergence (0.547 vs 0.506 per 1,000 UMI, BH p = 0.327).

The 775 unique neoantigen-forming mutations partition into three tiers: **93 shared** across both viral states, **467 SBS2-specific** for hot, immune-visible tumors, and **215 CNV-specific** for cold, immune-evasive tumors. Fifteen of the 93 shared mutations show patterns consistent with removal in CNV-HIGH cells, through fusion disruption or transcriptional silencing.

Candidates are ranked by clonal prevalence rather than expression, which surfaces one lead per tier: *COX4I1* p.Ala9Thr (CNV-specific, 24.9%), *SPRR1A* p.Val61Ile (shared, 16.9%), and *KRT6B* p.Glu342Lys (SBS2-specific, 7.9%). *KRT6B* arises in a clean T[C>T]W context and converts a non-binding wild-type peptide (IC50 27,793 nM) into a strong binder (81 nM). MHC-I binding is predicted across ten common HLA class I alleles covering 90.5% of the world population.

> Two claims from earlier drafts are retired. The expression-weighted ranking that put *ANXA1* on top has been replaced by clonal prevalence, because the *ANXA1* neoantigen is carried in only 1.1% of SBS2-HIGH cells despite near-universal expression. And genome-verified recomputation shows **no TCW-context enrichment** between the two populations under any of three definitions; the earlier enrichment signal traced to a SComatic annotation field and is not reproducible.

---

## Computational environment

All analyses ran on the titan/zeus HPC cluster at Texas Biomedical Research Institute under SLURM.

| Environment | Primary use |
|-------------|-------------|
| `ClusterCatcher` | Single-cell preprocessing, annotation, mutation calling, signatures (bundles its own per-rule sub-environments) |
| `NETWORK` | Co-expression networks, patient and lifecycle analysis, figure generation |
| `NEOANTIGEN` | SnpEff annotation, STAR chimeric alignment, MHCflurry binding prediction |
| `sc_pre` | Some single-cell diagnostics |

Key external tools: Cell Ranger, popV (Tabula Sapiens reference), Scrublet, SComatic, CytoTRACE2, inferCNVpy, Kraken2, minimap2, Leiden, SnpEff, STAR, MHCflurry, gseapy.

## Data and reference

| Resource | Detail |
|---|---|
| Single-cell RNA-seq | [GSE173468](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173468): 14 patients, 44 samples, 155,650 cells, 52,126 basal epithelial |
| Genome | GRCh38 (Cell Ranger `refdata-gex-GRCh38-2020-A`) |
| Proteome | Ensembl GRCh38 release 115, `pep.all.fa` (245,535 isoforms) |
| Mutational signatures | COSMIC v3.4, HNSCC-relevant subset, 15 signatures retained |
| Virus | HPV16 NC_001526.4 (7,906 bp) |
| A3 interactors | 174 genes compiled from published affinity-purification datasets |

Raw reads are retrieved with [SRAscraper](https://github.com/Diako-Lab/SRAscraper) and processed through the [ClusterCatcher](https://github.com/Diako-Lab/ClusterCatcher) single-cell pipeline, both released with this publication.

## Reproducing the analysis

Run in question order, except that Question 4 must precede Question 3.

```bash
# Question 1: preprocessing through to adata_final.h5ad
cd scripts/SINGLE_CELL/          && sbatch Run_Cluster_Catcher_Pipeline.sh

# Question 2: networks (compute, then figures)
cd scripts/NETWORK_SINGLE_CELL/  && sbatch RUN_THREE_NETWORK_PIPELINE.sh
                                 && bash MAKE_FIGURES.sh

# Question 4: HPV16 phases 1-4, then diagnostics and figure
cd scripts/HPV_ANALYSIS/         && sbatch RUN_HPV_ANALYSIS.sh

# Question 3: patient effects
cd scripts/PATIENT_SPECIFIC_EFFECTS/ && sbatch RUN_PATIENT_ANALYSIS.sh

# Question 5: neoantigens
cd scripts/NEOANTIGEN/           && sbatch RUN_NEOANTIGEN_PIEPLINE.sh
```

Several manuscript numbers are produced by assertion-based audit scripts rather than by the pipeline steps. These recompute each claim from on-disk output and report PASS or FAIL per number:

| Audit | Covers |
|---|---|
| `SINGLE_CELL/TROUBLESHOOTING/Diagnostic_section3_numbers.py` | Results 4.1 |
| `NETWORK_SINGLE_CELL/Diagnostic_A3_Interactor_Concordance.py` | Results 4.2 |
| `HPV_ANALYSIS/Diagnostic_Figure6_HostMarkers_and_IntegrationProxy.py` | Results 4.3 (run twice; see walkthrough) |
| `HPV_ANALYSIS/TROUBLESHOOTING/Diagnostic_Patient_Determinants_Table.py` | Results 4.4 |
| `NEOANTIGEN/diagnostic_section7_numbers.py` | Results 4.5 (run the group-aware diagnostic first) |

## Repository conventions

- `TROUBLESHOOTING/` directories hold diagnostics and prior drafts and are excluded from the walkthroughs, with the documented exceptions above.
- Figures are saved as PDF and PNG at 300 DPI, font sizes 28-34, with a consistent palette: SBS2 and A3A coral `#ed6a5a`, CNA and A3B mustard `#F6D155`, NORMAL blue `#4682b4`.
- `signature_weights_per_cell.txt` is written signatures x cells and requires transposition on load.
- Figure generation is a separate stage from compute, so panels can be revised without rebuilding networks.

## Citation

Lehle JD, Soleimanpour M, Haghjoo N, Rorex C, Li F, Mendez A, Rodriguez R, Sommer E, Chiang C-M, Ebrahimi D. *Deciphering the Network Architecture of APOBEC3-Driven Mutagenesis in HPV+ Head and Neck Cancers.* bioRxiv (2026).

Corresponding author: Diako Ebrahimi, debrahimi@txbiomed.org

Texas Biomedical Research Institute, San Antonio, TX; Department of Molecular Microbiology and Immunology, University of Texas at San Antonio; Simmons Comprehensive Cancer Center and Departments of Biochemistry and Pharmacology, University of Texas Southwestern Medical Center.
