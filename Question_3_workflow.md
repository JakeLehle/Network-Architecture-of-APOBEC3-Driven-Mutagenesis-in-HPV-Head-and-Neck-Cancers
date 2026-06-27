# Question 3 Workflow: Patient-Specific Effects (Figure 5)

## Question 3: Is the SBS2-HIGH population a conserved cross-patient program, or an artifact of a few individuals?

### Overview

Figure 5 asks whether the SBS2-HIGH basal population from the Figure 4 networks is a generalizable transcriptional program across patients or an artifact driven by one or two individuals, a well-known confounder in single-cell cancer studies. The analysis runs in three tiers of increasing stringency: expression profiling, somatic variant analysis, and network sensitivity testing.

The central finding is that although three patients (SC013, SC029, SC001) contribute 74.1% of the SBS2-HIGH cells, they are not transcriptionally distinct. Their HIGH cells intermingle with other patients in expression space (silhouette = 0.021), and the co-expression network survives leave-one-patient-out removal with the A3 wall intact at 100%. What distinguishes the high contributors is somatic variants in HPV16-associated genes rather than in the A3 program itself, several of which sit inside the Figure 4 network. This links the patient-level signal back to the virus and motivates the HPV16 lifecycle analysis in Figure 6.

---

## Directory Structure

```
scripts/PATIENT_SPECIFIC_EFFECTS/
├── patient_config.py                              # Shared configuration (paths, parameters, utilities)
├── RUN_PATIENT_ANALYSIS.sh                        # SLURM master script (runs all tiers)
│
├── Tier1A_Patient_Expression_GSEA.py              # Per-patient Wilcoxon DE + KEGG GSEA
├── Tier1B_HIGH_Cell_Transcriptional_Similarity.py # PCA, silhouette, profile correlations
├── Tier1C_Enriched_vs_Depleted_Patient_DE.py      # High-contributor vs other-patient DE
│
├── Tier2A_Expression_Haplotype_Proxies.py         # A3 ratios, activating-chain expression
├── Tier2B_SNP_Pattern_Analysis.py                 # SComatic variant sharing across patients
│
├── Tier3A_Leave_One_Patient_Out.py                # LOPO network sensitivity (-> Supp Fig 7)
├── Tier3B_Tumor_Only_Sensitivity.py               # Tumor-only control check (robustness, not a figure panel)
│
├── Generate_Supplemental_Patient_Effects.py       # Emits panels A-F (-> Fig 5 a-e + Supp Fig 7)
├── Analyze_SNP_Tier_Genes.py                      # Variant-to-gene mapping + KEGG per sharing tier
├── Diagnostic_HC_Network_HPV_Overlap.py           # HC-exclusive x network x HPV/immune overlap
├── Diagnostic_Recompute_Patient_Contributions.py  # Post-reselection contribution recompute
├── Diagnostic_Regenerate_Panel_A_C.py             # Post-reselection panel A-C regeneration
└── TROUBLESHOOTING/                               # Diagnostics (not documented here)
```

---

## Input Data

All inputs come from the Figure 4 single-cell network pipeline. `patient_config.py` defines the paths.

| Input | Path | Description |
|-------|------|-------------|
| AnnData | `data/FIG_4/00_input/adata_final.h5ad` | 155,650 cells, 27,736 genes (ClusterCatcher output) |
| Signature weights | `data/FIG_4/00_input/signature_weights_per_cell.txt` | Per-cell SBS signature weights |
| Three-group assignments | `data/FIG_4/01_group_selection/three_group_assignments.tsv` | 546 SBS2-HIGH, 546 CNV-HIGH, 546 NORMAL |
| DE genes | `data/FIG_4/NETWORK_SBS2_VS_NORMAL/02_differential_expression/` | scanpy Wilcoxon DE (FDR < 0.05) |
| Network partition | `data/FIG_4/NETWORK_SBS2_VS_NORMAL/04_communities/SC_best_partition.csv` | 2,948 genes in Leiden communities |
| A3 interactors | `data/FIG_4/00_input/Harris_A3_interactors.txt` | 174 catalogued A3 interactors |
| SComatic calls | `results_NMF_v0.1.1/all_samples.single_cell_genotype.filtered.tsv` | Per-cell somatic variant calls |
| GTF | `SC/ref/GRCh38/genes/genes_unzipped.gtf` | GRCh38 annotations for variant-to-gene mapping |

---

## Pipeline Execution

The tiered pipeline runs as a single SLURM job.

```bash
cd scripts/PATIENT_SPECIFIC_EFFECTS/
sbatch RUN_PATIENT_ANALYSIS.sh
```

Three phases run sequentially: Phase 1 (Tier 1A, 1B, 1C, 2A expression analyses), Phase 2 (Tier 2B variant pattern analysis), Phase 3 (Tier 3A LOPO across the three high contributors). Post-pipeline scripts run interactively:

```bash
conda run -n NETWORK python Generate_Supplemental_Patient_Effects.py
conda run -n NETWORK python Analyze_SNP_Tier_Genes.py
conda run -n NETWORK python Diagnostic_HC_Network_HPV_Overlap.py
```

---

## Script Descriptions

### Shared Configuration

**`patient_config.py`** defines all paths, parameters, and utilities. Key choices: three-group loading from `three_group_assignments.tsv` (546 each), with a backward-compatible wrapper returning (SBS2-HIGH, NORMAL) for the SBS2_VS_NORMAL analyses; network parameters matching the Figure 4 pipeline (scanpy DE at FDR < 0.05, DIFF threshold auto-selection by max fragmentation rate, full-network Leiden with resolution sweep and component-aware merge); the chain gene sets from the Figure 4 concordance analysis (the *RALY* / *HNRNPA2B1*-anchored activating chain, plus inhibiting-chain anchors including *ZNG1A*); and the style constants (font sizes 28-34, hex colors, PDF + PNG at 300 DPI).

### Tier 1: Expression

**`Tier1A_Patient_Expression_GSEA.py`** compares each patient's basal cells against all others (Wilcoxon), then runs GSEA prerank against KEGG, an unbiased screen for patient-specific pathway activity, and correlates per-patient A3 expression with SBS2-HIGH enrichment.

**`Tier1B_HIGH_Cell_Transcriptional_Similarity.py`** asks whether the 546 SBS2-HIGH cells cluster by patient or intermingle. PCA on the 3,907 DE genes (FDR < 0.05), silhouette score on patient label, and a per-patient mean-profile correlation matrix. The low silhouette (0.021) indicates a shared program, not patient-specific clusters.

**`Tier1C_Enriched_vs_Depleted_Patient_DE.py`** compares all basal cells of the three high contributors against the remaining patients, testing whether the high contributors are constitutively different independent of SBS2 status.

### Tier 2: Variants

**`Tier2A_Expression_Haplotype_Proxies.py`** measures expression proxies for germline A3 haplotype variation (A3A/A3B ratio, A3H expression) and the activating-chain expression program, testing whether high contributors carry an elevated baseline of the mutagenic program.

**`Tier2B_SNP_Pattern_Analysis.py`** maps the basal somatic variants to patients and classifies each by sharing breadth (universal, ten or more patients; broadly shared, five to nine; partially shared, two to four; patient-specific, one; and the HC-exclusive set carried by at least two of the three high contributors and no one else), with pairwise Jaccard similarity and per-patient composition profiles.

### Tier 3: Network Sensitivity

**`Tier3A_Leave_One_Patient_Out.py`** is the core test. For each high contributor it removes that patient from both SBS2-HIGH and NORMAL, then re-runs the full network pipeline (scanpy DE, Spearman matrices, max-fragmentation-rate threshold, full-network Leiden, component-aware merge), scoring four metrics: activating-chain recovery, community-structure ARI against the full network, gene-overlap Jaccard, and A3 wall integrity (percent negative DIFF edges on A3A). This feeds Supplemental Figure 7.

**`Tier3B_Tumor_Only_Sensitivity.py`** is a control check, not a figure panel. It re-runs the network two ways, once dropping the 25 normal-tissue HIGH cells while keeping the original controls, and once re-matching controls from tumor-only basal cells, to confirm that any A3A dropout traces to control reselection rather than to the loss of those normal-tissue HIGH cells.

### Post-Pipeline

**`Generate_Supplemental_Patient_Effects.py`** emits six panels. In the current figure layout these map as: Panel A (patient distribution) to Fig 5a, Panel B (per-patient A3A/A3B expression) to Fig 5b, Panel C (PCA + correlation) to Fig 5c, Panel E (five-tier variant sharing spectrum) to Fig 5d, Panel F (genome ideogram of variant positions by sharing tier) to Fig 5e, and Panel D (LOPO sensitivity) to Supplemental Figure 7. The `Supp_Panel_*` file names are a holdover from when this whole analysis was supplemental.

**`Analyze_SNP_Tier_Genes.py`** maps variants to genes via the GTF and cross-references each sharing tier against the network partition, the catalogued interactors, the A3 enzymes, and the activating and inhibiting chains, with KEGG enrichment per tier.

**`Diagnostic_HC_Network_HPV_Overlap.py`** closes the loop, identifying which HC-exclusive variant genes in significant HPV and immune KEGG pathways also sit inside the co-expression network, and which community they belong to.

**`Diagnostic_Recompute_Patient_Contributions.py`** and **`Diagnostic_Regenerate_Panel_A_C.py`** are post-reselection helpers that refresh the contribution numbers and regenerate panels A-C.

---

## Key Findings

1. **SBS2-HIGH is a shared program, not patient-specific.** The 546 cells come from 12 of 14 patients; three are disproportionately enriched (chi-squared p = 4.46e-309): SC013 (221 cells, 40.5%, 6.9-fold), SC029 (133, 24.4%, 2.1-fold), SC001 (50, 9.2%, 2.5-fold), together 74.1% of the population. But on PCA their cells intermingle with everyone else (silhouette = 0.021).

2. **A3 expression does not predict contribution.** The high contributors expressed A3A and A3B in the upper third of the cohort but were not the highest. SC005, the highest mean A3A (3.02), was under-enriched (0.7-fold) despite the largest basal pool in the dataset (9,502 cells).

3. **The network is robust to patient removal.** Removing any single high contributor shrank the network but did not collapse it: the A3 wall stayed 100% intact, activating-chain genes were largely preserved, and community structure stayed consistent (Supp Fig 7). Tier 3B confirmed this independently of control selection.

4. **High contributors are distinguished by HPV16-associated variants, not the A3 program.** Of 44,221 basal somatic variants, the HC-exclusive set was 532 variants across 383 genes. It contained no A3 enzyme and none of the activating or inhibiting chain genes; the A3 machinery and its cofactors (the activators *HNRNPA2B1* and *CCL20*, the inhibitory-chain member *ZNG1A*) were mutated across the cohort rather than specifically in the high contributors. The HC-exclusive set instead recovered HPV infection (KEGG adj p = 0.02).

5. **HC-exclusive HPV-pathway genes overlap the network.** Several genes behind that enrichment are members of the differential co-expression network, including the antigen-presentation gene *HLA-C*, the interferon-stimulated gene *MX1*, and the p53 regulator *MDM2*. Variants in these could let HPV-infected cells persist and keep accumulating A3-driven mutations, motivating the HPV16 lifecycle analysis in Figure 6.

---

## Figure 5 Panels

| Panel | Content | Source |
|-------|---------|--------|
| Fig 5a | Patient distribution of SBS2-HIGH cells | `Generate_Supplemental_Patient_Effects.py` (Panel A) |
| Fig 5b | Per-patient A3A and A3B expression | Panel B |
| Fig 5c | PCA of HIGH cells by patient + correlation heatmap | Panel C |
| Fig 5d | Five-tier variant sharing spectrum | Panel E |
| Fig 5e | Genome ideogram of variant positions by sharing tier | Panel F |
| Supp Fig 7 | LOPO network sensitivity (chain recovery, ARI, Jaccard, A3 wall) | Panel D |

---

## Output Structure

```
data/FIG_5/
├── 00_diagnostics/          patient enrichment recompute
├── 01_patient_expression/   GSEA, A3 summary, PCA, silhouette, Tier1C
├── 02_snp_haplotype/        haplotype proxies, variant sharing, gene_analysis/
├── 03_sensitivity/          LOPO_SC013/ LOPO_SC029/ LOPO_SC001/ (+ Tier3B)
└── FIGURE_5_PANELS/         Supp_Panel_A..F + variant_sharing_tiers.tsv
```

---

## Environment and Relationships

All scripts run in the `NETWORK` conda environment (scanpy, scipy, scikit-learn, igraph, leidenalg, networkx, statsmodels, gseapy). Figure 4 provides the network partition, DE genes, chains, and three-group assignments used throughout. Figure 5 establishes that the network is not a patient artifact and that the patient-level differentiator is HPV16-associated variants, several inside the network, which Figure 6 follows up through the viral lifecycle.
