# Figure 5 Walkthrough: Patient-Specific Effects

## Overview

Figure 5 asks whether the SBS2-HIGH epithelial cell population identified in Figure 4 is a generalizable transcriptional program across head and neck cancer patients or an artifact driven by one or two individuals. The analysis is organized into three tiers of increasing complexity: expression profiling, somatic variant analysis, and network sensitivity testing.

The central finding is that while three patients (SC029, SC013, SC001) contribute 67.6% of SBS2-HIGH cells, these patients are not transcriptionally distinct. Their HIGH cells intermingle with other patients on PCA (silhouette = 0.125), and the co-expression network survives leave-one-patient-out removal with the A3 wall intact at 100%. What does distinguish high-contributor patients is somatic mutations in HPV infection and antigen presentation genes (HLA-B, HLA-C, HLA-DRB1, TAP1) that overlap with the same immune pathways enriched in the Figure 4 network. This motivates the HPV lifecycle analysis in Figures 6-7.

---

## Directory Structure

```
scripts/PATIENT_SPECIFIC_EFFECTS/
├── patient_config.py                              # Shared configuration (paths, parameters, utilities)
├── RUN_PATIENT_ANALYSIS.sh                        # SLURM master script (runs all tiers overnight)
│
├── Tier1A_Patient_Expression_GSEA.py              # Per-patient Wilcoxon DE + KEGG GSEA
├── Tier1B_HIGH_Cell_Transcriptional_Similarity.py # PCA, silhouette, profile correlations
├── Tier1C_Enriched_vs_Depleted_Patient_DE.py      # High-contributor vs other patient DE
│
├── Tier2A_Expression_Haplotype_Proxies.py         # A3 ratios, activating chain expression
├── Tier2B_SNP_Pattern_Analysis.py                 # SComatic variant sharing across patients
│
├── Tier3A_Leave_One_Patient_Out.py                # LOPO network sensitivity (V4 pipeline)
│
├── Generate_Supplemental_Patient_Effects.py       # Supplemental figure panels A-F
├── Analyze_SNP_Tier_Genes.py                      # Variant-to-gene mapping + KEGG per tier
└── Diagnostic_HC_Network_HPV_Overlap.py           # HC-exclusive x network x HPV/antigen overlap
```

---

## Input Data

All inputs come from the Figure 4 single-cell network pipeline. The shared configuration file (`patient_config.py`) defines paths to these resources.

| Input | Path | Description |
|-------|------|-------------|
| AnnData | `data/FIG_4/00_input/adata_final.h5ad` | 155,650 cells, 27,736 genes (ClusterCatcher output) |
| Signature weights | `data/FIG_4/00_input/signature_weights_per_cell.txt` | Per-cell SBS signature weights (SigProfiler) |
| Three-group assignments | `data/FIG_4/01_group_selection/three_group_assignments.tsv` | 546 SBS2-HIGH, 546 CNV-HIGH, 546 NORMAL |
| DE genes | `data/FIG_4/NETWORK_SBS2_VS_NORMAL/02_differential_expression/` | Scanpy Wilcoxon DE results (FDR < 0.05) |
| Network partition | `data/FIG_4/NETWORK_SBS2_VS_NORMAL/04_communities/SC_best_partition.csv` | 2,948 genes in Leiden communities |
| Harris interactors | `data/FIG_4/00_input/Harris_A3_interactors.txt` | 175 known A3 protein interactors |
| SComatic calls | `results_NMF_v0.1.1/all_samples.single_cell_genotype.filtered.tsv` | Per-cell somatic variant calls |
| GTF | `SC/ref/GRCh38/genes/genes_unzipped.gtf` | GRCh38 gene annotations for variant mapping |

---

## Pipeline Execution

The entire pipeline runs as a single SLURM job. Typical runtime is 4-6 hours on 32 cores with 500 GB memory.

```bash
cd scripts/PATIENT_SPECIFIC_EFFECTS/
sbatch RUN_PATIENT_ANALYSIS.sh
```

The SLURM script runs three phases sequentially:

- **Phase 1** (~45 min): Tier 1A, 1B, 1C, 2A (expression analyses)
- **Phase 2** (~30-60 min): Tier 2B (SNP pattern analysis)
- **Phase 3** (~3-4 hrs): Tier 3A x3 (LOPO network sensitivity)

Post-pipeline scripts are run interactively:

```bash
conda run -n NETWORK python Generate_Supplemental_Patient_Effects.py
conda run -n NETWORK python Analyze_SNP_Tier_Genes.py
conda run -n NETWORK python Diagnostic_HC_Network_HPV_Overlap.py
```

---

## Script Descriptions

### Shared Configuration

**`patient_config.py`** defines all paths, parameters, and utility functions imported by every script.

Key design choices:
- Three-group loading via `load_three_groups()` returns SBS2-HIGH, CNV-HIGH, and NORMAL cell sets (546 each) from `three_group_assignments.tsv`
- Backward-compatible `load_groups()` wrapper returns (SBS2-HIGH, NORMAL) for scripts that only need the SBS2_VS_NORMAL comparison
- Network parameters match the V4 pipeline: scanpy DE at FDR < 0.05, DIFF threshold auto-selection via max fragmentation rate, full-network Leiden with resolution sweep and component-aware merge
- Gene sets from the Figure 4 concordance analysis: 9-gene activating chain (RALY, HNRNPA2B1, CCL20, KRT24, LCN2, LINC00278, RRAD, SMOX, UTY), 3 inhibiting chain anchors (SNHG3, THYN1, ZNG1A), 2 A3 interactor anchors (RALY, HNRNPA2B1)
- Style constants: font sizes 28-34, hex color codes, PDF + PNG at 300 DPI

---

### Tier 1: Expression Analyses

**`Tier1A_Patient_Expression_GSEA.py`** compares each patient's basal cells against all other patients' basal cells (Wilcoxon rank-sum), then runs GSEA prerank against KEGG pathways. This is an unbiased screen for patient-specific pathway activity. Also computes per-patient A3 family expression and correlates with SBS2-HIGH fold enrichment.

**`Tier1B_HIGH_Cell_Transcriptional_Similarity.py`** takes the 546 SBS2-HIGH cells and asks whether they cluster by patient or intermingle. Runs PCA on DE genes (3,907 at FDR < 0.05), computes silhouette score using patient as label, and builds a per-patient mean profile correlation matrix with hierarchical clustering. A low silhouette score (0.125) indicates a shared program rather than patient-specific clusters.

**`Tier1C_Enriched_vs_Depleted_Patient_DE.py`** compares all basal cells from the three high-contributing patients (SC029, SC013, SC001) against all basal cells from the remaining 11 patients. Asks whether high contributors are constitutively different at the transcriptional level, independent of SBS2 status.

---

### Tier 2: Variant and Haplotype Analysis

**`Tier2A_Expression_Haplotype_Proxies.py`** measures per-patient expression-level proxies for germline A3 haplotype variation: A3A/A3B ratio (deletion polymorphism proxy), A3H expression (haplotype stability proxy), and the activating chain gene expression program from Figure 4. Tests whether high-contributor patients show elevated baseline expression of the mutagenic co-expression program.

**`Tier2B_SNP_Pattern_Analysis.py`** maps SComatic single-cell somatic variant calls to patients and classifies variants into sharing tiers (Universal, Broadly shared, Partially shared, Patient-specific, HC-exclusive). Computes pairwise Jaccard similarity between patients, identifies germline candidates on the A3 locus, and generates per-patient variant composition profiles.

---

### Tier 3: Network Sensitivity

**`Tier3A_Leave_One_Patient_Out.py`** is the core sensitivity test. For each high-contributing patient (SC029, SC013, SC001), removes that patient's cells from both SBS2-HIGH and NORMAL groups, then re-runs the full V4 network pipeline: scanpy DE (Wilcoxon, BH-FDR < 0.05), Spearman correlation matrices, max fragmentation rate threshold auto-selection, full-network Leiden with resolution sweep (5 resolutions x 15 runs, composite score selection), and component-aware merge.

Four sensitivity metrics are computed per LOPO run:

| Metric | Definition | Result |
|--------|------------|--------|
| Activating chain recovery | Number of 9 concordant chain genes present in LOPO network | 8/9 (SC029, SC001), 3/9 (SC013) |
| Community structure ARI | Adjusted Rand index between LOPO and full community assignments | 0.21-0.34 |
| Gene overlap Jaccard | Jaccard index between LOPO and full node sets | 0.45-0.54 |
| A3 wall integrity | % negative DIFF edges for A3A | 100% (all runs) |

---

### Post-Pipeline Analysis

**`Generate_Supplemental_Patient_Effects.py`** produces the six-panel supplemental figure:

| Panel | Content |
|-------|---------|
| A | Patient distribution of SBS2-HIGH cells (sorted by HIGH count) |
| B | Per-patient A3A and A3B expression (bar charts) |
| C | PCA of HIGH cells colored by patient + correlation heatmap |
| D | LOPO sensitivity (4-metric bar charts: chain recovery, ARI, Jaccard, wall) |
| E | Variant sharing spectrum (5-tier stacked bar per patient) |
| F | Genome ideogram with variant positions colored by sharing tier |

**`Analyze_SNP_Tier_Genes.py`** maps variants to genes using the GRCh38 GTF (79.8% mapping rate) and cross-references per-tier gene lists against the SC network partition, Harris A3 interactors, A3 enzymes, activating chain genes, and inhibiting chain anchors. Runs KEGG enrichment per tier.

**`Diagnostic_HC_Network_HPV_Overlap.py`** closes the loop between the network and patient-level findings. Cross-references HC-exclusive variant genes that are in significant HPV/immune KEGG pathways against the SC network partition, identifying which specific pathway genes sit inside the co-expression network and which community they belong to.

---

## Output Structure

```
data/FIG_5/
├── 00_diagnostics/
│   └── patient_enrichment_SBS2_HIGH_v2.tsv
├── 01_patient_expression/
│   ├── patient_KEGG_NES_heatmap.pdf/png
│   ├── patient_A3_expression_summary.tsv
│   ├── patient_GSEA_top_pathways.tsv
│   ├── per_patient_GSEA_results/
│   ├── HIGH_cell_PCA_by_patient.pdf/png
│   ├── HIGH_cell_silhouette_score.txt
│   ├── HIGH_cell_patient_correlation_heatmap.pdf/png
│   ├── Tier1C_volcano_high_contrib_vs_other.pdf/png
│   └── Tier1C_high_contrib_KEGG_GSEA.tsv
├── 02_snp_haplotype/
│   ├── Tier2A_haplotype_proxies.tsv
│   ├── Tier2A_haplotype_proxies.pdf/png
│   ├── Tier2A_activating_chain_heatmap.pdf/png
│   ├── Tier2B_patient_variant_summary.tsv
│   ├── Tier2B_patient_variant_comparison.pdf/png
│   ├── Tier2B_germline_candidates.tsv
│   └── gene_analysis/
│       ├── SNP_tier_gene_analysis_report.tsv
│       ├── variant_to_gene_mapping.tsv
│       └── SNP_tier_{tier}_KEGG.tsv
├── 03_sensitivity/
│   ├── LOPO_SC029/
│   │   ├── LOPO_SC029_summary.tsv
│   │   ├── LOPO_SC029_communities.tsv
│   │   ├── LOPO_SC029_DE_all.tsv
│   │   ├── LOPO_SC029_DE_significant.tsv
│   │   ├── LOPO_SC029_threshold_sweep.tsv
│   │   ├── LOPO_SC029_resolution_sweep.tsv
│   │   └── LOPO_SC029_parameters.txt
│   ├── LOPO_SC013/  (same structure)
│   └── LOPO_SC001/  (same structure)
└── FIGURE_5_PANELS/
    ├── Supp_Panel_A_Patient_Distribution.pdf/png
    ├── Supp_Panel_B_A3_Expression.pdf/png
    ├── Supp_Panel_C_PCA_Correlation.pdf/png
    ├── Supp_Panel_D_LOPO_Sensitivity.pdf/png
    ├── Supp_Panel_E_Variant_Sharing.pdf/png
    ├── Supp_Panel_F_Ideogram.pdf/png
    └── variant_sharing_tiers.tsv
```

---

## Key Findings

1. **SBS2-HIGH cells are a shared program, not patient-specific.** Three patients contribute 67.6% of HIGH cells, but silhouette = 0.125 and mean inter-patient correlation = 0.775 indicate convergent transcriptional phenotypes.

2. **A3 expression does not predict patient contribution.** The highest A3A expressors (SC005, SC003) are not the highest HIGH-cell contributors, reinforcing that expression alone is insufficient.

3. **The network is robust to patient removal.** LOPO analysis preserves the A3 wall at 100%, recovers 8/9 activating chain genes (SC029, SC001), and maintains core community structure across all configurations.

4. **High-contributor patients carry mutations in HPV/immune genes, not A3 genes.** HC-exclusive variants hit zero A3 family members, zero activating chain genes, and zero known A3 interactors. Instead, they enrich for HPV infection (p = 0.007) and antigen presentation (p = 0.031).

5. **HC-exclusive HPV/immune pathway genes overlap with the SC network.** 15 of 33 HC-exclusive HPV/immune pathway genes are present in the Figure 4 network, including HLA-B, HLA-C (MHC class I, community 2), HLA-DRB1 (MHC class II, community 1), and TAP1 (peptide transporter, community 0). Mutations in antigen presentation machinery could impair immune clearance of HPV-infected cells, allowing them to persist and accumulate APOBEC-driven mutations.

---

## Environment

All scripts run in the `NETWORK` conda environment. Key dependencies: scanpy, scipy, scikit-learn, matplotlib, igraph, leidenalg, networkx, statsmodels, gseapy. The SLURM script checks all dependencies before execution and installs gseapy if missing.

---

## Relationship to Other Figures

- **Figure 4** provides the network partition, DE genes, concordant chains, and three-group assignments used throughout
- **Figure 5** validates that the network is not a patient-specific artifact and identifies HPV/immune mutations as the patient-level differentiator
- **Figures 6-7** follow up on the HPV infection signal, investigating the viral lifecycle framework and neoantigen landscape motivated by the HC-exclusive variant enrichment
