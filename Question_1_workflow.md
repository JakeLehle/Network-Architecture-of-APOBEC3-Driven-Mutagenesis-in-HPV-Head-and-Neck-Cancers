# Question 1 Workflow

## Question 1: What is the fundamental relationship between A3 expression and observed A3-induced mutations in HNSCC?

### Rationale

APOBEC3 (A3) cytidine deaminases with access to the nuclear compartment (A3A, A3B, A3C, A3H) are capable of editing genomic DNA, producing characteristic C>T and C>G mutations captured by the COSMIC SBS2 mutational signature. However, the presence of A3 enzyme expression alone does not guarantee that A3-driven mutagenesis will occur in a given tumor. This question establishes the foundational observation of the paper: **A3A and A3B expression is necessary but not sufficient for SBS2 mutagenesis in HNSCC**, implying the existence of unknown cofactors that regulate A3 enzymatic activity on genomic DNA. Additionally, this question dissects the individual and combinatorial contributions of A3A and A3B, establishing that A3B provides a baseline mutagenic state that saturates at high expression, while A3A continues to amplify mutagenesis without plateau.

Beyond the expression-SBS2 relationship, this question systematically excludes genetic explanations for the differential SBS2 burden. Germline analysis (pan-cancer, 7 cancer types) found no inherited variants in A3 coding regions or genome-wide that explain the HIGH/LOW split. Somatic analysis (HNSC, 53 vs 53 TMB-controlled) found that the mutations distinguishing HIGH from LOW tumors reflect consequences of high mutation burden (immune evasion, apoptosis resistance) and tumor subtype differences, not causal cofactors enabling A3 activity. Together, these negative results narrow the search space to **transcriptional regulation**, motivating the network analysis in Question 2.

### Key Results (v3, April 2026)

| Metric | Value |
|--------|-------|
| HNSC tumors matched (DIRECT only) | 426 |
| SBS2 range | 0 -- 728 mutations |
| SBS2 median | 18 |
| Tumors with SBS2 > 0 | 275 (64.6%) |
| Coral region (A3 + SBS2 HIGH) | ~213 |
| Cream region (A3 + SBS2 LOW) | ~213 |
| Teal region (no A3 + high SBS2) | 0 |
| Low-A3 outliers (A3A+A3B < 1.0, SBS2 > median) | 3 (0.7%) |
| Median-split quadrant medians | 5.5 / 18.0 / 24.0 / 20.5 |
| Both-high vs neither p-value | See diagnostic output |
| A3A-high-alone vs both-high p-value | ns (p=0.11, saturation signal) |
| SBS signature source | SigProfilerAssignment, COSMIC v3.4, exome mode |
| Germline enrichment (HNSC, 6938 tested) | 0 BH-significant variants |
| Pan-cancer germline (7 cancers, 728 patients) | 0 BH-significant A3 variants |
| Somatic enrichment (53 vs 53, TMB-adjusted) | CASP8 and HLA-A enriched in HIGH; 16 genes depleted |

### Data Sources

| Source | Description |
|--------|-------------|
| [GDC Data Portal](https://portal.gdc.cancer.gov/) | TCGA RNA-seq gene expression (STAR - Counts, FPKM-UQ) across 33 cancer types |
| [GDC Data Portal](https://portal.gdc.cancer.gov/) | TCGA MuTect2 Annotated Somatic Mutation MAFs across 33 cancer types (controlled access) |
| SigProfilerAssignment | COSMIC v3.4 SBS signature weights computed from PASS-filtered somatic SNVs |

### Directory Structure

```
scripts/TCGA/
├── DOWNLOAD_TCGA_COUNT_DATA.sh              # SLURM Job 1: RNA-seq data acquisition
├── Download_RNA-seq_Counts_TCGA.R           #   Step 1 -- Download from GDC
├── Organize_RNA-seq_Counts_TCGA.R           #   Step 2 -- Organize & build manifest
├── Master_TCGA_RNA-seq_Counts_Table.R       #   Step 3 -- Build master expression tables
│
├── ANALYZE_TCGA_COUNT_DATA.sh               # SLURM Job 2: Original analysis (Steps 4-8)
├── Prep_mutation_analysis_files.R           #   Step 4 -- Match A3 expression <-> SBS signatures
├── Patient_Level_HNSCC_TCGA_A3s_vs_SBS2.R  #   Step 5 (ORIGINAL) -- 6-panel Fig 1a
├── Diagnostic_A3C_A3H_bystander.R          #   Step 6 (ORIGINAL) -- ROC/AUC (Fig 1b)
├── Patient_Level_HNSCC_TCGA_3D_A3s.R       #   Step 7 (ORIGINAL) -- 3D scatter (Fig 1c)
├── A3A_A3B_additive_SBS2.R                 #   Step 8 (ORIGINAL) -- Additive model (Fig 1d-f)
│
├── Step05_Revised_HNSC_A3_vs_SBS2_v3.py    #   Step 5 (v3) -- 4-panel Fig 1 + supplemental
│                                            #     Uses new SigProfiler v3.4 SBS2 counts
│                                            #     DIRECT-only crosswalk (426 tumors)
│                                            #     Replaces Steps 5-8 R scripts for figure gen
│
├── Step05_Panel_1d_Saturation.py            #   Panel 1d -- A3B saturation fan + threshold sweep
│
├── HNSC_Somatic_Enrichment_Analysis.py      #   Somatic mutation enrichment (HIGH vs LOW)
│                                            #     Gene-level burden + TMB adjustment + KEGG
│                                            #     Uses Step03 groups from network pipeline
│
├── VCF/                                     # Branch B: Somatic mutation & SBS signature pipeline
│   ├── Diagnostic_TCGA_SNV_Availability.R   #   Step 0 -- Check GDC data availability
│   ├── Download_MuTect2_VCFs_TCGA.R         #   Step 1 -- Download MuTect2 MAFs (10,939 samples)
│   ├── Download_Pindel_VCFs_TCGA.R          #   Step 2 -- Download Pindel MAFs
│   ├── Organize_SNV_VCFs_TCGA.R             #   Step 3 -- Build unified manifest
│   ├── Consolidate_MAFs_TCGA.R              #   Step 4 -- Consolidate to pan-cancer master MAF
│   ├── Run_SigProfiler.py                   #   Step 5 -- Build SBS96 matrix + COSMIC assignment
│   ├── SIG_PROFILER_SCRIPT.sh               #   SLURM batch script for Steps 4-5
│   │
│   ├── Germline_SNP_Enrichment_Analysis.py          # HNSC-only germline (Tier 1+2)
│   ├── Germline_SNP_Enrichment_Analysis_Tier_1.py   # HNSC-only germline (Tier 1 only)
│   ├── Pan_Cancer_Germline_Feasibility_Diagnostic.py # Pan-cancer Phase 1: feasibility
│   └── Pan_Cancer_Germline_Enrichment.py             # Pan-cancer Phase 2: enrichment
│
└── TROUBLESHOOTING/
    ├── Diagnostic_Check_SBS2_Weight_Normalization.py   # Verified original SBS2 weights are
    │                                                    # absolute counts (CV=4.22), not ratios
    ├── Diagnostic_Compare_SBS_Weight_Sources.py        # Compared original vs new SigProfiler
    │                                                    # outputs (Spearman rho=0.90 for HNSC)
    ├── Diagnostic_Verify_Barcode_Matching.py           # Audited all 504 RNA-seq<->WES barcode
    │                                                    # pairs; verified 3 low-A3 outliers
    ├── Diagnostic_Figure1_Text_Numbers.py              # Extracted all text-ready numbers
    ├── Diagnostic_Verify_CONTEXT_Indexing.py           # Verified CONTEXT[5] == ref base (20/20)
    ├── Diagnostic_barcode_ambiguity_and_group_overlap.py # Compared 426 vs 502 group overlap
    ├── TCGA_troubleshooting.R
    ├── Compare_A3_Expression_Sources.R
    ├── FILTER_MASTER_TCGA_TSV.R
    ├── Global_TCGA_A3s_vs_SBS2.R
    ├── Patient_Level_HNSCC_TCGA_Network_Analysis.R
    ├── A3_Contribution_to_SBS2.R
    └── Diagnose_nonlinear_A3_Relationships.R
```

### Pipeline Overview

The analysis has three branches that converge:

```
 BRANCH A: RNA-seq Expression (SLURM Jobs 1-2)
 ──────────────────────────────────────────────
   Step 1-3: Download & build TCGA_master_FPKM_UQ.tsv
   Step 4: Extract A3 expression, match to mutation data
       |
       +---> TCGA_master_FPKM_UQ.tsv (11,505 samples)
            TCGA_sample_metadata_final.tsv (Project_ID)
            Mutation_Table_Tumors_TCGA.tsv (barcode crosswalk)


 BRANCH B: Somatic Mutation SBS Signatures (scripts/TCGA/VCF/)
 ──────────────────────────────────────────────────────────────
   Step 0: Diagnostic -- check GDC availability
   Step 1: Download MuTect2 MAFs (10,939 samples, 33 cancer types)
   Step 2: Download Pindel MAFs
   Step 3: Organize & verify files, build unified manifest
   Step 4: Consolidate per-sample MAFs -> pan-cancer master MAF
           + per-cancer MAFs + SigProfiler input files (PASS + ALL)
   Step 5: Build SBS96 matrix from CONTEXT column -> SigProfilerAssignment
       |
       +---> TCGA_SBS_signature_counts.tsv (absolute counts per sig)
            TCGA_SBS_signature_weights.tsv
            TCGA_MuTect2_master_manifest.tsv (Entity_ID -> Cancer_Type)
            per_cancer/TCGA-HNSC_mutations.maf.tsv (full annotated MAF)


 BRANCH C: Network Pipeline Group Selection (scripts/NETWORK/)
 ──────────────────────────────────────────────────────────────
   Step01: Clean expression matrix
   Step02: Merge expression + SBS signatures via crosswalk (UPDATED)
   Step03: Define HIGH/LOW groups (A3 median + SBS2 top/bottom 25%)
       |
       +---> TCGA-HNSC_SBS2_HIGH_group.pkl (53 tumors)
            TCGA-HNSC_SBS2_LOW_group.pkl (53 tumors)


 CONVERGENCE: Figure 1 + Enrichment Analyses
 ────────────────────────────────────────────
   Step05_Revised_HNSC_A3_vs_SBS2_v3.py
     Inputs: Branches A + B
     Output: 4-panel Figure 1, HNSC_A3_SBS2_matched_v3.tsv (426 tumors)

   Germline_SNP_Enrichment_Analysis.py (HNSC) +
   Pan_Cancer_Germline_Enrichment.py (7 cancers)
     Inputs: Branches A + B + per-cancer MAFs
     Output: 0 BH-significant germline variants

   HNSC_Somatic_Enrichment_Analysis.py
     Inputs: Branch C groups + per-cancer MAF + per-sample MAFs
     Output: CASP8/HLA-A enriched; NSD1/TP53 depleted; Apoptosis KEGG
       |
       v
   NARRATIVE: Genetic sequence (germline + somatic) cannot explain
   differential SBS2 -> cofactors must be transcriptional
   -> MOTIVATES QUESTION 2 (network analysis)
```

### Environments

| Environment | Scripts | Purpose |
|-------------|---------|---------|
| `RNA-seq_NovoGene` | Steps 1-4 (R), VCF Steps 0-4 (R) | TCGA data download, processing |
| `SComatic` | VCF Step 5 (Python) | SigProfilerAssignment |
| `NETWORK` | Revised Step 5, enrichment scripts, diagnostics | Figure generation, analysis |

---

### Revised Step 5 (v3): Figure 1 Generation (CURRENT)

**Script:** `Step05_Revised_HNSC_A3_vs_SBS2_v3.py`

The primary Figure 1 generation script. Replaces the original R scripts (Steps 5-8) with a single Python script. v3 uses DIRECT-only crosswalk matches (426 tumors), dropping all 76 CASE_ID matches which had ambiguous WES manifest entries.

**Change from v2 to v3:** After the barcode diagnostic (`Diagnostic_barcode_ambiguity_and_group_overlap.py`) showed all 76 CASE_ID matches have ambiguous WES manifest entries (duplicate rows) and were absent from the original analysis crosswalk file, v3 filters to DIRECT crosswalk matches only. This reduced the dataset from 502 to 426 tumors but provides verified RNA-to-WES barcode provenance for every sample.

**Dependencies:** `pandas`, `numpy`, `matplotlib`, `scipy`

**Input:**
- `TCGA_master_FPKM_UQ.tsv` -- full pan-cancer expression matrix (from Step 3)
- `TCGA_sample_metadata_final.tsv` -- Project_ID for cancer type filtering
- `Mutation_Table_Tumors_TCGA.tsv` -- barcode crosswalk (RNA-seq <-> WES), DIRECT matches only
- `TCGA_SBS_signature_counts.tsv` -- SigProfiler v3.4 absolute SBS counts (from VCF Step 5)
- `TCGA_MuTect2_master_manifest.tsv` -- HNSC WES Entity_IDs

**Output (-> `data/FIG_1/`):**
- `HNSC_A3_SBS2_matched_v3.tsv` -- 426 matched tumors (canonical dataset)
- `FIGURE_1_PANELS/Panel_1a_A3sum_vs_SBS2.pdf/.png`
- `FIGURE_1_PANELS/Panel_1b_A3A_vs_A3B_SBS2.pdf/.png`
- `FIGURE_1_PANELS/Panel_1c_Boxplot_Heatmap.pdf/.png`
- `FIGURE_1_PANELS/Supplemental_Low_A3_High_SBS2_Zoom.pdf/.png`
- `TROUBLESHOOTING/figure1_v3_pipeline_report.txt`

**Panel 1d (separate script):** `Step05_Panel_1d_Saturation.py` generates the A3B saturation fan chart and percentile threshold sweep, showing that A3B's contribution to SBS2 saturates (rho collapses around the 79th percentile) while A3A's does not (P90 keeps rising).

**v3 Key Numbers:**
- n = 426 tumors (DIRECT-only), quadrant medians: 5.5 / 18.0 / 24.0 / 20.5
- Both-high does NOT exceed A3A-high-alone (ns, p=0.11) -- the saturation signal
- A3B saturates, A3A does not

**Usage:**
```bash
conda run -n NETWORK python scripts/TCGA/Step05_Revised_HNSC_A3_vs_SBS2_v3.py
conda run -n NETWORK python scripts/TCGA/Step05_Panel_1d_Saturation.py
```

---

### Network Pipeline Step02 Update (April 2026)

**Script:** `scripts/NETWORK/Step02_Merge_SBS_Signatures.py` (UPDATED)

The network pipeline's Step02 was updated to use the new SigProfiler v3.4 counts with the DIRECT crosswalk merge strategy, replacing the old direct inner join on the legacy `Mutation_Table_Tumors_TCGA.tsv`.

**What changed:**
- Old: merged expression directly with `Mutation_Table_Tumors_TCGA.tsv` (unknown provenance, 65 COSMIC sigs, RNA-barcode-indexed) via inner join on Entity_ID
- New: loads SigProfiler v3.4 counts (WES-barcode-indexed), uses old file ONLY as crosswalk (RNA -> WES barcode mapping), DIRECT matches only

**Why:** Ensures the network pipeline (Steps 03-08) uses the same SBS2 values with documented provenance as Figure 1. The old pipeline got ~425 HNSC tumors with old SBS2 values; the updated pipeline gets 426 with v3.4 SBS2 values.

**Impact on downstream:** Steps 03-08 run unchanged. The output pkl has the same column structure (Entity_ID, Project_ID, SBS2, ENSG columns, A3 aliases). Two informational columns added (WES_Barcode, match_source) for audit trail.

---

### Germline SNP Enrichment Analysis

**Goal:** Test whether inherited germline variants in A3 coding regions or genome-wide explain why some A3-expressing tumors accumulate high SBS2 while others do not.

**HNSC-only scripts:**
- `Germline_SNP_Enrichment_Analysis.py` -- Tier 1 + Tier 2 combined
- `Germline_SNP_Enrichment_Analysis_Tier_1.py` -- Tier 1 only (high confidence)

**Pan-cancer scripts:**
- `Pan_Cancer_Germline_Feasibility_Diagnostic.py` -- Phase 1: feasibility across 33 cancers
- `Pan_Cancer_Germline_Enrichment.py` -- Phase 2: two-track enrichment across testable cancers

**Group selection:** Identical to network pipeline Step03 (above-median A3A+A3B, rank by SBS2, top/bottom 25%, equal-sized groups).

**Germline identification:** Two-tiered from MuTect2 FILTER column:
- Tier 1 (high confidence): `alt_allele_in_normal` flag
- Tier 2 (moderate confidence): `germline_risk` flag only

**Results (HNSC):**
- Track A (general): 6,938 recurrent germline SNPs tested, 93 nominal (p<0.05), 0 BH-significant
- Track B (A3-specific): 33 A3 gene variants cataloged, none significant

**Results (pan-cancer, 7 cancers, 728 patients):**
- Track A: 0 BH-significant variants in any cancer
- Track B: 187 unique A3 variants across 281 observations, none significant
- Suggestive patterns: A3H variants cluster in LOW (protective), A3D/A3B in HIGH (activity-enhancing), but underpowered

**Key limitation:** MuTect2 is a somatic caller; germline flags are incidental annotations, not a systematic germline survey. A definitive test would require dedicated germline calls (e.g., GATK HaplotypeCaller from MC3/dbGaP). This limitation is documented in the methods.

**Interpretation:** Common germline variation cannot explain differential SBS2 burden at detectable effect sizes.

---

### Somatic Mutation Enrichment Analysis (NEW, April 2026)

**Script:** `HNSC_Somatic_Enrichment_Analysis.py`

**Goal:** Test whether tumors with high SBS2 acquire different somatic mutations that could enable higher A3 activity. Uses MuTect2 for its intended purpose (somatic calling), addressing the primary limitation of the germline analysis.

**Groups:** Loaded from Step03 output pkls (53 HIGH + 53 LOW), ensuring identical patient sets with the network pipeline.

**Somatic variant selection:** FILTER == "PASS" from per-cancer MAF. Keeps SNPs and indels.

**Phases:**
1. TMB computation and comparison (HIGH vs LOW)
2. Track A: gene-level burden (non-silent, damaging/LOF, silent control) + TMB-adjusted logistic regression
3. Track B: A3-specific somatic variant catalog
4. APOBEC trinucleotide context from per-sample MAF CONTEXT column
5. Track C: gene set enrichment (Harris interactors, network communities, DDR, chromatin remodelers) + KEGG via Enrichr

**Key Results:**

*TMB confounding:* HIGH tumors have 1.5x more mutations (median 470 vs 232, p=4.7e-9). TMB adjustment is essential.

*Enriched in HIGH (TMB-adjusted):*
- CASP8 (15H/2L, Fisher p=0.001, LR p=0.008) -- apoptosis; loss allows heavily mutated cells to survive
- HLA-A (9H/2L, LR p=0.046) -- antigen presentation; loss enables immune evasion
- KEGG: Apoptosis only significant pathway (adj p=0.0006)

*Depleted in HIGH (TMB-adjusted, 16 genes):*
- NSD1 (4H/10L, LR p=0.001) -- histone H3K36 methyltransferase; loss defines a different HNSCC subtype
- TP53 (14H/25L in damaging track, p=0.04) -- more common in LOW, associated with NSD1-mutant subtype

*Silent control:* 1/311 nominal (AHNAK, giant gene). Confirms non-silent signal is not TMB artifact.

*A3 somatic:* 15 variants across 10 patients, scattered, no pattern. A3 genes are not somatic targets.

*APOBEC context:* 43.9% of all PASS SNPs at TCW motifs (expected for HNSC).

*Track C:* Community 10 nominally enriched (1 hit in 3 tested genes, p=0.026) but does not survive BH. No gene set reaches significance.

**Interpretation:** Somatic mutations distinguishing HIGH from LOW reflect consequences of high mutation burden (CASP8/HLA-A = survival and immune evasion) and tumor subtype differences (NSD1/TP53 depletion), NOT causal cofactors enabling A3 activity. This narrows the search to transcriptional regulation.

**Output (-> `data/FIG_1/SOMATIC_ENRICHMENT/`):**
- `HNSC_somatic_track_A_nonsilent.tsv` (with frac_apobec_context column)
- `HNSC_somatic_track_A_damaging.tsv`
- `HNSC_somatic_track_A_silent.tsv`
- `HNSC_somatic_track_A_tmb_adjusted.tsv` (with LR BH correction)
- `HNSC_somatic_track_B_a3_variants.tsv`
- `HNSC_somatic_track_C_geneset_enrichment.tsv`
- `HNSC_somatic_track_C_kegg_enrichr.tsv`
- `HNSC_somatic_tmb_comparison.tsv`
- `HNSC_somatic_apobec_context_by_gene.tsv`
- `HNSC_somatic_enrichment_report.txt`

**Usage:**
```bash
conda run -n NETWORK python scripts/TCGA/HNSC_Somatic_Enrichment_Analysis.py
```

---

### VCF Pipeline (scripts/TCGA/VCF/)

Independent pipeline that downloads all TCGA somatic mutation data and recomputes COSMIC SBS signature weights. This replaced a pre-processed signature file of unknown provenance.

**Scale:** 10,939 samples across 33 cancer types, 87.6 million total mutations, 6.58 million PASS-filtered somatic SNVs.

**Key output:** `TCGA_SBS_signature_counts.tsv` -- absolute mutation count attributions per COSMIC v3.4 signature per sample. Used as input to the revised Step 5 and the updated network pipeline Step02.

**CONTEXT column verification:** `Diagnostic_Verify_CONTEXT_Indexing.py` confirmed that the 11-mer CONTEXT column has the reference base at 0-indexed position [5] (1-indexed position 6). The trinucleotide extraction in `Run_SigProfiler.py` using positions [4],[5],[6] is correct (20/20 match in test).

**Comparison with original file (Diagnostic_Compare_SBS_Weight_Sources.py):**
- Original: 8,465 samples, 65 signatures (older COSMIC version), integer counts
- New: 10,936 samples, 86 signatures (COSMIC v3.4), integer counts
- HNSC-specific Spearman rho = 0.90 (strong rank agreement despite different decompositions)
- Decision: Use new counts for documented provenance and current COSMIC version

For full VCF pipeline documentation, see the `TCGA_Analysis` reference document.

---

### Original R Scripts (Steps 5-8, SUPERSEDED)

These R scripts produced the original 6-panel Figure 1. They are retained in the repository for reference but are superseded by `Step05_Revised_HNSC_A3_vs_SBS2_v3.py` for the current manuscript.

| Script | Original Panel | Status |
|--------|---------------|--------|
| `Patient_Level_HNSCC_TCGA_A3s_vs_SBS2.R` | Fig 1a (scatter) | Superseded by Step05 v3 Panel 1a |
| `Diagnostic_A3C_A3H_bystander.R` | Fig 1b (ROC) | Dropped per PI directive (4/2) |
| `Patient_Level_HNSCC_TCGA_3D_A3s.R` | Fig 1c (3D) | Dropped per PI directive (4/2) |
| `A3A_A3B_additive_SBS2.R` | Fig 1d-f (additive) | Superseded by Step05 v3 Panels 1b-c |

---

### Troubleshooting Scripts

| Script | Purpose | Key Finding |
|--------|---------|-------------|
| `Diagnostic_Check_SBS2_Weight_Normalization.py` | Check if original SBS2 weights are ratio-based | CV=4.22, absolute counts confirmed |
| `Diagnostic_Compare_SBS_Weight_Sources.py` | Compare original vs new SigProfiler outputs | HNSC rho=0.90, different COSMIC versions |
| `Diagnostic_Verify_Barcode_Matching.py` | Audit all RNA-seq<->WES barcode pairs | 502/502 pass, 3 outliers verified |
| `Diagnostic_Figure1_Text_Numbers.py` | Extract text-ready numbers for results section | All quadrant medians, p-values, region counts |
| `Diagnostic_Verify_CONTEXT_Indexing.py` | Verify CONTEXT column indexing in MAF files | Position [5] (0-indexed) == ref base, 20/20 confirmed |
| `Diagnostic_barcode_ambiguity_and_group_overlap.py` | Compare 426-only vs 502 groups, check CASE_ID ambiguity | All 76 CASE_ID matches ambiguous; HIGH Jaccard=0.82, LOW=0.80 |

---

### Figure 1 Summary (v3)

**Title:** A3A and A3B expression is necessary but not sufficient for SBS2 mutagenesis in HNSCC.

| Panel | Content | Script |
|-------|---------|--------|
| a | HNSC summed A3A+A3B vs SBS2 scatter with colored regions | `Step05_Revised_HNSC_A3_vs_SBS2_v3.py` |
| b | A3A vs A3B per tumor, colored by SBS2 (depth-sorted) | `Step05_Revised_HNSC_A3_vs_SBS2_v3.py` |
| c | Box-and-whisker + 2x2 heatmap of median-split quadrants | `Step05_Revised_HNSC_A3_vs_SBS2_v3.py` |
| d | A3B saturation fan + percentile threshold sweep | `Step05_Panel_1d_Saturation.py` |
| Supp | Zoomed view of low-A3/high-SBS2 boundary (4 tumors labeled) | `Step05_Revised_HNSC_A3_vs_SBS2_v3.py` |

**Narrative arc:**
1. A3A+A3B expression is necessary but not sufficient for SBS2 (Panel a)
2. Both A3A and A3B contribute, with SBS2 signal concentrated where both are high (Panel b)
3. A3B provides a baseline that A3A amplifies, with the highest burden requiring both (Panel c)
4. A3B saturates at high expression while A3A does not, explaining the "necessary but not sufficient" pattern (Panel d)
5. Germline variation (A3 coding + genome-wide) cannot explain differential SBS2 (germline analysis)
6. Somatic mutations in HIGH reflect consequences (CASP8/HLA-A) and subtype (NSD1/TP53), not causes (somatic analysis)
7. Unknown cofactors must regulate A3 activity at the transcriptional level -> **motivates Question 2 (network analysis)**
