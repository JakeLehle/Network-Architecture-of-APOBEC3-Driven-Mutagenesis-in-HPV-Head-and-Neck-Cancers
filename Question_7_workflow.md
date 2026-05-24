# NEOANTIGEN Pipeline Walkthrough

## Overview

This directory contains the Figure 7 neoantigen landscape analysis pipeline. The pipeline starts from SComatic single-cell somatic variant calls, annotates protein effects with SnpEff, predicts MHC-I binding with MHCflurry using real proteome context, ranks neoantigens by expression-weighted composite scores, analyzes RNA fusion events, and integrates immune escape evidence into a therapeutic tier classification.

The central finding: ANXA1 is the top neoantigen target in SBS2-HIGH cells (IC50 = 13.1 nM, 42 predicted neoantigen peptides, expressed in 99.5% of cells). The neoantigen landscape partitions into three clinically actionable therapeutic tiers based on immune escape evidence.

## Pipeline Flow

```
Step01 (NETWORK)          Step02 (NEOANTIGEN)        Step03 (NEOANTIGEN)
SComatic TSV              SnpEff annotation           MHCflurry binding
  + three-group             + germline                  + proteome-based
  assignments               subtraction                 peptide generation
       |                        |                           |
       v                        v                           v
  Per-group VCFs          Somatic protein-            Neoantigen calls
  + barcode lists         altering variants           + peptide details
  + config YAML                                       + binding summary
                                                           |
                               +---------------------------+
                               |                           |
                               v                           v
                    Step04 (NETWORK)              Step05 (NETWORK)
                    Expression-weighted           STAR chimeric
                    composite scoring             junction analysis
                    + group-unique rankings       + fusion enrichment
                               |                           |
                               +---------------------------+
                               |
                               v
                    Step06 (NETWORK)
                    Antigen loss + expression
                    silencing + immune escape
                    integration + vaccine scoring
                    + therapeutic tier classification
                               |
                               v
                    Generate_Figure7_Panels.py (NETWORK)
                    4-panel figure: burden, Venn,
                    ANXA1 peptides, ANXA1 gene track
```

## Conda Environments

| Environment | Scripts | Key packages |
|-------------|---------|-------------|
| NETWORK | Steps 01, 04, 06, Figure generation | scanpy, pandas, matplotlib, networkx, gseapy |
| NEOANTIGEN | Steps 02, 03 | SnpEff (java), mhcflurry, tensorflow |

The SLURM wrapper (`RUN_NEOANTIGEN_PIEPLINE.sh`) handles environment switching automatically between steps.

## Pipeline Scripts

### Step01_Prep_Neoantigen_Inputs.py

Loads the three-group cell assignments (`three_group_assignments.tsv`) and the master SComatic filtered TSV, splits somatic variants into per-group VCFs (SBS2_HIGH, CNV_HIGH, NORMAL), and generates barcode lists and a pipeline config YAML. Writes VCF with 1-based coordinates (SComatic `Start` is 0-based, Step01 adds +1 for VCF `POS`). Also computes per-group variant summaries including APOBEC TCW context counts.

Outputs to `data/FIG_7/01_neoantigen_inputs/`.

### Step02_SnpEff_Annotation.py

Runs SnpEff on per-group VCFs to classify variant effects (missense, synonymous, intronic, etc.). Parses the ANN field to extract gene, effect, HGVS protein notation, and transcript ID. Performs germline subtraction using NORMAL group variants as background: any variant present in NORMAL is flagged as germline and excluded from the disease group somatic calls. Produces filtered somatic protein-altering variant lists for each disease group.

Outputs to `data/FIG_7/02_snpeff_annotation/`.

### Step03_MHCflurry_Binding.py

The core neoantigen prediction step. For each somatic missense variant, generates mutant and wild-type peptides (lengths 8, 9, 10, 11) using actual protein sequence context from the Ensembl GRCh38 release 115 reference proteome (`pep.all.fa`, 245,535 isoforms). This replaces the earlier poly-alanine flanking approach, which padded peptides with alanine and gave systematically inaccurate IC50 values.

Protein sequence lookup uses a six-layer chain to maximize mapping success (98.6% for SBS2-HIGH):

1. ENST transcript ID (exact isoform from SnpEff annotation)
2. Gene symbol (canonical or longest isoform)
3. Gene alias (handles HUGO symbol updates: C4orf3 to C4orf33, TMEM199 to VMA12, SLC9A3R1 to NHERF1)
4. ENSG gene ID
5. Isoform scan (all isoforms at exact position)
6. Offset scan (positions +/-30, recovers signal peptide numbering shifts)

The offset scan is critical: ANXA1 was completely missing from the first two pipeline runs because SnpEff protein positions differ from the Ensembl proteome numbering by -30 due to signal peptide inclusion. The offset scan recovered ANXA1 and 76 other variants.

Peptides are scored against a 10-allele reference HLA panel (~80% population coverage) using MHCflurry Class1AffinityPredictor. Binders are defined as IC50 < 500 nM, strong binders as IC50 < 50 nM, and differential neoantigens as mutant IC50 < 500 nM with wild-type IC50 > 500 nM (new epitopes created by the mutation).

Outputs to `data/FIG_7/03_mhc_binding/`.

### Step04_Expression_Weighted_Ranking.py

Ranks neoantigens by an expression-weighted composite score: `composite = mean_expression * (500 / best_IC50)`, where expression is log1p-normalized from scanpy (no additional log transform). This naturally deprioritizes immunoglobulin genes that have high APOBEC motif density but low epithelial expression. Produces per-group rankings, group-unique neoantigen lists, and cross-group comparisons.

Outputs to `data/FIG_7/03_mhc_binding/` (expression-weighted files alongside binding results).

### Step05_Fusion_Analysis.py

Re-parses STAR chimeric junction output across all cells, filters to high-confidence fusions, and computes per-group fusion rates. Identifies group-exclusive fusion gene pairs and runs pathway enrichment (GO, KEGG, Reactome) using gseapy with exponential backoff retry logic for Enrichr 429 rate limiting. Cross-references fusion-disrupted genes with the neoantigen gene lists to identify asymmetric immune escape (e.g., SBS2-exclusive fusions hitting B2M and HLA-B antigen presentation genes, CNV-exclusive fusions hitting PSMB8/PSME2 proteasome genes).

Key finding: fusion rates are similar across groups (SBS2: 9.3/cell, CNV: 8.2/cell, NORMAL: 12.1/cell). Fusions do not differentiate groups; specific fusion partners do.

Outputs to `data/FIG_7/04_fusion_analysis/`.

### Step06_Integrated_Neoantigen_Analysis.py

Integrates neoantigen, fusion, and expression data into a unified immune escape analysis. Identifies antigen loss mutations (somatic variants in neoantigen genes that destroy the neoantigen), expression silencing (neoantigen genes > 2-fold downregulated in the opposite group), and multi-mechanism escape (genes with antigen loss + silencing + fusion disruption). Computes vaccine target scores with a 1.5x breadth bonus for shared targets and 25% per escape mechanism (these weights are interpretive, not statistical).

Classifies all neoantigen genes into four therapeutic tiers:

- Tier 1A (22 genes): shared neoantigens with CNV escape evidence (hot tumor targets, proven immune recognition)
- Tier 1B (276 genes): SBS2-only neoantigens (maintenance-phase targets, ANXA1 leads)
- Tier 2 (135 genes): CNV-only neoantigens (productive-phase targets + checkpoint blockade)
- Tier 3 (83 genes): shared neoantigens, no escape evidence (broad coverage vaccine backbone)

Outputs to `data/FIG_7/05_summary/` and `data/FIG_7/05_summary/THERAPEUTIC_TIERS/`.

### Generate_Figure7_Panels.py

Generates the four-panel Figure 7:

- Panel A: Neoantigen and RNA fusion burden comparison (SBS2-HIGH vs CNV-HIGH, total counts). Neoantigens show 1.8x enrichment in SBS2-HIGH; fusions are similar across groups.
- Panel B: Venn diagram of neoantigen gene overlap (381 SBS2, 240 CNV, 105 shared) with tier breakout annotations.
- Panel C: ANXA1 top neoantigen peptides grouped by two hotspot regions (pos 289 and pos 321-327), sorted by IC50 drop (WT minus mutant). Mutated amino acid highlighted in each peptide sequence.
- Panel D: ANXA1 protein gene track (full width) with domain annotations, lollipop mutations colored by hotspot region, and dashed highlight boxes connecting to Panel C.

Individual panels and combined 2x2 figure saved as PDF and PNG at 300 DPI.

### RUN_NEOANTIGEN_PIEPLINE.sh

SLURM wrapper that runs all six steps sequentially with conda environment switching, fail-fast error handling, and pre-flight checks (verifies input files exist before each step).

## Diagnostic Scripts (TROUBLE_SHOOTING/)

### Diagnostic_Proteome_Mapping.py

Investigated the initial 25% AA mismatch rate between SnpEff protein positions and the Ensembl proteome. Identified three causes: missing isoforms (fixed by isoform scan), gene alias mismatches (3 genes, fixed by alias table), and signal peptide numbering offsets (77 variants recovered by offset scan +/-30). Reduced mismatch rate from 25% to 0.5%.

### Diagnostic_Therapeutic_Tiers.py

Standalone tier classification from pipeline outputs. Used to validate tier assignments before integrating into Step06. Confirmed clean mapping: Tier 1A = shared genes with escape, Tier 1B = SBS2-only, Tier 2 = CNV-only, Tier 3 = shared without escape.

### Diagnostic_Figure7_Redesign.py

Pulled all data needed for the Figure 7 redesign: neoantigen counts per group, fusion rates per group, Venn diagram gene sets with tier cross-reference, ANXA1 peptide IC50 shifts by protein position, and mutation-to-domain mapping.

### Diagnostic_ANXA1_TCW_Context.py

Verified APOBEC trinucleotide context for all 8 ANXA1 somatic mutations. Identified a 0-based (SComatic) vs 1-based (VCF) coordinate offset that was causing incorrect trinucleotide lookups. With the correct offset applied, all 8 ANXA1 mutations are non-APOBEC (zero are in TCW context). Two of 8 have the right substitution type (C>G) but neither is in a TCW motif.

### Diagnostic_ANXA1_Provenance.py

Traced ANXA1 variant coordinates through the pipeline (SComatic TSV, per-group VCF, SnpEff annotation, somatic protein-altering output) to diagnose why 5 of 8 positions were not found in the SComatic master file. Confirmed the 0-based to 1-based coordinate conversion in Step01's `write_vcf` function (`pos = int(row['Start']) + 1`).

## Backup Scripts (BACKUP/)

Contains 14 archived scripts from the earlier Phase5B pipeline architecture before the six-step refactoring. Includes the original monolithic neoantigen pipeline, the initial STAR chimeric analysis scripts, and early troubleshooting scripts. Retained for provenance.

## Key Output Files

| File | Location | Description |
|------|----------|-------------|
| pipeline_config.yaml | data/FIG_7/01_neoantigen_inputs/ | All paths, parameters, group sizes |
| {group}.somatic_protein_altering.tsv | data/FIG_7/02_snpeff_annotation/ | Germline-subtracted missense variants |
| {group}_neoantigens.tsv | data/FIG_7/03_mhc_binding/ | All predicted binders with peptide details |
| {group}_all_peptide_results.tsv | data/FIG_7/03_mhc_binding/ | Full peptide prediction audit trail |
| {group}_expression_weighted_ranking.tsv | data/FIG_7/03_mhc_binding/ | Composite-scored rankings |
| all_filtered_junctions.tsv | data/FIG_7/04_fusion_analysis/ | Per-cell fusion events with gene pairs |
| per_group_junction_summary.tsv | data/FIG_7/04_fusion_analysis/ | Fusion rates per group |
| vaccine_target_scores.tsv | data/FIG_7/05_summary/ | Integrated scoring with escape bonuses |
| all_genes_tiered.tsv | data/FIG_7/05_summary/THERAPEUTIC_TIERS/ | Complete tier assignments |
| Figure7_Combined.pdf | data/FIG_7/FIGURE_7_PANELS/ | Combined 4-panel figure |

## Reference Data

- Ensembl GRCh38 release 115 proteome: `data/reference/Homo_sapiens.GRCh38.pep.all.fa` (245,535 protein isoforms, 19,885 gene symbols). The `pep.canonical.fa` file does not exist for release 115; `pep.all.fa` is the correct file. Selection is by longest isoform (Ensembl_canonical tag not present in release 115 headers).
