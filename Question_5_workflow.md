# Question 5 Workflow: Neoantigen Landscape (Figure 7)

## Question 5: What therapeutic neoantigen tiers fall out of the immune-visible and immune-evasive populations?

### Rationale

Question 4 (Fig 6) split the epithelial compartment into an immune-visible, maintenance-phase SBS2-HIGH state and an immune-evasive, productive-phase CNV-HIGH state. Because tumor-specific neoantigens are attractive targets for multivalent mRNA vaccines, and because how well that strategy works depends on the tumor's immune state, the neoantigens these two populations produce are relevant to different therapeutic settings. This question builds the neoantigen landscape from the single-cell somatic calls, predicts MHC-I binding, integrates immune-escape evidence, and ranks the result into therapeutic tiers matched to a tumor's position on the lifecycle axis.

The central finding: *ANXA1* is the top neoantigen target in SBS2-HIGH cells (best IC50 = 13.1 nM, 42 predicted peptides, expressed in 99.5% of cells), and the landscape partitions into three priority tiers by breadth of coverage and immune-escape evidence.

## Pipeline Flow

```
Step01 (NETWORK)          Step02 (NEOANTIGEN)        Step03 (NEOANTIGEN)
SComatic TSV              SnpEff annotation           MHCflurry binding
  + three-group             + germline                  + real-proteome
  assignments               subtraction                 peptide generation
       |                        |                           |
       v                        v                           v
  Per-group VCFs          Somatic protein-            Neoantigen calls
  + barcode lists         altering variants           + peptide details
                                                           |
                               +---------------------------+
                               |                           |
                               v                           v
                    Step04 (NETWORK)              Step05 (NETWORK)
                    Expression-weighted           STAR chimeric
                    composite scoring             fusion analysis
                               |                           |
                               +---------------------------+
                               |
                               v
                    Step06 (NETWORK)
                    Immune-escape integration
                    + therapeutic tier classification
                               |
                               v
                    Generate_Figure7_Panels.py (NETWORK)
                    4 panels: burden, Venn/tiers, ANXA1 peptides, ANXA1 track
```

## Conda Environments

| Environment | Scripts | Key packages |
|-------------|---------|--------------|
| NETWORK | Steps 01, 04, 05, 06, figure generation | scanpy, pandas, matplotlib, gseapy |
| NEOANTIGEN | Steps 02, 03 | SnpEff (java), mhcflurry, tensorflow |

`RUN_NEOANTIGEN_PIEPLINE.sh` switches environments automatically between steps.

## Pipeline Scripts (main directory)

### Step01_Prep_Neoantigen_Inputs.py
Loads the three-group assignments and the master SComatic filtered TSV, splits somatic variants into per-group VCFs (SBS2_HIGH, CNV_HIGH, NORMAL), and writes barcode lists and a config YAML. VCFs use 1-based coordinates (SComatic `Start` is 0-based, so `POS = Start + 1`). Also reports per-group APOBEC TCW-context counts. Outputs to `data/FIG_7/01_neoantigen_inputs/`.

### Step02_SnpEff_Annotation.py
Runs SnpEff on the per-group VCFs, parses the ANN field (gene, effect, HGVS protein notation, transcript), and performs germline subtraction using the NORMAL group as background: any variant present in NORMAL is dropped from the disease groups. Produces the germline-subtracted protein-altering variant lists. Outputs to `data/FIG_7/02_snpeff_annotation/`.

### Step03_MHCflurry_Binding.py
The core prediction step. For each somatic missense variant it builds mutant and wild-type peptides (lengths 8-11) from real protein context in the Ensembl GRCh38 release 115 proteome (`pep.all.fa`, 245,535 isoforms), replacing the earlier poly-alanine flanking that gave systematically wrong IC50 values. Protein lookup uses a six-layer chain (ENST, gene symbol, alias, ENSG, isoform scan, and a +/-30 offset scan), reaching 98.6-98.7% mapping. The offset scan matters: *ANXA1* was missing from the first two runs because SnpEff positions differ from the proteome by -30 due to signal-peptide numbering, and the scan recovered *ANXA1* and 76 other variants. Peptides are scored against a 10-allele HLA panel (~80% population coverage) with MHCflurry `Class1AffinityPredictor`; binders are IC50 < 500 nM, strong binders < 50 nM, and differential neoantigens have mutant < 500 nM with wild-type > 500 nM. Outputs to `data/FIG_7/03_mhc_binding/`.

### Step04_Expression_Weighted_Ranking.py
Ranks neoantigens by an expression-weighted composite, `composite = mean_expression * (500 / best_IC50)`, where expression is log1p-normalized. This deprioritizes immunoglobulin genes that have high APOBEC-motif density but low epithelial expression. Produces per-group and group-unique rankings. Outputs alongside the binding results in `data/FIG_7/03_mhc_binding/`.

### Step05_Fusion_Analysis.py
Re-parses STAR chimeric junctions across all cells, filters to high-confidence fusions, and computes per-group fusion rates. Identifies group-exclusive fusion pairs and runs pathway enrichment (gseapy with backoff for Enrichr rate limiting), then cross-references fusion-disrupted genes against the neoantigen lists to find asymmetric escape (for example SBS2-exclusive fusions hitting *B2M* and *HLA-B*, CNV-exclusive fusions hitting *PSMB8* and *PSME2*). Fusion rates are similar across groups, so fusions do not differentiate the populations; specific partners do. Outputs to `data/FIG_7/04_fusion_analysis/`.

### Step06_Integrated_Neoantigen_Analysis.py
Integrates neoantigen, fusion, and expression data into the immune-escape analysis: antigen-loss mutations (variants that destroy a neoantigen), expression silencing (neoantigen genes downregulated more than two-fold in the opposite group), and multi-mechanism escape. It computes vaccine target scores (a 1.5x breadth bonus for shared targets and 25% per escape mechanism, interpretive weights, not statistical) and classifies the 516 neoantigen-producing genes into four internal categories that the manuscript presents as three tiers (see below). Outputs to `data/FIG_7/05_summary/` and `.../THERAPEUTIC_TIERS/`.

### Generate_Figure7_Panels.py
Builds the four Figure 7 panels: A, neoantigen and fusion burden (SBS2-HIGH vs CNV-HIGH); B, the neoantigen-gene Venn with tier annotations; C, the *ANXA1* top peptides grouped by hotspot region and sorted by IC50 drop; D, the *ANXA1* gene track with domains and lollipop mutations colored by hotspot, linked to Panel C. Saved as PDF and PNG at 300 DPI.

### diagnostic_section7_numbers.py
Section 4.5 number-audit script: recomputes the burden, binder, TCW-context, fusion, tier, and *ANXA1* numbers from the pipeline outputs so each quantitative claim in the text is confirmable against one source.

### RUN_NEOANTIGEN_PIEPLINE.sh
SLURM wrapper that runs the six steps sequentially with environment switching, fail-fast handling, and pre-flight input checks.

## Therapeutic Tiers

Step06 computes four internal categories, which the manuscript collapses into three priority tiers:

| Manuscript tier | Genes | Internal categories | Lead | Setting |
|-----------------|-------|---------------------|------|---------|
| Tier 1 (highest) | 105 shared | 1A_hot_shared_escaped (22) + 3_broad_coverage (83) | *MDK* | Broadly conserved across tumor states |
| Tier 2 | 276 SBS2-specific | 1B_hot_sbs2_specific | *ANXA1* | Hot, immune-visible tumors |
| Tier 3 | 135 CNV-specific | 2_cold_cnv_specific | - | Cold, immune-evasive tumors |

Tier 1 is highest priority because its targets apply across the lifecycle axis; within it, the 22-gene escape subset (*KRT6A*, *PI3*, *KRT5*, *SERPINB2*) carries direct evidence of immune selection (in some CNV-HIGH cells the neoantigen-forming mutation has been removed by a correcting mutation or a fusion that deletes the junction, or the gene is downregulated), which de-risks immunogenicity. Tier 2 is immunostimulatory but unique to hot tumors. Tier 3 ranks last given the immunosuppressive state, but retains value because cold tumors otherwise lack targets. The three tiers sum to 516 genes; the Venn regions are 276 SBS2-specific / 105 shared / 135 CNV-specific (circle totals 381 SBS2, 240 CNV).

## Key Results

| Finding | Evidence |
|---------|----------|
| SBS2-HIGH carries more protein-altering variants | 1.13 vs 0.60 per cell |
| SBS2-HIGH produces more predicted neoantigens | 2,361 vs 1,331, 1.77-fold (Fig 7a) |
| SBS2-HIGH has more strong binders and differential neoantigens | 287 vs 164 (IC50 < 50 nM); 591 vs 352 |
| Neoantigen background is A3-skewed in SBS2-HIGH | TCW context 18.9% vs 9.0% (Fisher OR = 2.36, p = 0.0037); neoantigen-forming subset 17.8% vs 9.0% (p = 0.015) |
| Fusions do not track the immune divergence | CNV-HIGH 5,558 (10.2/cell), SBS2-HIGH 5,092 (9.3/cell), NORMAL 6,625 (12.1/cell) |
| Three tiers from 516 neoantigen genes | Tier 1 = 105 shared (*MDK*), Tier 2 = 276 SBS2-specific (*ANXA1*), Tier 3 = 135 CNV-specific |
| *ANXA1* is the top target | composite 324.9 (rank 1), 42 peptides, 8 strong binders, best IC50 = 13.1 nM, expressed in 99.5% of SBS2-HIGH cells (mean 8.48) |
| *ANXA1* mutations are not A3 products | seven somatic mutations across six positions, none in TCW context |

*ANXA1*'s mutations cluster in two hotspots within the annexin repeat domains: position 289 at the repeat 3/4 boundary (Phe289Leu, Phe289Ile) and positions 321 to 327 in repeat 4 (Asp321Glu, Ile322Val, Gln327Arg), with two more in repeat 3 (Asp258Glu, Met259Val). It ranks highly because it is recurrently mutated, near-uniformly expressed, and immune-visible, not because it is an A3 substrate.

## Figure 7 Panels

| Panel | Content | Manuscript |
|-------|---------|------------|
| A | Neoantigen and RNA fusion burden, SBS2-HIGH vs CNV-HIGH (1.77-fold neoantigen enrichment) | Fig 7a |
| B | Neoantigen-gene Venn (276 / 105 / 135) with tier breakout | Fig 7b |
| C | *ANXA1* top peptides by hotspot region, sorted by IC50 drop | Fig 7c |
| D | *ANXA1* gene track with domains and lollipop mutations | Fig 7d |

## Key Output Files

| File | Location | Description |
|------|----------|-------------|
| `pipeline_config.yaml` | `01_neoantigen_inputs/` | Paths, parameters, group sizes |
| `{group}.somatic_protein_altering.tsv` | `02_snpeff_annotation/` | Germline-subtracted missense variants |
| `{group}_neoantigens.tsv` | `03_mhc_binding/` | Predicted binders with peptide detail |
| `{group}_expression_weighted_ranking.tsv` | `03_mhc_binding/` | Composite-scored rankings |
| `per_group_junction_summary.tsv` | `04_fusion_analysis/` | Fusion rates per group |
| `vaccine_target_scores.tsv` | `05_summary/` | Integrated scoring with escape bonuses |
| `all_genes_tiered.tsv` | `05_summary/THERAPEUTIC_TIERS/` | Complete tier assignments |
| `Figure7_Combined.pdf` | `FIGURE_7_PANELS/` | Combined 4-panel figure |

## Reference Data

Ensembl GRCh38 release 115 proteome: `data/reference/Homo_sapiens.GRCh38.pep.all.fa` (245,535 isoforms, 19,885 gene symbols). The `pep.canonical.fa` file does not exist for release 115, so `pep.all.fa` is the correct source and selection is by longest isoform (the `Ensembl_canonical` tag is not present in release 115 headers).
