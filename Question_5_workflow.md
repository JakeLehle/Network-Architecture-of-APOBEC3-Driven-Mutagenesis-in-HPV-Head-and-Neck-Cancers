# Question 5 Workflow: Neoantigen Landscape (Figure 7)

## Question 5: What therapeutic neoantigens fall out of the immune-visible and immune-evasive populations, are they clonally prevalent, and are they APOBEC-driven?

### Rationale

Question 4 (Fig 6) split the epithelial compartment into an immune-visible, maintenance-phase SBS2-HIGH state and an immune-evasive, productive-phase CNV-HIGH state. Because tumor-specific neoantigens are attractive targets for multivalent mRNA vaccines, and because how well that strategy works depends on the tumor's immune state, the neoantigens these two populations produce are relevant to different therapeutic settings. This question builds the neoantigen landscape from the single-cell somatic calls, predicts MHC-I binding, compares neoantigen and RNA-fusion burden between the two populations on a per-cell footing, ranks the candidate neoantigens by how clonally prevalent they are, and traces the selected candidates down to their genome-verified mutational context and their position in the encoded protein.

The central population findings are unchanged: the SBS2-HIGH population produces more predicted neoantigens than CNV-HIGH (1.77-fold at the group level), and that excess is APOBEC-shaped, biased toward the TCW context, while RNA fusions do not track the immune divergence at all.

What changed is how the figure selects the genes it carries. An earlier expression-weighted composite put *ANXA1* on top because *ANXA1* is near-uniformly expressed, but its neoantigen-forming mutation is carried in only about 1% of tumor cells, which is too rare to be a strong vaccine target. The figure now selects candidates by clonal prevalence (carrier fraction among tumor cells), so it surfaces mutations that are actually present across the tumor. Three candidates are carried in the main figure, one per tier: *COX4I1* p.Ala9Thr (CNV-specific, carried in 24.9% of CNV-HIGH cells), *SPRR1A* p.Val61Ile (shared, carried in 16.9% of tumor cells pooled), and *KRT6B* p.Glu342Lys (SBS2-specific, carried in 7.9% of SBS2-HIGH cells and the one clean-TCW APOBEC hit of the three, with a five-log MHC-I binding gain).

### Note on units and numbering

- All tier and overlap counts in this question are at the **neoantigen (mutation) level**, not the gene level, because a vaccine encodes specific peptides. The neoantigen-level overlap is 467 SBS2-specific / 93 shared / 215 CNV-specific (SBS2 total 560, CNV total 308). The earlier gene-level Venn (272 / 82 / 143) is retired; a gene can be gene-level "shared" while its specific mutations differ between groups, so the two units genuinely disagree.
- The pipeline was renumbered when the expression-weighted composite step was retired. The old Step04 (expression-weighted ranking) is deleted; the old Step05 (fusion) is now Step04, and the old Step06 (integrated analysis) is now Step05. Step numbers below use the new scheme.

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
                    Step04a (shell) + Step04b   Analysis / verification stage
                    STAR chimeric align,        (TROUBLE_SHOOTING/, excluded from
                    then fusion analysis        the step-by-step) writes the locked
                               |                figure-input tables:
                               |                  - prevalence-weighted ranking
                               |                    (single source: full.tsv)
                               |                  - group-aware Panel A burden stats
                               |                  - genome-verified TCW table
                               |                  - native-frame protein-domain tracks
                               |                           |
                               +-------------+-------------+
                                             |
                                             v
                                   Step05 (NETWORK)
                                   Integrated analysis (prevalence-first):
                                   neoantigen-level tiers, shared-tier
                                   selection-evidence, per-gene aggregates
                                             |
                                             v
                          neoantigen_figure_utils.py  (shared loaders + renderers)
                               |                              |
                               v                              v
                    Generate_Figure7_Panels.py     Generate_Figure7_Individual_Panels.py
                    (main, four-row layout)        (every panel as its own file)

     Supplemental figure retired: the full ranking table (full.tsv) is the
     supplementary upload. The expression-weighted composite step is deleted.
```

## Conda Environments

| Environment | Scripts | Key packages |
|-------------|---------|--------------|
| NETWORK | Steps 01, 04b, 05, the prevalence ranking, all figure generation | scanpy, anndata, pandas, matplotlib, gseapy |
| NEOANTIGEN | Steps 02, 03 | SnpEff (java), mhcflurry, tensorflow |

`RUN_NEOANTIGEN_PIEPLINE.sh` switches environments automatically between steps.

## Pipeline Scripts (main directory)

### Step01_Prep_Neoantigen_Inputs.py
Loads the three-group assignments and the master SComatic filtered TSV, splits somatic variants into per-group VCFs (SBS2_HIGH, CNV_HIGH, NORMAL), and writes barcode lists and a config YAML. VCFs use 1-based coordinates (SComatic `Start` is 0-based, so `POS = Start + 1`). Outputs to `data/FIG_7/01_neoantigen_inputs/`.

### Step02_SnpEff_Annotation.py
Runs SnpEff on the per-group VCFs, parses the ANN field (gene, effect, HGVS protein notation, transcript), and performs germline subtraction using the NORMAL group as background: any variant present in NORMAL is dropped from the disease groups. Because SComatic calls variants within the basal-cell compartment only, a variant that reaches the disease groups is already somatic and cell-type restricted, and this NORMAL subtraction removes anything also present in the low-burden basal cells. Every candidate downstream is therefore somatic by construction; no separate germline filter is applied later. Outputs to `data/FIG_7/02_snpeff_annotation/`.

### Step03_MHCflurry_Binding.py
The core prediction step. For each somatic missense variant it builds mutant and wild-type peptides (lengths 8-11) from real protein context in the Ensembl GRCh38 release 115 proteome (`pep.all.fa`, 245,535 isoforms), replacing the earlier poly-alanine flanking that gave systematically wrong IC50 values. Protein lookup uses a six-layer chain (ENST, gene symbol, alias, ENSG, isoform scan, and a +/-30 offset scan), reaching 98.6-98.7% mapping. Peptides are scored against a 10-allele HLA panel (~80% population coverage) with MHCflurry `Class1AffinityPredictor`; binders are IC50 < 500 nM, strong binders < 50 nM, and differential neoantigens have mutant < 500 nM with wild-type > 500 nM. Writes `{group}_neoantigens.tsv` (binder calls, with `mut_position_in_peptide`) and `{group}_all_peptide_results.tsv` (full audit) to `data/FIG_7/03_mhc_binding/`. These feed the prevalence ranking, the neoantigen-level overlap (Panel B), and the mutated-residue red highlight.

### Step04a_STAR_Chimeric_Align.sh
Runs the STAR chimeric alignment across all cells to produce the per-cell chimeric-junction files that Step04b consumes. Kept as a shell step because it is the long, array-parallel alignment stage. (Was Step05a.)

### Step04b_Fusion_Analysis.py
Re-parses the STAR chimeric junctions, filters to high-confidence fusions, and computes per-group fusion burden. Identifies group-exclusive fusion pairs and runs pathway enrichment (gseapy with backoff for Enrichr rate limiting), then cross-references fusion-disrupted genes against the neoantigen lists to find asymmetric escape (for example SBS2-exclusive fusions hitting *B2M* and *HLA-B*, CNV-exclusive fusions hitting *PSMB8* and *PSME2*). The population-level conclusion is that fusion burden is similar across groups, so fusions do not differentiate the two states; only specific partners do. The cross-group overlap it writes feeds the fusion-disruption mechanism in Step05. Outputs to `data/FIG_7/04_fusion_analysis/`. (Was Step05b.)

### Step05_Integrated_Neoantigen_Analysis.py
Integrated analysis, rebuilt around the prevalence-weighted ranking. It reads neoantigen tier and per-niche prevalence from the single source (`neoantigen_prevalence_ranking_full.tsv`), the two `{group}_all_peptide_results.tsv` files, the fusion cross-reference, and adata for mean expression. It:

- Confirms the neoantigen-level overlap (467 / 93 / 215) independently from the `{group}_neoantigens.tsv` files and asserts it matches the ranking's tiers.
- Runs the selection-evidence analysis on the 93 shared neoantigens (see Selection-Evidence Methodology).
- Folds in the per-gene aggregates the retired Step04 used to provide (neoantigen and strong-binder counts per group), computed from `{group}_all_peptide_results.tsv`.
- Writes a prevalence-ordered annotated target table with no composite score.

Removed relative to the old Step06: the composite `vaccine_score` and its breadth/pressure multipliers, the hardcoded key-target list, the group-level neo:loss statistic, and dead imports. Outputs to `data/FIG_7/05_summary/`. (Was Step06.)

### diagnostic_section7_numbers.py
Section 4.5 number-audit script: recomputes the burden, binder, fusion, and exemplar numbers from the pipeline outputs so each quantitative claim in the text is confirmable against one source.

### RUN_NEOANTIGEN_PIEPLINE.sh
SLURM wrapper that runs the steps sequentially with environment switching, fail-fast handling, and pre-flight input checks.

## Analysis / Verification Stage

Read-only diagnostics in `TROUBLE_SHOOTING/`, excluded from the step-by-step, that produce the locked figure-input tables.

### Diagnostic_Prevalence_Weighted_Neoantigen_Ranking.py
Ranks every neoantigen mutation from both groups by clonal prevalence (primary) and MHC-I binding gain (secondary), and writes the single figure-input table. It unions the two `{group}_neoantigens.tsv` binder sets and collapses each mutation `(gene, hgvs_p)` to its highest-binding-gain peptide (`delta_IC50 = wt_IC50 - mut_IC50`); assigns tier from `in_sbs2`/`in_cnv` in `ref_tri_fasta.tsv`; counts carriers in both groups; and reports per-niche prevalence (`prevalence_sbs2`, `prevalence_cnv`, each over 546), the tier-conditional headline (`prevalence_tier`: specific over 546, shared over 1092), `prevalence_max`, and tier-conditional expression. Writes `neoantigen_prevalence_ranking_full.tsv` (the single figure source and the manuscript supplementary table) and a per-peptide `_long.tsv`.

### Group-aware burden diagnostic
Produces `panelA_burden_stats.tsv` (neoantigen per-cell per-UMI mean, SEM, BH-adjusted p) and the fusion burden stats; source for Panel A.

### Genome-verification / TCW diagnostic
Produces `ref_tri_fasta.tsv`: per-variant genome-verified reference base, trinucleotide context, TCW class (`is_tcw`, `is_tcw_ct`), tier membership (`in_sbs2`, `in_cnv`), and nucleotide alt. TCW class is read from the GRCh38 genome, not from SComatic's retired `REF_TRI`.

### Diagnostic_Fetch_Protein_Domains.py
Fetches UniProt features for the featured genes (`FIGURE_GENES = ['COX4I1', 'SPRR1A', 'KRT6B']`) and writes `protein_domains.tsv` in the native SnpEff-transcript frame. Must run on a node with outbound HTTPS to `rest.uniprot.org`.

## Figure Generation

The two figure scripts import a shared engine and read binding, prevalence, and expression from one locked table (`full.tsv`); nothing is recomputed at plot time.

### neoantigen_figure_utils.py
Single source for binding, prevalence, and tier-conditional expression is `full.tsv`; the burden bars read `panelA_burden_stats.tsv`, the tracks read `protein_domains.tsv`, and the overlap Venn is computed at the neoantigen level (unique gene+hgvs_p) so it matches the ranking's 467 / 93 / 215. The mutated residue in each displayed peptide is located from `mut_position_in_peptide` and drawn in red. Two color languages: the tier palette (CNV mustard #F6D155, shared purple #9B59B6, SBS2 coral #ed6a5a) for the prevalence and expression panels, and the TCW provenance palette (clean C>T coral #ed6a5a, C>G orange #E67E22, non-APOBEC gray #9AA0A6, disordered #B0BEC5) for the binding and track panels. Gene and protein names are drawn plain.

Note: `load_venn` counts unique (gene, hgvs_p) neoantigen mutations, not gene symbols, so Panel B reports 467 / 93 / 215.

### Generate_Figure7_Panels.py (main figure)
Four-row layout. Row 1: A burden, B neoantigen overlap (Venn, neoantigen level), C expression and carriage for the three featured genes. Row 2, D: MHC-I binding gain for the three featured mutations as horizontal stacked bars, ordered by prevalence, changed residue in red. Row 3, E: the *KRT6B* protein track full width. Row 4, F and G: the *COX4I1* and *SPRR1A* tracks paired. PDF and PNG at 300 DPI. Prerequisite: `protein_domains.tsv` must contain the three featured genes.

### Generate_Figure7_Individual_Panels.py
Every panel as its own standalone PDF and PNG for Illustrator: burden, Venn, expression-and-carriage, the three-mutation binding panel, and a track file per featured gene. Each carries its own legend.

## Figure Input Tables

| Table | Provides | Read by |
|-------|----------|---------|
| `neoantigen_prevalence_ranking_full.tsv` | per-mutation binding, tier, per-niche and tier-conditional prevalence, tier-conditional % expressing | binding (D), expression (C); also the manuscript supplementary table |
| `{group}_neoantigens.tsv` | neoantigen sets and `mut_position_in_peptide` | Venn (B), mutated-residue highlight |
| `panelA_burden_stats.tsv` | neoantigen per-cell per-UMI mean, SEM, BH-adjusted p | burden (A) |
| `ref_tri_fasta.tsv` | genome-verified TCW class, tier membership, nucleotide alt | tier and TCW color; consumed by the ranking |
| `protein_domains.tsv` + `{gene}_mutation_map.tsv` | native-frame domain boxes, disordered spans, lollipop positions | tracks (E, F, G) |

## Burden and Comparison Methodology

The neoantigen-versus-fusion comparison in Panel A is per cell and normalized per UMI so neither side is confounded by the shallower depth of SBS2-HIGH cells, and the two tests are corrected together by BH-FDR.

- Each burden is events per 1000 UMI per cell. The raw per-cell burden is depth-confounded and reported only as a cross-check.
- Fusion burden is additionally corrected for within-patient germline junctions (removes 57 junctions, 0.53%, across SC003, SC005, SC006; SBS2 mean 0.573 to 0.567).
- The group-level headline (equal 546-cell exposure) is a group-rate binomial test.

## Prevalence and Denominator Methodology

Carrier prevalence is counted from the single-cell genotype master and reported on tier-consistent denominators: SBS2-specific over 546 SBS2-HIGH cells, CNV-specific over 546 CNV-HIGH cells, shared over the combined 1092 (`prevalence_tier`). `prevalence_max`, the stronger niche, is reported alongside so a shared hit is not diluted by a weak second niche. Because scRNA dropout means a true carrier is only detected where the locus is covered, prevalence is a floor; relative ordering is the trustworthy read. Expression uses the same tier denominators so the expression bar and the carriage segment share one population.

## Selection-Evidence Methodology (shared tier)

For the 93 shared neoantigens, Step05 asks whether each shows active evidence of being selected against in CNV-HIGH. A shared neoantigen carries selection evidence if it shows at least one of two direct removal mechanisms:

- **fusion disruption**: the gene appears in the cross-group fusion overlap as a CNV fusion that removes the neoantigen junction.
- **expression silencing**: the gene's mean expression in CNV is below half its SBS2 mean (>2-fold down) while still expressed in SBS2.

Two design choices matter for the count and are reported transparently in Methods:

- An earlier antigen-loss signal (a gene-level, wild-type-binds / mutation-destroys flag on a *different* mutation than the neoantigen) was dropped as too weak and indirect; it inflated the count and pulled in CNV-enriched, direction-wrong hits.
- **HLA-A/B/C neoantigens are excluded from the count.** They remain in the binding analysis and the target table because their MHC-I gains are real, but those loci are hyperpolymorphic and mapping-artifact-prone, their apparent gains sit at germline-like prevalence, and HLA loss is its own escape story. They are flagged, not silently dropped, so the exclusion is auditable.

Per-niche prevalence (`prevalence_sbs2` vs `prevalence_cnv`) is reported alongside as the descriptive "minority of CNV cells still carry it" observation, with a depleted / flat / enriched breakdown of the counted set, but it is not the evidence trigger (depletion alone is confounded by A3A vs A3B generation rates).

## Key Results

| Finding | Evidence |
|---------|----------|
| SBS2-HIGH produces more predicted neoantigens (per-cell, depth-corrected) | 0.343 vs 0.184 per 1000 UMI, BH-adjusted p = 2.15e-06 (Panel A) |
| Group-level neoantigen excess | 4.34 vs 2.45 per cell; 2,370 vs 1,339 peptides, 1.77-fold, binomial p = 5.45e-65 |
| Holds at variant and gene level | neoantigen variants 560 vs 308 (1.82-fold, p = 9.50e-18); neoantigen genes 354 vs 225 (1.57-fold, p = 9.25e-08) |
| RNA fusions do not track the divergence | 0.567 vs 0.526 per 1000 UMI, BH-adjusted p = 0.327 (ns) |
| Neoantigen-level overlap (Panel B) | 467 SBS2-specific / 93 shared / 215 CNV-specific |
| Genome-verified | all 923 protein-altering variants match GRCh38 at their locus; TCW from the genome, not `REF_TRI` |
| Featured candidates span the tiers by prevalence | *COX4I1* p.Ala9Thr 24.9% (CNV), *SPRR1A* p.Val61Ile 16.9% (shared), *KRT6B* p.Glu342Lys 7.9% (SBS2) |
| The SBS2-specific prevalent candidate is APOBEC-driven | *KRT6B* p.Glu342Lys is clean TCW C>T, wt 27793 to mut 80.5 nM |

Featured candidates (Panel C), tier-conditional denominators:

| Gene | Tier | Prevalence (= carriage) | % expressing | Binding gain | TCW |
|------|------|-------------------------|--------------|--------------|-----|
| *COX4I1* p.Ala9Thr | CNV-specific | 24.9% (/546 CNV) | 99.8% | 996 to 321 nM | no |
| *SPRR1A* p.Val61Ile | shared | 16.9% (/1092) | 46.2% | 444 to 295 nM | no |
| *KRT6B* p.Glu342Lys | SBS2-specific | 7.9% (/546 SBS2) | 66.5% | 27793 to 80.5 nM | yes (C>T) |

## Exemplar Genes and Protein Tracks

Domains are drawn in the native SnpEff-transcript (HGVS) frame; where UniProt uses different numbering, the offset is solved and residue-verified before shifting.

| Gene | UniProt (native length) | Track content | Featured lollipop |
|------|-------------------------|---------------|-------------------|
| *KRT6B* | P04259 (564 aa) | keratin head, coils and rod with linkers, tail; disordered spans | p.Glu342Lys, clean TCW C>T (coral) |
| *COX4I1* | (short, ~169 aa) | short protein, backbone-dominant (exact boxes to confirm against `protein_domains.tsv`) | p.Ala9Thr, non-APOBEC (gray) |
| *SPRR1A* | (short, ~89 aa) | low-complexity proline-rich, backbone-dominant (exact boxes to confirm against `protein_domains.tsv`) | p.Val61Ile, non-APOBEC (gray) |

## Figure 7 Panels

Main figure (no separate supplemental figure; `full.tsv` is the supplementary upload):

| Panel | Content |
|-------|---------|
| A | per-cell mutational burden, neoantigen and RNA fusion, SBS2-HIGH vs CNV-HIGH |
| B | neoantigen overlap Venn (467 / 93 / 215, neoantigen level) |
| C | expression and neoantigen carriage for the three featured genes (bar = % expressing on the tier denominator; colored base = prevalence, colored by tier) |
| D | MHC-I binding gain of the three featured mutations (wild-type baseline + gain, changed residue in red) |
| E | *KRT6B* protein track |
| F | *COX4I1* protein track |
| G | *SPRR1A* protein track |

## Therapeutic Tiers

Neoantigens are grouped into three priority tiers at the neoantigen level: Tier 1 shared (93 neoantigens), broadly conserved across both viral states; Tier 2 SBS2-specific (467), for hot, immune-visible tumors; Tier 3 CNV-specific (215), for cold, immune-evasive tumors. Within Tier 1, the shared neoantigens that also carry selection evidence (fusion disruption or expression silencing in CNV, HLA-A/B/C excluded from the count) are the most compelling members, because active removal in the evasive state is independent evidence the immune system recognized them. The exact count of shared neoantigens with selection evidence comes from the Step05 run and replaces the earlier gene-level "22 genes" figure. Note that dropping the weak antigen-loss signal removes genes whose only evidence was antigen loss (for example *KRT5*), while genes with fusion or silencing evidence (for example *KRT6A*, *PI3*, *SERPINB2*, and the featured *SPRR1A*) remain.

## Key Output Files

| File | Location | Description |
|------|----------|-------------|
| `pipeline_config.yaml` | `01_neoantigen_inputs/` | Paths, parameters, group sizes |
| `{group}.somatic_protein_altering.tsv` | `02_snpeff_annotation/` | Germline-subtracted missense variants |
| `{group}_neoantigens.tsv` | `03_mhc_binding/` | Predicted binders with peptide detail and `mut_position_in_peptide` |
| `{group}_all_peptide_results.tsv` | `03_mhc_binding/` | Full per-peptide binding audit |
| `neoantigen_prevalence_ranking_full.tsv` | `06_prevalence_ranking/` | Single figure source and supplementary table (per mutation) |
| `neoantigen_prevalence_ranking_long.tsv` | `06_prevalence_ranking/` | Per-peptide/allele companion |
| `per_group_junction_summary.tsv` | `04_fusion_analysis/` | Fusion burden per group |
| `shared_neoantigen_selection_evidence.tsv` | `05_summary/` | Per shared neoantigen: mechanism flags, HLA flag, direction |
| `neoantigen_targets_annotated.tsv` | `05_summary/` | Prevalence-ordered neoantigen targets with escape annotations |
| `per_gene_neoantigen_aggregates.tsv` | `05_summary/` | Per-gene neoantigen and strong-binder counts (folded from old Step04) |
| `Figure7_main.pdf` | `figures/` | Combined main figure |
| `Figure7_{panel}.pdf` / `Figure7_{GENE}_track.pdf` | `figures/` | Standalone panels |

## Reference Data

Ensembl GRCh38 release 115 proteome: `data/reference/Homo_sapiens.GRCh38.pep.all.fa` (245,535 isoforms). `pep.canonical.fa` does not exist for release 115, so `pep.all.fa` is the source and selection is by longest isoform. Trinucleotide context for TCW is read from the GRCh38 genome FASTA, not from `REF_TRI`.

## Pending Reconciliation

1. The group-level TCW-context enrichment claim (previously 18.9% vs 9.0%, Fisher exact) was computed off the retired `REF_TRI` and should be recomputed from `ref_tri_fasta.tsv`.
2. Strong-binder and differential counts quoted in earlier drafts should be refreshed from the current `{group}_all_peptide_results.tsv` (now aggregated by Step05).
3. The exact protein-track content for *COX4I1* and *SPRR1A* should be read back from the domain-fetcher output and written into the Exemplar Genes table.
4. The shared-tier selection-evidence count and its escape-subset gene list should be taken from the Step05 run and written into the Tier 1 results paragraph.
