# Question 5 Workflow: Neoantigen Landscape (Figure 7)

**Maps to manuscript Results 4.5** — *Neoantigen landscape and therapeutic target identification.*

> **Script naming lags the manuscript.** The audit script is `diagnostic_section7_numbers.py` and its banner says "Section 4.5 / Figure 7"; the "section 7" in the filename dates from an earlier section order. Figure 7 is correct. This directory also uses `TROUBLE_SHOOTING` with an underscore, unlike `SINGLE_CELL/`, `NETWORK_SINGLE_CELL/` and `HPV_ANALYSIS/`.

> **CNA versus CNV.** The manuscript writes copy-number alteration as CNA; the code uses `CNV` throughout.

## Question 5: What therapeutic neoantigens fall out of the immune-visible and immune-evasive populations, are they clonally prevalent, and are they APOBEC-driven?

### Rationale

Question 4 (Fig 6) split the epithelial compartment into an immune-visible, maintenance-phase SBS2-HIGH state and an immune-evasive, productive-phase CNV-HIGH state. Because tumor-specific neoantigens are attractive targets for multivalent mRNA vaccines, and because how well that strategy works depends on the tumor's immune state, the neoantigens these two populations produce are relevant to different therapeutic settings. This question builds the neoantigen landscape from the single-cell somatic calls, predicts MHC-I binding, compares neoantigen and RNA-fusion burden between the two populations on a per-cell footing, ranks the candidate neoantigens by how clonally prevalent they are, and traces the selected candidates down to their genome-verified mutational context and their position in the encoded protein.

The central population finding is that the SBS2-HIGH population produces more predicted neoantigens than CNV-HIGH (1.82-fold at the variant level), an excess that survives normalization for sequencing depth, while RNA fusions do not track the immune divergence at all.

> **The TCW-enrichment claim is dead and must not be revived.** An earlier version of this document said the neoantigen excess "is APOBEC-shaped, biased toward the TCW context," resting on an 18.9% versus 9.0% Fisher comparison computed from SComatic's `REF_TRI` field. Genome-verified recomputation shows **no enrichment under any of three definitions**, all flat and non-significant (see Key Results). The excess is real; its TCW bias is not. The manuscript makes no TCW-enrichment claim, and the only APOBEC-context statement it retains is about the single *KRT6B* mutation.

An earlier expression-weighted composite put *ANXA1* on top because *ANXA1* is near-uniformly expressed (99.5% of SBS2-HIGH cells), but its neoantigen-forming mutation is carried in only 1.1% of them, too rare to be a strong vaccine target. The figure now selects candidates by clonal prevalence (carrier fraction among tumor cells), so it surfaces mutations actually present across the tumor. Three candidates are carried in the main figure, one per tier: *COX4I1* p.Ala9Thr (CNV-specific, 24.9% of CNV-HIGH cells), *SPRR1A* p.Val61Ile (shared, 16.9% of pooled tumor cells), and *KRT6B* p.Glu342Lys (SBS2-specific, 7.9% of SBS2-HIGH cells and the one clean-TCW hit of the three, with a large MHC-I binding gain).

### Note on units and numbering

- All tier and overlap counts are at the **neoantigen (mutation) level**, not the gene level, because a vaccine encodes specific peptides. The neoantigen-level overlap is 467 SBS2-specific / 93 shared / 215 CNV-specific (SBS2 total 560, CNV total 308, union 775). The earlier gene-level Venn (272 / 82 / 143) is retired; a gene can be gene-level "shared" while its specific mutations differ between groups, so the two units genuinely disagree.
- The pipeline was renumbered when the expression-weighted composite step was retired. The old Step04 (expression-weighted ranking) is deleted; the old Step05 (fusion) is now Step04, and the old Step06 (integrated analysis) is now Step05.

---

## Directory Structure

```
scripts/NEOANTIGEN/
├── RUN_NEOANTIGEN_PIEPLINE.sh              # SLURM wrapper (filename typo is as-shipped)
│
├── Step01_Prep_Neoantigen_Inputs.py        # per-group VCFs + barcode lists          [NETWORK]
├── Step02_SnpEff_Annotation.py             # SnpEff + germline subtraction           [NEOANTIGEN]
├── Step03_MHCflurry_Binding.py             # peptide build + MHC-I prediction        [NEOANTIGEN]
├── Step04a_STAR_Chimeric_Align.sh          # STAR chimeric alignment (array)
├── Step04b_Fusion_Analysis.py              # junction filtering + fusion burden      [NETWORK]
├── Step05_Integrated_Neoantigen_Analysis.py# tiers, selection evidence, aggregates   [NETWORK]
│
├── Phase5B_STAR_Chimeric_Pipeline.py       # predecessor of Step04a/b (superseded)
├── diagnostic_section7_numbers.py          # Results 4.5 number audit (assertion-based)
│
├── neoantigen_figure_utils.py              # shared loaders + renderers
├── Generate_Figure7_Panels.py              # main four-row figure
├── Generate_Figure7_Individual_Panels.py   # every panel as its own file
└── TROUBLE_SHOOTING/                       # verification stage; see below
```

`Phase5B_STAR_Chimeric_Pipeline.py` predates the Step04a/Step04b split and is retained for provenance. It is not called by the runner.

> **`TROUBLE_SHOOTING/` contains two nested `BACKUP/` trees** (`BEFORE_7-17-2026/` and `BEFORE_7-21-2026/`) each holding a full copy of the pipeline under the old step numbering, plus a third copy at the `BACKUP/` root. `Step03_MHCflurry_Binding.py` exists in four places. None but the main directory is current.

---

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
                    Step04a (shell) + Step04b   Verification stage
                    STAR chimeric align,        (TROUBLE_SHOOTING/) writes the
                    then fusion analysis        locked figure-input tables:
                               |                  - prevalence ranking (full.tsv)
                               |                  - Panel A burden stats
                               |                  - genome-verified TCW table
                               |                  - native-frame domain tracks
                               +-------------+-------------+
                                             |
                                             v
                                   Step05 (NETWORK)
                                   Integrated analysis (prevalence-first)
                                             |
                                             v
                          neoantigen_figure_utils.py
                               |                              |
                               v                              v
                    Generate_Figure7_Panels.py     Generate_Figure7_Individual_Panels.py

     No separate supplemental figure: full.tsv is the supplementary upload
     (manuscript Supp. File 5).
```

## Conda Environments

| Environment | Scripts | Key packages |
|-------------|---------|--------------|
| NETWORK | Steps 01, 04b, 05, the prevalence ranking, all figure generation | scanpy, anndata, pandas, matplotlib, gseapy |
| NEOANTIGEN | Steps 02, 03 | SnpEff (java), STAR 2.7.11b, mhcflurry, tensorflow, pysam |

`RUN_NEOANTIGEN_PIEPLINE.sh` switches environments automatically between steps.

## Pipeline Scripts (main directory)

### Step01_Prep_Neoantigen_Inputs.py
Loads the three-group assignments and the master SComatic filtered TSV, splits somatic variants into per-group VCFs (SBS2_HIGH, CNV_HIGH, NORMAL), and writes barcode lists and a config YAML. VCFs use 1-based coordinates (SComatic `Start` is 0-based, so `POS = Start + 1`). Outputs to `data/FIG_7/01_neoantigen_inputs/`.

### Step02_SnpEff_Annotation.py
Runs SnpEff (GRCh38.p14) on the per-group VCFs, parses the ANN field (gene, effect, HGVS protein notation, transcript), and performs germline subtraction using the NORMAL group as background: any variant present in NORMAL is dropped from the disease groups. Because SComatic calls variants within the basal-cell compartment only, a variant that reaches the disease groups is already somatic and cell-type restricted. Every candidate downstream is therefore somatic by construction; no separate germline filter is applied later. Outputs to `data/FIG_7/02_snpeff_annotation/`.

### Step03_MHCflurry_Binding.py
The core prediction step. For each somatic missense variant it builds mutant and wild-type peptides (lengths 8-11) from real protein context in the Ensembl GRCh38 release 115 proteome (`pep.all.fa`, 245,535 isoforms), replacing the earlier poly-alanine flanking that gave systematically wrong IC50 values. Protein lookup uses a six-layer chain (ENST, gene symbol, alias, ENSG, isoform scan, and a ±30 offset scan), reaching **98.6-98.7% mapping**. Peptides are scored against a 10-allele HLA panel with MHCflurry `Class1AffinityPredictor`; binders are IC50 < 500 nM, strong binders < 50 nM, and differential neoantigens have mutant < 500 nM with wild-type > 500 nM. Writes `{group}_neoantigens.tsv` (binder calls, with `mut_position_in_peptide`) and `{group}_all_peptide_results.tsv` (full audit) to `data/FIG_7/03_mhc_binding/`.

> **The HLA panel gives 90.5% world population coverage**, not the ~80% an earlier version of this document stated. The panel is A\*01:01, A\*02:01, A\*03:01, A\*24:02, B\*07:02, B\*08:01, B\*35:01, B\*44:02, C\*04:01, C\*07:01. The 90.5% comes from the IEDB Population Coverage tool, which is vendored at `TROUBLE_SHOOTING/population_coverage/` but is **not** re-run by the audit; it is a recorded value. Re-run it on the exact panel before any resubmission.

The ±30 offset scan exists because SnpEff protein positions differ from Ensembl GRCh38 r115 by up to ±30 residues where signal peptides are included.

### Step04a_STAR_Chimeric_Align.sh
Runs the STAR chimeric alignment across all cells to produce the per-cell chimeric-junction files that Step04b consumes. Kept as a shell step because it is the long, array-parallel alignment stage. Parsing uses `HPV_ANALYSIS/parse_chimeric_junctions.py`.

### Step04b_Fusion_Analysis.py
Re-parses the STAR chimeric junctions, filters to high-confidence fusions, and computes per-group fusion burden. Identifies group-exclusive fusion pairs and runs pathway enrichment (gseapy with backoff for Enrichr rate limiting), then cross-references fusion-disrupted genes against the neoantigen lists to find asymmetric escape. The population-level conclusion is that fusion burden is similar across groups, so fusions do not differentiate the two states; only specific partners do. The cross-group overlap it writes feeds the fusion-disruption mechanism in Step05. Outputs to `data/FIG_7/04_fusion_analysis/`.

### Step05_Integrated_Neoantigen_Analysis.py
Integrated analysis, built around the prevalence-weighted ranking. It reads neoantigen tier and per-niche prevalence from the single source (`neoantigen_prevalence_ranking_full.tsv`), the two `{group}_all_peptide_results.tsv` files, the fusion cross-reference, and adata for mean expression. It:

- Confirms the neoantigen-level overlap (467 / 93 / 215) independently and asserts it matches the ranking's tiers.
- Runs the selection-evidence analysis on the 93 shared neoantigens.
- Folds in the per-gene aggregates the retired Step04 used to provide.
- Writes a prevalence-ordered annotated target table with no composite score.

Removed relative to the old Step06: the composite `vaccine_score` and its multipliers, the hardcoded key-target list, the group-level neo:loss statistic, and dead imports. Outputs to `data/FIG_7/05_summary/`.

### diagnostic_section7_numbers.py
Assertion-based audit of every Results 4.5 number: PASS, **FAIL**, or SKIP per claim, with a tally. The September 2026 run returned **69 PASS, 0 FAIL**.

> **It depends on a TROUBLE_SHOOTING script.** CHECK A (group-rate counts, folds, binomial p) skips unless `panelA_grouprate_stats.tsv` exists, which is written by `Diagnostic_GroupAware_Expression_and_Carrier.py`. Run that first, then the audit, or the section's headline 1.82-fold and p = 9.50e-18 go unverified.

---

## Verification Stage (`TROUBLE_SHOOTING/`)

Read-only diagnostics that produce the locked figure-input tables.

### Diagnostic_Prevalence_Weighted_Neoantigen_Ranking.py
Ranks every neoantigen mutation from both groups by clonal prevalence (primary) and MHC-I binding gain (secondary), and writes the single figure-input table. It unions the two `{group}_neoantigens.tsv` binder sets and collapses each mutation `(gene, hgvs_p)` to its highest-binding-gain peptide (`delta_IC50 = wt_IC50 − mut_IC50`); assigns tier from `in_sbs2`/`in_cnv` in `ref_tri_fasta.tsv`; counts carriers in both groups; and reports per-niche prevalence, the tier-conditional headline (`prevalence_tier`), `prevalence_max`, and tier-conditional expression. Writes `neoantigen_prevalence_ranking_full.tsv` and a per-peptide `_long.tsv`.

### Diagnostic_GroupAware_Expression_and_Carrier.py
Produces `panelA_burden_stats.tsv` and `panelA_grouprate_stats.tsv`: the group-rate binomial tests, the per-cell per-UMI means with BH-adjusted p, and the tier-correct percent-expressing denominators. Source for Panel A and a prerequisite for the audit.

### TCW verification chain
`Diagnostic_Build_FASTA_Trinucleotide_Table.py`, `Diagnostic_FASTA_Arbiter.py`, `Diagnostic_TCW_Tri_vs_SComatic.py`, `Diagnostic_Differential_TCW_Audit.py` and `Diagnostic_FullSpectrum_TCW_Enrichment.py` together produce `ref_tri_fasta.tsv`: per-variant genome-verified reference base, trinucleotide context, TCW class (`is_tcw`, `is_tcw_ct`), tier membership, and nucleotide alt. TCW class is read from the GRCh38 genome, **not** from SComatic's `REF_TRI`, which produced a spurious enrichment signal.

### Diagnostic_Fetch_Protein_Domains.py
Fetches UniProt features for `FIGURE_GENES = ['COX4I1', 'SPRR1A', 'KRT6B']` and writes `protein_domains.tsv` in the native SnpEff-transcript frame. Must run on a node with outbound HTTPS to `rest.uniprot.org`.

### Others
`resolve_KRT6B_codon.py` and `resolve_KRT6B_HLA.py` (featured-candidate provenance), `Diagnostic_ANXA1_Provenance.py` and `Diagnostic_CAST_Provenance.py` (why the retired composite surfaced those genes), `Diagnostic_Proteome_Mapping.py`, `Diagnostic_HLA_panel.py`, the fusion audits, and `Collaborator_MHC_Peptide_Request.py` for the Ippolito validation request.

---

## Figure Generation

The two figure scripts import a shared engine and read binding, prevalence, and expression from one locked table; nothing is recomputed at plot time.

### neoantigen_figure_utils.py
Single source for binding, prevalence, and tier-conditional expression is `full.tsv`; the burden bars read `panelA_burden_stats.tsv`, the tracks read `protein_domains.tsv`, and the overlap Venn is computed at the neoantigen level (unique gene + hgvs_p) so it matches 467 / 93 / 215. The mutated residue in each displayed peptide is located from `mut_position_in_peptide` and drawn in red. Two color languages: the tier palette (CNV mustard `#F6D155`, shared purple `#9B59B6`, SBS2 coral `#ed6a5a`) for the prevalence and expression panels, and the TCW provenance palette (clean C>T coral `#ed6a5a`, C>G orange `#E67E22`, non-APOBEC gray `#9AA0A6`, disordered `#B0BEC5`) for the binding and track panels.

### Generate_Figure7_Panels.py (main figure)
Four-row layout. Row 1: A burden, B neoantigen overlap, C expression and carriage for the three featured genes. Row 2, D: MHC-I binding gain as horizontal stacked bars, ordered by prevalence, changed residue in red. Row 3, E: the *KRT6B* protein track full width. Row 4, F and G: the *COX4I1* and *SPRR1A* tracks paired. PDF and PNG at 300 DPI. Prerequisite: `protein_domains.tsv` must contain the three featured genes.

### Generate_Figure7_Individual_Panels.py
Every panel as its own standalone PDF and PNG for Illustrator, each with its own legend.

---

## Methodology

### Burden and comparison
Panel A is per cell and normalized per UMI so neither side is confounded by the shallower depth of SBS2-HIGH cells, and the two tests are corrected together by BH-FDR.

- Each burden is events per 1000 UMI per cell. The raw per-cell burden is depth-confounded and reported only as a cross-check (SBS2 4.203 vs CNV 3.799 carried variants per cell, MWU p = 0.331).
- Fusion burden is additionally corrected for within-patient germline junctions (removes 57 junctions, 0.53%, across SC003, SC005, SC006).
- The group-level headline uses a group-rate binomial test at equal 546-cell exposure.

### Prevalence and denominators
Carrier prevalence is counted from the single-cell genotype master on tier-consistent denominators: SBS2-specific over 546 SBS2-HIGH cells, CNV-specific over 546 CNV-HIGH cells, shared over the combined 1,092. `prevalence_max`, the stronger niche, is reported alongside so a shared hit is not diluted by a weak second niche. Because scRNA dropout means a true carrier is only detected where the locus is covered, prevalence is a floor; relative ordering is the trustworthy read. Expression uses the same tier denominators.

> **The SBS2-specific leader was decided on the tiebreak, not on prevalence.** *KRT6B* tied with *TACSTD2* at `prevalence_max` 0.0788 and won on binding gain (delta 27,713 versus 396). Methods 6.6 documents ranking "by clonal prevalence together with the MHC-I binding gain," so this is covered, but Results 4.5 describes prevalence alone and *KRT6B* is the candidate that needed the second criterion.

### Selection evidence (shared tier)
For the 93 shared neoantigens, Step05 asks whether each shows active evidence of being selected against in CNV-HIGH. A shared neoantigen carries selection evidence if it shows at least one of:

- **fusion disruption**: the gene appears in the cross-group fusion overlap as a CNV fusion that removes the neoantigen junction.
- **expression silencing**: the gene's mean expression in CNV is below half its SBS2 mean while still expressed in SBS2.

Two design choices, reported in Methods:

- An earlier antigen-loss signal (a gene-level flag on a *different* mutation than the neoantigen) was dropped as too weak and indirect.
- **HLA-A/B/C neoantigens are excluded from the count.** They remain in the binding analysis and the target table because their MHC-I gains are real, but those loci are hyperpolymorphic and mapping-artifact-prone, their apparent gains sit at germline-like prevalence, and HLA loss is its own escape story. Five HLA neoantigens had a mechanism and were excluded on this basis.

Per-niche prevalence is reported alongside as a descriptive breakdown but is not the evidence trigger, because depletion alone is confounded by A3A versus A3B generation rates.

---

## Key Results

All values verified September 2026 by `diagnostic_section7_numbers.py` (69 PASS, 0 FAIL) after running the group-aware diagnostic.

| Finding | Evidence |
|---------|----------|
| SBS2-HIGH produces more neoantigens (per-cell, depth-corrected) | 0.343 vs 0.184 per 1000 UMI, BH p = 2.15e-06 |
| Variant-level excess | 560 vs 308, **1.82-fold**, binomial p = 9.50e-18 |
| Peptide-level excess | 2,370 vs 1,339, 1.77-fold, p = 5.45e-65 (4.34 vs 2.45 peptides per cell) |
| Gene-level excess | 354 vs 225, 1.57-fold, p = 9.25e-08 |
| Strong binders (mut < 50 nM) | 264 SBS2 / 147 CNV |
| Differential binders (mut < 500, wt ≥ 500) | 615 SBS2 / 339 CNV |
| RNA fusions do not track the divergence | **0.547 vs 0.506** per 1000 UMI, BH p = 0.327 (ns) |
| Raw junction counts | SBS2 5,128 (9.4/cell), CNV 5,566 (10.2/cell), NORMAL 6,625 (12.1/cell) |
| Neoantigen-level overlap | 467 SBS2-specific / 93 shared / 215 CNV-specific (union 775) |
| Shared neoantigens with selection evidence | **15 of 93** |
| Selection-evidence breakdown | 12 depleted in CNV, 3 flat, 0 enriched; 10 by fusion, 6 by silencing; 5 HLA excluded |
| HLA panel coverage | 10 alleles, **90.5%** world population (IEDB) |
| Proteome mapping rate | 98.6-98.7% |

### TCW context: no enrichment, under any definition

Genome-verified from `ref_tri_fasta.tsv` across 672 SBS2-HIGH and 366 CNV-HIGH protein-altering variants:

| Definition | SBS2 | CNV | Fisher |
|---|---|---|---|
| A: among C>T protein-altering, fraction clean-TCW | 21/179 = 11.7% | 11/95 = 11.6% | OR 1.01, p = 1 |
| B: among all protein-altering, fraction TCW | 29/672 = 4.3% | 14/366 = 3.8% | OR 1.13, p = 0.75 |
| C: among neoantigen-forming variants, fraction TCW | 24/560 = 4.3% | 11/308 = 3.6% | OR 1.21, p = 0.72 |

All three flat and non-significant. This **resolves the pending reconciliation** that previously carried an 18.9% versus 9.0% enrichment from the retired `REF_TRI` field. The enrichment sentence stays removed.

### Featured candidates (Panel C), tier-conditional denominators

| Gene | Tier | Prevalence | % expressing | Binding gain | TCW |
|------|------|-----------|--------------|--------------|-----|
| *COX4I1* p.Ala9Thr | CNV-specific | 24.9% (/546 CNV) | 99.8% | 996 → 321 nM | no |
| *SPRR1A* p.Val61Ile | shared | 16.9% (/1092) | 46.2% | 444 → 295 nM | no |
| *KRT6B* p.Glu342Lys | SBS2-specific | 7.9% (/546 SBS2) | 66.5% | 27,793 → 80.5 nM | yes (C>T) |

For contrast, *ANXA1* is expressed in 99.5% of SBS2-HIGH cells but its neoantigen-forming mutation is carried in 1.1%. That gap is why the ranking moved from an expression-weighted composite to clonal prevalence.

---

## Exemplar Genes and Protein Tracks

Domains are drawn in the native SnpEff-transcript (HGVS) frame; where UniProt uses different numbering, the offset is solved and residue-verified before shifting.

| Gene | UniProt (native length) | Track content | Featured lollipop |
|------|-------------------------|---------------|-------------------|
| *KRT6B* | P04259 (564 aa) | keratin head, coils and rod with linkers, tail; disordered spans | p.Glu342Lys, clean TCW C>T (coral) |
| *COX4I1* | short, ~169 aa | short protein, backbone-dominant (**confirm boxes against `protein_domains.tsv`**) | p.Ala9Thr, non-APOBEC (gray) |
| *SPRR1A* | short, ~89 aa | low-complexity proline-rich, backbone-dominant (**confirm boxes**) | p.Val61Ile, non-APOBEC (gray) |

## Figure 7 Panels

No separate supplemental figure; `full.tsv` is manuscript Supplementary File 5.

| Panel | Content |
|-------|---------|
| A | per-cell mutational burden, neoantigen and RNA fusion, SBS2-HIGH vs CNV-HIGH |
| B | neoantigen overlap Venn (467 / 93 / 215, neoantigen level) |
| C | expression and neoantigen carriage for the three featured genes |
| D | MHC-I binding gain of the three featured mutations, changed residue in red |
| E | *KRT6B* protein track |
| F | *COX4I1* protein track |
| G | *SPRR1A* protein track |

## Therapeutic Tiers

Three priority tiers at the neoantigen level: Tier 1 shared (93), broadly conserved across both viral states; Tier 2 SBS2-specific (467), for hot, immune-visible tumors; Tier 3 CNV-specific (215), for cold, immune-evasive tumors.

Within Tier 1, the **15 of 93** shared neoantigens carrying selection evidence are the most compelling members, because active removal in the evasive state is independent evidence the immune system recognized them. The counted set decomposes as 12 depleted in CNV, 3 flat and 0 enriched, with 10 flagged by fusion disruption and 6 by expression silencing (one carries both). Five HLA-A/B/C neoantigens had a mechanism and were excluded from the count. Genes in the counted set include *SPRR1A*, *PI3*, *SERPINB2* and *KRT6A*.

> The manuscript describes the mechanism as "either disrupted by a CNA-HIGH RNA fusion ... or transcriptionally silenced." Because 10 + 6 = 16 against 15 mutations, one carries both, so "either / or" reads as more exclusive than the data. A minor wording point, noted rather than changed.

> **Retired:** the gene-level "22 genes" figure and the weak antigen-loss signal that supported it. Dropping that signal removes genes whose only evidence was antigen loss (for example *KRT5*).

## Key Output Files

| File | Location | Description |
|------|----------|-------------|
| `pipeline_config.yaml` | `01_neoantigen_inputs/` | Paths, parameters, group sizes |
| `{group}.somatic_protein_altering.tsv` | `02_snpeff_annotation/` | Germline-subtracted missense variants |
| `{group}_neoantigens.tsv` | `03_mhc_binding/` | Predicted binders with `mut_position_in_peptide` |
| `{group}_all_peptide_results.tsv` | `03_mhc_binding/` | Full per-peptide binding audit |
| `per_group_junction_summary.tsv` | `04_fusion_analysis/` | Fusion burden per group |
| `shared_neoantigen_selection_evidence.tsv` | `05_summary/` | Per shared neoantigen: mechanism flags, HLA flag, direction |
| `neoantigen_targets_annotated.tsv` | `05_summary/` | Prevalence-ordered targets with escape annotations |
| `per_gene_neoantigen_aggregates.tsv` | `05_summary/` | Per-gene neoantigen and strong-binder counts |
| `section45_number_audit.txt` | `05_summary/` | Audit report from `diagnostic_section7_numbers.py` |
| `neoantigen_prevalence_ranking_full.tsv` | `06_prevalence_ranking/` | Single figure source and Supp. File 5 |
| `neoantigen_prevalence_ranking_long.tsv` | `06_prevalence_ranking/` | Per-peptide/allele companion |
| `panelA_burden_stats.tsv`, `panelA_grouprate_stats.tsv` | `TROUBLESHOOTING/group_aware_expression/` | Panel A source; audit prerequisite |
| `ref_tri_fasta.tsv` | TROUBLE_SHOOTING outputs | Genome-verified TCW and tier membership |
| `Figure7_main.pdf` | `figures/` | Combined main figure |

## Reference Data

Ensembl GRCh38 release 115 proteome: `data/reference/Homo_sapiens.GRCh38.pep.all.fa` (245,535 isoforms). `pep.canonical.fa` does not exist for release 115, so `pep.all.fa` is the source and selection is by longest isoform. Trinucleotide context for TCW is read from the GRCh38 genome FASTA.

