
======================================================================
SECTION 2 TEXT NUMBERS — FIGURE 2
======================================================================

--- FIG_2 Top-Level Contents ---
  [DIR] 01_cleaned_expression/
  [DIR] 02_merged_with_SBS/
  [DIR] 03_differential_expression/
  [DIR] 04_correlation_networks/
  [DIR] 05_communities/
  [DIR] 06_centrality_metrics/
  [DIR] 08_pipeline_summary/
  [DIR] DIFF_threshold_sweep/
  [DIR] FIGURE_2_PANELS/

--- Network Group Sizes ---
  File: /master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_1/HNSC_network_groups_v3.tsv
  Columns: ['Entity_ID', 'Group', 'SBS2', 'A3A_plus_A3B']
  Groups (Group):
    HIGH: n=54
    LOW: n=57

--- Differential Expression ---
  DE directory: /master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_2/03_differential_expression/TCGA-HNSC
  Found: TCGA-HNSC_diffexpr_stats.csv
    Shape: (19101, 1)
    Columns: ['gene,gene_symbol,U_stat,p_value,logFC,mean_log1p_HIGH,mean_log1p_LOW,mean_raw_HIGH,mean_raw_LOW,p_adj']

--- Network Structure ---

  Community gene lists: /master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_2/05_communities/TCGA-HNSC/TCGA-HNSC_community_gene_lists.csv
  Columns: ['community', 'size', 'is_satellite', 'genes']
  Total communities: 44
  Total genes: 44
    Community 0: 1 genes
    Community 1: 1 genes
    Community 2: 1 genes
    Community 3: 1 genes
    Community 4: 1 genes
    Community 5: 1 genes
    Community 6: 1 genes
    Community 7: 1 genes
    Community 8: 1 genes
    Community 9: 1 genes
    Community 10: 1 genes
    Community 12: 1 genes
    Community 13: 1 genes
    Community 14: 1 genes
    Community 15: 1 genes
    Community 16: 1 genes
    Community 17: 1 genes
    Community 18: 1 genes
    Community 19: 1 genes
    Community 20: 1 genes
    Community 21: 1 genes
    Community 22: 1 genes
    Community 23: 1 genes
    Community 24: 1 genes
    Community 25: 1 genes
    Community 26: 1 genes
    Community 27: 1 genes
    Community 28: 1 genes
    Community 29: 1 genes
    Community 30: 1 genes
    Community 31: 1 genes
    Community 32: 1 genes
    Community 33: 1 genes
    Community 34: 1 genes
    Community 35: 1 genes
    Community 36: 1 genes
    Community 37: 1 genes
    Community 38: 1 genes
    Community 39: 1 genes
    Community 40: 1 genes
    Community 41: 1 genes
    Community 42: 1 genes
    Community 43: 1 genes
    Community 44: 1 genes

  Best partition: /master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_2/05_communities/TCGA-HNSC/TCGA-HNSC_best_partition.csv
  Shape: (1131, 3)
  Columns: ['gene', 'gene_symbol', 'community']

  Found summary file: /master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_2/05_communities/TCGA-HNSC/TCGA-HNSC_community_summary.txt
Community Summary: TCGA-HNSC
DIFF threshold: 0.65
Leiden resolution: 0.5
Clustering scope: full network (all 34 connected components)
Total genes: 1131

Community 0 (n=281)
  ABCC1, AC013394.1, ACOD1, ACTA1, ACTG1, ACTN2, ADCY2, ADRA2C, AFAP1L2, ANGPTL1
  ANKHD1, ANKRD2, AP2S1, AQP5, AQP9, ASB15, ASB5, ASPN, ATP1B4, AXDND1
  BAHD1, BRD3OS, BRINP1, BVES, C11orf68, C1QTNF3, C1orf105, CAMK1D, CAP2, CASQ2
  CBWD3, CCDC197, CCL21, CD14, CD163, CDC25A, CDCP1, CDRT15, CHEK1, CHST14
  CIART, CKM, CKMT2, CLEC2B, CLEC4E, CLEC7A, CLRN3, COL25A1, COX6A2, COX7A1
  CREG2, CRISP2, CRYAB, CSDC2, CSRP3, CTSL, CXCL9, CXXC1, CYP2W1, DDN
  DGKG, DIRAS1, DLL4, DMXL2, DNAJB4, DNAJC12, DNALI1, DTNA, ECEL1, EEF1A2
  EIF1AD, EIF2S3B, ELAPOR1, ENO3, ERBB4, ESAM, ETNK2, EYA4, FAM117A, FAM240C
  FAM78B, FANCA, FBXO27, FCGR1A, FCGR1B, FCGR2A, FCGR3A, FGF6, FIBCD1, FNDC5
  FREM2, FZD1, GALNT18, GBP1, GBP3, GBP5, GCH1, GDA, GDF3, GDF5
  GGACT, GPD1, GRK7, GTF3C3, HACD1, HJV, HK3, HMCN1, HNRNPM, HRC
  HSBP1L1, HSPB3, HSPB6, HSPB7, IDO1, IGF2, IGSF23, IL12RB2, IL17D, IL1RAPL1
  IMPG1, INPP5A, IP6K3, IQCM, IQGAP3, ITGAX, ITGB1BP2, ITIH3, KCMF1, KLHL13
  KLHL41, KLK4, KLRC3, LAYN, LDB3, LHFPL6, LHX8, LILRA6, LILRB3, LMOD2
  LRRC17, LRRC30, LRRN4CL, LTA4H, MAFA, MARCOL, MCHR2, MEFV, MEIS2, MGAM
  MLKL, MLLT11, MSR1, MTHFR, MYH1, MYH2, MYH3, MYL1, MYLPF, MYMK
  MYO18B, MYO7A, MYOG, MYOM2, MYOZ1, NACAD, NCKIPSD, NECAP1, NKX2-3, NKX3-1
  NLRC5, NLRP3, NOTUM, NRAP, NRK, ODF3, OMD, OR2B11, OR2K2, OR52L1
  P3H4, PC, PDCD1LG2, PGPEP1L, PHYHD1, PHYHIP, PITX3, PLCG2, PLEKHG6, PLP1
  PNMT, PNPLA7, POLE, POLM, POPDC3, PPFIBP1, PPP1R1B, PPP1R1C, PPP1R27, PPP1R3A
  PRKG1, PRKN, PRPH2, PRR27, PRSS2, PSEN2, PTBP1, PTGER3, PYGB, RABGAP1
  RASGEF1A, RASL10B, RBM19, RBM20, RBM24, RERGL, RETREG1, RGS17, RIC1, RNF208
  SCGB1D2, SCN4B, SCN5A, SDC3, SERPINA5, SERPING1, SGCA, SGPP1, SH2D3A, SHANK2
  SHISA4, SHOX2, SLC13A3, SLC13A5, SLC38A11, SLC38A4, SLTM, SMPX, SNX18, SORBS1
  SORCS2, SPARCL1, SQSTM1, SRPK3, SSTR5, STAR, S
  ... (truncated)

  Found summary file: /master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_2/05_communities/TCGA-HNSC/TCGA-HNSC_selected_parameters.txt
DIFF_THRESHOLD=0.65
THRESHOLD_METHOD=max_fragmentation_rate
CLUSTERING_SCOPE=full_network
N_COMPONENTS=34
LEIDEN_RESOLUTION=0.5
RESOLUTION_METHOD=composite_score_mod_x_ari_x_evenness
MODULARITY=0.4939
ARI=0.5221
EVENNESS=0.7515
COMPOSITE_SCORE=0.1938
A3_SEEDS=ENSG00000128383,ENSG00000179750


  Found summary file: /master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_2/08_pipeline_summary

  Found summary file: /master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_2/08_pipeline_summary/TCGA-HNSC_KEGG_enrichment_summary.csv
community,n_genes,n_mapped,top_kegg_term,top_kegg_pvalue,top_kegg_genes,n_kegg_significant
0,281,281,Cardiac muscle contraction,0.0485398637691769,TNNC1;TPM1;CASQ2;HRC;ATP1B4;COX6A2;COX7A1,1
1,257,257,Dilated cardiomyopathy,0.0469694905082206,CACNB1;SGCD;LAMA2;ITGAV;CACNG1;AGT;ADCY5,2
2,134,134,Coronavirus disease,0.7871730963168916,ACE2;RPS6;RPL14;RPSA;RPL29,0
3,130,130,Neutrophil extracellular trap formation,0.6059766665785621,HDAC5;HDAC3;CASP1;MAPK12;HDAC7,0
4,121,121,Gastric cancer,0.2505802971392994,CDKN2B;WNT5B;WNT3A;EGF;WNT16,0
5,44,44,Glycosaminoglycan biosynthesis,0.2428067423834619,HS3ST3B1;HS3ST3A1,0
6,33,33,Herpes simplex virus 1 infection,3.363598749136468e-05,ZNF571;ZNF175;ZNF461;ZNF383;ZFP30;ZNF181;ZNF600;ZFP28,1
7,23,23,Basal transcription factors,0.0182290815956702,ERCC3;TAF8,1
8,16,16,Epithelial cell signaling in Helicobacter pylori infection,0.0009173074065529596,CDC42;CXCL3;CXCL2,14
9,11,11,Synthesis and degradation of ketone bodies,0.0391376120777327,HMGCS2,6
10,6,6,Staphylococcus aureus infection,0.0588931127910019,KRT35,0
12,3,3,Systemic lupus erythematosus,0.0337009782356269,H2AC12,6
13,3,3,Endocrine and other factor-regulated calcium reabsorption,0.0253064523102616,DNM2,8
14,3,3,Pentose phosphate pathway,0.0197680671855171,ALDOA,8
15,3,3,No results,,,0
16,3,3,Olfactory transduction,0.0014275885134425,OR52W1;OR5M1,1
17,3,3,No results,,,0
18,3,3,Oxytocin signaling pathway,0.0372637172825571,RGS2,2
19,2,2,N/A (too few genes),,,0
20,2,2,N/A (too few genes),,,0
21,2,2,N/A (too few genes),,,0
22,2,2,N/A (too few genes),,,0
23,2,2,N/A (too few genes),,,0
24,2,2,N/A (too few genes),,,0
25,2,2,N/A (too few genes),,,0
26,2,2,N/A (too few genes),,,0
27,2,2,N/A (too few genes),,,0
28,2,2,N/A (too few genes),,,0
29,2,2,N/A (too few genes),,,0
30,2,2,N/A (too few genes),,,0
31,2,2,N/A (too few genes),,,0
32,2,2,N/A (too few genes),,,0
33,2,2,N/A (too few genes),,,0
34,2,2,N/A (too few genes),,,0
35,2,2,N/A (too few genes),,,0
36,2,2,N/A (too few genes),,,0
37,2,2,N/A (too few genes),,,0
38,2,2,N/A (too few genes),,,0
39,2,2,N/A (too few genes),,,0
40,2,2,N/A (too few genes),,,0
41,2,2,N/A (too few genes),,,0
42,2,2,N/A (too few genes),,,0
43,2,2,N/A (too few genes),,,0
44,4,4,Ribosome,0.0456020005475511,RPS4Y1,4


--- A3 Gene Centrality ---
  File: TCGA-HNSC_BOT_metrics.csv
           node  degree  betweenness  closeness  eigenvector  strength_abs gene_symbol
ENSG00000179750       0          0.0        0.0 1.343303e-10           0.0    APOBEC3B
ENSG00000128383       0          0.0        0.0 1.343303e-10           0.0    APOBEC3A
  File: TCGA-HNSC_DIFF_metrics.csv
           node  degree  betweenness  closeness  eigenvector  strength_abs gene_symbol
ENSG00000179750       1          0.0   0.058407 4.961808e-05      0.708126    APOBEC3B
ENSG00000128383       0          0.0   0.000000 3.903327e-33      0.000000    APOBEC3A
  File: TCGA-HNSC_DIFF_metrics_with_communities.csv
           node  degree  betweenness  closeness  eigenvector  strength_abs gene_symbol  community
ENSG00000179750       1          0.0   0.058407 4.961808e-05      0.708126    APOBEC3B        2.0
ENSG00000128383       0          0.0   0.000000 3.903327e-33      0.000000    APOBEC3A       43.0
  File: TCGA-HNSC_TOP_metrics.csv
           node  degree  betweenness  closeness  eigenvector  strength_abs gene_symbol
ENSG00000128383       0          0.0        0.0 1.128084e-24           0.0    APOBEC3A
ENSG00000179750       0          0.0        0.0 1.128084e-24           0.0    APOBEC3B

--- Diagnostic Audit Files ---

--- KEGG Enrichment ---

  TCGA-HNSC_KEGG_enrichment_summary.csv:
  Shape: (44, 7)
  Columns: ['community', 'n_genes', 'n_mapped', 'top_kegg_term', 'top_kegg_pvalue', 'top_kegg_genes', 'n_kegg_significant']

  TCGA-HNSC_community_00_KEGG.csv:
  Shape: (216, 10)
  Columns: ['Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Old P-value', 'Old Adjusted P-value', 'Odds Ratio', 'Combined Score', 'Genes']
  Significant (Adjusted P-value < 0.05): 1
       Gene_set                       Term Overlap  P-value  Adjusted P-value  Old P-value  Old Adjusted P-value  Odds Ratio  Combined Score                                     Genes
KEGG_2021_Human Cardiac muscle contraction    7/87 0.000225           0.04854            0                     0    6.271578       52.685325 TNNC1;TPM1;CASQ2;HRC;ATP1B4;COX6A2;COX7A1

  TCGA-HNSC_community_01_KEGG.csv:
  Shape: (200, 10)
  Columns: ['Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Old P-value', 'Old Adjusted P-value', 'Odds Ratio', 'Combined Score', 'Genes']
  Significant (Adjusted P-value < 0.05): 2
       Gene_set                                            Term Overlap  P-value  Adjusted P-value  Old P-value  Old Adjusted P-value  Odds Ratio  Combined Score                                    Genes
KEGG_2021_Human                          Dilated cardiomyopathy    7/96  0.00024          0.046969            0                     0    6.183281       51.530438 CACNB1;SGCD;LAMA2;ITGAV;CACNG1;AGT;ADCY5
KEGG_2021_Human Arrhythmogenic right ventricular cardiomyopathy    6/77  0.00047          0.046969            0                     0    6.623197       50.756391      CACNB1;SGCD;LAMA2;TCF7;ITGAV;CACNG1

  TCGA-HNSC_community_02_KEGG.csv:
  Shape: (137, 10)
  Columns: ['Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Old P-value', 'Old Adjusted P-value', 'Odds Ratio', 'Combined Score', 'Genes']
  Significant (Adjusted P-value < 0.05): 0

  TCGA-HNSC_community_03_KEGG.csv:
  Shape: (146, 10)
  Columns: ['Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Old P-value', 'Old Adjusted P-value', 'Odds Ratio', 'Combined Score', 'Genes']
  Significant (Adjusted P-value < 0.05): 0

  TCGA-HNSC_community_04_KEGG.csv:
  Shape: (146, 10)
  Columns: ['Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Old P-value', 'Old Adjusted P-value', 'Odds Ratio', 'Combined Score', 'Genes']
  Significant (Adjusted P-value < 0.05): 0

  TCGA-HNSC_community_05_KEGG.csv:
  Shape: (40, 10)
  Columns: ['Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Old P-value', 'Old Adjusted P-value', 'Odds Ratio', 'Combined Score', 'Genes']
  Significant (Adjusted P-value < 0.05): 0

  TCGA-HNSC_community_06_KEGG.csv:
  Shape: (30, 10)
  Columns: ['Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Old P-value', 'Old Adjusted P-value', 'Odds Ratio', 'Combined Score', 'Genes']
  Significant (Adjusted P-value < 0.05): 1
       Gene_set                             Term Overlap  P-value  Adjusted P-value  Old P-value  Old Adjusted P-value  Odds Ratio  Combined Score                                                 Genes
KEGG_2021_Human Herpes simplex virus 1 infection   8/498 0.000001          0.000034            0                     0   12.719673      174.273663 ZNF571;ZNF175;ZNF461;ZNF383;ZFP30;ZNF181;ZNF600;ZFP28

  TCGA-HNSC_community_07_KEGG.csv:
  Shape: (15, 10)
  Columns: ['Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Old P-value', 'Old Adjusted P-value', 'Odds Ratio', 'Combined Score', 'Genes']
  Significant (Adjusted P-value < 0.05): 1
       Gene_set                        Term Overlap  P-value  Adjusted P-value  Old P-value  Old Adjusted P-value  Odds Ratio  Combined Score      Genes
KEGG_2021_Human Basal transcription factors    2/45 0.001215          0.018229            0                     0   44.150609      296.373647 ERCC3;TAF8

  TCGA-HNSC_community_08_KEGG.csv:
  Shape: (50, 10)
  Columns: ['Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Old P-value', 'Old Adjusted P-value', 'Odds Ratio', 'Combined Score', 'Genes']
  Significant (Adjusted P-value < 0.05): 14
       Gene_set                                                          Term Overlap  P-value  Adjusted P-value  Old P-value  Old Adjusted P-value  Odds Ratio  Combined Score                 Genes
KEGG_2021_Human    Epithelial cell signaling in Helicobacter pylori infection    3/70 0.000022          0.000917            0                     0   68.600459      734.917073     CDC42;CXCL3;CXCL2
KEGG_2021_Human                                       IL-17 signaling pathway    3/94 0.000054          0.000917            0                     0   50.447168      495.852000      CSF3;CXCL3;CXCL2
KEGG_2021_Human Viral protein interaction with cytokine and cytokine receptor   3/100 0.000065          0.000917            0                     0   47.312450      456.304277      IL34;CXCL3;CXCL2
KEGG_2021_Human                        Cytokine-cytokine receptor interaction   4/295 0.000073          0.000917            0                     0   22.557847      214.746109 CSF3;IL34;CXCL3;CXCL2
KEGG_2021_Human                                   Chemokine signaling pathway   3/192 0.000445          0.004448            0                     0   24.169719      186.538209     CDC42;CXCL3;CXCL2
KEGG_2021_Human                                     Lipid and atherosclerosis   3/215 0.000619          0.005155            0                     0   21.522496      159.007431     CDC42;CXCL3;CXCL2
KEGG_2021_Human                                                 Legionellosis    2/57 0.000933          0.006460            0                     0   51.763636      361.140127           CXCL3;CXCL2
KEGG_2021_Human                                            Mineral absorption    2/60 0.001034          0.006460            0                     0   49.078818      337.398903             MT2A;MT1E
KEGG_2021_Human                                          Rheumatoid arthritis    2/93 0.002460          0.013668            0                     0   31.229199      187.609523           CXCL3;CXCL2
KEGG_2021_Human                                                    Amoebiasis   2/102 0.002950          0.013929            0                     0   28.405714      165.491408           CXCL3;CXCL2
KEGG_2021_Human                                  NF-kappa B signaling pathway   2/104 0.003064          0.013929            0                     0   27.845938      161.169364           CXCL3;CXCL2
KEGG_2021_Human                                         TNF signaling pathway   2/112 0.003543          0.014764            0                     0   25.810390      145.640570           CXCL3;CXCL2
KEGG_2021_Human                           NOD-like receptor signaling pathway   2/181 0.008992          0.034586            0                     0   15.806065       74.468559           CXCL3;CXCL2
KEGG_2021_Human               Kaposi sarcoma-associated herpesvirus infection   2/193 0.010171          0.036324            0                     0   14.804039       67.924456           CXCL3;CXCL2

  TCGA-HNSC_community_09_KEGG.csv:
  Shape: (9, 10)
  Columns: ['Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Old P-value', 'Old Adjusted P-value', 'Odds Ratio', 'Combined Score', 'Genes']
  Significant (Adjusted P-value < 0.05): 6
       Gene_set                                       Term Overlap  P-value  Adjusted P-value  Old P-value  Old Adjusted P-value  Odds Ratio  Combined Score   Genes
KEGG_2021_Human Synthesis and degradation of ketone bodies    1/10 0.005488          0.039138            0                     0  222.000000     1155.569978  HMGCS2
KEGG_2021_Human            Terpenoid backbone biosynthesis    1/22 0.012037          0.039138            0                     0   95.085714      420.260566  HMGCS2
KEGG_2021_Human                       Butanoate metabolism    1/28 0.015296          0.039138            0                     0   73.933333      309.051915  HMGCS2
KEGG_2021_Human           Mucin type O-glycan biosynthesis    1/36 0.019627          0.039138            0                     0   57.011429      224.101952 GALNTL6
KEGG_2021_Human       Other types of O-glycan biosynthesis    1/47 0.025555          0.039138            0                     0   43.354348      158.977812 GALNTL6
KEGG_2021_Human Valine, leucine and isoleucine degradation    1/48 0.026092          0.039138            0                     0   42.429787      154.704793  HMGCS2

  TCGA-HNSC_community_10_KEGG.csv:
  Shape: (4, 10)
  Columns: ['Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Old P-value', 'Old Adjusted P-value', 'Odds Ratio', 'Combined Score', 'Genes']
  Significant (Adjusted P-value < 0.05): 0

  TCGA-HNSC_community_12_KEGG.csv:
  Shape: (6, 10)
  Columns: ['Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Old P-value', 'Old Adjusted P-value', 'Odds Ratio', 'Combined Score', 'Genes']
  Significant (Adjusted P-value < 0.05): 6
       Gene_set                                    Term Overlap  P-value  Adjusted P-value  Old P-value  Old Adjusted P-value  Odds Ratio  Combined Score  Genes
KEGG_2021_Human            Systemic lupus erythematosus   1/135 0.020114          0.033701            0                     0   74.115672      289.519228 H2AC12
KEGG_2021_Human                                Ribosome   1/158 0.023514          0.033701            0                     0   63.184713      236.952026  RPL28
KEGG_2021_Human                             Necroptosis   1/159 0.023662          0.033701            0                     0   62.781646      235.047508 H2AC12
KEGG_2021_Human                              Alcoholism   1/186 0.027643          0.033701            0                     0   53.545946      192.144266 H2AC12
KEGG_2021_Human Neutrophil extracellular trap formation   1/189 0.028084          0.033701            0                     0   52.683511      188.214474 H2AC12
KEGG_2021_Human                     Coronavirus disease   1/232 0.034399          0.034399            0                     0   42.783550      144.168439  RPL28

  TCGA-HNSC_community_13_KEGG.csv:
  Shape: (8, 10)
  Columns: ['Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Old P-value', 'Old Adjusted P-value', 'Odds Ratio', 'Combined Score', 'Genes']
  Significant (Adjusted P-value < 0.05): 8
       Gene_set                                                      Term Overlap  P-value  Adjusted P-value  Old P-value  Old Adjusted P-value  Odds Ratio  Combined Score Genes
KEGG_2021_Human Endocrine and other factor-regulated calcium reabsorption    1/53 0.007929          0.025306            0                     0  191.778846      927.671855  DNM2
KEGG_2021_Human                    Bacterial invasion of epithelial cells    1/77 0.011506          0.025306            0                     0  131.059211      585.163932  DNM2
KEGG_2021_Human                                    Synaptic vesicle cycle    1/78 0.011655          0.025306            0                     0  129.350649      575.872814  DNM2
KEGG_2021_Human                          Fc gamma R-mediated phagocytosis    1/97 0.014480          0.025306            0                     0  103.651042      438.959797  DNM2
KEGG_2021_Human       Parathyroid hormone synthesis, secretion and action   1/106 0.015817          0.025306            0                     0   94.723810      392.791175 MMP25
KEGG_2021_Human                         Phospholipase D signaling pathway   1/148 0.022037          0.029383            0                     0   67.517007      257.579384  DNM2
KEGG_2021_Human                                      Salmonella infection   1/249 0.036889          0.037327            0                     0   39.816532      131.388757  DNM2
KEGG_2021_Human                                               Endocytosis   1/252 0.037327          0.037327            0                     0   39.334661      129.333497  DNM2

  TCGA-HNSC_community_14_KEGG.csv:
  Shape: (8, 10)
  Columns: ['Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Old P-value', 'Old Adjusted P-value', 'Odds Ratio', 'Combined Score', 'Genes']
  Significant (Adjusted P-value < 0.05): 8
       Gene_set                                                     Term Overlap  P-value  Adjusted P-value  Old P-value  Old Adjusted P-value  Odds Ratio  Combined Score Genes
KEGG_2021_Human                                Pentose phosphate pathway    1/30 0.004493          0.019768            0                     0  344.275862     1860.860041 ALDOA
KEGG_2021_Human                          Fructose and mannose metabolism    1/33 0.004942          0.019768            0                     0  311.953125     1656.465408 ALDOA
KEGG_2021_Human                             Glycolysis / Gluconeogenesis    1/67 0.010017          0.026711            0                     0  150.992424      695.092860 ALDOA
KEGG_2021_Human                                  HIF-1 signaling pathway   1/109 0.016262          0.026731            0                     0   92.078704      379.266738 ALDOA
KEGG_2021_Human                                    TNF signaling pathway   1/112 0.016707          0.026731            0                     0   89.576577      366.541985   LIF
KEGG_2021_Human Signaling pathways regulating pluripotency of stem cells   1/143 0.021298          0.027548            0                     0   69.911972      269.101444   LIF
KEGG_2021_Human                               JAK-STAT signaling pathway   1/162 0.024105          0.027548            0                     0   61.602484      229.490639   LIF
KEGG_2021_Human                   Cytokine-cytokine receptor interaction   1/295 0.043602          0.043602            0                     0   33.508503      104.970143   LIF

  TCGA-HNSC_community_16_KEGG.csv:
  Shape: (1, 10)
  Columns: ['Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Old P-value', 'Old Adjusted P-value', 'Odds Ratio', 'Combined Score', 'Genes']
  Significant (Adjusted P-value < 0.05): 1
       Gene_set                   Term Overlap  P-value  Adjusted P-value  Old P-value  Old Adjusted P-value  Odds Ratio  Combined Score        Genes
KEGG_2021_Human Olfactory transduction   2/440 0.001428          0.001428            0                     0   89.310502      585.141746 OR52W1;OR5M1

  TCGA-HNSC_community_18_KEGG.csv:
  Shape: (3, 10)
  Columns: ['Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Old P-value', 'Old Adjusted P-value', 'Odds Ratio', 'Combined Score', 'Genes']
  Significant (Adjusted P-value < 0.05): 2
       Gene_set                       Term Overlap  P-value  Adjusted P-value  Old P-value  Old Adjusted P-value  Odds Ratio  Combined Score Genes
KEGG_2021_Human Oxytocin signaling pathway   1/154 0.022924          0.037264            0                     0   64.849673      244.845776  RGS2
KEGG_2021_Human cGMP-PKG signaling pathway   1/167 0.024842          0.037264            0                     0   59.731928      220.721435  RGS2

  TCGA-HNSC_community_44_KEGG.csv:
  Shape: (4, 10)
  Columns: ['Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Old P-value', 'Old Adjusted P-value', 'Odds Ratio', 'Combined Score', 'Genes']
  Significant (Adjusted P-value < 0.05): 4
       Gene_set                             Term Overlap  P-value  Adjusted P-value  Old P-value  Old Adjusted P-value  Odds Ratio  Combined Score  Genes
KEGG_2021_Human                         Ribosome   1/158 0.031230          0.045602            0                     0   42.121019      146.007788 RPS4Y1
KEGG_2021_Human                    RNA transport   1/186 0.036687          0.045602            0                     0   35.695495      117.985708 GEMIN8
KEGG_2021_Human Regulation of actin cytoskeleton   1/218 0.042895          0.045602            0                     0   30.382488       95.674279 TMSB4Y
KEGG_2021_Human              Coronavirus disease   1/232 0.045602          0.045602            0                     0   28.520924       88.067013 RPS4Y1

======================================================================
DIAGNOSTIC COMPLETE
======================================================================

======================================================================
HARRIS INTERACTORS & CELL-TYPE MARKERS IN FIGURE 2 NETWORK
======================================================================

Network partition: /master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_2/05_communities/TCGA-HNSC/TCGA-HNSC_best_partition.csv
Total network genes: 1131
Columns: ['gene', 'gene_symbol', 'community']
Unique gene symbols in network: 1131

======================================================================
HARRIS A3 INTERACTORS
======================================================================

Loaded: /master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_4/00_input/Harris_A3_interactors.txt
Shape: (174, 6)
Columns: ['gene_symbol', 'A3_baits', 'source', 'R_loop_associated', 'confirmed_coIP', 'A3B_interactor']
Harris genes (gene_symbol): 174

Harris interactors in network: 9 / 174

Matched genes:
  CASC3  (Community 1)
  COA7  (Community 3)
  HEATR3  (Community 2)
  MOV10  (Community 3)
  NXF1  (Community 2)
  RBM3  (Community 1)
  SCO2  (Community 35)
  SRSF5  (Community 10)
  YTHDC2  (Community 2)

A3B-specific interactors in network: 3 / 51
  COA7  (Community 3)
  MOV10  (Community 3)
  RBM3  (Community 1)

======================================================================
CELL-TYPE MARKERS
======================================================================

Total markers checked: 26
Markers in network: 1 / 26

Matched markers:
  CD14  (Macrophage/Monocyte, Community 0)

Markers NOT in network:
  Basal epithelial: KRT5, KRT14, KRT15, TACSTD2, TP63
  Suprabasal epithelial: KRT4, KRT13, KRT19, IVL
  Fibroblast: COL1A1, COL3A1, DCN, LUM
  Endothelial: PECAM1, VWF, CDH5
  T cell: CD3E, CD3D, CD8A
  B cell: CD19, MS4A1
  Macrophage/Monocyte: CD68, CSF1R
  NK cell: NKG7, GNLY

======================================================================
SUMMARY FOR SECTION 2 TEXT
======================================================================

  Harris A3 interactors in network: 9 / 174
  Cell-type markers in network: 1 / 26
  A3B degree in DIFF network: 1
  A3A degree in DIFF network: 0
  A3B community: 2 (134 genes)
  A3B community KEGG: no significant terms

======================================================================
DIAGNOSTIC COMPLETE
======================================================================
