#!/usr/bin/env python3
"""
Summary_Figure4_Gene_Analysis.py
=================================

Comprehensive summary for Figure 4 results section preparation.

1. Catalogs all SC community genes with cross-references
2. Details the TCGA C2 (A3B+KRT5 basal) ↔ SC C0 (A3A+TACSTD2) overlap
3. Lists Known A3-Interactors found in SC communities with functional context
4. Generates a network plot of TCGA C2 highlighting SC C0 overlap genes
5. Produces a text summary for results section drafting

Usage:
  conda run -n NETWORK python Summary_Figure4_Gene_Analysis.py
"""

import os, json, pickle
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from datetime import datetime

# =============================================================================
# CONFIG
# =============================================================================
BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"

# SC data
SC_PARTITION   = os.path.join(BASE_DIR, "data/FIG_4/04_communities/SC_best_partition.csv")
SC_COMM_GENES  = os.path.join(BASE_DIR, "data/FIG_4/04_communities/SC_community_gene_lists.csv")
SC_DIFF_METRICS= os.path.join(BASE_DIR, "data/FIG_4/05_centrality_metrics/SC_DIFF_metrics.csv")
SC_GRAPH       = os.path.join(BASE_DIR, "data/FIG_4/04_communities/SC_G_comm.gpickle")
SC_DIFFEXPR    = os.path.join(BASE_DIR, "data/FIG_4/02_differential_expression/SC_diffexpr_stats.csv")
OVERLAP_DETAILS= os.path.join(BASE_DIR, "data/FIG_4/06_overlap_analysis/overlap_gene_details.csv")

# TCGA data
TCGA_PARTITION = os.path.join(BASE_DIR, "data/FIG_2/05_communities/TCGA-HNSC/TCGA-HNSC_best_partition.csv")
TCGA_GRAPH_DIR = os.path.join(BASE_DIR, "data/FIG_2/04_correlation_networks/TCGA-HNSC/graph_objects")
ENSG_TO_SYMBOL = os.path.join(BASE_DIR, "data/FIG_2/01_cleaned_expression/ensg_to_symbol.json")

# Harris interactors
HARRIS_ALL     = os.path.join(BASE_DIR, "data/FIG_4/00_input/Harris_A3_interactors.txt")
HARRIS_A3B     = os.path.join(BASE_DIR, "data/FIG_4/00_input/Harris_A3_interactors_A3B_only.txt")

# Output
OUTPUT_DIR     = os.path.join(BASE_DIR, "data/FIG_4/07_summary_analysis")
FIGURE_DIR     = os.path.join(BASE_DIR, "data/FIG_4/FIGURE_4_PANELS")

# Colors
COLOR_A3       = "#ed6a5a"
COLOR_OVERLAP  = "#f18f01"
COLOR_HARRIS   = "#fed766"
COLOR_BOTH     = "#ff6b6b"   # overlap + harris

A3_ENSG_TO_SYMBOL = {
    "ENSG00000128383": "APOBEC3A", "ENSG00000179750": "APOBEC3B",
    "ENSG00000244509": "APOBEC3C", "ENSG00000243811": "APOBEC3D",
    "ENSG00000128394": "APOBEC3F", "ENSG00000239713": "APOBEC3G",
    "ENSG00000100298": "APOBEC3H",
}

# Known functional context for key genes (for results section)
GENE_FUNCTIONS = {
    # Known A3-Interactors found in SC network
    "HNRNPA2B1": "RNA-binding protein; regulates mRNA splicing and transport; A3B interactor (McCann 2023 & Jang 2024); roles in RNA processing during transcription-associated mutagenesis",
    "HSPD1":     "Mitochondrial chaperonin (HSP60); protein folding; Jang 2024 A3 family interactor; involved in stress response and immune signaling",
    "RPL5":      "60S ribosomal protein L5; Jang 2024 A3 family interactor; also functions as MDM2 inhibitor activating p53 pathway",
    "TIMM8B":    "Mitochondrial import inner membrane translocase; Jang 2024 A3 family interactor; mitochondrial protein import",
    "RPL3":      "60S ribosomal protein L3; Jang 2024 A3 family interactor; ribosome biogenesis; p53-independent apoptosis regulator",
    # Key overlap genes (TCGA C2 ↔ SC C0)
    "TACSTD2":   "Trop-2; basal epithelial marker; overexpressed in HNSCC; therapeutic target (sacituzumab govitecan)",
    "KRT5":      "Keratin 5; canonical basal cell marker; Figure 2 anchor gene for basal community",
    "SCEL":      "Sciellin; cornified envelope protein; terminal differentiation marker for stratified epithelia",
    "CLIC3":     "Chloride intracellular channel 3; tumor microenvironment remodeling; prognostic in multiple cancers",
    "SULT2B1":   "Cholesterol sulfotransferase; lipid metabolism in epithelial differentiation",
    "AGR2":      "Anterior gradient 2; protein disulfide isomerase; promotes tumor growth in HNSCC",
    "S100A8":    "Calgranulin A; alarmin; innate immunity; overexpressed in HNSCC",
    "S100A9":    "Calgranulin B; alarmin; forms heterodimer with S100A8; tumor microenvironment",
    "IL1RN":     "IL-1 receptor antagonist; anti-inflammatory; modulates tumor immune response",
    "ALDH3B2":   "Aldehyde dehydrogenase 3B2; oxidative stress protection in epithelial cells",
    "KLK5":      "Kallikrein 5; serine protease; desquamation; epidermal barrier function",
    "KLK10":     "Kallikrein 10; serine protease; tumor suppressor in some contexts",
    "KLK11":     "Kallikrein 11; trypsin-like serine protease; biomarker in several cancers",
    "KLK13":     "Kallikrein 13; serine protease; epidermal homeostasis",
    "TGM1":      "Transglutaminase 1; crosslinks cornified envelope proteins; barrier function",
    "SPRR1A":    "Small proline-rich protein 1A; cornified envelope; epithelial stress response",
    "SPRR1B":    "Small proline-rich protein 1B; cornified envelope; squamous differentiation",
}


def banner(t):
    print(f"\n{'='*80}\n{t}\n{'='*80}", flush=True)

def log(m):
    print(m, flush=True)


def main():
    start = datetime.now()
    banner("FIGURE 4 — COMPREHENSIVE GENE ANALYSIS SUMMARY")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # =========================================================================
    # 1. Load all data
    # =========================================================================
    banner("1. LOADING DATA")

    # SC communities
    sc_part = pd.read_csv(SC_PARTITION)
    sc_g2c = dict(zip(sc_part["gene"], sc_part["community"]))
    sc_comms = {}
    for c in sorted(sc_part["community"].unique()):
        sc_comms[c] = set(sc_part[sc_part["community"]==c]["gene"])
    log(f"  SC communities: {len(sc_comms)}, total genes: {len(sc_g2c)}")

    # SC DE stats
    if os.path.exists(SC_DIFFEXPR):
        de_stats = pd.read_csv(SC_DIFFEXPR)
        gene_to_p = dict(zip(de_stats["gene"], de_stats["p_value"]))
        gene_to_fc = dict(zip(de_stats["gene"], de_stats["log2FC"]))
        log(f"  SC DE stats loaded: {len(de_stats)} genes")
    else:
        gene_to_p = {}; gene_to_fc = {}

    # SC centrality
    if os.path.exists(SC_DIFF_METRICS):
        cent = pd.read_csv(SC_DIFF_METRICS)
        gene_to_degree = dict(zip(cent["gene"], cent["degree"]))
        gene_to_hub = dict(zip(cent["gene"], cent["hub_rank"]))
        log(f"  SC centrality loaded: {len(cent)} genes")
    else:
        gene_to_degree = {}; gene_to_hub = {}

    # TCGA communities (ENSG → symbol)
    with open(ENSG_TO_SYMBOL) as f:
        e2s = json.load(f)
    tcga_part = pd.read_csv(TCGA_PARTITION)
    gc = "gene" if "gene" in tcga_part.columns else tcga_part.columns[0]
    tcga_comms = {}
    for c in sorted(tcga_part["community"].unique()):
        ensg_set = set(tcga_part[tcga_part["community"]==c][gc])
        symbols = set()
        for g in ensg_set:
            sym = e2s.get(g, A3_ENSG_TO_SYMBOL.get(g, g))
            symbols.add(sym)
        tcga_comms[c] = symbols
    log(f"  TCGA communities: {len(tcga_comms)}")
    for c, genes in sorted(tcga_comms.items()):
        log(f"    C{c}: {len(genes)} genes")

    # Harris interactors
    harris_all = set(); harris_a3b = set()
    for path, target in [(HARRIS_ALL, harris_all), (HARRIS_A3B, harris_a3b)]:
        if os.path.exists(path):
            with open(path) as f:
                for line in f:
                    g = line.strip().split("\t")[0].strip()
                    if g and not g.startswith("#"): target.add(g)
    log(f"  Harris ALL: {len(harris_all)}, A3B-only: {len(harris_a3b)}")

    # =========================================================================
    # 2. TCGA C2 ↔ SC C0 OVERLAP (key basal communities)
    # =========================================================================
    banner("2. TCGA C2 (A3B+KRT5) ↔ SC C0 (A3A+TACSTD2) OVERLAP")

    tcga_c2 = tcga_comms.get(2, set())
    sc_c0 = sc_comms.get(0, set())
    overlap_c2_c0 = tcga_c2 & sc_c0

    log(f"  TCGA C2: {len(tcga_c2)} genes (contains APOBEC3B, KRT5)")
    log(f"  SC C0:   {len(sc_c0)} genes (contains APOBEC3A, TACSTD2)")
    log(f"  Overlap: {len(overlap_c2_c0)} genes")
    log(f"\n  Overlapping genes:")
    for g in sorted(overlap_c2_c0):
        deg = gene_to_degree.get(g, "?")
        hub = gene_to_hub.get(g, "?")
        p = gene_to_p.get(g, None)
        fc = gene_to_fc.get(g, None)
        func = GENE_FUNCTIONS.get(g, "")
        p_str = f"p={p:.2e}" if p is not None else ""
        fc_str = f"FC={fc:.3f}" if fc is not None else ""
        in_harris = " [Known A3-Interactor]" if g in harris_all else ""
        log(f"    {g}: degree={deg}, hub_rank={hub}, {p_str}, {fc_str}{in_harris}")
        if func:
            log(f"      → {func}")

    # =========================================================================
    # 3. ALL TCGA ↔ SC OVERLAPS BY COMMUNITY PAIR
    # =========================================================================
    banner("3. ALL CROSS-COMMUNITY OVERLAPS")

    all_overlaps = []
    for sc_id, sc_genes in sorted(sc_comms.items()):
        for tcga_id, tcga_genes in sorted(tcga_comms.items()):
            overlap = sc_genes & tcga_genes
            if overlap:
                all_overlaps.append({
                    "SC_comm": sc_id, "TCGA_comm": tcga_id,
                    "overlap": len(overlap),
                    "sc_size": len(sc_genes), "tcga_size": len(tcga_genes),
                    "genes": sorted(overlap)
                })

    log(f"  Total community pairs with overlap: {len(all_overlaps)}")
    for ov in sorted(all_overlaps, key=lambda x: -x["overlap"]):
        log(f"  SC C{ov['SC_comm']} ({ov['sc_size']}) × TCGA C{ov['TCGA_comm']} ({ov['tcga_size']}): "
            f"{ov['overlap']} genes — {', '.join(ov['genes'][:10])}"
            f"{'...' if len(ov['genes'])>10 else ''}")

    # =========================================================================
    # 4. KNOWN A3-INTERACTORS IN SC NETWORK
    # =========================================================================
    banner("4. KNOWN A3-INTERACTORS IN SC COMMUNITIES")

    all_sc_genes = set()
    for g in sc_comms.values(): all_sc_genes |= g

    harris_in_sc = harris_all & all_sc_genes
    log(f"  Known A3-Interactors in SC network: {len(harris_in_sc)} / {len(harris_all)}")

    for g in sorted(harris_in_sc):
        comm = sc_g2c.get(g, "?")
        deg = gene_to_degree.get(g, "?")
        in_a3b = "A3B-specific" if g in harris_a3b else "A3-family"
        func = GENE_FUNCTIONS.get(g, "No annotation")
        log(f"  {g}: Community C{comm}, degree={deg}, {in_a3b}")
        log(f"    → {func}")

    # =========================================================================
    # 5. SC C0 DETAILED GENE TABLE
    # =========================================================================
    banner("5. SC COMMUNITY 0 — DETAILED GENE TABLE")

    c0_rows = []
    for g in sorted(sc_c0):
        in_tcga_c2 = g in tcga_c2
        in_any_tcga = any(g in tg for tg in tcga_comms.values())
        which_tcga = [f"C{c}" for c, tg in tcga_comms.items() if g in tg]
        in_harris_all = g in harris_all
        in_harris_a3b = g in harris_a3b

        c0_rows.append({
            "gene": g,
            "degree": gene_to_degree.get(g, 0),
            "hub_rank": gene_to_hub.get(g, 9999),
            "p_value": gene_to_p.get(g, 1.0),
            "log2FC": gene_to_fc.get(g, 0.0),
            "in_TCGA_C2": in_tcga_c2,
            "in_any_TCGA": in_any_tcga,
            "TCGA_communities": ";".join(which_tcga),
            "known_A3_interactor_all": in_harris_all,
            "known_A3_interactor_A3B": in_harris_a3b,
        })

    c0_df = pd.DataFrame(c0_rows).sort_values("hub_rank")
    c0_path = os.path.join(OUTPUT_DIR, "SC_C0_detailed_gene_table.csv")
    c0_df.to_csv(c0_path, index=False)
    log(f"  Saved: {c0_path}")
    log(f"  SC C0 genes in TCGA C2: {c0_df['in_TCGA_C2'].sum()}")
    log(f"  SC C0 genes in any TCGA community: {c0_df['in_any_TCGA'].sum()}")
    log(f"  SC C0 Known A3-Interactors: {c0_df['known_A3_interactor_all'].sum()}")

    # =========================================================================
    # 6. FULL CROSS-REFERENCE TABLE (all SC communities)
    # =========================================================================
    banner("6. FULL CROSS-REFERENCE TABLE (ALL SC COMMUNITIES)")

    full_rows = []
    for g, c in sorted(sc_g2c.items(), key=lambda x: (x[1], x[0])):
        which_tcga = [f"C{tc}" for tc, tg in tcga_comms.items() if g in tg]
        full_rows.append({
            "gene": g, "sc_community": c,
            "degree": gene_to_degree.get(g, 0),
            "hub_rank": gene_to_hub.get(g, 9999),
            "p_value": gene_to_p.get(g, 1.0),
            "log2FC": gene_to_fc.get(g, 0.0),
            "in_TCGA": len(which_tcga) > 0,
            "TCGA_communities": ";".join(which_tcga),
            "known_A3_interactor": g in harris_all,
            "A3_gene": g.startswith("APOBEC3"),
        })

    full_df = pd.DataFrame(full_rows)
    full_path = os.path.join(OUTPUT_DIR, "SC_all_communities_cross_reference.csv")
    full_df.to_csv(full_path, index=False)
    log(f"  Saved: {full_path}")

    # Summary stats
    log(f"\n  Cross-reference summary:")
    log(f"    Total SC community genes: {len(full_df)}")
    log(f"    In any TCGA community: {full_df['in_TCGA'].sum()} ({100*full_df['in_TCGA'].mean():.1f}%)")
    log(f"    Known A3-Interactors: {full_df['known_A3_interactor'].sum()}")
    log(f"    A3 genes: {full_df['A3_gene'].sum()}")

    for sc_c in sorted(full_df["sc_community"].unique()):
        sub = full_df[full_df["sc_community"]==sc_c]
        n_tcga = sub["in_TCGA"].sum()
        n_harris = sub["known_A3_interactor"].sum()
        n_a3 = sub["A3_gene"].sum()
        extras = []
        if n_tcga > 0: extras.append(f"{n_tcga} TCGA overlap")
        if n_harris > 0: extras.append(f"{n_harris} A3-interactor")
        if n_a3 > 0: extras.append(f"{n_a3} A3")
        extra_str = f" — {', '.join(extras)}" if extras else ""
        log(f"    SC C{sc_c}: {len(sub)} genes{extra_str}")

    # =========================================================================
    # 7. TCGA C2 NETWORK PLOT WITH SC C0 OVERLAP HIGHLIGHTED
    # =========================================================================
    banner("7. TCGA C2 NETWORK PLOT (highlighting SC C0 overlap)")

    # Try to load the TCGA DIFF graph
    tcga_diff_path = os.path.join(TCGA_GRAPH_DIR, "TCGA-HNSC_G_diff_noiso.gpickle")
    if os.path.exists(tcga_diff_path):
        with open(tcga_diff_path, "rb") as f:
            G_tcga = pickle.load(f)
        log(f"  Loaded TCGA DIFF graph: {G_tcga.number_of_nodes()} nodes, {G_tcga.number_of_edges()} edges")

        # Get TCGA C2 genes as ENSG
        tcga_c2_ensg = set(tcga_part[tcga_part["community"]==2][gc])
        # Subgraph of C2
        c2_in_graph = [n for n in tcga_c2_ensg if n in G_tcga.nodes()]
        if len(c2_in_graph) > 2:
            G_c2 = G_tcga.subgraph(c2_in_graph).copy()
            log(f"  TCGA C2 subgraph: {G_c2.number_of_nodes()} nodes, {G_c2.number_of_edges()} edges")

            # Map ENSG → symbol for display
            def sym(ensg):
                return e2s.get(ensg, A3_ENSG_TO_SYMBOL.get(ensg, ensg))

            # Identify special nodes
            overlap_ensg = set()
            harris_ensg = set()
            a3_ensg = set()
            for ensg in c2_in_graph:
                s = sym(ensg)
                if s in overlap_c2_c0:
                    overlap_ensg.add(ensg)
                if s in harris_all:
                    harris_ensg.add(ensg)
                if ensg in A3_ENSG_TO_SYMBOL:
                    a3_ensg.add(ensg)

            log(f"  In SC C0 overlap: {len(overlap_ensg)}")
            log(f"  Known A3-Interactors: {len(harris_ensg)}")
            log(f"  A3 genes: {len(a3_ensg)}")

            # Layout
            pos = nx.spring_layout(G_c2, seed=42, weight="abs_weight",
                                   k=2.0/np.sqrt(max(G_c2.number_of_nodes(),1)),
                                   iterations=200)

            fig, ax = plt.subplots(figsize=(18, 16))

            # Edges
            for u, v, d in G_c2.edges(data=True):
                w = d.get("weight", 0)
                color = "firebrick" if w > 0 else "steelblue"
                ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                        color=color, alpha=0.2, linewidth=0.5 + 2.0*abs(w))

            # Base nodes
            deg = dict(G_c2.degree())
            md = max(deg.values()) if deg else 1
            all_nodes = list(G_c2.nodes())
            non_special = [n for n in all_nodes if n not in overlap_ensg and n not in a3_ensg]
            special = [n for n in all_nodes if n in overlap_ensg or n in a3_ensg]

            # Non-overlapping nodes (gray)
            if non_special:
                ns = [80 + 400*(deg[n]/md) for n in non_special]
                nx.draw_networkx_nodes(G_c2, pos, nodelist=non_special, ax=ax,
                                       node_size=ns, node_color="#b0b0b0",
                                       alpha=0.6, edgecolors="black", linewidths=0.5)

            # Overlap nodes (orange)
            overlap_only = [n for n in overlap_ensg if n not in a3_ensg]
            if overlap_only:
                ns = [120 + 600*(deg[n]/md) for n in overlap_only]
                nx.draw_networkx_nodes(G_c2, pos, nodelist=overlap_only, ax=ax,
                                       node_size=ns, node_color=COLOR_OVERLAP,
                                       alpha=0.9, edgecolors="black", linewidths=1.5)

            # A3 nodes (red)
            if a3_ensg:
                a3_list = list(a3_ensg)
                ns = [200 + 800*(deg[n]/md) for n in a3_list]
                nx.draw_networkx_nodes(G_c2, pos, nodelist=a3_list, ax=ax,
                                       node_size=ns, node_color=COLOR_A3,
                                       alpha=1.0, edgecolors="black", linewidths=2.0)

            # Harris rings on overlap nodes
            harris_overlap = [n for n in overlap_ensg if sym(n) in harris_all]
            if harris_overlap:
                ns = [3.0*(120 + 600*(deg[n]/md)) for n in harris_overlap]
                nx.draw_networkx_nodes(G_c2, pos, nodelist=harris_overlap, ax=ax,
                                       node_size=ns, node_color="none",
                                       edgecolors=COLOR_HARRIS, linewidths=4.0)

            # Labels — all overlap + A3 + top 5 hubs
            labels = {}
            for n in list(overlap_ensg) + list(a3_ensg):
                labels[n] = sym(n)
            top_hubs = sorted(deg.items(), key=lambda x: -x[1])[:5]
            for n, d in top_hubs:
                if n not in labels:
                    labels[n] = sym(n)

            # Color labels
            label_colors = {}
            for n in labels:
                if n in a3_ensg:
                    label_colors[n] = COLOR_A3
                elif n in overlap_ensg:
                    label_colors[n] = "#d35400"
                else:
                    label_colors[n] = "black"

            for n, txt in labels.items():
                x, y = pos[n]
                col = label_colors.get(n, "black")
                ax.annotate(txt, (x, y), fontsize=12, fontweight="bold",
                            color=col, ha="center", va="bottom",
                            xytext=(0, 8), textcoords="offset points",
                            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8))

            # Legend
            legend_el = [
                Patch(fc="#b0b0b0", ec="black", lw=0.5, label=f"TCGA C2 only ({len(non_special)})"),
                Patch(fc=COLOR_OVERLAP, ec="black", lw=1.5, label=f"Also in SC C0 ({len(overlap_only)})"),
                Patch(fc=COLOR_A3, ec="black", lw=2, label=f"APOBEC3 ({len(a3_ensg)})"),
                Patch(fc="none", ec=COLOR_HARRIS, lw=3, label="Known A3-Interactor"),
            ]
            ax.legend(handles=legend_el, loc="upper left", fontsize=14, framealpha=0.9)
            ax.set_title(f"TCGA Bulk Community 2 (160 genes) — SC C0 Overlap Highlighted\n"
                         f"A3B + KRT5 basal cofactor network",
                         fontsize=18, pad=15)
            ax.axis("off")
            plt.tight_layout()

            for ext in ["pdf", "png"]:
                path = os.path.join(FIGURE_DIR, f"Panel_4f_TCGA_C2_with_SC_overlap.{ext}")
                plt.savefig(path, dpi=300, bbox_inches="tight")
            plt.close()
            log(f"  [SAVE] TCGA C2 overlap plot → {FIGURE_DIR}")
        else:
            log(f"  WARNING: Not enough TCGA C2 genes in graph ({len(c2_in_graph)})")
    else:
        log(f"  WARNING: TCGA DIFF graph not found at {tcga_diff_path}")

    # =========================================================================
    # 8. RESULTS SECTION SUMMARY
    # =========================================================================
    banner("8. RESULTS SECTION SUMMARY DATA")

    summary = []
    summary.append("=" * 70)
    summary.append("FIGURE 4 — RESULTS SECTION SUMMARY DATA")
    summary.append("=" * 70)
    summary.append("")
    summary.append("PIPELINE PARAMETERS:")
    summary.append(f"  Cell selection: L-method elbow (SBS2 >= 1.999), 546 HIGH / 546 LOW basal cells")
    summary.append(f"  DE genes: 9,244 (raw p < 0.05, 16,343 tested)")
    summary.append(f"  DIFF threshold: |delta-rho| >= 0.40")
    summary.append(f"  DIFF network: 1,254 nodes, 2,806 edges, LCC = 1,202 (95.9%)")
    summary.append(f"  Communities: 14 (Leiden, resolution 1.0, modularity 0.72)")
    summary.append(f"  Total community genes: 1,202")
    summary.append("")
    summary.append("KEY FINDINGS:")
    summary.append(f"  1. APOBEC3A in SC Community 0 (299 genes, hub_rank={gene_to_hub.get('APOBEC3A','?')}, degree={gene_to_degree.get('APOBEC3A','?')})")
    summary.append(f"     - Co-expressed with epithelial differentiation genes: SCEL, SULT2B1, CLIC3, KLK family")
    summary.append(f"     - Community also contains TACSTD2 (Trop-2, basal marker)")
    summary.append(f"  2. APOBEC3H in SC Community 8 (67 genes, hub_rank={gene_to_hub.get('APOBEC3H','?')}, degree={gene_to_degree.get('APOBEC3H','?')})")
    summary.append(f"  3. APOBEC3B max |delta-rho| = 0.26 — absent from communities (cofactors are post-transcriptional)")
    summary.append("")
    summary.append(f"CROSS-RESOLUTION OVERLAP:")
    summary.append(f"  Total SC genes in any TCGA community: {full_df['in_TCGA'].sum()} / {len(full_df)}")
    summary.append(f"  TCGA C2 (A3B+KRT5 basal) ↔ SC C0 (A3A+TACSTD2): {len(overlap_c2_c0)} shared genes")
    if overlap_c2_c0:
        summary.append(f"    Shared genes: {', '.join(sorted(overlap_c2_c0))}")
    summary.append("")
    summary.append(f"KNOWN A3-INTERACTORS IN SC NETWORK:")
    summary.append(f"  Total: {len(harris_in_sc)} / {len(harris_all)} ({100*len(harris_in_sc)/max(len(harris_all),1):.1f}%)")
    for g in sorted(harris_in_sc):
        comm = sc_g2c.get(g, "?")
        summary.append(f"    {g} — SC C{comm} — {GENE_FUNCTIONS.get(g, '')}")
    summary.append("")
    summary.append("BIOLOGICAL INTERPRETATION:")
    summary.append("  - SC C0 recapitulates key elements of TCGA C2 basal cofactor network at single-cell resolution")
    summary.append("  - A3A (not A3B) anchors the SC basal community, suggesting cell-type-specific A3 usage")
    summary.append("  - A3B's cofactor relationships are post-transcriptional (protein-protein, not co-expression)")
    summary.append("  - Known A3-interactors (HNRNPA2B1, RPL5, RPL3, HSPD1, TIMM8B) cluster in SC C0,")
    summary.append("    consistent with a shared A3-associated regulatory module in basal epithelial cells")
    summary.append("  - The KLK protease family (KLK5, KLK10, KLK11, KLK13) and cornified envelope genes")
    summary.append("    (SCEL, TGM1, SPRR1A/B) in SC C0 point to epithelial barrier disruption programs")
    summary.append("    differentially active in high-SBS2 basal cells")

    summary_text = "\n".join(summary)
    summary_path = os.path.join(OUTPUT_DIR, "Figure4_results_summary.txt")
    with open(summary_path, "w") as f:
        f.write(summary_text)
    print(summary_text)
    log(f"\n  [SAVE] Summary → {summary_path}")

    elapsed = datetime.now() - start
    banner(f"COMPLETE | Elapsed: {elapsed}")


if __name__ == "__main__":
    main()
