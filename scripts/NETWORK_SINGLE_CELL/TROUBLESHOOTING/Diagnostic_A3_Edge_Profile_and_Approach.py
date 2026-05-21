#!/usr/bin/env python3
"""
Diagnostic_A3_Edge_Profile_and_Approach.py
============================================

Three-part diagnostic for understanding A3 co-expression decoupling:

Q1: A3 Edge Profile
    For A3A and A3B in each network, report positive vs negative DIFF
    edges, and specifically which Harris interactors they connect to
    directly and by what sign. Determines whether the concordance wall
    at A3 is real or artifactual.

Q2: Concordant Chain Mapping with Closest Approach
    From each Harris interactor, flood-fill through concordant edges
    (up-in-tumor + positive DIFF, or up-in-normal + negative DIFF).
    Map the full extent of each concordant component. Then measure the
    shortest distance from the component boundary to A3A/A3B across
    the full (unfiltered) graph. Reports how close concordant chains
    get before the concordance breaks.

Q3: A3 Boundary Analysis
    Identifies the specific genes that sit at the boundary between
    concordant chains and A3, meaning the last concordant node
    before the edge sign or DE direction flips. These are the
    decoupling points.

Usage:
    conda run -n NETWORK python Diagnostic_A3_Edge_Profile_and_Approach.py
"""

import os
import pickle
import numpy as np
import pandas as pd
import networkx as nx
from collections import defaultdict

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG4_ROOT = os.path.join(BASE_DIR, "data/FIG_4")
HARRIS_PATH = os.path.join(FIG4_ROOT, "00_input/Harris_A3_interactors.txt")
OUTPUT_DIR = os.path.join(FIG4_ROOT, "DIAGNOSTIC_CONCORDANCE")
os.makedirs(OUTPUT_DIR, exist_ok=True)

NETWORKS = [
    {
        "name": "SBS2_VS_NORMAL",
        "label": "SBS2-HIGH vs NORMAL",
        "dir": os.path.join(FIG4_ROOT, "NETWORK_SBS2_VS_NORMAL"),
    },
    {
        "name": "CNV_VS_NORMAL",
        "label": "CNV-HIGH vs NORMAL",
        "dir": os.path.join(FIG4_ROOT, "NETWORK_CNV_VS_NORMAL"),
    },
]

A3_SYMBOLS = {"APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
              "APOBEC3F", "APOBEC3G", "APOBEC3H"}
A3_ALIAS = {"APOBEC3A": "A3A", "APOBEC3B": "A3B", "APOBEC3C": "A3C",
            "APOBEC3D": "A3D", "APOBEC3F": "A3F", "APOBEC3G": "A3G",
            "APOBEC3H": "A3H"}
A3_TARGETS = ["APOBEC3A", "APOBEC3B"]


def log(msg):
    print(msg, flush=True)

def banner(title, char="="):
    print(f"\n{char * 80}\n  {title}\n{char * 80}", flush=True)

def section(title):
    print(f"\n  --- {title} ---", flush=True)


# =============================================================================
# LOAD DATA
# =============================================================================

def load_harris():
    harris_df = pd.read_csv(HARRIS_PATH, sep="\t")
    harris_genes = set(harris_df["gene_symbol"].values)
    log(f"Harris A3 interactors: {len(harris_genes)} genes")
    return harris_genes


def load_network(net_config):
    net_dir = net_config["dir"]
    result = {"name": net_config["name"], "label": net_config["label"]}

    part_df = pd.read_csv(
        os.path.join(net_dir, "04_communities/SC_best_partition.csv"))
    result["gene_to_comm"] = dict(zip(part_df["gene"], part_df["community"]))

    with open(os.path.join(net_dir,
              "04_communities/SC_G_comm.gpickle"), "rb") as f:
        result["G_comm"] = pickle.load(f)

    de_df = pd.read_csv(
        os.path.join(net_dir,
                     "02_differential_expression/SC_diffexpr_stats.csv"))
    result["de_log2fc"] = dict(zip(de_df["gene"], de_df["log2FC"]))

    return result


def get_a3_community(G_comm, gene_to_comm):
    """Return subgraph of the A3A/A3B community, isolates removed."""
    a3_comm = None
    for g in A3_TARGETS:
        if g in gene_to_comm:
            a3_comm = gene_to_comm[g]
            break
    if a3_comm is None:
        return None, None

    comm_nodes = [g for g, c in gene_to_comm.items() if c == a3_comm]
    G_sub = G_comm.subgraph(comm_nodes).copy()
    isolates = [n for n in G_sub.nodes() if G_sub.degree(n) == 0]
    G_sub.remove_nodes_from(isolates)
    return G_sub, a3_comm


# =============================================================================
# Q1: A3 EDGE PROFILE
# =============================================================================

def q1_a3_edge_profile(G_sub, de_log2fc, harris_genes, name):
    """Detailed edge profile for A3A and A3B."""
    banner(f"Q1: A3 EDGE PROFILE ({name})")

    rows = []

    for a3_gene in A3_TARGETS:
        if a3_gene not in G_sub:
            log(f"  {A3_ALIAS[a3_gene]}: NOT in community")
            continue

        alias = A3_ALIAS[a3_gene]
        a3_fc = de_log2fc.get(a3_gene, 0)
        a3_de = "up_tumor" if a3_fc > 0 else "up_normal"

        neighbors = list(G_sub.neighbors(a3_gene))
        n_total = len(neighbors)

        pos_edges = []
        neg_edges = []
        for nb in neighbors:
            w = G_sub[a3_gene][nb].get("weight", 0)
            nb_fc = de_log2fc.get(nb, 0)
            nb_de = "up_tumor" if nb_fc > 0 else "up_normal"
            is_harris = nb in harris_genes
            is_a3 = nb in A3_SYMBOLS

            entry = {
                "neighbor": nb, "edge_weight": w,
                "nb_log2FC": nb_fc, "nb_de": nb_de,
                "is_harris": is_harris, "is_a3": is_a3,
            }
            if w > 0:
                pos_edges.append(entry)
            else:
                neg_edges.append(entry)

            rows.append({
                "a3_gene": alias,
                "a3_log2FC": round(a3_fc, 4),
                "neighbor": nb,
                "edge_weight": round(w, 4),
                "edge_direction": "positive" if w > 0 else "negative",
                "neighbor_log2FC": round(nb_fc, 4),
                "neighbor_de": nb_de,
                "is_harris": is_harris,
                "is_a3": is_a3,
            })

        # Summary
        log(f"\n  {alias}: log2FC={a3_fc:.4f} ({a3_de})")
        log(f"    Total edges: {n_total}")
        log(f"    Positive DIFF (gained in tumor): {len(pos_edges)}")
        log(f"    Negative DIFF (gained in normal): {len(neg_edges)}")
        pct_pos = len(pos_edges) / n_total * 100 if n_total > 0 else 0
        log(f"    Percent positive: {pct_pos:.1f}%")

        # Harris neighbors breakdown
        harris_pos = [e for e in pos_edges if e["is_harris"]]
        harris_neg = [e for e in neg_edges if e["is_harris"]]
        log(f"\n    Harris interactor neighbors: "
            f"{len(harris_pos) + len(harris_neg)} total")

        if harris_pos:
            log(f"    Harris + positive DIFF (concordant activator edge):")
            for e in sorted(harris_pos, key=lambda x: -x["edge_weight"]):
                log(f"      {e['neighbor']}: w={e['edge_weight']:.4f}, "
                    f"log2FC={e['nb_log2FC']:.4f} ({e['nb_de']})")

        if harris_neg:
            log(f"    Harris + negative DIFF (concordant repressor edge):")
            for e in sorted(harris_neg, key=lambda x: x["edge_weight"]):
                log(f"      {e['neighbor']}: w={e['edge_weight']:.4f}, "
                    f"log2FC={e['nb_log2FC']:.4f} ({e['nb_de']})")

        # Neighbor DE breakdown by edge sign
        section(f"{alias} edge-vs-DE cross-tabulation")
        pos_up = sum(1 for e in pos_edges if e["nb_de"] == "up_tumor")
        pos_dn = sum(1 for e in pos_edges if e["nb_de"] == "up_normal")
        neg_up = sum(1 for e in neg_edges if e["nb_de"] == "up_tumor")
        neg_dn = sum(1 for e in neg_edges if e["nb_de"] == "up_normal")
        log(f"    Positive edge + up_tumor:  {pos_up}  "
            f"<-- concordant activator")
        log(f"    Positive edge + up_normal: {pos_dn}  "
            f"<-- discordant")
        log(f"    Negative edge + up_tumor:  {neg_up}  "
            f"<-- discordant")
        log(f"    Negative edge + up_normal: {neg_dn}  "
            f"<-- concordant repressor")

    # Save full edge table
    if rows:
        df = pd.DataFrame(rows)
        out_path = os.path.join(OUTPUT_DIR,
                                f"{name}_A3_edge_profile.tsv")
        df.to_csv(out_path, sep="\t", index=False)
        log(f"\n  [SAVE] A3 edge profile -> {out_path}")


# =============================================================================
# Q2: CONCORDANT CHAIN MAPPING + CLOSEST APPROACH
# =============================================================================

def build_concordant_subgraph(G_sub, de_log2fc, mode="activator"):
    """Build concordance-filtered subgraph.

    Activator: positive DIFF edges between up-in-tumor nodes.
    Repressor: negative DIFF edges between up-in-normal nodes.
    A3 genes included in both as reachability targets.
    """
    a3_nodes = set(n for n in G_sub.nodes() if n in A3_SYMBOLS)

    if mode == "activator":
        eligible = set(n for n in G_sub.nodes()
                       if de_log2fc.get(n, 0) > 0) | a3_nodes
        G_conc = nx.Graph()
        for u, v, d in G_sub.edges(data=True):
            w = d.get("weight", 0)
            if w > 0 and u in eligible and v in eligible:
                G_conc.add_edge(u, v, weight=w,
                                distance=1.0 / max(w, 0.001))
    else:
        eligible = set(n for n in G_sub.nodes()
                       if de_log2fc.get(n, 0) <= 0) | a3_nodes
        G_conc = nx.Graph()
        for u, v, d in G_sub.edges(data=True):
            w = d.get("weight", 0)
            if w < 0 and u in eligible and v in eligible:
                G_conc.add_edge(u, v, weight=w,
                                distance=1.0 / max(abs(w), 0.001))

    # Add isolated eligible nodes so they show as "in subgraph but
    # no concordant edges"
    for n in eligible:
        if n in G_sub.nodes():
            G_conc.add_node(n)

    return G_conc


def q2_concordant_chains(G_sub, de_log2fc, harris_genes, name):
    """Map concordant chains from every Harris interactor and measure
    closest approach to A3A/A3B."""
    banner(f"Q2: CONCORDANT CHAIN MAPPING ({name})")

    G_act = build_concordant_subgraph(G_sub, de_log2fc, "activator")
    G_rep = build_concordant_subgraph(G_sub, de_log2fc, "repressor")

    log(f"  Activator subgraph: {G_act.number_of_nodes()} nodes, "
        f"{G_act.number_of_edges()} edges")
    log(f"  Repressor subgraph: {G_rep.number_of_nodes()} nodes, "
        f"{G_rep.number_of_edges()} edges")

    harris_in_comm = [g for g in harris_genes
                      if g in G_sub.nodes() and g not in A3_SYMBOLS]

    # Build distance graph on full (unfiltered) community for gap
    # measurement: distance = 1/|weight|
    G_full_dist = G_sub.copy()
    for u, v, d in G_full_dist.edges(data=True):
        d["distance"] = 1.0 / max(abs(d.get("weight", 0.001)), 0.001)

    all_rows = []

    for mode_name, G_conc, de_check in [
        ("activator", G_act, lambda fc: fc > 0),
        ("repressor", G_rep, lambda fc: fc <= 0),
    ]:
        section(f"{mode_name.upper()} chains")

        # Which Harris interactors qualify for this mode?
        qualified = [g for g in harris_in_comm
                     if de_check(de_log2fc.get(g, 0))]
        log(f"  Harris interactors qualifying as {mode_name}: "
            f"{len(qualified)}")

        for gene in sorted(qualified):
            fc = de_log2fc.get(gene, 0)

            # Is this gene connected in the concordant subgraph?
            if gene not in G_conc or G_conc.degree(gene) == 0:
                all_rows.append({
                    "gene": gene, "mode": mode_name,
                    "log2FC": round(fc, 4),
                    "concordant_edges": 0,
                    "concordant_component_size": 0,
                    "reaches_a3a": False, "reaches_a3b": False,
                    "path_to_a3a": None, "path_to_a3b": None,
                    "hops_to_a3a": None, "hops_to_a3b": None,
                    "closest_approach_gene": None,
                    "gap_hops_to_a3": None,
                })
                continue

            # Concordant component from this gene
            conc_component = set(nx.node_connected_component(G_conc, gene))
            conc_edges = G_conc.degree(gene)

            # Can we reach A3A or A3B within the concordant subgraph?
            reaches = {}
            paths = {}
            for target in A3_TARGETS:
                if target in conc_component:
                    try:
                        p = nx.shortest_path(G_conc, gene, target,
                                             weight="distance")
                        reaches[target] = True
                        paths[target] = p
                    except nx.NetworkXNoPath:
                        reaches[target] = False
                else:
                    reaches[target] = False

            # If not reachable: find closest approach
            closest_gene = None
            gap_hops = None

            if not any(reaches.values()):
                # From each node in the concordant component, find the
                # shortest path to any A3 target in the FULL graph
                best_gap = float("inf")
                best_boundary = None

                for boundary_node in conc_component:
                    for target in A3_TARGETS:
                        if target not in G_full_dist:
                            continue
                        try:
                            gap = nx.shortest_path_length(
                                G_full_dist, boundary_node, target,
                                weight=None)  # unweighted hops
                            if gap < best_gap:
                                best_gap = gap
                                best_boundary = boundary_node
                        except nx.NetworkXNoPath:
                            continue

                if best_boundary is not None and best_gap < float("inf"):
                    closest_gene = best_boundary
                    gap_hops = best_gap

            row = {
                "gene": gene, "mode": mode_name,
                "log2FC": round(fc, 4),
                "concordant_edges": conc_edges,
                "concordant_component_size": len(conc_component),
                "reaches_a3a": reaches.get("APOBEC3A", False),
                "reaches_a3b": reaches.get("APOBEC3B", False),
                "path_to_a3a": (" -> ".join(paths["APOBEC3A"])
                                if "APOBEC3A" in paths else None),
                "path_to_a3b": (" -> ".join(paths["APOBEC3B"])
                                if "APOBEC3B" in paths else None),
                "hops_to_a3a": (len(paths["APOBEC3A"]) - 1
                                if "APOBEC3A" in paths else None),
                "hops_to_a3b": (len(paths["APOBEC3B"]) - 1
                                if "APOBEC3B" in paths else None),
                "closest_approach_gene": closest_gene,
                "gap_hops_to_a3": gap_hops,
            }
            all_rows.append(row)

            # Log interesting cases
            if any(reaches.values()):
                targets_hit = [A3_ALIAS[t] for t, v in reaches.items() if v]
                log(f"\n  {gene} ({mode_name}): REACHES "
                    f"{', '.join(targets_hit)}")
                for t, p in paths.items():
                    log(f"    -> {A3_ALIAS[t]}: {' -> '.join(p)} "
                        f"({len(p)-1} hops)")
            elif conc_edges > 0:
                log(f"\n  {gene} ({mode_name}): {len(conc_component)} "
                    f"concordant nodes, {conc_edges} concordant edges")
                if closest_gene:
                    log(f"    Closest approach: {closest_gene} "
                        f"({gap_hops} hops to nearest A3 in full graph)")

    # Save
    if all_rows:
        df = pd.DataFrame(all_rows)
        df.sort_values(["mode", "concordant_component_size"],
                       ascending=[True, False], inplace=True)
        out_path = os.path.join(OUTPUT_DIR,
                                f"{name}_concordant_chains.tsv")
        df.to_csv(out_path, sep="\t", index=False)
        log(f"\n  [SAVE] Concordant chains -> {out_path}")

    return all_rows


# =============================================================================
# Q3: BOUNDARY ANALYSIS
# =============================================================================

def q3_boundary_analysis(G_sub, de_log2fc, harris_genes, name):
    """Identify decoupling boundary genes: nodes that sit between
    concordant chains and A3, where concordance breaks."""
    banner(f"Q3: A3 BOUNDARY / DECOUPLING ANALYSIS ({name})")

    for a3_gene in A3_TARGETS:
        if a3_gene not in G_sub:
            continue

        alias = A3_ALIAS[a3_gene]
        a3_fc = de_log2fc.get(a3_gene, 0)

        section(f"Direct neighbors of {alias}")

        neighbors = list(G_sub.neighbors(a3_gene))
        log(f"  {alias} has {len(neighbors)} direct neighbors")

        # For each neighbor, classify the edge and the neighbor
        boundary_genes = []
        for nb in neighbors:
            w = G_sub[a3_gene][nb].get("weight", 0)
            nb_fc = de_log2fc.get(nb, 0)
            nb_de = "up_tumor" if nb_fc > 0 else "up_normal"
            edge_dir = "positive" if w > 0 else "negative"
            is_harris = nb in harris_genes

            # Check if this neighbor is in a concordant chain
            # (does it have concordant edges to OTHER genes beyond A3?)
            nb_neighbors = [nn for nn in G_sub.neighbors(nb)
                            if nn != a3_gene and nn not in A3_SYMBOLS]
            concordant_beyond = 0
            for nn in nb_neighbors:
                nn_w = G_sub[nb][nn].get("weight", 0)
                nn_fc = de_log2fc.get(nn, 0)
                nn_de = "up_tumor" if nn_fc > 0 else "up_normal"

                if nb_de == "up_tumor" and nn_de == "up_tumor" and nn_w > 0:
                    concordant_beyond += 1
                elif nb_de == "up_normal" and nn_de == "up_normal" and nn_w < 0:
                    concordant_beyond += 1

            boundary_genes.append({
                "a3_gene": alias,
                "neighbor": nb,
                "edge_weight": round(w, 4),
                "edge_direction": edge_dir,
                "neighbor_log2FC": round(nb_fc, 4),
                "neighbor_de": nb_de,
                "is_harris": is_harris,
                "concordant_edges_beyond": concordant_beyond,
                "total_edges_beyond": len(nb_neighbors),
            })

        # Sort by concordant edges beyond (most connected = most likely
        # to be a real boundary into a concordant chain)
        boundary_genes.sort(key=lambda x: -x["concordant_edges_beyond"])

        # Report the boundary genes that connect A3 to concordant chains
        log(f"\n  Genes bridging {alias} to concordant chains:")
        log(f"  (sorted by concordant edges beyond {alias})\n")

        active_bridges = [b for b in boundary_genes
                          if b["concordant_edges_beyond"] > 0]
        if active_bridges:
            for b in active_bridges[:20]:
                harris_mark = " [Harris]" if b["is_harris"] else ""
                log(f"    {b['neighbor']}: "
                    f"edge to {alias}={b['edge_direction']} "
                    f"(w={b['edge_weight']:.4f}), "
                    f"DE={b['neighbor_de']} "
                    f"(log2FC={b['neighbor_log2FC']:.4f}), "
                    f"concordant beyond={b['concordant_edges_beyond']}/"
                    f"{b['total_edges_beyond']}{harris_mark}")

                # Is the edge to A3 concordant or discordant?
                if b["neighbor_de"] == "up_tumor" and \
                   b["edge_direction"] == "positive":
                    log(f"      Edge to {alias}: CONCORDANT "
                        f"(up_tumor + positive DIFF)")
                elif b["neighbor_de"] == "up_normal" and \
                     b["edge_direction"] == "negative":
                    log(f"      Edge to {alias}: CONCORDANT "
                        f"(up_normal + negative DIFF)")
                else:
                    log(f"      Edge to {alias}: DISCORDANT "
                        f"({b['neighbor_de']} + {b['edge_direction']} DIFF) "
                        f"<-- DECOUPLING POINT")
        else:
            log(f"    No neighbors with concordant edges beyond {alias}")

        # Save boundary table
        if boundary_genes:
            df = pd.DataFrame(boundary_genes)
            out_path = os.path.join(
                OUTPUT_DIR,
                f"{name}_{alias}_boundary_genes.tsv")
            df.to_csv(out_path, sep="\t", index=False)
            log(f"\n  [SAVE] {alias} boundary -> {out_path}")

        # Summary statistics
        section(f"{alias} boundary summary")
        n_conc_bridges = len(active_bridges)
        n_conc_to_a3 = sum(
            1 for b in active_bridges
            if (b["neighbor_de"] == "up_tumor" and
                b["edge_direction"] == "positive") or
               (b["neighbor_de"] == "up_normal" and
                b["edge_direction"] == "negative"))
        n_decouple = n_conc_bridges - n_conc_to_a3

        log(f"  Total neighbors: {len(neighbors)}")
        log(f"  Neighbors with concordant chains beyond: {n_conc_bridges}")
        log(f"    Edge to {alias} is also concordant: {n_conc_to_a3} "
            f"(full concordant bridge)")
        log(f"    Edge to {alias} is discordant: {n_decouple} "
            f"(decoupling point)")
        log(f"  Neighbors with NO concordant chains beyond: "
            f"{len(boundary_genes) - n_conc_bridges}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    banner("A3 EDGE PROFILE + CONCORDANT APPROACH DIAGNOSTIC")

    harris_genes = load_harris()

    for net_config in NETWORKS:
        name = net_config["name"]
        log(f"\nLoading {name}...")
        net_data = load_network(net_config)

        G_sub, a3_comm = get_a3_community(
            net_data["G_comm"], net_data["gene_to_comm"])
        if G_sub is None:
            log("  [SKIP] No A3 community"); continue

        de_log2fc = net_data["de_log2fc"]

        log(f"  A3 community C{a3_comm}: {G_sub.number_of_nodes()} nodes, "
            f"{G_sub.number_of_edges()} edges")

        # Q1
        q1_a3_edge_profile(G_sub, de_log2fc, harris_genes, name)

        # Q2
        q2_concordant_chains(G_sub, de_log2fc, harris_genes, name)

        # Q3
        q3_boundary_analysis(G_sub, de_log2fc, harris_genes, name)

    banner("DIAGNOSTIC COMPLETE")
    log(f"\nOutput: {OUTPUT_DIR}")
    for f in sorted(os.listdir(OUTPUT_DIR)):
        log(f"  {f}")


if __name__ == "__main__":
    main()
