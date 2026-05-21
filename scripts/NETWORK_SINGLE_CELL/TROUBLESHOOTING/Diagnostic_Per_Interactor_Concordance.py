#!/usr/bin/env python3
"""
Diagnostic_Per_Interactor_Concordance.py
==========================================

For each Harris A3 interactor in the A3 community, independently flood
through concordant edges and report:
  1. Concordant component size and members
  2. Nearest node to A3A (hop count in full graph)
  3. Nearest node to A3B (hop count in full graph)
  4. Whether multiple interactors share a component

Separately reports A3-boundary-anchored chains (genes that extend from
A3 neighbors outward through concordant edges but are NOT anchored by
a Harris interactor).

Usage:
    conda run -n NETWORK python Diagnostic_Per_Interactor_Concordance.py
"""

import os
import pickle
import numpy as np
import pandas as pd
import networkx as nx

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
    return set(harris_df["gene_symbol"].values)


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
# BUILD CONCORDANCE SUBGRAPHS
# =============================================================================

def build_concordance_subgraphs(G_sub, de_log2fc):
    a3_nodes = set(n for n in G_sub.nodes() if n in A3_SYMBOLS)

    # Activator: positive DIFF edges between up-in-tumor nodes (+ A3)
    G_act = nx.Graph()
    act_eligible = set(n for n in G_sub.nodes()
                       if de_log2fc.get(n, 0) > 0) | a3_nodes
    for u, v, d in G_sub.edges(data=True):
        w = d.get("weight", 0)
        if w > 0 and u in act_eligible and v in act_eligible:
            G_act.add_edge(u, v, weight=w)
    for n in act_eligible:
        if n in G_sub.nodes():
            G_act.add_node(n)

    # Repressor: negative DIFF edges between up-in-normal nodes (+ A3)
    G_rep = nx.Graph()
    rep_eligible = set(n for n in G_sub.nodes()
                       if de_log2fc.get(n, 0) <= 0) | a3_nodes
    for u, v, d in G_sub.edges(data=True):
        w = d.get("weight", 0)
        if w < 0 and u in rep_eligible and v in rep_eligible:
            G_rep.add_edge(u, v, weight=w)
    for n in rep_eligible:
        if n in G_sub.nodes():
            G_rep.add_node(n)

    return G_act, G_rep


# =============================================================================
# FIND NEAREST APPROACH TO A3
# =============================================================================

def find_nearest_a3(component_nodes, G_full, a3_targets_in_graph):
    """From a set of nodes, find the nearest node to each A3 target
    using unweighted hop count in the full community graph."""
    results = {}
    for target in a3_targets_in_graph:
        best_node = None
        best_hops = float("inf")
        for node in component_nodes:
            if node == target:
                best_node = target
                best_hops = 0
                break
            try:
                hops = nx.shortest_path_length(G_full, node, target,
                                               weight=None)
                if hops < best_hops:
                    best_hops = hops
                    best_node = node
            except nx.NetworkXNoPath:
                continue
        alias = A3_ALIAS.get(target, target)
        results[alias] = {
            "nearest_node": best_node,
            "hops": best_hops if best_hops < float("inf") else None,
        }
    return results


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def analyze_network(net_data, harris_genes):
    name = net_data["name"]
    label = net_data["label"]
    G_comm = net_data["G_comm"]
    gene_to_comm = net_data["gene_to_comm"]
    de_log2fc = net_data["de_log2fc"]

    banner(f"PER-INTERACTOR CONCORDANCE: {label}")

    G_sub, a3_comm = get_a3_community(G_comm, gene_to_comm)
    if G_sub is None:
        log("  [SKIP] No A3 community"); return

    log(f"  A3 community C{a3_comm}: {G_sub.number_of_nodes()} nodes, "
        f"{G_sub.number_of_edges()} edges")

    a3_in_graph = [g for g in A3_TARGETS if g in G_sub.nodes()]
    harris_in_comm = sorted(g for g in harris_genes
                            if g in G_sub.nodes() and g not in A3_SYMBOLS)
    log(f"  A3 genes in community: {', '.join(A3_ALIAS[g] for g in a3_in_graph)}")
    log(f"  Harris interactors in community: {len(harris_in_comm)}")

    G_act, G_rep = build_concordance_subgraphs(G_sub, de_log2fc)
    log(f"  Activator subgraph: {G_act.number_of_edges()} edges")
    log(f"  Repressor subgraph: {G_rep.number_of_edges()} edges")

    # =====================================================================
    # PART 1: Per-interactor flooding
    # =====================================================================
    section("PART 1: Per-Harris-interactor concordance flooding")

    # Track which components we've already seen (to detect sharing)
    seen_components = {}  # frozenset -> list of interactors

    interactor_rows = []

    for gene in harris_in_comm:
        fc = de_log2fc.get(gene, 0)
        de_dir = "up_tumor" if fc > 0 else "up_normal"

        # Pick matching concordance subgraph
        if de_dir == "up_tumor":
            G_conc = G_act
            mode = "activator"
        else:
            G_conc = G_rep
            mode = "repressor"

        # Check if gene has concordant edges
        if gene not in G_conc or G_conc.degree(gene) == 0:
            interactor_rows.append({
                "gene": gene,
                "log2FC": round(fc, 4),
                "de_direction": de_dir,
                "mode": mode,
                "concordant_edges": 0,
                "component_size": 0,
                "component_members": "",
                "shared_with": "",
                "nearest_to_A3A": None,
                "hops_to_A3A": None,
                "nearest_to_A3B": None,
                "hops_to_A3B": None,
            })
            continue

        # Flood fill
        component = frozenset(nx.node_connected_component(G_conc, gene))
        conc_edges = G_conc.degree(gene)

        # Check for shared component
        if component in seen_components:
            seen_components[component].append(gene)
        else:
            seen_components[component] = [gene]

        # Remove A3 from member list for reporting (they're targets)
        members = sorted(component - set(A3_TARGETS))
        members_str = ", ".join(members)

        # Nearest approach to A3
        a3_dist = find_nearest_a3(component, G_sub, a3_in_graph)

        row = {
            "gene": gene,
            "log2FC": round(fc, 4),
            "de_direction": de_dir,
            "mode": mode,
            "concordant_edges": conc_edges,
            "component_size": len(component),
            "component_members": members_str,
            "shared_with": "",  # filled after loop
        }

        for a3g in A3_TARGETS:
            alias = A3_ALIAS[a3g]
            info = a3_dist.get(alias, {})
            row[f"nearest_to_{alias}"] = info.get("nearest_node")
            row[f"hops_to_{alias}"] = info.get("hops")

        interactor_rows.append(row)

        # Log
        log(f"\n  {gene} ({mode}, log2FC={fc:.3f}):")
        log(f"    Concordant edges from {gene}: {conc_edges}")
        log(f"    Component size: {len(component)}")
        log(f"    Members: {members_str[:120]}"
            f"{'...' if len(members_str) > 120 else ''}")

        for a3g in A3_TARGETS:
            alias = A3_ALIAS[a3g]
            info = a3_dist.get(alias, {})
            if info.get("hops") is not None:
                if info["hops"] == 0:
                    log(f"    -> {alias}: REACHED (in component)")
                else:
                    log(f"    -> {alias}: nearest node = "
                        f"{info['nearest_node']} "
                        f"({info['hops']} hops in full graph)")
            else:
                log(f"    -> {alias}: no path")

    # Fill in shared_with
    for row in interactor_rows:
        gene = row["gene"]
        for comp, genes_list in seen_components.items():
            if gene in genes_list and len(genes_list) > 1:
                others = [g for g in genes_list if g != gene]
                row["shared_with"] = ", ".join(others)

    # Summary: shared components
    section("Component sharing summary")
    for comp, genes_list in seen_components.items():
        if len(genes_list) > 1:
            log(f"  Shared component ({len(comp)} nodes): "
                f"{', '.join(genes_list)}")
        elif len(genes_list) == 1:
            log(f"  Independent component ({len(comp)} nodes): "
                f"{genes_list[0]}")

    # Interactors with no concordant edges
    no_conc = [r for r in interactor_rows if r["concordant_edges"] == 0]
    if no_conc:
        section("Interactors with NO concordant edges")
        for r in no_conc:
            log(f"  {r['gene']} ({r['mode']}, log2FC={r['log2FC']})")

    # Save interactor table
    df_int = pd.DataFrame(interactor_rows)
    df_int.sort_values(["mode", "component_size"], ascending=[True, False],
                       inplace=True)
    int_path = os.path.join(OUTPUT_DIR,
                            f"{name}_per_interactor_concordance.tsv")
    df_int.to_csv(int_path, sep="\t", index=False)
    log(f"\n  [SAVE] Per-interactor table -> {int_path}")

    # =====================================================================
    # PART 2: A3-boundary-anchored chains
    # =====================================================================
    section("PART 2: A3-boundary-anchored concordance chains")

    boundary_rows = []
    a3_set = set(a3_in_graph)

    for a3g in a3_in_graph:
        alias = A3_ALIAS[a3g]
        for nb in G_sub.neighbors(a3g):
            if nb in A3_SYMBOLS:
                continue

            nb_fc = de_log2fc.get(nb, 0)
            nb_de = "up_tumor" if nb_fc > 0 else "up_normal"
            edge_w = G_sub[a3g][nb].get("weight", 0)
            edge_dir = "positive" if edge_w > 0 else "negative"
            is_harris = nb in harris_genes

            # Is the edge to A3 concordant?
            edge_concordant = ((nb_de == "up_tumor" and edge_w > 0) or
                               (nb_de == "up_normal" and edge_w < 0))

            # Check concordant edges beyond A3
            if nb_de == "up_tumor":
                G_conc = G_act
                mode = "activator"
            else:
                G_conc = G_rep
                mode = "repressor"

            if nb in G_conc and G_conc.degree(nb) > 0:
                comp = set(nx.node_connected_component(G_conc, nb))
                comp_no_a3 = comp - a3_set
                conc_edges_beyond = sum(
                    1 for nn in G_conc.neighbors(nb)
                    if nn != a3g and nn not in A3_SYMBOLS)
            else:
                comp_no_a3 = set()
                conc_edges_beyond = 0

            if conc_edges_beyond == 0:
                continue  # not a boundary gene with chains beyond

            # Is this chain already covered by a Harris interactor?
            harris_in_comp = comp_no_a3 & harris_genes
            anchored_by_harris = len(harris_in_comp) > 0

            members = sorted(comp_no_a3)
            boundary_rows.append({
                "a3_gene": alias,
                "boundary_gene": nb,
                "boundary_log2FC": round(nb_fc, 4),
                "boundary_de": nb_de,
                "edge_to_a3_weight": round(edge_w, 4),
                "edge_to_a3_concordant": edge_concordant,
                "mode": mode,
                "concordant_edges_beyond": conc_edges_beyond,
                "chain_size": len(comp_no_a3),
                "chain_members": ", ".join(members[:20]),
                "is_harris": is_harris,
                "chain_has_harris": anchored_by_harris,
                "harris_in_chain": ", ".join(sorted(harris_in_comp)),
            })

            tag = "[Harris]" if is_harris else ""
            conc_tag = "CONCORDANT" if edge_concordant else "DISCORDANT"
            harris_tag = (f" (contains Harris: {', '.join(sorted(harris_in_comp))})"
                          if anchored_by_harris else " (no Harris)")
            log(f"\n  {alias} -> {nb} {tag}: "
                f"edge {edge_dir} (w={edge_w:.4f}), {conc_tag}")
            log(f"    {nb} is {nb_de} (log2FC={nb_fc:.3f}), "
                f"{mode} mode")
            log(f"    Chain beyond: {len(comp_no_a3)} nodes, "
                f"{conc_edges_beyond} concordant edges{harris_tag}")
            if members:
                log(f"    Members: {', '.join(members[:15])}"
                    f"{'...' if len(members) > 15 else ''}")

    if boundary_rows:
        df_bnd = pd.DataFrame(boundary_rows)
        df_bnd.sort_values(["a3_gene", "chain_size"], ascending=[True, False],
                           inplace=True)
        bnd_path = os.path.join(OUTPUT_DIR,
                                f"{name}_boundary_chain_detail.tsv")
        df_bnd.to_csv(bnd_path, sep="\t", index=False)
        log(f"\n  [SAVE] Boundary chain table -> {bnd_path}")

        # Summary
        section("Boundary chain summary")
        n_with_harris = sum(1 for r in boundary_rows if r["chain_has_harris"])
        n_without = sum(1 for r in boundary_rows if not r["chain_has_harris"])
        n_concordant = sum(1 for r in boundary_rows
                          if r["edge_to_a3_concordant"])
        n_discordant = sum(1 for r in boundary_rows
                          if not r["edge_to_a3_concordant"])
        log(f"  Total boundary chains: {len(boundary_rows)}")
        log(f"    With Harris interactor: {n_with_harris}")
        log(f"    Without Harris (boundary-only): {n_without}")
        log(f"    Concordant edge to A3: {n_concordant}")
        log(f"    Discordant edge to A3 (wall): {n_discordant}")

    # =====================================================================
    # OVERALL SUMMARY
    # =====================================================================
    banner(f"SUMMARY: {label}")

    active_interactors = [r for r in interactor_rows
                          if r["concordant_edges"] > 0]
    log(f"  Harris interactors with concordant chains: "
        f"{len(active_interactors)}/{len(interactor_rows)}")

    for r in sorted(active_interactors, key=lambda x: -x["component_size"]):
        shared = f" (shared with {r['shared_with']})" if r['shared_with'] else ""
        a3a_dist = r.get("hops_to_A3A")
        a3b_dist = r.get("hops_to_A3B")
        a3a_str = (f"A3A: {r['nearest_to_A3A']} @ {a3a_dist} hops"
                   if a3a_dist is not None else "A3A: no path")
        a3b_str = (f"A3B: {r['nearest_to_A3B']} @ {a3b_dist} hops"
                   if a3b_dist is not None else "A3B: no path")
        log(f"  {r['gene']} ({r['mode']}): "
            f"{r['component_size']} nodes{shared}")
        log(f"    {a3a_str}")
        log(f"    {a3b_str}")

    no_chain = [r for r in interactor_rows if r["concordant_edges"] == 0]
    if no_chain:
        log(f"\n  Harris interactors with NO concordant chain: {len(no_chain)}")
        for r in no_chain:
            log(f"    {r['gene']} ({r['mode']}, log2FC={r['log2FC']})")


# =============================================================================
# MAIN
# =============================================================================

def main():
    banner("PER-INTERACTOR CONCORDANCE DIAGNOSTIC")

    harris_genes = load_harris()
    log(f"Harris A3 interactors: {len(harris_genes)}")

    for net_config in NETWORKS:
        log(f"\nLoading {net_config['name']}...")
        net_data = load_network(net_config)
        analyze_network(net_data, harris_genes)

    banner("DIAGNOSTIC COMPLETE")
    log(f"\nOutput: {OUTPUT_DIR}")
    for f in sorted(os.listdir(OUTPUT_DIR)):
        if "per_interactor" in f or "boundary_chain" in f:
            log(f"  {f}")


if __name__ == "__main__":
    main()
