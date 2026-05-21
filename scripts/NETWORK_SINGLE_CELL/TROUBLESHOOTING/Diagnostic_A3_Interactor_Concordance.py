#!/usr/bin/env python3
"""
Diagnostic_A3_Interactor_Concordance.py
========================================

Concordance-first path tracing for Harris A3 interactors.

Two concordance subgraphs per network:
  1. Activator subgraph: positive DIFF edges only, nodes with log2FC > 0
     (up in tumor, co-expression gained in tumor)
  2. Repressor subgraph: negative DIFF edges only, nodes with log2FC < 0
     (up in normal, co-expression maintained in normal)

From each concordant Harris interactor, flood-fills through the matching
concordance subgraph. If A3A or A3B is reachable via unbroken concordant
edges, runs Dijkstra to find the shortest concordant path. If not
reachable, reports the longest concordant chain and nearest approach.

Usage:
    conda run -n NETWORK python Diagnostic_A3_Interactor_Concordance.py
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


def log(msg):
    print(msg, flush=True)

def banner(title, char="="):
    print(f"\n{char * 80}\n  {title}\n{char * 80}", flush=True)


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

    # Partition
    part_path = os.path.join(net_dir, "04_communities/SC_best_partition.csv")
    part_df = pd.read_csv(part_path)
    result["gene_to_comm"] = dict(zip(part_df["gene"], part_df["community"]))

    # Community graph
    graph_path = os.path.join(net_dir, "04_communities/SC_G_comm.gpickle")
    with open(graph_path, "rb") as f:
        result["G_comm"] = pickle.load(f)

    # DE stats
    de_path = os.path.join(net_dir,
                           "02_differential_expression/SC_diffexpr_stats.csv")
    de_df = pd.read_csv(de_path)
    result["de_log2fc"] = dict(zip(de_df["gene"], de_df["log2FC"]))

    return result


# =============================================================================
# BUILD CONCORDANCE SUBGRAPHS
# =============================================================================

def build_concordance_subgraphs(G_sub, de_log2fc):
    """Build activator and repressor concordance subgraphs.

    Activator: only positive DIFF edges between nodes with log2FC > 0
    Repressor: only negative DIFF edges between nodes with log2FC <= 0

    A3 genes are included in both subgraphs regardless of DE direction
    (they are the targets we're trying to reach).
    """
    up_tumor = set()
    up_normal = set()
    for n in G_sub.nodes():
        fc = de_log2fc.get(n, 0)
        if fc > 0:
            up_tumor.add(n)
        else:
            up_normal.add(n)

    # A3 genes go into both (they are reachability targets)
    a3_in_graph = set(n for n in G_sub.nodes() if n in A3_SYMBOLS)

    # Activator subgraph: positive edges, up-in-tumor nodes + A3 targets
    G_act = nx.Graph()
    act_nodes = up_tumor | a3_in_graph
    for u, v, d in G_sub.edges(data=True):
        w = d.get("weight", 0)
        if w > 0 and u in act_nodes and v in act_nodes:
            G_act.add_edge(u, v, weight=w, distance=1.0 / max(w, 0.001))

    # Repressor subgraph: negative edges, up-in-normal nodes + A3 targets
    G_rep = nx.Graph()
    rep_nodes = up_normal | a3_in_graph
    for u, v, d in G_sub.edges(data=True):
        w = d.get("weight", 0)
        if w < 0 and u in rep_nodes and v in rep_nodes:
            G_rep.add_edge(u, v, weight=w,
                           distance=1.0 / max(abs(w), 0.001))

    return G_act, G_rep


# =============================================================================
# CLASSIFY INTERACTORS
# =============================================================================

def classify_interactor(gene, G_sub, de_log2fc):
    """Classify a Harris interactor by DE + edge concordance."""
    log2fc = de_log2fc.get(gene, 0)
    de_direction = "up_tumor" if log2fc > 0 else "up_normal"

    if gene not in G_sub:
        return None

    neighbors = list(G_sub.neighbors(gene))
    if len(neighbors) == 0:
        return {
            "gene": gene, "log2FC": log2fc, "de_direction": de_direction,
            "n_edges": 0, "n_pos_edges": 0, "n_neg_edges": 0,
            "pct_concordant_edges": 0.0, "classification": "isolated",
        }

    n_pos = sum(1 for nb in neighbors
                if G_sub[gene][nb].get("weight", 0) > 0)
    n_neg = len(neighbors) - n_pos
    n_total = len(neighbors)

    if de_direction == "up_tumor":
        pct_conc = n_pos / n_total
    else:
        pct_conc = n_neg / n_total

    if pct_conc >= 0.6:
        classification = ("concordant_activator" if de_direction == "up_tumor"
                          else "concordant_repressor")
    elif pct_conc <= 0.4:
        classification = "discordant"
    else:
        classification = "mixed"

    return {
        "gene": gene, "log2FC": log2fc, "de_direction": de_direction,
        "n_edges": n_total, "n_pos_edges": n_pos, "n_neg_edges": n_neg,
        "pct_concordant_edges": round(pct_conc, 3),
        "classification": classification,
    }


# =============================================================================
# CONCORDANT FLOOD FILL + DIJKSTRA
# =============================================================================

def trace_concordant_path(gene, G_conc, de_log2fc, harris_genes):
    """From a concordant Harris interactor, flood-fill through the
    concordance subgraph and find the shortest concordant path to
    A3A or A3B if reachable.

    Returns path info or None.
    """
    if gene not in G_conc:
        return None

    a3_targets = [g for g in ["APOBEC3A", "APOBEC3B"] if g in G_conc]
    if not a3_targets:
        return None

    # Flood fill: find all reachable nodes through concordant edges
    reachable = set(nx.node_connected_component(G_conc, gene))

    # Which A3 targets are reachable?
    reachable_a3 = [t for t in a3_targets if t in reachable]

    if reachable_a3:
        # Dijkstra to find shortest concordant path to nearest A3
        best_path = None
        best_dist = float("inf")
        best_target = None

        for target in reachable_a3:
            try:
                path = nx.shortest_path(G_conc, source=gene,
                                        target=target, weight="distance")
                dist = nx.shortest_path_length(G_conc, source=gene,
                                               target=target,
                                               weight="distance")
                if dist < best_dist:
                    best_path = path
                    best_dist = dist
                    best_target = target
            except nx.NetworkXNoPath:
                continue

        if best_path is None:
            return None

        # Annotate each step
        steps = []
        for i in range(len(best_path) - 1):
            u, v = best_path[i], best_path[i + 1]
            edge_w = G_conc[u][v].get("weight", 0)
            v_fc = de_log2fc.get(v, 0)
            steps.append({
                "from": u,
                "to": v,
                "edge_weight": edge_w,
                "to_log2FC": v_fc,
                "to_de_direction": "up_tumor" if v_fc > 0 else "up_normal",
                "to_is_a3": v in A3_SYMBOLS,
                "to_is_harris": v in harris_genes,
            })

        # Strongest edge along path (bottleneck)
        min_abs_w = min(abs(s["edge_weight"]) for s in steps)
        # Total path strength (product of weights)
        path_strength = 1.0
        for s in steps:
            path_strength *= abs(s["edge_weight"])

        return {
            "source": gene,
            "target": best_target,
            "target_alias": A3_ALIAS.get(best_target, best_target),
            "path": best_path,
            "path_length": len(best_path) - 1,
            "total_distance": best_dist,
            "steps": steps,
            "reaches_a3": True,
            "reachable_component_size": len(reachable),
            "min_edge_weight": min_abs_w,
            "path_strength": path_strength,
        }

    else:
        # A3 not reachable. Report the concordant component.
        # Find the farthest reachable node (longest concordant chain)
        # using BFS depth from the source gene
        depth = nx.single_source_shortest_path_length(G_conc, gene)
        max_depth = max(depth.values())
        farthest = [n for n, d in depth.items() if d == max_depth]

        return {
            "source": gene,
            "target": None,
            "target_alias": None,
            "path": None,
            "path_length": 0,
            "total_distance": None,
            "steps": [],
            "reaches_a3": False,
            "reachable_component_size": len(reachable),
            "max_concordant_depth": max_depth,
            "farthest_genes": farthest[:5],
        }


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def analyze_network(net_data, harris_genes):
    name = net_data["name"]
    label = net_data["label"]
    G_comm = net_data["G_comm"]
    gene_to_comm = net_data["gene_to_comm"]
    de_log2fc = net_data["de_log2fc"]

    banner(f"CONCORDANCE ANALYSIS: {label}")

    # Find A3 community
    a3_comm = None
    for g in ["APOBEC3A", "APOBEC3B"]:
        if g in gene_to_comm:
            a3_comm = gene_to_comm[g]
            break
    if a3_comm is None:
        log("  [SKIP] A3A/A3B not in partition"); return

    comm_nodes = [g for g, c in gene_to_comm.items() if c == a3_comm]
    G_sub = G_comm.subgraph(comm_nodes).copy()
    isolates = [n for n in G_sub.nodes() if G_sub.degree(n) == 0]
    G_sub.remove_nodes_from(isolates)

    log(f"  A3 community C{a3_comm}: {G_sub.number_of_nodes()} nodes, "
        f"{G_sub.number_of_edges()} edges")

    # Report A3 genes in community and their DE direction
    a3_in_comm = [g for g in A3_SYMBOLS if g in G_sub.nodes()]
    log(f"\n  A3 genes in community:")
    for g in sorted(a3_in_comm):
        fc = de_log2fc.get(g, 0)
        direction = "up_tumor" if fc > 0 else "up_normal"
        log(f"    {A3_ALIAS.get(g, g)}: log2FC={fc:.4f} ({direction})")

    # Build concordance subgraphs
    G_act, G_rep = build_concordance_subgraphs(G_sub, de_log2fc)
    log(f"\n  Activator subgraph: {G_act.number_of_nodes()} nodes, "
        f"{G_act.number_of_edges()} edges")
    log(f"  Repressor subgraph: {G_rep.number_of_nodes()} nodes, "
        f"{G_rep.number_of_edges()} edges")

    # Check A3 reachability in each subgraph
    for sg_name, sg in [("Activator", G_act), ("Repressor", G_rep)]:
        a3_in_sg = [g for g in ["APOBEC3A", "APOBEC3B"] if g in sg]
        if a3_in_sg:
            for a3g in a3_in_sg:
                comp = nx.node_connected_component(sg, a3g)
                log(f"    {sg_name}: {A3_ALIAS[a3g]} in component of "
                    f"{len(comp)} nodes")
        else:
            log(f"    {sg_name}: neither A3A nor A3B present")

    # Classify all Harris interactors
    harris_in_comm = [g for g in harris_genes
                      if g in G_sub.nodes() and g not in A3_SYMBOLS]
    log(f"\n  Harris interactors in community: {len(harris_in_comm)}")

    results = []
    for gene in sorted(harris_in_comm):
        result = classify_interactor(gene, G_sub, de_log2fc)
        if result is not None:
            results.append(result)

    # Summary
    classifications = {}
    for r in results:
        c = r["classification"]
        classifications[c] = classifications.get(c, 0) + 1

    log(f"\n  === CLASSIFICATION SUMMARY ===")
    for c in ["concordant_activator", "concordant_repressor",
              "mixed", "discordant", "isolated"]:
        count = classifications.get(c, 0)
        if count > 0:
            genes = [r["gene"] for r in results if r["classification"] == c]
            log(f"    {c}: {count}")
            log(f"      {', '.join(genes)}")

    # Save concordance table
    df_conc = pd.DataFrame([{
        "gene": r["gene"],
        "log2FC": round(r["log2FC"], 4),
        "de_direction": r["de_direction"],
        "n_edges": r["n_edges"],
        "n_pos_edges": r["n_pos_edges"],
        "n_neg_edges": r["n_neg_edges"],
        "pct_concordant_edges": r["pct_concordant_edges"],
        "classification": r["classification"],
    } for r in results])
    df_conc.sort_values("pct_concordant_edges", ascending=False,
                        inplace=True)
    conc_path = os.path.join(OUTPUT_DIR, f"{name}_harris_concordance.tsv")
    df_conc.to_csv(conc_path, sep="\t", index=False)
    log(f"\n  [SAVE] Concordance table -> {conc_path}")

    # =================================================================
    # CONCORDANT PATH TRACING (Dijkstra on filtered subgraphs)
    # =================================================================
    banner(f"CONCORDANT PATH TRACING: {label}")

    path_rows = []

    # --- Concordant activators -> activator subgraph ---
    activators = [r for r in results
                  if r["classification"] == "concordant_activator"]
    log(f"\n  Tracing {len(activators)} concordant activators through "
        f"activator subgraph...")

    for r in activators:
        gene = r["gene"]
        p = trace_concordant_path(gene, G_act, de_log2fc, harris_genes)
        if p is None:
            log(f"\n  {gene}: not in activator subgraph")
            continue

        if p["reaches_a3"]:
            path_str = " -> ".join(p["path"])
            log(f"\n  {gene} -> {p['target_alias']}: "
                f"CONCORDANT PATH FOUND ({p['path_length']} hops)")
            log(f"    Path: {path_str}")
            log(f"    Component size: {p['reachable_component_size']}, "
                f"min |edge|: {p['min_edge_weight']:.4f}, "
                f"path strength: {p['path_strength']:.6f}")

            for step in p["steps"]:
                markers = []
                if step["to_is_a3"]:
                    markers.append("A3")
                if step["to_is_harris"]:
                    markers.append("Harris")
                mk = f" [{', '.join(markers)}]" if markers else ""
                log(f"      {step['from']} --"
                    f"(w={step['edge_weight']:.4f})--> "
                    f"{step['to']} "
                    f"(log2FC={step['to_log2FC']:.4f}){mk}")

            path_rows.append({
                "source_gene": gene,
                "classification": "concordant_activator",
                "source_log2FC": round(r["log2FC"], 4),
                "target": p["target_alias"],
                "path_length": p["path_length"],
                "path": " -> ".join(p["path"]),
                "reaches_a3": True,
                "component_size": p["reachable_component_size"],
                "min_abs_edge": round(p["min_edge_weight"], 4),
                "path_strength": round(p["path_strength"], 6),
            })
        else:
            log(f"\n  {gene}: A3A/A3B NOT reachable in activator subgraph")
            log(f"    Concordant component: "
                f"{p['reachable_component_size']} nodes, "
                f"max depth: {p['max_concordant_depth']}")
            log(f"    Farthest: {', '.join(p['farthest_genes'])}")

            path_rows.append({
                "source_gene": gene,
                "classification": "concordant_activator",
                "source_log2FC": round(r["log2FC"], 4),
                "target": None,
                "path_length": 0,
                "path": None,
                "reaches_a3": False,
                "component_size": p["reachable_component_size"],
                "min_abs_edge": None,
                "path_strength": None,
            })

    # --- Concordant repressors -> repressor subgraph ---
    repressors = [r for r in results
                  if r["classification"] == "concordant_repressor"]
    log(f"\n  Tracing {len(repressors)} concordant repressors through "
        f"repressor subgraph...")

    for r in repressors:
        gene = r["gene"]
        p = trace_concordant_path(gene, G_rep, de_log2fc, harris_genes)
        if p is None:
            log(f"\n  {gene}: not in repressor subgraph")
            continue

        if p["reaches_a3"]:
            path_str = " -> ".join(p["path"])
            log(f"\n  {gene} -> {p['target_alias']}: "
                f"CONCORDANT PATH FOUND ({p['path_length']} hops)")
            log(f"    Path: {path_str}")
            log(f"    Component size: {p['reachable_component_size']}, "
                f"min |edge|: {p['min_edge_weight']:.4f}, "
                f"path strength: {p['path_strength']:.6f}")

            for step in p["steps"]:
                markers = []
                if step["to_is_a3"]:
                    markers.append("A3")
                if step["to_is_harris"]:
                    markers.append("Harris")
                mk = f" [{', '.join(markers)}]" if markers else ""
                log(f"      {step['from']} --"
                    f"(w={step['edge_weight']:.4f})--> "
                    f"{step['to']} "
                    f"(log2FC={step['to_log2FC']:.4f}){mk}")

            path_rows.append({
                "source_gene": gene,
                "classification": "concordant_repressor",
                "source_log2FC": round(r["log2FC"], 4),
                "target": p["target_alias"],
                "path_length": p["path_length"],
                "path": " -> ".join(p["path"]),
                "reaches_a3": True,
                "component_size": p["reachable_component_size"],
                "min_abs_edge": round(p["min_edge_weight"], 4),
                "path_strength": round(p["path_strength"], 6),
            })
        else:
            log(f"\n  {gene}: A3A/A3B NOT reachable in repressor subgraph")
            log(f"    Concordant component: "
                f"{p['reachable_component_size']} nodes, "
                f"max depth: {p['max_concordant_depth']}")
            log(f"    Farthest: {', '.join(p['farthest_genes'])}")

            path_rows.append({
                "source_gene": gene,
                "classification": "concordant_repressor",
                "source_log2FC": round(r["log2FC"], 4),
                "target": None,
                "path_length": 0,
                "path": None,
                "reaches_a3": False,
                "component_size": p["reachable_component_size"],
                "min_abs_edge": None,
                "path_strength": None,
            })

    # Save path results
    if path_rows:
        df_paths = pd.DataFrame(path_rows)
        df_paths.sort_values(["reaches_a3", "path_strength"],
                             ascending=[False, False], inplace=True)
        path_out = os.path.join(OUTPUT_DIR,
                                f"{name}_concordant_paths.tsv")
        df_paths.to_csv(path_out, sep="\t", index=False)
        log(f"\n  [SAVE] Path table -> {path_out}")

    # =================================================================
    # SUMMARY
    # =================================================================
    banner(f"SUMMARY: {label}")

    n_reach_a3 = sum(1 for p in path_rows if p["reaches_a3"])
    n_no_reach = sum(1 for p in path_rows if not p["reaches_a3"])
    log(f"  Total concordant interactors traced: {len(path_rows)}")
    log(f"    Reach A3A or A3B: {n_reach_a3}")
    log(f"    Do NOT reach A3A/A3B: {n_no_reach}")

    if n_reach_a3 > 0:
        reaching = [p for p in path_rows if p["reaches_a3"]]
        log(f"\n  === CONCORDANT PATHS TO A3 ===")
        for p in reaching:
            log(f"    [{p['classification']}] {p['source_gene']} -> "
                f"{p['target']}: {p['path_length']} hops, "
                f"strength={p['path_strength']}")
            log(f"      {p['path']}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    banner("A3 INTERACTOR CONCORDANCE DIAGNOSTIC V2")
    banner("Concordance-first Dijkstra on filtered subgraphs")

    harris_genes = load_harris()

    for net_config in NETWORKS:
        log(f"\nLoading {net_config['name']}...")
        net_data = load_network(net_config)
        analyze_network(net_data, harris_genes)

    banner("DIAGNOSTIC COMPLETE")
    log(f"\nOutput: {OUTPUT_DIR}")
    for f in sorted(os.listdir(OUTPUT_DIR)):
        log(f"  {f}")


if __name__ == "__main__":
    main()
