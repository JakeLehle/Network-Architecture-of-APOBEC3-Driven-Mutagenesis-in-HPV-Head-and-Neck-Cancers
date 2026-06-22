#!/usr/bin/env python3
"""
Diagnostic_A3_Interactor_Concordance.py  (V3)
=============================================

Concordance-first path tracing for Harris A3 interactors.

V3 changes (2026-06-21):
  1. MULTI-COMMUNITY. Analyzes EVERY community that contains A3A or A3B,
     not just the first one found. In CNV-HIGH vs NORMAL the two enzymes
     split (A3A and A3B land in different communities), so the old
     single-community pass silently skipped A3B's community (the
     productive-phase driver, where most Harris interactors sit).

  2. SELF-CONSISTENT CLASSIFICATION. A Harris interactor is called
     "concordant" only on edges that would also survive into the
     concordance subgraph: the edge sign matches the gene's DE direction
     AND the partner shares that direction (or the partner is an A3
     target). Previously the classifier counted any sign-matched edge
     while the subgraph additionally required the partner to match
     direction, so a gene could be classed concordant_activator yet be
     absent from the activator subgraph ("not in activator subgraph",
     traced 0). The classified concordant set is now guaranteed to be a
     subset of the traceable subgraph nodes.

  3. EXPLICIT A3 WALL VERDICT. For each A3 gene in each community, reports
     whether it appears in the activator and/or repressor subgraph, the
     size of its concordant component in each, how many concordant
     interactors reach it, and a verdict. An up-in-tumor A3 gene that is
     absent from the activator subgraph (its co-expression sits only on
     the negative-DIFF / repressor side) is flagged WALLED.

Concordance subgraphs per community:
  Activator: positive DIFF edges between up-in-tumor nodes (+ A3 targets)
  Repressor: negative DIFF edges between up-in-normal nodes (+ A3 targets)

From each concordant Harris interactor, flood-fills through the matching
concordance subgraph. If A3A or A3B is reachable via unbroken concordant
edges, runs Dijkstra to find the shortest concordant path; otherwise
reports the concordant component reached.

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

# A3 genes we trace toward (the DE-significant seeds in the partition)
A3_TARGETS = ["APOBEC3A", "APOBEC3B"]


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

    part_path = os.path.join(net_dir, "04_communities/SC_best_partition.csv")
    part_df = pd.read_csv(part_path)
    result["gene_to_comm"] = dict(zip(part_df["gene"], part_df["community"]))

    graph_path = os.path.join(net_dir, "04_communities/SC_G_comm.gpickle")
    with open(graph_path, "rb") as f:
        result["G_comm"] = pickle.load(f)

    de_path = os.path.join(net_dir,
                           "02_differential_expression/SC_diffexpr_stats.csv")
    de_df = pd.read_csv(de_path)
    result["de_log2fc"] = dict(zip(de_df["gene"], de_df["log2FC"]))

    return result


# =============================================================================
# DIRECTION + CONCORDANCE SUBGRAPHS
# =============================================================================

def node_directions(G_sub, de_log2fc):
    """Partition community nodes into up-in-tumor / up-in-normal by log2FC."""
    up_tumor, up_normal = set(), set()
    for n in G_sub.nodes():
        if de_log2fc.get(n, 0) > 0:
            up_tumor.add(n)
        else:
            up_normal.add(n)
    return up_tumor, up_normal


def build_concordance_subgraphs(G_sub, up_tumor, up_normal):
    """Build activator and repressor concordance subgraphs.

    Activator: positive DIFF edges between up-in-tumor nodes (+ A3 targets)
    Repressor: negative DIFF edges between up-in-normal nodes (+ A3 targets)

    A3 genes are included as endpoints in both subgraphs (they are the
    reachability targets), so a sign-matched edge into A3 is always kept.
    """
    a3_in_graph = set(n for n in G_sub.nodes() if n in A3_SYMBOLS)
    act_nodes = up_tumor | a3_in_graph
    rep_nodes = up_normal | a3_in_graph

    G_act = nx.Graph()
    for u, v, d in G_sub.edges(data=True):
        w = d.get("weight", 0)
        if w > 0 and u in act_nodes and v in act_nodes:
            G_act.add_edge(u, v, weight=w, distance=1.0 / max(w, 0.001))

    G_rep = nx.Graph()
    for u, v, d in G_sub.edges(data=True):
        w = d.get("weight", 0)
        if w < 0 and u in rep_nodes and v in rep_nodes:
            G_rep.add_edge(u, v, weight=w, distance=1.0 / max(abs(w), 0.001))

    return G_act, G_rep


# =============================================================================
# CLASSIFY INTERACTORS  (consistent with the concordance subgraphs)
# =============================================================================

def classify_interactor(gene, G_sub, de_log2fc, up_tumor, up_normal):
    """Classify a Harris interactor by DE direction + concordant edges.

    A "concordant edge" uses the SAME rule as the concordance subgraphs:
      up_tumor gene  -> positive edge to an up_tumor partner (or A3)
      up_normal gene -> negative edge to an up_normal partner (or A3)
    so any gene called concordant has >=1 edge in the matching subgraph.
    """
    log2fc = de_log2fc.get(gene, 0)
    de_direction = "up_tumor" if log2fc > 0 else "up_normal"

    if gene not in G_sub:
        return None

    neighbors = list(G_sub.neighbors(gene))
    n_total = len(neighbors)
    if n_total == 0:
        return {
            "gene": gene, "log2FC": log2fc, "de_direction": de_direction,
            "n_edges": 0, "n_pos_edges": 0, "n_neg_edges": 0,
            "n_concordant_edges": 0, "pct_concordant_edges": 0.0,
            "classification": "isolated",
        }

    n_pos = sum(1 for nb in neighbors
                if G_sub[gene][nb].get("weight", 0) > 0)
    n_neg = n_total - n_pos

    # Concordant edges: sign matches DE direction AND partner shares
    # direction (or partner is an A3 target). Matches subgraph membership.
    if de_direction == "up_tumor":
        n_conc = sum(1 for nb in neighbors
                     if G_sub[gene][nb].get("weight", 0) > 0
                     and (nb in up_tumor or nb in A3_SYMBOLS))
    else:
        n_conc = sum(1 for nb in neighbors
                     if G_sub[gene][nb].get("weight", 0) < 0
                     and (nb in up_normal or nb in A3_SYMBOLS))

    pct_conc = n_conc / n_total

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
        "n_concordant_edges": n_conc,
        "pct_concordant_edges": round(pct_conc, 3),
        "classification": classification,
    }


# =============================================================================
# CONCORDANT FLOOD FILL + DIJKSTRA
# =============================================================================

def trace_concordant_path(gene, G_conc, de_log2fc, harris_genes):
    """Flood-fill from a concordant interactor through the concordance
    subgraph; shortest concordant path to A3A/A3B if reachable."""
    if gene not in G_conc:
        return None

    a3_targets = [g for g in A3_TARGETS if g in G_conc]
    if not a3_targets:
        return None

    reachable = set(nx.node_connected_component(G_conc, gene))
    reachable_a3 = [t for t in a3_targets if t in reachable]

    if reachable_a3:
        best_path, best_dist, best_target = None, float("inf"), None
        for target in reachable_a3:
            try:
                path = nx.shortest_path(G_conc, source=gene,
                                        target=target, weight="distance")
                dist = nx.shortest_path_length(G_conc, source=gene,
                                               target=target,
                                               weight="distance")
                if dist < best_dist:
                    best_path, best_dist, best_target = path, dist, target
            except nx.NetworkXNoPath:
                continue
        if best_path is None:
            return None

        steps = []
        for i in range(len(best_path) - 1):
            u, v = best_path[i], best_path[i + 1]
            edge_w = G_conc[u][v].get("weight", 0)
            v_fc = de_log2fc.get(v, 0)
            steps.append({
                "from": u, "to": v, "edge_weight": edge_w, "to_log2FC": v_fc,
                "to_de_direction": "up_tumor" if v_fc > 0 else "up_normal",
                "to_is_a3": v in A3_SYMBOLS, "to_is_harris": v in harris_genes,
            })

        min_abs_w = min(abs(s["edge_weight"]) for s in steps)
        path_strength = 1.0
        for s in steps:
            path_strength *= abs(s["edge_weight"])

        return {
            "source": gene, "target": best_target,
            "target_alias": A3_ALIAS.get(best_target, best_target),
            "path": best_path, "path_length": len(best_path) - 1,
            "total_distance": best_dist, "steps": steps, "reaches_a3": True,
            "reachable_component_size": len(reachable),
            "min_edge_weight": min_abs_w, "path_strength": path_strength,
        }

    else:
        depth = nx.single_source_shortest_path_length(G_conc, gene)
        max_depth = max(depth.values())
        farthest = [n for n, d in depth.items() if d == max_depth]
        return {
            "source": gene, "target": None, "target_alias": None, "path": None,
            "path_length": 0, "total_distance": None, "steps": [],
            "reaches_a3": False, "reachable_component_size": len(reachable),
            "max_concordant_depth": max_depth, "farthest_genes": farthest[:5],
        }


# =============================================================================
# PER-COMMUNITY ANALYSIS
# =============================================================================

def _trace_group(results, classification, G_conc, subgraph_label,
                 comm, de_log2fc, harris_genes):
    """Trace one concordance class through its subgraph; return path rows."""
    group = [r for r in results if r["classification"] == classification]
    log(f"\n  Tracing {len(group)} {classification} through "
        f"{subgraph_label} subgraph...")
    rows = []
    for r in group:
        gene = r["gene"]
        p = trace_concordant_path(gene, G_conc, de_log2fc, harris_genes)
        if p is None:
            # With the V3 classifier this should not happen for concordant
            # genes; kept as a guard.
            log(f"\n  {gene}: not in {subgraph_label} subgraph (unexpected)")
            continue

        if p["reaches_a3"]:
            log(f"\n  {gene} -> {p['target_alias']}: CONCORDANT PATH FOUND "
                f"({p['path_length']} hops)")
            log(f"    Path: {' -> '.join(p['path'])}")
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
                log(f"      {step['from']} --(w={step['edge_weight']:.4f})--> "
                    f"{step['to']} (log2FC={step['to_log2FC']:.4f}){mk}")
            rows.append({
                "community": comm, "source_gene": gene,
                "classification": classification,
                "source_log2FC": round(r["log2FC"], 4),
                "target": p["target_alias"], "path_length": p["path_length"],
                "path": " -> ".join(p["path"]), "reaches_a3": True,
                "component_size": p["reachable_component_size"],
                "min_abs_edge": round(p["min_edge_weight"], 4),
                "path_strength": round(p["path_strength"], 6),
            })
        else:
            log(f"\n  {gene}: A3 NOT reachable in {subgraph_label} subgraph "
                f"(component {p['reachable_component_size']} nodes, "
                f"depth {p['max_concordant_depth']}, "
                f"farthest: {', '.join(p['farthest_genes'])})")
            rows.append({
                "community": comm, "source_gene": gene,
                "classification": classification,
                "source_log2FC": round(r["log2FC"], 4),
                "target": None, "path_length": 0, "path": None,
                "reaches_a3": False,
                "component_size": p["reachable_component_size"],
                "min_abs_edge": None, "path_strength": None,
            })
    return rows


def analyze_a3_community(net_data, comm, harris_genes):
    """Run the concordance analysis for a single A3-containing community.

    Returns (conc_rows, path_rows, wall_rows).
    """
    label = net_data["label"]
    G_comm = net_data["G_comm"]
    gene_to_comm = net_data["gene_to_comm"]
    de_log2fc = net_data["de_log2fc"]

    comm_nodes = [g for g, c in gene_to_comm.items() if c == comm]
    G_sub = G_comm.subgraph(comm_nodes).copy()
    isolates = [n for n in G_sub.nodes() if G_sub.degree(n) == 0]
    G_sub.remove_nodes_from(isolates)

    a3_in_comm = sorted(g for g in A3_TARGETS if g in G_sub.nodes())
    a3_lbl = ", ".join(A3_ALIAS[g] for g in a3_in_comm) if a3_in_comm else "none"

    banner(f"CONCORDANCE: {label} | C{comm} (A3: {a3_lbl})")
    log(f"  Community C{comm}: {G_sub.number_of_nodes()} nodes, "
        f"{G_sub.number_of_edges()} edges")

    up_tumor, up_normal = node_directions(G_sub, de_log2fc)

    log(f"\n  A3 genes in community:")
    for g in a3_in_comm:
        fc = de_log2fc.get(g, 0)
        d = "up_tumor" if fc > 0 else "up_normal"
        log(f"    {A3_ALIAS[g]}: log2FC={fc:.4f} ({d})")

    G_act, G_rep = build_concordance_subgraphs(G_sub, up_tumor, up_normal)
    log(f"\n  Activator subgraph: {G_act.number_of_nodes()} nodes, "
        f"{G_act.number_of_edges()} edges")
    log(f"  Repressor subgraph: {G_rep.number_of_nodes()} nodes, "
        f"{G_rep.number_of_edges()} edges")

    for sg_name, sg in [("Activator", G_act), ("Repressor", G_rep)]:
        present = [g for g in A3_TARGETS if g in sg]
        if present:
            for a3g in present:
                comp = nx.node_connected_component(sg, a3g)
                log(f"    {sg_name}: {A3_ALIAS[a3g]} in component of "
                    f"{len(comp)} nodes")
        else:
            log(f"    {sg_name}: neither A3A nor A3B present")

    # --- classify Harris interactors in this community ---
    harris_in_comm = [g for g in harris_genes
                      if g in G_sub.nodes() and g not in A3_SYMBOLS]
    log(f"\n  Harris interactors in community: {len(harris_in_comm)}")

    results = []
    for gene in sorted(harris_in_comm):
        r = classify_interactor(gene, G_sub, de_log2fc, up_tumor, up_normal)
        if r is not None:
            r["community"] = comm
            results.append(r)

    classifications = {}
    for r in results:
        classifications[r["classification"]] = \
            classifications.get(r["classification"], 0) + 1
    log(f"\n  === CLASSIFICATION SUMMARY ===")
    for c in ["concordant_activator", "concordant_repressor",
              "mixed", "discordant", "isolated"]:
        count = classifications.get(c, 0)
        if count > 0:
            genes = [r["gene"] for r in results if r["classification"] == c]
            log(f"    {c}: {count}")
            log(f"      {', '.join(genes)}")

    # --- trace concordant classes through their subgraphs ---
    banner(f"PATH TRACING: {label} | C{comm}", char="-")
    path_rows = []
    path_rows += _trace_group(results, "concordant_activator", G_act,
                              "activator", comm, de_log2fc, harris_genes)
    path_rows += _trace_group(results, "concordant_repressor", G_rep,
                              "repressor", comm, de_log2fc, harris_genes)

    # --- A3 wall verdict per A3 gene in this community ---
    wall_rows = []
    log(f"\n  === A3 WALL STATUS (C{comm}) ===")
    for a3g in a3_in_comm:
        fc = de_log2fc.get(a3g, 0)
        d = "up_tumor" if fc > 0 else "up_normal"
        in_act = a3g in G_act
        act_comp = len(nx.node_connected_component(G_act, a3g)) if in_act else 0
        in_rep = a3g in G_rep
        rep_comp = len(nx.node_connected_component(G_rep, a3g)) if in_rep else 0
        reached_by = [p["source_gene"] for p in path_rows
                      if p["reaches_a3"] and p["target"] == A3_ALIAS[a3g]]

        if d == "up_tumor":
            if in_act and reached_by:
                verdict = "COUPLED (activating chain reaches it)"
            elif in_act:
                verdict = "in activator subgraph; no concordant interactor path"
            else:
                verdict = ("WALLED (up in tumor but absent from activator "
                           "subgraph; co-expression only on negative-DIFF side)")
        else:
            if in_rep and reached_by:
                verdict = "COUPLED (repressor chain reaches it)"
            elif in_rep:
                verdict = "in repressor subgraph; no concordant interactor path"
            else:
                verdict = ("WALLED (up in normal but absent from repressor "
                           "subgraph)")

        log(f"    {A3_ALIAS[a3g]} ({d}): {verdict}")
        log(f"      activator: {'yes' if in_act else 'no'} "
            f"(component {act_comp}); "
            f"repressor: {'yes' if in_rep else 'no'} "
            f"(component {rep_comp}); "
            f"concordant interactors reaching it: {len(reached_by)}")
        if reached_by:
            log(f"      reached by: {', '.join(sorted(reached_by))}")

        wall_rows.append({
            "community": comm, "a3_gene": A3_ALIAS[a3g], "log2FC": round(fc, 4),
            "direction": d, "in_activator": in_act,
            "activator_component": act_comp, "in_repressor": in_rep,
            "repressor_component": rep_comp,
            "n_concordant_reaching": len(reached_by),
            "reached_by": ";".join(sorted(reached_by)) if reached_by else "",
            "verdict": verdict,
        })

    return results, path_rows, wall_rows


# =============================================================================
# PER-NETWORK ORCHESTRATION
# =============================================================================

def analyze_network(net_data, harris_genes):
    name = net_data["name"]
    label = net_data["label"]
    gene_to_comm = net_data["gene_to_comm"]

    banner(f"NETWORK: {label}")

    placement = {g: gene_to_comm[g] for g in A3_TARGETS if g in gene_to_comm}
    if not placement:
        log("  [SKIP] A3A/A3B not in partition"); return

    log("  A3 placement: " +
        ", ".join(f"{A3_ALIAS[g]}=C{c}" for g, c in placement.items()))
    a3_comms = sorted(set(placement.values()))
    if len(a3_comms) > 1:
        log("  NOTE: A3A and A3B are in different communities (split); "
            "analyzing each.")

    all_conc, all_paths, all_wall = [], [], []
    for comm in a3_comms:
        c_rows, p_rows, w_rows = analyze_a3_community(
            net_data, comm, harris_genes)
        all_conc += c_rows
        all_paths += p_rows
        all_wall += w_rows

    # --- save tables (community-tagged) ---
    if all_conc:
        df_conc = pd.DataFrame([{
            "community": r["community"], "gene": r["gene"],
            "log2FC": round(r["log2FC"], 4), "de_direction": r["de_direction"],
            "n_edges": r["n_edges"], "n_pos_edges": r["n_pos_edges"],
            "n_neg_edges": r["n_neg_edges"],
            "n_concordant_edges": r["n_concordant_edges"],
            "pct_concordant_edges": r["pct_concordant_edges"],
            "classification": r["classification"],
        } for r in all_conc])
        df_conc.sort_values(["community", "pct_concordant_edges"],
                            ascending=[True, False], inplace=True)
        p = os.path.join(OUTPUT_DIR, f"{name}_harris_concordance.tsv")
        df_conc.to_csv(p, sep="\t", index=False)
        log(f"\n  [SAVE] Concordance table -> {p}")

    if all_paths:
        df_paths = pd.DataFrame(all_paths)
        df_paths.sort_values(["community", "reaches_a3", "path_strength"],
                             ascending=[True, False, False], inplace=True)
        p = os.path.join(OUTPUT_DIR, f"{name}_concordant_paths.tsv")
        df_paths.to_csv(p, sep="\t", index=False)
        log(f"  [SAVE] Path table -> {p}")

    if all_wall:
        df_wall = pd.DataFrame(all_wall)
        p = os.path.join(OUTPUT_DIR, f"{name}_a3_wall_status.tsv")
        df_wall.to_csv(p, sep="\t", index=False)
        log(f"  [SAVE] A3 wall status -> {p}")

    # --- network summary ---
    banner(f"SUMMARY: {label}")
    n_reach = sum(1 for p in all_paths if p["reaches_a3"])
    log(f"  Communities analyzed: {len(a3_comms)} "
        f"({', '.join('C' + str(c) for c in a3_comms)})")
    log(f"  Concordant interactors traced: {len(all_paths)} "
        f"(reach A3: {n_reach})")
    log(f"  A3 wall status:")
    for w in all_wall:
        log(f"    C{w['community']} {w['a3_gene']} ({w['direction']}): "
            f"{w['verdict']}")
    if n_reach > 0:
        log(f"\n  === CONCORDANT PATHS TO A3 ===")
        for p in all_paths:
            if p["reaches_a3"]:
                log(f"    C{p['community']} [{p['classification']}] "
                    f"{p['source_gene']} -> {p['target']}: "
                    f"{p['path_length']} hops, strength={p['path_strength']}")
                log(f"      {p['path']}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    banner("A3 INTERACTOR CONCORDANCE DIAGNOSTIC V3")
    banner("Multi-community, self-consistent, with explicit wall verdict")

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
