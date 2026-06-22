#!/usr/bin/env python3
"""
Diagnostic_A3_Network_Structure.py > Still called Diagnostic_A3_Interactor_Concordance.py but this is the updated version
==================================

Structural map of the A3 neighborhood in each differential co-expression
network. Supersedes three earlier diagnostics by folding their useful
parts into one multi-community pass:

  - Diagnostic_A3_Interactor_Concordance.py  (V3; toward-A3 tracing)
  - Diagnostic_A3_Edge_Profile_and_Approach.py  (edge cross-tab, boundary)
  - Diagnostic_Per_Interactor_Concordance.py  (per-interactor / A3-anchored flood)

Why a rewrite. Tracing from each Harris interactor toward A3 always
reports zero when A3 is walled, which tells us nothing about the actual
structure. Instead this script ENUMERATES the concordant structure and
annotates it, and quantifies the wall directly from A3's edge composition.

MULTI-COMMUNITY. Every community containing A3A or A3B is analyzed
separately. In CNV-HIGH vs NORMAL the enzymes split (A3A in one
community, A3B in another), so each is profiled in its own community.

Concordance definitions (DIFF = HIGH_corr - LOW_corr):
  positive DIFF edge = co-expression gained in tumor (HIGH)
  negative DIFF edge = co-expression retained in normal (LOW)
  activator-concordant edge   = positive DIFF between two up-in-tumor genes
  repressor-concordant edge   = negative DIFF between two up-in-normal genes
  WALL edge                   = negative DIFF between two up-in-tumor genes
                                (co-induced but co-expression lost)

Per community it reports:
  PART A  community edge composition cross-tab, then per A3 gene an edge
          profile (sign x partner direction), direct Harris neighbors by
          sign, and a wall verdict.
  PART B  enumeration of every activator and repressor concordant
          subnetwork (size >= 2), each annotated with Harris members,
          A3 membership, top hubs, and hop-gap to the nearest A3 gene in
          the full community graph. Plus a per-Harris index.
  PART C  A3 boundary / decoupling: A3's direct neighbors classified as
          concordant bridges vs discordant decoupling points, and whether
          each neighbor anchors a concordant chain beyond A3.

Usage:
    conda run -n NETWORK python Diagnostic_A3_Network_Structure.py
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
    {"name": "SBS2_VS_NORMAL", "label": "SBS2-HIGH vs NORMAL",
     "dir": os.path.join(FIG4_ROOT, "NETWORK_SBS2_VS_NORMAL")},
    {"name": "CNV_VS_NORMAL", "label": "CNV-HIGH vs NORMAL",
     "dir": os.path.join(FIG4_ROOT, "NETWORK_CNV_VS_NORMAL")},
]

A3_SYMBOLS = {"APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
              "APOBEC3F", "APOBEC3G", "APOBEC3H"}
A3_ALIAS = {"APOBEC3A": "A3A", "APOBEC3B": "A3B", "APOBEC3C": "A3C",
            "APOBEC3D": "A3D", "APOBEC3F": "A3F", "APOBEC3G": "A3G",
            "APOBEC3H": "A3H"}
A3_TARGETS = ["APOBEC3A", "APOBEC3B"]

MIN_SUBNET_SIZE = 2     # smallest concordant subnetwork to enumerate
TOP_SUBNETS_PRINTED = 12


def log(msg):
    print(msg, flush=True)

def banner(title, char="="):
    print(f"\n{char * 80}\n  {title}\n{char * 80}", flush=True)

def section(title):
    print(f"\n  --- {title} ---", flush=True)


# =============================================================================
# LOAD
# =============================================================================

def load_harris():
    df = pd.read_csv(HARRIS_PATH, sep="\t")
    genes = set(df["gene_symbol"].values)
    log(f"Harris A3 interactors: {len(genes)} genes")
    return genes


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


# =============================================================================
# DIRECTIONS + CONCORDANCE SUBGRAPHS
# =============================================================================

def node_directions(G_sub, de_log2fc):
    up_tumor, up_normal = set(), set()
    for n in G_sub.nodes():
        if de_log2fc.get(n, 0) > 0:
            up_tumor.add(n)
        else:
            up_normal.add(n)
    return up_tumor, up_normal


def build_concordance_subgraphs(G_sub, up_tumor, up_normal):
    """Edge-only concordance subgraphs (A3 included as endpoints).

    Activator: positive DIFF edges between up-in-tumor nodes (+ A3).
    Repressor: negative DIFF edges between up-in-normal nodes (+ A3).
    A node appears only if it has >=1 qualifying edge.
    """
    a3 = set(n for n in G_sub.nodes() if n in A3_SYMBOLS)
    act_nodes = up_tumor | a3
    rep_nodes = up_normal | a3

    G_act = nx.Graph()
    G_rep = nx.Graph()
    for u, v, d in G_sub.edges(data=True):
        w = d.get("weight", 0)
        if w > 0 and u in act_nodes and v in act_nodes:
            G_act.add_edge(u, v, weight=w)
        elif w < 0 and u in rep_nodes and v in rep_nodes:
            G_rep.add_edge(u, v, weight=w)
    return G_act, G_rep


def edge_class(w, du, dv):
    """Classify a community edge by sign and the two endpoint directions."""
    if w > 0 and du == "up_tumor" and dv == "up_tumor":
        return "activator_concordant"
    if w < 0 and du == "up_normal" and dv == "up_normal":
        return "repressor_concordant"
    if w < 0 and du == "up_tumor" and dv == "up_tumor":
        return "wall"                      # co-induced, co-expression lost
    if w > 0 and du == "up_normal" and dv == "up_normal":
        return "anti_normal"               # co-normal, gained in tumor
    return "cross_pos" if w > 0 else "cross_neg"


# =============================================================================
# PART A: COMMUNITY EDGE COMPOSITION + A3 EDGE PROFILE
# =============================================================================

def part_a_edge_profile(G_sub, up_tumor, up_normal, de_log2fc,
                        harris, a3_in_comm, comm, name):
    dirof = lambda g: "up_tumor" if g in up_tumor else "up_normal"

    # --- community-level edge composition cross-tab ---
    section(f"Community C{comm} edge composition")
    counts = {}
    for u, v, d in G_sub.edges(data=True):
        c = edge_class(d.get("weight", 0), dirof(u), dirof(v))
        counts[c] = counts.get(c, 0) + 1
    total = sum(counts.values())
    log(f"    nodes: {G_sub.number_of_nodes()}, edges: {total}")
    log(f"    up_tumor nodes: {len(up_tumor & set(G_sub.nodes()))}, "
        f"up_normal nodes: {len(up_normal & set(G_sub.nodes()))}")
    for c in ["activator_concordant", "repressor_concordant", "wall",
              "anti_normal", "cross_pos", "cross_neg"]:
        n = counts.get(c, 0)
        pct = 100 * n / total if total else 0
        log(f"    {c:22s}: {n:7d}  ({pct:5.1f}%)")
    conc = counts.get("activator_concordant", 0) + \
        counts.get("repressor_concordant", 0)
    log(f"    -> concordant fraction: {100*conc/total:.1f}%  "
        f"(low fraction = wall-dominated community)")

    # --- per-A3-gene edge profile ---
    rows = []
    wall_rows = []
    for a3g in a3_in_comm:
        alias = A3_ALIAS[a3g]
        a3fc = de_log2fc.get(a3g, 0)
        a3dir = dirof(a3g)
        nbrs = list(G_sub.neighbors(a3g))
        ntot = len(nbrs)

        pos_up = pos_dn = neg_up = neg_dn = 0
        harris_pos, harris_neg = [], []
        for nb in nbrs:
            w = G_sub[a3g][nb].get("weight", 0)
            nbdir = dirof(nb)
            if w > 0 and nbdir == "up_tumor":
                pos_up += 1
            elif w > 0:
                pos_dn += 1
            elif w < 0 and nbdir == "up_tumor":
                neg_up += 1
            else:
                neg_dn += 1
            if nb in harris and nb not in A3_SYMBOLS:
                (harris_pos if w > 0 else harris_neg).append(
                    (nb, round(w, 3), round(de_log2fc.get(nb, 0), 2), nbdir))
            rows.append({
                "community": comm, "a3_gene": alias,
                "a3_log2FC": round(a3fc, 4), "neighbor": nb,
                "edge_weight": round(w, 4),
                "edge_direction": "positive" if w > 0 else "negative",
                "neighbor_log2FC": round(de_log2fc.get(nb, 0), 4),
                "neighbor_de": nbdir,
                "is_harris": (nb in harris and nb not in A3_SYMBOLS),
                "edge_class": edge_class(w, a3dir, nbdir),
            })

        section(f"{alias} edge profile (C{comm}, log2FC={a3fc:.3f}, {a3dir})")
        log(f"    total edges: {ntot}")
        log(f"    positive + up_tumor  : {pos_up:4d}  concordant activator")
        log(f"    negative + up_tumor  : {neg_up:4d}  WALL (co-induced, decohered)")
        log(f"    negative + up_normal : {neg_dn:4d}  concordant repressor")
        log(f"    positive + up_normal : {pos_dn:4d}  discordant")
        if harris_pos:
            log(f"    Harris (positive edge): "
                + ", ".join(f"{g}(w={w},fc={fc},{d})"
                            for g, w, fc, d in harris_pos))
        if harris_neg:
            log(f"    Harris (negative edge): "
                + ", ".join(f"{g}(w={w},fc={fc},{d})"
                            for g, w, fc, d in harris_neg))

        # wall verdict
        if a3dir == "up_tumor":
            if pos_up == 0:
                verdict = (f"WALLED: zero gained co-expression with the "
                           f"up-program; {neg_up} wall edges")
            elif pos_up < neg_up:
                verdict = (f"MOSTLY WALLED: {pos_up} concordant-activator vs "
                           f"{neg_up} wall edges")
            else:
                verdict = (f"PARTIALLY COUPLED: {pos_up} concordant-activator "
                           f"edges")
        else:
            verdict = "up_normal A3 (unexpected in tumor comparison)"
        log(f"    verdict: {verdict}")

        wall_rows.append({
            "community": comm, "a3_gene": alias, "log2FC": round(a3fc, 4),
            "direction": a3dir, "total_edges": ntot,
            "concordant_activator_edges": pos_up, "wall_edges": neg_up,
            "concordant_repressor_edges": neg_dn,
            "discordant_pos_edges": pos_dn,
            "n_harris_neighbors": len(harris_pos) + len(harris_neg),
            "verdict": verdict,
        })

    return rows, wall_rows


# =============================================================================
# PART B: SUBNETWORK ENUMERATION
# =============================================================================

def part_b_subnetworks(G_sub, G_act, G_rep, de_log2fc, harris,
                       a3_in_comm, comm, name):
    # hop distance from each A3 gene to every node in the full community
    a3_hops = {a3g: nx.single_source_shortest_path_length(G_sub, a3g)
               for a3g in a3_in_comm if a3g in G_sub}

    def gap_to_a3(component):
        present = [A3_ALIAS[a3g] for a3g in a3_in_comm if a3g in component]
        if present:
            return 0, present
        best = None
        for a3g, dist in a3_hops.items():
            for node in component:
                h = dist.get(node)
                if h is not None and (best is None or h < best):
                    best = h
        return best, []

    subnet_rows = []
    harris_index = []

    for mode, G_conc in [("activator", G_act), ("repressor", G_rep)]:
        comps = [c for c in nx.connected_components(G_conc)
                 if len(c) >= MIN_SUBNET_SIZE]
        comps.sort(key=len, reverse=True)
        section(f"{mode.upper()} subnetworks in C{comm}: {len(comps)} "
                f"(>= {MIN_SUBNET_SIZE} nodes)")

        for i, comp in enumerate(comps):
            sub = G_conc.subgraph(comp)
            hubs = sorted(sub.degree(), key=lambda x: -x[1])[:5]
            hub_str = ", ".join(f"{g}({d})" for g, d in hubs)
            h_members = sorted((comp & harris) - A3_SYMBOLS)
            a3_members = sorted(A3_ALIAS[g] for g in comp if g in A3_TARGETS)
            gap, contains = gap_to_a3(comp)

            tag = ""
            if a3_members:
                tag = f"  <-- contains {', '.join(a3_members)}"
            elif gap is not None:
                tag = f"  (nearest A3 {gap} hop{'s' if gap != 1 else ''} away)"
            if i < TOP_SUBNETS_PRINTED or a3_members or h_members:
                log(f"    [{mode[:3]}-{i}] {len(comp)} nodes, "
                    f"{sub.number_of_edges()} edges, "
                    f"Harris={len(h_members)}{tag}")
                log(f"         hubs: {hub_str}")
                if h_members:
                    log(f"         Harris: {', '.join(h_members)}")

            subnet_rows.append({
                "community": comm, "mode": mode, "subnet_id": f"{mode[:3]}-{i}",
                "size": len(comp), "n_edges": sub.number_of_edges(),
                "n_harris": len(h_members),
                "harris_members": ", ".join(h_members),
                "a3_members": ", ".join(a3_members),
                "top_hubs": hub_str,
                "gap_hops_to_a3": gap,
            })
            for hg in h_members:
                harris_index.append({
                    "community": comm, "harris_gene": hg, "mode": mode,
                    "subnet_id": f"{mode[:3]}-{i}", "subnet_size": len(comp),
                    "contains_a3": ", ".join(a3_members) if a3_members else "",
                    "gap_hops_to_a3": gap,
                })

    return subnet_rows, harris_index


# =============================================================================
# PART C: A3 BOUNDARY / DECOUPLING
# =============================================================================

def part_c_boundary(G_sub, G_act, G_rep, up_tumor, up_normal,
                    de_log2fc, harris, a3_in_comm, comm, name):
    dirof = lambda g: "up_tumor" if g in up_tumor else "up_normal"
    rows = []

    for a3g in a3_in_comm:
        alias = A3_ALIAS[a3g]
        a3dir = dirof(a3g)
        section(f"{alias} boundary (C{comm})")
        nbrs = [n for n in G_sub.neighbors(a3g) if n not in A3_SYMBOLS]

        n_bridge = n_decouple = 0
        for nb in nbrs:
            w = G_sub[a3g][nb].get("weight", 0)
            nbdir = dirof(nb)
            edge_conc = ((nbdir == "up_tumor" and w > 0) or
                         (nbdir == "up_normal" and w < 0))
            # does nb anchor a concordant chain beyond A3?
            G_conc = G_act if nbdir == "up_tumor" else G_rep
            beyond = 0
            if nb in G_conc:
                beyond = sum(1 for x in G_conc.neighbors(nb)
                             if x != a3g and x not in A3_SYMBOLS)
            if beyond == 0:
                continue
            if edge_conc:
                n_bridge += 1
                kind = "concordant bridge"
            else:
                n_decouple += 1
                kind = "DECOUPLING POINT (wall edge into a real chain)"
            rows.append({
                "community": comm, "a3_gene": alias, "boundary_gene": nb,
                "edge_weight": round(w, 4), "neighbor_de": nbdir,
                "neighbor_log2FC": round(de_log2fc.get(nb, 0), 4),
                "edge_to_a3_concordant": edge_conc,
                "concordant_edges_beyond": beyond,
                "is_harris": (nb in harris and nb not in A3_SYMBOLS),
                "kind": kind,
            })

        log(f"    neighbors anchoring concordant chains: "
            f"{n_bridge + n_decouple} "
            f"({n_bridge} concordant bridges, {n_decouple} decoupling points)")
        chain_rows = [r for r in rows
                      if r["a3_gene"] == alias and r["community"] == comm]
        for r in sorted(chain_rows,
                        key=lambda x: -x["concordant_edges_beyond"])[:15]:
            hm = " [Harris]" if r["is_harris"] else ""
            log(f"      {r['boundary_gene']}{hm}: edge w={r['edge_weight']}, "
                f"{r['neighbor_de']}, chains beyond="
                f"{r['concordant_edges_beyond']} -> {r['kind']}")
    return rows


# =============================================================================
# PER-COMMUNITY + PER-NETWORK
# =============================================================================

def analyze_a3_community(net_data, comm, harris):
    G_comm = net_data["G_comm"]
    gene_to_comm = net_data["gene_to_comm"]
    de_log2fc = net_data["de_log2fc"]
    name, label = net_data["name"], net_data["label"]

    comm_nodes = [g for g, c in gene_to_comm.items() if c == comm]
    G_sub = G_comm.subgraph(comm_nodes).copy()
    G_sub.remove_nodes_from([n for n in G_sub.nodes()
                             if G_sub.degree(n) == 0])

    a3_in_comm = sorted(g for g in A3_TARGETS if g in G_sub.nodes())
    a3_lbl = ", ".join(A3_ALIAS[g] for g in a3_in_comm) or "none"
    banner(f"{label} | C{comm} (A3: {a3_lbl})")
    log(f"  Community C{comm}: {G_sub.number_of_nodes()} nodes, "
        f"{G_sub.number_of_edges()} edges")

    up_tumor, up_normal = node_directions(G_sub, de_log2fc)
    G_act, G_rep = build_concordance_subgraphs(G_sub, up_tumor, up_normal)

    a_rows, wall_rows = part_a_edge_profile(
        G_sub, up_tumor, up_normal, de_log2fc, harris, a3_in_comm, comm, name)
    b_rows, h_index = part_b_subnetworks(
        G_sub, G_act, G_rep, de_log2fc, harris, a3_in_comm, comm, name)
    c_rows = part_c_boundary(
        G_sub, G_act, G_rep, up_tumor, up_normal, de_log2fc, harris,
        a3_in_comm, comm, name)

    return {"edges": a_rows, "wall": wall_rows, "subnets": b_rows,
            "harris_index": h_index, "boundary": c_rows}


def analyze_network(net_data, harris):
    name, label = net_data["name"], net_data["label"]
    gene_to_comm = net_data["gene_to_comm"]
    banner(f"NETWORK: {label}")

    placement = {g: gene_to_comm[g] for g in A3_TARGETS if g in gene_to_comm}
    if not placement:
        log("  [SKIP] A3A/A3B not in partition"); return
    log("  A3 placement: " +
        ", ".join(f"{A3_ALIAS[g]}=C{c}" for g, c in placement.items()))
    a3_comms = sorted(set(placement.values()))
    if len(a3_comms) > 1:
        log("  NOTE: A3A and A3B split across communities; analyzing each.")

    agg = {k: [] for k in
           ["edges", "wall", "subnets", "harris_index", "boundary"]}
    for comm in a3_comms:
        out = analyze_a3_community(net_data, comm, harris)
        for k in agg:
            agg[k] += out[k]

    # save tables
    saves = [
        ("edges", f"{name}_A3_edge_profile.tsv"),
        ("wall", f"{name}_A3_wall_status.tsv"),
        ("subnets", f"{name}_concordant_subnetworks.tsv"),
        ("harris_index", f"{name}_harris_subnetwork_index.tsv"),
        ("boundary", f"{name}_A3_boundary.tsv"),
    ]
    for key, fname in saves:
        if agg[key]:
            pd.DataFrame(agg[key]).to_csv(
                os.path.join(OUTPUT_DIR, fname), sep="\t", index=False)

    # network summary
    banner(f"SUMMARY: {label}")
    log("  A3 wall status:")
    for w in agg["wall"]:
        log(f"    C{w['community']} {w['a3_gene']} ({w['direction']}): "
            f"{w['verdict']}")
        log(f"      edges {w['total_edges']}: "
            f"act-concordant {w['concordant_activator_edges']}, "
            f"wall {w['wall_edges']}, "
            f"rep-concordant {w['concordant_repressor_edges']}")
    log("\n  Largest concordant subnetworks:")
    for mode in ["activator", "repressor"]:
        ms = [s for s in agg["subnets"] if s["mode"] == mode]
        ms.sort(key=lambda x: -x["size"])
        if ms:
            top = ms[0]
            log(f"    {mode}: C{top['community']} {top['size']} nodes, "
                f"{top['n_harris']} Harris, "
                f"A3={top['a3_members'] or 'none'}, "
                f"gap_to_A3={top['gap_hops_to_a3']}")
    n_h = len(set((h["community"], h["harris_gene"])
                  for h in agg["harris_index"]))
    log(f"\n  Harris interactors placed in a concordant subnetwork: {n_h}")


def main():
    banner("A3 NETWORK STRUCTURE DIAGNOSTIC")
    banner("Multi-community: edge composition, subnetwork enumeration, wall")
    harris = load_harris()
    for net_config in NETWORKS:
        log(f"\nLoading {net_config['name']}...")
        analyze_network(load_network(net_config), harris)
    banner("DIAGNOSTIC COMPLETE")
    log(f"\nOutput: {OUTPUT_DIR}")
    for f in sorted(os.listdir(OUTPUT_DIR)):
        log(f"  {f}")


if __name__ == "__main__":
    main()
