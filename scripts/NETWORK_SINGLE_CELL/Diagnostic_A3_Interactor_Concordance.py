#!/usr/bin/env python3
"""
Diagnostic_A3_Network_Structure_and_Numbers.py > Still called Diagnostic_A3_Interactor_Concordance.py but this is the updated version
=============================================================================
Merge of two diagnostics into a single pass so the structural map and every
Section 4.2 / methods number come out of one run with full cross-context:

  (1) BUILD PROVENANCE HARVEST  (folded in from Harvest_Section4.2_numbers.py)
      Per network, reads the on-disk build outputs and prints the values that
      fill the bracketed numbers in the results and methods, each tagged with
      the sentence it fills and flagged against the confirmed draft value:
        - 04_communities/SC_selected_parameters.txt   thr, resolution,
          modularity, ARI, NMI, evenness, N_GENES, N_COMMUNITIES, N_SATELLITE
        - 02_differential_expression/SC_diffexpr_stats.csv  A3A/A3B fdr + log2FC
        - 04_communities/SC_threshold_sweep.csv        fragmentation margin
        - 04_communities/SC_G_comm.gpickle             A3A/A3B full-graph degree
        - per-group wall fraction recomputed from the graph

  (2) A3 NETWORK STRUCTURE  (unchanged core of Diagnostic_A3_Network_Structure.py)
      Per A3-containing community: edge-composition cross-tab + per-A3 edge
      profile and wall verdict (PART A), full activator/repressor concordant
      subnetwork enumeration with Harris/A3 annotation and hop-gap (PART B),
      and A3 boundary / decoupling-point classification (PART C). Writes the
      same five TSVs to DIAGNOSTIC_CONCORDANCE/ as before.

  (3) CONSOLIDATED BRACKET REPORT
      After both, a per-network report listing each Section 4.2 number next to
      the confirmed value and a flag, so manuscript drift is loud.

SCOPE: SBS2-HIGH vs NORMAL and CNV-HIGH vs NORMAL only (SBS2-vs-CNV is spun off).

Confirmed values baked into the flags are the on-disk truth as of job 129065
(verified 2026-06-22). KNOWN DRAFT FIX: CNV-HIGH resolution is 0.80, the prose
still says 0.70.

Usage:
    conda run -n NETWORK python Diagnostic_A3_Network_Structure_and_Numbers.py

Author: Jake Lehle / Texas Biomedical Research Institute
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
     "dir": os.path.join(FIG4_ROOT, "NETWORK_SBS2_VS_NORMAL"),
     # confirmed on-disk values for the flag (NOT the draft, which is wrong on res)
     "exp": {"genes": 2948, "groups": 23, "thr": 0.40, "res": 0.70,
             "recovery": 54, "deg_A3A": 55, "deg_A3B": 16}},
    {"name": "CNV_VS_NORMAL", "label": "CNV-HIGH vs NORMAL",
     "dir": os.path.join(FIG4_ROOT, "NETWORK_CNV_VS_NORMAL"),
     "exp": {"genes": 4886, "groups": 38, "thr": 0.45, "res": 0.80,
             "recovery": 109, "deg_A3A": 15, "deg_A3B": 136}},
]

# Manuscript prose values that differ from on-disk truth (printed as fixes).
DRAFT_FIXES = {
    "CNV_VS_NORMAL": [("Leiden resolution", "prose says 0.70", "disk = 0.80")],
}

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

def flag(found, expected, tol=0):
    if expected is None:
        return "(no value)"
    try:
        ok = abs(float(found) - float(expected)) <= tol
    except (TypeError, ValueError):
        ok = str(found) == str(expected)
    return "OK" if ok else f"FLAG (confirmed {expected})"


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
    result = {"name": net_config["name"], "label": net_config["label"],
              "dir": net_dir, "exp": net_config["exp"]}

    part_df = pd.read_csv(
        os.path.join(net_dir, "04_communities/SC_best_partition.csv"))
    result["gene_to_comm"] = dict(zip(part_df["gene"], part_df["community"]))
    result["partition"] = part_df

    with open(os.path.join(net_dir,
              "04_communities/SC_G_comm.gpickle"), "rb") as f:
        result["G_comm"] = pickle.load(f)

    de_df = pd.read_csv(
        os.path.join(net_dir,
                     "02_differential_expression/SC_diffexpr_stats.csv"))
    result["de_df"] = de_df
    result["de_log2fc"] = dict(zip(de_df["gene"], de_df["log2FC"]))
    return result


# =============================================================================
# (1) BUILD PROVENANCE HARVEST  (folded in from the standalone harvester)
# =============================================================================

def read_params(net_dir):
    path = os.path.join(net_dir, "04_communities", "SC_selected_parameters.txt")
    if not os.path.exists(path):
        return None
    d = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if "=" in line and not line.startswith("#"):
                k, v = line.split("=", 1)
                d[k.strip()] = v.strip()
    return d


def fragmentation_margin(net_dir, selected_thr):
    """Largest positive delta-component step at an A3-valid upper threshold.
    Mirrors Step03's max-fragmentation-rate selection."""
    path = os.path.join(net_dir, "04_communities", "SC_threshold_sweep.csv")
    if not os.path.exists(path):
        return None, None
    df = pd.read_csv(path).sort_values("threshold").reset_index(drop=True)
    df["d_comp"] = df["components"].diff()
    valid = df[(df.get("deg_A3A", 1) >= 1) & (df.get("deg_A3B", 1) >= 1)]
    margin_row = None
    if "d_comp" in valid:
        pos = valid[valid["d_comp"] > 0]
        if len(pos):
            margin_row = pos.loc[pos["d_comp"].idxmax()]
    return df, margin_row


def harvest_build(net_data):
    """Print the build numbers that fill the methods + results-para-3 brackets."""
    net_dir, exp = net_data["dir"], net_data["exp"]
    G = net_data["G_comm"]
    de = net_data["de_df"]

    banner(f"[BUILD HARVEST] {net_data['label']}", "-")
    params = read_params(net_dir)
    harvested = {}

    if params:
        ng = int(params.get("N_GENES", -1))
        nc = int(params.get("N_COMMUNITIES", -1))
        thr = float(params.get("DIFF_THRESHOLD", -1))
        res = float(params.get("LEIDEN_RESOLUTION", -1))
        harvested.update(genes=ng, groups=nc, thr=thr, res=res,
                         mod=params.get("MODULARITY"), ari=params.get("ARI"),
                         nmi=params.get("NMI"), evenness=params.get("EVENNESS"))
        log(f"  N_GENES        = {ng}   {flag(ng, exp['genes'])}   -> results 'N genes'")
        log(f"  N_COMMUNITIES  = {nc}   {flag(nc, exp['groups'])}   -> results 'N gene groups'")
        log(f"  DIFF_THRESHOLD = {thr}   {flag(thr, exp['thr'])}   -> methods threshold")
        log(f"  LEIDEN_RES     = {res}   {flag(res, exp['res'])}   -> methods resolution")
        log(f"  MODULARITY     = {params.get('MODULARITY')}   ARI = {params.get('ARI')}   "
            f"NMI = {params.get('NMI')}   EVENNESS = {params.get('EVENNESS')}   -> methods [NN]")
    else:
        log("  [MISSING] SC_selected_parameters.txt")

    # A3 differential expression (four adj-p)
    for a3, lab in ((("APOBEC3A"), "A3A"), (("APOBEC3B"), "A3B")):
        row = de[de["gene"] == a3]
        if len(row):
            r = row.iloc[0]
            harvested[f"{lab}_fdr"] = r["fdr"]
            harvested[f"{lab}_log2fc"] = float(r["log2FC"])
            log(f"  {lab} adj p (fdr) = {r['fdr']:.3e}   log2FC = {float(r['log2FC']):.4f}"
                f"   -> methods DE para")
        else:
            log(f"  {lab}: ABSENT from SC_diffexpr_stats.csv")

    # A3 full-graph degree
    for a3, lab in (("APOBEC3A", "A3A"), ("APOBEC3B", "A3B")):
        d = G.degree(a3) if a3 in G else "ABSENT"
        harvested[f"{lab}_deg"] = d
        log(f"  {lab} degree = {d}   {flag(d, exp.get('deg_'+lab))}   -> results para 3")

    # fragmentation margin
    if params:
        _, mrow = fragmentation_margin(net_dir, harvested.get("thr"))
        if mrow is not None:
            log(f"  fragmentation margin: +{int(mrow['d_comp'])} components into "
                f"threshold {mrow['threshold']}   -> methods 'fragmentation margins'")
            harvested["frag_margin"] = int(mrow["d_comp"])

    return harvested


def group_wall_fraction(G, gene_to_comm, de_log2fc):
    """Per A3-community wall fraction recomputed from the graph (results para 7)."""
    a3_comms = sorted({gene_to_comm.get(a) for a in A3_TARGETS
                       if gene_to_comm.get(a) is not None})
    out = []
    for c in a3_comms:
        nodes = [g for g, cc in gene_to_comm.items() if cc == c and g in G]
        sub = G.subgraph(nodes)
        ne = sub.number_of_edges()
        wall = pos = rep = 0
        for u, v, d in sub.edges(data=True):
            w = d.get("weight", 0.0)
            uu, vv = de_log2fc.get(u, 0) > 0, de_log2fc.get(v, 0) > 0
            if w < 0 and uu and vv:
                wall += 1
            elif w < 0 and not uu and not vv:
                rep += 1
            elif w > 0:
                pos += 1
        a3_here = [A3_ALIAS[a] for a in A3_TARGETS if gene_to_comm.get(a) == c]
        out.append({"community": c, "a3": ",".join(a3_here),
                    "nodes": len(nodes), "edges": ne,
                    "wall_pct": round(100 * wall / ne, 1) if ne else 0.0})
    return out


# =============================================================================
# (2) A3 NETWORK STRUCTURE  -- directions, concordance subgraphs, edge class
# =============================================================================

def node_directions(G_sub, de_log2fc):
    up_tumor, up_normal = set(), set()
    for n in G_sub.nodes():
        (up_tumor if de_log2fc.get(n, 0) > 0 else up_normal).add(n)
    return up_tumor, up_normal


def build_concordance_subgraphs(G_sub, up_tumor, up_normal):
    a3 = set(n for n in G_sub.nodes() if n in A3_SYMBOLS)
    act_nodes = up_tumor | a3
    rep_nodes = up_normal | a3
    G_act, G_rep = nx.Graph(), nx.Graph()
    for u, v, d in G_sub.edges(data=True):
        w = d.get("weight", 0)
        if w > 0 and u in act_nodes and v in act_nodes:
            G_act.add_edge(u, v, weight=w)
        elif w < 0 and u in rep_nodes and v in rep_nodes:
            G_rep.add_edge(u, v, weight=w)
    return G_act, G_rep


def edge_class(w, du, dv):
    if w > 0 and du == "up_tumor" and dv == "up_tumor":
        return "activator_concordant"
    if w < 0 and du == "up_normal" and dv == "up_normal":
        return "repressor_concordant"
    if w < 0 and du == "up_tumor" and dv == "up_tumor":
        return "wall"
    if w > 0 and du == "up_normal" and dv == "up_normal":
        return "anti_normal"
    return "cross_pos" if w > 0 else "cross_neg"


# ---- PART A -----------------------------------------------------------------
def part_a_edge_profile(G_sub, up_tumor, up_normal, de_log2fc,
                        harris, a3_in_comm, comm, name):
    dirof = lambda g: "up_tumor" if g in up_tumor else "up_normal"

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
    conc = counts.get("activator_concordant", 0) + counts.get("repressor_concordant", 0)
    log(f"    -> concordant fraction: {100*conc/total:.1f}%  "
        f"(low fraction = wall-dominated community)")

    rows, wall_rows = [], []
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
            log("    Harris (positive edge): " +
                ", ".join(f"{g}(w={w},fc={fc},{d})" for g, w, fc, d in harris_pos))
        if harris_neg:
            log("    Harris (negative edge): " +
                ", ".join(f"{g}(w={w},fc={fc},{d})" for g, w, fc, d in harris_neg))

        if a3dir == "up_tumor":
            if pos_up == 0:
                verdict = f"WALLED: zero gained co-expression with the up-program; {neg_up} wall edges"
            elif pos_up < neg_up:
                verdict = f"MOSTLY WALLED: {pos_up} concordant-activator vs {neg_up} wall edges"
            else:
                verdict = f"PARTIALLY COUPLED: {pos_up} concordant-activator edges"
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


# ---- PART B -----------------------------------------------------------------
def part_b_subnetworks(G_sub, G_act, G_rep, de_log2fc, harris,
                       a3_in_comm, comm, name):
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

    subnet_rows, harris_index = [], []
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
            gap, _ = gap_to_a3(comp)
            tag = ""
            if a3_members:
                tag = f"  <-- contains {', '.join(a3_members)}"
            elif gap is not None:
                tag = f"  (nearest A3 {gap} hop{'s' if gap != 1 else ''} away)"
            if i < TOP_SUBNETS_PRINTED or a3_members or h_members:
                log(f"    [{mode[:3]}-{i}] {len(comp)} nodes, "
                    f"{sub.number_of_edges()} edges, Harris={len(h_members)}{tag}")
                log(f"         hubs: {hub_str}")
                if h_members:
                    log(f"         Harris: {', '.join(h_members)}")
            subnet_rows.append({
                "community": comm, "mode": mode, "subnet_id": f"{mode[:3]}-{i}",
                "size": len(comp), "n_edges": sub.number_of_edges(),
                "n_harris": len(h_members), "harris_members": ", ".join(h_members),
                "a3_members": ", ".join(a3_members), "top_hubs": hub_str,
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


# ---- PART C -----------------------------------------------------------------
def part_c_boundary(G_sub, G_act, G_rep, up_tumor, up_normal,
                    de_log2fc, harris, a3_in_comm, comm, name):
    dirof = lambda g: "up_tumor" if g in up_tumor else "up_normal"
    rows = []
    for a3g in a3_in_comm:
        alias = A3_ALIAS[a3g]
        section(f"{alias} boundary (C{comm})")
        nbrs = [n for n in G_sub.neighbors(a3g) if n not in A3_SYMBOLS]
        n_bridge = n_decouple = 0
        for nb in nbrs:
            w = G_sub[a3g][nb].get("weight", 0)
            nbdir = dirof(nb)
            edge_conc = ((nbdir == "up_tumor" and w > 0) or
                         (nbdir == "up_normal" and w < 0))
            G_conc = G_act if nbdir == "up_tumor" else G_rep
            beyond = 0
            if nb in G_conc:
                beyond = sum(1 for x in G_conc.neighbors(nb)
                             if x != a3g and x not in A3_SYMBOLS)
            if beyond == 0:
                continue
            if edge_conc:
                n_bridge += 1; kind = "concordant bridge"
            else:
                n_decouple += 1; kind = "DECOUPLING POINT (wall edge into a real chain)"
            rows.append({
                "community": comm, "a3_gene": alias, "boundary_gene": nb,
                "edge_weight": round(w, 4), "neighbor_de": nbdir,
                "neighbor_log2FC": round(de_log2fc.get(nb, 0), 4),
                "edge_to_a3_concordant": edge_conc,
                "concordant_edges_beyond": beyond,
                "is_harris": (nb in harris and nb not in A3_SYMBOLS),
                "kind": kind,
            })
        log(f"    neighbors anchoring concordant chains: {n_bridge + n_decouple} "
            f"({n_bridge} concordant bridges, {n_decouple} decoupling points)")
        chain_rows = [r for r in rows
                      if r["a3_gene"] == alias and r["community"] == comm]
        for r in sorted(chain_rows,
                        key=lambda x: -x["concordant_edges_beyond"])[:15]:
            hm = " [Harris]" if r["is_harris"] else ""
            log(f"      {r['boundary_gene']}{hm}: edge w={r['edge_weight']}, "
                f"{r['neighbor_de']}, chains beyond={r['concordant_edges_beyond']} "
                f"-> {r['kind']}")
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
    G_sub.remove_nodes_from([n for n in G_sub.nodes() if G_sub.degree(n) == 0])

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


def consolidated_report(net_data, harvested, agg):
    """Per-network Section 4.2 bracket fill, drawing build numbers from the
    harvest and structural numbers from the enumeration."""
    name, label, exp = net_data["name"], net_data["label"], net_data["exp"]
    gene_to_comm, G, de_log2fc = (net_data["gene_to_comm"], net_data["G_comm"],
                                  net_data["de_log2fc"])
    banner(f"[SECTION 4.2 BRACKET REPORT] {label}", "#")

    log("  Methods / results para 3:")
    log(f"    network size   : {harvested.get('genes')} genes, "
        f"{harvested.get('groups')} gene groups")
    log(f"    DIFF threshold : {harvested.get('thr')}  (margin +{harvested.get('frag_margin','?')} components)")
    log(f"    Leiden res     : {harvested.get('res')}   modularity {harvested.get('mod')}   ARI {harvested.get('ari')}")
    recovery = len(set(net_data["partition"]["gene"].astype(str)) & load_harris.cache)
    log(f"    recovery       : {recovery} / 174   {flag(recovery, exp['recovery'])}")

    log("  Methods DE para (adj p, log2FC):")
    log(f"    A3A fdr {harvested.get('A3A_fdr'):.3e}  log2FC {harvested.get('A3A_log2fc'):.4f}")
    log(f"    A3B fdr {harvested.get('A3B_fdr'):.3e}  log2FC {harvested.get('A3B_log2fc'):.4f}")

    log("  Results para 3 (A3 degree, placement):")
    a3a_c, a3b_c = gene_to_comm.get("APOBEC3A"), gene_to_comm.get("APOBEC3B")
    log(f"    A3A degree {harvested.get('A3A_deg')}  (C{a3a_c}) | "
        f"A3B degree {harvested.get('A3B_deg')}  (C{a3b_c}) | "
        f"{'SAME group' if a3a_c == a3b_c else 'SPLIT'}")

    log("  Results para 7 (wall):")
    for w in agg["wall"]:
        pos = w["concordant_activator_edges"] + w["discordant_pos_edges"]
        neg = w["wall_edges"] + w["concordant_repressor_edges"]
        log(f"    C{w['community']} {w['a3_gene']}: {w['total_edges']} edges, "
            f"{pos} positive, {neg} negative "
            f"({w['wall_edges']} wall + {w['concordant_repressor_edges']} repressor)")
    gw = group_wall_fraction(G, gene_to_comm, de_log2fc)
    log("    group wall fraction: " +
        ", ".join(f"C{r['community']}[{r['a3']}] {r['wall_pct']}%" for r in gw)
        + f"   -> range {min(r['wall_pct'] for r in gw)}-{max(r['wall_pct'] for r in gw)}%")

    log("  Results para 5/6 (concordant chains):")
    for mode in ("activator", "repressor"):
        ms = sorted((s for s in agg["subnets"] if s["mode"] == mode),
                    key=lambda x: -x["size"])
        if ms:
            t = ms[0]
            log(f"    largest {mode}: C{t['community']} size {t['size']}, "
                f"Harris [{t['harris_members'] or 'none'}], A3 [{t['a3_members'] or 'none'}]")
            log(f"        hubs: {t['top_hubs']}")

    for fixlabel, says, disk in DRAFT_FIXES.get(name, []):
        log(f"  *** DRAFT FIX: {fixlabel} -- {says}, {disk} ***")


def analyze_network(net_data, harris):
    name, label = net_data["name"], net_data["label"]
    gene_to_comm = net_data["gene_to_comm"]
    banner(f"NETWORK: {label}")

    # (1) build provenance harvest
    harvested = harvest_build(net_data)

    placement = {g: gene_to_comm[g] for g in A3_TARGETS if g in gene_to_comm}
    if not placement:
        log("  [SKIP] A3A/A3B not in partition"); return
    log("\n  A3 placement: " +
        ", ".join(f"{A3_ALIAS[g]}=C{c}" for g, c in placement.items()))
    a3_comms = sorted(set(placement.values()))
    if len(a3_comms) > 1:
        log("  NOTE: A3A and A3B split across communities; analyzing each.")

    # (2) structural enumeration
    agg = {k: [] for k in ["edges", "wall", "subnets", "harris_index", "boundary"]}
    for comm in a3_comms:
        out = analyze_a3_community(net_data, comm, harris)
        for k in agg:
            agg[k] += out[k]

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

    # (3) consolidated bracket report
    consolidated_report(net_data, harvested, agg)


def main():
    banner("A3 NETWORK STRUCTURE + SECTION 4.2 NUMBERS")
    harris = load_harris()
    load_harris.cache = harris  # for recovery count in the report
    for net_config in NETWORKS:
        log(f"\nLoading {net_config['name']}...")
        analyze_network(load_network(net_config), harris)
    banner("DIAGNOSTIC COMPLETE")
    log(f"\nOutput TSVs: {OUTPUT_DIR}")
    for f in sorted(os.listdir(OUTPUT_DIR)):
        log(f"  {f}")


if __name__ == "__main__":
    main()
