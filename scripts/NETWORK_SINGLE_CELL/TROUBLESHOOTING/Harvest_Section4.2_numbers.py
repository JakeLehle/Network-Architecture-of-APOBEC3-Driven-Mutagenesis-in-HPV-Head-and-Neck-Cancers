#!/usr/bin/env python3
"""
Harvest_Section4.2_numbers.py
=============================================================================
Reads the CURRENT on-disk Figure 4 network outputs (post CNV-HIGH reselection,
job 129065) and prints every bracketed number in the Section 4.2 results and
the network-analysis methods, mapped to the exact sentence it fills.

WHY THIS EXISTS
  The older Diagnostic_section4_numbers.py walks the whole FIG_4 tree and still
  matches legacy NETWORK_A/B/C names; its saved output.md reports the PRE-
  reselection CNV-HIGH run (DIFF 0.40, 5,174 genes, 16 groups, 116 Harris, A3A
  and A3B in the SAME community). This harvester instead targets the two
  retained networks by fixed path and recomputes the group-level wall fraction
  directly from the graph, so nothing depends on a stale log.

SCOPE
  Only SBS2-HIGH vs NORMAL and CNV-HIGH vs NORMAL are read. SBS2-HIGH vs CNV-HIGH
  is intentionally excluded (spun off to the tumor-vs-tumor paper).

SOURCES (per network, under data/FIG_4/NETWORK_<name>/)
  04_communities/SC_selected_parameters.txt   threshold, resolution, modularity,
                                               ARI, NMI, evenness, N_GENES,
                                               N_COMMUNITIES, N_SATELLITE
  04_communities/SC_best_partition.csv         gene -> community (+ a3_alias,
                                               is_satellite); A3 placement and
                                               Harris recovery
  04_communities/SC_G_comm.gpickle             thresholded DIFF graph w/ signed
                                               edge weights -> A3 degrees and
                                               per-group wall fraction
  04_communities/SC_threshold_sweep.csv        fragmentation margin at selection
  02_differential_expression/SC_diffexpr_stats.csv
                                               per-gene fdr / log2FC (cols:
                                               gene, p_value, fdr, log2FC, ...)
  ../DIAGNOSTIC_CONCORDANCE/<name>_A3_wall_status.tsv
  ../DIAGNOSTIC_CONCORDANCE/<name>_concordant_subnetworks.tsv

Run in the NETWORK conda env (needs pandas, networkx):
  conda run -n NETWORK python Harvest_Section4.2_numbers.py

Author: Jake Lehle / Texas Biomedical Research Institute
"""

import os
import sys
import pickle
import pandas as pd
import networkx as nx

# =============================================================================
# CONFIG (edit only if the data root moves)
# =============================================================================
FIG4 = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_4"
CONCORD = os.path.join(FIG4, "DIAGNOSTIC_CONCORDANCE")
HARRIS_PATH = os.path.join(FIG4, "00_input", "Harris_A3_interactors.txt")

A3A, A3B = "APOBEC3A", "APOBEC3B"

# Retained networks only. expected_* are the current draft values; FLAG fires
# on mismatch so any drift is loud.
NETWORKS = [
    {"key": "SBS2_VS_NORMAL", "dir": "NETWORK_SBS2_VS_NORMAL",
     "label": "SBS2-HIGH vs NORMAL",
     "exp_genes": 2948, "exp_groups": 23, "exp_thr": 0.40,
     "exp_recovery": 54, "exp_deg": (55, 16)},
    {"key": "CNV_VS_NORMAL", "dir": "NETWORK_CNV_VS_NORMAL",
     "label": "CNV-HIGH vs NORMAL",
     "exp_genes": 4886, "exp_groups": 38, "exp_thr": 0.45,
     "exp_recovery": 109, "exp_deg": (None, None)},
]

SEP = "=" * 74


def flag(found, expected, tol=0):
    """Return a short agreement marker for found-vs-expected."""
    if expected is None:
        return "(no draft value)"
    try:
        ok = abs(float(found) - float(expected)) <= tol
    except (TypeError, ValueError):
        ok = str(found) == str(expected)
    return "OK matches draft" if ok else f"FLAG draft says {expected}"


# =============================================================================
# LOADERS
# =============================================================================
def load_params(net_dir):
    """Parse SC_selected_parameters.txt (KEY=VALUE per line) into a dict."""
    path = os.path.join(net_dir, "04_communities", "SC_selected_parameters.txt")
    if not os.path.exists(path):
        return None, path
    d = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if "=" in line and not line.startswith("#"):
                k, v = line.split("=", 1)
                d[k.strip()] = v.strip()
    return d, path


def load_graph(net_dir):
    """Load SC_G_comm.gpickle (pickle.dump under the hood for any nx version)."""
    path = os.path.join(net_dir, "04_communities", "SC_G_comm.gpickle")
    if not os.path.exists(path):
        return None, path
    try:
        with open(path, "rb") as f:
            G = pickle.load(f)
    except Exception:
        G = nx.read_gpickle(path)  # fallback for older nx
    return G, path


def load_partition(net_dir):
    path = os.path.join(net_dir, "04_communities", "SC_best_partition.csv")
    if not os.path.exists(path):
        return None, path
    return pd.read_csv(path), path


def load_de(net_dir):
    path = os.path.join(net_dir, "02_differential_expression",
                        "SC_diffexpr_stats.csv")
    if not os.path.exists(path):
        return None, path
    return pd.read_csv(path), path


def load_harris():
    if not os.path.exists(HARRIS_PATH):
        return None
    df = pd.read_csv(HARRIS_PATH, sep="\t")
    col = "gene_symbol" if "gene_symbol" in df.columns else df.columns[0]
    return set(df[col].astype(str))


def load_tsv(name, suffix):
    path = os.path.join(CONCORD, f"{name}{suffix}")
    if not os.path.exists(path):
        return None, path
    return pd.read_csv(path, sep="\t"), path


# =============================================================================
# ANALYSIS
# =============================================================================
def a3_degrees(G):
    """Full-graph degree of A3A and A3B (the 'degree = N' values in results)."""
    out = {}
    for a3 in (A3A, A3B):
        out[a3] = G.degree(a3) if (G is not None and a3 in G) else "ABSENT"
    return out


def a3_communities(part):
    """Community id of A3A / A3B from the partition (same group vs split)."""
    out = {}
    if part is None:
        return out
    for a3 in (A3A, A3B):
        row = part[part["gene"] == a3]
        out[a3] = int(row["community"].iloc[0]) if len(row) else "ABSENT"
    return out


def harris_recovery(part, harris):
    if part is None or harris is None:
        return None
    return len(set(part["gene"].astype(str)) & harris)


def group_wall_fraction(G, part, de):
    """For each community holding A3A or A3B, recompute the wall fraction
    directly from the graph: an edge is a wall edge when its DIFF weight is
    negative AND both endpoints are up in tumor (log2FC > 0). Returns the
    84-94% range quoted in results para 7."""
    if G is None or part is None or de is None:
        return []
    log2fc = dict(zip(de["gene"].astype(str), de["log2FC"]))
    comm_of = dict(zip(part["gene"].astype(str), part["community"]))
    a3_comms = sorted({comm_of.get(a3) for a3 in (A3A, A3B)
                       if comm_of.get(a3) is not None})
    results = []
    for c in a3_comms:
        nodes = [g for g, cc in comm_of.items() if cc == c and g in G]
        sub = G.subgraph(nodes)
        n_edges = sub.number_of_edges()
        wall = pos = rep = 0
        for u, v, d in sub.edges(data=True):
            w = d.get("weight", 0.0)
            uu, vv = log2fc.get(u, 0.0) > 0, log2fc.get(v, 0.0) > 0
            if w < 0 and uu and vv:
                wall += 1
            elif w < 0 and (not uu) and (not vv):
                rep += 1
            elif w > 0:
                pos += 1
        a3_here = [a for a in (A3A, A3B) if comm_of.get(a) == c]
        results.append({
            "community": c, "a3_genes": ",".join(a3_here) or "(none)",
            "nodes": len(nodes), "edges": n_edges,
            "wall_pct": round(100 * wall / n_edges, 1) if n_edges else 0.0,
            "pos_pct": round(100 * pos / n_edges, 1) if n_edges else 0.0,
            "rep_pct": round(100 * rep / n_edges, 1) if n_edges else 0.0,
        })
    return results


def fragmentation_margin(net_dir, selected_thr):
    """Print the rows of SC_threshold_sweep.csv around the selected threshold
    so the fragmentation margin (delta-component) can be read off."""
    path = os.path.join(net_dir, "04_communities", "SC_threshold_sweep.csv")
    if not os.path.exists(path):
        print(f"    [missing] {path}")
        return
    df = pd.read_csv(path)
    print(f"    columns: {list(df.columns)}")
    thr_col = next((c for c in df.columns if "thr" in c.lower()), df.columns[0])
    try:
        df["_d"] = (df[thr_col].astype(float) - float(selected_thr)).abs()
        window = df.sort_values("_d").head(5).drop(columns="_d")
        print(window.to_string(index=False))
    except Exception:
        print(df.head(10).to_string(index=False))


# =============================================================================
# MAIN
# =============================================================================
def main():
    harris = load_harris()
    print(SEP)
    print("SECTION 4.2 NUMBER HARVEST  (retained networks only)")
    print(f"Harris list: {'loaded ' + str(len(harris)) if harris else 'MISSING'} interactors")
    print(SEP)

    for net in NETWORKS:
        net_dir = os.path.join(FIG4, net["dir"])
        print(f"\n{SEP}\nNETWORK: {net['label']}   [{net['dir']}]\n{SEP}")
        if not os.path.isdir(net_dir):
            print(f"  [MISSING DIR] {net_dir}")
            continue

        params, ppath = load_params(net_dir)
        part, _ = load_partition(net_dir)
        G, gpath = load_graph(net_dir)
        de, _ = load_de(net_dir)

        # ---- A) Network build (methods + results para 3) ----------------
        print("\n[A] NETWORK BUILD  -> methods (threshold/Leiden) + results para 3")
        if params:
            ng = int(params.get("N_GENES", -1))
            nc = int(params.get("N_COMMUNITIES", -1))
            thr = float(params.get("DIFF_THRESHOLD", -1))
            res = params.get("LEIDEN_RESOLUTION", "?")
            mod = params.get("MODULARITY", "?")
            ari = params.get("ARI", "?")
            print(f"    N_GENES        = {ng:>6}   {flag(ng, net['exp_genes'])}   (results: '{net['exp_genes']} genes')")
            print(f"    N_COMMUNITIES  = {nc:>6}   {flag(nc, net['exp_groups'])}   (results: '{net['exp_groups']} gene groups')")
            print(f"    DIFF_THRESHOLD = {thr:>6}   {flag(thr, net['exp_thr'])}   (methods threshold)")
            print(f"    LEIDEN_RES     = {res:>6}   {flag(res, 0.70)}   (methods resolution [0.70])")
            print(f"    MODULARITY     = {mod}   ARI = {ari}   NMI = {params.get('NMI','?')}   EVENNESS = {params.get('EVENNESS','?')}")
            print(f"    N_SATELLITE    = {params.get('N_SATELLITE','?')}   N_COMPONENTS = {params.get('N_COMPONENTS','?')}")
            print(f"    source: {ppath}")
        else:
            print(f"    [MISSING] {ppath}")

        print("\n    Fragmentation margin (methods 'fragmentation margins'):")
        if params:
            fragmentation_margin(net_dir, params.get("DIFF_THRESHOLD", net["exp_thr"]))

        # ---- B) A3 adjusted-p + log2FC (methods DE para) ----------------
        print("\n[B] A3 DIFFERENTIAL EXPRESSION  -> methods DE para (four adj p)")
        if de is not None:
            for a3, lab in ((A3A, "A3A"), (A3B, "A3B")):
                row = de[de["gene"] == a3]
                if len(row):
                    r = row.iloc[0]
                    print(f"    {lab:3} ({a3}):  adj p (fdr) = {r['fdr']:.3e}   "
                          f"log2FC = {float(r['log2FC']):.4f}   raw p = {r['p_value']:.3e}")
                else:
                    print(f"    {lab}: ABSENT from SC_diffexpr_stats.csv")
        else:
            print("    [MISSING] SC_diffexpr_stats.csv")

        # ---- C) Recovery, A3 placement, A3 degree (results para 3) ------
        print("\n[C] RECOVERY / A3 PLACEMENT / A3 DEGREE  -> results para 3")
        rec = harris_recovery(part, harris)
        print(f"    Harris recovery = {rec} / 174   {flag(rec, net['exp_recovery'])}")
        comms = a3_communities(part)
        print(f"    A3A community = {comms.get(A3A)}   A3B community = {comms.get(A3B)}   "
              f"({'SAME group' if comms.get(A3A) == comms.get(A3B) else 'SPLIT into different groups'})")
        deg = a3_degrees(G)
        eda, edb = net["exp_deg"]
        print(f"    A3A degree = {deg.get(A3A)}   {flag(deg.get(A3A), eda)}")
        print(f"    A3B degree = {deg.get(A3B)}   {flag(deg.get(A3B), edb)}")
        print(f"    graph source: {gpath}")

        # ---- D) Per-group wall fraction (results para 7, 84-94%) --------
        print("\n[D] PER-GROUP WALL FRACTION (recomputed from graph)  -> results para 7")
        gw = group_wall_fraction(G, part, de)
        if gw:
            for r in gw:
                print(f"    C{r['community']} [{r['a3_genes']}]: {r['nodes']} nodes, "
                      f"{r['edges']} edges, wall {r['wall_pct']}%  "
                      f"(pos {r['pos_pct']}%, repressor {r['rep_pct']}%)")
            print(f"    -> quoted range = {min(r['wall_pct'] for r in gw)}"
                  f" to {max(r['wall_pct'] for r in gw)} percent")
        else:
            print("    [could not compute - missing graph/partition/DE]")

        # ---- E) Per-A3 edge counts / wall (results para 7) --------------
        print("\n[E] PER-A3 EDGE COUNTS  -> results para 7 ('all NN ... were negative')")
        wall, wpath = load_tsv(net["key"], "_A3_wall_status.tsv")
        if wall is not None:
            for _, r in wall.iterrows():
                pos = int(r.get("concordant_activator_edges", 0)) + int(r.get("discordant_pos_edges", 0))
                neg = int(r.get("wall_edges", 0)) + int(r.get("concordant_repressor_edges", 0))
                print(f"    C{r['community']} {r['a3_gene']}: total {int(r['total_edges'])} edges, "
                      f"{pos} positive, {neg} negative "
                      f"({int(r.get('wall_edges',0))} wall + "
                      f"{int(r.get('concordant_repressor_edges',0))} repressor)")
            print(f"    source: {wpath}")
        else:
            print(f"    [MISSING] {wpath}")

        # ---- F) Concordant subnetworks (results para 5 + 6) -------------
        print("\n[F] CONCORDANT SUBNETWORKS  -> results para 5 (activator) + 6 (70-node inhibitor)")
        sub, spath = load_tsv(net["key"], "_concordant_subnetworks.tsv")
        if sub is not None:
            for mode in ("activator", "repressor"):
                ms = sub[sub["mode"] == mode].sort_values("size", ascending=False)
                if len(ms):
                    top = ms.iloc[0]
                    extra = {k: top[k] for k in top.index
                             if k not in ("mode", "size", "community")}
                    print(f"    largest {mode}: C{top['community']}, size {int(top['size'])}")
                    print(f"        {extra}")
                else:
                    print(f"    largest {mode}: none")
            print(f"    source: {spath}")
        else:
            print(f"    [MISSING] {spath}")

    print(f"\n{SEP}\nHARVEST COMPLETE\n{SEP}")


if __name__ == "__main__":
    main()
