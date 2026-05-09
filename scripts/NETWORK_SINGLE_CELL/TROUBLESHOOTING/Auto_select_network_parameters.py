#!/usr/bin/env python3
"""
Auto_Select_Network_Parameters.py
====================================

Automatically selects optimal DIFF threshold and Leiden resolution
for the SC differential co-expression network pipeline.

DIFF THRESHOLD SELECTION:
  Criterion: Peak connected components with LCC < 300 genes.
  Logic:
    1. Find the threshold where # connected components is maximized
    2. If LCC at that threshold is < 300, use it
    3. If LCC > 300, step up to the next higher threshold where LCC < 300
    4. Verify the selected threshold still has reasonable component count
       (within 80% of the peak)

  Rationale: The component peak marks where the network fragments along
  natural biological boundaries. The LCC constraint ensures the largest
  community is pathway-scale (< 300 genes) and biologically interpretable.

LEIDEN RESOLUTION SELECTION:
  Criterion: Largest modularity jump with ARI stability floor.
  Logic:
    1. Compute delta-Q between consecutive resolutions
    2. Find the resolution with the largest single-step Q increase
       (the "structural jump" where Leiden discovers real communities)
    3. Select that resolution if ARI > 0.65 (stability floor)
    4. If ARI fails, step back to the previous resolution
    5. Report efficiency (Q gained per ARI lost) for the next step
       to confirm diminishing returns

  Rationale: The trivial solution (2 communities at low resolution) never
  produces a large Q jump. The structural jump identifies where the
  algorithm found genuine community boundaries. The stability floor
  (matching TCGA bulk tolerance of ARI ~0.63) ensures reproducibility.

USAGE:
  # Select DIFF threshold from sweep data:
  python Auto_Select_Network_Parameters.py --mode threshold \
      --sweep_csv path/to/sweep_detailed.csv

  # Select Leiden resolution from resolution sweep:
  python Auto_Select_Network_Parameters.py --mode resolution \
      --sweep_csv path/to/SC_resolution_sweep.csv

  # Select both (provide both sweep files):
  python Auto_Select_Network_Parameters.py --mode both \
      --threshold_csv path/to/sweep_detailed.csv \
      --resolution_csv path/to/SC_resolution_sweep.csv

  # Write selected parameters to a config override file:
  python Auto_Select_Network_Parameters.py --mode both \
      --threshold_csv path/to/sweep_detailed.csv \
      --resolution_csv path/to/SC_resolution_sweep.csv \
      --output_config path/to/selected_params.txt

OUTPUT:
  Prints selected parameters to stdout.
  Optionally writes a key=value config file for the bash wrapper to source.

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd


# =============================================================================
# DIFF THRESHOLD SELECTION
# =============================================================================

def select_diff_threshold(sweep_csv, max_lcc=300, comp_fraction=0.80):
    """
    Select optimal DIFF threshold from sweep data.

    Parameters:
      sweep_csv:      Path to sweep_detailed.csv from Step02.1
      max_lcc:        Maximum LCC size for biological interpretability
      comp_fraction:  Minimum fraction of peak components to accept

    Returns:
      dict with selected threshold and supporting metrics
    """
    df = pd.read_csv(sweep_csv)

    # Standardize column names (handle both formats)
    col_map = {}
    for col in df.columns:
        cl = col.strip().lower()
        if 'threshold' in cl or 'thresh' in cl:
            col_map[col] = 'threshold'
        elif 'node' in cl:
            col_map[col] = 'nodes'
        elif 'edge' in cl:
            col_map[col] = 'edges'
        elif 'avgdeg' in cl or 'avg_deg' in cl or 'avg deg' in cl:
            col_map[col] = 'avg_deg'
        elif 'density' in cl:
            col_map[col] = 'density'
        elif cl in ['comp', 'components', 'n_components']:
            col_map[col] = 'components'
        elif cl in ['lcc', 'lcc_size']:
            col_map[col] = 'lcc'
    df = df.rename(columns=col_map)

    # Handle column names from the detailed sweep format
    # The sweep CSV might have columns like: Threshold, Nodes, Edges, AvgDeg, Density, Comp, LCC, ...
    # Try alternate column detection
    if 'components' not in df.columns:
        for col in df.columns:
            if col.strip() in ['Comp', 'comp', 'Components']:
                df = df.rename(columns={col: 'components'})
    if 'lcc' not in df.columns:
        for col in df.columns:
            if col.strip() in ['LCC', 'lcc', 'LCC_size']:
                df = df.rename(columns={col: 'lcc'})
    if 'threshold' not in df.columns:
        for col in df.columns:
            if col.strip() in ['Threshold', 'threshold']:
                df = df.rename(columns={col: 'threshold'})

    # Ensure numeric
    for col in ['threshold', 'nodes', 'edges', 'avg_deg', 'density', 'components', 'lcc']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    print(f"\n  DIFF Threshold Selection")
    print(f"  {'='*60}")
    print(f"  Input: {sweep_csv}")
    print(f"  Rows: {len(df)}")
    print(f"  Max LCC: {max_lcc}")
    print(f"  Min comp fraction: {comp_fraction}")

    # Find component peak
    peak_idx = df['components'].idxmax()
    peak_row = df.loc[peak_idx]
    peak_comp = peak_row['components']
    peak_thresh = peak_row['threshold']
    peak_lcc = peak_row['lcc']

    print(f"\n  Component peak: threshold={peak_thresh}, "
          f"components={int(peak_comp)}, LCC={int(peak_lcc)}")

    # Check if LCC at peak is acceptable
    if peak_lcc <= max_lcc:
        selected = peak_row
        reason = "Component peak with LCC within limit"
    else:
        # Step up to find threshold where LCC < max_lcc
        candidates = df[(df['threshold'] > peak_thresh) & (df['lcc'] <= max_lcc)]
        if len(candidates) > 0:
            # Among candidates, pick the one closest to peak threshold
            # (most permissive that still satisfies LCC constraint)
            selected = candidates.iloc[0]  # Already sorted by threshold ascending
            reason = f"LCC at peak ({int(peak_lcc)}) > {max_lcc}, stepped up to next qualifying threshold"
        else:
            # No threshold gives LCC < max_lcc, use the one closest
            above_peak = df[df['threshold'] > peak_thresh].copy()
            if len(above_peak) > 0:
                above_peak['lcc_diff'] = above_peak['lcc'] - max_lcc
                best = above_peak.loc[above_peak['lcc'].idxmin()]
                selected = best
                reason = f"No threshold gives LCC < {max_lcc}, using smallest LCC above peak"
            else:
                selected = peak_row
                reason = "Fallback to peak (no better option)"

    # Verify component count is reasonable
    selected_comp = selected['components']
    if selected_comp < peak_comp * comp_fraction:
        print(f"  WARNING: Selected threshold has {int(selected_comp)} components "
              f"({selected_comp/peak_comp:.0%} of peak {int(peak_comp)})")

    sel_thresh = selected['threshold']
    sel_nodes = int(selected['nodes'])
    sel_edges = int(selected['edges'])
    sel_avgdeg = selected['avg_deg']
    sel_density = selected['density']
    sel_comp = int(selected['components'])
    sel_lcc = int(selected['lcc'])

    print(f"\n  SELECTED: threshold = {sel_thresh}")
    print(f"    Reason:     {reason}")
    print(f"    Nodes:      {sel_nodes:,}")
    print(f"    Edges:      {sel_edges:,}")
    print(f"    Avg degree: {sel_avgdeg:.1f}")
    print(f"    Density:    {sel_density:.6f}")
    print(f"    Components: {sel_comp}")
    print(f"    LCC:        {sel_lcc}")

    # Show context: thresholds around the selection
    print(f"\n  Context (thresholds near selection):")
    print(f"  {'Threshold':>10s} {'Nodes':>7s} {'Edges':>8s} {'AvgDeg':>7s} "
          f"{'Comp':>5s} {'LCC':>6s} {'Selected':>10s}")
    for _, row in df.iterrows():
        marker = "  <<<" if row['threshold'] == sel_thresh else ""
        print(f"  {row['threshold']:>10.2f} {int(row['nodes']):>7,} {int(row['edges']):>8,} "
              f"{row['avg_deg']:>7.1f} {int(row['components']):>5} "
              f"{int(row['lcc']):>6}{marker}")

    return {
        'threshold': sel_thresh,
        'nodes': sel_nodes,
        'edges': sel_edges,
        'avg_deg': sel_avgdeg,
        'density': sel_density,
        'components': sel_comp,
        'lcc': sel_lcc,
        'reason': reason,
    }


# =============================================================================
# LEIDEN RESOLUTION SELECTION
# =============================================================================

def select_leiden_resolution(resolution_csv, min_ari=0.65):
    """
    Select optimal Leiden resolution using modularity jump detection.

    Method:
      1. Compute delta-Q between consecutive resolutions
      2. Find the resolution where the largest delta-Q occurs
         (the "structural jump" where Leiden finds real communities)
      3. Select that resolution if ARI > stability floor
      4. If ARI fails, step back to previous resolution

    This avoids the trivial-solution trap (r=0.2 with 2 communities)
    because the trivial solution never produces a large Q jump.

    The stability floor (ARI > 0.65) matches the TCGA bulk analysis
    where r=0.80 was selected with ARI=0.63.

    Parameters:
      resolution_csv:  Path to SC_resolution_sweep.csv from Step03
      min_ari:         Minimum ARI stability floor (default: 0.65)

    Returns:
      dict with selected resolution and supporting metrics
    """
    df = pd.read_csv(resolution_csv)

    # Standardize column names
    col_map = {}
    for col in df.columns:
        cl = col.strip().lower()
        if 'resolution' in cl:
            col_map[col] = 'resolution'
        elif 'modularity' in cl and 'std' not in cl:
            col_map[col] = 'modularity'
        elif 'modularity_std' in cl or 'mod_std' in cl:
            col_map[col] = 'modularity_std'
        elif 'ari' in cl and 'std' not in cl:
            col_map[col] = 'ari'
        elif 'nmi' in cl and 'std' not in cl:
            col_map[col] = 'nmi'
        elif 'communit' in cl or 'n_communities' in cl:
            col_map[col] = 'communities'
    df = df.rename(columns=col_map)

    # Ensure numeric and sort by resolution
    for col in ['resolution', 'modularity', 'ari', 'nmi', 'communities']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    df = df.sort_values('resolution').reset_index(drop=True)

    print(f"\n  Leiden Resolution Selection (Modularity Jump Detection)")
    print(f"  {'='*60}")
    print(f"  Input: {resolution_csv}")
    print(f"  Rows: {len(df)}")
    print(f"  ARI stability floor: {min_ari}")

    # Compute delta-Q between consecutive resolutions
    df['delta_q'] = df['modularity'].diff().fillna(0)

    # Compute stability cost (ARI drop)
    df['delta_ari'] = -df['ari'].diff().fillna(0)  # positive = stability lost

    # Efficiency: modularity gained per stability lost
    df['efficiency'] = df['delta_q'] / df['delta_ari'].clip(lower=0.001)

    print(f"\n  Resolution sweep with modularity jumps:")
    print(f"  {'Res':>5s} {'Comms':>6s} {'Q':>7s} {'dQ':>7s} {'ARI':>7s} "
          f"{'dARI':>7s} {'NMI':>7s} {'Effic':>7s} {'Status':>8s}")
    print(f"  {'-'*5} {'-'*6} {'-'*7} {'-'*7} {'-'*7} "
          f"{'-'*7} {'-'*7} {'-'*7} {'-'*8}")

    for _, row in df.iterrows():
        passes = "OK" if row['ari'] >= min_ari else "low ARI"
        print(f"  {row['resolution']:>5.1f} {int(row['communities']):>6} "
              f"{row['modularity']:>7.4f} {row['delta_q']:>+7.4f} "
              f"{row['ari']:>7.4f} {row['delta_ari']:>+7.4f} "
              f"{row['nmi']:>7.4f} {row['efficiency']:>7.2f} {passes:>8s}")

    # Find the largest modularity jump
    # Skip the first row (no delta)
    jump_candidates = df.iloc[1:].copy()
    if len(jump_candidates) == 0:
        print(f"\n  ERROR: Only one resolution tested")
        return {'resolution': df.iloc[0]['resolution']}

    best_jump_idx = jump_candidates['delta_q'].idxmax()
    best_jump = df.loc[best_jump_idx]
    best_jump_dq = best_jump['delta_q']

    print(f"\n  Largest modularity jump: dQ={best_jump_dq:+.4f} at resolution={best_jump['resolution']}")
    print(f"    Q={best_jump['modularity']:.4f}, ARI={best_jump['ari']:.4f}, "
          f"NMI={best_jump['nmi']:.4f}, Communities={int(best_jump['communities'])}")

    # Check if the resolution after the biggest jump passes the stability floor
    if best_jump['ari'] >= min_ari:
        selected = best_jump
        reason = f"Largest modularity jump (dQ={best_jump_dq:+.4f}), ARI passes floor"
    else:
        # Step back to the previous resolution
        prev_idx = best_jump_idx - 1
        if prev_idx >= 0:
            prev = df.loc[prev_idx]
            selected = prev
            reason = (f"Largest jump at r={best_jump['resolution']} "
                      f"(ARI={best_jump['ari']:.3f} < {min_ari}), "
                      f"stepped back to r={prev['resolution']}")
            print(f"  ARI at jump ({best_jump['ari']:.4f}) < floor ({min_ari})")
            print(f"  Stepping back to resolution={prev['resolution']}")
        else:
            selected = best_jump
            reason = "Largest jump (no earlier resolution available)"

    # Also report: what happens one step beyond the jump?
    next_idx = best_jump_idx + 1
    if next_idx < len(df):
        beyond = df.loc[next_idx]
        marginal_q = beyond['modularity'] - selected['modularity']
        marginal_ari = selected['ari'] - beyond['ari']
        print(f"\n  Beyond selected (r={beyond['resolution']}): "
              f"dQ={marginal_q:+.4f}, ARI cost={marginal_ari:+.4f}")
        if marginal_q > 0 and marginal_ari > 0:
            print(f"    Efficiency: {marginal_q/max(marginal_ari,0.001):.2f}x "
                  f"({'worth it' if marginal_q/max(marginal_ari,0.001) > 0.5 else 'diminishing returns'})")

    sel_res = selected['resolution']
    sel_comms = int(selected['communities'])
    sel_q = selected['modularity']
    sel_ari = selected['ari']
    sel_nmi = selected['nmi']

    print(f"\n  SELECTED: resolution = {sel_res}")
    print(f"    Reason:      {reason}")
    print(f"    Communities: {sel_comms}")
    print(f"    Modularity:  {sel_q:.4f}")
    print(f"    ARI:         {sel_ari:.4f}")
    print(f"    NMI:         {sel_nmi:.4f}")

    return {
        'resolution': sel_res,
        'communities': sel_comms,
        'modularity': sel_q,
        'ari': sel_ari,
        'nmi': sel_nmi,
        'reason': reason,
    }


# =============================================================================
# WRITE CONFIG OVERRIDE
# =============================================================================

def write_config(threshold_result, resolution_result, output_path):
    """Write selected parameters to a sourceable config file."""
    with open(output_path, 'w') as f:
        f.write("# Auto-selected network parameters\n")
        f.write(f"# Generated by Auto_Select_Network_Parameters.py\n\n")

        if threshold_result:
            f.write(f"DIFF_THRESHOLD={threshold_result['threshold']}\n")
            f.write(f"DIFF_NODES={threshold_result['nodes']}\n")
            f.write(f"DIFF_EDGES={threshold_result['edges']}\n")
            f.write(f"DIFF_LCC={threshold_result['lcc']}\n")
            f.write(f"DIFF_COMPONENTS={threshold_result['components']}\n")

        if resolution_result:
            f.write(f"LEIDEN_RESOLUTION={resolution_result['resolution']}\n")
            f.write(f"LEIDEN_COMMUNITIES={resolution_result['communities']}\n")
            f.write(f"LEIDEN_MODULARITY={resolution_result['modularity']:.4f}\n")
            f.write(f"LEIDEN_ARI={resolution_result['ari']:.4f}\n")

    print(f"\n  Config written: {output_path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Auto-select DIFF threshold and Leiden resolution")
    parser.add_argument('--mode', choices=['threshold', 'resolution', 'both'],
                        default='both')
    parser.add_argument('--threshold_csv', '--sweep_csv',
                        help='Path to threshold sweep CSV (from Step02.1)')
    parser.add_argument('--resolution_csv',
                        help='Path to resolution sweep CSV (from Step03)')
    parser.add_argument('--output_config',
                        help='Path to write config override file')
    parser.add_argument('--max_lcc', type=int, default=300,
                        help='Maximum LCC size (default: 300)')
    parser.add_argument('--min_ari', type=float, default=0.65,
                        help='Minimum ARI stability floor (default: 0.65)')

    args = parser.parse_args()

    threshold_result = None
    resolution_result = None

    if args.mode in ['threshold', 'both']:
        if not args.threshold_csv:
            print("ERROR: --threshold_csv required for threshold selection")
            sys.exit(1)
        threshold_result = select_diff_threshold(
            args.threshold_csv, max_lcc=args.max_lcc)

    if args.mode in ['resolution', 'both']:
        if not args.resolution_csv:
            print("ERROR: --resolution_csv required for resolution selection")
            sys.exit(1)
        resolution_result = select_leiden_resolution(
            args.resolution_csv,
            min_ari=args.min_ari)

    if args.output_config:
        write_config(threshold_result, resolution_result, args.output_config)

    # Print summary for bash wrapper to capture
    print(f"\n  {'='*60}")
    print(f"  SUMMARY")
    print(f"  {'='*60}")
    if threshold_result:
        print(f"  DIFF_THRESHOLD = {threshold_result['threshold']}")
    if resolution_result:
        print(f"  LEIDEN_RESOLUTION = {resolution_result['resolution']}")


if __name__ == "__main__":
    main()
