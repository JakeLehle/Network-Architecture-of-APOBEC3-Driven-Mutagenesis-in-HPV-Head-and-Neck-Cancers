#!/usr/bin/env python3
"""
Diagnostic_Fetch_Protein_Domains.py
===================================
Fetch protein domain / repeat / zinc-finger features from UniProt for the five
Figure 7 genes and cache them as protein_domains.tsv for the gene track.

Everything is rendered in the NATIVE (HGVS / SnpEff-transcript) frame so lollipops
sit at their labeled positions (the numbers used everywhere else in the paper).
For genes whose SnpEff transcript is N-terminally shifted from the UniProt entry
(ANXA1 +11, CAST +83), the domain boxes are shifted by the solved offset K so they
line up with the native lollipop coordinates; the backbone then runs 1..(UniProt
length + K), the leading K residues being the transcript-specific extension.

Feature handling
----------------
- A full feature dump is printed per gene so we can see everything UniProt offers.
- Specific features (Domain, Repeat, Zinc finger, DNA binding, Motif, Coiled coil)
  and named non-disordered Regions are drawn as functional domain boxes; on overlap,
  the more specific type wins (so calpastatin's inhibitory domains are not shadowed
  by the big "Disordered" Regions).
- "Disordered" Regions are emitted as a separate category (distinct muted style),
  never as functional boxes.

Frame safety: mutations are numbered in the SnpEff transcript frame, so we verify
each mutation's reference residue against the UniProt sequence (directly, or after
solving a constant offset). Domains are trusted only when every residue maps.

Inputs (read-only): data/FIG_7/02_snpeff_annotation/{group}.somatic_protein_altering.tsv
Outputs: data/FIG_7/protein_domains/protein_domains.tsv (+ <gene>_mutation_map.tsv)
Env: NETWORK (pandas) on a login node with outbound HTTPS to rest.uniprot.org.
Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import re
import json
import time
import urllib.request
import urllib.parse
from datetime import datetime
import pandas as pd

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
ANNOTATION_DIR = os.path.join(BASE_DIR, "data/FIG_7/02_snpeff_annotation")
OUTPUT_DIR = os.path.join(BASE_DIR, "data/FIG_7/protein_domains")
OUTPUT_TSV = os.path.join(OUTPUT_DIR, "protein_domains.tsv")

FIGURE_GENES = ['KLF3', 'CAST', 'SERPINB2', 'KRT6B', 'ANXA1']
GROUPS = ['SBS2_HIGH', 'CNV_HIGH']

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{acc}.json"
REQUEST_TIMEOUT = 30
RETRIES = 3
RETRY_DELAY = 3

# functional feature types we draw as boxes; Region is included but a Region whose
# description contains "disorder" is split off into its own layer (never a box).
KEEP_TYPES = {'Domain', 'Repeat', 'Zinc finger', 'DNA binding', 'DNA-binding region',
              'Region', 'Motif', 'Coiled coil'}

PALETTE = ['#AED6F1', '#A9DFBF', '#F9E79F', '#D7BDE2', '#F5B7B1',
           '#A3E4D7', '#FAD7A0', '#D5F5E3', '#FADBD8', '#D6EAF8']
GAP_COLOR = '#E8E8E8'         # unannotated backbone
DISORDER_COLOR = '#B0BEC5'    # disordered regions: distinct muted blue-gray

CURATED = {
    'ANXA1': {
        'acc': 'P04083', 'length': 346,
        'domains': [
            (1, 41, 'N-term', '#F5B7B1'),
            (42, 109, 'Repeat 1', '#AED6F1'),
            (110, 122, 'Linker 1-2', '#D5DBDB'),
            (123, 196, 'Repeat 2', '#A9DFBF'),
            (197, 220, 'Linker 2-3', '#D5DBDB'),
            (221, 287, 'Repeat 3', '#F9E79F'),
            (288, 290, 'Linker 3-4', '#D5DBDB'),
            (291, 346, 'Repeat 4', '#D7BDE2'),
        ],
    },
}

AA3TO1 = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q',
          'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
          'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
          'Tyr': 'Y', 'Val': 'V', 'Ter': '*', 'Sec': 'U'}

report = []
def log(msg=""):
    print(msg, flush=True)
    report.append(str(msg))
def banner(title, char="="):
    log(""); log(char * 80); log(f"  {title}"); log(char * 80)


# =============================================================================
# Mutations from the annotation
# =============================================================================
def load_gene_mutations():
    frames = []
    for group in GROUPS:
        p = os.path.join(ANNOTATION_DIR, f"{group}.somatic_protein_altering.tsv")
        if not os.path.exists(p):
            continue
        df = pd.read_csv(p, sep='\t')
        lm = {c.lower(): c for c in df.columns}
        cc = lm.get('chrom', lm.get('#chrom'))
        pc, rc, ac = lm.get('pos'), lm.get('ref'), lm.get('alt')
        gc, hc = lm.get('gene'), lm.get('hgvs_p')
        if not all([cc, pc, gc, hc]):
            continue
        sub = df.rename(columns={cc: 'chrom', pc: 'gpos', rc: 'ref', ac: 'alt',
                                 gc: 'gene', hc: 'hgvs_p'})
        frames.append(sub[['gene', 'chrom', 'gpos', 'ref', 'alt', 'hgvs_p']])
    if not frames:
        return {}
    allm = pd.concat(frames, ignore_index=True)
    xm_pos, wt_aa = [], []
    for h in allm['hgvs_p'].astype(str):
        m = re.match(r'p\.([A-Za-z]{3})(\d+)', h)
        xm_pos.append(int(m.group(2)) if m else None)
        wt_aa.append(AA3TO1.get(m.group(1)) if m else None)
    allm['xm_pos'] = xm_pos
    allm['wt_aa'] = wt_aa
    out = {}
    for gene in FIGURE_GENES:
        g = allm[(allm['gene'] == gene) & allm['xm_pos'].notna() & allm['wt_aa'].notna()]
        g = g.drop_duplicates(subset=['chrom', 'gpos', 'alt']).reset_index(drop=True)
        if len(g):
            out[gene] = g
    return out


# =============================================================================
# UniProt fetch
# =============================================================================
def _get(url):
    req = urllib.request.Request(url, headers={'User-Agent': 'NMF-paper-figure7/1.0'})
    last = None
    for attempt in range(1, RETRIES + 1):
        try:
            with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT) as r:
                return json.loads(r.read().decode('utf-8'))
        except Exception as e:
            last = e
            log(f"    [retry {attempt}/{RETRIES}] {e}")
            time.sleep(RETRY_DELAY)
    log(f"    [ERROR] fetch failed ({last})")
    return None


def uniprot_search(gene, reviewed=True, size=10):
    q = f"gene_exact:{gene} AND organism_id:9606"
    if reviewed:
        q += " AND reviewed:true"
    params = urllib.parse.urlencode({'query': q, 'format': 'json', 'size': size})
    return (_get(f"{UNIPROT_SEARCH}?{params}") or {}).get('results', [])


def uniprot_by_accession(acc):
    return _get(UNIPROT_ENTRY.format(acc=acc))


def parse_entry(entry):
    """Return (acc, seq, length, all_features) with EVERY feature (unfiltered)."""
    acc = entry.get('primaryAccession', '')
    sq = entry.get('sequence', {}) or {}
    seq = sq.get('value', '')
    length = int(sq.get('length', len(seq)))
    feats = []
    for f in entry.get('features', []):
        loc = f.get('location', {})
        try:
            s = int(loc['start']['value']); e = int(loc['end']['value'])
        except (KeyError, TypeError, ValueError):
            continue
        feats.append({'type': f.get('type', ''), 'name': f.get('description', '') or f.get('type', ''),
                      'start': s, 'end': e})
    return acc, seq, length, feats


def pick_best_entry(results, max_mut_pos):
    parsed = [parse_entry(r) for r in results]
    parsed = [p for p in parsed if p[1]]
    if not parsed:
        return None
    spanning = [p for p in parsed if max_mut_pos and p[2] >= max_mut_pos]
    pool = spanning if spanning else parsed
    pool.sort(key=lambda p: p[2], reverse=True)
    return pool[0]


def dump_features(all_features):
    log(f"  UniProt features returned: {len(all_features)} (full list)")
    for f in sorted(all_features, key=lambda x: (x['start'], x['end'])):
        if f['type'] in KEEP_TYPES:
            log(f"    {f['type']:14s} {f['start']:4d}-{f['end']:<4d}  {f['name']}")


def split_features(all_features):
    """functional boxes vs disordered regions (from the drawable KEEP_TYPES)."""
    functional, disordered = [], []
    for f in all_features:
        if f['type'] not in KEEP_TYPES:
            continue
        if f['type'] == 'Region' and 'disorder' in f['name'].lower():
            disordered.append(f)
        else:
            functional.append(f)
    return functional, disordered


def pick_nonoverlapping(feats):
    """Keep a non-overlapping set by position (earliest start, then shortest span).
    This shows finer sub-features (KLF3 repressor domain, KRT6B coils) rather than
    letting a larger overlapping feature (an IF-rod Domain, an inactive motif) swallow
    them. Disordered regions are split off beforehand, so they never compete here,
    which is what lets CAST's inhibitory domains through."""
    kept = []
    for f in sorted(feats, key=lambda f: (f['start'], f['end'])):
        if all(f['end'] < k['start'] or f['start'] > k['end'] for k in kept):
            kept.append(f)
    return kept


# =============================================================================
# Frame check + offset mapping
# =============================================================================
def check_frame(seq, ref_residues):
    nchk = nmat = 0
    for pos, aa in ref_residues:
        if 1 <= pos <= len(seq):
            nchk += 1
            if seq[pos - 1] == aa:
                nmat += 1
    return nchk, nmat


def solve_offset(seq, ref_residues):
    if not ref_residues:
        return None, 0, 0, []
    xs = [p for p, _ in ref_residues]
    lo, hi = min(xs) - len(seq), max(xs) - 1
    bestK, bestN = None, -1
    for K in range(lo, hi + 1):
        n = sum(1 for pos, aa in ref_residues
                if 1 <= pos - K <= len(seq) and seq[pos - K - 1] == aa)
        if n > bestN:
            bestN, bestK = n, K
    per = []
    for pos, aa in ref_residues:
        up = pos - bestK
        got = seq[up - 1] if 1 <= up <= len(seq) else '?'
        per.append((pos, up, aa, got, got == aa))
    return bestK, bestN, len(ref_residues), per


# =============================================================================
# Native-frame track builder
# =============================================================================
def build_track_native(functional, disordered, uniprot_length, K):
    """Emit segments in NATIVE coords (uniprot coord + K). Functional features form
    a contiguous domain/gap partition; disordered regions are overlapping underlay
    rows with a distinct category. Returns (segments, native_length)."""
    native_len = uniprot_length + K
    fn = sorted([{'start': f['start'] + K, 'end': f['end'] + K, 'name': f['name'],
                  'type': f['type']} for f in functional], key=lambda f: (f['start'], f['end']))
    segs, cur, ci = [], 1, 0
    for f in fn:
        if f['start'] > cur:
            segs.append({'start': cur, 'end': f['start'] - 1, 'name': '', 'type': 'unannotated',
                         'category': 'unannotated', 'color': GAP_COLOR})
        segs.append({'start': f['start'], 'end': f['end'], 'name': f['name'], 'type': f['type'],
                     'category': 'domain', 'color': PALETTE[ci % len(PALETTE)]})
        ci += 1; cur = f['end'] + 1
    if cur <= native_len:
        segs.append({'start': cur, 'end': native_len, 'name': '', 'type': 'unannotated',
                     'category': 'unannotated', 'color': GAP_COLOR})
    for d in disordered:
        segs.append({'start': d['start'] + K, 'end': d['end'] + K, 'name': d['name'],
                     'type': 'Disordered', 'category': 'disordered', 'color': DISORDER_COLOR})
    return segs, native_len


def curated_native(spec, K):
    """Curated contiguous domains shifted into native coords, plus the N-terminal
    transcript extension (native 1..K) as unannotated backbone."""
    segs, native_len = [], spec['length'] + K
    if K and K > 0:
        segs.append({'start': 1, 'end': K, 'name': 'N-term extension', 'type': 'unannotated',
                     'category': 'unannotated', 'color': GAP_COLOR})
    for (s, e, nm, col) in spec['domains']:
        cat = 'unannotated' if nm.startswith(('Linker', 'N-term', 'C-term')) else 'domain'
        segs.append({'start': s + K, 'end': e + K, 'name': nm,
                     'type': 'Repeat' if cat == 'domain' else 'unannotated',
                     'category': cat, 'color': col})
    return segs, native_len


def write_mutation_map(gene, gm, per, seq):
    m2 = {int(x): up for x, up, _, _, _ in per}
    mm = gm.copy()
    mm['native_pos'] = mm['xm_pos'].astype(int)      # plot coordinate (HGVS frame)
    mm['uniprot_pos'] = mm['xm_pos'].astype(int).map(m2)
    mm['uniprot_residue'] = [seq[int(u) - 1] if pd.notna(u) and 1 <= int(u) <= len(seq) else ''
                             for u in mm['uniprot_pos']]
    cols = ['gene', 'chrom', 'gpos', 'ref', 'alt', 'hgvs_p', 'xm_pos', 'native_pos',
            'uniprot_pos', 'wt_aa', 'uniprot_residue']
    mm[cols].to_csv(os.path.join(OUTPUT_DIR, f"{gene}_mutation_map.tsv"), sep='\t', index=False)
    return len(mm)


def emit(rows, gene, acc, uniprot_len, native_len, K, numbering_ok, source, segs):
    for i, s in enumerate(segs):
        rows.append({'gene': gene, 'uniprot_acc': acc, 'uniprot_length': uniprot_len,
                     'native_length': native_len, 'frame_offset_K': K, 'numbering_ok': numbering_ok,
                     'source': source, 'seg_index': i, 'category': s['category'],
                     'start': s['start'], 'end': s['end'], 'name': s['name'],
                     'type': s['type'], 'color': s['color']})


def log_segments(segs):
    dom = [s for s in segs if s['category'] == 'domain']
    dis = [s for s in segs if s['category'] == 'disordered']
    log(f"  emitting {len(segs)} native segments ({len(dom)} domain boxes, {len(dis)} disordered):")
    for s in dom:
        log(f"    domain      {s['start']:4d}-{s['end']:<4d}  {s['name']}")
    for s in dis:
        log(f"    disordered  {s['start']:4d}-{s['end']:<4d}  {s['name']}")


# =============================================================================
# MAIN
# =============================================================================
def main():
    banner("FETCH PROTEIN DOMAINS FROM UNIPROT (native frame; specific-over-disordered)")
    log(f"  {datetime.now().isoformat(timespec='seconds')}")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    muts = load_gene_mutations()
    rows = []

    for gene in FIGURE_GENES:
        banner(gene, "-")
        gm = muts.get(gene)
        ref_res = sorted({(int(r.xm_pos), r.wt_aa) for r in gm.itertuples()}) if gm is not None else []
        max_mut = max((p for p, _ in ref_res), default=0)

        # ---- CURATED PATH (ANXA1)
        if gene in CURATED:
            spec = CURATED[gene]
            log(f"  curated override: UniProt {spec['acc']} ({spec['length']} aa)")
            entry = uniprot_by_accession(spec['acc'])
            seq = (entry or {}).get('sequence', {}).get('value', '') if entry else ''
            K, nmat, ntot, per = solve_offset(seq, ref_res) if seq else (None, 0, len(ref_res), [])
            all_ok = (ntot > 0 and nmat == ntot and all(1 <= up <= spec['length'] for _, up, _, _, _ in per))
            log(f"  frame-offset solve: K={K}  (uniprot_pos = xm_pos - K)  matched {nmat}/{ntot}")
            for xm, up, aa, got, okm in per:
                log(f"    {xm:4d} -> {str(up):>4}   {aa} vs {got}   {'OK' if okm else 'MISMATCH'}")
            if all_ok:
                log(f"  -> mapping verified; curated domains in native frame (shift +{K}).")
                if gm is not None:
                    n = write_mutation_map(gene, gm, per, seq)
                    log(f"  wrote {gene}_mutation_map.tsv ({n} lollipops, native positions)")
                segs, native_len = curated_native(spec, K)
                log_segments(segs)
                emit(rows, gene, spec['acc'], spec['length'], native_len, K, True, 'curated', segs)
            else:
                bblen = max(spec['length'], max_mut)
                log(f"  [FALLBACK] mapping not verified; plain backbone at length {bblen}")
                emit(rows, gene, spec['acc'], spec['length'], bblen, K, False, 'curated_fallback',
                     [{'start': 1, 'end': bblen, 'name': '', 'type': 'backbone',
                       'category': 'backbone', 'color': GAP_COLOR}])
            continue

        # ---- AUTOMATED PATH
        results = uniprot_search(gene, reviewed=True) or uniprot_search(gene, reviewed=False)
        best = pick_best_entry(results, max_mut) if results else None
        if not best:
            log(f"  [ERROR] no usable UniProt entry; plain backbone at length {max_mut}")
            emit(rows, gene, '', 0, max_mut, None, False, 'no_entry',
                 [{'start': 1, 'end': max_mut, 'name': '', 'type': 'backbone',
                   'category': 'backbone', 'color': GAP_COLOR}])
            continue

        acc, seq, length, all_features = best
        log(f"  UniProt {acc}  length={length}  (mutations reach {max_mut})")
        dump_features(all_features)
        functional, disordered = split_features(all_features)
        picked = pick_nonoverlapping(functional)
        log(f"  {len(functional)} functional + {len(disordered)} disordered features; "
            f"{len(picked)} functional kept after overlap resolution")

        nchk, nmat = check_frame(seq, ref_res)
        if nchk > 0 and nmat == nchk:
            K, source, ok = 0, 'uniprot', True
            log(f"  frame check: {nmat}/{nchk} match -> native frame = UniProt frame (K=0)")
        else:
            Ks, noff, ntoff, per = solve_offset(seq, ref_res)
            ok = (ntoff > 0 and noff == ntoff and all(1 <= up <= length for _, up, _, _, _ in per))
            log(f"  identity failed; offset solve: K={Ks}  matched {noff}/{ntoff}")
            for xm, up, aa, got, okm in per:
                log(f"    {xm:4d} -> {str(up):>4}   {aa} vs {got}   {'OK' if okm else 'MISMATCH'}")
            if ok:
                K, source = Ks, 'uniprot_offset'
                log(f"  -> colinear with {acc} at K={K}; domains shifted +{K} into native frame")
                if gm is not None:
                    n = write_mutation_map(gene, gm, per, seq)
                    log(f"  wrote {gene}_mutation_map.tsv ({n} lollipops, native positions)")
            else:
                bblen = max(length, max_mut)
                log(f"  -> no clean colinear mapping (best {noff}/{ntoff}): isoform divergence; "
                    f"plain backbone at {bblen}")
                emit(rows, gene, acc, length, bblen, None, False, 'uniprot_fallback',
                     [{'start': 1, 'end': bblen, 'name': '', 'type': 'backbone',
                       'category': 'backbone', 'color': GAP_COLOR}])
                continue

        segs, native_len = build_track_native(picked, disordered, length, K)
        log_segments(segs)
        emit(rows, gene, acc, length, native_len, K, ok, source, segs)

    out = pd.DataFrame(rows)
    out.to_csv(OUTPUT_TSV, sep='\t', index=False)

    banner("SUMMARY")
    for gene in FIGURE_GENES:
        g = out[out['gene'] == gene]
        if len(g) == 0:
            continue
        r0 = g.iloc[0]
        n_dom = int((g['category'] == 'domain').sum())
        n_dis = int((g['category'] == 'disordered').sum())
        status = {'curated': 'curated domains (native, mapping verified)',
                  'uniprot': 'UniProt domains (native)',
                  'uniprot_offset': f"UniProt domains (native, +{int(r0['frame_offset_K'])})",
                  'uniprot_fallback': 'plain backbone (isoform divergence)',
                  'curated_fallback': 'plain backbone (mapping failed)',
                  'no_entry': 'plain backbone (no entry)'}.get(r0['source'], r0['source'])
        log(f"  {gene:9s} {str(r0['uniprot_acc']):8s} native_len={int(r0['native_length']):<5d} "
            f"domains={n_dom:<2d} disordered={n_dis:<2d}  {status}")
    log(f"\n  Wrote {len(out)} rows -> {OUTPUT_TSV}")
    log("  Gene track: x-axis = native HGVS positions; draw domain boxes, underlay")
    log("  'disordered' rows in a distinct muted style, lollipops from <gene>_mutation_map.tsv.")
    with open(os.path.join(OUTPUT_DIR, "fetch_protein_domains_report.txt"), 'w') as fh:
        fh.write('\n'.join(report))


if __name__ == "__main__":
    main()
