#!/usr/bin/env python3
"""
enrichr_kegg_wrapper.py  (proposed helper for the DEG -> KEGG stage)
===================================================================

Standardized Enrichr/KEGG call for the single-cell DEG enrichment. It merges
the two rate-limit patterns already in the repo:

  - preemptive fixed delay before every call
    (from TROUBLESHOOTING/KEGG_A3_Neighborhood_Diagnostic.py)
  - exponential backoff retry on HTTP 429
    (from scripts/NETWORK/Step08_Pipeline_Summary.py)

Why it differs from the diagnostic's swallow-and-return-empty:
  The A3 neighborhood diagnostic fires dozens of calls and can afford to drop
  a single failed community and keep going. The DEG stage makes only a handful
  of calls (2 comparisons, optionally split up/down), so a swallowed 429 would
  silently blank an entire comparison's KEGG and read as "no enrichment" on the
  figure. That is the same silent-masking failure mode as the stale SnpEff
  cache. So here: retry on 429, then FAIL LOUD (raise) rather than return empty.
  An empty frame is returned ONLY for the legitimate "too few genes to test"
  case, and that is logged explicitly.

Standardized on the canonical 15 s spacing (the diagnostic used 16 s).

Intended to be embedded in the DEG KEGG compute script (matching how
KEGG_A3_Neighborhood_Diagnostic.py and Step08 each embed their own wrapper),
or factored into a shared helper if preferred.
"""

import time
import pandas as pd

KEGG_LIBRARY = "KEGG_2021_Human"

# Rate-limit conventions (canonical 15 s spacing, matching Step08).
ENRICHR_PRE_DELAY   = 15    # seconds slept BEFORE every call (preemptive throttle)
ENRICHR_MAX_RETRIES = 5     # attempts total on HTTP 429
ENRICHR_BACKOFF     = 30    # first backoff (s); doubles each retry: 30/60/120/240
MIN_GENES_ENRICHR   = 5     # below this, nothing to test


def run_enrichr_kegg(gene_list, description="query",
                     background=None,
                     pre_delay=ENRICHR_PRE_DELAY,
                     max_retries=ENRICHR_MAX_RETRIES,
                     backoff=ENRICHR_BACKOFF,
                     logfn=print):
    """
    KEGG pathway enrichment via gseapy/Enrichr, throttled and retried.

    Parameters
    ----------
    gene_list : iterable of str
        Gene symbols. ENSG-looking IDs are dropped (unmapped).
    description : str
        Label for logging (e.g. "SBS2_VS_NORMAL_up").
    background : iterable of str or None
        Optional custom background (e.g. all tested genes). None keeps Enrichr's
        default genome-wide background, matching the existing repo wrappers.
        NOTE: with some gseapy versions, passing a background switches enrichr to
        a local hypergeometric computation rather than the web service; confirm
        behavior on the HPC gseapy version before relying on it (decision 2).
    logfn : callable
        Logger (defaults to print; pass the script's `log` to capture in report).

    Returns
    -------
    pd.DataFrame
        Enrichr results sorted by Adjusted P-value. Empty ONLY when the input is
        below MIN_GENES_ENRICHR (logged).

    Raises
    ------
    RuntimeError
        If a non-429 error occurs, or all 429 retries are exhausted. Loud on
        purpose so a failed comparison is never mistaken for "no enrichment".
    """
    import gseapy as gp

    gene_list = [g for g in gene_list
                 if isinstance(g, str) and not g.startswith("ENSG")]
    if len(gene_list) < MIN_GENES_ENRICHR:
        logfn(f"    [{description}] only {len(gene_list)} genes "
              f"(< {MIN_GENES_ENRICHR}); skipping enrichment.")
        return pd.DataFrame()

    kwargs = dict(
        gene_list=list(gene_list),
        gene_sets=KEGG_LIBRARY,
        organism="human",
        outdir=None,
        no_plot=True,
        verbose=False,
    )
    if background is not None:
        kwargs["background"] = list(background)

    delay = backoff
    for attempt in range(1, max_retries + 1):
        time.sleep(pre_delay)          # preemptive throttle before EVERY call
        try:
            enr = gp.enrichr(**kwargs)
            df = enr.results
            if len(df) > 0:
                df = df.sort_values("Adjusted P-value").reset_index(drop=True)
            logfn(f"    [{description}] {len(df)} pathways tested, "
                  f"{int((df['Adjusted P-value'] < 0.05).sum()) if len(df) else 0} "
                  f"sig (FDR<0.05).")
            return df
        except Exception as e:
            is_429 = "429" in str(e) or "Too Many Requests" in str(e)
            if is_429 and attempt < max_retries:
                logfn(f"    [{description}] 429 on attempt {attempt}/"
                      f"{max_retries}; backing off {delay}s")
                time.sleep(delay)
                delay *= 2
                continue
            raise RuntimeError(
                f"Enrichr KEGG failed for '{description}' after "
                f"{attempt} attempt(s): {e}"
            )
