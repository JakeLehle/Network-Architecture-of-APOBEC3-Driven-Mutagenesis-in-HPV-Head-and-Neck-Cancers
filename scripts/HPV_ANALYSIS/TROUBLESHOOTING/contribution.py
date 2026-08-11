#!/usr/bin/env python3
"""
contribution.py
===============
Single source of truth for the patient-level contribution denominator.

WHY THIS MODULE EXISTS
----------------------
Six scripts compute patient-level fold enrichment for the SBS2-HIGH and
CNV-HIGH groups. Each one used to hardcode its own denominator, which is
exactly the silent-drift failure mode that has bitten this project before.
Every fold in the pipeline now flows through the two functions here, and the
active setting is announced in every log.

THE SETTING
-----------
`CONTRIBUTION_DENOMINATOR` is read from patient_config.py. Two values:

  'all_basal'  expected share = the patient's share of ALL basal cells.
               Asks whether something about the PATIENT drives contribution.
               A3A capability looks partly constitutive (the determinants
               diagnostic finds normal-adjacent tissue already carrying A3A),
               so conditioning the denominator on tumor tissue would partly
               condition on the exposure being measured.

  'tumor'      expected share = the patient's share of TUMOR basal only.
               Conditions on selection eligibility, since Step00B seeds the
               tumor groups from tumor basal exclusively.

Every script computes and logs BOTH folds regardless of the setting. The
setting only decides which one is the headline column and which reference set
the contribution chi-square uses.

VIRUS-DERIVED MEASURES ARE UNAFFECTED. Viral load, lifecycle phase, and every
other virus quantity stay restricted to tumor cells under either setting,
because normal-adjacent basal would dilute them toward zero by construction.

COLUMN CONVENTION
-----------------
Every table written by these scripts carries, for each fate:

    fold_<fate>_all_basal    always computed
    fold_<fate>_tumor        always computed
    fold_<fate>              the ACTIVE one, copied from whichever is selected
    denominator              the string, so a stale table is self-identifying

Downstream scripts read the suffixed column matching the active setting and
fall back to the bare name with a warning.

NO HARDCODED RESULTS
--------------------
`derive_contributors` derives the high-contributor set from the fold column at
runtime. Lists in patient_config are treated as EXPECTATIONS only: a mismatch
is logged loudly and the derived set is what gets used.

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import numpy as np

# -----------------------------------------------------------------------------
# Read the setting from patient_config, with safe fallbacks so a partially
# updated config cannot silently break every script at import time.
# -----------------------------------------------------------------------------
try:
    from patient_config import CONTRIBUTION_DENOMINATOR
except ImportError:
    CONTRIBUTION_DENOMINATOR = 'all_basal'
    print("[contribution] WARNING: CONTRIBUTION_DENOMINATOR not found in "
          "patient_config.py; falling back to 'all_basal'.", flush=True)

try:
    from patient_config import HC_THRESHOLD
except ImportError:
    HC_THRESHOLD = 2.0

VALID_DENOMINATORS = ('all_basal', 'tumor')
if CONTRIBUTION_DENOMINATOR not in VALID_DENOMINATORS:
    raise ValueError(
        f"CONTRIBUTION_DENOMINATOR must be one of {VALID_DENOMINATORS}, "
        f"got {CONTRIBUTION_DENOMINATOR!r} (check patient_config.py)")


def short(p):
    """Strip the 'Patient ' prefix for compact logging."""
    return str(p).replace('Patient ', '')


def fold_enrichment(n_group_p, n_group_total, n_ref_p, n_ref_total):
    """
    Observed share of the group divided by the expected share from the
    reference set. Returns 0.0 when undefined rather than NaN, matching the
    original scripts so no downstream comparison changes behaviour.
    """
    if n_group_total <= 0 or n_ref_total <= 0 or n_ref_p <= 0:
        return 0.0
    return (n_group_p / n_group_total) / (n_ref_p / n_ref_total)


def both_folds(n_group_p, n_group_total,
               n_basal_p, n_basal_total,
               n_tumor_p, n_tumor_total):
    """
    Compute the fold under BOTH denominators plus the active one.

    Returns a dict with keys 'all_basal', 'tumor', 'active'.
    """
    f_all = fold_enrichment(n_group_p, n_group_total, n_basal_p, n_basal_total)
    f_tum = fold_enrichment(n_group_p, n_group_total, n_tumor_p, n_tumor_total)
    return {'all_basal': f_all,
            'tumor': f_tum,
            'active': f_tum if CONTRIBUTION_DENOMINATOR == 'tumor' else f_all}


def reference_count(n_basal_p, n_tumor_p):
    """Per-patient reference count under the active denominator."""
    return int(n_tumor_p) if CONTRIBUTION_DENOMINATOR == 'tumor' else int(n_basal_p)


def reference_total(n_basal_total, n_tumor_total):
    """Cohort reference total under the active denominator."""
    return (int(n_tumor_total) if CONTRIBUTION_DENOMINATOR == 'tumor'
            else int(n_basal_total))


def active_col(prefix):
    """
    Column name of the fold under the active denominator.
    active_col('fold_sbs2') -> 'fold_sbs2_all_basal' or 'fold_sbs2_tumor'
    """
    return f"{prefix}_{CONTRIBUTION_DENOMINATOR}"


def attach_folds(df, prefix):
    """
    Given a DataFrame carrying '<prefix>_all_basal' and '<prefix>_tumor',
    add the bare active column and the denominator tag. Returns the DataFrame.
    """
    need = [f"{prefix}_all_basal", f"{prefix}_tumor"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise KeyError(f"attach_folds({prefix!r}): missing {missing}")
    df[prefix] = df[active_col(prefix)]
    df['denominator'] = CONTRIBUTION_DENOMINATOR
    return df


def read_fold(df, prefix, logger=print, label=''):
    """
    Pull the fold column matching the active denominator out of a table that
    was written by an upstream diagnostic.

    Prefers '<prefix>_<denominator>'. Falls back to the bare '<prefix>' with a
    loud warning, because a bare column means the upstream table predates the
    toggle and its denominator cannot be verified.
    """
    suffixed = active_col(prefix)
    if suffixed in df.columns:
        return df[suffixed]
    if prefix in df.columns:
        tag = df['denominator'].iloc[0] if 'denominator' in df.columns else 'UNKNOWN'
        logger(f"  [WARN] {label}: '{suffixed}' not present; using legacy "
               f"'{prefix}' (table reports denominator = {tag}).")
        logger(f"         Re-run the upstream diagnostic so the denominator is "
               f"explicit before trusting these folds.")
        return df[prefix]
    raise KeyError(f"read_fold: neither '{suffixed}' nor '{prefix}' in "
                   f"{label or 'table'} (have: {list(df.columns)})")


def derive_contributors(df, fold_col, patient_col='patient',
                        threshold=None, expected=None,
                        label='contributors', logger=print):
    """
    Derive a high-contributor set from a fold column at runtime.

    `expected` is an EXPECTATION, never an operative definition. A mismatch is
    logged loudly and the DERIVED set is returned regardless.
    """
    thr = HC_THRESHOLD if threshold is None else threshold
    if fold_col not in df.columns:
        raise KeyError(f"derive_contributors: '{fold_col}' not in table "
                       f"(have: {list(df.columns)})")
    vals = df[fold_col].replace([np.inf, -np.inf], np.nan)
    derived = set(df.loc[vals >= thr, patient_col])
    logger(f"  {label}: derived from '{fold_col}' at >= {thr:.1f}x "
           f"[denominator: {CONTRIBUTION_DENOMINATOR}]")
    logger(f"    -> {sorted(short(p) for p in derived)}")
    if expected is not None:
        exp = set(expected)
        if derived != exp:
            logger(f"    WARNING: disagrees with the patient_config expectation.")
            logger(f"      derived only : {sorted(short(p) for p in derived - exp)}")
            logger(f"      expected only: {sorted(short(p) for p in exp - derived)}")
            logger(f"      The DERIVED set is in use. Update patient_config if "
                   f"this change is intended.")
        else:
            logger(f"    matches the patient_config expectation.")
    return derived


def announce(logger, n_basal_total, n_tumor_total):
    """Standard header block. Every script prints this so no log is ambiguous."""
    pct_norm = (100.0 * (n_basal_total - n_tumor_total) / n_basal_total
                if n_basal_total else 0.0)
    logger(f"  CONTRIBUTION_DENOMINATOR = '{CONTRIBUTION_DENOMINATOR}'")
    logger(f"    all basal   = {n_basal_total:,}")
    logger(f"    tumor basal = {n_tumor_total:,}  "
           f"({pct_norm:.1f}% of basal is normal-adjacent)")
    logger(f"    Both folds are computed and printed; this setting picks the "
           f"headline column and the chi-square reference set.")
    logger(f"    Virus-derived measures (load, lifecycle phase) remain "
           f"tumor-restricted under either setting.")
