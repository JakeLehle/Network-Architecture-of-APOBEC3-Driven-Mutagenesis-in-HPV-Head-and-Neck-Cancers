#!/usr/bin/env bash
# ============================================================================
# fig7_timestamp_audit.sh
# Trace the CNV-HIGH neoantigen chain that feeds the Figure 7 Venn and flag
# any file older than the reselection in three_group_assignments.tsv.
#
# Logic: the diagnostic computes the Venn live from 03_mhc_binding/*_neoantigens.tsv,
# so the chain feeding it is:
#   three_group_assignments.tsv  (canonical selection, the reference clock)
#     -> 01_neoantigen_inputs/barcodes_CNV_HIGH.tsv
#       -> 01_neoantigen_inputs/scomatic_CNV_HIGH.vcf
#         -> 02_snpeff_annotation/CNV_HIGH.somatic_protein_altering.tsv
#           -> 03_mhc_binding/CNV_HIGH_neoantigens.tsv   (what the Venn reads)
# The FIRST file in this list older than the reference is where the rebuild stopped.
# ============================================================================

set -u

BASE="/master/jlehle/WORKING/2026_NMF_PAPER"
FIG7="$BASE/data/FIG_7"
REF_FILE="$BASE/data/FIG_4/01_group_selection/three_group_assignments.tsv"

# Chain in dependency order (earliest input -> final Venn input)
CHAIN=(
  "$FIG7/01_neoantigen_inputs/barcodes_CNV_HIGH.tsv"
  "$FIG7/01_neoantigen_inputs/scomatic_CNV_HIGH.vcf"
  "$FIG7/02_snpeff_annotation/CNV_HIGH.somatic_protein_altering.tsv"
  "$FIG7/03_mhc_binding/CNV_HIGH_neoantigens.tsv"
)

# For cross-reference: the SBS2 and NORMAL equivalents (should be untouched by
# reselection, so these are a control - they can legitimately predate the ref).
CONTROL=(
  "$FIG7/01_neoantigen_inputs/barcodes_SBS2_HIGH.tsv"
  "$FIG7/01_neoantigen_inputs/barcodes_NORMAL.tsv"
  "$FIG7/03_mhc_binding/SBS2_HIGH_neoantigens.tsv"
)

echo "============================================================================"
echo "  FIGURE 7 CNV-HIGH CHAIN: TIMESTAMP AUDIT"
echo "============================================================================"

# --- Reference clock ---------------------------------------------------------
if [[ ! -f "$REF_FILE" ]]; then
  echo "[FATAL] Reference selection file not found:"
  echo "        $REF_FILE"
  echo "        Cannot establish the reselection clock. Aborting."
  exit 1
fi

REF_EPOCH=$(stat -c %Y "$REF_FILE")
REF_HUMAN=$(stat -c '%y' "$REF_FILE")
echo
echo "REFERENCE (reselection clock):"
echo "  three_group_assignments.tsv"
echo "    modified: $REF_HUMAN"
echo

# --- Helper: report one file relative to the reference ----------------------
report_file () {
  local f="$1"
  local tag="$2"   # "CHAIN" or "CONTROL"
  if [[ ! -f "$f" ]]; then
    printf "  [MISSING] %s\n" "$f"
    return
  fi
  local epoch human delta status
  epoch=$(stat -c %Y "$f")
  human=$(stat -c '%y' "$f")
  delta=$(( epoch - REF_EPOCH ))

  if (( epoch < REF_EPOCH )); then
    status="STALE  (older than reselection by $(( -delta / 3600 ))h $(( (-delta % 3600) / 60 ))m)"
  else
    status="fresh  (newer than reselection by $(( delta / 3600 ))h $(( (delta % 3600) / 60 ))m)"
  fi

  printf "  [%s] %-55s\n" "$status" "$(basename "$f")"
  printf "         modified: %s\n" "$human"
  printf "         path:     %s\n" "$f"
}

# --- Walk the chain ----------------------------------------------------------
echo "----------------------------------------------------------------------------"
echo "  CNV-HIGH CHAIN (dependency order; first STALE = where rebuild stopped)"
echo "----------------------------------------------------------------------------"
FIRST_STALE=""
for f in "${CHAIN[@]}"; do
  report_file "$f" "CHAIN"
  if [[ -f "$f" ]]; then
    e=$(stat -c %Y "$f")
    if (( e < REF_EPOCH )) && [[ -z "$FIRST_STALE" ]]; then
      FIRST_STALE="$f"
    fi
  fi
  echo
done

# --- Controls (SBS2 / NORMAL, expected to be unchanged by reselection) -------
echo "----------------------------------------------------------------------------"
echo "  CONTROL FILES (SBS2 / NORMAL - may legitimately predate the reselection)"
echo "----------------------------------------------------------------------------"
for f in "${CONTROL[@]}"; do
  report_file "$f" "CONTROL"
  echo
done

# --- Verdict -----------------------------------------------------------------
echo "============================================================================"
echo "  VERDICT"
echo "============================================================================"
if [[ -n "$FIRST_STALE" ]]; then
  echo "  First stale file in the CNV-HIGH chain:"
  echo "    $FIRST_STALE"
  echo
  echo "  -> The rebuild stopped before this step. Everything from here downstream"
  echo "     (including 03_mhc_binding/CNV_HIGH_neoantigens.tsv that the Venn reads)"
  echo "     is carrying pre-reselection CNV-HIGH membership."
  echo "  -> Rerun the NEOANTIGEN chain from the Step that produces this file."
else
  echo "  No file in the CNV-HIGH chain is older than the reselection."
  echo "  -> Timestamps do NOT explain a stale Venn. Next: membership check on"
  echo "     barcodes_CNV_HIGH.tsv vs the CNV_HIGH rows of three_group_assignments.tsv,"
  echo "     then the hardcode/recompute check. The file may have been re-touched"
  echo "     without its contents actually being rebuilt (timestamp newer, content old)."
fi
echo "============================================================================"
