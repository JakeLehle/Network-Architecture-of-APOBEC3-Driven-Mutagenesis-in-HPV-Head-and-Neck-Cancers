#!/usr/bin/env bash
# ============================================================================
# fig7_membership_check.sh
# Confirm the CNV_HIGH barcodes used by the FIG_7 neoantigen pipeline match
# the reselected CNV_HIGH 546 in the canonical FIG_4 selection file.
#
# Decisive test for:
#   (b) rebuild ran on OLD membership  -> barcodes will NOT match
#   (c) Venn genuinely robust          -> barcodes WILL match (the new 546)
# ============================================================================

set -u

BASE="/master/jlehle/WORKING/2026_NMF_PAPER"
SELECTION="$BASE/data/FIG_4/01_group_selection/three_group_assignments.tsv"
BARCODES="$BASE/data/FIG_7/01_neoantigen_inputs/barcodes_CNV_HIGH.tsv"

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

echo "============================================================================"
echo "  FIG_7 CNV_HIGH MEMBERSHIP CHECK"
echo "============================================================================"

for f in "$SELECTION" "$BARCODES"; do
  if [[ ! -f "$f" ]]; then
    echo "[FATAL] Missing required file: $f"
    exit 1
  fi
done

# --- Inspect the selection file header so we parse the right columns ---------
echo
echo "Selection file header (first line) and a sample row:"
head -n 1 "$SELECTION" | sed 's/\t/ | /g'
echo "  ---"
sed -n '2p' "$SELECTION" | sed 's/\t/ | /g'
echo

# --- Locate the barcode column and the group/label column -------------------
# We do NOT assume positions. Find a column whose header looks like a barcode
# id and one whose header looks like the group assignment.
HEADER=$(head -n 1 "$SELECTION")

# Candidate names for the barcode column and the group column (case-insensitive)
bc_col=$(awk -F'\t' 'NR==1{
  for(i=1;i<=NF;i++){
    h=tolower($i)
    if(h=="barcode"||h=="cell"||h=="cell_id"||h=="cellid"||h=="index"||h=="cell_barcode"){print i; exit}
  }
}' "$SELECTION")

grp_col=$(awk -F'\t' 'NR==1{
  for(i=1;i<=NF;i++){
    h=tolower($i)
    if(h=="group"||h=="assignment"||h=="label"||h=="tier"||h=="population"||h=="class"){print i; exit}
  }
}' "$SELECTION")

if [[ -z "${bc_col:-}" || -z "${grp_col:-}" ]]; then
  echo "[WARN] Could not auto-detect columns from header."
  echo "       Detected barcode col: '${bc_col:-NONE}', group col: '${grp_col:-NONE}'"
  echo "       Columns in file:"
  echo "$HEADER" | tr '\t' '\n' | nl
  echo
  echo "       Edit bc_col / grp_col manually and re-run. Aborting to avoid a"
  echo "       false match on the wrong columns."
  exit 2
fi

echo "Auto-detected: barcode column = $bc_col, group column = $grp_col"
echo

# --- What group labels exist, and which is CNV_HIGH? ------------------------
echo "Group labels present in selection file (with counts):"
awk -F'\t' -v g="$grp_col" 'NR>1{print $g}' "$SELECTION" | sort | uniq -c | sort -rn
echo

# Match CNV_HIGH tolerantly (CNV_HIGH, CNV-HIGH, CNV HIGH, cnv_high...)
CNV_LABEL=$(awk -F'\t' -v g="$grp_col" 'NR>1{print $g}' "$SELECTION" \
  | sort -u \
  | grep -iE 'cnv[ _-]?high' | head -n1)

if [[ -z "$CNV_LABEL" ]]; then
  echo "[FATAL] No CNV_HIGH-like label found in group column. Labels above."
  exit 3
fi
echo "Using CNV-HIGH label: '$CNV_LABEL'"
echo

# --- Normalize barcodes on both sides ---------------------------------------
# Strip a trailing '-1' (or '-N') Cell Ranger suffix and any sample prefix
# ('SAMPLE_AAAC...') so the two files compare on the raw 16bp where possible.
# We keep BOTH raw and normalized comparisons so you can see if a mismatch is
# real or just a formatting artifact.

normalize () {
  # reads barcodes on stdin, writes normalized to stdout
  sed -E 's/-[0-9]+$//' | sed -E 's/^.*_([ACGT]{16})$/\1/'
}

# CNV_HIGH barcodes from the canonical selection
awk -F'\t' -v g="$grp_col" -v b="$bc_col" -v lab="$CNV_LABEL" \
  'NR>1 && $g==lab {print $b}' "$SELECTION" | sort -u > "$WORK/sel_raw.txt"

# Barcodes the pipeline actually used (file may be headerless; handle both)
# If the first line looks like a header (non-ACGT), drop it.
first="$(head -n1 "$BARCODES")"
if [[ "$first" =~ ^[ACGTacgt]{8,} || "$first" =~ [ACGT]{16} ]]; then
  cat "$BARCODES" | sort -u > "$WORK/bc_raw.txt"
else
  tail -n +2 "$BARCODES" | sort -u > "$WORK/bc_raw.txt"
fi

# Normalized versions
normalize < "$WORK/sel_raw.txt" | sort -u > "$WORK/sel_norm.txt"
normalize < "$WORK/bc_raw.txt"  | sort -u > "$WORK/bc_norm.txt"

n_sel_raw=$(wc -l < "$WORK/sel_raw.txt")
n_bc_raw=$(wc -l < "$WORK/bc_raw.txt")
n_sel_norm=$(wc -l < "$WORK/sel_norm.txt")
n_bc_norm=$(wc -l < "$WORK/bc_norm.txt")

echo "----------------------------------------------------------------------------"
echo "  COUNTS"
echo "----------------------------------------------------------------------------"
printf "  Selection CNV_HIGH barcodes : %5d raw / %5d normalized\n" "$n_sel_raw" "$n_sel_norm"
printf "  Pipeline  CNV_HIGH barcodes : %5d raw / %5d normalized\n" "$n_bc_raw" "$n_bc_norm"
echo "  (expect 546 on each side)"
echo

# --- Overlap on raw, then normalized ----------------------------------------
report_overlap () {
  local sel="$1" bc="$2" tag="$3"
  local shared only_sel only_bc
  shared=$(comm -12 "$sel" "$bc" | wc -l)
  only_sel=$(comm -23 "$sel" "$bc" | wc -l)
  only_bc=$(comm -13 "$sel" "$bc" | wc -l)
  echo "----------------------------------------------------------------------------"
  echo "  OVERLAP ($tag)"
  echo "----------------------------------------------------------------------------"
  printf "  Shared (in both)              : %5d\n" "$shared"
  printf "  Only in FIG_4 selection       : %5d\n" "$only_sel"
  printf "  Only in FIG_7 pipeline file   : %5d\n" "$only_bc"
  echo
  # stash for verdict
  echo "$shared $only_sel $only_bc" > "$WORK/ov_${tag}.txt"
}

report_overlap "$WORK/sel_raw.txt" "$WORK/bc_raw.txt" "raw"
report_overlap "$WORK/sel_norm.txt" "$WORK/bc_norm.txt" "normalized"

# --- Verdict ----------------------------------------------------------------
read -r sh os ob < "$WORK/ov_normalized.txt"
echo "============================================================================"
echo "  VERDICT"
echo "============================================================================"
if (( os == 0 && ob == 0 && sh > 0 )); then
  echo "  PERFECT MATCH (normalized): all $sh CNV_HIGH barcodes identical."
  echo "  -> Case (c): the FIG_7 pipeline used the RESELECTED CNV_HIGH 546."
  echo "     The unchanged Venn (276/105/135) is GENUINE robustness, not a stale"
  echo "     artifact. You are in the clear. Consider noting in the text that the"
  echo "     neoantigen gene partition is stable to CNV-HIGH reselection."
elif (( sh == 0 )); then
  echo "  ZERO overlap. The two files share no barcodes at all."
  echo "  -> Likely a formatting/namespace mismatch (suffix, prefix, or the"
  echo "     selection uses a different cell-id scheme). Inspect the sample rows"
  echo "     printed above before concluding. Do NOT interpret as biology yet."
else
  echo "  PARTIAL / NON-MATCH (normalized): shared=$sh, only_FIG4=$os, only_FIG7=$ob"
  echo "  -> Case (b): the rebuild ran on a DIFFERENT CNV_HIGH set than the"
  echo "     corrected FIG_4 selection. The Venn is stale-by-wrong-input despite"
  echo "     fresh timestamps. Trace which selection source Step01 actually read"
  echo "     (grep the NEOANTIGEN Step01 script for the path it loads), then rerun"
  echo "     the chain from Step01 against three_group_assignments.tsv."
  echo
  echo "  Sample of barcodes ONLY in FIG_4 selection (missing from pipeline):"
  comm -23 "$WORK/sel_norm.txt" "$WORK/bc_norm.txt" | head -n 5 | sed 's/^/    /'
  echo "  Sample of barcodes ONLY in FIG_7 pipeline (extra / from old set):"
  comm -13 "$WORK/sel_norm.txt" "$WORK/bc_norm.txt" | head -n 5 | sed 's/^/    /'
fi
echo "============================================================================"
