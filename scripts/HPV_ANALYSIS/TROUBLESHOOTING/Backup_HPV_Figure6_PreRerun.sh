#!/usr/bin/env bash
###############################################################################
# Backup_HPV_Figure6_PreRerun.sh
#
# Backs up and clears ONLY the population-dependent (now stale) Figure 6
# outputs before regenerating them against the new CNV-HIGH group. Everything
# moved is preserved in a timestamped backup dir (nothing is deleted). The
# per-cell Phase 1-4 outputs are population-independent and are left untouched.
#
# WHAT THIS MOVES (stale, built from the OLD CNV-HIGH):
#   - data/FIG_6/FIGURE_6_PANELS/            (Generate_Figure6_Lifecycle_Panels.py output)
#   - the lifecycle diagnostic output dir    (Diagnostic_HPV_Lifecycle_Markers.py output)
#
# WHAT THIS KEEPS (per-cell, group-independent, do NOT regenerate):
#   - data/FIG_6/00_input/adata_v_pp.h5ad
#   - data/FIG_6/01_raw_hpv16_counts/        (master tables, raw counts, summaries)
#   - data/FIG_6/03_hpv16_genome/            (per_cell_hpv16_gene_counts.tsv, alignment, ref)
#
# USAGE (two steps, verify-first):
#   1) Preview (default, moves nothing):
#        bash Backup_HPV_Figure6_PreRerun.sh
#   2) Execute once the preview looks right:
#        DRY_RUN=false bash Backup_HPV_Figure6_PreRerun.sh
#
# BEFORE RUNNING: confirm the diagnostic's real output dir (one line):
#   grep -n 'OUTPUT_DIR' \
#     scripts/HPV_ANALYSIS/TROUBLESHOOTING/Diagnostic_HPV_Lifecycle_Markers.py
#   then set DIAG_OUTPUT_DIR below to whatever path it prints.
###############################################################################

set -euo pipefail

# =============================================================================
# CONFIG
# =============================================================================
PROJECT_ROOT="/master/jlehle/WORKING/2026_NMF_PAPER"
FIG6="${PROJECT_ROOT}/data/FIG_6"
HPV_SCRIPTS="${PROJECT_ROOT}/scripts/HPV_ANALYSIS"
GROUP_FILE="${PROJECT_ROOT}/data/FIG_4/01_group_selection/three_group_assignments.tsv"

# Stale outputs to back up + clear ------------------------------------------
FIG6_PANELS="${FIG6}/FIGURE_6_PANELS"

# !!! VERIFY THIS PATH (see header) before executing. Best-guess default. !!!
DIAG_OUTPUT_DIR="${FIG6}/TROUBLESHOOTING"

# Dry run by default; set DRY_RUN=false to actually move files.
DRY_RUN="${DRY_RUN:-true}"

TS="$(date +%Y%m%d_%H%M%S)"
BACKUP_DIR="${FIG6}/_BACKUP_pre_CNV_rerun_${TS}"
MANIFEST="${BACKUP_DIR}/MANIFEST.txt"

# =============================================================================
# HELPERS
# =============================================================================
log() { printf '%s\n' "$*"; }

run() {  # execute, or just print in dry-run
    if [ "${DRY_RUN}" = "true" ]; then
        log "    [dry-run] $*"
    else
        eval "$@"
    fi
}

ensure_backup_dir() {
    if [ "${DRY_RUN}" != "true" ] && [ ! -d "${BACKUP_DIR}" ]; then
        mkdir -p "${BACKUP_DIR}"
        {
            echo "Backup created: $(date '+%Y-%m-%d %H:%M:%S')"
            echo "Reason: pre-rerun of Figure 6 / lifecycle diagnostic with new CNV-HIGH group"
            echo "Group file: ${GROUP_FILE}"
            echo "Moved items:"
        } > "${MANIFEST}"
    fi
}

backup_and_clear() {  # $1 = path to move into the backup dir
    local src="$1"
    if [ -e "${src}" ]; then
        log "  BACK UP + CLEAR -> ${src}"
        ensure_backup_dir
        run "mv '${src}' '${BACKUP_DIR}/'"
        if [ "${DRY_RUN}" != "true" ]; then
            echo "  ${src}" >> "${MANIFEST}"
        fi
    else
        log "  [skip, not present] ${src}"
    fi
}

# =============================================================================
# GATE: confirm the group file is the regenerated three-group file
# =============================================================================
log "============================================================"
log " PRE-RERUN BACKUP  (DRY_RUN=${DRY_RUN})"
log "============================================================"

if [ ! -f "${GROUP_FILE}" ]; then
    log "ERROR: group file not found: ${GROUP_FILE}"
    exit 1
fi

log ""
log "Group file gate:"
log "  path:     ${GROUP_FILE}"
log "  modified: $(stat -c '%y' "${GROUP_FILE}" 2>/dev/null || date -r "${GROUP_FILE}" '+%Y-%m-%d %H:%M:%S')"
log "  group counts (expect 546 CNV_HIGH, 546 NORMAL, 546 SBS2_HIGH):"
tail -n +2 "${GROUP_FILE}" | cut -f2 | sort | uniq -c | sed 's/^/    /'
TOTAL=$(($(wc -l < "${GROUP_FILE}") - 1))
log "  total cells: ${TOTAL}"
if [ "${TOTAL}" -ne 1638 ]; then
    log "  WARNING: total is not 1638 (3 x 546). Confirm this is the regenerated file before proceeding."
fi

# =============================================================================
# BACK UP + CLEAR the stale population-dependent outputs
# =============================================================================
log ""
log "Backing up stale, population-dependent outputs:"
backup_and_clear "${FIG6_PANELS}"
backup_and_clear "${DIAG_OUTPUT_DIR}"

# =============================================================================
# REASSURANCE: list what is being KEPT (per-cell, group-independent)
# =============================================================================
log ""
log "Keeping (NOT touched, do not regenerate):"
for keep in \
    "${FIG6}/00_input/adata_v_pp.h5ad" \
    "${FIG6}/01_raw_hpv16_counts" \
    "${FIG6}/03_hpv16_genome"
do
    if [ -e "${keep}" ]; then
        log "  [keep] ${keep}"
    else
        log "  [keep, not present] ${keep}"
    fi
done

# =============================================================================
# OPTIONAL: archive clearly superseded Phase5A (old two-population) outputs.
# These are NOT used by the current figure. Enable only if you want them out
# of the way. Set ARCHIVE_SUPERSEDED=true to include them.
# =============================================================================
ARCHIVE_SUPERSEDED="${ARCHIVE_SUPERSEDED:-false}"
if [ "${ARCHIVE_SUPERSEDED}" = "true" ]; then
    log ""
    log "Archiving superseded Phase5A outputs:"
    SUPERSEDED_DIR="${FIG6}/_ARCHIVE_superseded_phase5A_${TS}"
    for s in \
        "${FIG6}/04_population_profiles_v2" \
        "${FIG6}/03_hpv16_genome/hpv16_gene_by_population.tsv"
    do
        if [ -e "${s}" ]; then
            log "  ARCHIVE -> ${s}"
            run "mkdir -p '${SUPERSEDED_DIR}'"
            run "mv '${s}' '${SUPERSEDED_DIR}/'"
        else
            log "  [skip, not present] ${s}"
        fi
    done
fi

# =============================================================================
# DONE
# =============================================================================
log ""
log "============================================================"
if [ "${DRY_RUN}" = "true" ]; then
    log " DRY RUN complete. Nothing moved."
    log " If the list above is correct, run:"
    log "   DRY_RUN=false bash $(basename "$0")"
else
    log " Backup complete."
    log " Backup dir: ${BACKUP_DIR}"
    [ -f "${MANIFEST}" ] && log " Manifest:   ${MANIFEST}"
fi
log "============================================================"
log ""
log "Next, regenerate against the new CNV-HIGH group:"
log "  cd ${HPV_SCRIPTS}"
log "  conda run -n NETWORK python TROUBLESHOOTING/Diagnostic_HPV_Lifecycle_Markers.py"
log "  conda run -n NETWORK python Generate_Figure6_Lifecycle_Panels.py"
