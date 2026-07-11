#!/bin/bash
#SBATCH --job-name=STAR_chimeric
#SBATCH --output=/master/jlehle/WORKING/LOGS/STAR_chimeric_%A_%a.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/STAR_chimeric_%A_%a.err
#SBATCH --time=6-00:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=10
#SBATCH --partition=normal
# NOTE: the array range is supplied by the master runner at submit time
#       (sbatch --wait --array=1-N%4 Step05a_STAR_Chimeric_Align.sh),
#       so it is intentionally NOT hardcoded here.
#
# Step05a_STAR_Chimeric_Align.sh
# ================================
# Figure 7: Neoantigen Landscape - STAR chimeric alignment (GENERATOR).
#
# Re-aligns one sample's raw FASTQs with STARsolo in chimeric-detection mode to
# produce Chimeric.out.junction. Cell Ranger BAMs strip supplementary
# alignments, so chimeric reads must be detected by re-aligning from FASTQ.
# This is the missing generation stage: the junctions consumed by Step05b are
# built here at runtime rather than imported from a pre-existing location.
#
# The alignment is GROUP-AGNOSTIC. It aligns every read in the sample against
# the full 10x whitelist; group membership (SBS2_HIGH/CNV_HIGH/NORMAL) is
# applied later, at count time, in Step05b. So these junction files are valid
# regardless of how the three groups are drawn.
#
# Self-contained: paths are defined here, the STAR version is asserted, and a
# sample is skipped only if its own non-empty Chimeric.out.junction already
# exists (idempotent resume). The genome index and whitelist are prepared once,
# before the array, by the master runner and Step01 respectively.
#
# Env: NEOANTIGEN (STAR 2.7.11b). Invoked per array task; SRR is line
#      SLURM_ARRAY_TASK_ID of sample_list.txt, or pass an SRR id as $1.

set -euo pipefail

# =============================================================================
# CONFIGURATION (canonical paths; keep in sync with pipeline_config.yaml)
# =============================================================================
PROJECT_ROOT="/master/jlehle/WORKING/2026_NMF_PAPER"
FASTQ_BASE="/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/fastq/GSE173468"
STAR_CHIMERIC_DIR="${PROJECT_ROOT}/data/FIG_7/04_fusion_analysis/star_chimeric"
INDEX_DIR="${STAR_CHIMERIC_DIR}/genome_index_2.7.11b"
WHITELIST="${STAR_CHIMERIC_DIR}/3M-february-2018.txt"
SAMPLE_LIST="${STAR_CHIMERIC_DIR}/sample_list.txt"

REQUIRED_STAR_VERSION="2.7.11b"
THREADS="${SLURM_CPUS_PER_TASK:-10}"

# Chimeric-detection parameters (identical to the original validated run)
CHIM_SEGMENT_MIN=10
CHIM_JUNCTION_OVERHANG_MIN=10

# =============================================================================
# ENVIRONMENT
# =============================================================================
set +u   # conda's activate scripts reference unbound vars; re-enable after
source ~/anaconda3/bin/activate 2>/dev/null || true
conda activate NEOANTIGEN
set -u

# --- Hard version assert: junction format + barcode recovery depend on it ----
STAR_VERSION="$(STAR --version 2>/dev/null || echo 'MISSING')"
if [ "${STAR_VERSION}" != "${REQUIRED_STAR_VERSION}" ]; then
    echo "FATAL: STAR ${REQUIRED_STAR_VERSION} required, found '${STAR_VERSION}'."
    echo "       The Chimeric.out.junction column layout and CB recovery that"
    echo "       Step05b parses are version-specific. Refusing to run."
    exit 1
fi

# =============================================================================
# RESOLVE SAMPLE
# =============================================================================
if [ -n "${1:-}" ]; then
    SRR_ID="$1"
elif [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    SRR_ID="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")"
else
    echo "FATAL: no SRR provided (pass as \$1 or run as a SLURM array)."
    exit 1
fi

if [ -z "${SRR_ID}" ]; then
    echo "FATAL: empty SRR id (SLURM_ARRAY_TASK_ID out of range for ${SAMPLE_LIST}?)"
    exit 1
fi

echo "=============================================="
echo "STAR chimeric alignment"
echo "  Sample: ${SRR_ID}   STAR: ${STAR_VERSION}   Threads: ${THREADS}"
echo "  Date:   $(date)"
echo "=============================================="

R1="${FASTQ_BASE}/${SRR_ID}/${SRR_ID}_S1_L001_R1_001.fastq.gz"
R2="${FASTQ_BASE}/${SRR_ID}/${SRR_ID}_S1_L001_R2_001.fastq.gz"
OUT_DIR="${STAR_CHIMERIC_DIR}/${SRR_ID}"
JUNCTION="${OUT_DIR}/Chimeric.out.junction"

# --- Preconditions -----------------------------------------------------------
for f in "${R1}" "${R2}"; do
    [ -f "${f}" ] || { echo "FATAL: missing FASTQ ${f}"; exit 1; }
done
[ -f "${INDEX_DIR}/SA" ]        || { echo "FATAL: STAR index missing at ${INDEX_DIR} (master builds it before the array)"; exit 1; }
[ -f "${WHITELIST}" ]           || { echo "FATAL: whitelist missing at ${WHITELIST} (Step01 prepares it)"; exit 1; }

# --- Idempotent resume: skip if a non-empty junction file already exists -----
if [ -s "${JUNCTION}" ]; then
    echo "SKIP: ${JUNCTION} already exists and is non-empty ($(wc -l < "${JUNCTION}") lines)."
    exit 0
fi

# Clean any partial previous attempt so STAR starts fresh
rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}"

echo "R1: ${R1}"
echo "R2: ${R2}"
echo "Output: ${OUT_DIR}"

# =============================================================================
# STAR chimeric alignment
# =============================================================================
STAR \
    --runThreadN "${THREADS}" \
    --genomeDir "${INDEX_DIR}" \
    --readFilesIn "${R2}" "${R1}" \
    --readFilesCommand zcat \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "${WHITELIST}" \
    --soloCBstart 1 \
    --soloCBlen 16 \
    --soloUMIstart 17 \
    --soloUMIlen 12 \
    --soloBarcodeReadLength 0 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
    --chimSegmentMin "${CHIM_SEGMENT_MIN}" \
    --chimJunctionOverhangMin "${CHIM_JUNCTION_OVERHANG_MIN}" \
    --chimOutType Junctions WithinBAM SoftClip \
    --chimOutJunctionFormat 1 \
    --chimMultimapNmax 10 \
    --outFileNamePrefix "${OUT_DIR}/" \
    --limitBAMsortRAM 30000000000

# =============================================================================
# Verify output (so `sbatch --wait` sees a non-zero exit on failure)
# =============================================================================
if [ ! -s "${JUNCTION}" ]; then
    echo "FATAL: STAR finished but produced no Chimeric.out.junction for ${SRR_ID}"
    exit 1
fi

echo "OK: ${SRR_ID} -> ${JUNCTION} ($(wc -l < "${JUNCTION}") lines) at $(date)"
