#!/bin/bash
# run_star_chimeric.sh — Run STARsolo with chimeric detection for one sample
# Usage: bash run_star_chimeric.sh SRR_ID

SRR_ID=$1
if [ -z "$SRR_ID" ]; then
    echo "Usage: $0 SRR_ID"
    exit 1
fi

echo "Processing $SRR_ID at $(date)"

FASTQ_DIR="/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/fastq/GSE173468/$SRR_ID"
R1="$FASTQ_DIR/${SRR_ID}_S1_L001_R1_001.fastq.gz"
R2="$FASTQ_DIR/${SRR_ID}_S1_L001_R2_001.fastq.gz"
OUT_DIR="/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_6/05_neoantigen/fusion_detection/star_chimeric/$SRR_ID"

mkdir -p "$OUT_DIR"

# Check inputs
if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "ERROR: FASTQs not found for $SRR_ID"
    exit 1
fi

echo "R1: $R1"
echo "R2: $R2"
echo "Output: $OUT_DIR"

# Run STARsolo with chimeric detection
STAR \
    --runThreadN 8 \
    --genomeDir /master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_6/05_neoantigen/fusion_detection/star_chimeric/genome_index_2.7.11b \
    --readFilesIn "$R2" "$R1" \
    --readFilesCommand zcat \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist /master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_6/05_neoantigen/fusion_detection/star_chimeric/3M-february-2018.txt \
    --soloCBstart 1 \
    --soloCBlen 16 \
    --soloUMIstart 17 \
    --soloUMIlen 12 \
    --soloBarcodeReadLength 0 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
    --chimSegmentMin 10 \
    --chimJunctionOverhangMin 10 \
    --chimOutType Junctions WithinBAM SoftClip \
    --chimOutJunctionFormat 1 \
    --chimMultimapNmax 10 \
    --outFileNamePrefix "$OUT_DIR/" \
    --limitBAMsortRAM 30000000000

echo "STAR complete for $SRR_ID at $(date)"

# Check output
if [ -f "$OUT_DIR/Chimeric.out.junction" ]; then
    N_CHIM=$(wc -l < "$OUT_DIR/Chimeric.out.junction")
    echo "Chimeric junctions: $N_CHIM"
else
    echo "WARNING: No Chimeric.out.junction produced"
fi
