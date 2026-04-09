#!/bin/bash
# Troubleshoot_SnpEff_STAR.sh
# Test fixes for both issues. Run interactively in NEOANTIGEN env.

echo "=============================================="
echo "Fix 1: SnpEff Java heap space"
echo "=============================================="

# Test with 200GB heap on first 10 lines of a VCF
VCF="/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_6/05_neoantigen/inputs/scomatic_SBS2_HIGH.vcf"
TEST_VCF="/tmp/test_snpeff_10lines.vcf"

# Extract header + first 10 data lines
grep "^#" "$VCF" > "$TEST_VCF"
grep -v "^#" "$VCF" | head -10 >> "$TEST_VCF"

echo "Test VCF: $(wc -l < $TEST_VCF) lines"
echo "Running: snpEff -Xmx200g ann GRCh38.p14 ..."

snpEff -Xmx200g ann -noStats -canon GRCh38.p14 "$TEST_VCF" > /tmp/test_snpeff_out.vcf 2>/tmp/test_snpeff_err.txt
EXITCODE=$?

echo "Exit code: $EXITCODE"
if [ $EXITCODE -eq 0 ]; then
    echo "SUCCESS! SnpEff works with -Xmx200g"
    echo "Output lines: $(wc -l < /tmp/test_snpeff_out.vcf)"
    echo "First annotation:"
    grep -v "^#" /tmp/test_snpeff_out.vcf | head -1 | cut -f8 | tr ';' '\n' | grep "^ANN=" | head -1 | cut -c1-200
else
    echo "FAILED. Error:"
    cat /tmp/test_snpeff_err.txt | tail -20
fi

echo ""
echo "=============================================="
echo "Fix 2: STAR whitelist (decompress .gz)"  
echo "=============================================="

WL_GZ="/master/jlehle/cellranger-8.0.1/lib/python/cellranger/barcodes/3M-february-2018.txt.gz"
WL_OUT="/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_6/05_neoantigen/fusion_detection/star_chimeric/3M-february-2018.txt"

if [ -f "$WL_OUT" ]; then
    echo "Decompressed whitelist already exists: $WL_OUT"
else
    echo "Decompressing: $WL_GZ"
    zcat "$WL_GZ" > "$WL_OUT"
    echo "Done. Lines: $(wc -l < $WL_OUT)"
fi

echo "First 3 barcodes:"
head -3 "$WL_OUT"
echo "Total barcodes: $(wc -l < $WL_OUT)"

# Quick STAR test with decompressed whitelist
echo ""
echo "Testing STAR with decompressed whitelist..."

# Just check that STAR accepts the whitelist (dry run with tiny input)
TEST_R1="/tmp/test_r1.fq"
TEST_R2="/tmp/test_r2.fq"

# Create minimal test FASTQs
echo -e "@read1\nACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIII" > "$TEST_R1"
echo -e "@read1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" > "$TEST_R2"

STAR \
    --runThreadN 1 \
    --genomeDir /master/jlehle/WORKING/SC/ref/GRCh38/star \
    --readFilesIn "$TEST_R2" "$TEST_R1" \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "$WL_OUT" \
    --soloCBstart 1 --soloCBlen 16 \
    --soloUMIstart 17 --soloUMIlen 12 \
    --chimSegmentMin 10 \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix /tmp/star_test/ \
    --runMode alignReads \
    2>/tmp/star_test_err.txt

STAR_EXIT=$?
echo "STAR exit code: $STAR_EXIT"
if [ $STAR_EXIT -eq 0 ] || grep -q "completed successfully" /tmp/star_test_err.txt 2>/dev/null; then
    echo "SUCCESS! STAR accepts the decompressed whitelist"
else
    echo "Error output:"
    tail -5 /tmp/star_test_err.txt
fi

# Clean up
rm -rf /tmp/star_test/ /tmp/test_r1.fq /tmp/test_r2.fq /tmp/test_snpeff_*.vcf /tmp/test_snpeff_*.txt

echo ""
echo "=============================================="
echo "Summary"
echo "=============================================="
echo "SnpEff fix: add -Xmx200g flag to snpEff command"
echo "STAR fix:   use decompressed whitelist at $WL_OUT"
echo ""
echo "If both passed, update scripts and re-run."
