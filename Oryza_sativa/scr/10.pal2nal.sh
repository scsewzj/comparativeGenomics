#!/usr/bin/env bash

echo "==============================================="
echo "Minh Ngoc VU - M2 GENIOMHE Comparative genomics"
echo "==============================================="
echo " "

WORK_DIR="/home/ngoc/comparativeGenomics/Oryza_sativa/paml-workflow"
export WORK_DIR


echo "===== $(date): Running PAL2NAL to align CDS on the low-stringent dataset..."

mkdir -p $WORK_DIR/low/pal2nal_out

run_pal2nal() {
    cds_file=$1
    base=$(basename "$cds_file" .cds.fa)
    $WORK_DIR/pal2nal.pl \
        $WORK_DIR/low/clustal_out/${base}.aln \
        $cds_file \
        -output paml > $WORK_DIR/low/pal2nal_out/${base}.phy
}

export -f run_pal2nal

find $WORK_DIR/low/cds -name "*.cds.fa" | parallel -j 16 --bar run_pal2nal

PHY_COUNT=$(find $WORK_DIR/low/pal2nal_out -name "*.phy" 2>/dev/null | wc -l | tr -d ' ')
echo "===== $(date): PAL2NAL completed for the low-stringent dataset with $PHY_COUNT alignments."


echo "===== $(date): Running PAL2NAL to align CDS on the high-stringent dataset..."

mkdir -p $WORK_DIR/high/pal2nal_out

run_pal2nal() {
    cds_file=$1
    base=$(basename "$cds_file" .cds.fa)
    $WORK_DIR/pal2nal.pl \
        $WORK_DIR/high/clustal_out/${base}.aln \
        $cds_file \
        -output paml > $WORK_DIR/high/pal2nal_out/${base}.phy
}

find $WORK_DIR/high/cds -name "*.cds.fa" | parallel -j 16 --bar run_pal2nal

PHY_COUNT=$(find $WORK_DIR/high/pal2nal_out -name "*.phy" 2>/dev/null | wc -l | tr -d ' ')
echo "===== $(date): PAL2NAL completed for the high-stringent dataset with $PHY_COUNT alignments."