#!/usr/bin/env bash

echo "==============================================="
echo "Minh Ngoc VU - M2 GENIOMHE Comparative genomics"
echo "==============================================="
echo " "

WORK_DIR="/home/ngoc/comparativeGenomics/Oryza_sativa"


echo "===== $(date): Running CLUSTALW2 on the low-stringent dataset"

mkdir -p $WORK_DIR/paml-workflow/low/clustal_out

find $WORK_DIR/paml-workflow/low/aa -name "*.pep.fa" | \
  parallel -j 16 --bar --eta \
  'base=$(basename {} .pep.fa); clustalw2 -quiet -align -infile={} -outfile='"$WORK_DIR"'/paml-workflow/low/clustal_out/${base}.aln'

ALN_COUNT=$(find $WORK_DIR/paml-workflow/low/clustal_out -name "*.aln" -type f 2>/dev/null | wc -l | tr -d ' ')
echo "===== $(date): CLUSTALW completed with $ALN_COUNT alignments"



echo "===== $(date): Running CLUSTALW2 on the high-stringent dataset"

mkdir -p $WORK_DIR/paml-workflow/high/clustal_out

find $WORK_DIR/paml-workflow/high/aa -name "*.pep.fa" | \
  parallel -j 16 --bar --eta \
  'base=$(basename {} .pep.fa); clustalw2 -quiet -align -infile={} -outfile='"$WORK_DIR"'/paml-workflow/high/clustal_out/${base}.aln'

ALN_COUNT=$(find $WORK_DIR/paml-workflow/high/clustal_out -name "*.aln" -type f 2>/dev/null | wc -l | tr -d ' ')
echo "===== $(date): CLUSTALW completed with $ALN_COUNT alignments"