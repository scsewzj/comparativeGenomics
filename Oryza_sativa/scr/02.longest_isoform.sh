#!/usr/bin/env bash

echo "==============================================="
echo "Minh Ngoc VU - M2 GENIOMHE Comparative genomics"
echo "==============================================="
echo " "

echo "===== $(date): Keeping the longest isoform per protein"

WORK_DIR="/home/ngoc/comparativeGenomics/Oryza_sativa"

echo "Number of protein sequences:"
grep -c "^>" $WORK_DIR/data/osativa.pep.filtered.fa

echo "===== $(date): Processing the PEP file"
python scr/02.longest_isoform_pep.py $WORK_DIR/data/osativa.pep.filtered.fa > $WORK_DIR/data/osativa.pep.lgiso.fa

echo "===== $(date): Processing the corresponding CDS file"
python scr/02.longest_isoform_cds.py $WORK_DIR/data/osativa.pep.lgiso.fa $WORK_DIR/data/Oryza_sativa.IRGSP-1.0.cds.all.fa > $WORK_DIR/data/osativa.cds.lgiso.fa

echo "Number of longest isoforms (post-filtering):"
echo "PEP:"
grep -c "^>" $WORK_DIR/data/osativa.pep.lgiso.fa
echo "CDS:"
grep -c "^>" $WORK_DIR/data/osativa.cds.lgiso.fa

echo "===== $(date): Longest isoforms are kept, saved at data/osativa.pep.lgiso.fa and data/osativa.cds.lgiso.fa"

