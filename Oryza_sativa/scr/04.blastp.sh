#!/usr/bin/env bash

echo "==============================================="
echo "Minh Ngoc VU - M2 GENIOMHE Comparative genomics"
echo "==============================================="
echo " "


WORK_DIR="/home/ngoc/comparativeGenomics/Oryza_sativa"


echo "===== $(date): Creating BLAST protein database..."

makeblastdb \
  -in $WORK_DIR/data/osativa.pep.lgiso.fa \
  -dbtype prot \
  -out osativa \
  -title osativa_IRGSP1

mkdir -p $WORK_DIR/data/db
mv osativa* $WORK_DIR/data/db/


echo "===== $(date): Performing BLASTP all-against-all, eval=1e-3, 3 threads, no max result constraints, outfmt 6..."

blastp \
  -query $WORK_DIR/data/osativa.pep.lgiso.fa \
  -db $WORK_DIR/data/db/osativa \
  -outfmt 6 \
  -out $WORK_DIR/data/osativa_blastp_results.txt \
  -evalue 1e-3 \
  -num_threads 3

echo "===== $(date): BLASTP completed, output saved at data/osativa_blastp_results.txt"