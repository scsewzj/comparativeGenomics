#!/usr/bin/env bash

echo "==============================================="
echo "Minh Ngoc VU - M2 GENIOMHE Comparative genomics"
echo "==============================================="
echo " "

echo "===== $(date): Filtering the proteome"

WORK_DIR="/home/ngoc/comparativeGenomics/Oryza_sativa"

echo "===== $(date): Initial assessment:"

echo "Number of sequences:"
grep -c "^>" $WORK_DIR/data/Oryza_sativa.IRGSP-1.0.pep.all.fa

echo "Number of sequences per chromosome:"
awk '/^>/' $WORK_DIR/data/Oryza_sativa.IRGSP-1.0.pep.all.fa | \
  awk '{print $3}' | \
  awk -F':' '{print $2":"$3}' | \
  sort -d | \
  uniq -c

echo "Number of sequences per transcript_biotype:"
awk '/^>/' $WORK_DIR/data/Oryza_sativa.IRGSP-1.0.pep.all.fa | \
  awk '{print $7}' | \
  awk -F':' '{print $2}' | \
  sort -d | \
  uniq -c


echo "===== $(date): Filtering nontranslating and Mt Pt sequences..."
awk '
/^>/ {
    keep = ($3 !~ /:(Mt|Pt):/ && $7 !~ /nontranslating_CDS/)
}
keep
' data/Oryza_sativa.IRGSP-1.0.pep.all.fa > data/osativa.pep.filtered.fa

echo "===== $(date): Post-filtering assessment:"
echo "Number of sequences:"
grep -c "^>" $WORK_DIR/data/osativa.pep.filtered.fa

echo "Number of sequences per chromosome:"
awk '/^>/' $WORK_DIR/data/osativa.pep.filtered.fa | \
  awk '{print $3}' | \
  awk -F':' '{print $2":"$3}' | \
  sort -d | \
  uniq -c

echo "Number of sequences per transcript_biotype:"
awk '/^>/' $WORK_DIR/data/osativa.pep.filtered.fa | \
  awk '{print $7}' | \
  awk -F':' '{print $2}' | \
  sort -d | \
  uniq -c

echo "===== $(date): Filtering completed, saved at data/osativa.pep.filtered.fa."

