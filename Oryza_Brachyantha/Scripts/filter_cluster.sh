#!/bin/bash

# -----------------------------
# Usage
# -----------------------------
usage() {
  cat <<EOF
Usage: $(basename "$0") -i INPUT [-s STRICT] [-m MODERATE] [-r RELAXED]

Description:
  Filter a BLAST TSV file into strict, moderate, and relaxed homologs, 
  and prepare clustering-ready files.

Options:
  -i FILE    Input BLAST TSV file (mandatory)
  -s FILE    Output strict homologs file (default: strict.tsv)
  -m FILE    Output moderate homologs file (default: moderate.tsv)
  -r FILE    Output relaxed homologs file (default: relaxed.tsv)
  -h         Show this help message

Examples:
  # Default filenames
  $(basename "$0") -i treated_blast.tsv

  # Custom filenames
  $(basename "$0") -i treated_blast.tsv -s strict_hits.tsv -m moderate_hits.tsv -r relaxed_hits.tsv
EOF
}

# -----------------------------
# Defaults
# -----------------------------
strict_out="strict.tsv"
moderate_out="moderate.tsv"
relaxed_out="relaxed.tsv"

strict_cluster="strict_clustering.tsv"
moderate_cluster="moderate_clustering.tsv"
relaxed_cluster="relaxed_clustering.tsv"

# -----------------------------
# Parse arguments
# -----------------------------
while getopts "i:s:m:r:h" opt; do
  case $opt in
    i) input="$OPTARG" ;;
    s) strict_out="$OPTARG" ;;
    m) moderate_out="$OPTARG" ;;
    r) relaxed_out="$OPTARG" ;;
    h) usage; exit 0 ;;
    *) usage; exit 1 ;;
  esac
done

# Check mandatory input
if [[ -z "$input" ]]; then
  echo "Error: Input file is mandatory."
  usage
  exit 1
fi

# -----------------------------
# Process all three filters in one awk
# -----------------------------
awk -v strict="$strict_out" -v moderate="$moderate_out" -v relaxed="$relaxed_out" '
{
    qcov = ($8-$7+1)/$13*100
    scov = ($10-$9+1)/$14*100

    # Strict: %identity >=50, qcov>=70, scov>=70, evalue<=1e-10
    if ($3>=50 && qcov>=70 && scov>=70 && $11<=1e-10) print > strict

    # Moderate: %identity >=30, qcov>=50, scov>=50, evalue<=1e-10
    if ($3>=30 && qcov>=50 && scov>=50 && $11<=1e-10) print > moderate

    # Relaxed: %identity >=20, qcov>=30 OR scov>=30, evalue<=1e-5
    if ($3>=20 && (qcov>=30 || scov>=30) && $11<=1e-5) print > relaxed
}
' "$input"

echo "Homolog filtering done!"

# -----------------------------
# Prepare clustering files
# -----------------------------
for cat in strict moderate relaxed; do
    out_file="${cat}.tsv"   # e.g., strict.tsv
    cluster_file="${cat}_pairs.tsv"  # e.g., strict_clustering.tsv

    awk '{print $1 "\t" $2 "\t" $12}' "$out_file" > "$cluster_file"
    echo "Clustering-ready file for $cat homologs -> $cluster_file"
done
