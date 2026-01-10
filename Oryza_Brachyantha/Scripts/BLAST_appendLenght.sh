#!/bin/bash

# -----------------------------
# Usage
# -----------------------------
usage() {
  cat <<EOF
Usage: $(basename "$0") -i ALL_VS_ALL -f FASTA [-o OUTPUT]

Description:
  Process all-vs-all BLAST hits and peptide FASTA to remove self-hits,
  compute sequence lengths, annotate hits with lengths, and keep the best hit per pair.

Options:
  -i FILE    Input all-vs-all TSV file (mandatory)
  -f FILE    Input peptide FASTA file (mandatory)
  -o FILE    Output TSV file (default: treated_blast.tsv)
  -h         Show this help message

Example:
  $(basename "$0") -i all_vs_all.tsv -f filtered_pep.fasta -o final.tsv
EOF
}

# -----------------------------
# Parse arguments
# -----------------------------
output_file="treated_blast.tsv"

while getopts "i:f:o:h" opt; do
  case $opt in
    i) all_vs_all="$OPTARG" ;;
    f) fasta_file="$OPTARG" ;;
    o) output_file="$OPTARG" ;;
    h) usage; exit 0 ;;
    *) usage; exit 1 ;;
  esac
done

# Check mandatory arguments
if [[ -z "$all_vs_all" || -z "$fasta_file" ]]; then
  echo "Error: Input files are mandatory."
  usage
  exit 1
fi

# -----------------------------
# Main processing
# -----------------------------
awk -v fasta="$fasta_file" '
BEGIN {
    # -----------------------------
    # Step 1: Load sequence lengths
    # -----------------------------
    while ((getline line < fasta) > 0) {
        if (line ~ /^>/) {
            if (seq != "") lengths[id] = length(seq)
            split(line, a, " ")
            id = substr(a[1], 2)
            seq = ""
        } else {
            seq = seq line
        }
    }
    if (seq != "") lengths[id] = length(seq)
    close(fasta)
}

# -----------------------------
# Step 2: Process all-vs-all hits
# -----------------------------
{
    # Normalize pair order FIRST
    if ($1 < $2) {
        g1 = $1
        g2 = $2
    } else {
        g1 = $2
        g2 = $1
    }

    # Skip self-hits (after normalization)
    if (g1 == g2) next

    pair = g1 FS g2

    # Add sequence lengths
    len1 = lengths[g1]
    len2 = lengths[g2]

    # Keep best hit based on column 4
    if (pair in maxlen) {
        if ($4 > maxlen[pair]) {
            maxlen[pair] = $4
            best[pair] = $0 FS len1 FS len2
        }
    } else {
        maxlen[pair] = $4
        best[pair] = $0 FS len1 FS len2
    }
}

END {
    for (p in best) print best[p]
}

' "$all_vs_all" > "$output_file"

echo "Done! Output written to $output_file"
