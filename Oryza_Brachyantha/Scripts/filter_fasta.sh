#!/bin/bash

usage() {
  cat <<EOF
Usage:
  $(basename "$0") [OPTIONS]

Description:
  Filter a FASTA file by skipping sequences whose headers match
  selected conditions. Each skip condition can be enabled independently.

Options:
  --input FILE               Input FASTA file
                             (default: Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa)

  --output FILE              Output FASTA file
                             (default: filtered_pep.fasta)

  --skip-superscaffold       Skip entries with 'superscaffold' in the header
  --skip-mitochond           Skip entries with 'mitochond' in the header
  --skip-chlorop             Skip entries with 'chlorop' in the header
  --skip-arabidopsis         Skip entries with
                             'Source:Projected from Arabidopsis thaliana'

  -h, --help                 Show this help message and exit

Examples:
  # Skip superscaffold entries only
  $(basename "$0") --skip-superscaffold

  # Skip mitochondrial and chloroplast entries
  $(basename "$0") --skip-mitochond --skip-chlorop

  # Use custom input/output and skip all non-nuclear entries
  $(basename "$0") --input input.fa --output output.fa \\
    --skip-superscaffold --skip-mitochond --skip-chlorop --skip-arabidopsis
EOF
}

# Defaults
input_fasta="Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa"
output_fasta="filtered_pep.fasta"

skip_superscaffold=0
skip_mitochond=0
skip_chlorop=0
skip_arabidopsis=0

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input)
      input_fasta="$2"
      shift 2
      ;;
    --output)
      output_fasta="$2"
      shift 2
      ;;
    --skip-superscaffold)
      skip_superscaffold=1
      shift
      ;;
    --skip-mitochond)
      skip_mitochond=1
      shift
      ;;
    --skip-chlorop)
      skip_chlorop=1
      shift
      ;;
    --skip-arabidopsis)
      skip_arabidopsis=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Error: Unknown option '$1'" >&2
      echo
      usage
      exit 1
      ;;
  esac
done

awk \
  -v ss="$skip_superscaffold" \
  -v mito="$skip_mitochond" \
  -v chloro="$skip_chlorop" \
  -v arab="$skip_arabidopsis" '
  /^>/ {
    skip = 0
    if (ss && $0 ~ /superscaffold/) skip = 1
    if (mito && $0 ~ /mitochond/) skip = 1
    if (chloro && $0 ~ /chlorop/) skip = 1
    if (arab && $0 ~ /Source:Projected from Arabidopsis thaliana/) skip = 1
  }
  skip == 0 { print }
' "$input_fasta" > "$output_fasta"
