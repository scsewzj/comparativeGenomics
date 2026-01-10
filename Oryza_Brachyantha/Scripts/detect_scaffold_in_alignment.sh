#!/bin/bash

# Base directory containing clusters
base_dir="kaks_input/peptides"

# Iterate over all FASTA files in cluster_* directories
find "$base_dir" -type f -name "*.fasta" | while read fasta; do
    # Check if any header contains 'scaffold'
    if grep -q "^>.*scaffold" "$fasta"; then
        echo "$fasta"
    fi
done
