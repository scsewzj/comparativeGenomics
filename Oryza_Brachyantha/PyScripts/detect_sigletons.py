#!/usr/bin/env python3

import argparse
import os

# -----------------------------
# Parse command-line arguments
# -----------------------------
parser = argparse.ArgumentParser(
    description="Compare two TSV files and split duplicates and singletons using genes from the first two columns."
)
parser.add_argument("-filtered", required=True, help="Path to the filtered TSV file")
parser.add_argument("-reference", required=True, help="Path to the reference TSV file")
parser.add_argument("-out_duplicates", default=None, help="Output file for duplicates")
parser.add_argument("-out_singletons", default=None, help="Output file for singletons")
parser.add_argument("-delimiter", default="\t", help="Delimiter (default: tab)")

args = parser.parse_args()

# -----------------------------
# Set default output names
# -----------------------------
filtered_base = os.path.splitext(os.path.basename(args.filtered))[0]

if args.out_duplicates is None:
    args.out_duplicates = f"{filtered_base}.duplicate.tsv"

if args.out_singletons is None:
    args.out_singletons = f"{filtered_base}.singleton.tsv"

# -----------------------------
# Function to read first two columns as unique genes
# -----------------------------
def read_first_two_columns(filepath, delimiter):
    genes = set()
    with open(filepath, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split(delimiter)
            for col in parts[:2]:  # take only first two columns
                if col:
                    genes.add(col)
    return genes

# -----------------------------
# Read filtered and reference genes
# -----------------------------
filtered_genes = read_first_two_columns(args.filtered, args.delimiter)
reference_genes = read_first_two_columns(args.reference, args.delimiter)

# -----------------------------
# Compare genes
# -----------------------------
duplicates = sorted(filtered_genes & reference_genes)
singletons = sorted(reference_genes - filtered_genes)

# -----------------------------
# Write output files
# -----------------------------
with open(args.out_duplicates, "w") as f:
    for gene in duplicates:
        f.write(f"{gene}\n")

with open(args.out_singletons, "w") as f:
    for gene in singletons:
        f.write(f"{gene}\n")

print(f"Found {len(duplicates)} duplicates and {len(singletons)} singletons.")
print(f"Duplicates written to: {args.out_duplicates}")
print(f"Singletons written to: {args.out_singletons}")
