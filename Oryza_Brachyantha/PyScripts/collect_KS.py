#!/usr/bin/env python3
import os
import glob
import csv
import argparse

# -----------------------------
# Argument parsing
# -----------------------------
parser = argparse.ArgumentParser(
    description="Extract dS, dN, omega, t, kappa and SE values from 2YN output files"
)
parser.add_argument(
    "base_dir",
    help="Input directory containing gene pair files or subdirectories"
)
parser.add_argument(
    "-o", "--out_file",
    default="KaKs_all.tsv",
    help="Output TSV file (default: KaKs_all.tsv)"
)

args = parser.parse_args()

base_dir = args.base_dir
out_file = args.out_file

# -----------------------------
# Check input directory
# -----------------------------
if not os.path.isdir(base_dir):
    print(f"❌ Error: {base_dir} is not a valid directory")
    exit(1)

# -----------------------------
# Open TSV for writing
# -----------------------------
with open(out_file, "w", newline="") as tsvfile:
    writer = csv.writer(tsvfile, delimiter="\t")
    # New header
    writer.writerow(["pair", "dS", "dN", "t", "kappa", "omega", "dN_SE", "dS_SE"])

    # Loop over files or directories
    for pair_path in sorted(glob.glob(os.path.join(base_dir, "*"))):
        pair_name = os.path.basename(pair_path)

        # Assume rst file is either the file itself (if it's a file) or inside directory
        if os.path.isfile(pair_path):
            rst_file = pair_path
        else:
            rst_file = os.path.join(pair_path, "rst")
            if not os.path.isfile(rst_file):
                continue

        # Read the rst file
        with open(rst_file) as f:
            lines = [line.strip() for line in f if line.strip()]

        # Take last line with data
        data_line = lines[-1]
        fields = data_line.split()

        if len(fields) < 9:
            print(f"⚠️ Skipping {rst_file}: not enough columns")
            continue

        # Extract values
        t = fields[1]
        kappa = fields[2]
        omega = fields[3]

        dN, dN_SE = fields[4], fields[6]
        dS, dS_SE = fields[7], fields[9]

        writer.writerow([pair_name, dS, dN, t, kappa, omega, dN_SE, dS_SE])

print(f"✅ dS/dN/omega extracted to {out_file}")
