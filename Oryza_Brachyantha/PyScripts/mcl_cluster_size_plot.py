#!/usr/bin/env python3

import os
import argparse
from collections import Counter
import matplotlib.pyplot as plt

# =========================
# Argument parsing
# =========================
parser = argparse.ArgumentParser(
    description="Plot MCL cluster size distribution"
)

parser.add_argument(
    "-i", "--input",
    required=True,
    help="MCL cluster file (one cluster per line)"
)

parser.add_argument(
    "-o", "--outdir",
    default="plots",
    help="Output directory (default: plots/)"
)

args = parser.parse_args()

infile = args.input
outdir = args.outdir

# =========================
# Prepare output path
# =========================
os.makedirs(outdir, exist_ok=True)

base = os.path.basename(infile)
name = os.path.splitext(base)[0]
outfile = os.path.join(outdir, f"cluster_distribution_{name}.png")

# =========================
# Read MCL clusters
# =========================
cluster_sizes = []

with open(infile) as f:
    for line in f:
        if line.strip():
            genes = line.strip().split()
            cluster_sizes.append(len(genes))

# =========================
# Count clusters per size
# =========================
size_counts = Counter(cluster_sizes)

sizes = sorted(size_counts)
counts = [size_counts[s] for s in sizes]

# =========================
# Plot
# =========================
plt.figure(figsize=(7, 5))
plt.bar(sizes, counts)
plt.xlabel("Cluster size (number of genes)")
plt.ylabel("Number of clusters")
plt.title("MCL Cluster Size Distribution")

plt.tight_layout()
plt.savefig(outfile, dpi=300)
plt.close()

print(f"Plot saved to: {outfile}")
