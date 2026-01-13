#!/usr/bin/env python3

import os
import argparse
from collections import Counter
import matplotlib.pyplot as plt

# =========================
# Argument parsing
# =========================
parser = argparse.ArgumentParser(
    description="Plot MCL cluster size distribution and count genes per category"
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
cluster_genes = []

with open(infile) as f:
    for line in f:
        if line.strip():
            genes = line.strip().split()
            cluster_sizes.append(len(genes))
            cluster_genes.append(genes)

# =========================
# Count clusters per size
# =========================
size_counts = Counter(cluster_sizes)

sizes = sorted(size_counts)
counts = [size_counts[s] for s in sizes]

# =========================
# Count clusters and genes per category
# =========================
def count_categories_with_genes(clusters):
    categories = {
        "size_2": {"clusters": 0, "genes": set()},
        "size_3_9": {"clusters": 0, "genes": set()},
        "size_10_19": {"clusters": 0, "genes": set()},
        "size_20_plus": {"clusters": 0, "genes": set()}
    }

    for genes in clusters:
        size = len(genes)
        if size == 2:
            categories["size_2"]["clusters"] += 1
            categories["size_2"]["genes"].update(genes)
        elif 3 <= size <= 9:
            categories["size_3_9"]["clusters"] += 1
            categories["size_3_9"]["genes"].update(genes)
        elif 10 <= size <= 19:
            categories["size_10_19"]["clusters"] += 1
            categories["size_10_19"]["genes"].update(genes)
        elif size >= 20:
            categories["size_20_plus"]["clusters"] += 1
            categories["size_20_plus"]["genes"].update(genes)

    # Convert gene sets to counts
    for cat in categories:
        categories[cat]["genes"] = len(categories[cat]["genes"])

    return categories

category_counts = count_categories_with_genes(cluster_genes)

print("Cluster and gene counts per category:")
for cat, info in category_counts.items():
    print(f"  {cat}: {info['clusters']} clusters, {info['genes']} unique genes")

# =========================
# Plot
# =========================
plt.figure(figsize=(7, 5))
plt.bar(sizes, counts)
plt.xlabel("Cluster size")
plt.ylabel("Frequency")
plt.title("() O. Brachyantha - Cluster Sizes Distribution")

plt.tight_layout()
plt.savefig(outfile, dpi=300)
plt.close()

print(f"Plot saved to: {outfile}")
