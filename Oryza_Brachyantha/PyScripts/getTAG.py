#!/usr/bin/env python3

import sys
import os
import argparse

import argparse

parser = argparse.ArgumentParser(
    description="Identify Tandemly Arrayed Genes (TAGs)"
)

parser.add_argument(
    "-b", "--blast",
    default="clusters/strict.tsv",
    help="BLAST tabular file (default: clusters/strict.tsv)"
)

parser.add_argument(
    "-g", "--bed",
    default="gene_pos.bed",
    help="BED file sorted by gene start position (default: gene_pos.bed)"
)

parser.add_argument(
    "-m", "--mcl",
    default="clusters/strict_mcl.tabular",
    help="MCL cluster file (default: clusters/strict_mcl.tabular)"
)

parser.add_argument(
    "-s", "--max_spacer",
    type=int,
    default=10,
    help="Maximum number of spacer genes allowed (default: 10)"
)

parser.add_argument(
    "-o", "--output",
    default="TAGs.tsv",
    help="Output TAG file (default: TAGs.tsv)"
)

args = parser.parse_args()

BLAST_FILE   = args.blast
BED_FILE     = args.bed
CLUSTER_FILE = args.mcl
MAX_SPACERS  = args.max_spacer
OUTPUT_FILE  = args.output

############################
# 1. LOAD BED FILE
############################
bed_genes = []
gene_index = {}

with open(BED_FILE) as bed:
    for line in bed:
        if line.strip() == "" or line.startswith("#"):
            continue
        fields = line.strip().split()
        gene = fields[0]
        gene_index[gene] = len(bed_genes)
        bed_genes.append(gene)

############################
# 2. LOAD CLUSTERS (list of sets)
############################
clusters = []

with open(CLUSTER_FILE) as cl:
    for line in cl:
        if line.strip() == "" or line.startswith("#"):
            continue
        genes = line.strip().split()
        clusters.append(set(genes))

############################
# 3. FAST LOOKUP: gene → clusters
############################
gene_to_clusters = {}

for i, cluster in enumerate(clusters):
    for gene in cluster:
        gene_to_clusters.setdefault(gene, set()).add(i)

############################
# 4. PROCESS BLAST & IDENTIFY TAGs
############################
tags = {}  
# key: (gene1, gene2)
# value: number of spacer genes

with open(BLAST_FILE) as blast:
    for line in blast:
        if line.strip() == "" or line.startswith("#"):
            continue

        g1, g2 = line.strip().split()[:2]

        # both genes must be in BED and clusters
        if g1 not in gene_index or g2 not in gene_index:
            continue
        
         # genes must exist in clusters
        if g1 not in gene_to_clusters or g2 not in gene_to_clusters:
            continue

        # same cluster?
        if gene_to_clusters[g1].isdisjoint(gene_to_clusters[g2]):
            continue

        # distance in BED list
        idx1 = gene_index[g1]
        idx2 = gene_index[g2]

        spacer_count = abs(idx1 - idx2) - 1

        if spacer_count <= MAX_SPACERS:
            pair = tuple(sorted([g1, g2]))
            tags[pair] = spacer_count

############################
# 4. OUTPUT
############################

# print("gene1\tgene2\tspacer_genes")

# for (g1, g2), spacers in tags.items():
#     print(f"{g1}\t{g2}\t{spacers}")

with open(OUTPUT_FILE, "w") as out:
    out.write("gene1\tgene2\tspacer_genes\n")
    for (g1, g2), spacers in sorted(tags.items()):
        out.write(f"{g1}\t{g2}\t{spacers}\n")

print(f"The results have been written at {OUTPUT_FILE}.")

