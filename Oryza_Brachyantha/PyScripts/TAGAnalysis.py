#!/usr/bin/env python3

import os
import argparse
from collections import defaultdict, deque
from itertools import combinations

############################
# ARGUMENTS
############################

parser = argparse.ArgumentParser(
    description="Identify TAG arrays and compute statistics"
)

parser.add_argument("-b", "--blast", required=True,
                    help="BLAST tabular file")
parser.add_argument("-g", "--bed", required=True,
                    help="BED file sorted by gene order")
parser.add_argument("-m", "--mcl", required=True,
                    help="MCL cluster file")

args = parser.parse_args()

BLAST_FILE = args.blast
BED_FILE = args.bed
CLUSTER_FILE = args.mcl

SPACER_RUNS = [(0,0), (1,5), (6,10)]

TAG0_FILE = "TAGs_0_spacers.txt"
NONTAG_FILE = "non_TAGs.txt"
BASE_DIR = "TAGS"

# =========================
# Prepare output directory path
# =========================
os.makedirs(BASE_DIR, exist_ok=True)

############################
# 1. LOAD BED FILE
############################

gene_index = {}
gene_chrom = {}
chrom_genes = defaultdict(list)

with open(BED_FILE) as bed:
    for line in bed:
        if not line.strip() or line.startswith("#"):
            continue
        f = line.strip().split()
        gene, chrom = f[0], f[1]
        gene_chrom[gene] = chrom
        gene_index[gene] = len(chrom_genes[chrom])
        chrom_genes[chrom].append(gene)

############################
# 2. LOAD MCL CLUSTERS
############################

gene_to_clusters = defaultdict(set)

with open(CLUSTER_FILE) as cl:
    for i, line in enumerate(cl):
        if not line.strip() or line.startswith("#"):
            continue
        genes = line.strip().split()
        for g in genes:
            gene_to_clusters[g].add(i)

############################
# 3. LOAD BLAST PAIRS ONCE
############################

blast_pairs = []

with open(BLAST_FILE) as blast:
    for line in blast:
        if not line.strip() or line.startswith("#"):
            continue
        g1, g2 = line.strip().split()[:2]
        blast_pairs.append((g1, g2))

############################
# 4. MAIN LOOP OVER SPACER VALUES
############################
genes_in_any_tag = set()


print(
    "Spacer_range\tTAG_pairs\tTAG_arrays\tMax_array\t"
    "Arrays_size2\tArrays_3_5\tArrays_6_9\tArrays_10plus"
)

for MAX_SPACERS in SPACER_RUNS:

    ############################
    # 4A. FIND TAG PAIRS
    ############################

    tag_pairs = set()

    for g1, g2 in blast_pairs:

        if g1 not in gene_index or g2 not in gene_index:
            continue
        if gene_chrom[g1] != gene_chrom[g2]:
            continue
        if gene_to_clusters[g1].isdisjoint(gene_to_clusters[g2]):
            continue

        idx1, idx2 = gene_index[g1], gene_index[g2]
        spacers = abs(idx1 - idx2) - 1

        if spacers >= MAX_SPACERS[0] and spacers <= MAX_SPACERS[1]:
            pair = tuple(sorted((g1, g2)))
            tag_pairs.add(pair)
        
            # Track TAG genes globally
            genes_in_any_tag.update(pair)

        

    ############################
    # 4B. BUILD TAG ARRAYS (single-linkage)
    ############################

    chrom_graph = defaultdict(lambda: defaultdict(set))

    for g1, g2 in tag_pairs:
        chrom = gene_chrom[g1]
        chrom_graph[chrom][g1].add(g2)
        chrom_graph[chrom][g2].add(g1)

    tag_arrays = []

    for chrom, graph in chrom_graph.items():
        visited = set()
        for gene in graph:
            if gene in visited:
                continue
            queue = deque([gene])
            comp = set()
            while queue:
                g = queue.popleft()
                if g in visited:
                    continue
                visited.add(g)
                comp.add(g)
                queue.extend(graph[g] - visited)
            if len(comp) >= 2:
                tag_arrays.append(comp)
                
    ############################
    # 4Bbis. WRITE TAGs WITH 0 SPACERS
    ############################

    if MAX_SPACERS == (0, 0):
        OUT_FILE0 = os.path.join(BASE_DIR, TAG0_FILE)
        with open(OUT_FILE0, "w") as out:
            for arr in tag_arrays:
                out.write("\t".join(sorted(arr)) + "\n")
    print(f"TAGs-0 saved to: {OUT_FILE0}")


    ############################
    # 4C. STATISTICS
    ############################

    num_arrays = len(tag_arrays)
    num_pairs = sum(len(a) * (len(a) - 1) // 2 for a in tag_arrays)
    max_array = max((len(a) for a in tag_arrays), default=0)

    size2 = size3_5 = size6_9 = size10p = 0

    for arr in tag_arrays:
        n = len(arr)
        if n == 2:
            size2 += 1
        elif 3 <= n <= 5:
            size3_5 += 1
        elif 6 <= n <= 9:
            size6_9 += 1
        else:
            size10p += 1

    ############################
    # 4D. PRINT SUMMARY LINE
    ############################

    print(
        f"{MAX_SPACERS}\t{num_pairs}\t{num_arrays}\t{max_array}\t"
        f"{size2}\t{size3_5}\t{size6_9}\t{size10p}"
    )
    
############################
# 5. WRITE NON-TAG GENES
############################

all_genes = set(gene_index.keys())
non_tag_genes = sorted(all_genes - genes_in_any_tag)

OUT_FILE = os.path.join(BASE_DIR, NONTAG_FILE)
with open(OUT_FILE, "w") as out:
    for g in non_tag_genes:
        out.write(g + "\n")
print(f"Non-TAGs saved to: {OUT_FILE}")

