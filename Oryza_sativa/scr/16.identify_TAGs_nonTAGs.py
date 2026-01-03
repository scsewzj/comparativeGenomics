#!/usr/bin/env python3

print("===============================================")
print("Minh Ngoc VU - M2 GENIOMHE Comparative genomics")
print("===============================================")
print(" ")

from collections import defaultdict
from datetime import datetime
import sys
import os

if len(sys.argv) < 3:
    print("Usage: python script.py <MCL_clusters_file> <dataset>")
    sys.exit(1)

WORK_DIR = "/home/ngoc/comparativeGenomics/Oryza_sativa/"
OUT_DIR = os.path.join(WORK_DIR, "mcl/TAGs")
MCL_FILE = os.path.join(WORK_DIR, "mcl/output", sys.argv[1])
MCL_DATASET = sys.argv[2]
GFF_FILE = os.path.join(WORK_DIR, "data/Oryza_sativa.IRGSP-1.0.41.gff3")

SPACER_CLASSES = {"TAG_0": (0, 0),
                  "TAG_1_5": (1, 5),
                  "TAG_6_10": (6, 10)}

# ----- parse the annotation file
print(f"===== {datetime.now().strftime('%c')}: Parsing the annotation file {GFF_FILE}...")
gene_pos = {}
genes_by_chr = defaultdict(list)

with open(GFF_FILE) as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.rstrip().split("\t")
        if len(fields) < 9 or fields[2] != "mRNA":
            continue

        chrom = fields[0]
        start = int(fields[3])
        attrs = fields[8]

        gene_id = None
        for x in attrs.split(";"):
            if x.startswith("transcript_id="):
                gene_id = x.replace("transcript_id=", "").strip()
                break

        if gene_id:
            gene_pos[gene_id] = (chrom, start)
            genes_by_chr[chrom].append((start, gene_id))


# sort genes by position and create index lookup
gene_index_by_chr = {}
for chrom in genes_by_chr:
    genes_by_chr[chrom].sort()
    gene_index_by_chr[chrom] = {g: i for i, (_, g) in enumerate(genes_by_chr[chrom])}


# ----- read MCL clusters
print(f"===== {datetime.now().strftime('%c')}: Reading MCL clusters from {MCL_FILE}...")
clusters = []

with open(MCL_FILE) as f:
    for line in f:
        fam = line.strip().split("\t")
        if fam:
            clusters.append(set(fam))
print(f"===== {datetime.now().strftime('%c')}: Completed reading {len(clusters)} clusters.")



# ----- spacer count between consecutive family members
def spacer_count(g1, g2, family, chrom):
    """Count non-family genes between g1 and g2"""
    i1 = gene_index_by_chr[chrom][g1]
    i2 = gene_index_by_chr[chrom][g2]
    
    lo, hi = sorted([i1, i2])
    ordered = [g for _, g in genes_by_chr[chrom]]
    
    return sum(1 for g in ordered[lo+1:hi] if g not in family)


# TAGs: pairs and arrays
print(f"===== {datetime.now().strftime('%c')}: Detecting TAG pairs and arrays...")

pair_TAGs = defaultdict(list)
array_TAGs = defaultdict(list)
non_TAG_pairs = []  # Add this to track non-TAG pairs

for fam_id, family in enumerate(clusters, 1):
    by_chr = defaultdict(list)

    # group family genes by chromosome
    for g in family:
        if g in gene_pos:
            chrom, start = gene_pos[g]
            by_chr[chrom].append((start, g))

    # process chromosomes
    for chrom, genes in by_chr.items():
        genes.sort()
        ordered_genes = [g for _, g in genes]
        
        # check each consecutive pair
        for i in range(1, len(ordered_genes)):
            prev = ordered_genes[i-1]
            curr = ordered_genes[i]
            spacers = spacer_count(prev, curr, family, chrom)
            
            # check if this pair is a TAG in ANY spacer class
            is_TAG = False
            for cls, (lo, hi) in SPACER_CLASSES.items():
                if lo <= spacers <= hi:
                    is_TAG = True
                    pair_TAGs[cls].append((fam_id, chrom, prev, curr, spacers))
            
            # if not a TAG in any class, record as non-TAG
            if not is_TAG:
                non_TAG_pairs.append((fam_id, chrom, prev, curr, spacers))
        
        # build arrays for each spacer class
        for cls, (lo, hi) in SPACER_CLASSES.items():
            current_array = [ordered_genes[0]]
            current_spacers = []

            for i in range(1, len(ordered_genes)):
                prev = ordered_genes[i-1]
                curr = ordered_genes[i]
                spacers = spacer_count(prev, curr, family, chrom)
                
                if lo <= spacers <= hi:
                    current_array.append(curr)
                    current_spacers.append(spacers)
                else:
                    if len(current_array) >= 2:
                        array_TAGs[cls].append((fam_id, chrom, current_array.copy(), current_spacers.copy()))
                    current_array = [curr]
                    current_spacers = []

            # last array
            if len(current_array) >= 2:
                array_TAGs[cls].append((fam_id, chrom, current_array.copy(), current_spacers.copy()))


for cls in SPACER_CLASSES:
    with open(os.path.join(OUT_DIR, f"{MCL_DATASET}_{cls}_pairs.tsv"), "w") as f:
        f.write("family_id\tchr\tgene1\tgene2\tspacers\n")
        for r in pair_TAGs[cls]:
            f.write("\t".join(map(str, r)) + "\n")

    with open(os.path.join(OUT_DIR, f"{MCL_DATASET}_{cls}_arrays.tsv"), "w") as f:
        f.write("family_id\tchr\tarray_size\tgenes\tspacers\n")
        for fam_id, chrom, genes, spacers in array_TAGs[cls]:
            array_size = len(genes)
            f.write(f"{fam_id}\t{chrom}\t{array_size}\t{','.join(genes)}\t{','.join(map(str, spacers))}\n")

# non-TAG pairs
with open(os.path.join(OUT_DIR, f"{MCL_DATASET}_non_TAG_pairs.tsv"), "w") as f:
    f.write("family_id\tchr\tgene1\tgene2\tspacers\n")
    for r in non_TAG_pairs:
        f.write("\t".join(map(str, r)) + "\n")

print(f"===== {datetime.now().strftime('%c')}: TAG pairs and arrays detection complete, output saved in {OUT_DIR}.")
print(f"===== {datetime.now().strftime('%c')}: Found {len(non_TAG_pairs)} non-TAG pairs (spacers > 10).")