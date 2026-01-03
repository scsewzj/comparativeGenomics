#!/usr/bin/env python3

from Bio import SeqIO
import sys

longest_per_gene = {}

for rec in SeqIO.parse(sys.argv[1], "fasta"):
    # Os08t0254300-00 → Os08t0254300
    gene = rec.id.split("-")[0]

    seq_len = len(rec)
    if gene not in longest_per_gene or seq_len > len(longest_per_gene[gene]):
        longest_per_gene[gene] = rec

# write results
SeqIO.write(longest_per_gene.values(), sys.stdout, "fasta")