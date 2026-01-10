#!/usr/bin/env python3
import sys
from itertools import combinations
from Bio import SeqIO
import os

# ----- CONFIG -----
mcl_file = sys.argv[1]        # e.g., clusters.txt
cds_fasta = sys.argv[2]       # e.g., cds.fasta
prot_fasta = sys.argv[3]      # e.g., protein.fasta
out_dir = sys.argv[4]         # e.g., "pairs_output"
os.makedirs(out_dir, exist_ok=True)

# ----- LOAD CDS AND PROTEIN SEQUENCES -----
cds_dict = SeqIO.to_dict(SeqIO.parse(cds_fasta, "fasta"))
prot_dict = SeqIO.to_dict(SeqIO.parse(prot_fasta, "fasta"))

# ----- LOAD MCL CLUSTERS -----
clusters = {}
with open(mcl_file) as f:
    for idx, line in enumerate(f):
        line = line.strip()
        if not line:
            continue
        genes = line.split()
        cluster_id = f"cluster_{idx+1}"
        clusters[cluster_id] = genes

# ----- GENERATE PAIRS AND WRITE FASTA FILES IN PAIR DIRECTORIES -----
total_pairs = 0
for cluster_id, genes in clusters.items():
    if len(genes) < 2:
        continue  # skip singletons
    
    for gene1, gene2 in combinations(genes, 2):
        total_pairs += 1
        pair_dir_name = f"{gene1}_{gene2}"
        pair_dir_path = os.path.join(out_dir, pair_dir_name)
        os.makedirs(pair_dir_path, exist_ok=True)

        # Prepare CDS file
        cds_out_file = os.path.join(pair_dir_path, f"{pair_dir_name}_cds.fasta")
        cds_seqs = []
        if gene1 in cds_dict and gene2 in cds_dict:
            cds_seqs = [cds_dict[gene1], cds_dict[gene2]]
            SeqIO.write(cds_seqs, cds_out_file, "fasta")
        else:
            print(f"Warning: {gene1} or {gene2} not found in CDS FASTA")
        
        # Prepare Protein file
        prot_out_file = os.path.join(pair_dir_path, f"{pair_dir_name}_prot.fasta")
        prot_seqs = []
        if gene1 in prot_dict and gene2 in prot_dict:
            prot_seqs = [prot_dict[gene1], prot_dict[gene2]]
            SeqIO.write(prot_seqs, prot_out_file, "fasta")
        else:
            print(f"Warning: {gene1} or {gene2} not found in protein FASTA")

print(f"Prepared {total_pairs} gene pair directories with CDS and protein FASTA files in {out_dir}.")
