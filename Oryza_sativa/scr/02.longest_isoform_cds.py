#!/usr/bin/env python3

from Bio import SeqIO
import sys

pep_fa = sys.argv[1]   # filtered PEP fasta
cds_fa = sys.argv[2]   # original CDS fasta

# 1. collect transcript IDs from PEP
kept_transcripts = set()

for rec in SeqIO.parse(pep_fa, "fasta"):
    # >Os08t0254300-00 pep ...
    tx = rec.id
    kept_transcripts.add(tx)

# 2. filter CDS
kept_cds = []

for rec in SeqIO.parse(cds_fa, "fasta"):
    # >Os08t0254300-00 cds ...
    if rec.id in kept_transcripts:
        kept_cds.append(rec)

# 3. write output
SeqIO.write(kept_cds, sys.stdout, "fasta")