#!/usr/bin/env python3

import os
import re
from Bio import SeqIO
import pandas as pd

WORK_DIR = "/home/ngoc/comparativeGenomics/Oryza_sativa"
pepfile = os.path.join(WORK_DIR, "data/osativa.pep.lgiso.fa")
outfile = os.path.join(WORK_DIR, "data/protein_info.txt")

def extract_field(pattern, line, default=""):
    match = re.search(pattern, line)
    return match.group(1) if match else default

def parse_fasta_header(line):
    line = line[1:].strip()  # rm ">"
    
    fields = {
        'gene_id': extract_field(r"gene:([^ ]+)", line),
        'gene_symbol': extract_field(r"gene_symbol:([^ ]+)", line),
        'gene_biotype': extract_field(r"gene_biotype:([^ ]+)", line),
        'transcript': extract_field(r"transcript:([^ ]+)", line),
        'transcript_biotype': extract_field(r"transcript_biotype:([^ ]+)", line),
        'description': extract_field(r"description:(.+)$", line),
    }
    
    coord_match = re.search(r"chromosome:[^:]+:([^:]+):([^:]+):([^:]+):([^ ]+)", line)
    if coord_match:
        chrom, start, end, strand = coord_match.groups()
        fields.update({
            'chromosome': chrom,
            'start': start,
            'end': end,
            'strand': strand,
            'gene_length': str(abs(int(end) - int(start)))
        })
    else:
        fields.update({'chromosome': '', 'start': '', 'end': '', 'strand': '', 'gene_length': ''})
    
    return fields


protein_data = []

for record in SeqIO.parse(pepfile, "fasta"):
    info = parse_fasta_header(record.description)
    
    info['protein_id'] = record.id
    info['protein_length'] = len(str(record.seq).replace('*', ''))
    
    protein_data.append(info)


df = pd.DataFrame(protein_data)

column_order = ['protein_id', 'gene_id', 'gene_symbol', 'gene_biotype', 
                'transcript', 'transcript_biotype', 'chromosome', 'start', 
                'end', 'gene_length', 'protein_length', 'strand', 'description']

df[column_order].to_csv(outfile, sep='\t', index=False)

print(f"Processed {len(df)} proteins")
print(f"Saved to: {outfile}")
print(df[column_order].head())