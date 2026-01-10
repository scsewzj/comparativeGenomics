#!/usr/bin/env python3
import sys
import re

def read_fasta_genes(fasta_file):
    """Read gene IDs from a FASTA file (header lines only)."""
    genes = set()
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                gene_id = line[1:].split()[0]  # take first word after >
                genes.add(gene_id)
    return genes

def extract_gene_positions(gff_file, filtered_genes, output_file):
    """Extract gene positions from GFF3 filtered by gene list."""
    with open(gff_file, 'r') as gff, open(output_file, 'w') as out:
        for line in gff:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            seqid, source, feature, start, end, score, strand, phase, attributes = parts
            if feature != "gene":
                continue
            # Extract gene ID robustly
            match = re.search(r"ID=([^;]+)", attributes) # <re.Match object; span=(0, 18), match='ID=gene:OB09G27280'>
            # print(match.group(1)) # gene:OB09G27220
            if match:
                gene_id = match.group(1).split(':')[1] + ".1"
                if gene_id in filtered_genes:
                    out.write(f"{gene_id}\t{seqid}\t{start}\t{end}\n")

def main():
    if len(sys.argv) < 2:
        print("Usage: ./makeBedFile.py genome.gff3 [pep.fasta]")
        sys.exit(1)

    gff_file = sys.argv[1]
    fasta_file = sys.argv[2] if len(sys.argv) > 2 else "pep_filtered.fa"
    output_file = "gene_pos.bed"

    print(f"Extracting gene positions from {gff_file} filtered by {fasta_file}...")

    filtered_genes = read_fasta_genes(fasta_file)
    # print(filtered_genes)
    extract_gene_positions(gff_file, filtered_genes, output_file)

    print(f"Done. Output saved in {output_file}")

if __name__ == "__main__":
    main()
