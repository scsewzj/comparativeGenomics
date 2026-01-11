import argparse
from Bio import SeqIO

# Set up argument parser
parser = argparse.ArgumentParser(description="Parse FASTA headers into structured TSV.")
parser.add_argument("input_fasta", help="Input FASTA file")
parser.add_argument("output_tsv", help="Output TSV file")
args = parser.parse_args()

# Open output file and write header
with open(args.output_tsv, "w") as out_f:
    out_f.write("\t".join([
        "protein_id", "gene_id", "gene_biotype", "transcript",
        "transcript_biotype", "chromosome", "start", "end",
        "gene_length", "protein_length", "strand"
    ]) + "\n")

    # Parse FASTA sequences
    for record in SeqIO.parse(args.input_fasta, "fasta"):
        header = record.description
        seq = record.seq

        # Initialize dictionary to hold info
        info = {
            "protein_id": record.id,
            "gene_id": "",
            "gene_biotype": "",
            "transcript": "",
            "transcript_biotype": "",
            "chromosome": "",
            "start": "",
            "end": "",
            "strand": ""
        }

        # Split header into key:value parts
        for part in header.split():
            if ":" in part:
                key, value = part.split(":", 1)
                if key == "gene":
                    info["gene_id"] = value
                elif key == "gene_biotype":
                    info["gene_biotype"] = value
                elif key == "transcript":
                    info["transcript"] = value
                elif key == "transcript_biotype":
                    info["transcript_biotype"] = value
                elif key == "chromosome":
                    info["chromosome"] = value
                elif key == "start":
                    info["start"] = int(value)
                elif key == "end":
                    info["end"] = int(value)
                elif key == "strand":
                    info["strand"] = value

        # Calculate gene length
        if info["start"] and info["end"]:
            gene_length = abs(info["end"] - info["start"]) + 1
        else:
            gene_length = ""

        # Protein length from sequence
        protein_length = len(seq)

        # Write line to output
        out_f.write("\t".join([
            info["protein_id"],
            info["gene_id"],
            info["gene_biotype"],
            info["transcript"],
            info["transcript_biotype"],
            info["chromosome"],
            str(info["start"]),
            str(info["end"]),
            str(gene_length),
            str(protein_length),
            info["strand"]
        ]) + "\n")

print(f"Output saved to {args.output_tsv}")
