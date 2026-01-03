#!/usr/bin/env python3

import os
from datetime import datetime

print("===============================================")
print("Minh Ngoc VU - M2 GENIOMHE Comparative genomics")
print("===============================================")
print(" ")
print(f"===== {datetime.now().strftime('%c')}: Extracting AA and CDS for each hit of BLASTP output")

WORK_DIR = "/home/ngoc/comparativeGenomics/Oryza_sativa"

BLAST_OUTPUT_H = os.path.join(WORK_DIR, "data/filtered/osativa_H.txt")
BLAST_OUTPUT_L = os.path.join(WORK_DIR, "data/filtered/osativa_L.txt")

CDS = os.path.join(WORK_DIR, "data/osativa.cds.lgiso.fa")
PEP = os.path.join(WORK_DIR, "data/osativa.pep.lgiso.fa")

OUT_AA_H  = os.path.join(WORK_DIR, "paml-workflow/high/aa")
OUT_CDS_H = os.path.join(WORK_DIR, "paml-workflow/high/cds")
os.makedirs(OUT_AA_H, exist_ok=True)
os.makedirs(OUT_CDS_H, exist_ok=True)

OUT_AA_L  = os.path.join(WORK_DIR, "paml-workflow/low/aa")
OUT_CDS_L = os.path.join(WORK_DIR, "paml-workflow/low/cds")
os.makedirs(OUT_AA_L, exist_ok=True)
os.makedirs(OUT_CDS_L, exist_ok=True)


def fasta2dict(filepath):
    """
    Reads a FASTA file into a dictionary {ID: Sequence}
    """
    seqs = {}
    current_id = None
    current_seq = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            
            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = "".join(current_seq)
                
                # >Os08t0254300-00 cds ... -> Os08t0254300-00
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
                
        if current_id:
            seqs[current_id] = "".join(current_seq)
            
    return seqs


def extract_hits(blast_file, cds_db, pep_db, out_aa_dir, out_cds_dir):
    """
    Extract AA and CDS sequences for each BLAST hit
    """
    with open(blast_file, 'r') as f:
        for line in f:
            if not line.strip():
                continue

            cols = line.rstrip().split('\t')
            q_id = cols[0]
            s_id = cols[1]

            hit = f"{q_id}_{s_id}"

            if (q_id in cds_db and s_id in cds_db and
                q_id in pep_db and s_id in pep_db):

                with open(os.path.join(out_aa_dir, f"{hit}.pep.fa"), 'w') as out_aa:
                    out_aa.write(f">{q_id}\n{pep_db[q_id]}\n")
                    out_aa.write(f">{s_id}\n{pep_db[s_id]}\n")

                with open(os.path.join(out_cds_dir, f"{hit}.cds.fa"), 'w') as out_cds:
                    out_cds.write(f">{q_id}\n{cds_db[q_id]}\n")
                    out_cds.write(f">{s_id}\n{cds_db[s_id]}\n")

            else:
                print(f"CDS or AA sequence missing for query {q_id} and subject {s_id}")



def main():
    cds_db = fasta2dict(CDS)
    pep_db = fasta2dict(PEP)

    extract_hits(BLAST_OUTPUT_H,
                 cds_db,
                 pep_db,
                 OUT_AA_H,
                 OUT_CDS_H)

    extract_hits(BLAST_OUTPUT_L,
                 cds_db,
                 pep_db,
                 OUT_AA_L,
                 OUT_CDS_L)

if __name__ == "__main__":
    main()
    print(f"===== {datetime.now().strftime('%c')}:  Extracting completed.")