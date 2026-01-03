#!/usr/bin/env python3

import os
import pandas as pd

DATA_DIR = "/home/ngoc/comparativeGenomics/Oryza_sativa/data"

# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
blast_cols = ["qseqid","sseqid","pident","length","mismatch","gapopen",
              "qstart","qend","sstart","send","evalue","bitscore"]

dt = pd.read_csv(os.path.join(DATA_DIR, "osativa_blastp_results_dedup.txt"), 
                 sep="\t", names=blast_cols, header=None)
length_dt = pd.read_csv(os.path.join(DATA_DIR, "protein_info.txt"), 
                      sep="\t")
length_dt = length_dt.sort_values(["protein_id", "transcript"])

length_dict = length_dt.set_index("transcript")["protein_length"].to_dict()

dt["query_length"] = dt["qseqid"].map(length_dict)
dt["subject_length"] = dt["sseqid"].map(length_dict)

print("Number of missing query length: ", sum(dt["query_length"].isna()))
print("Number of missing subject length: ", sum(dt["subject_length"].isna()))

dt.to_csv(os.path.join(DATA_DIR, "osativa_blastp_results_dedup_lengths.txt"), sep="\t", index=False, header=False)

