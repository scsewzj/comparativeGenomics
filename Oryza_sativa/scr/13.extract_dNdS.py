#!/usr/bin/env python3

from datetime import datetime
import re
import os
import pandas as pd

print("===============================================")
print("Minh Ngoc VU - M2 GENIOMHE Comparative genomics")
print("===============================================")
print(" ")

FLOAT = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"

header_pattern = re.compile(
    r"seq\.\s+seq\.\s+S\s+N\s+t\s+kappa\s+omega\s+dN\s+\+-\s+SE\s+dS\s+\+-\s+SE"
)

data_pattern = re.compile(
    rf"\s*(\d+)\s+(\d+)\s+"              # seq1 seq2
    rf"({FLOAT})\s+({FLOAT})\s+"         # S N
    rf"({FLOAT})\s+({FLOAT})\s+({FLOAT})\s+"  # t kappa omega
    rf"({FLOAT})\s*\+\-\s*({FLOAT})\s+"       # dN dN_SE
    rf"({FLOAT})\s*\+\-\s*({FLOAT})"          # dS dS_SE
)

def extract_dnds(yn00_dir, output_csv):
    """
    Parse all .out files in a YN00 directory and save as a CSV.
    """
    rows = []

    for fname in os.listdir(yn00_dir):
        if not fname.endswith(".out"):
            continue

        file_path = os.path.join(yn00_dir, fname)
        s1 = fname.split("_")[0]
        s2 = fname.split("_")[1].split(".")[0]

        with open(file_path, "r") as f:
            look_for_data_line = False

            for line in f:
                if header_pattern.search(line):
                    look_for_data_line = True
                    continue

                if look_for_data_line:
                    if line.strip() == "":
                        continue  # skip blank lines

                    match = data_pattern.match(line)
                    if match:
                        seq1, seq2, S, N, t, kappa, omega, dN, dN_se, dS, dS_se = match.groups()

                        rows.append({
                            "hit": fname.split(".")[0],
                            "seq1": s1,
                            "seq2": s2,
                            "S": float(S),
                            "N": float(N),
                            "t": float(t),
                            "kappa": float(kappa),
                            "omega": float(omega),
                            "dN": float(dN),
                            "dN_SE": float(dN_se),
                            "dS": float(dS),
                            "dS_SE": float(dS_se),
                        })
                    else:
                        print(f"WARNING: Could not parse line in {fname}:")
                        print(line.strip())

                    look_for_data_line = False

    df = pd.DataFrame(rows)
    df.to_csv(output_csv, index=False)
    print(f"===== {datetime.now().strftime('%c')}: Extracted dN/dS for {df.shape[0]} hits from {yn00_dir}")


PAML_DIR = "/home/ngoc/comparativeGenomics/Oryza_sativa/paml-workflow"
YN00_L = os.path.join(PAML_DIR, "low/yn00_out")
YN00_H = os.path.join(PAML_DIR, "high/yn00_out")


print(f"===== {datetime.now().strftime('%c')}: Extracting dN dS (Yang and Nielsen (2000)) from YN00 output of the low stringent dataset...")
extract_dnds(YN00_L, os.path.join(PAML_DIR, "low/dNdS.csv"))

print(f"===== {datetime.now().strftime('%c')}: Extracting dN dS (Yang and Nielsen (2000)) from YN00 output of the high stringent dataset...")
extract_dnds(YN00_H, os.path.join(PAML_DIR, "high/dNdS.csv"))
