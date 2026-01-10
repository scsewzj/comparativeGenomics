#!/usr/bin/env python3

import sys
import os
from itertools import combinations

# --------- INPUTS ---------
mcl_file = sys.argv[1]     # e.g. clusters.mcl
out_dir = sys.argv[2]      # e.g. cluster_pairs

os.makedirs(out_dir, exist_ok=True)

# --------- PROCESS MCL FILE ---------
with open(mcl_file) as f:
    for idx, line in enumerate(f, start=1):
        line = line.strip()
        if not line:
            continue

        genes = line.split("\t")
        if len(genes) < 2:
            continue  # skip singletons

        cluster_name = f"cluster_{idx}"
        out_file = os.path.join(out_dir, f"{cluster_name}.pairs.txt")

        with open(out_file, "w") as out:
            for g1, g2 in combinations(genes, 2):
                out.write(f"{g1}\t{g2}\n")

print(f"✅ Pair files generated in: {out_dir}")
