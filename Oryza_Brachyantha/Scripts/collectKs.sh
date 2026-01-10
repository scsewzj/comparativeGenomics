#!/bin/bash

output="Ks_all_clusters.tsv"
echo -e "Cluster\tPair\tKs" > $output

for cluster in cluster_*; do
    [ -d "$cluster" ] || continue

    for pairdir in "$cluster"/*; do
        [ -d "$pairdir" ] || continue

        # detect yn00 output file
        f=""
        for candidate in yn00.out yn yn00_output.txt; do
            if [[ -f "$pairdir/$candidate" ]]; then
                f="$pairdir/$candidate"
            fi
        done
        
        [ -z "$f" ] && continue

        pair=$(basename "$pairdir")

        # Extract Ks: look for line containing "t=" (Ks field in YN00 summary)
        ks=$(grep -E "t\s*=" "$f" | head -1 | sed -E 's/.*t *= *([0-9.]+).*/\1/')

        # skip if extraction failed
        [[ -z "$ks" ]] && continue

        echo -e "${cluster}\t${pair}\t${ks}" >> $output
    done
done

echo "Saved Ks table to $output"
