#!/usr/bin/env bash

out="KaKs_all.tsv"
echo -e "Cluster\tPair\tKa\tKs\tomega" > "$out"

# Loop over all clusters
for cluster in cluster_3526; do
    [ -d "$cluster" ] || continue

    # Loop over all gene pairs
    for pair_dir in "$cluster"/*; do
        [ -d "$pair_dir" ] || continue

        yn_file="$pair_dir/yn"
        [ -f "$yn_file" ] || continue

        # Get pair name
        pair=$(basename "$pair_dir")

        # Extract last line after 'seq. seq.' header, ignoring empty lines
        vals=$(awk 'NF && /seq\. seq\./{getline; print; exit}' "$yn_file")

        if [ -z "$vals" ]; then
            # Skip if no data found
            continue
        fi

        # Extract Ka (dN), Ks (dS), omega
        omega=$(echo "$vals" | awk '{print $7}')
        Ka=$(echo "$vals" | awk '{print $8}')
        Ks=$(echo "$vals" | awk '{print $10}')

        echo -e "${cluster}\t${pair}\t${Ka}\t${Ks}\t${omega}" >> "$out"
    done
done

echo "Ka/Ks/omega extraction completed. Output: $out"
