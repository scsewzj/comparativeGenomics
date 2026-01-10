#!/usr/bin/env bash

# Default input and output files
INPUT="${1:-blast_with_lengths.tsv}"
OUTPUT="${2:-longest_pairs.tsv}"

awk '
{
    if ($1 < $2)
        pair = $1 FS $2
    else
        pair = $2 FS $1

    if (pair in maxlen) {
        if ($4 > maxlen[pair]) {
            maxlen[pair] = $4
            best[pair] = $0
        }
    } else {
        maxlen[pair] = $4
        best[pair] = $0
    }
}
END {
    for (p in best)
        print best[p]
}
' "$INPUT" > "$OUTPUT"
