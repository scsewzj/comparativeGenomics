#!/bin/bash


# -----------------------------
# Usage
# -----------------------------
usage() {
  cat <<EOF

  # Usage:
  #   ./makeOverRep.sh [input_file] [output_file]
  #
  # Arguments:
  #   input_file   Optional. Default: TAGS/TAGs_0_spacers.txt
  #   output_file  Optional. Default: TAGS/TAGs_OB2GO_Overrep.txt
  #
  #   -h         Show this help message
  #
  # Description:
  #   Processes a list of TAGs, joins with OB_to_Os_IDs and Os_Prot2Gen mappings,
  #   counts occurrences, filters by minimum count (>=10), and outputs
  #   the list of TAGs for overrepresentation analysis.

  Example:
    # Using default input/output
    $./makeOverRep.sh

    # Using custom input and output files
    $./makeOverRep.sh my_custom_spacers.txt my_custom_OB2GO.txt
EOF
}

# -----------------------------
# Parse arguments
# -----------------------------

while getopts ":h" opt; do
  case "${opt}" in
    h) usage; exit 0 ;;
    *) usage; exit 1 ;;
  esac
done

# Default input/output
DEFAULT_INPUT="TAGS/TAGs_0_spacers.txt"

# Automatically create output file name by adding suffix
# Extract the base name without directory and extension
INPUT="${1:-$DEFAULT_INPUT}"

BASENAME="$(basename "$INPUT" .txt)"
OUTPUT="${2:-TAGS/${BASENAME}_OB2GO_Overrep.txt}"

# Temporary files
TMP1="TAGS/TAGs_0_spacers_IDs.txt"
TMP2="TAGS/temp0.txt"
TMP3="TAGS/OB2GO_temp.txt"
TMP4="TAGS/OB2GO_temp.clean.txt"
TMP5="TAGS/OB2GO_temp.uniq.txt"

# Step 1: Extract first PEP ID and TAGs Array Size, then sort
awk '{print $1, NF}' "$INPUT" | sort > "$TMP1"

# Step 2: Join with homologous O. Sativa PEP IDs.txt
join -1 1 -2 1 OB_to_Os_IDs.txt "$TMP1" > "$TMP2"

# Step 3: Sort O.B and join with Os_Prot2Gen.txt
sort -k 2 "$TMP2" | join -1 2 -2 1 - Data/Os_Prot2Gen.txt > "$TMP3"

# Cleanup temporary file
rm "$TMP2"

# Step 4: Keep O.S and TAGs Array Size
awk '{print $4 "\t" $3}' "$TMP3" > "$TMP4"

# Cleanup temporary file
rm "$TMP3"

# Step 5: Count unique occurrences, sort, and reorder
sort -k 1 "$TMP4" | uniq -c | sort -nr -k 1 | awk '{print $2 "\t" $3}' | sort -k 2 -nr > "$TMP5"

# Cleanup temporary file
rm "$TMP4"

# Step 6: Keep only "important" TAGs Arrays (Size >= 10)
awk '{if($2 >= 10){print $1}}' "$TMP5" > "$OUTPUT"

# Cleanup temporary file
rm "$TMP5"

echo "Done! Output saved to $OUTPUT"
