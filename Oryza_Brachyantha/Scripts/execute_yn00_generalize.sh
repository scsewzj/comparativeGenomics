#!/bin/bash

usage() {
  cat <<EOF
Usage:
  $(basename "$0") INPUT_DIR YN00_CTL_TEMPLATE

Arguments:
  INPUT_DIR            Directory containing subdirectories with *.ali.phy files
                       (e.g. ks_filtered_input)

  YN00_CTL_TEMPLATE    yn00 control file template
                       (e.g. yn00.ctl_master)

Example:
  $(basename "$0") ks_filtered_input yn00.ctl_master
EOF
}

# -----------------------------
# Check arguments
# -----------------------------
if [[ $# -ne 2 ]]; then
  echo "Error: 2 arguments required."
  usage
  exit 1
fi

INPUT_DIR="$1"
CTL_TEMPLATE="$2"

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Error: INPUT_DIR is not a directory: $INPUT_DIR" >&2
  exit 1
fi

if [[ ! -f "$CTL_TEMPLATE" ]]; then
  echo "Error: Control template not found: $CTL_TEMPLATE" >&2
  exit 1
fi

# -----------------------------
# Main loop
# -----------------------------
for dir in "$INPUT_DIR"/*; do
    # Skip non-directories
    [[ -d "$dir" ]] || continue

    for d in "$dir"/*.ali.phy; do
        # Skip if no .ali.phy files exist
        [[ -e "$d" ]] || continue

        a=$(basename "$d")
        pair="${a%.ali.phy}"

        awk -v file="$d" -v out="$dir/yn" '
            {
                gsub("XXXXX", file)
                gsub("outfile *= *yn", "outfile = " out)
                print
            }
        ' "$CTL_TEMPLATE" > "$dir/yn00.ctl"

        yn00 "$dir/yn00.ctl"

        mv -i 2YN.dN 2YN.dS 2YN.t "$dir"
        mv -i rst rst1 rub "$dir"
    done
done
