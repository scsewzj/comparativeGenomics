#!/usr/bin/env bash

echo "==============================================="
echo "Minh Ngoc VU - M2 GENIOMHE Comparative genomics"
echo "==============================================="
echo " "

WORK_DIR="/home/ngoc/comparativeGenomics/Oryza_sativa/paml-workflow/"


echo "===== $(date): Running PAML YN00 for the low-stringent dataset..."

CONTROL_DIR="$WORK_DIR/low/controls"
OUT_DIR="$WORK_DIR/low/yn00_out"

mkdir -p "$OUT_DIR"

echo "Found $(find "$CONTROL_DIR" -name "*.ctl" | wc -l | tr -d ' ') control files."

find "$CONTROL_DIR" -name "*.ctl" | \
  parallel -j 16 --bar --eta \
  '
  tmpdir=$(mktemp -d)
  cp {} "$tmpdir/yn00.ctl"
  (cd "$tmpdir" && yn00 > /dev/null 2>&1)
  rm -rf "$tmpdir"
  '

# remove unnecessary result files
# rm 2YN* rst* rub

YN00_COUNT=$(find "$OUT_DIR" -name "*.out" | wc -l | tr -d ' ')
echo "===== $(date): PAML YN00 completed for $YN00_COUNT pairs in the low-stringent dataset."


echo "===== $(date): Running PAML YN00 for the high-stringent dataset..."

CONTROL_DIR="$WORK_DIR/high/controls"
OUT_DIR="$WORK_DIR/high/yn00_out"

mkdir -p "$OUT_DIR"

echo "Found $(find "$CONTROL_DIR" -name "*.ctl" | wc -l | tr -d ' ') control files."

find "$CONTROL_DIR" -name "*.ctl" | \
  parallel -j 16 --bar --eta \
  '
  tmpdir=$(mktemp -d)
  cp {} "$tmpdir/yn00.ctl"
  (cd "$tmpdir" && yn00 > /dev/null 2>&1)
  rm -rf "$tmpdir"
  '

YN00_COUNT=$(find "$OUT_DIR" -name "*.out" | wc -l | tr -d ' ')
echo "===== $(date): PAML YN00 completed for $YN00_COUNT pairs for the high-stringent dataset."