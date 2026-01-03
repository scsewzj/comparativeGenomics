#!/usr/bin/env bash

echo "==============================================="
echo "Minh Ngoc VU - M2 GENIOMHE Comparative genomics"
echo "==============================================="
echo " "

WORK_DIR="/home/ngoc/comparativeGenomics/Oryza_sativa/paml-workflow"

echo "===== $(date): Creating control files for each pair (low-stringent dataset) to run PAML..."

PAL2NAL_DIR="$WORK_DIR/low/pal2nal_out"
CONTROL_MASTER="$WORK_DIR/low/yncontrol.txt"
CTL_OUT_DIR="$WORK_DIR/low/controls"
YN00_OUT_DIR="$WORK_DIR/low/yn00_out"

mkdir -p "$CTL_OUT_DIR"


find "$PAL2NAL_DIR" -name "*.phy" -type f | \
  parallel -j 16 --bar --eta \
  'base=$(basename {} .phy); sed "s|XXXXX|{}|g; s|yn|'"$YN00_OUT_DIR"'/${base}.out|g" '"$CONTROL_MASTER"' > '"$CTL_OUT_DIR"'/${base}.ctl'

CTL_COUNT=$(find "$CTL_OUT_DIR" -name "*.ctl" -type f 2>/dev/null | wc -l | tr -d ' ')
echo "===== $(date): Creating controls completed (low-stringent dataset) for $CTL_COUNT files."


echo "===== $(date): Creating control files for each pair (high-stringent dataset) to run PAML..."

PAL2NAL_DIR="$WORK_DIR/high/pal2nal_out"
CONTROL_MASTER="$WORK_DIR/high/yncontrol.txt"
CTL_OUT_DIR="$WORK_DIR/high/controls"
YN00_OUT_DIR="$WORK_DIR/high/yn00_out"

mkdir -p "$CTL_OUT_DIR"

find "$PAL2NAL_DIR" -name "*.phy" -type f | \
  parallel -j 16 --bar --eta \
  'base=$(basename {} .phy); sed "s|XXXXX|{}|g; s|yn|'"$YN00_OUT_DIR"'/${base}.out|g" '"$CONTROL_MASTER"' > '"$CTL_OUT_DIR"'/${base}.ctl'


CTL_COUNT=$(find "$CTL_OUT_DIR" -name "*.ctl" -type f 2>/dev/null | wc -l | tr -d ' ')
echo "===== $(date): Creating controls completed (high-stringent dataset) for $CTL_COUNT files."