#!/bin/bash

X=13.79
Y=-6.75
Z=-9.33

W=20
H=20
D=20

INPUT="data/target/7G0B_cleaned.pdb"
OUTPUT="data/target/7G0B"

mk_prepare_receptor.py -i "$INPUT" \
  -o "$OUTPUT" \
  -p -v \
  --box_center "$X" "$Y" "$Z" \
  --box_size "$W" "$H" "$D" \
  --default_altloc A
