#!/bin/bash
set -euo pipefail
shopt -s nullglob

VINA="./bin/vina"
RECEPTOR="data/target/7G0B.pdbqt"
CONFIG="data/target/7G0B.box.txt"
LIGANDS="data/ligands/pdbqt"
OUTPUT="results/vina"

mkdir -p "$OUTPUT"

ligand_files=("$LIGANDS"/*.pdbqt)
if [ ${#ligand_files[@]} -eq 0 ]; then
  echo "No .pdbqt files found in $LIGANDS" >&2
  exit 1
fi

for ligand in "${ligand_files[@]}"; do
  name=$(basename "$ligand" .pdbqt)
  echo "Docking $name..."

  "$VINA" \
    --receptor "$RECEPTOR" \
    --ligand "$ligand" \
    --config "$CONFIG" \
    --exhaustiveness 16 \
    --num_modes 9 \
    --out "$OUTPUT/${name}_poses.pdbqt" \
    | tee "$OUTPUT/${name}.log"
done
