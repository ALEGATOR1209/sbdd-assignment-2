#!/usr/bin/env python3

from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

LIGANDS = Path("data/ligands")
SMILES = LIGANDS / "smiles"
MOL = LIGANDS / "mol"
PDBQT = LIGANDS / "pdbqt"


def prepare_mols():
  for ligand in SMILES.iterdir():
    mol = Chem.MolFromSmiles(ligand.read_text())
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)

    output = MOL / (ligand.name.removesuffix(".smiles") + ".mol")
    output.parent.mkdir(exist_ok=True)

    Chem.MolToMolFile(mol, str(output))


def prepare_pbqt():
  for ligand in MOL.iterdir():
    preparator = MoleculePreparation()
    supplier = Chem.SDMolSupplier(str(ligand), removeHs=False)

    for i, mol in enumerate(supplier, start=1):
      if mol is None:
        print(f"Record {i}: RDKit couldn't parse it, skipping")
        continue

    for setup in preparator.prepare(mol):
      pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
      if is_ok:
        output = PDBQT / (ligand.name.removesuffix(".mol") + ".pdbqt")
        output.parent.mkdir(exist_ok=True)

        with open(output, "w") as f:
          f.write(pdbqt_string)
      else:
        print(f"Record {i}: {error_msg}")


def prepare_ligands():
  prepare_mols()
  prepare_pbqt()

if __name__ == "__main__":
  prepare_ligands()
