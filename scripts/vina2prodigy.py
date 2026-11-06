#!/usr/bin/env python3
"""
Build a PRODIGY-LIG-ready complex from a receptor PDB and a Vina PDBQT.

    python3 vina2prodigy.py receptor.pdb poses.pdbqt

Writes complex_pose1.pdb ... and prints the prodigy_lig command to run.
Options:
    --pose N        only this pose (1-based). Default: all models in the file.
    --chain X       chain for the ligand (default A, i.e. the receptor's chain)
    --resname XXX   ligand residue name (default LIG)
    --resnum N      ligand residue number (default 900)
    --run           run prodigy_lig directly and print the table

PRODIGY-LIG ignores hydrogens and does not care whether the ligand is
ATOM or HETATM, so no protonation or element-column fixing is needed.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys


def read_models(pdbqt_path):
    """Yield (pose_index, vina_score, [atom_lines]) for each MODEL."""
    models = []
    cur, score, idx = None, None, 0
    with open(pdbqt_path) as fh:
        for line in fh:
            if line.startswith("MODEL"):
                idx += 1
                cur, score = [], None
            elif line.startswith("ENDMDL"):
                if cur is not None:
                    models.append((idx, score, cur))
                cur, score = None, None
            elif line.startswith("REMARK") and "VINA RESULT" in line:
                m = re.search(r"VINA RESULT:\s*(-?\d+\.\d+)", line)
                if m:
                    score = float(m.group(1))
            elif line.startswith(("ATOM", "HETATM")):
                if cur is None:          # single-model file, no MODEL record
                    idx = 1
                    cur, score = [], score
                cur.append(line.rstrip("\n"))
    if cur:                              # file ended without ENDMDL
        models.append((idx if idx else 1, score, cur))
    return models


def to_pdb_lines(atom_lines, chain, resname, resnum):
    """PDBQT atom line -> PDB HETATM line.

    A PDBQT record is a PDB record with partial charge and AutoDock type
    appended past column 70, so truncating at 66 and rewriting the
    residue fields is sufficient.
    """
    out = []
    for line in atom_lines:
        line = line.ljust(80)
        name = line[12:16]
        x, y, z = line[30:38], line[38:46], line[46:54]
        # element: trust the atom name's leading letters, not the AD type
        raw = name.strip()
        el = raw[0] if raw[0].isalpha() else (raw[1] if len(raw) > 1 else "C")
        if raw[:2].capitalize() in ("Cl", "Br", "Fe", "Zn", "Mg", "Mn"):
            el = raw[:2].capitalize()
        out.append("HETATM%5d %-4s %3s %1s%4d    %8s%8s%8s  1.00 20.00          %2s"
                   % (len(out) + 1, name.strip()[:4], resname[:3], chain,
                      resnum, x, y, z, el.rjust(2)))
    return out


def receptor_lines(pdb_path):
    """Protein ATOM/HETATM records, dropping waters, END and model markers."""
    keep = []
    for line in open(pdb_path):
        rec = line[:6]
        if rec.startswith(("ATOM", "HETATM")):
            if line[17:20].strip() in ("HOH", "WAT", "DOD"):
                continue
            keep.append(line.rstrip("\n"))
        elif rec.startswith("TER"):
            keep.append(line.rstrip("\n"))
    return keep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("receptor")
    ap.add_argument("pdbqt")
    ap.add_argument("--pose", type=int, default=None)
    ap.add_argument("--chain", default="A")
    ap.add_argument("--resname", default="LIG")
    ap.add_argument("--resnum", type=int, default=900)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--run", action="store_true")
    args = ap.parse_args()

    rec = receptor_lines(args.receptor)
    if not rec:
        sys.exit("No atoms read from %s" % args.receptor)

    models = read_models(args.pdbqt)
    if not models:
        sys.exit("No poses read from %s" % args.pdbqt)
    if args.pose:
        models = [m for m in models if m[0] == args.pose]
        if not models:
            sys.exit("Pose %d not in file (%d poses present)."
                     % (args.pose, len(read_models(args.pdbqt))))

    os.makedirs(args.outdir, exist_ok=True)
    written = []
    for idx, score, atoms in models:
        lig = to_pdb_lines(atoms, args.chain, args.resname, args.resnum)
        path = os.path.join(args.outdir, "complex_pose%d.pdb" % idx)
        with open(path, "w") as fh:
            fh.write("\n".join(rec) + "\n")
            fh.write("\n".join(lig) + "\n")
            fh.write("END\n")
        written.append((path, idx, score, len(lig)))
        print("pose %-2d  vina %s kcal/mol  %2d ligand atoms  -> %s"
              % (idx, ("%.2f" % score) if score is not None else "  n/a",
                 len(lig), path))

    cmd_tmpl = "prodigy_lig -c %s %s:%s -i %%s -v" % (
        args.chain, args.chain, args.resname)
    print("\nRun:\n  " + cmd_tmpl % written[0][0])

    if args.run:
        exe = shutil.which("prodigy_lig")
        if not exe:
            sys.exit("\nprodigy_lig not on PATH -- pip install prodigy-lig")
        print("\npose\tvina\tprodigy_dG")
        for path, idx, score, _ in written:
            res = subprocess.run(
                [exe, "-c", args.chain, "%s:%s" % (args.chain, args.resname),
                 "-i", path, "-v"],
                capture_output=True, text=True)
            if res.returncode != 0:
                print("%d\t%s\tFAILED: %s"
                      % (idx, score, res.stderr.strip().splitlines()[-1:]))
                continue
            rows = [r for r in res.stdout.strip().splitlines()
                    if not r.startswith("Job name")]
            dg = rows[0].split("\t")[1].strip() if rows else "?"
            print("%d\t%s\t%s"
                  % (idx, ("%.2f" % score) if score is not None else "n/a", dg))


if __name__ == "__main__":
    main()
