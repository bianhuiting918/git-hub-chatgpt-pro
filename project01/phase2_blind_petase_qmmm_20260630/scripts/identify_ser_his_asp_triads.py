from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path


def parse_pdb_atoms(path: Path):
    residues = defaultdict(dict)
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            altloc = line[16].strip()
            if altloc not in ("", "A"):
                continue
            resname = line[17:20].strip()
            chain = line[21].strip() or "_"
            try:
                resseq = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            icode = line[26].strip()
            atom = line[12:16].strip()
            residues[(chain, resseq, icode, resname)][atom] = (x, y, z)
    return residues


def dist(a, b) -> float:
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def min_named_distance(residue_a, names_a, residue_b, names_b):
    best = None
    best_names = None
    for name_a in names_a:
        if name_a not in residue_a:
            continue
        for name_b in names_b:
            if name_b not in residue_b:
                continue
            value = dist(residue_a[name_a], residue_b[name_b])
            if best is None or value < best:
                best = value
                best_names = (name_a, name_b)
    return best, best_names


def residue_label(key) -> str:
    chain, resseq, icode, resname = key
    suffix = icode if icode else ""
    return f"{chain}:{resname}{resseq}{suffix}"


def find_triads(path: Path):
    residues = parse_pdb_atoms(path)
    ser = [(key, atoms) for key, atoms in residues.items() if key[3] == "SER" and "OG" in atoms]
    his = [
        (key, atoms)
        for key, atoms in residues.items()
        if key[3] in {"HIS", "HID", "HIE", "HIP"} and ("ND1" in atoms or "NE2" in atoms)
    ]
    asp = [(key, atoms) for key, atoms in residues.items() if key[3] == "ASP" and ("OD1" in atoms or "OD2" in atoms)]

    rows = []
    for ser_key, ser_atoms in ser:
        for his_key, his_atoms in his:
            sh_dist, sh_atoms = min_named_distance(ser_atoms, ["OG"], his_atoms, ["ND1", "NE2"])
            if sh_dist is None or sh_dist > 4.0:
                continue
            for asp_key, asp_atoms in asp:
                hd_dist, hd_atoms = min_named_distance(his_atoms, ["ND1", "NE2"], asp_atoms, ["OD1", "OD2"])
                if hd_dist is None or hd_dist > 4.0:
                    continue
                rows.append(
                    {
                        "pdb": path.stem.upper(),
                        "ser": residue_label(ser_key),
                        "his": residue_label(his_key),
                        "asp": residue_label(asp_key),
                        "ser_his_min_A": f"{sh_dist:.2f}",
                        "ser_his_atoms": "-".join(sh_atoms),
                        "his_asp_min_A": f"{hd_dist:.2f}",
                        "his_asp_atoms": "-".join(hd_atoms),
                    }
                )
    rows.sort(key=lambda row: (row["pdb"], float(row["ser_his_min_A"]) + float(row["his_asp_min_A"])))
    return rows


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: identify_ser_his_asp_triads.py <pdb-dir> [out-tsv]", file=sys.stderr)
        return 2
    pdb_dir = Path(sys.argv[1])
    out_path = Path(sys.argv[2]) if len(sys.argv) > 2 else None

    rows = []
    for path in sorted(pdb_dir.glob("*.pdb")):
        rows.extend(find_triads(path))

    fieldnames = ["pdb", "ser", "his", "asp", "ser_his_min_A", "ser_his_atoms", "his_asp_min_A", "his_asp_atoms"]
    if out_path:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    else:
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
