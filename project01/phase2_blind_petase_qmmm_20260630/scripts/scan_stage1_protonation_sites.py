from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

TRIAD_BY_PDB = {
    "5XJH": ("A", 160, "SER", 237, "HIS", 206, "ASP"),
    "5YFE": ("A", 134, "SER", 211, "HIS", 180, "ASP"),
    "6EQD": ("A", 160, "SER", 237, "HIS", 206, "ASP"),
    "6EQE": ("A", 160, "SER", 237, "HIS", 206, "ASP"),
    "6EQF": ("A", 160, "SER", 237, "HIS", 206, "ASP"),
    "6EQG": ("A", 160, "SER", 237, "HIS", 206, "ASP"),
    "6EQH": ("A", 160, "SER", 237, "HIS", 206, "ASP"),
    "6ILW": ("A", 160, "SER", 237, "HIS", 206, "ASP"),
    "6QGC": ("A", 160, "SER", 237, "HIS", 206, "ASP"),
}
TITRATABLE = {"ASP", "GLU", "HIS", "LYS", "ARG", "CYS", "TYR"}
SIDECHAIN_ATOMS = {
    "ASP": {"OD1", "OD2"},
    "GLU": {"OE1", "OE2"},
    "HIS": {"ND1", "NE2"},
    "LYS": {"NZ"},
    "ARG": {"NE", "NH1", "NH2"},
    "CYS": {"SG"},
    "TYR": {"OH"},
    "SER": {"OG"},
}
DEFAULT_STATE = {
    "ASP": "deprotonated_minus1",
    "GLU": "deprotonated_minus1",
    "HIS": "neutral_tautomer_to_assign",
    "LYS": "protonated_plus1",
    "ARG": "protonated_plus1",
    "CYS": "neutral_thiol",
    "TYR": "neutral_phenol",
}


def dist(a, b) -> float:
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def parse_pdb(path: Path):
    residues = defaultdict(dict)
    waters = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom = line[12:16].strip()
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
        record = line[:6].strip()
        key = (chain, resseq, icode, resname)
        if record == "HETATM" and resname in {"HOH", "WAT", "DOD"} and atom in {"O", "OW"}:
            waters.append((f"{chain}:{resname}{resseq}{icode}", (x, y, z)))
        elif record == "ATOM":
            residues[key][atom] = (x, y, z)
    return residues, waters


def label(key) -> str:
    chain, resseq, icode, resname = key
    return f"{chain}:{resname}{resseq}{icode}"


def min_distance_to_points(atoms: dict[str, tuple[float, float, float]], names: set[str], points):
    best = None
    best_atom = ""
    best_target = ""
    for name in names:
        if name not in atoms:
            continue
        for target_label, target_point in points:
            value = dist(atoms[name], target_point)
            if best is None or value < best:
                best = value
                best_atom = name
                best_target = target_label
    return best, best_atom, best_target


def analyze_histidine(atoms: dict[str, tuple[float, float, float]], environment_points) -> tuple[str, str]:
    nd1_best, _, nd1_target = min_distance_to_points(atoms, {"ND1"}, environment_points)
    ne2_best, _, ne2_target = min_distance_to_points(atoms, {"NE2"}, environment_points)
    notes = []
    if nd1_best is not None:
        notes.append(f"ND1_nearest={nd1_best:.2f}A_to_{nd1_target}")
    if ne2_best is not None:
        notes.append(f"NE2_nearest={ne2_best:.2f}A_to_{ne2_target}")
    if nd1_best is None or ne2_best is None:
        return "neutral_tautomer_to_assign", ";".join(notes)
    if ne2_best + 0.35 < nd1_best:
        return "HID_candidate_delta_protonated_NE2_acceptor", ";".join(notes)
    if nd1_best + 0.35 < ne2_best:
        return "HIE_candidate_epsilon_protonated_ND1_acceptor", ";".join(notes)
    return "HIS_tautomer_ambiguous_sensitivity_required", ";".join(notes)


def scan_one(path: Path):
    pdb = path.stem.split("_")[0].upper()
    residues, waters = parse_pdb(path)
    triad = TRIAD_BY_PDB[pdb]
    chain, ser_i, ser_name, his_i, his_name, asp_i, asp_name = triad
    triad_keys = [(chain, ser_i, "", ser_name), (chain, his_i, "", his_name), (chain, asp_i, "", asp_name)]
    triad_points = []
    for key in triad_keys:
        atoms = residues.get(key, {})
        for atom_name in SIDECHAIN_ATOMS.get(key[3], set()):
            if atom_name in atoms:
                triad_points.append((f"{label(key)}:{atom_name}", atoms[atom_name]))
    water_points = [(f"{water_label}:O", point) for water_label, point in waters]
    environment_points = triad_points + water_points
    rows = []
    for key, atoms in sorted(residues.items(), key=lambda item: (item[0][0], item[0][1], item[0][2], item[0][3])):
        resname = key[3]
        if resname not in TITRATABLE:
            continue
        distance_to_triad, atom_name, target = min_distance_to_points(atoms, SIDECHAIN_ATOMS[resname], triad_points)
        _, _, env_target = min_distance_to_points(atoms, SIDECHAIN_ATOMS[resname], environment_points)
        include = False
        reason = []
        if key in triad_keys:
            include = True
            reason.append("catalytic_triad_member")
        if distance_to_triad is not None and distance_to_triad <= 8.0:
            include = True
            reason.append("within_8A_of_triad")
        if resname == "HIS":
            include = True
            reason.append("all_histidines_require_tautomer_assignment")
        if not include:
            continue
        proposed = DEFAULT_STATE[resname]
        note = ""
        sensitivity = "no"
        if key == (chain, his_i, "", his_name):
            proposed = "catalytic_HIS_neutral_tautomer_sensitivity_required"
            sensitivity = "yes"
            note = "test_HID_vs_HIE_vs_HIP_only_if_pH_or_network_supports"
        elif key == (chain, asp_i, "", asp_name):
            proposed = "catalytic_ASP_deprotonated_primary"
            sensitivity = "yes"
            note = "also_test_neutral_ASP_if_His_Asp_Hbond_network_forces_it"
        elif resname == "HIS":
            proposed, note = analyze_histidine(atoms, environment_points)
            sensitivity = "yes" if "ambiguous" in proposed or "candidate" in proposed else "no"
        elif distance_to_triad is not None and distance_to_triad <= 4.0 and resname in {"ASP", "GLU", "CYS", "TYR"}:
            sensitivity = "yes"
            note = "near_catalytic_center_check_pKa_shift"
        rows.append({
            "pdb": pdb,
            "prepared_pdb": str(path),
            "residue": label(key),
            "resname": resname,
            "default_pH7_state": DEFAULT_STATE[resname],
            "proposed_stage1_state": proposed,
            "sensitivity_required": sensitivity,
            "nearest_triad_distance_A": "" if distance_to_triad is None else f"{distance_to_triad:.2f}",
            "nearest_triad_atom": atom_name,
            "nearest_triad_target": target,
            "nearest_environment_target": env_target,
            "selection_reason": ";".join(reason),
            "note": note,
        })
    return rows


def write_tsv(path: Path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: scan_stage1_protonation_sites.py <prepared-pdb-dir> <out-dir>", file=sys.stderr)
        return 2
    pdb_dir = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    rows = []
    for path in sorted(pdb_dir.glob("*_chainA_initial_clean_v2.pdb")):
        rows.extend(scan_one(path))
    fieldnames = [
        "pdb", "prepared_pdb", "residue", "resname", "default_pH7_state", "proposed_stage1_state",
        "sensitivity_required", "nearest_triad_distance_A", "nearest_triad_atom", "nearest_triad_target",
        "nearest_environment_target", "selection_reason", "note",
    ]
    write_tsv(out_dir / "protonation_site_scan.tsv", rows, fieldnames)
    primary = [row for row in rows if row["pdb"] == "6EQE"]
    hypothesis_rows = [row for row in primary if row["sensitivity_required"] == "yes" or "catalytic" in row["selection_reason"]]
    write_tsv(out_dir / "protonation_hypothesis_manifest.tsv", hypothesis_rows, fieldnames)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
