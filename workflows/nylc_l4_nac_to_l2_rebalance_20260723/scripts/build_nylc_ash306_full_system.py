#!/usr/bin/env python3
"""Build a neutral full-system NylC-C18 ASH306 topology preflight."""
import argparse
import hashlib
import json
import math
import pathlib
import re
import shutil

import numpy as np

SOURCE_GRO_SHA256 = "182f68a41c72c490e0046323c9b8a0f8514ccf7af7dc4d578dcb87342c1c8e2d"
PROBE_GRO_SHA256 = "a3998b6b5f5df5a40c78861dfd4b65983dff5d96bba1cd7ed238c7d57528ea44"
PROBE_TOP_SHA256 = "630b45eb20713b1f604ae11a193ded778edcfb52f44b70c9377d6858bd7690f9"
CHAIN_H_FIRST = 8949
CHAIN_H_LAST = 10272
L2_FIRST = 10273
L2_LAST = 10351
SOURCE_ATOMS = 133589
FINAL_ATOMS = 133589
GATE_RESIDUES = tuple(range(261, 267))  # Thr267 is excluded from the gate group.
HEAVY_SHIFT_TOLERANCE_NM = 0.0011
NAC_DISTANCE_MAX_NM = 0.35
NAC_ANGLE_MIN_DEG = 95.0
NAC_ANGLE_MAX_DEG = 115.0
MIN_HEAVY_CONTACT_NM = 0.20


def sha256(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_gro(path):
    lines = pathlib.Path(path).read_text().splitlines()
    count = int(lines[1])
    atom_lines = lines[2:2 + count]
    if len(atom_lines) != count or len(lines) < count + 3:
        raise ValueError(f"malformed GRO: {path}")
    atoms = []
    for source_index, line in enumerate(atom_lines, 1):
        atoms.append({
            "source_index": source_index,
            "resid": int(line[:5]),
            "resname": line[5:10].strip(),
            "atomname": line[10:15].strip(),
            "xyz_nm": tuple(float(line[a:b]) for a, b in ((20, 28), (28, 36), (36, 44))),
            "line": line,
        })
    fields = [float(value) for value in lines[2 + count].split()]
    if len(fields) == 3:
        cell = ((fields[0], 0.0, 0.0), (0.0, fields[1], 0.0), (0.0, 0.0, fields[2]))
    elif len(fields) == 9:
        cell = (
            (fields[0], fields[3], fields[4]),
            (fields[5], fields[1], fields[6]),
            (fields[7], fields[8], fields[2]),
        )
    else:
        raise ValueError(f"unsupported GRO box with {len(fields)} fields")
    return lines[0], atoms, lines[2 + count], cell


def is_heavy(atom):
    return not atom["atomname"].upper().startswith("H")


def minimum_periodic_distances(point, coordinates, cell):
    point_array = np.asarray(point, dtype=float)
    coordinates_array = np.asarray(coordinates, dtype=float)
    cell_array = np.asarray(cell, dtype=float)
    inverse = np.linalg.inv(cell_array)
    fractional = (coordinates_array - point_array) @ inverse
    fractional -= np.rint(fractional)
    shifts = np.asarray(
        [(i, j, k) for i in (-1.0, 0.0, 1.0) for j in (-1.0, 0.0, 1.0) for k in (-1.0, 0.0, 1.0)]
    )
    cartesian = (fractional[:, None, :] + shifts[None, :, :]) @ cell_array
    return np.linalg.norm(cartesian, axis=2).min(axis=1)


def select_farthest_sodium(sodiums, solute_coordinates, cell):
    if not sodiums or not solute_coordinates:
        raise ValueError("sodium and solute coordinates are required")
    ranked = []
    for sodium in sodiums:
        minimum = float(minimum_periodic_distances(sodium["xyz_nm"], solute_coordinates, cell).min())
        ranked.append({**sodium, "minimum_solute_distance_nm": minimum})
    return max(ranked, key=lambda item: (item["minimum_solute_distance_nm"], -item["source_index"]))


def update_na_molecule_count(text):
    pattern = re.compile(r"(?m)^(\s*NA\s+)144(\s*(?:;.*)?)$")
    updated, count = pattern.subn(r"\g<1>143\g<2>", text)
    if count != 1:
        raise ValueError(f"expected one NA 144 molecule line, observed {count}")
    return updated


def extract_chain_itp(text):
    lines = text.splitlines()
    start = next((index for index, line in enumerate(lines) if re.match(r"^\s*\[\s*moleculetype\s*\]", line)), None)
    if start is None:
        raise ValueError("chain topology has no moleculetype section")
    stop = len(lines)
    for index in range(start + 1, len(lines)):
        stripped = lines[index].strip()
        if stripped.startswith("#ifdef POSRES") or stripped == "; Include Position restraint file":
            stop = index
            break
        if re.match(r"^\[\s*system\s*\]", stripped):
            stop = index
            break
    result = "\n".join(lines[start:stop]).rstrip() + "\n"
    if "Protein_chain_H" not in result or "[ atoms ]" not in result:
        raise ValueError("extracted chain topology is incomplete")
    return result


def max_chain_heavy_shift(source_atoms, probe_atoms):
    source = {(atom["resid"], atom["atomname"]): atom["xyz_nm"] for atom in source_atoms if is_heavy(atom)}
    probe = {(atom["resid"], atom["atomname"]): atom["xyz_nm"] for atom in probe_atoms if is_heavy(atom)}
    if set(source) != set(probe):
        raise ValueError("chain heavy-atom identity mismatch")
    terminal = {(355, "OC1"), (355, "OC2")}
    ordinary = set(source) - terminal
    values = [
        max(abs(a - b) for a, b in zip(source[key], probe[key]))
        for key in ordinary
    ]
    direct = max(
        max(abs(a - b) for a, b in zip(source[key], probe[key]))
        for key in terminal
    )
    swapped = max(
        max(abs(a - b) for a, b in zip(source[(355, left)], probe[(355, right)]))
        for left, right in (("OC1", "OC2"), ("OC2", "OC1"))
    )
    values.append(min(direct, swapped))
    return max(values, default=0.0), "minimum_of_direct_and_OC1_OC2_swapped"


def rewrite_atom_numbers(atoms):
    lines = []
    for new_index, atom in enumerate(atoms, 1):
        number = new_index % 100000
        lines.append(atom["line"][:15] + f"{number:5d}" + atom["line"][20:])
        atom["new_index"] = new_index
    return lines


def minimum_image_vector(vector, cell):
    cell_array = np.asarray(cell, dtype=float)
    fractional = np.asarray(vector, dtype=float) @ np.linalg.inv(cell_array)
    fractional -= np.rint(fractional)
    shifts = np.asarray(
        [(i, j, k) for i in (-1.0, 0.0, 1.0) for j in (-1.0, 0.0, 1.0) for k in (-1.0, 0.0, 1.0)]
    )
    candidates = (fractional[None, :] + shifts) @ cell_array
    return candidates[np.argmin(np.linalg.norm(candidates, axis=1))]


def angle_deg(left, right):
    left_norm = np.linalg.norm(left)
    right_norm = np.linalg.norm(right)
    cosine = float(np.dot(left, right) / (left_norm * right_norm))
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def unique(atoms, resid, resname, atomname):
    matches = [
        atom for atom in atoms
        if atom["resid"] == resid and atom["resname"] == resname and atom["atomname"] == atomname
    ]
    if len(matches) != 1:
        raise ValueError(f"expected one {resid}:{resname}:{atomname}, observed {len(matches)}")
    return matches[0]


def copy_topology_inputs(build_dir, output_dir):
    for path in build_dir.glob("*.itp"):
        if path.name != "topol_Protein_chain_H.itp":
            shutil.copy2(path, output_dir / path.name)
    mdp = build_dir / "em_preflight.processed.mdp"
    if not mdp.is_file():
        raise ValueError("stable em_preflight.processed.mdp is missing")
    shutil.copy2(mdp, output_dir / "grompp_preflight.mdp")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-gro", type=pathlib.Path, required=True)
    parser.add_argument("--source-build-dir", type=pathlib.Path, required=True)
    parser.add_argument("--probe-dir", type=pathlib.Path, required=True)
    parser.add_argument("--output-dir", type=pathlib.Path, required=True)
    args = parser.parse_args()

    source_gro = args.source_gro.resolve()
    build_dir = args.source_build_dir.resolve()
    probe_dir = args.probe_dir.resolve()
    output = args.output_dir.resolve()
    if output.exists():
        raise ValueError(f"refusing to overwrite existing output: {output}")
    if sha256(source_gro) != SOURCE_GRO_SHA256:
        raise ValueError("immutable source GRO SHA mismatch")
    probe_gro = probe_dir / "chainH_ASH306.gro"
    probe_top = probe_dir / "chainH.top"
    if sha256(probe_gro) != PROBE_GRO_SHA256 or sha256(probe_top) != PROBE_TOP_SHA256:
        raise ValueError("authoritative ASH306 probe SHA mismatch")
    probe_pass = json.loads((probe_dir / "PASS.json").read_text())
    if probe_pass.get("status") != "PROBE_PASS_ASH306_CHAIN":
        raise ValueError("authoritative ASH306 probe is not PASS")

    title, source_atoms, box_line, cell = read_gro(source_gro)
    _, probe_atoms, _, _ = read_gro(probe_gro)
    if len(source_atoms) != SOURCE_ATOMS:
        raise ValueError(f"source atom count changed: {len(source_atoms)}")
    if len(probe_atoms) != CHAIN_H_LAST - CHAIN_H_FIRST + 2:
        raise ValueError(f"probe chain atom count changed: {len(probe_atoms)}")
    chain_source = source_atoms[CHAIN_H_FIRST - 1:CHAIN_H_LAST]
    maximum_shift, symmetry_rule = max_chain_heavy_shift(chain_source, probe_atoms)
    if maximum_shift > HEAVY_SHIFT_TOLERANCE_NM:
        raise ValueError(f"probe heavy coordinates moved by {maximum_shift} nm")

    ash_atoms = [atom for atom in probe_atoms if atom["resid"] == 306 and atom["resname"] == "ASH"]
    thr_atoms = [atom for atom in probe_atoms if atom["resid"] == 267 and atom["resname"] == "THR"]
    if len(ash_atoms) != 13 or "HD2" not in {atom["atomname"] for atom in ash_atoms}:
        raise ValueError("ASH306 is incomplete")
    if not {"H1", "H2", "H3", "OG1", "HG1"}.issubset({atom["atomname"] for atom in thr_atoms}):
        raise ValueError("N-terminal Thr267 is incomplete")

    solute = [
        atom["xyz_nm"] for atom in source_atoms[:L2_LAST]
        if is_heavy(atom)
    ]
    sodiums = [atom for atom in source_atoms if atom["resname"] == "NA" and atom["atomname"] == "NA"]
    if len(sodiums) != 144:
        raise ValueError(f"expected 144 Na atoms, observed {len(sodiums)}")
    removed = select_farthest_sodium(sodiums, solute, cell)

    inserted = []
    for atom in probe_atoms:
        inserted.append({**atom, "source_index": None, "origin": "ash306_probe"})
    merged = [
        {**atom, "origin": "source"} for atom in source_atoms[:CHAIN_H_FIRST - 1]
    ] + inserted + [
        {**atom, "origin": "source"} for atom in source_atoms[CHAIN_H_LAST:]
        if atom["source_index"] != removed["source_index"]
    ]
    if len(merged) != FINAL_ATOMS:
        raise ValueError(f"final atom count mismatch: {len(merged)}")

    output.mkdir(parents=True, exist_ok=False)
    copy_topology_inputs(build_dir, output)
    chain_itp = extract_chain_itp(probe_top.read_text())
    (output / "topol_Protein_chain_H.itp").write_text(chain_itp)
    topol = update_na_molecule_count((build_dir / "topol.top").read_text())
    (output / "topol.top").write_text(topol)
    gro_lines = rewrite_atom_numbers(merged)
    (output / "system.gro").write_text(
        f"{title} | neutral ASH306 full-system preflight\n{FINAL_ATOMS}\n"
        + "\n".join(gro_lines) + "\n" + box_line + "\n"
    )

    active_chain = merged[CHAIN_H_FIRST - 1:CHAIN_H_FIRST - 1 + len(probe_atoms)]
    ligand = merged[CHAIN_H_FIRST - 1 + len(probe_atoms):CHAIN_H_FIRST - 1 + len(probe_atoms) + 79]
    thr = unique(active_chain, 267, "THR", "OG1")
    carbon = unique(ligand, 356, "L2", "C12")
    oxygen = unique(ligand, 356, "L2", "O2")
    attack = minimum_image_vector(np.asarray(thr["xyz_nm"]) - np.asarray(carbon["xyz_nm"]), cell)
    carbonyl = minimum_image_vector(np.asarray(oxygen["xyz_nm"]) - np.asarray(carbon["xyz_nm"]), cell)
    nac_distance = float(np.linalg.norm(attack))
    nac_angle = angle_deg(attack, carbonyl)
    protein_heavy = [atom["xyz_nm"] for atom in merged[:CHAIN_H_FIRST - 1 + len(probe_atoms)] if is_heavy(atom)]
    ligand_heavy = [atom["xyz_nm"] for atom in ligand if is_heavy(atom)]
    minimum_contact = min(
        float(minimum_periodic_distances(point, protein_heavy, cell).min())
        for point in ligand_heavy
    )
    if not (nac_distance <= NAC_DISTANCE_MAX_NM and NAC_ANGLE_MIN_DEG <= nac_angle <= NAC_ANGLE_MAX_DEG):
        raise ValueError(f"NAC geometry changed: {nac_distance} nm, {nac_angle} deg")
    if minimum_contact < MIN_HEAVY_CONTACT_NM:
        raise ValueError(f"ligand-protein heavy clash: {minimum_contact} nm")

    audit = {
        "schema_version": 1,
        "candidate_id": "nylc_C18_trueT267_freeGS",
        "status": "READY_FOR_GROMPP_ASH306_FULL_SYSTEM",
        "source_gro_sha256": SOURCE_GRO_SHA256,
        "probe_job_id": 61717760,
        "probe_gro_sha256": PROBE_GRO_SHA256,
        "probe_top_sha256": PROBE_TOP_SHA256,
        "source_atom_count": SOURCE_ATOMS,
        "final_atom_count": FINAL_ATOMS,
        "chain_h_atom_count": len(probe_atoms),
        "ash306_atom_count": len(ash_atoms),
        "ash306_has_hd2": "HD2" in {atom["atomname"] for atom in ash_atoms},
        "thr267_terminal_atoms_present": sorted({"H1", "H2", "H3", "OG1", "HG1"} & {atom["atomname"] for atom in thr_atoms}),
        "gate_residues": list(GATE_RESIDUES),
        "thr267_excluded_from_gate_group": True,
        "maximum_chain_heavy_coordinate_shift_nm": maximum_shift,
        "terminal_oxygen_symmetry_rule": symmetry_rule,
        "removed_sodium": {
            "source_global_atom": removed["source_index"],
            "residue_number": removed["resid"],
            "coordinate_nm": list(removed["xyz_nm"]),
            "minimum_solute_distance_nm": removed["minimum_solute_distance_nm"],
        },
        "molecule_counts": {
            "Protein_chain_A": 1,
            "Protein_chain_E": 1,
            "Protein_chain_D": 1,
            "Protein_chain_H": 1,
            "PA66_L2": 1,
            "SOL": 40990,
            "NA": 143,
            "CL": 124,
        },
        "active_atom_mapping": {
            "thr267_og1_new": thr["new_index"],
            "l2_c12_new": carbon["new_index"],
            "l2_o2_new": oxygen["new_index"],
            "l2_first_new": ligand[0]["new_index"],
            "l2_last_new": ligand[-1]["new_index"],
        },
        "nac_distance_nm": nac_distance,
        "nac_angle_deg": nac_angle,
        "minimum_ligand_protein_heavy_distance_nm": minimum_contact,
        "position_restraints_transferred": False,
        "next_gate": "GROMACS grompp plus ParmEd atom-count and total-charge audit",
        "interpretation": "Neutral ASH306 full-system build only; not a TS, RC, PMF, or barrier.",
    }
    (output / "build_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    print(json.dumps(audit, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
