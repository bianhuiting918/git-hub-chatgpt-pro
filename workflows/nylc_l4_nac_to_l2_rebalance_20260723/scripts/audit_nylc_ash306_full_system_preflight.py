#!/usr/bin/env python3
"""Audit GROMACS/ParmEd entry for the neutral NylC ASH306 full system."""
import argparse
import hashlib
import json
import os
import pathlib
import re

import parmed as pmd

EXPECTED_ATOMS = 133589
CHARGE_TOLERANCE = 1.0e-4
PROBE_STATUS = "PROBE_PASS_ASH306_CHAIN"


def sha256(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=pathlib.Path, required=True)
    parser.add_argument("--probe-pass", type=pathlib.Path, required=True)
    parser.add_argument("--grompp-log", type=pathlib.Path, required=True)
    parser.add_argument("--amber-output-dir", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()

    build = args.build_dir.resolve()
    amber_output = args.amber_output_dir.resolve()
    build_audit = json.loads((build / "build_audit.json").read_text())
    probe = json.loads(args.probe_pass.read_text())
    grompp_text = args.grompp_log.read_text(errors="replace")
    if probe.get("status") != PROBE_STATUS:
        raise ValueError("probe is not PROBE_PASS_ASH306_CHAIN")
    if build_audit.get("status") != "READY_FOR_GROMPP_ASH306_FULL_SYSTEM":
        raise ValueError("builder did not reach grompp gate")
    if re.search(r"fatal error|atom count mismatch|number of coordinates.*does not match", grompp_text, re.I):
        raise ValueError("grompp log contains a fatal topology/coordinate error")

    previous = pathlib.Path.cwd()
    os.chdir(build)
    try:
        structure = pmd.load_file("topol.top", xyz="system.gro")
    finally:
        os.chdir(previous)
    atom_count = len(structure.atoms)
    total_charge = float(sum(atom.charge for atom in structure.atoms))
    ash_residues = [residue for residue in structure.residues if residue.name == "ASH"]
    ash306_atom_names = sorted(atom.name for residue in ash_residues for atom in residue.atoms)
    ash306_charge = float(sum(atom.charge for residue in ash_residues for atom in residue.atoms))
    l2_residues = [residue for residue in structure.residues if residue.name == "L2"]
    l2_atoms = sum(len(residue.atoms) for residue in l2_residues)

    passed = (
        atom_count == EXPECTED_ATOMS
        and abs(total_charge) <= CHARGE_TOLERANCE
        and len(ash_residues) == 1
        and len(ash306_atom_names) == 13
        and "HD2" in ash306_atom_names
        and abs(ash306_charge) <= CHARGE_TOLERANCE
        and len(l2_residues) == 1
        and l2_atoms == 79
        and build_audit["nac_distance_nm"] <= 0.35
        and 95.0 <= build_audit["nac_angle_deg"] <= 115.0
        and build_audit["minimum_ligand_protein_heavy_distance_nm"] >= 0.20
    )

    result = {
        "schema_version": 1,
        "candidate_id": "nylc_C18_trueT267_freeGS",
        "status": "PASS_ASH306_FULL_SYSTEM_PREFLIGHT" if passed else "FAIL_ASH306_FULL_SYSTEM_PREFLIGHT",
        "probe_status": probe.get("status"),
        "atom_count": atom_count,
        "expected_atom_count": EXPECTED_ATOMS,
        "total_topology_charge_e": total_charge,
        "ash306_residue_count": len(ash_residues),
        "ash306_atom_names": ash306_atom_names,
        "ash306_topology_charge_e": ash306_charge,
        "l2_residue_count": len(l2_residues),
        "l2_atom_count": l2_atoms,
        "nac_distance_nm": build_audit["nac_distance_nm"],
        "nac_angle_deg": build_audit["nac_angle_deg"],
        "minimum_ligand_protein_heavy_distance_nm": build_audit["minimum_ligand_protein_heavy_distance_nm"],
        "removed_sodium": build_audit["removed_sodium"],
        "active_atom_mapping": build_audit["active_atom_mapping"],
        "grompp_log_sha256": sha256(args.grompp_log),
        "system_gro_sha256": sha256(build / "system.gro"),
        "topol_top_sha256": sha256(build / "topol.top"),
        "chain_h_itp_sha256": sha256(build / "topol_Protein_chain_H.itp"),
        "interpretation": "Neutral ASH306 GROMACS/ParmEd entry preflight only; not a TS, RC, PMF, or barrier.",
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))
    if not passed:
        return 1

    amber_output.mkdir(parents=True, exist_ok=True)
    structure.save(str(amber_output / "system.prmtop"), overwrite=True)
    structure.save(str(amber_output / "system.rst7"), overwrite=True)
    result["system_prmtop_sha256"] = sha256(amber_output / "system.prmtop")
    result["system_rst7_sha256"] = sha256(amber_output / "system.rst7")
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
