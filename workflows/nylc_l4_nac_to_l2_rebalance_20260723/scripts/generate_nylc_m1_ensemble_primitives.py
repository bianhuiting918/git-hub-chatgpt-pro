#!/usr/bin/env python3
"""Generate small, auditable metrics from one fully unrestrained M1 replica."""

from __future__ import annotations

import argparse
import json
import math
import pathlib

import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import calc_bonds, distance_array, minimize_vectors

import analyze_nylc_m1_proton_geometry as proton

GATE_AXIS = np.asarray(
    (0.4904295935403325, 0.813080787402625, -0.3136533866174433),
    dtype=float,
)
GATE_BASELINE_NM = -1.51109
THR267_OG1_INDEX1 = 8960
L2_C_INDEX1 = 10287
L2_O_INDEX1 = 10288
L2_N_INDEX1 = 10289


def parse_ndx(path: pathlib.Path) -> dict[str, list[int]]:
    groups: dict[str, list[int]] = {}
    current = None
    for raw in pathlib.Path(path).read_text(errors="strict").splitlines():
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith("[") and line.endswith("]"):
            current = line[1:-1].strip()
            if not current or current in groups:
                raise ValueError(f"invalid or duplicate NDX group: {line}")
            groups[current] = []
            continue
        if current is None:
            raise ValueError("NDX atom membership appears before a group")
        groups[current].extend(int(value) for value in line.split())
    if not groups:
        raise ValueError(f"{path}: no NDX groups")
    return groups


def gate_opening_nm(vector_nm) -> float:
    vector = np.asarray(vector_nm, dtype=float)
    if vector.shape != (3,) or not np.all(np.isfinite(vector)):
        raise ValueError("gate vector must contain three finite values")
    return float(np.dot(vector, GATE_AXIS) - GATE_BASELINE_NM)


def pocket_retained(contact_count: int, com_nm: float) -> bool:
    return int(contact_count) >= 3 and float(com_nm) <= 1.2


def _identity(atom, index1: int, resname: str, name: str) -> None:
    if atom.index + 1 != index1 or atom.resname != resname or atom.name != name:
        raise ValueError(
            f"atom {index1} maps to {atom.resname}:{atom.name} index1={atom.index + 1}, "
            f"expected {resname}:{name}"
        )


def _write_xvg(path: pathlib.Path, rows) -> None:
    with path.open("w") as handle:
        for row in rows:
            handle.write(" ".join(f"{float(value):.10g}" for value in row) + "\n")


def _minimum_center_vector(left, right, box) -> np.ndarray:
    vector = np.asarray(right.center_of_mass() - left.center_of_mass(), dtype=float)
    return np.asarray(minimize_vectors(vector, box), dtype=float)


def _heavy(group):
    return group.select_atoms("not name H*")


def _local_resname_matches_m1(record: dict, residue) -> bool:
    expected = str(record["resname"])
    observed = str(residue.resname)
    if observed == expected:
        return True
    return (
        expected == "ASP"
        and observed == "ASH"
        and int(record.get("resid", -1)) == 306
        and int(residue.resid) == 306
    )


def _local_groups(universe, selection: dict, thr_atom):
    local_indices = []
    for record in selection["local_rmsd"]["local_residues"]:
        residue_instance = int(record["residue_instance"])
        residue = universe.residues[residue_instance - 1]
        if not _local_resname_matches_m1(record, residue):
            raise ValueError(
                f"residue instance {residue_instance} is {residue.resname}, "
                f"expected {record['resname']}"
            )
        local_indices.extend(_heavy(residue.atoms).indices.tolist())
    ligand_heavy = _heavy(universe.select_atoms("resname L2"))
    local_indices.extend(ligand_heavy.indices.tolist())
    local = universe.atoms[np.asarray(local_indices, dtype=int)]
    expected = int(selection["local_rmsd"]["local_total_heavy_atom_count"])
    if len(local) != expected:
        raise ValueError(f"local heavy atom count {len(local)} != {expected}")

    segment_residues = list(thr_atom.segment.residues)
    thr_position = next(
        index for index, residue in enumerate(segment_residues)
        if residue.ix == thr_atom.residue.ix
    )
    active_residues = segment_residues[thr_position:thr_position + 89]
    alignment_indices = []
    for residue in active_residues:
        selected = residue.atoms.select_atoms("name N CA C")
        if len(selected) != 3:
            raise ValueError(f"alignment residue {residue.resname}{residue.resid} lacks N/CA/C")
        alignment_indices.extend(selected.indices.tolist())
    alignment = universe.atoms[np.asarray(alignment_indices, dtype=int)]
    if len(alignment) != 267:
        raise ValueError(f"alignment atom count {len(alignment)} != 267")
    return alignment, local


def generate(
    tpr: pathlib.Path,
    xtc: pathlib.Path,
    source_cycle_ndx: pathlib.Path,
    selection_manifest: pathlib.Path,
    candidate_id: str,
    velocity_seed: int,
    output_dir: pathlib.Path,
) -> dict:
    output_dir.mkdir(parents=True, exist_ok=False)
    selection = json.loads(selection_manifest.read_text())
    universe = mda.Universe(str(tpr), str(xtc))
    og = universe.atoms[THR267_OG1_INDEX1 - 1]
    carbon = universe.atoms[L2_C_INDEX1 - 1]
    oxygen = universe.atoms[L2_O_INDEX1 - 1]
    nitrogen = universe.atoms[L2_N_INDEX1 - 1]
    _identity(og, THR267_OG1_INDEX1, "THR", "OG1")
    _identity(carbon, L2_C_INDEX1, "L2", "C12")
    _identity(oxygen, L2_O_INDEX1, "L2", "O2")
    _identity(nitrogen, L2_N_INDEX1, "L2", "N3")

    ligand = universe.select_atoms("resname L2")
    ligand_heavy = _heavy(ligand)
    if len(ligand) != 79 or len(ligand_heavy) != 33:
        raise ValueError(f"L2 atom counts are {len(ligand)}/{len(ligand_heavy)}, expected 79/33")
    protein_heavy = _heavy(universe.select_atoms("protein"))
    water_o = universe.select_atoms("resname SOL and name OW")
    if len(water_o) == 0:
        water_o = universe.select_atoms("(resname WAT or resname HOH) and (name O or name OH2)")
    if len(water_o) == 0:
        raise ValueError("no explicit water oxygen atoms")

    groups = parse_ndx(source_cycle_ndx)
    for required in ("Core", "Gate"):
        if required not in groups:
            raise ValueError(f"source NDX lacks {required}")
    core = universe.atoms[np.asarray(groups["Core"], dtype=int) - 1]
    gate = universe.atoms[np.asarray(groups["Gate"], dtype=int) - 1]
    alignment, local = _local_groups(universe, selection, og)

    thr = og.residue
    nalpha = thr.atoms.select_atoms("name N")
    hg = thr.atoms.select_atoms("name HG1")
    if len(nalpha) != 1 or len(hg) != 1:
        raise ValueError("Thr267 Nalpha/HG1 mapping failed")
    segment_residues = list(og.segment.residues)
    thr_position = next(i for i, residue in enumerate(segment_residues) if residue.ix == thr.ix)
    asp306 = segment_residues[thr_position + 39]
    asp308 = segment_residues[thr_position + 41]
    if asp306.resname != "ASH" or asp308.resname != "ASP":
        raise ValueError(f"expected ASH306/ASP308, got {asp306.resname}/{asp308.resname}")
    asp306_o = asp306.atoms.select_atoms("name OD1 OD2")
    asp308_o = asp308.atoms.select_atoms("name OD1 OD2")
    water_h_by_o = {
        int(ow.index): [
            atom for atom in ow.residue.atoms if atom.name.upper().startswith("H")
        ]
        for ow in water_o
    }

    local_residue_indices = sorted(
        {int(item["residue_instance"]) - 1 for item in selection["local_rmsd"]["local_residues"]}
    )
    pocket_atom_indices = []
    for residue_index in local_residue_indices:
        residue = universe.residues[residue_index]
        if residue.resname != "L2":
            pocket_atom_indices.extend(_heavy(residue.atoms).indices.tolist())
    pocket = universe.atoms[np.asarray(pocket_atom_indices, dtype=int)]
    if len(pocket) != int(selection["local_rmsd"]["local_protein_heavy_atom_count"]):
        raise ValueError("pocket heavy atom universe differs from frozen selection")

    distance_rows = []
    angle_rows = []
    gate_rows = []
    pocket_rows = []
    proton_rows = []
    fingerprint_times = []
    alignment_positions = []
    local_positions = []
    minimum_protein_A = math.inf
    minimum_water_A = math.inf

    for ts in universe.trajectory:
        time = float(ts.time)
        box = ts.dimensions
        ogp = og.position.copy()
        hgp = hg[0].position.copy()
        np_ = nalpha[0].position.copy()
        cp = carbon.position.copy()
        op = oxygen.position.copy()
        attack_A = proton.scalar_bond(ogp, cp, box)
        attack_angle = proton.scalar_angle(op, cp, ogp, box)
        nac = bool(attack_A <= 3.5 and 95.0 <= attack_angle <= 115.0)
        distance_rows.append((time, attack_A / 10.0))
        angle_rows.append((time, attack_angle))

        gate_vector_A = _minimum_center_vector(core, gate, box)
        gate_rows.append((time, gate_opening_nm(gate_vector_A / 10.0)))

        pair_distances = distance_array(ligand_heavy.positions, pocket.positions, box=box)
        contact_count = int(np.count_nonzero(pair_distances <= 4.5))
        ligand_pocket_A = float(np.linalg.norm(_minimum_center_vector(pocket, ligand_heavy, box)))
        pocket_rows.append((time, contact_count, ligand_pocket_A / 10.0))
        minimum_protein_A = min(
            minimum_protein_A,
            float(np.min(distance_array(ligand_heavy.positions, protein_heavy.positions, box=box))),
        )
        minimum_water_A = min(
            minimum_water_A,
            float(np.min(distance_array(ligand_heavy.positions, water_o.positions, box=box))),
        )

        routes = {
            "OgH_to_Nalpha_preorg": proton.best_direct(
                ogp, [hgp], np.asarray([np_]), box
            ),
            "OgH_to_Asp306": proton.best_direct(ogp, [hgp], asp306_o.positions, box),
            "OgH_to_Asp308": proton.best_direct(ogp, [hgp], asp308_o.positions, box),
        }
        routes["OgH_to_water"], _ = proton.best_water_acceptor(
            ogp, [hgp], water_o, box
        )
        routes["OgH_via_water_to_Asp306"] = proton.water_bridge(
            ogp, hgp, asp306_o, water_o, water_h_by_o, box
        )
        routes["OgH_via_water_to_Asp308"] = proton.water_bridge(
            ogp, hgp, asp308_o, water_o, water_h_by_o, box
        )
        proton_rows.append(
            proton.sanitize({"time_ps": time, "nac": nac, "channels": routes})
        )
        if nac:
            fingerprint_times.append(time)
            alignment_positions.append(alignment.positions.copy())
            local_positions.append(local.positions.copy())

    _write_xvg(output_dir / "nac_distance.xvg", distance_rows)
    _write_xvg(output_dir / "nac_angle.xvg", angle_rows)
    _write_xvg(output_dir / "gate_opening.xvg", gate_rows)
    _write_xvg(output_dir / "pocket_state.xvg", pocket_rows)
    with (output_dir / "proton_network.jsonl").open("w") as handle:
        for row in proton_rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    (output_dir / "minimum_contact.json").write_text(
        json.dumps(
            {
                "minimum_ligand_protein_heavy_nm": minimum_protein_A / 10.0,
                "minimum_ligand_water_heavy_nm": minimum_water_A / 10.0,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    np.savez_compressed(
        output_dir / "nac_fingerprints.npz",
        time_ps=np.asarray(fingerprint_times, dtype=float),
        alignment_A=np.asarray(alignment_positions, dtype=float).reshape(
            (-1, len(alignment), 3)
        ),
        local_A=np.asarray(local_positions, dtype=float).reshape((-1, len(local), 3)),
    )
    manifest = {
        "schema_version": 1,
        "candidate_id": candidate_id,
        "velocity_seed": int(velocity_seed),
        "microstate": "M1_NalphaH2_Asp306H",
        "gate_definition": "NylC residues 261-266; Thr267 excluded",
        "fully_unrestrained": True,
        "source": {
            "tpr": str(tpr),
            "xtc": str(xtc),
            "tpr_sha256": proton.sha256(tpr),
            "xtc_sha256": proton.sha256(xtc),
            "selection_manifest": str(selection_manifest),
            "selection_manifest_sha256": proton.sha256(selection_manifest),
        },
        "atom_mapping_index1": {
            "thr267_og1": THR267_OG1_INDEX1,
            "l2_c": L2_C_INDEX1,
            "l2_o": L2_O_INDEX1,
            "l2_n": L2_N_INDEX1,
        },
        "frame_count": len(distance_rows),
        "nac_frame_count": len(fingerprint_times),
        "pocket_heavy_atom_count": len(pocket),
        "local_heavy_atom_count": len(local),
        "alignment_atom_count": len(alignment),
        "scientific_scope": "fixed_topology_preorganization_not_proton_transfer",
    }
    (output_dir / "replica_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tpr", type=pathlib.Path, required=True)
    parser.add_argument("--xtc", type=pathlib.Path, required=True)
    parser.add_argument("--source-cycle-ndx", type=pathlib.Path, required=True)
    parser.add_argument("--selection-manifest", type=pathlib.Path, required=True)
    parser.add_argument("--candidate-id", required=True)
    parser.add_argument("--velocity-seed", type=int, required=True)
    parser.add_argument("--output-dir", type=pathlib.Path, required=True)
    args = parser.parse_args()
    manifest = generate(
        args.tpr,
        args.xtc,
        args.source_cycle_ndx,
        args.selection_manifest,
        args.candidate_id,
        args.velocity_seed,
        args.output_dir,
    )
    print(json.dumps(manifest, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
