#!/usr/bin/env python3
"""Fail-closed audit of the selected full-system NylC A1 NAC frame."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import pathlib
from types import SimpleNamespace
from typing import Callable, Iterable, Mapping

import numpy as np

SELECTED_TIME_PS = 354.000
TIME_TOLERANCE_PS = 0.001
COORDINATE_TOLERANCE_NM = 0.0015
NAC_DISTANCE_MAX_NM = 0.35
NAC_ANGLE_MIN_DEG = 95.0
NAC_ANGLE_MAX_DEG = 115.0
SEVERE_CLASH_CUTOFF_NM = 0.18
EXPECTED_ATOM_COUNT = 133589
EXPECTED_L2_ATOM_COUNT = 79
EXPECTED_L2_HEAVY_COUNT = 33
EXPECTED_SOURCE_HASHES = {
    "tpr": "c60078a92c2ace51facde4ef64e453f690177fc88b4d6363427f935944fa2e43",
    "xtc": "1a54f1b5b9f139b746c22d9e0f7e9a4a94eb8154bf2b881888986eedca933d89",
}
REACTIVE_INDEX1 = {
    "thr267_og1": 8960,
    "l2_c": 10287,
    "l2_o": 10288,
    "l2_n": 10289,
}
SCIENTIFIC_SCOPE = "classical_fixed_topology_preorganization_not_proton_transfer"


class AuditError(RuntimeError):
    """A frozen representative-frame contract was not satisfied."""

    def __init__(
        self,
        message: str,
        *,
        stage: str = "audit",
        gate: str = "audit",
    ):
        super().__init__(message)
        self.stage = str(stage)
        self.gate = str(gate)


def _require_frozen_time(time_ps: float) -> float:
    requested = float(time_ps)
    if (
        not math.isfinite(requested)
        or abs(requested - SELECTED_TIME_PS) > TIME_TOLERANCE_PS
    ):
        raise AuditError(
            f"requested time {requested!r} ps differs from frozen selected time "
            f"{SELECTED_TIME_PS:.3f} ps by more than {TIME_TOLERANCE_PS:.3f} ps",
            stage="source_time",
            gate="source_time_unique",
        )
    return SELECTED_TIME_PS


def _sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_source_hashes(
    tpr: pathlib.Path,
    xtc: pathlib.Path,
    declared: Mapping[str, str],
    *,
    hash_file: Callable[[pathlib.Path], str] = _sha256,
) -> dict[str, str]:
    declared_copy = dict(declared)
    if declared_copy != EXPECTED_SOURCE_HASHES:
        raise AuditError(
            "declared source hashes do not equal the frozen TPR/XTC authority"
        )
    observed = {
        "tpr": hash_file(pathlib.Path(tpr)),
        "xtc": hash_file(pathlib.Path(xtc)),
    }
    if observed != EXPECTED_SOURCE_HASHES:
        raise AuditError(
            "observed source hashes do not equal the frozen TPR/XTC authority"
        )
    return observed


def _select_unique_time(
    frames: Iterable[object],
    target_time_ps: float,
    tolerance_ps: float = TIME_TOLERANCE_PS,
):
    matches = [
        frame
        for frame in frames
        if abs(float(frame.time) - float(target_time_ps)) <= float(tolerance_ps)
    ]
    if len(matches) != 1:
        raise AuditError(
            f"expected exactly one frame within {tolerance_ps:.6g} ps of "
            f"{target_time_ps:.6f} ps; found {len(matches)}"
        )
    return matches[0]


def _validated_box(box) -> np.ndarray:
    dimensions = np.asarray(box, dtype=float)
    if dimensions.shape != (6,) or not np.all(np.isfinite(dimensions)):
        raise AuditError("box must contain six finite GROMACS dimensions")
    if np.any(dimensions[:3] <= 0.0):
        raise AuditError("box lengths must be positive")
    if np.any(dimensions[3:] <= 0.0) or np.any(dimensions[3:] >= 180.0):
        raise AuditError("box angles must lie strictly between 0 and 180 degrees")
    return dimensions


def _box_vectors_A(box) -> np.ndarray:
    lx, ly, lz, alpha_deg, beta_deg, gamma_deg = _validated_box(box)
    alpha, beta, gamma = np.deg2rad((alpha_deg, beta_deg, gamma_deg))
    sin_gamma = math.sin(gamma)
    if abs(sin_gamma) < 1.0e-12:
        raise AuditError("box gamma angle is singular")
    a = np.asarray((lx, 0.0, 0.0), dtype=float)
    b = np.asarray((ly * math.cos(gamma), ly * sin_gamma, 0.0), dtype=float)
    cx = lz * math.cos(beta)
    cy = lz * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / sin_gamma
    cz2 = lz * lz - cx * cx - cy * cy
    if cz2 <= 0.0:
        raise AuditError("box vectors are singular")
    c = np.asarray((cx, cy, math.sqrt(cz2)), dtype=float)
    return np.asarray((a, b, c), dtype=float)


def _minimum_image_A(vectors, box) -> np.ndarray:
    values = np.asarray(vectors, dtype=float)
    if values.shape[-1:] != (3,) or not np.all(np.isfinite(values)):
        raise AuditError("coordinate vector must end in three finite values")
    matrix = _box_vectors_A(box)
    fractional = values @ np.linalg.inv(matrix)
    fractional -= np.rint(fractional)
    return fractional @ matrix


def _maximum_pbc_displacement_nm(source_positions_A, extracted_positions_A, box) -> float:
    source = np.asarray(source_positions_A, dtype=float)
    extracted = np.asarray(extracted_positions_A, dtype=float)
    if source.shape != extracted.shape or source.ndim != 2 or source.shape[1] != 3:
        raise AuditError("coordinate identity requires matching N by 3 arrays")
    if not np.all(np.isfinite(source)) or not np.all(np.isfinite(extracted)):
        raise AuditError("coordinate identity requires finite coordinates")
    if len(source) == 0:
        raise AuditError("coordinate identity cannot audit an empty system")
    displacement_A = _minimum_image_A(extracted - source, box)
    return float(np.max(np.linalg.norm(displacement_A, axis=1)) / 10.0)


def _validate_coordinate_identity(source_positions_A, extracted_positions_A, box) -> float:
    maximum_nm = _maximum_pbc_displacement_nm(
        source_positions_A, extracted_positions_A, box
    )
    if maximum_nm > COORDINATE_TOLERANCE_NM:
        raise AuditError(
            f"coordinate identity maximum displacement {maximum_nm:.9g} nm "
            f"exceeds {COORDINATE_TOLERANCE_NM:.9g} nm"
        )
    return maximum_nm


def _validate_box_identity(source_box, extracted_box) -> dict[str, list[float]]:
    source = _validated_box(source_box)
    extracted = _validated_box(extracted_box)
    length_delta_nm = np.max(np.abs(source[:3] - extracted[:3])) / 10.0
    angle_delta_deg = np.max(np.abs(source[3:] - extracted[3:]))
    if length_delta_nm > COORDINATE_TOLERANCE_NM or angle_delta_deg > 0.001:
        raise AuditError(
            "extracted/source box identity exceeds length or angle tolerance"
        )
    return {
        "source_A_deg": [float(value) for value in source],
        "extracted_A_deg": [float(value) for value in extracted],
    }


def _require_identity(
    atom,
    index1: int,
    resname: str,
    name: str,
    resid: int | None = None,
) -> dict[str, object]:
    observed_index1 = int(atom.index) + 1
    observed_resname = str(atom.resname)
    observed_name = str(atom.name)
    observed_resid = int(atom.resid)
    expected_label = f"{resname}{resid if resid is not None else ''}:{name}"
    if (
        observed_index1 != int(index1)
        or observed_resname != resname
        or observed_name != name
        or (resid is not None and observed_resid != int(resid))
    ):
        raise AuditError(
            f"atom index1={observed_index1} is "
            f"{observed_resname}{observed_resid}:{observed_name}; "
            f"expected {expected_label} at index1={index1}"
        )
    return {
        "index1": observed_index1,
        "segid": str(getattr(atom, "segid", "")),
        "resid": observed_resid,
        "resname": observed_resname,
        "name": observed_name,
    }


def _validate_a1_bonds(nalpha, og1) -> dict[str, list[str]]:
    n_bonded_names = [str(value) for value in nalpha.bonded_atoms.names]
    n_hydrogens = sorted(name for name in n_bonded_names if name.upper().startswith("H"))
    if n_hydrogens != ["H1", "H2", "HG1"]:
        raise AuditError(
            f"Thr267 N bonded hydrogens are {n_hydrogens}, expected H1/H2/HG1"
        )
    og1_bonded = sorted(str(value) for value in og1.bonded_atoms.names)
    if og1_bonded != ["CB"]:
        raise AuditError(
            f"Thr267 OG1 bonded atoms are {og1_bonded}, expected only CB "
            "(and therefore no OG1-HG1 bond)"
        )
    return {
        "n_bonded_hydrogens": n_hydrogens,
        "og1_bonded_atoms": og1_bonded,
    }


def _joint_nac(oxygen_A, carbon_A, og1_A, box) -> dict[str, object]:
    oxygen_vector_A = _minimum_image_A(
        np.asarray(oxygen_A, dtype=float) - np.asarray(carbon_A, dtype=float), box
    )
    attack_vector_A = _minimum_image_A(
        np.asarray(og1_A, dtype=float) - np.asarray(carbon_A, dtype=float), box
    )
    oxygen_norm = float(np.linalg.norm(oxygen_vector_A))
    attack_norm = float(np.linalg.norm(attack_vector_A))
    if oxygen_norm <= 0.0 or attack_norm <= 0.0:
        raise AuditError("joint NAC angle is undefined for a zero-length vector")
    cosine = float(
        np.clip(
            np.dot(oxygen_vector_A, attack_vector_A)
            / (oxygen_norm * attack_norm),
            -1.0,
            1.0,
        )
    )
    angle_deg = float(np.degrees(np.arccos(cosine)))
    distance_nm = attack_norm / 10.0
    return {
        "c_og1_distance_nm": distance_nm,
        "o_c_og1_angle_deg": angle_deg,
        "distance_max_nm": NAC_DISTANCE_MAX_NM,
        "angle_range_deg": [NAC_ANGLE_MIN_DEG, NAC_ANGLE_MAX_DEG],
        "joint_pass": bool(
            distance_nm <= NAC_DISTANCE_MAX_NM
            and NAC_ANGLE_MIN_DEG <= angle_deg <= NAC_ANGLE_MAX_DEG
        ),
    }


def _validate_gate_membership(gate) -> list[int]:
    residues = sorted({int(value) for value in gate.resids})
    expected = list(range(261, 267))
    if residues != expected or 267 in residues:
        raise AuditError(
            f"Gate residues are {residues}, expected 261-266 with Thr267 excluded"
        )
    return residues


def _gate_opening_record(
    core,
    gate,
    box,
    *,
    gate_opening: Callable[[np.ndarray], float],
) -> dict[str, object]:
    if len(core) == 0 or len(gate) == 0:
        raise AuditError("Core and Gate groups must both contain atoms")
    gate_resids = _validate_gate_membership(gate)
    vector_A = _minimum_image_A(
        np.asarray(gate.center_of_mass(), dtype=float)
        - np.asarray(core.center_of_mass(), dtype=float),
        box,
    )
    opening_nm = float(gate_opening(vector_A / 10.0))
    if not math.isfinite(opening_nm):
        raise AuditError("gate opening is not finite")
    return {
        "core_atom_count": len(core),
        "gate_atom_count": len(gate),
        "gate_resids": gate_resids,
        "thr267_excluded": True,
        "opening_nm": opening_nm,
    }


def _validate_counts(system_atom_count: int, ligand) -> dict[str, int]:
    if int(system_atom_count) != EXPECTED_ATOM_COUNT:
        raise AuditError(
            f"system atom count {system_atom_count} != {EXPECTED_ATOM_COUNT}"
        )
    l2_atoms = len(ligand)
    l2_heavy = sum(
        1 for name in ligand.names if not str(name).upper().startswith("H")
    )
    if l2_atoms != EXPECTED_L2_ATOM_COUNT or l2_heavy != EXPECTED_L2_HEAVY_COUNT:
        raise AuditError(
            f"L2 atom counts {l2_atoms}/{l2_heavy} != "
            f"{EXPECTED_L2_ATOM_COUNT}/{EXPECTED_L2_HEAVY_COUNT}"
        )
    return {
        "system_atoms": int(system_atom_count),
        "l2_atoms": l2_atoms,
        "l2_heavy_atoms": l2_heavy,
    }


def _atom_record(atom) -> dict[str, object]:
    return {
        "index1": int(atom.index) + 1,
        "segid": str(getattr(atom, "segid", "")),
        "resid": int(atom.resid),
        "resname": str(atom.resname),
        "name": str(atom.name),
    }


def _minimum_contact(ligand, partners, box) -> dict[str, object]:
    if len(ligand) == 0 or len(partners) == 0:
        raise AuditError("minimum contact groups must both contain atoms")
    ligand_positions = np.asarray(ligand.positions, dtype=float)
    partner_positions = np.asarray(partners.positions, dtype=float)
    vectors_A = partner_positions[None, :, :] - ligand_positions[:, None, :]
    distances_A = np.linalg.norm(_minimum_image_A(vectors_A, box), axis=2)
    if not np.all(np.isfinite(distances_A)):
        raise AuditError("minimum contact distances are not finite")
    flat_index = int(np.argmin(distances_A))
    ligand_index, partner_index = np.unravel_index(flat_index, distances_A.shape)
    distance_nm = float(distances_A[ligand_index, partner_index] / 10.0)
    return {
        "distance_nm": distance_nm,
        "severe_clash_cutoff_nm": SEVERE_CLASH_CUTOFF_NM,
        "contact_pass": bool(distance_nm >= SEVERE_CLASH_CUTOFF_NM),
        "ligand": _atom_record(ligand[ligand_index]),
        "partner": _atom_record(partners[partner_index]),
    }


def _validate_minimum_contacts(
    contacts: Mapping[str, Mapping[str, object]],
) -> bool:
    required = ("ligand_protein_heavy", "ligand_water_heavy")
    for label in required:
        if label not in contacts:
            raise AuditError(
                f"minimum-contact result lacks {label}",
                stage="minimum_contacts",
                gate="minimum_contacts",
            )
        record = contacts[label]
        distance_nm = float(record.get("distance_nm", float("nan")))
        cutoff_nm = float(record.get("severe_clash_cutoff_nm", float("nan")))
        passed = bool(record.get("contact_pass", False))
        if (
            not math.isfinite(distance_nm)
            or cutoff_nm != SEVERE_CLASH_CUTOFF_NM
            or not passed
            or distance_nm < SEVERE_CLASH_CUTOFF_NM
        ):
            raise AuditError(
                f"{label} minimum {distance_nm:.9g} nm fails the frozen "
                f"{SEVERE_CLASH_CUTOFF_NM:.9g} nm severe-clash cutoff",
                stage="minimum_contacts",
                gate="minimum_contacts",
            )
    return True


def _parse_ndx(path: pathlib.Path) -> dict[str, list[int]]:
    groups: dict[str, list[int]] = {}
    current: str | None = None
    for raw in pathlib.Path(path).read_text(encoding="utf-8", errors="strict").splitlines():
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith("[") and line.endswith("]"):
            current = line[1:-1].strip()
            if not current or current in groups:
                raise AuditError(f"invalid or duplicate NDX group: {line}")
            groups[current] = []
        elif current is None:
            raise AuditError("NDX atom membership appears before a group")
        else:
            try:
                groups[current].extend(int(value) for value in line.split())
            except ValueError as exc:
                raise AuditError(f"invalid integer in NDX group {current}") from exc
    for required in ("Core", "Gate"):
        if required not in groups or not groups[required]:
            raise AuditError(f"source NDX lacks nonempty {required} group")
    return groups


def _group_from_index1(universe, index1: Iterable[int], label: str):
    values = np.asarray(list(index1), dtype=int)
    if (
        values.ndim != 1
        or len(values) == 0
        or np.any(values < 1)
        or np.any(values > len(universe.atoms))
    ):
        raise AuditError(f"{label} NDX indices fall outside the system")
    return universe.atoms[values - 1]


def _load_existing_gate_opening() -> Callable[[np.ndarray], float]:
    source = pathlib.Path(__file__).with_name("generate_nylc_m1_ensemble_primitives.py")
    spec = importlib.util.spec_from_file_location(
        "_nylc_existing_gate_primitive", source
    )
    if spec is None or spec.loader is None:
        raise AuditError(f"cannot load existing gate primitive from {source}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.gate_opening_nm


def _scientific_status(gates: Mapping[str, bool]) -> str:
    if gates and all(bool(value) for value in gates.values()):
        return "PASS_A1_REPRESENTATIVE_NAC_FRAME"
    return "FAIL_A1_REPRESENTATIVE_NAC_FRAME"


def _heavy(group):
    return group.select_atoms("not name H*")


def audit_frame(
    tpr,
    xtc,
    gro,
    ndx,
    time_ps: float,
    source_hashes: Mapping[str, str],
) -> dict[str, object]:
    frozen_time_ps = _require_frozen_time(time_ps)
    tpr_path = pathlib.Path(tpr)
    xtc_path = pathlib.Path(xtc)
    gro_path = pathlib.Path(gro)
    ndx_path = pathlib.Path(ndx)

    observed_hashes = _validate_source_hashes(
        tpr_path, xtc_path, source_hashes
    )

    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise AuditError("MDAnalysis is required for representative-frame audit") from exc

    source = mda.Universe(str(tpr_path), str(xtc_path))
    frame_records = [
        SimpleNamespace(frame=int(ts.frame), time=float(ts.time))
        for ts in source.trajectory
    ]
    selected = _select_unique_time(frame_records, frozen_time_ps)
    source.trajectory[int(selected.frame)]
    selected_time = float(source.trajectory.ts.time)
    source_positions_A = source.atoms.positions.copy()
    source_box = _validated_box(source.trajectory.ts.dimensions).copy()

    extracted = mda.Universe(str(tpr_path), str(gro_path))
    extracted_positions_A = extracted.atoms.positions.copy()
    extracted_box = _validated_box(extracted.trajectory.ts.dimensions).copy()

    if len(source.atoms) != EXPECTED_ATOM_COUNT:
        raise AuditError(
            f"source system atom count {len(source.atoms)} != {EXPECTED_ATOM_COUNT}"
        )
    if len(extracted.atoms) != len(source.atoms):
        raise AuditError(
            "extracted/source atom counts differ before coordinate identity audit"
        )
    box_record = _validate_box_identity(source_box, extracted_box)
    maximum_displacement_nm = _validate_coordinate_identity(
        source_positions_A, extracted_positions_A, source_box
    )

    expected_identities = {
        "thr267_og1": ("THR", "OG1", 267),
        "l2_c": ("L2", "C12", None),
        "l2_o": ("L2", "O2", None),
        "l2_n": ("L2", "N3", None),
    }
    atom_mapping: dict[str, dict[str, object]] = {}
    for key, index1 in REACTIVE_INDEX1.items():
        resname, name, resid = expected_identities[key]
        source_atom = source.atoms[index1 - 1]
        extracted_atom = extracted.atoms[index1 - 1]
        source_record = _require_identity(source_atom, index1, resname, name, resid)
        extracted_record = _require_identity(
            extracted_atom, index1, resname, name, resid
        )
        if source_record != extracted_record:
            raise AuditError(f"source/extracted identity differs for {key}")
        atom_mapping[key] = source_record

    thr267 = source.atoms[REACTIVE_INDEX1["thr267_og1"] - 1].residue
    nalpha = thr267.atoms.select_atoms("name N")
    if len(nalpha) != 1:
        raise AuditError(f"Thr267 N-alpha selection returned {len(nalpha)} atoms")
    a1_bonds = _validate_a1_bonds(
        nalpha[0], source.atoms[REACTIVE_INDEX1["thr267_og1"] - 1]
    )

    ligand = extracted.select_atoms("resname L2")
    counts = _validate_counts(len(extracted.atoms), ligand)
    ligand_heavy = _heavy(ligand)
    protein_heavy = _heavy(extracted.select_atoms("protein"))
    water_oxygen = extracted.select_atoms("resname SOL and name OW")
    if len(water_oxygen) == 0:
        water_oxygen = extracted.select_atoms(
            "(resname WAT or resname HOH) and (name O or name OH2)"
        )
    if len(water_oxygen) == 0:
        raise AuditError("no explicit water oxygen atoms were found")

    og1 = extracted.atoms[REACTIVE_INDEX1["thr267_og1"] - 1]
    carbon = extracted.atoms[REACTIVE_INDEX1["l2_c"] - 1]
    oxygen = extracted.atoms[REACTIVE_INDEX1["l2_o"] - 1]
    nac = _joint_nac(oxygen.position, carbon.position, og1.position, extracted_box)
    if not nac["joint_pass"]:
        raise AuditError(
            "selected frame fails the joint NAC distance-and-angle requirement"
        )

    ndx_groups = _parse_ndx(ndx_path)
    core = _group_from_index1(extracted, ndx_groups["Core"], "Core")
    gate = _group_from_index1(extracted, ndx_groups["Gate"], "Gate")
    gate_record = _gate_opening_record(
        core,
        gate,
        extracted_box,
        gate_opening=_load_existing_gate_opening(),
    )

    minimum_contacts = {
        "ligand_protein_heavy": _minimum_contact(
            ligand_heavy, protein_heavy, extracted_box
        ),
        "ligand_water_heavy": _minimum_contact(
            ligand_heavy, water_oxygen, extracted_box
        ),
    }
    contacts_pass = _validate_minimum_contacts(minimum_contacts)

    gates = {
        "source_hashes": True,
        "source_time_unique": True,
        "coordinate_identity": True,
        "atom_mapping": True,
        "a1_bonds": True,
        "joint_nac": bool(nac["joint_pass"]),
        "gate_definition": bool(gate_record["thr267_excluded"]),
        "counts_and_box": True,
        "minimum_contacts": contacts_pass,
    }
    result = {
        "schema_version": 1,
        "scientific_status": _scientific_status(gates),
        "scientific_scope": SCIENTIFIC_SCOPE,
        "source": {
            "tpr": str(tpr_path),
            "xtc": str(xtc_path),
            "sha256": observed_hashes,
            "requested_time_ps": float(time_ps),
            "frozen_selected_time_ps": frozen_time_ps,
            "selected_time_ps": selected_time,
            "time_tolerance_ps": TIME_TOLERANCE_PS,
        },
        "extracted_gro": str(gro_path),
        "coordinate_identity": {
            "maximum_pbc_displacement_nm": maximum_displacement_nm,
            "tolerance_nm": COORDINATE_TOLERANCE_NM,
        },
        "box": box_record,
        "counts": counts,
        "atom_mapping": atom_mapping,
        "a1_bonds": a1_bonds,
        "joint_nac": nac,
        "gate": gate_record,
        "minimum_contacts": minimum_contacts,
        "gates": gates,
    }
    if result["scientific_status"] != "PASS_A1_REPRESENTATIVE_NAC_FRAME":
        raise AuditError("representative frame did not satisfy every scientific gate")
    return result


def _failure_result(error: AuditError) -> dict[str, object]:
    return {
        "schema_version": 1,
        "scientific_status": "FAIL_A1_REPRESENTATIVE_NAC_FRAME",
        "scientific_scope": SCIENTIFIC_SCOPE,
        "failed_stage": error.stage,
        "failed_gate": error.gate,
        "error": str(error),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tpr", type=pathlib.Path, required=True)
    parser.add_argument("--xtc", type=pathlib.Path, required=True)
    parser.add_argument("--gro", type=pathlib.Path, required=True)
    parser.add_argument("--ndx", type=pathlib.Path, required=True)
    parser.add_argument("--time-ps", type=float, default=SELECTED_TIME_PS)
    parser.add_argument(
        "--tpr-sha256", default=EXPECTED_SOURCE_HASHES["tpr"]
    )
    parser.add_argument(
        "--xtc-sha256", default=EXPECTED_SOURCE_HASHES["xtc"]
    )
    parser.add_argument("--output", type=pathlib.Path)
    args = parser.parse_args(argv)
    return_code = 0
    try:
        result = audit_frame(
            args.tpr,
            args.xtc,
            args.gro,
            args.ndx,
            args.time_ps,
            {"tpr": args.tpr_sha256, "xtc": args.xtc_sha256},
        )
    except AuditError as error:
        result = _failure_result(error)
        return_code = 1
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered, encoding="utf-8")
    return return_code


if __name__ == "__main__":
    raise SystemExit(main())
