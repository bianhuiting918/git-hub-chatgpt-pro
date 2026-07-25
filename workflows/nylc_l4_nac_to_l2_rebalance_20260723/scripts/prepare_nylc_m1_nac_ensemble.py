#!/usr/bin/env python3
"""Build one immutable NylC M1 candidate from the frozen real-NAC selection."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

from prepare_nylc_nalpha_h2_ash306 import (
    drop_gro_atom,
    export_chain_pdb,
    extract_molecule_itp,
    replace_gro_slice,
    rewrite_chain_itp,
)
from select_nylc_m1_nac_ensemble import read_gro


DEFAULT_REACTIVE_IDENTITIES = {
    "thr267_og1": (267, "THR", "OG1"),
    "l2_c": (356, "L2", "C12"),
    "l2_o": (356, "L2", "O2"),
    "l2_n": (356, "L2", "N3"),
}
MICROSTATE = {
    "Thr267_Nalpha": "NH2_neutral_proxy",
    "Thr267_Ogamma": "OH",
    "Asp306": "ASH_HD2_neutral",
    "Asp308": "ASP_minus",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _parse_itp_atoms(text):
    section = ""
    atoms = []
    for line in text.splitlines():
        match = re.match(r"\s*\[\s*([^]]+?)\s*\]", line)
        if match:
            section = match.group(1).strip().lower()
            continue
        if section != "atoms" or not line.strip() or line.lstrip().startswith((";", "#")):
            continue
        fields = line.split(";", 1)[0].split()
        if len(fields) >= 8 and fields[0].isdigit():
            atoms.append(
                {
                    "nr": int(fields[0]),
                    "type": fields[1],
                    "resid": int(fields[2]),
                    "resname": fields[3],
                    "name": fields[4],
                    "charge": float(fields[6]),
                }
            )
    return atoms


def _gro_atoms_from_text(text):
    lines = text.splitlines()
    atom_count = int(lines[1])
    if len(lines) != atom_count + 3:
        raise ValueError("GRO atom-count mismatch")
    atoms = []
    for index, line in enumerate(lines[2 : 2 + atom_count], 1):
        atoms.append(
            {
                "nr": index,
                "resid": int(line[:5]),
                "resname": line[5:10].strip(),
                "name": line[10:15].strip(),
                "xyz": tuple(
                    float(line[start:stop])
                    for start, stop in ((20, 28), (28, 36), (36, 44))
                ),
            }
        )
    box = tuple(float(value) for value in lines[-1].split())
    return atoms, box


def _is_hydrogen(name):
    letters = "".join(character for character in name if character.isalpha()).upper()
    return bool(letters) and letters[0] == "H"


def _distance(left, right):
    return math.dist(left, right)


def _angle(left, center, right):
    u = np.asarray(left, dtype=float) - np.asarray(center, dtype=float)
    v = np.asarray(right, dtype=float) - np.asarray(center, dtype=float)
    cosine = float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def _source_geometry(atoms, indices, gate):
    og1 = atoms[indices["thr267_og1"] - 1]["coord_nm"]
    carbon = atoms[indices["l2_c"] - 1]["coord_nm"]
    oxygen = atoms[indices["l2_o"] - 1]["coord_nm"]
    distance = _distance(og1, carbon)
    angle = _angle(og1, carbon, oxygen)
    return {
        "attack_distance_nm": distance,
        "attack_angle_deg": angle,
        "joint_nac": (
            distance <= gate["distance_max_nm"]
            and gate["angle_min_deg"] <= angle <= gate["angle_max_deg"]
        ),
    }


def _validate_source(source_gro, authority):
    atoms, _ = read_gro(source_gro)
    identities = authority.get("reactive_atom_identity_m0", DEFAULT_REACTIVE_IDENTITIES)
    observed = {}
    for label, global_index1 in authority["reactive_global_index1_m0"].items():
        atom = atoms[int(global_index1) - 1]
        identity = (atom["resid"], atom["resname"], atom["atom_name"])
        expected = tuple(identities[label])
        if identity != expected:
            raise ValueError(
                f"reactive identity mismatch for {label}: expected {expected}, observed {identity}"
            )
        observed[label] = {
            "global_index1": int(global_index1),
            "resid": atom["resid"],
            "resname": atom["resname"],
            "atom_name": atom["atom_name"],
        }
    geometry = _source_geometry(
        atoms, authority["reactive_global_index1_m0"], authority["nac_gate"]
    )
    if not geometry["joint_nac"]:
        raise ValueError("selected source no longer satisfies the frozen NAC gate")
    return observed, geometry


def _selected_record(source_gro, authority):
    selection_path = Path(authority["ensemble_selection"])
    selection = json.loads(selection_path.read_text())
    source_resolved = source_gro.resolve()
    matches = [
        item
        for item in selection["selected"]
        if Path(item["source_gro"]).resolve() == source_resolved
    ]
    if len(matches) != 1:
        raise ValueError("source GRO is not present in the frozen selection")
    selected = matches[0]
    observed_sha = sha256_file(source_gro)
    if observed_sha != selected["source_gro_sha256"]:
        raise ValueError("source GRO SHA256 disagrees with the frozen selection")
    return selected, observed_sha, selection_path


def _run(command, cwd, stdout_path, stderr_path, input_text=None):
    result = subprocess.run(
        [str(item) for item in command],
        cwd=cwd,
        input=input_text,
        text=True,
        capture_output=True,
    )
    stdout_path.write_text(result.stdout)
    stderr_path.write_text(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed with exit code {result.returncode}: {' '.join(map(str, command))}"
        )
    return result


def _copy_legacy_assets(legacy_build_root, build_dir):
    required = [
        "topol.top",
        "PA66_L2_GMX.itp",
        "source_cycle.ndx",
        "topol_Protein_chain_A.itp",
        "topol_Protein_chain_D.itp",
        "topol_Protein_chain_E.itp",
    ]
    for name in required:
        source = legacy_build_root / name
        if not source.is_file():
            raise FileNotFoundError(f"missing audited legacy build asset: {source}")
        shutil.copy2(source, build_dir / name)
    patterns = [
        "posre_Protein_chain_A*.itp",
        "posre_Protein_chain_D*.itp",
        "posre_Protein_chain_E*.itp",
        "posre_l2_*.itp",
    ]
    for pattern in patterns:
        matches = list(legacy_build_root.glob(pattern))
        if not matches:
            raise FileNotFoundError(f"missing audited legacy build assets: {pattern}")
        for source in matches:
            shutil.copy2(source, build_dir / source.name)


def _box_matrix(box):
    if len(box) == 3:
        return np.diag(np.asarray(box, dtype=float))
    return np.asarray(
        [
            [box[0], box[3], box[4]],
            [box[5], box[1], box[6]],
            [box[7], box[8], box[2]],
        ],
        dtype=float,
    )


def _minimum_cross_distance(left, right, box, chunk_size=16):
    matrix = _box_matrix(box)
    inverse = np.linalg.inv(matrix)
    minimum = math.inf
    for start in range(0, len(left), chunk_size):
        delta = left[start : start + chunk_size, None, :] - right[None, :, :]
        fractional = delta @ inverse
        fractional -= np.round(fractional)
        delta = fractional @ matrix
        block_minimum = float(np.sqrt(np.sum(delta * delta, axis=2)).min())
        minimum = min(minimum, block_minimum)
    return minimum


def _cross_contact_audit(full_atoms, box, chain_start, chain_end):
    coords = np.asarray([atom["xyz"] for atom in full_atoms], dtype=float)
    chain_idx0 = np.arange(chain_start - 1, chain_end, dtype=int)
    rest_idx0 = np.concatenate(
        (np.arange(0, chain_start - 1, dtype=int), np.arange(chain_end, len(full_atoms), dtype=int))
    )
    minimum_all = _minimum_cross_distance(coords[chain_idx0], coords[rest_idx0], box)
    chain_heavy = np.asarray(
        [index for index in chain_idx0 if not _is_hydrogen(full_atoms[index]["name"])],
        dtype=int,
    )
    rest_heavy = np.asarray(
        [index for index in rest_idx0 if not _is_hydrogen(full_atoms[index]["name"])],
        dtype=int,
    )
    minimum_heavy = _minimum_cross_distance(
        coords[chain_heavy], coords[rest_heavy], box, chunk_size=32
    )
    return minimum_all, minimum_heavy


def _position_restraint_text(atoms, force_constant):
    lines = [
        "; generated for NylC M1 NalphaH2/Asp306H active chain",
        "[ position_restraints ]",
        "; atom type fx fy fz",
    ]
    lines.extend(
        f"{atom['nr']:6d} 1 {force_constant} {force_constant} {force_constant}"
        for atom in atoms
        if not _is_hydrogen(atom["name"])
    )
    return "\n".join(lines) + "\n"


def _materialize_candidate(source_gro, output_dir, authority, selected):
    chain_gen = output_dir / "chain_gen"
    build_dir = output_dir / "build"
    chain_gen.mkdir()
    build_dir.mkdir()

    active_start, active_end = authority["active_chain_global_index1_range_m0"]
    source_text = source_gro.read_text()
    pdb_text, export_audit = export_chain_pdb(
        source_text, active_start, active_end, ash_resid=306
    )
    if export_audit["chain_atom_count"] != 1324:
        raise ValueError("active-chain source atom count is not 1324")
    chain_pdb = chain_gen / "chainH_ASH306.pdb"
    chain_pdb.write_text(pdb_text)
    (chain_gen / "export_audit.json").write_text(
        json.dumps(export_audit, indent=2, sort_keys=True) + "\n"
    )

    gmx = authority["gmx"]
    _run(
        [
            gmx,
            "pdb2gmx",
            "-f",
            "chainH_ASH306.pdb",
            "-o",
            "chainH_ASH306.gro",
            "-p",
            "chainH.top",
            "-i",
            "posre_chainH.itp",
            "-ff",
            "amber99sb-ildn",
            "-water",
            "tip3p",
            "-ignh",
        ],
        chain_gen,
        chain_gen / "pdb2gmx.stdout",
        chain_gen / "pdb2gmx.stderr",
    )

    legacy_build_root = Path(authority["legacy_build_root"])
    _copy_legacy_assets(legacy_build_root, build_dir)

    rewritten, itp_rewrite_audit = rewrite_chain_itp(
        (chain_gen / "chainH.top").read_text()
    )
    chain_itp = extract_molecule_itp(rewritten).replace(
        "posre_chainH.itp", "posre_Protein_chain_H.itp"
    )
    if "POSRES_WEAK" not in chain_itp:
        chain_itp += (
            "\n#ifdef POSRES_WEAK\n"
            '#include "posre_Protein_chain_H_weak.itp"\n'
            "#endif\n"
        )
    chain_text, gro_rewrite_audit = drop_gro_atom(
        (chain_gen / "chainH_ASH306.gro").read_text()
    )
    full_text, merge_audit = replace_gro_slice(
        source_text, chain_text, active_start, active_end
    )

    new_atoms = _parse_itp_atoms(chain_itp)
    old_atoms = _parse_itp_atoms(
        (legacy_build_root / "topol_Protein_chain_H.itp").read_text()
    )
    chain_atoms, _ = _gro_atoms_from_text(chain_text)
    if len(new_atoms) != len(chain_atoms) or len(new_atoms) != 1324:
        raise ValueError("M1 active-chain atom count is not 1324")
    if [
        (atom["nr"], atom["resid"], atom["resname"], atom["name"])
        for atom in new_atoms
    ] != [
        (atom["nr"], atom["resid"], atom["resname"], atom["name"])
        for atom in chain_atoms
    ]:
        raise ValueError("M1 active-chain topology/GRO atom identities differ")

    old_charge = sum(atom["charge"] for atom in old_atoms)
    new_charge = sum(atom["charge"] for atom in new_atoms)
    if abs(old_charge + 4.0) > 1e-5 or abs(new_charge + 4.0) > 1e-5:
        raise ValueError("active-chain charge is not -4.0 e")
    charge_delta = new_charge - old_charge
    if abs(charge_delta) > 1e-6:
        raise ValueError("M1 charge transfer changed the full-system charge")

    ordinary_old = {
        (atom["resid"], atom["name"]): (atom["type"], atom["charge"])
        for atom in old_atoms
        if atom["resid"] not in (267, 306)
    }
    ordinary_new = {
        (atom["resid"], atom["name"]): (atom["type"], atom["charge"])
        for atom in new_atoms
        if atom["resid"] not in (267, 306)
    }
    if ordinary_old != ordinary_new:
        raise ValueError("ordinary active-chain topology atoms changed")

    thr267 = [atom for atom in new_atoms if atom["resid"] == 267]
    asp306 = [atom for atom in new_atoms if atom["resid"] == 306]
    asp308 = [atom for atom in new_atoms if atom["resid"] == 308]
    thr_names = {atom["name"] for atom in thr267}
    asp306_names = {atom["name"] for atom in asp306}
    asp308_names = {atom["name"] for atom in asp308}
    if not {"N", "H1", "H2", "OG1", "HG1"}.issubset(thr_names) or "H3" in thr_names:
        raise ValueError("Thr267 is not the NalphaH2/OgammaH microstate")
    if not {"OD1", "OD2", "HD2"}.issubset(asp306_names) or {
        atom["resname"] for atom in asp306
    } != {"ASH"}:
        raise ValueError("Asp306 is not ASH with HD2")
    if not {"OD1", "OD2"}.issubset(asp308_names) or "HD2" in asp308_names:
        raise ValueError("Asp308 is not the deprotonated ASP control state")

    source_atoms, _ = _gro_atoms_from_text(source_text)
    merged_atoms, merged_box = _gro_atoms_from_text(full_text)
    if len(source_atoms) != len(merged_atoms):
        raise ValueError("full-system atom count changed")
    for source_atom, merged_atom in zip(
        source_atoms[: active_start - 1], merged_atoms[: active_start - 1]
    ):
        if source_atom != merged_atom:
            raise ValueError("atoms before active-chain slice changed")
    for source_atom, merged_atom in zip(
        source_atoms[active_end:], merged_atoms[active_end:]
    ):
        if source_atom != merged_atom:
            raise ValueError("atoms after active-chain slice changed")

    source_heavy = {
        (atom["resid"], atom["name"]): atom["xyz"]
        for atom in source_atoms[active_start - 1 : active_end]
        if not _is_hydrogen(atom["name"])
    }
    new_heavy = {
        (atom["resid"], atom["name"]): atom["xyz"]
        for atom in chain_atoms
        if not _is_hydrogen(atom["name"])
    }
    ordinary_keys = set(source_heavy) & set(new_heavy) - {(355, "OC1"), (355, "OC2")}
    maximum_heavy_displacement = max(
        max(
            abs(source_value - new_value)
            for source_value, new_value in zip(source_heavy[key], new_heavy[key])
        )
        for key in ordinary_keys
    )
    if maximum_heavy_displacement > 0.0011:
        raise ValueError("ordinary active-chain heavy atoms moved by more than 0.0011 nm")

    def active_global_index(resid, name):
        atom = next(
            atom for atom in new_atoms if atom["resid"] == resid and atom["name"] == name
        )
        return active_start - 1 + atom["nr"]

    built_indices = {
        "thr267_og1": active_global_index(267, "OG1"),
        "l2_c": authority["reactive_global_index1_m0"]["l2_c"],
        "l2_o": authority["reactive_global_index1_m0"]["l2_o"],
        "l2_n": authority["reactive_global_index1_m0"]["l2_n"],
    }
    built_geometry = {
        "attack_distance_nm": _distance(
            merged_atoms[built_indices["thr267_og1"] - 1]["xyz"],
            merged_atoms[built_indices["l2_c"] - 1]["xyz"],
        ),
        "attack_angle_deg": _angle(
            merged_atoms[built_indices["thr267_og1"] - 1]["xyz"],
            merged_atoms[built_indices["l2_c"] - 1]["xyz"],
            merged_atoms[built_indices["l2_o"] - 1]["xyz"],
        ),
    }
    gate = authority["nac_gate"]
    built_geometry["joint_nac"] = (
        built_geometry["attack_distance_nm"] <= gate["distance_max_nm"]
        and gate["angle_min_deg"]
        <= built_geometry["attack_angle_deg"]
        <= gate["angle_max_deg"]
    )
    if not built_geometry["joint_nac"]:
        raise ValueError("built M1 source no longer satisfies the NAC gate")

    (build_dir / "topol_Protein_chain_H.itp").write_text(chain_itp)
    (build_dir / "system_M1.gro").write_text(full_text)
    (build_dir / "posre_Protein_chain_H.itp").write_text(
        _position_restraint_text(new_atoms, 1000)
    )
    (build_dir / "posre_Protein_chain_H_weak.itp").write_text(
        _position_restraint_text(new_atoms, 100)
    )

    minimum_all, minimum_heavy = _cross_contact_audit(
        merged_atoms, merged_box, active_start, active_end
    )
    if minimum_all < 0.08:
        raise ValueError("minimum active-chain/rest distance is below 0.08 nm")
    if minimum_heavy < 0.18:
        raise ValueError("minimum active-chain/rest heavy-atom contact is below 0.18 nm")

    _run(
        [
            gmx,
            "grompp",
            "-f",
            authority["em_mdp"],
            "-c",
            "system_M1.gro",
            "-r",
            "system_M1.gro",
            "-p",
            "topol.top",
            "-o",
            "preflight.tpr",
            "-po",
            "preflight.processed.mdp",
            "-maxwarn",
            "0",
        ],
        build_dir,
        build_dir / "grompp.stdout",
        build_dir / "grompp.stderr",
    )

    return {
        "charge_delta_e": round(charge_delta, 8),
        "old_chain_charge_e": round(old_charge, 8),
        "new_chain_charge_e": round(new_charge, 8),
        "active_chain_atom_count": len(new_atoms),
        "full_atom_count": len(merged_atoms),
        "full_atom_count_unchanged": len(source_atoms) == len(merged_atoms),
        "reaction_geometry": built_geometry,
        "source_selected_time_ps": selected["time_ps"],
        "atom_mapping_global_index1": built_indices,
        "minimum_chain_rest_distance_nm": minimum_all,
        "minimum_heavy_atom_contact_nm": minimum_heavy,
        "maximum_ordinary_heavy_displacement_nm": maximum_heavy_displacement,
        "grompp_maxwarn": 0,
        "itp_rewrite": itp_rewrite_audit,
        "gro_rewrite": gro_rewrite_audit,
        "merge": merge_audit,
        "sha256": {
            "system_M1_gro": sha256_file(build_dir / "system_M1.gro"),
            "chainH_itp": sha256_file(build_dir / "topol_Protein_chain_H.itp"),
            "preflight_tpr": sha256_file(build_dir / "preflight.tpr"),
        },
    }


def build_candidate(source_gro: Path, candidate_dir: Path, authority: dict) -> dict:
    source_gro = Path(source_gro)
    candidate_dir = Path(candidate_dir)
    if candidate_dir.exists():
        raise FileExistsError(f"refusing to overwrite candidate directory: {candidate_dir}")
    selected, source_sha, selection_path = _selected_record(source_gro, authority)
    observed_identities, source_geometry = _validate_source(source_gro, authority)

    candidate_dir.mkdir(parents=True)
    try:
        materialized = _materialize_candidate(
            source_gro, candidate_dir, authority, selected
        )
        if abs(materialized["charge_delta_e"]) > 1e-6:
            raise ValueError("candidate build charge delta is not zero")
        if not materialized["full_atom_count_unchanged"]:
            raise ValueError("candidate build changed the full atom count")
        if not materialized["reaction_geometry"]["joint_nac"]:
            raise ValueError("candidate build does not retain source NAC geometry")
        if materialized["minimum_chain_rest_distance_nm"] < 0.08:
            raise ValueError("candidate build fails the chain/rest contact gate")
        if materialized["minimum_heavy_atom_contact_nm"] < 0.18:
            raise ValueError("candidate build fails the heavy-atom contact gate")
        if materialized["grompp_maxwarn"] != 0:
            raise ValueError("candidate build did not pass grompp -maxwarn 0")

        audit = {
            "schema_version": 1,
            "status": "PASS_TECHNICAL_M1_BUILD_GROMPP",
            "scientific_status": "M1_SOURCE_BUILD_READY_FOR_EM_NOT_MD_VALIDATED",
            "candidate_id": selected["candidate_id"],
            "event_id": selected["event_id"],
            "source_selected_time_ps": selected["time_ps"],
            "source_gro": str(source_gro),
            "source_gro_sha256": source_sha,
            "selection_manifest": str(selection_path),
            "selection_manifest_sha256": sha256_file(selection_path),
            "microstate": MICROSTATE,
            "gate_definition": "NylC residues 261-266; Thr267 excluded",
            "source_reactive_atom_identities": observed_identities,
            "source_reaction_geometry": source_geometry,
            **materialized,
        }
        (candidate_dir / "build_audit.json").write_text(
            json.dumps(audit, indent=2, sort_keys=True) + "\n"
        )
        (candidate_dir / "PASS.json").write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "status": "PASS_TECHNICAL_M1_BUILD_GROMPP",
                    "candidate_id": selected["candidate_id"],
                    "build_audit": str(candidate_dir / "build_audit.json"),
                    "next_gate": "double-precision flexible-water EM",
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        return audit
    except Exception as error:
        (candidate_dir / "FAIL.json").write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "status": "FAIL_TECHNICAL_BUILD",
                    "candidate_id": selected["candidate_id"],
                    "error_type": type(error).__name__,
                    "reason": str(error),
                    "outputs_preserved": True,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        raise


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--authority", type=Path, required=True)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--candidate-index", type=int, required=True)
    parser.add_argument("--candidate-dir", type=Path, required=True)
    parser.add_argument("--legacy-build-root", type=Path, required=True)
    parser.add_argument("--gmx", required=True)
    parser.add_argument("--em-mdp", type=Path, required=True)
    args = parser.parse_args()

    authority = json.loads(args.authority.read_text())
    authority["ensemble_selection"] = str(args.selection)
    authority["legacy_build_root"] = str(args.legacy_build_root)
    authority["gmx"] = args.gmx
    authority["em_mdp"] = str(args.em_mdp)
    authority.setdefault("reactive_atom_identity_m0", DEFAULT_REACTIVE_IDENTITIES)

    selection = json.loads(args.selection.read_text())
    selected = selection["selected"]
    if args.candidate_index < 0 or args.candidate_index >= len(selected):
        print("candidate index is outside the frozen selection", file=sys.stderr)
        return 2
    record = selected[args.candidate_index]
    try:
        audit = build_candidate(
            Path(record["source_gro"]), args.candidate_dir, authority
        )
    except (FileExistsError, FileNotFoundError, RuntimeError, ValueError) as error:
        print(str(error), file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "status": audit["status"],
                "candidate_id": audit["candidate_id"],
                "candidate_dir": str(args.candidate_dir),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
