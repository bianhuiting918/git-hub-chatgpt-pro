#!/usr/bin/env python3
"""Prepare and audit the minimal NylC A1 Step2 QM-water A2 preflight.

Scope: deterministic attacking-water selection, fixed Step2 QM contract, and
two independent unrestrained A2 legs.  This file never constructs a hydrolysis
product, transition state, path, PMF, barrier, or mechanism claim.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import pathlib
import re
import warnings
from typing import Any, Iterable, Mapping, Sequence

HERE = pathlib.Path(__file__).resolve().parent
ACYL_PATH = HERE / "prepare_audit_nylc_a1_step1_acyl_endpoint_stability.py"
_SPEC = importlib.util.spec_from_file_location("_nylc_a1_acyl_authority", ACYL_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"cannot load frozen acyl authority {ACYL_PATH}")
ACYL = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(ACYL)

TASK_ROOT = pathlib.Path(
    "/work/home/acshdt1dks/nylon_pa66_scnet_20260708/"
    "l4_nac_to_l2_rebalance_20260723"
)
PRMTOP = ACYL.PRMTOP
EXPECTED_PRMTOP_SHA256 = "a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0"
AUTHORITY_COMMIT = "19f327b422fd5081b119e4ff4c1883edefd85ba4"
SENSITIVITY_AUDIT_RELATIVE = pathlib.Path(
    "audit/nylc_a1_acyl_endpoint_connectivity_sensitivity_62216380_v1.json"
)
EXPECTED_SENSITIVITY_AUDIT_BLOB_SHA1 = "ca4fd713540e2621391e2f5219e023d7aa3fafdf"
SOURCE_ROOT = (
    TASK_ROOT
    / "a1_activated_nac_20260726/qmmm/"
    "a1_step1_acyl_release_md_continuation"
)
SOURCES = (
    {
        "seed": "seed26723",
        "attempt": "attempt_62216380_0",
        "restart_sha256": "5d8f76d2c90e3e8c707b640c55a93938f53e18dc30d6f93d54adda467da25f41",
        "velocity_seeds": (2672301, 2672302),
    },
    {
        "seed": "seed26737",
        "attempt": "attempt_62216380_1",
        "restart_sha256": "4cc60ad4d7be099b3f76040f1ee8c49b91deecebb20ed98570bc192b3d511132",
        "velocity_seeds": (2673701, 2673702),
    },
)
REACTIVE = {
    "nalpha": 8949,
    "og1": 8960,
    "hg1": 8961,
    "c11": 10286,
    "c12": 10287,
    "o2": 10288,
    "n3": 10289,
}
WATER_NAMES = {"SOL", "WAT", "HOH", "TIP3", "TIP3P"}
WATER_FILTER = {
    "c12_ow_A": (2.40, 3.50),
    "o2_c12_ow_deg": (95.0, 125.0),
    "h_nalpha_A_max": 2.40,
    "ow_h_nalpha_deg_min": 135.0,
    "ow_nonbonded_heavy_A_min": 1.80,
    "h_nonbonded_heavy_A_min": 1.10,
}
RANK_TARGET = {
    "c12_ow_A": (2.90, 0.45),
    "o2_c12_ow_deg": (107.0, 12.0),
    "h_nalpha_A": (1.90, 0.35),
    "ow_h_nalpha_deg": (180.0, 20.0),
}
EXPECTED_CONTRACT = {
    "qm_atom_count": 149,
    "qmcharge": 0,
    "spin": 1,
    "link_atom_count": 6,
    "electron_count": 518,
}
MD_STEPS = 1000
MD_DT_PS = 0.0005
MD_NTWX = 20
EXPECTED_FRAMES = MD_STEPS // MD_NTWX
TAIL_FRAMES = 20
A2_GATE = {
    "acyl_sensitivity_occupancy_min": 0.80,
    "neutral_water_occupancy_min": 0.80,
    "no_addition_occupancy_min": 1.00,
    "attack_residence_occupancy_min": 0.50,
    "ow_h_A_max": 1.25,
    "c12_ow_no_addition_A_min": 1.80,
    "resident_c12_ow_A_max": 3.80,
    "resident_angle_deg": (85.0, 130.0),
    "resident_h_nalpha_A_max": 2.70,
    "resident_hbond_angle_deg_min": 120.0,
}
NEXT = "STEP2_PRODUCT_ENDPOINT_BLOCKED_PENDING_A2_PASS"
SCIENTIFIC_SCOPE = (
    "EXPLORATORY_A2_QM_WATER_PREFLIGHT_ONLY_NOT_PRODUCT_TS_PATH_PMF_BARRIER_OR_MECHANISM"
)
HARD_PATTERNS = {
    "sander_bomb": r"SANDER BOMB",
    "segmentation": r"segmentation",
    "forrtl": r"forrtl",
    "nan": r"\bnan\b",
    "fatal": r"FATAL",
}


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_blob_sha1(path: pathlib.Path) -> str:
    payload = path.read_bytes()
    header = f"blob {len(payload)}\0".encode("ascii")
    return hashlib.sha1(header + payload).hexdigest()


def write_json(path: pathlib.Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def source_for_index(index: int) -> dict[str, Any]:
    if index not in (0, 1):
        raise ValueError("seed index must be 0 or 1")
    source = dict(SOURCES[index])
    root = SOURCE_ROOT / source["attempt"]
    source.update(
        {
            "root": root,
            "restart": root / "release_md_endpoint.rst7",
            "manifest": root / "ENDPOINT_MANIFEST.json",
            "result": root / "RESULT.json",
        }
    )
    return source


def terminal_payload(
    seed: str, classification: str, reason: str, detail: Mapping[str, Any] | None = None
) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "status": "NOT_EVALUATED_STEP2_A2_QMWATER_PREFLIGHT",
        "scientific_scope": SCIENTIFIC_SCOPE,
        "seed": seed,
        "classification": classification,
        "reason": reason,
        "detail": dict(detail or {}),
        "NEXT": NEXT,
        "product_endpoint_implemented": False,
    }


def write_not_evaluated(
    output: pathlib.Path,
    seed: str,
    classification: str,
    reason: str,
    detail: Mapping[str, Any] | None = None,
) -> None:
    payload = terminal_payload(seed, classification, reason, detail)
    write_json(output / "RESULT.json", payload)
    write_json(output / "NOT_EVALUATED.json", payload)


def _xyz(value: Any) -> tuple[float, float, float]:
    if hasattr(value, "xx"):
        return float(value.xx), float(value.xy), float(value.xz)
    return tuple(float(component) for component in value)  # type: ignore[return-value]


def _box_lengths(structure: Any) -> tuple[float, float, float]:
    if structure.box is None or len(structure.box) < 6:
        raise ValueError("source restart lacks a periodic box")
    box = tuple(float(value) for value in structure.box[:6])
    if any(abs(angle - 90.0) > 1.0e-3 for angle in box[3:]):
        raise ValueError(f"only orthorhombic minimum-image selection is supported: {box}")
    if any(length <= 0.0 for length in box[:3]):
        raise ValueError(f"invalid box lengths: {box[:3]}")
    return box[0], box[1], box[2]


def _delta(
    left: Any, right: Any, box: Sequence[float]
) -> tuple[float, float, float]:
    a, b = _xyz(left), _xyz(right)
    values = []
    for x, y, length in zip(a, b, box):
        value = x - y
        value -= round(value / length) * length
        values.append(value)
    return values[0], values[1], values[2]


def _norm(vector: Sequence[float]) -> float:
    return math.sqrt(sum(value * value for value in vector))


def _distance(left: Any, right: Any, box: Sequence[float]) -> float:
    return _norm(_delta(left, right, box))


def _angle_from_vectors(left: Sequence[float], right: Sequence[float]) -> float:
    denominator = _norm(left) * _norm(right)
    if denominator == 0.0:
        raise ValueError("undefined angle from coincident atoms")
    cosine = sum(x * y for x, y in zip(left, right)) / denominator
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def _angle(first: Any, center: Any, third: Any, box: Sequence[float]) -> float:
    return _angle_from_vectors(
        _delta(first, center, box), _delta(third, center, box)
    )


def _atomic_number(atom: Any) -> int:
    return int(getattr(atom, "atomic_number", 0) or 0)


def _complete_h2o(residue: Any) -> tuple[Any, tuple[Any, Any]] | None:
    if str(residue.name).upper() not in WATER_NAMES or len(residue.atoms) != 3:
        return None
    oxygen = [atom for atom in residue.atoms if _atomic_number(atom) == 8]
    hydrogens = [atom for atom in residue.atoms if _atomic_number(atom) == 1]
    if len(oxygen) != 1 or len(hydrogens) != 2:
        return None
    partners = {atom.idx for atom in oxygen[0].bond_partners}
    if {atom.idx for atom in hydrogens} != partners:
        return None
    return oxygen[0], tuple(sorted(hydrogens, key=lambda atom: atom.idx))


def _water_hbond(
    oxygen: Any, hydrogens: Sequence[Any], nalpha: Any, box: Sequence[float]
) -> tuple[Any, float, float]:
    measured = []
    for hydrogen in hydrogens:
        distance = _distance(hydrogen, nalpha, box)
        angle = _angle(oxygen, hydrogen, nalpha, box)
        measured.append((distance, -angle, hydrogen.idx, hydrogen, angle))
    measured.sort(key=lambda item: (item[0], item[1], item[2]))
    selected = measured[0]
    return selected[3], selected[0], selected[4]


def _no_clash(
    structure: Any,
    oxygen: Any,
    hydrogens: Sequence[Any],
    box: Sequence[float],
) -> tuple[bool, float, float]:
    own = {oxygen.idx, *(atom.idx for atom in hydrogens)}
    heavy = [
        atom
        for atom in structure.atoms
        if atom.idx not in own and _atomic_number(atom) > 1
    ]
    min_oxygen = min(_distance(oxygen, atom, box) for atom in heavy)
    min_hydrogen = min(
        _distance(hydrogen, atom, box)
        for hydrogen in hydrogens
        for atom in heavy
    )
    return (
        min_oxygen >= WATER_FILTER["ow_nonbonded_heavy_A_min"]
        and min_hydrogen >= WATER_FILTER["h_nonbonded_heavy_A_min"],
        min_oxygen,
        min_hydrogen,
    )


def select_water(structure: Any) -> tuple[dict[str, int], list[dict[str, Any]]]:
    box = _box_lengths(structure)
    c12 = structure.atoms[REACTIVE["c12"] - 1]
    o2 = structure.atoms[REACTIVE["o2"] - 1]
    nalpha = structure.atoms[REACTIVE["nalpha"] - 1]
    funnel = {
        "all_residues": len(structure.residues),
        "water_named_residues": 0,
        "complete_h2o": 0,
        "c12_ow_distance": 0,
        "buergi_dunitz_angle": 0,
        "nalpha_hbond": 0,
        "no_hard_clash": 0,
    }
    selected: list[dict[str, Any]] = []
    for residue in structure.residues:
        if str(residue.name).upper() not in WATER_NAMES:
            continue
        funnel["water_named_residues"] += 1
        parsed = _complete_h2o(residue)
        if parsed is None:
            continue
        funnel["complete_h2o"] += 1
        oxygen, hydrogens = parsed
        r = _distance(c12, oxygen, box)
        if not WATER_FILTER["c12_ow_A"][0] <= r <= WATER_FILTER["c12_ow_A"][1]:
            continue
        funnel["c12_ow_distance"] += 1
        theta = _angle(o2, c12, oxygen, box)
        if not (
            WATER_FILTER["o2_c12_ow_deg"][0]
            <= theta
            <= WATER_FILTER["o2_c12_ow_deg"][1]
        ):
            continue
        funnel["buergi_dunitz_angle"] += 1
        donor_h, h_n, phi = _water_hbond(oxygen, hydrogens, nalpha, box)
        if (
            h_n > WATER_FILTER["h_nalpha_A_max"]
            or phi < WATER_FILTER["ow_h_nalpha_deg_min"]
        ):
            continue
        funnel["nalpha_hbond"] += 1
        clear, min_ow, min_hw = _no_clash(
            structure, oxygen, hydrogens, box
        )
        if not clear:
            continue
        funnel["no_hard_clash"] += 1
        score = (
            ((r - RANK_TARGET["c12_ow_A"][0]) / RANK_TARGET["c12_ow_A"][1]) ** 2
            + (
                (theta - RANK_TARGET["o2_c12_ow_deg"][0])
                / RANK_TARGET["o2_c12_ow_deg"][1]
            )
            ** 2
            + (
                (h_n - RANK_TARGET["h_nalpha_A"][0])
                / RANK_TARGET["h_nalpha_A"][1]
            )
            ** 2
            + (
                (phi - RANK_TARGET["ow_h_nalpha_deg"][0])
                / RANK_TARGET["ow_h_nalpha_deg"][1]
            )
            ** 2
        )
        selected.append(
            {
                "residue_index1": residue.idx + 1,
                "residue_name": residue.name,
                "oxygen_index1": oxygen.idx + 1,
                "hydrogen_indices1": [atom.idx + 1 for atom in hydrogens],
                "donor_h_index1": donor_h.idx + 1,
                "c12_ow_A": r,
                "o2_c12_ow_deg": theta,
                "h_nalpha_A": h_n,
                "ow_h_nalpha_deg": phi,
                "minimum_ow_other_heavy_A": min_ow,
                "minimum_hw_other_heavy_A": min_hw,
                "score": score,
            }
        )
    selected.sort(
        key=lambda item: (
            item["score"],
            item["residue_index1"],
            item["oxygen_index1"],
        )
    )
    return funnel, selected


def _parse_qmmask(qmmask: str, expected_count: int = 146) -> list[int]:
    values: list[int] = []
    for token in qmmask.split(","):
        token = token.strip()
        if not token:
            continue
        if token.startswith("@"):
            token = token[1:]
        if not token.isdigit():
            raise ValueError(f"Step1 qmmask is not an explicit atom list: {token!r}")
        values.append(int(token))
    if len(values) != expected_count or len(set(values)) != expected_count:
        raise ValueError(
            f"qmmask is not an explicit unique {expected_count}-atom list"
        )
    return values


def _boundary_bonds(structure: Any, qm_indices1: Iterable[int]) -> list[dict[str, int]]:
    qm = {int(index) - 1 for index in qm_indices1}
    boundaries = []
    for bond in structure.bonds:
        left, right = bond.atom1.idx, bond.atom2.idx
        if (left in qm) != (right in qm):
            q = left if left in qm else right
            m = right if left in qm else left
            boundaries.append({"qm_atom": q + 1, "mm_atom": m + 1})
    return sorted(boundaries, key=lambda item: (item["qm_atom"], item["mm_atom"]))


def _load_json(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def validate_authority(
    source: Mapping[str, Any], code_root: pathlib.Path
) -> tuple[dict[str, Any], dict[str, Any]]:
    for key in ("restart", "manifest", "result"):
        path = pathlib.Path(source[key])
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"missing source {key}: {path}")
    restart = pathlib.Path(source["restart"])
    if sha256(restart) != source["restart_sha256"]:
        raise ValueError("fixed Step1 endpoint restart SHA256 mismatch")
    if sha256(PRMTOP) != EXPECTED_PRMTOP_SHA256:
        raise ValueError("frozen system.prmtop SHA256 mismatch")
    manifest = _load_json(pathlib.Path(source["manifest"]))
    result = _load_json(pathlib.Path(source["result"]))
    if result.get("status") != "PASS_TECHNICAL_A1_ACYL_ENDPOINT_STABILITY":
        raise ValueError("Step1 continuation lacks technical PASS")
    if result.get("seed") != source["seed"]:
        raise ValueError("Step1 result seed mismatch")
    contract = manifest.get("qm_contract", {})
    expected_step1 = {
        "qm_atom_count": 146,
        "qmcharge": 0,
        "electron_count_including_link_h": 510,
        "link_atom_count": 6,
        "step1_qm_water_count": 0,
    }
    for key, expected in expected_step1.items():
        if contract.get(key) != expected:
            raise ValueError(f"Step1 QM authority mismatch for {key}")
    stages = manifest.get("stages", [])
    if (
        not stages
        or stages[-1].get("stage") != "release_md"
        or stages[-1].get("technical_pass") is not True
        or stages[-1].get("restart_sha256") != source["restart_sha256"]
    ):
        raise ValueError("Step1 source manifest does not bind the fixed release restart")
    audit_path = code_root / SENSITIVITY_AUDIT_RELATIVE
    if (
        not audit_path.is_file()
        or git_blob_sha1(audit_path) != EXPECTED_SENSITIVITY_AUDIT_BLOB_SHA1
    ):
        raise ValueError("versioned acyl connectivity sensitivity authority changed")
    sensitivity = _load_json(audit_path)
    if (
        sensitivity.get("result", {}).get("status")
        != "PASS_ACYL_CONNECTIVITY_SENSITIVITY_V1_REPRODUCED"
        or sensitivity.get("result", {}).get("provisional_scientific_gate")
        != "PROVISIONAL_ACYL_ENDPOINT_SUPPORTED_FOR_EXPLORATORY_ENDPOINT_WORK"
    ):
        raise ValueError("acyl sensitivity authority is not eligible for exploratory A2")
    matching = [
        item
        for item in sensitivity.get("per_seed", [])
        if item.get("seed") == source["seed"]
    ]
    if len(matching) != 1 or matching[0].get("sensitivity_product_like_frames") != "20/20":
        raise ValueError("seed lacks 20/20 acyl sensitivity support")
    return manifest, sensitivity


def qmmm_block(qmmask: str) -> str:
    return f"""&qmmm
  qmmask='{qmmask}',
  qmcharge=0,
  spin=1,
  qm_theory='DFTB3',
  dftb_telec=200.0,
  qmshake=0,
/
"""


def md_input(seed: str, leg: int, velocity_seed: int, qmmask: str) -> str:
    return f"""NylC A1 Step2 A2 QM-water preflight {seed} leg{leg}
&cntrl
  imin=0, irest=0, ntx=1, nstlim={MD_STEPS}, dt={MD_DT_PS},
  ntb=1, cut=10.0, ntt=3, gamma_ln=2.0,
  temp0=300.0, tempi=300.0, ig={velocity_seed},
  ntc=1, ntf=1, ntwx={MD_NTWX}, ntpr=20,
  ntxo=1, ifqnt=1, ntr=0, nmropt=0,
/
{qmmm_block(qmmask)}
"""


def prepare(
    seed_index: int,
    output: pathlib.Path,
    scratch: pathlib.Path,
    code_root: pathlib.Path,
    github_commit: str,
) -> None:
    source = source_for_index(seed_index)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    output.mkdir(parents=True)
    if not re.fullmatch(r"[0-9a-f]{40}", github_commit):
        write_not_evaluated(
            output,
            source["seed"],
            "NOT_EVALUATED_CODE_AUTHORITY",
            "runtime code snapshot is not bound to a full 40-hex Git commit",
            {
                "observed": github_commit,
                "scientific_authority_commit": AUTHORITY_COMMIT,
            },
        )
        return
    try:
        source_manifest, sensitivity = validate_authority(source, code_root)
        import parmed as pmd

        structure = pmd.load_file(str(PRMTOP), xyz=str(source["restart"]))
        if (
            len(structure.atoms) != ACYL.EXPECTED_SYSTEM_ATOMS
            or structure.box is None
        ):
            raise ValueError("full-system restart atom count or box changed")
        ACYL._validate_structure(structure)
        funnel, candidates = select_water(structure)
    except Exception as error:
        write_not_evaluated(
            output,
            source["seed"],
            "NOT_EVALUATED_A2_PREPARATION_AUTHORITY",
            f"{type(error).__name__}: {error}",
        )
        return
    write_json(
        output / "WATER_SELECTION.json",
        {
            "schema_version": 1,
            "seed": source["seed"],
            "candidate_universe": "all named full-system water residues in fixed final restart",
            "filter": WATER_FILTER,
            "rank_target": RANK_TARGET,
            "funnel_denominators": funnel,
            "candidate_count": len(candidates),
            "candidates": candidates,
            "selection_policy": "lowest deterministic score; ties by residue and atom index",
        },
    )
    with (output / "WATER_CANDIDATES.tsv").open("w", encoding="utf-8") as handle:
        handle.write(
            "rank\tresidue_index1\tresidue_name\toxygen_index1\thydrogen_indices1"
            "\tdonor_h_index1\tc12_ow_A\to2_c12_ow_deg\th_nalpha_A"
            "\tow_h_nalpha_deg\tscore\n"
        )
        for rank, item in enumerate(candidates, 1):
            handle.write(
                f"{rank}\t{item['residue_index1']}\t{item['residue_name']}\t"
                f"{item['oxygen_index1']}\t"
                f"{','.join(str(value) for value in item['hydrogen_indices1'])}\t"
                f"{item['donor_h_index1']}\t{item['c12_ow_A']:.8f}\t"
                f"{item['o2_c12_ow_deg']:.8f}\t{item['h_nalpha_A']:.8f}\t"
                f"{item['ow_h_nalpha_deg']:.8f}\t{item['score']:.10f}\n"
            )
    if not candidates:
        write_not_evaluated(
            output,
            source["seed"],
            "NOT_EVALUATED_STEP2_NO_ATTACK_WATER",
            "no complete H2O passed the predeclared geometry filter",
            {"funnel_denominators": funnel},
        )
        return
    chosen = candidates[0]
    try:
        base_mask = str(source_manifest["qm_contract"]["qmmask"])
        base_indices = _parse_qmmask(base_mask)
        water_indices = [
            int(chosen["oxygen_index1"]),
            *[int(value) for value in chosen["hydrogen_indices1"]],
        ]
        if set(base_indices).intersection(water_indices):
            raise ValueError("selected water already overlaps Step1 QM region")
        qm_indices = base_indices + water_indices
        if len(qm_indices) != 149 or len(set(qm_indices)) != 149:
            raise ValueError("Step2 QM mask is not 149 unique atoms")
        boundaries = _boundary_bonds(structure, qm_indices)
        derived_contract = {
            "qm_atom_count": len(qm_indices),
            "qmcharge": 0,
            "spin": 1,
            "link_atom_count": len(boundaries),
            "electron_count": 510 + 8,
            "electron_derivation": (
                "frozen Step1 510 plus neutral H2O DFTB3 valence electrons "
                "(O=6,H=1,H=1)"
            ),
        }
        if any(
            derived_contract[key] != expected
            for key, expected in EXPECTED_CONTRACT.items()
        ):
            raise ValueError(
                f"independently derived Step2 contract {derived_contract} "
                f"!= {EXPECTED_CONTRACT}"
            )
        qmmask = ",".join(f"@{index}" for index in qm_indices)
        box = _box_lengths(structure)
        scratch.mkdir(parents=True, exist_ok=False)
        for leg, velocity_seed in enumerate(source["velocity_seeds"]):
            stage = scratch / f"a2_leg{leg}"
            stage.mkdir()
            (stage / "stage.in").write_text(
                md_input(source["seed"], leg, velocity_seed, qmmask),
                encoding="utf-8",
            )
            write_json(
                stage / "PREPARED.json",
                {
                    "seed": source["seed"],
                    "leg": leg,
                    "velocity_seed": velocity_seed,
                    "input_restart": str(source["restart"]),
                    "input_restart_sha256": source["restart_sha256"],
                    "expected_contract": EXPECTED_CONTRACT,
                },
            )
    except Exception as error:
        write_not_evaluated(
            output,
            source["seed"],
            "NOT_EVALUATED_STEP2_QM_CONTRACT_PREPARATION",
            f"{type(error).__name__}: {error}",
        )
        return
    manifest = {
        "schema_version": 1,
        "status": "READY_STEP2_A2_QMWATER_PREFLIGHT",
        "scientific_scope": SCIENTIFIC_SCOPE,
        "github_commit": github_commit,
        "authority_commit": AUTHORITY_COMMIT,
        "seed_index": seed_index,
        "seed": source["seed"],
        "source": {
            "restart": str(source["restart"]),
            "restart_sha256": source["restart_sha256"],
            "manifest": str(source["manifest"]),
            "result": str(source["result"]),
            "sensitivity_audit_sha256": EXPECTED_SENSITIVITY_AUDIT_BLOB_SHA1,
        },
        "prmtop": str(PRMTOP),
        "prmtop_sha256": EXPECTED_PRMTOP_SHA256,
        "selected_water": chosen,
        "qm_contract": {
            "expected": EXPECTED_CONTRACT,
            "derived": derived_contract,
            "boundary_bonds": boundaries,
            "qmmask": qmmask,
            "base_qm_atom_count": 146,
            "selected_complete_water_atom_indices1": water_indices,
        },
        "box_lengths_A": list(box),
        "protocol": {
            "legs": 2,
            "execution": "PARALLEL_TWO_INDEPENDENT_VELOCITY_LEGS_2X8_RANKS_NOT_PAIRED",
            "nstlim": MD_STEPS,
            "dt_ps": MD_DT_PS,
            "duration_ps": MD_STEPS * MD_DT_PS,
            "ntwx": MD_NTWX,
            "expected_frames_per_leg": EXPECTED_FRAMES,
            "tail_frames": TAIL_FRAMES,
            "reactive_restraints": False,
            "position_restraints": False,
            "trajectory_policy": "scratch_only",
        },
        "A2_gate": A2_GATE,
        "NEXT": NEXT,
        "product_endpoint_implemented": False,
    }
    write_json(output / "A2_MANIFEST.json", manifest)
    write_json(
        output / "READY.json",
        {
            "status": "READY_STEP2_A2_QMWATER_PREFLIGHT",
            "seed": source["seed"],
            "selected_water_oxygen_index1": chosen["oxygen_index1"],
            "NEXT": NEXT,
        },
    )


def _last_nstep(text: str) -> int | None:
    matches = re.findall(r"\bNSTEP\s*=\s*(\d+)", text, re.I)
    return int(matches[-1]) if matches else None


def engine_banner_region(text: str) -> str:
    """Return only engine QM initialization text, excluding the mdin echo."""
    marker = re.search(
        r"(?im)^\s*(?:\|\s*)?(?:QMMM:|QM/MM\b|QM-MM\b)",
        text,
    )
    if marker is None:
        return ""
    first_step = re.search(r"(?im)^\s*NSTEP\s*=", text[marker.start():])
    stop = (
        marker.start() + first_step.start()
        if first_step is not None
        else len(text)
    )
    return text[marker.start():stop]


def parse_banner(text: str) -> dict[str, list[int]]:
    region = engine_banner_region(text)
    patterns = {
        "qm_atom_count": (
            r"\bnquant\s*[=:]\s*(\d+)",
            r"number\s+of\s+(?:QM|quantum)\s+atoms\s*[=:]\s*(\d+)",
            r"(?:contains|there\s+are)\s+(\d+)\s+(?:QM|quantum)\s+atoms",
        ),
        "qmcharge": (
            r"(?:QMMM:|QM/MM|QM-MM)[^\n]{0,100}?\bcharge\b\s*[=:]\s*(-?\d+)",
            r"\bQM\s+charge\s*[=:]\s*(-?\d+)",
        ),
        "link_atom_count": (
            r"(?:number\s+of\s+(?:QM/MM\s+)?link\s+atoms|\bnlink)\s*[=:]\s*(\d+)",
            r"there\s+(?:are|is)\s+(\d+)\s+(?:QM/MM\s+)?link\s+atoms",
            r"link\s+atom\s+count\s*[=:]\s*(\d+)",
        ),
        "electron_count": (
            r"(?:total\s+)?number\s+of\s+(?:QM\s+|quantum\s+)?electrons\s*[=:]\s*(\d+)",
            r"\b(?:nelectrons|nelec)\s*[=:]\s*(\d+)",
            r"(?:specified\s+with|there\s+are)\s+(\d+)\s+(?:QM\s+|quantum\s+)?electrons",
        ),
    }
    observed: dict[str, list[int]] = {}
    for key, alternatives in patterns.items():
        values: list[int] = []
        for pattern in alternatives:
            values.extend(int(value) for value in re.findall(pattern, region, re.I))
        observed[key] = sorted(set(values))
    return observed


def banner_pass(observed: Mapping[str, Sequence[int]]) -> bool:
    return all(
        list(observed.get(key, [])) == [expected]
        for key, expected in EXPECTED_CONTRACT.items()
        if key != "spin"
    )


def _read_md_frames(
    trajectory: pathlib.Path, atom_count: int
) -> tuple[list[Any], list[str]]:
    import parmed as pmd

    captured: list[str] = []
    reader = None
    frames: list[Any] = []
    with warnings.catch_warnings(record=True) as observed:
        warnings.simplefilter("always")
        try:
            reader = pmd.amber.AmberMdcrd(
                str(trajectory), atom_count, hasbox=True, mode="r"
            )
            frames = list(reader.coordinates)
        finally:
            if reader is not None and hasattr(reader, "close"):
                reader.close()
        captured.extend(str(item.message) for item in observed)
    return frames, captured


def _coordinates(frame: Any) -> dict[int, Any]:
    return {index + 1: frame[index] for index in range(len(frame))}


def _acyl_sensitivity_like(g: Mapping[str, Any]) -> bool:
    gate = ACYL.PRODUCT_GATE
    return bool(
        g["attack_A"] <= gate["attack_A"][1]
        and gate["hg1_n3_A"][0] <= g["hg1_n3_A"] <= gate["hg1_n3_A"][1]
        and g["nalpha_hg1_A"] >= gate["nalpha_hg1_A_min"]
        and g["qPT_A"] >= gate["qPT_A_min"]
        and g["c12_n3_A"] >= gate["c12_n3_A_min"]
        and gate["c12_o2_A"][0] <= g["c12_o2_A"] <= gate["c12_o2_A"][1]
        and g["product_out_of_plane_A"] <= gate["product_out_of_plane_A_max"]
        and g["product_angle_sum_deg"] >= gate["product_angle_sum_deg_min"]
        and g["hg1_nearest_qm_heavy_atom"] == REACTIVE["n3"]
    )


def frame_geometry(
    frame: Any, manifest: Mapping[str, Any], topology: Any
) -> dict[str, Any]:
    coords = _coordinates(frame)
    box = tuple(float(value) for value in manifest["box_lengths_A"])
    water = manifest["selected_water"]
    ow = int(water["oxygen_index1"])
    hydrogens = [int(value) for value in water["hydrogen_indices1"]]
    heavy = [
        index
        for index in _parse_qmmask(manifest["qm_contract"]["qmmask"], 149)
        if _atomic_number(topology.atoms[index - 1]) > 1
    ]
    acyl = ACYL.geometry_from_coordinates(coords, heavy)
    r_add = _distance(coords[REACTIVE["c12"]], coords[ow], box)
    attack_angle = _angle(
        coords[REACTIVE["o2"]], coords[REACTIVE["c12"]], coords[ow], box
    )
    oh = [_distance(coords[ow], coords[index], box) for index in hydrogens]
    hbonds = []
    for index in hydrogens:
        hbonds.append(
            (
                _distance(coords[index], coords[REACTIVE["nalpha"]], box),
                _angle(
                    coords[ow],
                    coords[index],
                    coords[REACTIVE["nalpha"]],
                    box,
                ),
                index,
            )
        )
    hbonds.sort(key=lambda item: (item[0], -item[1], item[2]))
    h_nalpha, hbond_angle, donor_h = hbonds[0]
    result = dict(acyl)
    result.update(
        {
            "water_c12_ow_A": r_add,
            "water_attack_angle_deg": attack_angle,
            "water_oh_A": oh,
            "water_h_nalpha_A": h_nalpha,
            "water_hbond_angle_deg": hbond_angle,
            "water_donor_h_index1": donor_h,
        }
    )
    result["acyl_sensitivity_like"] = _acyl_sensitivity_like(result)
    result["neutral_water_like"] = all(
        distance <= A2_GATE["ow_h_A_max"] for distance in oh
    )
    result["no_water_addition"] = (
        r_add >= A2_GATE["c12_ow_no_addition_A_min"]
    )
    result["attack_resident"] = bool(
        r_add <= A2_GATE["resident_c12_ow_A_max"]
        and A2_GATE["resident_angle_deg"][0]
        <= attack_angle
        <= A2_GATE["resident_angle_deg"][1]
        and h_nalpha <= A2_GATE["resident_h_nalpha_A_max"]
        and hbond_angle >= A2_GATE["resident_hbond_angle_deg_min"]
    )
    return result


def audit_leg(output: pathlib.Path, scratch: pathlib.Path, leg: int) -> None:
    manifest = _load_json(output / "A2_MANIFEST.json")
    stage = scratch / f"a2_leg{leg}"
    text = (
        (stage / "stage.out").read_text(encoding="utf-8", errors="replace")
        if (stage / "stage.out").is_file()
        else ""
    )
    engine_rc = (
        int((stage / "engine.rc").read_text(encoding="utf-8").strip())
        if (stage / "engine.rc").is_file()
        else 999
    )
    banner = parse_banner(text)
    hard = {
        name: len(re.findall(pattern, text, re.I))
        for name, pattern in HARD_PATTERNS.items()
    }
    scc = len(re.findall(r"Convergence could not be achieved", text, re.I))
    vlimit = len(re.findall(r"vlimit\s+exceeded", text, re.I))
    overflow = len(re.findall(r"BOND\s*=\s*\*+", text, re.I))
    trajectory = stage / "a2.mdcrd"
    restart = stage / "stage.rst7"
    frames: list[Any] = []
    parse_warnings: list[str] = []
    if trajectory.is_file():
        try:
            frames, parse_warnings = _read_md_frames(
                trajectory, ACYL.EXPECTED_SYSTEM_ATOMS
            )
        except Exception as error:
            parse_warnings = [f"{type(error).__name__}: {error}"]
    technical = bool(
        engine_rc == 0
        and _last_nstep(text) == MD_STEPS
        and restart.is_file()
        and restart.stat().st_size > 0
        and len(frames) == EXPECTED_FRAMES
        and not parse_warnings
        and scc == 0
        and vlimit == 0
        and overflow == 0
        and sum(hard.values()) == 0
        and banner_pass(banner)
    )
    geometries: list[dict[str, Any]] = []
    geometry_error = None
    if technical:
        try:
            import parmed as pmd

            topology = pmd.load_file(str(PRMTOP))
            geometries = [
                frame_geometry(frame, manifest, topology) for frame in frames
            ]
            if not all(
                math.isfinite(float(value))
                for geometry in geometries
                for key, value in geometry.items()
                if isinstance(value, (int, float)) and not isinstance(value, bool)
            ):
                raise ValueError("non-finite frame geometry")
        except Exception as error:
            geometry_error = f"{type(error).__name__}: {error}"
            geometries = []
            technical = False
    tail = geometries[-TAIL_FRAMES:] if len(geometries) >= TAIL_FRAMES else []
    occupancies = {}
    for key in (
        "acyl_sensitivity_like",
        "neutral_water_like",
        "no_water_addition",
        "attack_resident",
    ):
        occupancies[key] = (
            sum(bool(frame[key]) for frame in tail) / len(tail) if tail else 0.0
        )
    scientific_pass = bool(
        technical
        and len(tail) == TAIL_FRAMES
        and occupancies["acyl_sensitivity_like"]
        >= A2_GATE["acyl_sensitivity_occupancy_min"]
        and occupancies["neutral_water_like"]
        >= A2_GATE["neutral_water_occupancy_min"]
        and occupancies["no_water_addition"]
        >= A2_GATE["no_addition_occupancy_min"]
        and occupancies["attack_resident"]
        >= A2_GATE["attack_residence_occupancy_min"]
    )
    result = {
        "schema_version": 1,
        "seed": manifest["seed"],
        "leg": leg,
        "engine_exit_code": engine_rc,
        "technical_pass": technical,
        "banner_contract_pass": banner_pass(banner),
        "banner_observed": banner,
        "banner_expected": EXPECTED_CONTRACT,
        "banner_region_found": bool(engine_banner_region(text)),
        "banner_region_sha256": hashlib.sha256(
            engine_banner_region(text).encode("utf-8", errors="replace")
        ).hexdigest(),
        "scientific_A2_leg_pass": scientific_pass,
        "classification": (
            "PASS_EXPLORATORY_A2_LEG"
            if scientific_pass
            else (
                "NOT_EVALUATED_A2_QM_CONTRACT"
                if not banner_pass(banner)
                else (
                    "NOT_EVALUATED_A2_TECHNICAL"
                    if not technical
                    else "A2_LEG_NOT_REVALIDATED"
                )
            )
        ),
        "frame_count": len(frames),
        "tail_denominator": len(tail),
        "tail_occupancies": occupancies,
        "final_geometry": geometries[-1] if geometries else {},
        "diagnostics": {
            "last_nstep": _last_nstep(text),
            "scc_warnings": scc,
            "vlimit_warnings": vlimit,
            "bond_overflow": overflow,
            "hard_error_hits": hard,
            "trajectory_parse_warnings": parse_warnings,
            "geometry_error": geometry_error,
        },
        "restart_path": str(restart),
        "restart_sha256": sha256(restart) if restart.is_file() else "",
        "trajectory_persisted": False,
        "NEXT": NEXT,
    }
    write_json(output / f"A2_LEG_{leg}.json", result)


def finalize(output: pathlib.Path) -> None:
    manifest_path = output / "A2_MANIFEST.json"
    if not manifest_path.is_file():
        if (
            (output / "RESULT.json").is_file()
            and (output / "NOT_EVALUATED.json").is_file()
        ):
            return
        raise ValueError("A2 manifest is absent and no terminal NOT_EVALUATED result exists")
    manifest = _load_json(manifest_path)
    paths = [output / f"A2_LEG_{leg}.json" for leg in (0, 1)]
    if not all(path.is_file() for path in paths):
        write_not_evaluated(
            output,
            manifest["seed"],
            "NOT_EVALUATED_A2_INCOMPLETE_DENOMINATOR",
            "both expected A2 legs are not auditable",
            {"present": [path.is_file() for path in paths]},
        )
        return
    legs = [_load_json(path) for path in paths]
    contract_ok = all(leg.get("banner_contract_pass") is True for leg in legs)
    technical = all(leg.get("technical_pass") is True for leg in legs)
    scientific = all(leg.get("scientific_A2_leg_pass") is True for leg in legs)
    if scientific:
        status = "PASS_EXPLORATORY_STEP2_A2_ACYL_BASIN_REVALIDATED"
        classification = status
    elif not contract_ok:
        status = "NOT_EVALUATED_STEP2_A2_QM_CONTRACT"
        classification = status
    elif not technical:
        status = "NOT_EVALUATED_STEP2_A2_TECHNICAL"
        classification = status
    else:
        status = "FAIL_STEP2_A2_NOT_REVALIDATED"
        classification = status
    result = {
        "schema_version": 1,
        "status": status,
        "scientific_scope": SCIENTIFIC_SCOPE,
        "seed": manifest["seed"],
        "classification": classification,
        "denominator": {
            "expected_legs": 2,
            "audited_legs": len(legs),
            "legs_are_independent_velocity_initializations": True,
            "legs_are_antithetic_pairs": False,
            "execution_was_serial": False,
        },
        "selected_water": manifest["selected_water"],
        "qm_contract": manifest["qm_contract"],
        "A2_gate": A2_GATE,
        "per_leg": legs,
        "interpretation": (
            "A2 endpoint preflight only. A PASS supports short-window persistence "
            "of the Step1 acyl connectivity with one fixed QM water; it is not a "
            "Step2 product, TS, path, PMF, barrier, or mechanism result."
        ),
        "NEXT": NEXT,
        "product_endpoint_implemented": False,
    }
    for name in ("PASS.json", "NOT_EVALUATED.json", "FAIL.json"):
        path = output / name
        if path.exists():
            path.unlink()
    write_json(output / "RESULT.json", result)
    if scientific:
        write_json(output / "PASS.json", result)
    elif not contract_ok or not technical:
        write_json(output / "NOT_EVALUATED.json", result)
    else:
        write_json(output / "FAIL.json", result)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode", choices=("prepare", "audit-leg", "finalize"), required=True
    )
    parser.add_argument("--seed-index", type=int)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--scratch", type=pathlib.Path)
    parser.add_argument("--code-root", type=pathlib.Path)
    parser.add_argument("--github-commit", default="unknown")
    parser.add_argument("--leg", type=int, choices=(0, 1))
    args = parser.parse_args()
    output = args.output.resolve()
    if args.mode == "prepare":
        if args.seed_index is None or args.scratch is None or args.code_root is None:
            parser.error("prepare requires --seed-index, --scratch and --code-root")
        prepare(
            args.seed_index,
            output,
            args.scratch.resolve(),
            args.code_root.resolve(),
            args.github_commit,
        )
    elif args.mode == "audit-leg":
        if args.scratch is None or args.leg is None:
            parser.error("audit-leg requires --scratch and --leg")
        audit_leg(output, args.scratch.resolve(), args.leg)
    else:
        finalize(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
