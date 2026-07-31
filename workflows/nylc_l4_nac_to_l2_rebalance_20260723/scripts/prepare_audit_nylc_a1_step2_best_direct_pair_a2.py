#!/usr/bin/env python3
"""Two-source adapter for the best seed26737 direct-water A2 events.

The established 149-QM A2 engine and audit remain unchanged. Array task 0
uses task4/frame375 and array task 1 uses task6/frame328. Each source is
SHA-frozen and uses the same OW13046/H13047/H13048 direct water.
"""
from __future__ import annotations

import importlib.util
import itertools
import json
import os
import pathlib
import re
from typing import Any, Mapping

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
BASE_PATH = HERE / "prepare_audit_nylc_a1_step2_qmwater_endpoint.py"
_SPEC = importlib.util.spec_from_file_location("_nylc_a2_base", BASE_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"cannot load frozen A2 base driver {BASE_PATH}")
BASE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(BASE)

TASK_ROOT = pathlib.Path(
    "/work/home/acshdt1dks/nylon_pa66_scnet_20260708/"
    "l4_nac_to_l2_rebalance_20260723"
)
SOURCE_TASKS = (4, 6)
SOURCE_RESTART_SHA256S = (
    "cbe9fbddb0ea8c3a6f12c9ae03d19864b6f4d6e2582901fe35ffa828d0f0c1ca",
    "5c4227dc8752f1a181e5b17d2890e6bce5a8a790f3ed5f91a181445e72006efe",
)
VELOCITY_SEEDS_BY_SLOT = (
    (26737421, 26737422),
    (26737631, 26737632),
)
SOURCE_SLOT = int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))
if SOURCE_SLOT not in (0, 1):
    raise ValueError("best-direct pair requires SLURM_ARRAY_TASK_ID 0 or 1")
SOURCE_TASK = SOURCE_TASKS[SOURCE_SLOT]
ATTEMPT_ROOT = (
    TASK_ROOT
    / "a1_activated_nac_20260726/qmmm/"
    f"a1_step2_water_network_nearmiss_continuation/attempt_62507122_{SOURCE_TASK}"
)
SOURCE_RESTART = ATTEMPT_ROOT / "direct_event_0.rst7"
SOURCE_RESTART_SHA256 = SOURCE_RESTART_SHA256S[SOURCE_SLOT]
SOURCE_RESULT = ATTEMPT_ROOT / "RESULT.json"
SOURCE_SAMPLING_MANIFEST = ATTEMPT_ROOT / "SOURCE_MANIFEST.json"
EXPECTED_SAMPLING_SOURCE_RESTART_SHA256 = (
    "42c46ebe61ad3016c86a91880ac1d6f94d7f7c7fdcda385bd89935277d378583"
)
SELECTED_WATER_ATOMS1 = (13046, 13047, 13048)
SELECTED_DONOR_H1 = 13048
ALLOWED_REACTED_TOPOLOGY_BONDS1 = {
    (8949, 8961),
    (10287, 10289),
}
VELOCITY_SEEDS = VELOCITY_SEEDS_BY_SLOT[SOURCE_SLOT]

_ORIGINAL_SOURCE_FOR_INDEX = BASE.source_for_index
_ORIGINAL_VALIDATE_AUTHORITY = BASE.validate_authority
_ORIGINAL_SELECT_WATER = BASE.select_water
_ORIGINAL_PREPARE = BASE.prepare
_ORIGINAL_AUDIT_LEG = BASE.audit_leg
_ORIGINAL_RESOLVE_LEG_BANNER = BASE.resolve_leg_banner_contract
EXPECTED_DFTB_DOUBLY_OCCUPIED = 194


def _load_json(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def source_for_index(index: int) -> dict[str, Any]:
    if index != 0:
        raise ValueError("best-direct pair A2 uses seed-index 0 within each array task")
    old = _ORIGINAL_SOURCE_FOR_INDEX(1)
    return {
        "seed": "seed26737",
        "attempt": f"attempt_62507122_{SOURCE_TASK}/direct_event_0",
        "restart_sha256": SOURCE_RESTART_SHA256,
        "velocity_seeds": VELOCITY_SEEDS,
        "root": ATTEMPT_ROOT,
        "restart": SOURCE_RESTART,
        "manifest": old["manifest"],
        "result": SOURCE_RESULT,
        "sampling_manifest": SOURCE_SAMPLING_MANIFEST,
    }


def _triclinic_cell(box: Any) -> tuple[np.ndarray, np.ndarray]:
    a, b, c, alpha, beta, gamma = (float(value) for value in box[:6])
    alpha, beta, gamma = np.deg2rad((alpha, beta, gamma))
    sin_gamma = float(np.sin(gamma))
    if min(a, b, c) <= 0.0 or abs(sin_gamma) <= 1.0e-12:
        raise ValueError("invalid triclinic periodic box")
    cell = np.array(
        [
            [a, b * np.cos(gamma), c * np.cos(beta)],
            [
                0.0,
                b * sin_gamma,
                c
                * (np.cos(alpha) - np.cos(beta) * np.cos(gamma))
                / sin_gamma,
            ],
            [0.0, 0.0, 0.0],
        ],
        dtype=float,
    )
    z2 = c * c - cell[0, 2] ** 2 - cell[1, 2] ** 2
    if not np.isfinite(z2) or z2 <= 1.0e-12:
        raise ValueError("degenerate triclinic periodic box")
    cell[2, 2] = np.sqrt(z2)
    return cell, np.linalg.inv(cell)


def audit_bonded_pbc_integrity(structure: Any) -> dict[str, Any]:
    cell, inverse = _triclinic_cell(structure.box)
    coordinates = np.array(
        [[atom.xx, atom.xy, atom.xz] for atom in structure.atoms], dtype=float
    )
    pairs0 = np.array(
        [(bond.atom1.idx, bond.atom2.idx) for bond in structure.bonds], dtype=int
    )
    displacement = coordinates[pairs0[:, 0]] - coordinates[pairs0[:, 1]]
    raw = np.linalg.norm(displacement, axis=1)
    fractional = displacement @ inverse.T
    centers = np.rint(fractional)
    best_squared = np.full(len(displacement), np.inf)
    for shift in itertools.product((-1, 0, 1), repeat=3):
        candidate = displacement - (centers + np.array(shift)) @ cell.T
        best_squared = np.minimum(
            best_squared, np.einsum("ij,ij->i", candidate, candidate)
        )
    minimum_image = np.sqrt(best_squared)
    long_indices = np.flatnonzero(raw > 3.0)
    long_pairs = {
        tuple(sorted((int(pairs0[i, 0]) + 1, int(pairs0[i, 1]) + 1)))
        for i in long_indices
    }
    pbc_split = [
        int(i)
        for i in long_indices
        if float(minimum_image[i]) <= 3.0
    ]
    if pbc_split:
        raise ValueError(
            f"PBC-split bonded fragments remain in selected direct-event restart: {len(pbc_split)}"
        )
    if long_pairs != ALLOWED_REACTED_TOPOLOGY_BONDS1:
        raise ValueError(
            f"unexpected long topology bonds {sorted(long_pairs)}; "
            f"expected only reacted QM bonds "
            f"{sorted(ALLOWED_REACTED_TOPOLOGY_BONDS1)}"
        )
    other = [
        float(raw[i])
        for i in range(len(raw))
        if tuple(
            sorted((int(pairs0[i, 0]) + 1, int(pairs0[i, 1]) + 1))
        )
        not in ALLOWED_REACTED_TOPOLOGY_BONDS1
    ]
    return {
        "bond_count": int(len(raw)),
        "raw_bonds_gt3_A": int(len(long_indices)),
        "pbc_split_bonds": 0,
        "allowed_reacted_topology_bonds1": [
            list(pair) for pair in sorted(ALLOWED_REACTED_TOPOLOGY_BONDS1)
        ],
        "max_nonreacted_raw_bond_A": max(other),
    }


def validate_best_direct_pair_authority(
    source: Mapping[str, Any], code_root: pathlib.Path
) -> tuple[dict[str, Any], dict[str, Any]]:
    old_source = _ORIGINAL_SOURCE_FOR_INDEX(1)
    old_manifest, sensitivity = _ORIGINAL_VALIDATE_AUTHORITY(
        old_source, code_root
    )
    if BASE.sha256(SOURCE_RESTART) != SOURCE_RESTART_SHA256:
        raise ValueError("selected direct-event restart SHA256 mismatch")
    result = _load_json(SOURCE_RESULT)
    sampling = _load_json(SOURCE_SAMPLING_MANIFEST)
    if (
        result.get("status") != "PASS_TECHNICAL_STEP2_WATER_NETWORK_SAMPLING"
        or result.get("technical_complete") is not True
        or result.get("seed") != "seed26737"
        or result.get("task_index") != SOURCE_TASK
        or result.get("frame_count") != 500
        or result.get("engine", {}).get("banner_contract_pass") is not True
    ):
        raise ValueError("selected water-network technical authority changed")
    if (
        sampling.get("source", {}).get("restart_sha256")
        != EXPECTED_SAMPLING_SOURCE_RESTART_SHA256
    ):
        raise ValueError("selected sampling lineage no longer binds seed26737 endpoint")
    matches = [
        event
        for event in result.get("direct_events", [])
        if event.get("selected_restart_sha256") == SOURCE_RESTART_SHA256
        and event.get("route") == "direct"
        and event.get("consecutive_frames", 0) >= 5
    ]
    if len(matches) != 1:
        raise ValueError("selected fixed direct event is absent or ambiguous")
    geometry = matches[0].get("representative_geometry", {})
    if (
        geometry.get("attack_water_oxygen_index1") != SELECTED_WATER_ATOMS1[0]
        or tuple(geometry.get("attack_water_hydrogen_indices1", ()))
        != SELECTED_WATER_ATOMS1[1:]
        or geometry.get("attack_donor_h_index1") != SELECTED_DONOR_H1
    ):
        raise ValueError("selected direct-event water identity changed")
    import parmed as pmd

    structure = pmd.load_file(str(BASE.PRMTOP), xyz=str(SOURCE_RESTART))
    audit = audit_bonded_pbc_integrity(structure)
    source["pbc_integrity_audit"] = audit
    return old_manifest, sensitivity


def fixed_select_water(structure: Any):
    funnel, candidates = _ORIGINAL_SELECT_WATER(structure)
    matches = [
        candidate
        for candidate in candidates
        if (
            int(candidate["oxygen_index1"]),
            *[int(value) for value in candidate["hydrogen_indices1"]],
        )
        == SELECTED_WATER_ATOMS1
        and int(candidate["donor_h_index1"]) == SELECTED_DONOR_H1
    ]
    if len(matches) != 1:
        raise ValueError(
            "fixed OW13046/H13047/H13048 direct water does not pass "
            "the unchanged 149-QM A2 filter"
        )
    return funnel, matches


def normalize_amber_atom_mask(qmmask: str) -> str:
    tokens = [token.strip() for token in qmmask.split(",") if token.strip()]
    if not tokens:
        raise ValueError("empty Amber atom mask")
    values = [token[1:] if token.startswith("@") else token for token in tokens]
    if any(not value.isdigit() for value in values):
        raise ValueError(f"non-explicit Amber atom mask: {qmmask!r}")
    return "@" + ",".join(values)


def prepare_with_amber18_mask(*args: Any, **kwargs: Any) -> None:
    _ORIGINAL_PREPARE(*args, **kwargs)
    output = pathlib.Path(kwargs.get("output", args[1] if len(args) > 1 else ""))
    manifest_path = output / "A2_MANIFEST.json"
    if not manifest_path.is_file():
        return
    manifest = _load_json(manifest_path)
    old_mask = str(manifest["qm_contract"]["qmmask"])
    new_mask = normalize_amber_atom_mask(old_mask)
    if BASE._parse_qmmask(new_mask, 149) != BASE._parse_qmmask(old_mask, 149):
        raise ValueError("Amber mask normalization changed atom identities")
    manifest["qm_contract"]["qmmask"] = new_mask
    BASE.write_json(manifest_path, manifest)
    scratch = pathlib.Path(kwargs.get("scratch", args[2] if len(args) > 2 else ""))
    needle = f"qmmask='{old_mask}'"
    replacement = f"qmmask='{new_mask}'"
    for leg in (0, 1):
        stage_input = scratch / f"a2_leg{leg}" / "stage.in"
        text = stage_input.read_text(encoding="utf-8")
        if text.count(needle) != 1:
            raise ValueError(f"unexpected qmmask occurrence count in {stage_input}")
        stage_input.write_text(text.replace(needle, replacement), encoding="utf-8")


def parse_direct_engine_banner(text: str) -> dict[str, list[int]]:
    """Parse Amber18 QMMM options plus DFTB valence-level authority."""
    start = re.search(r"(?im)^\s*QMMM\s+options:\s*$", text)
    if start is None:
        return {
            "qm_atom_count": [],
            "qmcharge": [],
            "link_atom_count": [],
            "electron_count": [],
            "doubly_occupied_levels": [],
        }
    first_step = re.search(r"(?im)^\s*NSTEP\s*=", text[start.start():])
    stop = start.start() + first_step.start() if first_step else len(text)
    region = text[start.start():stop]
    patterns = {
        "qm_atom_count": r"\bnquant\s*[=:]\s*(\d+)",
        "qmcharge": r"\bqmcharge\s*[=:]\s*(-?\d+)",
        "link_atom_count": r"\bnlink\s*[=:]\s*(\d+)",
        "doubly_occupied_levels": (
            r"NO\.\s+OF\s+DOUBLY\s+OCCUPIED\s+LEVELS\s*=\s*(\d+)"
        ),
    }
    observed = {
        key: sorted({int(value) for value in re.findall(pattern, region, re.I)})
        for key, pattern in patterns.items()
    }
    observed["electron_count"] = []
    return observed


def resolve_direct_leg_banner_contract(
    manifest: Mapping[str, Any],
    prepared: Mapping[str, Any],
    leg_banner: Mapping[str, Any],
) -> tuple[dict[str, list[int]], str]:
    expected = dict(BASE.EXPECTED_CONTRACT)
    raw_ok = (
        list(leg_banner.get("qm_atom_count", [])) == [149]
        and list(leg_banner.get("qmcharge", [])) == [0]
        and list(leg_banner.get("link_atom_count", [])) == [6]
        and list(leg_banner.get("doubly_occupied_levels", []))
        == [EXPECTED_DFTB_DOUBLY_OCCUPIED]
    )
    try:
        if not raw_ok:
            raise ValueError("direct leg Amber18/DFTB banner mismatch")
        if dict(manifest["qm_contract"]["expected"]) != expected:
            raise ValueError("manifest Step2 contract changed")
        if dict(prepared["expected_contract"]) != expected:
            raise ValueError("prepared Step2 contract changed")
        BASE._parse_qmmask(manifest["qm_contract"]["qmmask"], 149)
        effective = {
            "qm_atom_count": [149],
            "qmcharge": [0],
            "link_atom_count": [6],
            "electron_count": [518],
        }
        return effective, "LEG_ENGINE_DFTB_VALENCE_PLUS_FIXED_COMPOSITION"
    except (KeyError, TypeError, ValueError):
        return _ORIGINAL_RESOLVE_LEG_BANNER(
            manifest, prepared, leg_banner
        )


def _minimum_image_distance(
    coordinates: np.ndarray, left: int, right: int, box: Any
) -> float:
    cell, inverse = _triclinic_cell(box)
    displacement = coordinates[left - 1] - coordinates[right - 1]
    fractional = displacement @ inverse.T
    center = np.rint(fractional)
    best = np.inf
    for shift in itertools.product((-1, 0, 1), repeat=3):
        candidate = displacement - (center + np.array(shift)) @ cell.T
        best = min(best, float(np.linalg.norm(candidate)))
    return best


def audit_thr267_heavy_skeleton(restart: pathlib.Path) -> dict[str, Any]:
    import parmed as pmd

    structure = pmd.load_file(str(BASE.PRMTOP), xyz=str(restart))
    coordinates = np.array(
        [[atom.xx, atom.xy, atom.xz] for atom in structure.atoms], dtype=float
    )
    pairs = {
        "Nalpha-CA": (8949, 8952, 1.85),
        "CA-CB": (8952, 8954, 1.95),
        "CB-OG1": (8954, 8960, 1.85),
        "C267-N268": (8962, 8964, 1.85),
    }
    distances = {
        name: _minimum_image_distance(
            coordinates, left, right, structure.box
        )
        for name, (left, right, _limit) in pairs.items()
    }
    checks = {
        name: distances[name] <= limit
        for name, (_left, _right, limit) in pairs.items()
    }
    return {
        "pass": all(checks.values()),
        "distances_A": distances,
        "checks": checks,
    }


def audit_leg_with_persisted_thr_integrity(
    output: pathlib.Path, scratch: pathlib.Path, leg: int
) -> None:
    _ORIGINAL_AUDIT_LEG(output, scratch, leg)
    result_path = output / f"A2_LEG_{leg}.json"
    result = _load_json(result_path)
    restart = scratch / f"a2_leg{leg}" / "stage.rst7"
    integrity = (
        audit_thr267_heavy_skeleton(restart)
        if restart.is_file()
        else {"pass": False, "distances_A": {}, "checks": {}}
    )
    result["thr267_heavy_skeleton"] = integrity
    if not integrity["pass"]:
        result["technical_pass"] = False
        result["scientific_A2_leg_pass"] = False
        result["classification"] = "NOT_EVALUATED_A2_CHEMICAL_INTEGRITY"
    BASE.write_json(result_path, result)


BASE.source_for_index = source_for_index
BASE.validate_authority = validate_best_direct_pair_authority
BASE.select_water = fixed_select_water
BASE.prepare = prepare_with_amber18_mask
BASE.parse_banner = parse_direct_engine_banner
BASE.resolve_leg_banner_contract = resolve_direct_leg_banner_contract
BASE.audit_leg = audit_leg_with_persisted_thr_integrity
BASE.SCIENTIFIC_SCOPE = (
    "EXPLORATORY_BEST_DIRECT_PAIR_149QM_A2_PREFLIGHT_ONLY_"
    "NOT_PRODUCT_TS_PATH_PMF_BARRIER_OR_MECHANISM"
)


def main() -> int:
    return BASE.main()


if __name__ == "__main__":
    raise SystemExit(main())
