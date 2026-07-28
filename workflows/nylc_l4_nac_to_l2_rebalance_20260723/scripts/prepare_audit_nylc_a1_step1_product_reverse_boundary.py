#!/usr/bin/env python3
"""NylC A1 Step1 product-side reverse-boundary scout.

Restrained windows are geometry probes only. A future shooting candidate can be
written only from a later fully reactive-restraint-free release trajectory.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import importlib.util
import json
import math
import os
import pathlib
import re
import shutil
from typing import Any, Mapping, Sequence


WORKFLOW_ROOT = pathlib.Path(__file__).resolve().parents[1]
AC_PATH = pathlib.Path(__file__).with_name(
    "prepare_audit_nylc_a1_step1_acyl_endpoint_stability.py"
)
TETRA_PATH = pathlib.Path(__file__).with_name(
    "prepare_audit_nylc_a1_step1_tetra_parallel_scout.py"
)


def _load_module(name: str, path: pathlib.Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import required module {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


AC = _load_module("_a1_acyl_endpoint_base", AC_PATH)
TETRA = _load_module("_a1_tetra_scout_base", TETRA_PATH)
BASE = AC.BASE

TASK_ROOT = pathlib.Path(
    "/work/home/acshdt1dks/nylon_pa66_scnet_20260708/"
    "l4_nac_to_l2_rebalance_20260723"
)
SOURCE_ROOT = (
    TASK_ROOT
    / "a1_activated_nac_20260726/qmmm/"
    "a1_step1_acyl_release_md_continuation"
)
AUTHORITY_AUDIT = (
    WORKFLOW_ROOT
    / "audit/nylc_a1_acyl_endpoint_connectivity_sensitivity_62216380_v1.json"
)
AUTHORITY_COMMIT = "19f327b422fd5081b119e4ff4c1883edefd85ba4"
SOURCE_JOB = "62216380"
PRMTOP = AC.PRMTOP

SOURCES = (
    {
        "seed_index": 0,
        "seed": "seed26723",
        "attempt": "attempt_62216380_0",
        "restart": SOURCE_ROOT / "attempt_62216380_0/release_md_endpoint.rst7",
        "manifest": SOURCE_ROOT / "attempt_62216380_0/ENDPOINT_MANIFEST.json",
        "result": SOURCE_ROOT / "attempt_62216380_0/RESULT.json",
        "restart_sha256": (
            "5d8f76d2c90e3e8c707b640c55a93938f53e18dc30d6f93d54adda467da25f41"
        ),
    },
    {
        "seed_index": 1,
        "seed": "seed26737",
        "attempt": "attempt_62216380_1",
        "restart": SOURCE_ROOT / "attempt_62216380_1/release_md_endpoint.rst7",
        "manifest": SOURCE_ROOT / "attempt_62216380_1/ENDPOINT_MANIFEST.json",
        "result": SOURCE_ROOT / "attempt_62216380_1/RESULT.json",
        "restart_sha256": (
            "4cc60ad4d7be099b3f76040f1ee8c49b91deecebb20ed98570bc192b3d511132"
        ),
    },
)

WINDOWS = (
    {
        "window_index": 0,
        "targets_A": {"attack": 1.55, "cn": 2.05, "nalpha_hg1": 1.45},
        "forces_kcal_mol_A2": {"attack": 10.0, "cn": 10.0, "nalpha_hg1": 5.0},
        "maxcyc": 300,
        "ncyc": 100,
    },
    {
        "window_index": 1,
        "targets_A": {"attack": 1.70, "cn": 1.90, "nalpha_hg1": 1.35},
        "forces_kcal_mol_A2": {"attack": 15.0, "cn": 15.0, "nalpha_hg1": 7.5},
        "maxcyc": 400,
        "ncyc": 100,
    },
    {
        "window_index": 2,
        "targets_A": {"attack": 1.85, "cn": 1.75, "nalpha_hg1": 1.25},
        "forces_kcal_mol_A2": {"attack": 20.0, "cn": 20.0, "nalpha_hg1": 10.0},
        "maxcyc": 500,
        "ncyc": 150,
    },
    {
        "window_index": 3,
        "targets_A": {"attack": 2.00, "cn": 1.60, "nalpha_hg1": 1.10},
        "forces_kcal_mol_A2": {"attack": 25.0, "cn": 25.0, "nalpha_hg1": 15.0},
        "maxcyc": 600,
        "ncyc": 150,
    },
)

RELEASE_STEPS = 500
RELEASE_DT_PS = 0.0005
RELEASE_NTWX = 10
RELEASE_FRAME_COUNT = RELEASE_STEPS // RELEASE_NTWX
RELEASE_EXCLUDED_INITIAL_PS = 0.05
MIN_BOUNDARY_RUN_FRAMES = 3
MAX_CANDIDATES_PER_SEED = 2

SCIENTIFIC_BOUNDARY = "NOT_EVALUATED_TS_COMMITTOR_PMF_BARRIER_MECHANISM"
SCOPE = (
    "step1_product_side_reverse_boundary_scout_not_ts_committor_pmf_barrier_or_mechanism"
)
TECHNICAL_PASS = "PASS_TECHNICAL_A1_PRODUCT_REVERSE_BOUNDARY_SCOUT"
TECHNICAL_FAIL = "NOT_EVALUATED_TECHNICAL_A1_PRODUCT_REVERSE_BOUNDARY_SCOUT"

REACTIVE = AC.REACTIVE_ATOMS
ALLOWED_HG1_NEAREST = {REACTIVE["nalpha"], REACTIVE["n3"]}


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: pathlib.Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def require_nonempty(path: pathlib.Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing required nonempty file {path}")


def source_from_index(seed_index: int) -> dict[str, Any]:
    index = int(seed_index)
    if index not in (0, 1):
        raise ValueError("seed index must be 0 or 1")
    return dict(SOURCES[index])


def window_spec(window_index: int) -> dict[str, Any]:
    index = int(window_index)
    if not 0 <= index < len(WINDOWS):
        raise ValueError(f"window index {index} is outside 0..{len(WINDOWS) - 1}")
    return dict(WINDOWS[index])


def _finite_geometry(geometry: Mapping[str, Any]) -> bool:
    if not geometry:
        return False
    return all(
        isinstance(value, (int, float)) and math.isfinite(float(value))
        for key, value in geometry.items()
        if key != "hg1_nearest_qm_heavy_atom"
    )


def product_side_sensitivity_gate(geometry: Mapping[str, Any]) -> bool:
    """The versioned v1 gate removes only the formal attack lower bound."""
    gate = AC.PRODUCT_GATE
    return bool(
        geometry
        and geometry["attack_A"] <= gate["attack_A"][1]
        and gate["hg1_n3_A"][0]
        <= geometry["hg1_n3_A"]
        <= gate["hg1_n3_A"][1]
        and geometry["nalpha_hg1_A"] >= gate["nalpha_hg1_A_min"]
        and geometry["qPT_A"] >= gate["qPT_A_min"]
        and geometry["c12_n3_A"] >= gate["c12_n3_A_min"]
        and gate["c12_o2_A"][0]
        <= geometry["c12_o2_A"]
        <= gate["c12_o2_A"][1]
        and geometry["product_out_of_plane_A"]
        <= gate["product_out_of_plane_A_max"]
        and geometry["product_angle_sum_deg"]
        >= gate["product_angle_sum_deg_min"]
        and geometry["hg1_nearest_qm_heavy_atom"] == REACTIVE["n3"]
    )


def validate_authority(seed_index: int) -> tuple[dict[str, Any], str, dict[str, Any]]:
    source = source_from_index(seed_index)
    for path in (
        AUTHORITY_AUDIT,
        source["restart"],
        source["manifest"],
        source["result"],
        PRMTOP,
    ):
        require_nonempty(path)

    audit = json.loads(AUTHORITY_AUDIT.read_text(encoding="utf-8"))
    if audit.get("result", {}).get("status") != (
        "PASS_ACYL_CONNECTIVITY_SENSITIVITY_V1_REPRODUCED"
    ):
        raise ValueError("connectivity sensitivity authority is not PASS")
    if str(audit.get("source", {}).get("job_id")) != SOURCE_JOB:
        raise ValueError("connectivity sensitivity source job changed")
    if "ALLOW_PRODUCT_SIDE_REVERSE_BOUNDARY_SCOUT" not in str(
        audit.get("result", {}).get("next_action", "")
    ):
        raise ValueError("authority does not allow the product-side reverse scout")
    if sha256(source["restart"]) != source["restart_sha256"]:
        raise ValueError("source release_md endpoint SHA256 changed")
    if sha256(PRMTOP) != BASE.EXPECTED_PRMTOP_SHA256:
        raise ValueError("frozen Step1 prmtop SHA256 changed")

    source_manifest = json.loads(source["manifest"].read_text(encoding="utf-8"))
    source_result = json.loads(source["result"].read_text(encoding="utf-8"))
    if source_manifest.get("seed") != source["seed"]:
        raise ValueError("source manifest seed changed")
    if int(source_manifest.get("seed_index", -1)) != seed_index:
        raise ValueError("source manifest seed index changed")
    if source_manifest.get("status") != AC.TECHNICAL_PASS:
        raise ValueError("source manifest is not technical PASS")
    if source_result.get("seed") != source["seed"]:
        raise ValueError("source result seed changed")
    if int(source_result.get("seed_index", -1)) != seed_index:
        raise ValueError("source result seed index changed")
    if source_result.get("status") != AC.TECHNICAL_PASS:
        raise ValueError("source result status is not technical PASS")
    if source_result.get("technical_complete") is not True:
        raise ValueError("source result is not technically complete")
    stages = source_manifest.get("stages", [])
    if not stages or stages[-1].get("stage") != "release_md":
        raise ValueError("source manifest lacks terminal release_md stage")
    if stages[-1].get("technical_pass") is not True:
        raise ValueError("source release_md stage is not technical PASS")
    if stages[-1].get("restart_sha256") != source["restart_sha256"]:
        raise ValueError("source manifest terminal restart SHA256 changed")
    expected_chain = [
        {
            "stage": stage.get("stage"),
            "input_restart_sha256": stage.get("input_restart_sha256"),
            "output_restart_sha256": stage.get("restart_sha256"),
        }
        for stage in stages
    ]
    if source_result.get("stage_sha_chain") != expected_chain:
        raise ValueError("source result SHA chain differs from source manifest")
    if expected_chain[-1]["output_restart_sha256"] != source["restart_sha256"]:
        raise ValueError("source result terminal restart SHA256 changed")

    qm_contract = source_manifest.get("qm_contract", {})
    expected_qm = {
        "qm_atom_count": BASE.EXPECTED_QM_ATOMS,
        "qmcharge": BASE.QMCHARGE,
        "electron_count_including_link_h": BASE.EXPECTED_ELECTRONS,
        "link_atom_count": BASE.EXPECTED_LINK_ATOMS,
        "step1_qm_water_count": 0,
    }
    for key, expected in expected_qm.items():
        if qm_contract.get(key) != expected:
            raise ValueError(f"source QM contract changed for {key}")

    authority = AC.validate_source_authority(AC.SOURCES[seed_index])
    if not isinstance(authority, tuple) or len(authority) != 2:
        raise ValueError("base authority did not return qmmask and source coordinate")
    qmmask = authority[0]
    if qm_contract.get("qmmask") != qmmask:
        raise ValueError("source qmmask differs from frozen base authority")

    import parmed as pmd

    structure = pmd.load_file(str(PRMTOP), xyz=str(source["restart"]))
    AC._validate_structure(structure)
    geometry = AC.geometry_from_coordinates(
        AC._coordinate_map(structure), qm_contract["qm_heavy_atom_indices"]
    )
    if not _finite_geometry(geometry):
        raise ValueError("source endpoint geometry is not finite")
    if not product_side_sensitivity_gate(geometry):
        raise ValueError(
            "source final restart does not pass the versioned product-side sensitivity gate"
        )
    return source_manifest, qmmask, geometry


def initialize(
    seed_index: int,
    root: pathlib.Path,
    start_rst7: pathlib.Path,
    github_commit: str,
) -> None:
    if root.exists():
        raise FileExistsError(f"attempt root already exists: {root}")
    if start_rst7.exists() or not start_rst7.parent.is_dir():
        raise ValueError("scratch start restart target must be fresh")
    source_manifest, qmmask, geometry = validate_authority(seed_index)
    source = source_from_index(seed_index)
    shutil.copy2(source["restart"], start_rst7)
    if sha256(start_rst7) != source["restart_sha256"]:
        raise ValueError("scratch copy of source restart failed SHA256 verification")
    root.mkdir(parents=True, exist_ok=False)
    payload = {
        "schema_version": 1,
        "status": "READY_A1_PRODUCT_REVERSE_BOUNDARY_SCOUT",
        "scientific_status": SCIENTIFIC_BOUNDARY,
        "scope": SCOPE,
        "github_commit": github_commit,
        "authority": {
            "audit": str(AUTHORITY_AUDIT),
            "authority_commit": AUTHORITY_COMMIT,
            "authority_status": (
                "PASS_ACYL_CONNECTIVITY_SENSITIVITY_V1_REPRODUCED"
            ),
            "formal_endpoint_result_preserved": (
                "FAIL_NO_RELEASE_STABLE_A1_ACYL_PRODUCT_ENDPOINT"
            ),
        },
        "seed_index": seed_index,
        "seed": source["seed"],
        "source": {
            "attempt": source["attempt"],
            "restart": str(source["restart"]),
            "restart_sha256": source["restart_sha256"],
            "manifest": str(source["manifest"]),
            "result": str(source["result"]),
            "source_job": SOURCE_JOB,
        },
        "prmtop": str(PRMTOP),
        "prmtop_sha256": sha256(PRMTOP),
        "qm_contract": dict(source_manifest["qm_contract"]),
        "source_geometry": geometry,
        "start_restart": str(start_rst7),
        "start_restart_sha256": sha256(start_rst7),
        "window_specs": list(WINDOWS),
        "windows": [],
        "releases": [],
        "candidates": [],
        "candidate_limit": MAX_CANDIDATES_PER_SEED,
        "restraint_policy": {
            "restrained_coordinates": [
                "C12_OG1_breaking",
                "C12_N3_forming",
                "Nalpha_HG1_return",
            ],
            "unrestrained_response_coordinates": [
                "HG1_N3",
                "C12_O2",
                "C12_out_of_plane",
                "attack_angle",
            ],
            "restrained_endpoint_shooting_forbidden": True,
        },
    }
    write_json(root / "SOURCE_MANIFEST.json", payload)


def reverse_restraints(spec: Mapping[str, Any]) -> str:
    targets = spec["targets_A"]
    forces = spec["forces_kcal_mol_A2"]
    return "".join(
        (
            AC.distance_restraint(
                REACTIVE["og1"], REACTIVE["c12"], targets["attack"], forces["attack"]
            ),
            AC.distance_restraint(
                REACTIVE["c12"], REACTIVE["n3"], targets["cn"], forces["cn"]
            ),
            AC.distance_restraint(
                REACTIVE["nalpha"],
                REACTIVE["hg1"],
                targets["nalpha_hg1"],
                forces["nalpha_hg1"],
            ),
        )
    )


def minimization_input(seed: str, spec: Mapping[str, Any], qmmask: str) -> str:
    return f"""NylC A1 product reverse boundary {seed} window {spec['window_index']}
&cntrl
  imin=1, ntmin=2, maxcyc={spec['maxcyc']}, ncyc={spec['ncyc']}, dx0=0.005,
  ntb=1, cut=10.0, ntpr=5, ntxo=1, ifqnt=1,
  ntr=1, nmropt=1, restraint_wt=1.0,
  restraintmask='{BASE.NON_QM_SOLUTE_HEAVY_MASK}', drms=0.10,
/
{BASE.qmmm_block(qmmask)}&wt type='END' /
DISANG=restraints.RST
DUMPAVE=restraint.dat
"""


def release_seed(seed: str, window_index: int) -> int:
    base = 126723 if seed == "seed26723" else 126737
    return base + 1000 * int(window_index)


def release_input(seed: str, window_index: int, qmmask: str) -> str:
    ig = release_seed(seed, window_index)
    return f"""NylC A1 product reverse boundary unrestrained release {seed} window {window_index}
&cntrl
  imin=0, irest=0, ntx=1, nstlim={RELEASE_STEPS}, dt={RELEASE_DT_PS},
  ntb=1, cut=10.0, ntt=3, gamma_ln=2.0,
  tempi=300.0, temp0=300.0, ig={ig},
  ntc=1, ntf=1, ntwx={RELEASE_NTWX}, ntpr=10,
  ntxo=1, ifqnt=1, ntr=0, nmropt=0,
/
{BASE.qmmm_block(qmmask)}
"""


def prepare_window(
    window_index: int,
    root: pathlib.Path,
    output: pathlib.Path,
    scratch: pathlib.Path,
    input_rst7: pathlib.Path,
) -> None:
    if output.exists() or scratch.exists():
        raise FileExistsError("window output or scratch directory already exists")
    require_nonempty(input_rst7)
    manifest = json.loads((root / "SOURCE_MANIFEST.json").read_text(encoding="utf-8"))
    spec = window_spec(window_index)
    if len(manifest["windows"]) != window_index:
        raise ValueError("window order is not strict serial 0..3")
    expected_input_sha = (
        manifest["start_restart_sha256"]
        if window_index == 0
        else manifest["windows"][-1]["output_restart_sha256"]
    )
    if sha256(input_rst7) != expected_input_sha:
        raise ValueError("window input restart breaks SHA256 inheritance")

    output.mkdir(parents=True, exist_ok=False)
    scratch.mkdir(parents=True, exist_ok=False)
    (scratch / "stage.in").write_text(
        minimization_input(manifest["seed"], spec, manifest["qm_contract"]["qmmask"]),
        encoding="utf-8",
    )
    (scratch / "restraints.RST").write_text(
        reverse_restraints(spec), encoding="utf-8"
    )
    prepared = {
        "schema_version": 1,
        "status": "READY_RESTRAINED_REVERSE_WINDOW",
        "seed": manifest["seed"],
        "window_spec": spec,
        "input_restart": str(input_rst7),
        "input_restart_sha256": expected_input_sha,
        "restrained": True,
        "shooting_forbidden": True,
        "unrestrained_response_coordinates": ["C12_O2", "C12_OOP", "HG1_N3"],
    }
    write_json(output / "WINDOW_MANIFEST.json", prepared)
    write_json(scratch / "PREPARED.json", prepared)


def chemical_guard(geometry: Mapping[str, Any]) -> dict[str, bool]:
    nearest = geometry.get("hg1_nearest_qm_heavy_atom")
    checks = {
        "finite_geometry": _finite_geometry(geometry),
        "c12_o2_1p15_to_1p55_A": bool(
            geometry and 1.15 <= geometry["c12_o2_A"] <= 1.55
        ),
        "reactant_plane_oop_le_0p60_A": bool(
            geometry
            and geometry["c12_reactant_plane_out_of_plane_A"] <= 0.60
        ),
        "attack_not_overcompressed": bool(
            geometry and geometry["attack_A"] >= 1.25
        ),
        "cn_not_overcompressed": bool(
            geometry and geometry["c12_n3_A"] >= 1.25
        ),
        "proton_not_detached": bool(
            geometry
            and min(geometry["nalpha_hg1_A"], geometry["hg1_n3_A"]) <= 1.60
        ),
        "proton_not_misrouted": nearest in ALLOWED_HG1_NEAREST,
        "attack_angle_80_to_140_deg": bool(
            geometry and 80.0 <= geometry["attack_angle_deg"] <= 140.0
        ),
    }
    checks["all"] = all(checks.values())
    return checks


def opposite_response(
    previous: Mapping[str, Any], current: Mapping[str, Any]
) -> dict[str, bool]:
    checks = {
        "attack_not_reverse_by_gt_0p15_A": bool(
            current["attack_A"] >= previous["attack_A"] - 0.15
        ),
        "cn_not_reverse_by_gt_0p15_A": bool(
            current["c12_n3_A"] <= previous["c12_n3_A"] + 0.15
        ),
        "nalpha_hg1_not_reverse_by_gt_0p15_A": bool(
            current["nalpha_hg1_A"] <= previous["nalpha_hg1_A"] + 0.15
        ),
    }
    checks["all"] = all(checks.values())
    return checks


def restrained_hint(geometry: Mapping[str, Any]) -> dict[str, bool]:
    checks = {
        "attack_1p45_to_2p20_A": bool(
            geometry and 1.45 <= geometry["attack_A"] <= 2.20
        ),
        "cn_1p45_to_2p20_A": bool(
            geometry and 1.45 <= geometry["c12_n3_A"] <= 2.20
        ),
        "bond_balance_le_0p60_A": bool(
            geometry
            and abs(geometry["attack_A"] - geometry["c12_n3_A"]) <= 0.60
        ),
        "abs_qpt_le_0p50_A": bool(
            geometry and abs(geometry["qPT_A"]) <= 0.50
        ),
        "attack_angle_90_to_130_deg": bool(
            geometry and 90.0 <= geometry["attack_angle_deg"] <= 130.0
        ),
    }
    checks["all"] = all(checks.values())
    return checks


def audit_window(
    root: pathlib.Path, output: pathlib.Path, scratch: pathlib.Path
) -> dict[str, Any]:
    manifest_path = root / "SOURCE_MANIFEST.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    prepared = json.loads((output / "WINDOW_MANIFEST.json").read_text(encoding="utf-8"))
    restart = scratch / "stage.rst7"
    technical, geometry, diagnostics = TETRA._technical(
        scratch / "stage.out", restart, manifest
    )
    previous = (
        manifest["source_geometry"]
        if not manifest["windows"]
        else manifest["windows"][-1]["geometry"]
    )
    guard = chemical_guard(geometry) if technical else {"all": False}
    direction = (
        opposite_response(previous, geometry) if technical else {"all": False}
    )
    hint = restrained_hint(geometry) if technical else {"all": False}
    guard_stop = bool(technical and (not guard["all"] or not direction["all"]))
    hint_pass = bool(technical and guard["all"] and direction["all"] and hint["all"])
    classification = (
        "NOT_EVALUATED_TECHNICAL_FAILURE"
        if not technical
        else "OUT_OF_SCOPE_CHEMISTRY_OR_OPPOSITE_RESPONSE"
        if guard_stop
        else "REVERSE_BOUNDARY_RESTRAINED_HINT"
        if hint_pass
        else "REVERSE_SCAN_NO_HINT"
    )
    result = {
        "schema_version": 1,
        "status": (
            "PASS_TECHNICAL_A1_PRODUCT_REVERSE_WINDOW"
            if technical
            else "NOT_EVALUATED_TECHNICAL_A1_PRODUCT_REVERSE_WINDOW"
        ),
        "scientific_status": SCIENTIFIC_BOUNDARY,
        "classification": classification,
        "technical_pass": technical,
        "shooting_forbidden": True,
        "guard_stop": guard_stop,
        "restrained_hint": hint_pass,
        "chemical_guard": guard,
        "direction_guard": direction,
        "hint_gate": hint,
        "geometry": geometry,
        "diagnostics": diagnostics,
        "window_spec": prepared["window_spec"],
        "input_restart_sha256": prepared["input_restart_sha256"],
        "output_restart": str(restart),
        "output_restart_sha256": (
            sha256(restart) if restart.is_file() and restart.stat().st_size else None
        ),
    }
    write_json(output / "RESULT.json", result)
    write_json(
        output / ("PASS.json" if technical else "NOT_EVALUATED.json"), result
    )
    manifest["windows"].append(result)
    write_json(manifest_path, manifest)
    if not technical:
        raise RuntimeError("restrained reverse window did not technically PASS")
    return result


def prepare_release(
    root: pathlib.Path,
    window_index: int,
    output: pathlib.Path,
    scratch: pathlib.Path,
    input_rst7: pathlib.Path,
) -> None:
    if output.exists() or scratch.exists():
        raise FileExistsError("release output or scratch directory already exists")
    require_nonempty(input_rst7)
    manifest = json.loads((root / "SOURCE_MANIFEST.json").read_text(encoding="utf-8"))
    matching = [
        item for item in manifest["windows"] if item["window_spec"]["window_index"] == window_index
    ]
    if len(matching) != 1:
        raise ValueError("release does not map to exactly one audited window")
    window = matching[0]
    if window.get("technical_pass") is not True or window.get("restrained_hint") is not True:
        raise ValueError("release requires a technical-PASS restrained hint")
    if window.get("guard_stop"):
        raise ValueError("guard-stopped window cannot enter release")
    if sha256(input_rst7) != window["output_restart_sha256"]:
        raise ValueError("release input does not match the restrained window restart")
    if len(manifest["candidates"]) >= MAX_CANDIDATES_PER_SEED:
        raise ValueError("candidate limit already reached")

    output.mkdir(parents=True, exist_ok=False)
    scratch.mkdir(parents=True, exist_ok=False)
    (scratch / "stage.in").write_text(
        release_input(
            manifest["seed"], window_index, manifest["qm_contract"]["qmmask"]
        ),
        encoding="utf-8",
    )
    prepared = {
        "schema_version": 1,
        "status": "READY_FULLY_UNRESTRAINED_REVERSE_RELEASE",
        "seed": manifest["seed"],
        "window_index": window_index,
        "input_restart": str(input_rst7),
        "input_restart_sha256": window["output_restart_sha256"],
        "ig": release_seed(manifest["seed"], window_index),
        "nstlim": RELEASE_STEPS,
        "dt_ps": RELEASE_DT_PS,
        "ntwx": RELEASE_NTWX,
        "ntr": 0,
        "nmropt": 0,
        "disang": None,
        "reactive_restraints": [],
        "environment_restraint": False,
        "trajectory_policy": "scratch_only",
        "initial_excluded_ps": RELEASE_EXCLUDED_INITIAL_PS,
    }
    write_json(output / "RELEASE_MANIFEST.json", prepared)
    write_json(scratch / "PREPARED.json", prepared)


def boundary_checks(geometry: Mapping[str, Any]) -> dict[str, bool]:
    nearest = geometry.get("hg1_nearest_qm_heavy_atom")
    tetra_response = bool(
        geometry
        and (
            geometry["c12_o2_A"] >= 1.28
            or geometry["c12_reactant_plane_out_of_plane_A"] >= 0.15
        )
    )
    checks = {
        "attack_1p55_to_2p15_A": bool(
            geometry and 1.55 <= geometry["attack_A"] <= 2.15
        ),
        "cn_1p50_to_2p10_A": bool(
            geometry and 1.50 <= geometry["c12_n3_A"] <= 2.10
        ),
        "bond_balance_le_0p45_A": bool(
            geometry
            and abs(geometry["attack_A"] - geometry["c12_n3_A"]) <= 0.45
        ),
        "nalpha_hg1_1p10_to_1p55_A": bool(
            geometry and 1.10 <= geometry["nalpha_hg1_A"] <= 1.55
        ),
        "hg1_n3_1p10_to_1p55_A": bool(
            geometry and 1.10 <= geometry["hg1_n3_A"] <= 1.55
        ),
        "abs_qpt_le_0p30_A": bool(
            geometry and abs(geometry["qPT_A"]) <= 0.30
        ),
        "attack_angle_90_to_130_deg": bool(
            geometry and 90.0 <= geometry["attack_angle_deg"] <= 130.0
        ),
        "c12_o2_1p20_to_1p50_A": bool(
            geometry and 1.20 <= geometry["c12_o2_A"] <= 1.50
        ),
        "reactant_plane_oop_le_0p60_A": bool(
            geometry
            and geometry["c12_reactant_plane_out_of_plane_A"] <= 0.60
        ),
        "tetrahedral_response": tetra_response,
        "proton_not_misrouted": nearest in ALLOWED_HG1_NEAREST,
    }
    checks["all"] = all(checks.values())
    return checks


def qualifying_runs(flags: Sequence[bool]) -> list[tuple[int, int]]:
    runs: list[tuple[int, int]] = []
    start: int | None = None
    for index, flag in enumerate(flags):
        if flag and start is None:
            start = index
        if start is not None and (not flag or index == len(flags) - 1):
            end = index if flag and index == len(flags) - 1 else index - 1
            if end - start + 1 >= MIN_BOUNDARY_RUN_FRAMES:
                runs.append((start, end))
            start = None
    return runs


def _last_nstep(text: str) -> int | None:
    matches = re.findall(r"\bNSTEP\s*=\s*(\d+)", text, re.I)
    return int(matches[-1]) if matches else None


def write_candidate_frame(
    trajectory: pathlib.Path, frame_index: int, target: pathlib.Path
) -> None:
    if target.exists():
        raise FileExistsError(target)
    import parmed as pmd

    structure = pmd.load_file(str(PRMTOP))
    reader = pmd.amber.AmberMdcrd(
        str(trajectory), len(structure.atoms), hasbox=True, mode="r"
    )
    coordinate_frames = reader.coordinates
    if not 0 <= frame_index < len(coordinate_frames):
        raise ValueError("candidate frame index outside trajectory")
    structure.coordinates = coordinate_frames[frame_index]
    box_frames = getattr(reader, "box", None)
    if box_frames is not None and len(box_frames) > frame_index:
        observed_box = list(box_frames[frame_index])
        if len(observed_box) == 3:
            structure.box = observed_box + [90.0, 90.0, 90.0]
        elif len(observed_box) >= 6:
            structure.box = observed_box[:6]
    structure.save(str(target), overwrite=False)


def audit_release(
    root: pathlib.Path, output: pathlib.Path, scratch: pathlib.Path
) -> dict[str, Any]:
    manifest_path = root / "SOURCE_MANIFEST.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    prepared = json.loads((output / "RELEASE_MANIFEST.json").read_text(encoding="utf-8"))
    stage_out = scratch / "stage.out"
    restart = scratch / "stage.rst7"
    trajectory = scratch / "release.mdcrd"
    text = (
        stage_out.read_text(encoding="utf-8", errors="replace")
        if stage_out.is_file()
        else ""
    )
    hard = AC._hard_hits(text)
    scc = len(re.findall(r"Convergence could not be achieved", text, re.I))
    vlimit = len(re.findall(r"vlimit\s+exceeded", text, re.I))
    overflow = len(re.findall(r"BOND\s*=\s*\*+", text, re.I))
    complete = _last_nstep(text) == RELEASE_STEPS
    restart_geometry: dict[str, Any] = {}
    if restart.is_file() and restart.stat().st_size:
        restart_geometry = AC._read_structure_geometry(restart, manifest)

    frames: list[dict[str, Any]] = []
    parse_warnings: list[str] = []
    if trajectory.is_file() and trajectory.stat().st_size:
        try:
            observed = AC._read_md_frames(trajectory, manifest)
            frames = list(observed[0])
            parse_warnings = [str(item) for item in observed[1]]
        except Exception as error:
            parse_warnings = [f"{type(error).__name__}: {error}"]
    technical = bool(
        complete
        and restart.is_file()
        and restart.stat().st_size
        and _finite_geometry(restart_geometry)
        and len(frames) == RELEASE_FRAME_COUNT
        and not parse_warnings
        and scc == 0
        and vlimit == 0
        and overflow == 0
        and sum(hard.values()) == 0
        and all(_finite_geometry(frame) for frame in frames)
    )

    frame_trace = []
    boundary_flags: list[bool] = []
    for frame_index, geometry in enumerate(frames):
        nstep = (frame_index + 1) * RELEASE_NTWX
        time_ps = nstep * RELEASE_DT_PS
        checks = boundary_checks(geometry)
        eligible_time = time_ps > RELEASE_EXCLUDED_INITIAL_PS
        boundary = bool(technical and eligible_time and checks["all"])
        boundary_flags.append(boundary)
        frame_trace.append(
            {
                "frame_index": frame_index,
                "nstep": nstep,
                "time_ps": time_ps,
                "after_excluded_initial_window": eligible_time,
                "boundary": boundary,
            }
        )

    runs = qualifying_runs(boundary_flags) if technical else []
    candidate_record: dict[str, Any] | None = None
    if runs and len(manifest["candidates"]) < MAX_CANDIDATES_PER_SEED:
        run_start, run_end = runs[0]
        selected_frame = (run_start + run_end) // 2
        candidate_number = len(manifest["candidates"]) + 1
        candidate_path = root / f"candidate_{candidate_number:02d}.rst7"
        write_candidate_frame(trajectory, selected_frame, candidate_path)
        candidate_geometry = AC._read_structure_geometry(candidate_path, manifest)
        candidate_gate = boundary_checks(candidate_geometry)
        if not candidate_gate["all"]:
            raise ValueError("extracted unrestrained candidate fails its source-frame gate")
        selected_nstep = (selected_frame + 1) * RELEASE_NTWX
        selected_time = selected_nstep * RELEASE_DT_PS
        candidate_record = {
            "candidate_number": candidate_number,
            "path": str(candidate_path),
            "sha256": sha256(candidate_path),
            "source_kind": "FULLY_UNRESTRAINED_RELEASE_PROBE_FRAME",
            "source_window_index": prepared["window_index"],
            "source_trajectory": str(trajectory),
            "trajectory_persistence": "scratch_only_not_persisted",
            "frame_index": selected_frame,
            "nstep": selected_nstep,
            "time_ps": selected_time,
            "qualifying_run": {
                "start_frame": run_start,
                "end_frame": run_end,
                "frame_count": run_end - run_start + 1,
            },
            "geometry": candidate_geometry,
            "gate": candidate_gate,
            "restraint_proof": {
                "ntr": 0,
                "nmropt": 0,
                "disang": None,
                "reactive_restraints": [],
                "environment_restraint": False,
                "excluded_initial_ps": RELEASE_EXCLUDED_INITIAL_PS,
            },
            "eligible_for_future_paired_velocity_shooting_preflight": True,
            "scientific_status": SCIENTIFIC_BOUNDARY,
        }
        manifest["candidates"].append(candidate_record)

    diagnostics = {
        "complete_to_nstep": _last_nstep(text),
        "expected_nstep": RELEASE_STEPS,
        "frame_count": len(frames),
        "expected_frame_count": RELEASE_FRAME_COUNT,
        "parse_warnings": parse_warnings,
        "scc_warnings": scc,
        "vlimit_warnings": vlimit,
        "bond_overflow": overflow,
        "hard_error_hits": hard,
        "trajectory_persisted": False,
    }
    result = {
        "schema_version": 1,
        "status": (
            "PASS_TECHNICAL_A1_PRODUCT_REVERSE_RELEASE"
            if technical
            else "NOT_EVALUATED_TECHNICAL_A1_PRODUCT_REVERSE_RELEASE"
        ),
        "scientific_status": SCIENTIFIC_BOUNDARY,
        "technical_pass": technical,
        "scientific_gate": (
            "PASS_UNRESTRAINED_PRODUCT_SIDE_REVERSE_BOUNDARY_CANDIDATE"
            if candidate_record
            else "NO_UNRESTRAINED_REVERSE_BOUNDARY_CANDIDATE_IN_RELEASE"
            if technical
            else "NOT_EVALUATED_TECHNICAL_REVERSE_RELEASE"
        ),
        "window_index": prepared["window_index"],
        "input_restart_sha256": prepared["input_restart_sha256"],
        "restart_sha256": (
            sha256(restart) if restart.is_file() and restart.stat().st_size else None
        ),
        "restart_geometry": restart_geometry,
        "candidate": candidate_record,
        "qualifying_runs": [
            {
                "start_frame": start,
                "end_frame": end,
                "frame_count": end - start + 1,
            }
            for start, end in runs
        ],
        "frame_gate_trace": frame_trace,
        "diagnostics": diagnostics,
        "shooting_forbidden_inputs": {
            "restrained_window_endpoint": True,
            "release_initial_frame": True,
        },
    }
    write_json(output / "RESULT.json", result)
    write_json(
        output / ("PASS.json" if technical else "NOT_EVALUATED.json"), result
    )
    manifest["releases"].append(result)
    write_json(manifest_path, manifest)
    if not technical:
        raise RuntimeError("unrestrained reverse release did not technically PASS")
    return result


def _inheritance_ok(manifest: Mapping[str, Any]) -> bool:
    previous = manifest["start_restart_sha256"]
    for item in manifest.get("windows", []):
        if item.get("input_restart_sha256") != previous:
            return False
        output_sha = item.get("output_restart_sha256")
        if not isinstance(output_sha, str) or not re.fullmatch(r"[0-9a-f]{64}", output_sha):
            return False
        previous = output_sha
    return True


def _write_hashes(root: pathlib.Path) -> None:
    records = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.name == "SHA256.tsv":
            continue
        if path.suffix == ".json" or path.suffix == ".rst7":
            records.append(f"{sha256(path)}  {path.relative_to(root)}")
    (root / "SHA256.tsv").write_text(
        "\n".join(records) + ("\n" if records else ""), encoding="utf-8"
    )


def finalize(
    root: pathlib.Path,
    stop_reason: str,
    technical_failure: bool = False,
) -> dict[str, Any]:
    manifest_path = root / "SOURCE_MANIFEST.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    windows = manifest.get("windows", [])
    releases = manifest.get("releases", [])
    candidates = manifest.get("candidates", [])
    windows_technical = bool(windows) and all(
        item.get("technical_pass") is True for item in windows
    )
    window_indices = [
        item.get("window_spec", {}).get("window_index") for item in windows
    ]
    window_indices_strict = window_indices == list(range(len(windows)))
    expected_release_windows = {
        index
        for index, item in zip(window_indices, windows)
        if item.get("restrained_hint") is True
    }
    observed_release_windows = [item.get("window_index") for item in releases]
    observed_release_indices_valid = all(
        isinstance(index, int)
        and not isinstance(index, bool)
        and 0 <= index < len(WINDOWS)
        for index in observed_release_windows
    )
    release_windows_unique = (
        observed_release_indices_valid
        and len(observed_release_windows) == len(set(observed_release_windows))
    )
    release_window_set_exact = (
        observed_release_indices_valid
        and set(observed_release_windows) == expected_release_windows
    )
    window_by_index = dict(zip(window_indices, windows))
    release_inputs_match_windows = (
        observed_release_indices_valid
        and all(
            item.get("input_restart_sha256")
            == window_by_index.get(item.get("window_index"), {}).get(
                "output_restart_sha256"
            )
            for item in releases
        )
    )
    releases_technical = all(
        item.get("technical_pass") is True for item in releases
    )
    release_contract = bool(
        window_indices_strict
        and release_windows_unique
        and release_window_set_exact
        and release_inputs_match_windows
        and releases_technical
    )
    inheritance = _inheritance_ok(manifest)
    terminal_complete = bool(
        (stop_reason == "COMPLETED_ALL_WINDOWS" and len(windows) == len(WINDOWS))
        or (
            stop_reason == "CANDIDATE_LIMIT_REACHED"
            and len(candidates) == MAX_CANDIDATES_PER_SEED
        )
        or (
            stop_reason == "GUARD_STOP"
            and windows
            and windows[-1].get("guard_stop") is True
        )
    )
    candidate_limit_ok = len(candidates) <= MAX_CANDIDATES_PER_SEED
    technical = bool(
        not technical_failure
        and windows_technical
        and releases_technical
        and release_contract
        and inheritance
        and terminal_complete
        and candidate_limit_ok
    )
    if not technical:
        scientific_gate = "NOT_EVALUATED_TECHNICAL_PRODUCT_SIDE_REVERSE_BOUNDARY_SCOUT"
    elif candidates:
        scientific_gate = "PASS_UNRESTRAINED_PRODUCT_SIDE_REVERSE_BOUNDARY_CANDIDATE"
    elif stop_reason == "GUARD_STOP":
        scientific_gate = "NO_CANDIDATE_CHEMICAL_GUARD_STOP"
    else:
        scientific_gate = "FAIL_NO_UNRESTRAINED_PRODUCT_SIDE_REVERSE_BOUNDARY_CANDIDATE"

    result = {
        "schema_version": 1,
        "status": TECHNICAL_PASS if technical else TECHNICAL_FAIL,
        "technical_complete": technical,
        "scientific_gate": scientific_gate,
        "scientific_status": SCIENTIFIC_BOUNDARY,
        "scope": SCOPE,
        "seed": manifest["seed"],
        "seed_index": manifest["seed_index"],
        "stop_reason": stop_reason,
        "window_count": len(windows),
        "release_count": len(releases),
        "candidate_count": len(candidates),
        "candidate_limit": MAX_CANDIDATES_PER_SEED,
        "candidates": candidates,
        "inheritance_sha256_verified": inheritance,
        "expected_release_windows": sorted(expected_release_windows),
        "observed_release_windows": observed_release_windows,
        "release_windows_unique": release_windows_unique,
        "release_window_set_exact": release_window_set_exact,
        "release_inputs_match_windows": release_inputs_match_windows,
        "release_contract_verified": release_contract,
        "candidate_limit_verified": candidate_limit_ok,
        "restrained_endpoints_shooting_forbidden": True,
        "automatic_downstream_action": "NONE",
    }
    write_json(root / "RESULT.json", result)
    selected = root / ("PASS.json" if technical else "NOT_EVALUATED.json")
    opposite = root / ("NOT_EVALUATED.json" if technical else "PASS.json")
    write_json(selected, result)
    if opposite.exists():
        opposite.unlink()
    _write_hashes(root)
    return result


def merge_if_ready(output_root: pathlib.Path, array_job: str) -> bool:
    audit = output_root / "audit" / f"nylc_a1_product_reverse_boundary_{array_job}.json"
    audit.parent.mkdir(parents=True, exist_ok=True)
    lock_path = audit.with_suffix(audit.suffix + ".lock")
    with lock_path.open("a+", encoding="utf-8") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            attempts = [
                output_root / f"attempt_{array_job}_{seed_index}"
                for seed_index in (0, 1)
            ]
            result_paths = [attempt / "RESULT.json" for attempt in attempts]
            if not all(path.is_file() for path in result_paths):
                return False
            results = [
                json.loads(path.read_text(encoding="utf-8"))
                for path in result_paths
            ]
            if {int(item.get("seed_index", -1)) for item in results} != {0, 1}:
                raise ValueError("cross-seed result lacks exact seed denominator 2")
            ordered = sorted(results, key=lambda item: item["seed_index"])
            if any(item.get("technical_complete") is not True for item in ordered):
                status = (
                    "NOT_EVALUATED_TECHNICAL_PRODUCT_SIDE_REVERSE_BOUNDARY_SCOUT"
                )
            else:
                positive = sum(
                    int(item.get("candidate_count", 0) >= 1) for item in ordered
                )
                if positive == 2:
                    status = "PASS_REPRODUCED_PRODUCT_SIDE_REVERSE_BOUNDARY_SCOUT"
                elif positive == 1:
                    status = "NOT_REPRODUCED_PRODUCT_SIDE_REVERSE_BOUNDARY_SCOUT"
                else:
                    status = (
                        "FAIL_NO_UNRESTRAINED_PRODUCT_SIDE_REVERSE_BOUNDARY_CANDIDATE"
                    )
            payload = {
                "schema_version": 1,
                "denominator_seeds": 2,
                "status": status,
                "scientific_status": SCIENTIFIC_BOUNDARY,
                "scope": SCOPE,
                "per_seed": ordered,
                "automatic_downstream_action": "NONE",
                "step2_boundary": (
                    "REVALIDATE_BOTH_BASINS_UNDER_FINAL_STEP2_QM_WATER_HAMILTONIAN"
                ),
            }
            if audit.exists():
                existing = json.loads(audit.read_text(encoding="utf-8"))
                if existing != payload:
                    raise FileExistsError(
                        f"cross-seed audit differs from existing {audit}"
                    )
                return True
            write_json(audit, payload)
            return True
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "source_authority_commit": AUTHORITY_COMMIT,
        "source_job": SOURCE_JOB,
        "seed_count": len(SOURCES),
        "window_count_per_seed": len(WINDOWS),
        "candidate_limit_per_seed": MAX_CANDIDATES_PER_SEED,
        "release_steps": RELEASE_STEPS,
        "release_dt_ps": RELEASE_DT_PS,
        "release_trajectory_policy": "scratch_only",
        "restrained_endpoint_shooting_forbidden": True,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        required=True,
        choices=(
            "describe",
            "initialize",
            "prepare-window",
            "audit-window",
            "prepare-release",
            "audit-release",
            "finalize",
            "merge-if-ready",
        ),
    )
    parser.add_argument("--seed-index", type=int)
    parser.add_argument("--window-index", type=int)
    parser.add_argument("--root", type=pathlib.Path)
    parser.add_argument("--output", type=pathlib.Path)
    parser.add_argument("--output-root", type=pathlib.Path)
    parser.add_argument("--scratch", type=pathlib.Path)
    parser.add_argument("--start-rst7", type=pathlib.Path)
    parser.add_argument("--input-rst7", type=pathlib.Path)
    parser.add_argument("--github-commit", default="unknown")
    parser.add_argument("--stop-reason", default="UNKNOWN")
    parser.add_argument("--array-job")
    parser.add_argument("--technical-failure", action="store_true")
    args = parser.parse_args()

    if args.mode == "describe":
        print(json.dumps(describe(), sort_keys=True))
    elif args.mode == "initialize":
        initialize(
            args.seed_index,
            args.root.resolve(),
            args.start_rst7.resolve(),
            args.github_commit,
        )
    elif args.mode == "prepare-window":
        prepare_window(
            args.window_index,
            args.root.resolve(),
            args.output.resolve(),
            args.scratch.resolve(),
            args.input_rst7.resolve(),
        )
    elif args.mode == "audit-window":
        audit_window(
            args.root.resolve(), args.output.resolve(), args.scratch.resolve()
        )
    elif args.mode == "prepare-release":
        prepare_release(
            args.root.resolve(),
            args.window_index,
            args.output.resolve(),
            args.scratch.resolve(),
            args.input_rst7.resolve(),
        )
    elif args.mode == "audit-release":
        audit_release(
            args.root.resolve(), args.output.resolve(), args.scratch.resolve()
        )
    elif args.mode == "finalize":
        finalize(
            args.root.resolve(), args.stop_reason, args.technical_failure
        )
    else:
        merge_if_ready(args.output_root.resolve(), args.array_job)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
