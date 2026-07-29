#!/usr/bin/env python3
"""NylC A1 Step1 balanced four-distance bridge scout.

The restrained bridge is four independent Amber distance restraints.  It is not
a linear-combination collective-variable bias.  Restrained endpoints are never
shooting inputs; only a fully reaction-coordinate-unrestrained release frame
can be recorded as a future paired-shooting preflight candidate.
"""

from __future__ import annotations

import argparse
import fcntl
import importlib.util
import json
import os
import pathlib
import shutil
from typing import Any, Mapping


WORKFLOW_ROOT = pathlib.Path(__file__).resolve().parents[1]
REVERSE_PATH = pathlib.Path(__file__).with_name(
    "prepare_audit_nylc_a1_step1_product_reverse_boundary.py"
)


def _load_module(name: str, path: pathlib.Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import required module {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


REVERSE = _load_module("_a1_product_reverse_boundary", REVERSE_PATH)
AC = REVERSE.AC
BASE = REVERSE.BASE
REACTIVE = AC.REACTIVE_ATOMS
SOURCES = REVERSE.SOURCES
PRMTOP = REVERSE.PRMTOP

BRANCHES = (
    {
        "branch_index": 0,
        "window_index": 0,
        "target_xi_bond_A": -0.25,
        "attack_A": 1.675,
        "cn_A": 1.925,
        "nalpha_hg1_A": 1.30,
        "hg1_n3_A": 1.30,
        "force_attack_cn": 15.0,
        "force_proton": 10.0,
        "maxcyc": 500,
        "ncyc": 150,
    },
    {
        "branch_index": 1,
        "window_index": 1,
        "target_xi_bond_A": 0.00,
        "attack_A": 1.800,
        "cn_A": 1.800,
        "nalpha_hg1_A": 1.30,
        "hg1_n3_A": 1.30,
        "force_attack_cn": 15.0,
        "force_proton": 10.0,
        "maxcyc": 500,
        "ncyc": 150,
    },
    {
        "branch_index": 2,
        "window_index": 2,
        "target_xi_bond_A": 0.25,
        "attack_A": 1.925,
        "cn_A": 1.675,
        "nalpha_hg1_A": 1.30,
        "hg1_n3_A": 1.30,
        "force_attack_cn": 15.0,
        "force_proton": 10.0,
        "maxcyc": 500,
        "ncyc": 150,
    },
)

MAX_CANDIDATES_PER_SEED = 1
MIN_BOUNDARY_RUN_FRAMES = 3
SCIENTIFIC_BOUNDARY = "NOT_EVALUATED_TS_COMMITTOR_PMF_BARRIER_MECHANISM"
SCOPE = "step1_balanced_four_distance_bridge_scout_not_a_linear_combination_cv"
TECHNICAL_PASS = "PASS_TECHNICAL_A1_BALANCED_FOUR_DISTANCE_BRIDGE_SCOUT"
TECHNICAL_FAIL = "NOT_EVALUATED_TECHNICAL_A1_BALANCED_FOUR_DISTANCE_BRIDGE_SCOUT"

# Reuse the frozen source/QM/proton guards and all release-frame gates exactly.
REVERSE.MAX_CANDIDATES_PER_SEED = MAX_CANDIDATES_PER_SEED
REVERSE.MIN_BOUNDARY_RUN_FRAMES = MIN_BOUNDARY_RUN_FRAMES
REVERSE.SCIENTIFIC_BOUNDARY = SCIENTIFIC_BOUNDARY
REVERSE.SCOPE = SCOPE
REVERSE.TECHNICAL_PASS = TECHNICAL_PASS
REVERSE.TECHNICAL_FAIL = TECHNICAL_FAIL


def write_json(path: pathlib.Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def sha256(path: pathlib.Path) -> str:
    return REVERSE.sha256(path)


def require_nonempty(path: pathlib.Path) -> None:
    REVERSE.require_nonempty(path)


def source_from_index(seed_index: int) -> dict[str, Any]:
    return REVERSE.source_from_index(seed_index)


def branch_spec(branch_index: int) -> dict[str, Any]:
    if not 0 <= int(branch_index) < len(BRANCHES):
        raise ValueError("branch index must be 0, 1, or 2")
    return dict(BRANCHES[int(branch_index)])


def bridge_restraints(spec: Mapping[str, Any]) -> str:
    """Return exactly four independent distance restraints, never a joint CV."""
    return "".join(
        (
            AC.distance_restraint(
                REACTIVE["og1"], REACTIVE["c12"], float(spec["attack_A"]),
                float(spec["force_attack_cn"]),
            ),
            AC.distance_restraint(
                REACTIVE["c12"], REACTIVE["n3"], float(spec["cn_A"]),
                float(spec["force_attack_cn"]),
            ),
            AC.distance_restraint(
                REACTIVE["nalpha"], REACTIVE["hg1"], float(spec["nalpha_hg1_A"]),
                float(spec["force_proton"]),
            ),
            AC.distance_restraint(
                REACTIVE["hg1"], REACTIVE["n3"], float(spec["hg1_n3_A"]),
                float(spec["force_proton"]),
            ),
        )
    )


def minimization_input(seed: str, spec: Mapping[str, Any], qmmask: str) -> str:
    return f"""NylC A1 balanced four-distance bridge {seed} branch {spec['branch_index']}
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


def initialize(
    seed_index: int, branch_index: int, root: pathlib.Path,
    start_rst7: pathlib.Path, github_commit: str,
) -> None:
    if root.exists():
        raise FileExistsError(f"attempt root already exists: {root}")
    if start_rst7.exists() or not start_rst7.parent.is_dir():
        raise ValueError("scratch start restart target must be fresh")
    source_manifest, qmmask, geometry = REVERSE.validate_authority(seed_index)
    source = source_from_index(seed_index)
    spec = branch_spec(branch_index)
    shutil.copy2(source["restart"], start_rst7)
    if sha256(start_rst7) != source["restart_sha256"]:
        raise ValueError("scratch copy of frozen product restart changed")
    root.mkdir(parents=True, exist_ok=False)
    write_json(root / "SOURCE_MANIFEST.json", {
        "schema_version": 1,
        "status": "READY_A1_BALANCED_FOUR_DISTANCE_BRIDGE_SCOUT",
        "scientific_status": SCIENTIFIC_BOUNDARY,
        "scope": SCOPE,
        "github_commit": github_commit,
        "seed_index": seed_index,
        "seed": source["seed"],
        "branch_index": branch_index,
        "source": dict(source),
        "prmtop": str(PRMTOP),
        "prmtop_sha256": sha256(PRMTOP),
        "qm_contract": dict(source_manifest["qm_contract"]),
        "source_geometry": geometry,
        "start_restart": str(start_rst7),
        "start_restart_sha256": sha256(start_rst7),
        "branch_spec": spec,
        "windows": [],
        "releases": [],
        "candidates": [],
        "candidate_limit_per_seed": MAX_CANDIDATES_PER_SEED,
        "bias_representation": "four_independent_distance_restraints_not_a_linear_combination_cv",
        "bias_distances": (
            "OG1-C12, C12-N3, Nalpha-HG1, HG1-N3"
        ),
        "target_qpt_A_approximation": (
            "two independent proton-distance targets Nalpha-HG1=HG1-N3=1.30 A"
        ),
        "restrained_endpoint_shooting_forbidden": True,
        "automatic_downstream_action": "NONE",
    })


def prepare_window(
    root: pathlib.Path, output: pathlib.Path, scratch: pathlib.Path,
    input_rst7: pathlib.Path,
) -> None:
    if output.exists() or scratch.exists():
        raise FileExistsError("bridge output or scratch directory already exists")
    require_nonempty(input_rst7)
    manifest = json.loads((root / "SOURCE_MANIFEST.json").read_text(encoding="utf-8"))
    if manifest["windows"]:
        raise ValueError("each independent task permits exactly one bridge window")
    if sha256(input_rst7) != manifest["start_restart_sha256"]:
        raise ValueError("bridge input restart differs from its fixed source restart")
    spec = manifest["branch_spec"]
    output.mkdir(parents=True, exist_ok=False)
    scratch.mkdir(parents=True, exist_ok=False)
    (scratch / "stage.in").write_text(
        minimization_input(manifest["seed"], spec, manifest["qm_contract"]["qmmask"]),
        encoding="utf-8",
    )
    (scratch / "restraints.RST").write_text(bridge_restraints(spec), encoding="utf-8")
    prepared = {
        "schema_version": 1,
        "status": "READY_RESTRAINED_BALANCED_FOUR_DISTANCE_BRIDGE",
        "seed": manifest["seed"],
        "branch_index": manifest["branch_index"],
        "window_spec": spec,
        "input_restart": str(input_rst7),
        "input_restart_sha256": manifest["start_restart_sha256"],
        "restrained": True,
        "shooting_forbidden": True,
        "bias_representation": manifest["bias_representation"],
    }
    write_json(output / "WINDOW_MANIFEST.json", prepared)
    write_json(scratch / "PREPARED.json", prepared)


def audit_window(root: pathlib.Path, output: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    # Existing technical, chemical, direction, and restrained-hint gates are reused unchanged.
    return REVERSE.audit_window(root, output, scratch)


def prepare_release(
    root: pathlib.Path, output: pathlib.Path, scratch: pathlib.Path,
    input_rst7: pathlib.Path,
) -> None:
    manifest = json.loads((root / "SOURCE_MANIFEST.json").read_text(encoding="utf-8"))
    REVERSE.prepare_release(
        root, int(manifest["branch_index"]), output, scratch, input_rst7
    )


def audit_release(root: pathlib.Path, output: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    # Existing 0.25 ps fully unrestrained release, boundary_checks, and 3-frame run gate.
    return REVERSE.audit_release(root, output, scratch)


def finalize(root: pathlib.Path, stop_reason: str, technical_failure: bool = False) -> dict[str, Any]:
    manifest_path = root / "SOURCE_MANIFEST.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    windows = manifest.get("windows", [])
    releases = manifest.get("releases", [])
    candidates = manifest.get("candidates", [])
    window = windows[0] if len(windows) == 1 else {}
    hint = window.get("restrained_hint") is True
    release_contract = (
        len(releases) == (1 if hint else 0)
        and all(item.get("technical_pass") is True for item in releases)
    )
    technical = bool(
        not technical_failure
        and len(windows) == 1
        and window.get("technical_pass") is True
        and release_contract
        and len(candidates) <= 1
        and stop_reason in {"COMPLETED_ONE_BRANCH", "CANDIDATE_LIMIT_REACHED", "GUARD_STOP"}
    )
    result = {
        "schema_version": 1,
        "status": TECHNICAL_PASS if technical else TECHNICAL_FAIL,
        "technical_complete": technical,
        "scientific_status": SCIENTIFIC_BOUNDARY,
        "scope": SCOPE,
        "seed": manifest["seed"],
        "seed_index": manifest["seed_index"],
        "branch_index": manifest["branch_index"],
        "stop_reason": stop_reason,
        "restrained_hint": hint,
        "candidate": candidates[0] if candidates else None,
        "candidate_count_for_this_independent_branch": len(candidates),
        "candidate_limit_per_seed": MAX_CANDIDATES_PER_SEED,
        "bias_representation": manifest["bias_representation"],
        "restrained_endpoints_shooting_forbidden": True,
        "automatic_downstream_action": "NONE",
    }
    write_json(root / "RESULT.json", result)
    selected = root / ("PASS.json" if technical else "NOT_EVALUATED.json")
    opposite = root / ("NOT_EVALUATED.json" if technical else "PASS.json")
    write_json(selected, result)
    if opposite.exists():
        opposite.unlink()
    REVERSE._write_hashes(root)
    return result


def merge_if_ready(output_root: pathlib.Path, array_job: str, seed_index: int) -> bool:
    """Write one deterministic per-seed selection only after all three branches finish."""
    audit = output_root / "audit" / (
        f"nylc_a1_balanced_four_distance_bridge_{array_job}_seed{seed_index}.json"
    )
    audit.parent.mkdir(parents=True, exist_ok=True)
    with (audit.with_suffix(audit.suffix + ".lock")).open("a+", encoding="utf-8") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        try:
            task_indices = [seed_index * 3 + branch for branch in range(3)]
            paths = [output_root / f"attempt_{array_job}_{index}" / "RESULT.json" for index in task_indices]
            if not all(path.is_file() for path in paths):
                return False
            results = [json.loads(path.read_text(encoding="utf-8")) for path in paths]
            if [item.get("branch_index") for item in sorted(results, key=lambda item: item["branch_index"])] != [0, 1, 2]:
                raise ValueError("per-seed bridge denominator is not exactly three branches")
            if any(item.get("seed_index") != seed_index for item in results):
                raise ValueError("per-seed bridge result has a mismatched seed")
            technical = all(item.get("technical_complete") is True for item in results)
            proposals = [item for item in sorted(results, key=lambda item: item["branch_index"]) if item.get("candidate")]
            selected = proposals[0].get("candidate") if technical and proposals else None
            payload = {
                "schema_version": 1,
                "status": (
                    "NOT_EVALUATED_TECHNICAL_A1_BALANCED_FOUR_DISTANCE_BRIDGE_SCOUT"
                    if not technical else
                    "PASS_ONE_TIME_DECORRELATED_BRIDGE_CANDIDATE"
                    if selected else
                    "FAIL_NO_UNRESTRAINED_BRIDGE_CANDIDATE"
                ),
                "denominator_branches": 3,
                "seed_index": seed_index,
                "scientific_status": SCIENTIFIC_BOUNDARY,
                "scope": SCOPE,
                "bias_representation": "four_independent_distance_restraints_not_a_linear_combination_cv",
                "per_branch": sorted(results, key=lambda item: item["branch_index"]),
                "selected_candidate": selected,
                "candidate_count_for_seed": 1 if selected else 0,
                "selection_policy": "first branch-indexed proposal after existing continuous-3-frame release gate",
                "automatic_downstream_action": "NONE",
            }
            if audit.exists():
                if json.loads(audit.read_text(encoding="utf-8")) != payload:
                    raise FileExistsError(f"existing per-seed audit differs: {audit}")
                return True
            write_json(audit, payload)
            return True
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "seed_count": 2,
        "branch_count": len(BRANCHES),
        "array_task_count": 6,
        "array_mapping": "task_index=seed_index*3+branch_index",
        "branch_specs": list(BRANCHES),
        "candidate_limit_per_seed": MAX_CANDIDATES_PER_SEED,
        "continuous_frame_requirement": MIN_BOUNDARY_RUN_FRAMES,
        "release_steps": REVERSE.RELEASE_STEPS,
        "release_dt_ps": REVERSE.RELEASE_DT_PS,
        "bias_representation": "four_independent_distance_restraints_not_a_linear_combination_cv",
        "restrained_endpoint_shooting_forbidden": True,
        "automatic_downstream_action": "NONE",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=(
        "describe", "initialize", "prepare-window", "audit-window",
        "prepare-release", "audit-release", "finalize", "merge-if-ready",
    ))
    parser.add_argument("--seed-index", type=int)
    parser.add_argument("--branch-index", type=int)
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
        initialize(args.seed_index, args.branch_index, args.root.resolve(), args.start_rst7.resolve(), args.github_commit)
    elif args.mode == "prepare-window":
        prepare_window(args.root.resolve(), args.output.resolve(), args.scratch.resolve(), args.input_rst7.resolve())
    elif args.mode == "audit-window":
        audit_window(args.root.resolve(), args.output.resolve(), args.scratch.resolve())
    elif args.mode == "prepare-release":
        prepare_release(args.root.resolve(), args.output.resolve(), args.scratch.resolve(), args.input_rst7.resolve())
    elif args.mode == "audit-release":
        audit_release(args.root.resolve(), args.output.resolve(), args.scratch.resolve())
    elif args.mode == "finalize":
        finalize(args.root.resolve(), args.stop_reason, args.technical_failure)
    else:
        merge_if_ready(args.output_root.resolve(), args.array_job, args.seed_index)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
