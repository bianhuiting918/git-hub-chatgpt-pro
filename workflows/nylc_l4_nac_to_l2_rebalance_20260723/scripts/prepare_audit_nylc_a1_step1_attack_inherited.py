#!/usr/bin/env python3
"""Inherited, attack-only restrained A1 minimisation chain; not a PMF protocol."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import pathlib
import re
from typing import Any

_BASE_PATH = pathlib.Path(__file__).with_name("prepare_audit_nylc_a1_step1_pt2_cn_scout.py")
_spec = importlib.util.spec_from_file_location("_a1_pt2_cn_base", _BASE_PATH)
if _spec is None or _spec.loader is None:
    raise RuntimeError(f"cannot import immutable base {_BASE_PATH}")
BASE = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(BASE)

ATTACK_TARGETS_A = (3.0, 2.8, 2.6, 2.4, 2.2, 2.0, 1.8, 1.65, 1.50)
ATTACK_FORCE_KCAL_MOL_A2 = (10.0, 10.0, 10.0, 10.0, 25.0, 25.0, 25.0, 50.0, 50.0)
MAXCYC = (300, 300, 400, 400, 500, 600, 700, 800, 1000)
NCYC = (75, 75, 100, 100, 125, 150, 175, 200, 250)
DRMS = 0.10
PT2_AND_CN_ARE_AUDIT_ONLY = True
PASS_STATUS = "PASS_TECHNICAL_A1_ATTACK_INHERITED_CHAIN"
FAIL_STATUS = "NOT_EVALUATED_A1_ATTACK_INHERITED_CHAIN"
SCIENTIFIC_STATUS = "NOT_EVALUATED_TS_PMF_BARRIER_MECHANISM"
SCIENTIFIC_ACTION = "DO_NOT_START_PMF"
SCOPE = "inherited_restrained_attack_closure_not_ts_pmf_barrier_or_mechanism"

# Immutable atom mapping inherited from the frozen Step-1 preflight authority.
THR267_N = BASE.THR267_N
THR267_OG1 = BASE.THR267_OG1
TRANSFERRED_HG1 = BASE.TRANSFERRED_HG1
L2_C12 = BASE.L2_C12
L2_O2 = BASE.L2_O2
L2_N3 = BASE.L2_N3
C12_C11 = BASE.C12_C11


def write_json(path: pathlib.Path, payload: dict[str, Any]) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def seed_from_index(index: int) -> dict[str, Any]:
    if not 0 <= int(index) < len(BASE.SOURCES):
        raise ValueError(f"seed array index {index} is outside 0..1")
    return dict(BASE.SOURCES[int(index)])


def window_from_index(seed_index: int, window_index: int) -> dict[str, Any]:
    source = seed_from_index(seed_index)
    if not 0 <= int(window_index) < len(ATTACK_TARGETS_A):
        raise ValueError(f"window index {window_index} is outside 0..8")
    i = int(window_index)
    return {
        "seed_index": int(seed_index),
        "seed": source["seed"],
        "candidate": source["candidate"],
        "window_index": i,
        "attack_target_A": ATTACK_TARGETS_A[i],
        "attack_force_kcal_mol_A2": ATTACK_FORCE_KCAL_MOL_A2[i],
        "maxcyc": MAXCYC[i],
        "ncyc": NCYC[i],
        "tag": f"{source['seed']}_attack_{i:02d}_{ATTACK_TARGETS_A[i]:.2f}A",
    }


def _xyz(structure: Any, atom_index: int) -> tuple[float, float, float]:
    atom = structure.atoms[atom_index - 1]
    return atom.xx, atom.xy, atom.xz


def _distance(structure: Any, first: int, second: int) -> float:
    left, right = _xyz(structure, first), _xyz(structure, second)
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right)))


def _angle(structure: Any, first: int, center: int, third: int) -> float:
    a, b, c = (_xyz(structure, i) for i in (first, center, third))
    u = tuple(x - y for x, y in zip(a, b))
    v = tuple(x - y for x, y in zip(c, b))
    denom = math.sqrt(sum(x * x for x in u)) * math.sqrt(sum(x * x for x in v))
    if denom == 0:
        return float("nan")
    cosine = sum(x * y for x, y in zip(u, v)) / denom
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def _point_plane_distance(point: tuple[float, float, float],
                          a: tuple[float, float, float],
                          b: tuple[float, float, float],
                          c: tuple[float, float, float]) -> float:
    ab = tuple(y - x for x, y in zip(a, b))
    ac = tuple(y - x for x, y in zip(a, c))
    normal = (
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    )
    norm = math.sqrt(sum(x * x for x in normal))
    if norm == 0:
        return float("nan")
    return abs(sum(n * (p - x) for n, p, x in zip(normal, point, a))) / norm


def geometry(structure: Any) -> dict[str, float]:
    n_h = _distance(structure, THR267_N, TRANSFERRED_HG1)
    h_n3 = _distance(structure, TRANSFERRED_HG1, L2_N3)
    c12 = _xyz(structure, L2_C12)
    c11 = _xyz(structure, C12_C11)
    o2 = _xyz(structure, L2_O2)
    n3 = _xyz(structure, L2_N3)
    return {
        "attack_A": _distance(structure, THR267_OG1, L2_C12),
        "nalpha_hg1_A": n_h,
        "hg1_n3_A": h_n3,
        "qPT_A": n_h - h_n3,
        "qCN_A": _distance(structure, L2_C12, L2_N3),
        "c12_o2_A": _distance(structure, L2_C12, L2_O2),
        "attack_angle_deg": _angle(structure, L2_O2, L2_C12, THR267_OG1),
        "c12_out_of_plane_A": _point_plane_distance(c12, o2, n3, c11),
    }


def attack_restraint(target: float, force: float) -> str:
    """The only reaction-coordinate restraint in this protocol."""
    return (
        f"&rst iat=8960,10287, r1={max(0.1, target - 0.35):.3f}, "
        f"r2={target - 0.05:.3f}, r3={target + 0.05:.3f}, r4=4.500, "
        f"rk2={force:.1f}, rk3={force:.1f}, /\n"
    )


def amber_input(title: str, qmmask: str, window: dict[str, Any]) -> str:
    return f"""{title}
&cntrl
  imin=1, ntmin=2, maxcyc={window['maxcyc']}, ncyc={window['ncyc']}, dx0=0.005,
  ntb=1, cut=10.0, ntpr=5, ntxo=1,
  ifqnt=1, nmropt=1, ntr=1,
  restraint_wt=1.0,
  restraintmask='{BASE.NON_QM_SOLUTE_HEAVY_MASK}',
  drms={DRMS:.2f},
/
{BASE.qmmm_block(qmmask)}&wt type='END' /
DISANG=restraints.RST
DUMPAVE=restraint.dat
"""


def _load_geometry(prmtop: pathlib.Path, rst7: pathlib.Path) -> dict[str, float]:
    import parmed as pmd
    return geometry(pmd.load_file(str(prmtop), xyz=str(rst7)))


def initialize(seed_index: int, start_rst7: pathlib.Path, chain_state: pathlib.Path,
               github_commit: str) -> None:
    source = seed_from_index(seed_index)
    qmmask, source_gro = BASE.validate_authority(source)
    import parmed as pmd
    structure = pmd.load_file(str(BASE.PRMTOP), xyz=str(source_gro))
    if len(structure.atoms) != BASE.EXPECTED_SYSTEM_ATOMS or structure.box is None:
        raise ValueError("immutable coordinate transplant changed atom count or box")
    for atom_index, resname, name in (
        (THR267_N, "THR", "N"), (THR267_OG1, "THR", "OG1"),
        (TRANSFERRED_HG1, "THR", "HG1"), (L2_C12, "L2", "C12"),
        (L2_O2, "L2", "O2"), (L2_N3, "L2", "N3"),
    ):
        atom = structure.atoms[atom_index - 1]
        if atom.residue.name != resname or atom.name != name:
            raise ValueError(f"frozen atom-order identity changed at {atom_index}")
    if start_rst7.exists():
        raise FileExistsError(start_rst7)
    start_rst7.parent.mkdir(parents=True, exist_ok=True)
    structure.save(str(start_rst7), overwrite=False)
    start_geometry = geometry(structure)
    write_json(chain_state, {
        "schema_version": 1,
        "status": "READY_A1_ATTACK_INHERITED_CHAIN",
        "github_commit": github_commit,
        "source": {
            "seed": source["seed"], "candidate": source["candidate"],
            "source_gro": str(source_gro), "source_gro_sha256": sha256(source_gro),
            "frozen_prmtop": str(BASE.PRMTOP),
            "frozen_prmtop_sha256": sha256(BASE.PRMTOP),
            "coordinate_transplant": "ParMed frozen prmtop plus same-order source GRO",
        },
        "qm_contract": {
            "qm_atom_count": BASE.EXPECTED_QM_ATOMS, "qmcharge": BASE.QMCHARGE,
            "electron_count_including_link_h": BASE.EXPECTED_ELECTRONS,
            "link_atom_count": BASE.EXPECTED_LINK_ATOMS, "step1_qm_water_count": 0,
            "qmmask": qmmask,
        },
        "baseline_3p0_geometry": start_geometry,
        "scientific_status": SCIENTIFIC_STATUS, "next_action": SCIENTIFIC_ACTION,
        "scope": SCOPE,
    })


def prepare(seed_index: int, window_index: int, output: pathlib.Path,
            start_rst7: pathlib.Path, chain_state: pathlib.Path,
            github_commit: str) -> None:
    if output.exists():
        raise FileExistsError(output)
    if not start_rst7.is_file() or start_rst7.stat().st_size == 0:
        raise ValueError(f"missing inherited restart {start_rst7}")
    state = json.loads(chain_state.read_text(encoding="utf-8"))
    if state.get("status") != "READY_A1_ATTACK_INHERITED_CHAIN":
        raise ValueError("chain state is not READY")
    window = window_from_index(seed_index, window_index)
    output.mkdir(parents=True)
    (output / "stage.in").write_text(
        amber_input(f"NylC A1 inherited attack-only {window['tag']}",
                    state["qm_contract"]["qmmask"], window), encoding="utf-8")
    (output / "restraints.RST").write_text(
        attack_restraint(window["attack_target_A"], window["attack_force_kcal_mol_A2"]),
        encoding="utf-8")
    manifest = {
        "schema_version": 1, "status": "READY_A1_ATTACK_INHERITED_WINDOW",
        "github_commit": github_commit, "window": window,
        "inherited_restart": str(start_rst7),
        "frozen_prmtop": state["source"]["frozen_prmtop"],
        "frozen_prmtop_sha256": state["source"]["frozen_prmtop_sha256"],
        "baseline_3p0_geometry": state["baseline_3p0_geometry"],
        "protocol": {
            "reaction_restraint": "attack_only_8960_10287",
            "attack_angle_restrained": False,
            "pt2_restrained": False, "cn_restrained": False,
            "PT2_AND_CN_ARE_AUDIT_ONLY": PT2_AND_CN_ARE_AUDIT_ONLY,
            "non_qm_solute_heavy_position_restraint": BASE.NON_QM_SOLUTE_HEAVY_MASK,
            "restart_inheritance": "strict_serial_previous_window_restart",
        },
        "scientific_status": SCIENTIFIC_STATUS, "next_action": SCIENTIFIC_ACTION,
        "scope": SCOPE,
    }
    write_json(output / "WINDOW_MANIFEST.json", manifest)


def inspect_stage(prmtop: pathlib.Path, directory: pathlib.Path) -> dict[str, Any]:
    output = directory / "stage.out"
    restart = directory / "stage.rst7"
    text = output.read_text(encoding="utf-8", errors="replace") if output.is_file() else ""
    hard = {name: len(re.findall(pattern, text, re.I))
            for name, pattern in BASE.HARD_PATTERNS.items()}
    complete = "FINAL RESULTS" in text and bool(re.search(r"Run\s+done", text))
    scc = len(re.findall(r"Convergence could not be achieved", text, re.I))
    vlimit = len(re.findall(r"vlimit\s+exceeded", text, re.I))
    bond_overflow = len(re.findall(r"BOND\s*=\s*\*+", text, re.I))
    records = re.findall(
        r"^\s*(\d+)\s+(-?\d+\.\d+E[+-]\d+)\s+(\d+\.\d+E[+-]\d+)\s+"
        r"(\d+\.\d+E[+-]\d+)", text, re.M)
    observed: dict[str, float] = {}
    if restart.is_file() and restart.stat().st_size:
        observed = _load_geometry(prmtop, restart)
    finite = bool(observed) and all(math.isfinite(value) for value in observed.values())
    technical = bool(complete and scc == 0 and vlimit == 0 and bond_overflow == 0
                     and sum(hard.values()) == 0 and finite)
    return {
        "complete": complete, "scc_warnings": scc, "vlimit_warnings": vlimit,
        "bond_overflow": bond_overflow, "hard_error_hits": hard,
        "geometry": observed, "technical_pass": technical,
        "last_minimization_record": ({"nstep": int(records[-1][0]),
            "energy_kcal_mol": float(records[-1][1]), "rms": float(records[-1][2]),
            "gmax": float(records[-1][3])} if records else None),
    }


def audit(output: pathlib.Path, scratch: pathlib.Path) -> None:
    manifest = json.loads((output / "WINDOW_MANIFEST.json").read_text(encoding="utf-8"))
    if manifest.get("status") != "READY_A1_ATTACK_INHERITED_WINDOW":
        raise ValueError("window manifest is not READY")
    stage = inspect_stage(pathlib.Path(manifest["frozen_prmtop"]), scratch)
    window = manifest["window"]
    final = stage["geometry"]
    technical = stage["technical_pass"]
    residual = abs(final["attack_A"] - window["attack_target_A"]) if technical else None
    baseline_co = manifest["baseline_3p0_geometry"]["c12_o2_A"]
    elongation = final["c12_o2_A"] - baseline_co if technical else None
    reactant_guard = bool(technical and final["nalpha_hg1_A"] <= 1.20
                          and final["qPT_A"] <= -0.40 and final["qCN_A"] <= 1.50)
    guard_reason = ("OBSERVED_COUPLED_RESPONSE_OUTSIDE_REACTANT_GUARD"
                    if technical and not reactant_guard else None)
    restrained_close = bool(technical and final["attack_A"] <= 1.70
                            and residual is not None and residual <= 0.20)
    tetrahedral_like = bool(
        restrained_close and elongation is not None and elongation >= 0.06
        and 1.28 <= final["c12_o2_A"] <= 1.45
        and final["c12_out_of_plane_A"] >= 0.20
        and 1.30 <= final["qCN_A"] <= 1.60
        and 90.0 <= final["attack_angle_deg"] <= 130.0 and reactant_guard)
    forced_close = bool(restrained_close and not tetrahedral_like)
    classification = ("TETRAHEDRAL_LIKE_RESTRAINED" if tetrahedral_like else
                      "FORCED_CLOSE_CONTACT" if forced_close else
                      guard_reason if guard_reason else "NO_TETRAHEDRAL_LIKE_RESPONSE")
    result = {
        "schema_version": 1, "status": PASS_STATUS if technical else FAIL_STATUS,
        "scientific_status": SCIENTIFIC_STATUS, "next_action": SCIENTIFIC_ACTION,
        "scope": SCOPE, "window": window, "classification": classification,
        "stopped_by_reactant_guard": bool(guard_reason),
        "gates": {
            "technical": technical, "target_residual_le_0p20_A": bool(residual is not None and residual <= 0.20),
            "attack_le_1p70_A": bool(technical and final["attack_A"] <= 1.70),
            "carbonyl_elongation_from_baseline_A": elongation,
            "c12_o2_1p28_to_1p45_A": bool(technical and 1.28 <= final["c12_o2_A"] <= 1.45),
            "c12_out_of_plane_ge_0p20_A": bool(technical and final["c12_out_of_plane_A"] >= 0.20),
            "qCN_1p30_to_1p60_A": bool(technical and 1.30 <= final["qCN_A"] <= 1.60),
            "attack_angle_90_to_130_deg": bool(technical and 90.0 <= final["attack_angle_deg"] <= 130.0),
            "reactant_guard": reactant_guard,
            "TETRAHEDRAL_LIKE_RESTRAINED": tetrahedral_like,
            "FORCED_CLOSE_CONTACT": forced_close,
        },
        "target_residual_A": residual, "observed": final, "stage": stage,
    }
    write_json(output / "RESULT.json", result)
    write_json(output / ("PASS.json" if technical else "NOT_EVALUATED.json"), result)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=("initialize", "prepare", "audit"))
    parser.add_argument("--seed-index", type=int, required=True)
    parser.add_argument("--window-index", type=int)
    parser.add_argument("--output", type=pathlib.Path)
    parser.add_argument("--scratch", type=pathlib.Path)
    parser.add_argument("--start-rst7", type=pathlib.Path)
    parser.add_argument("--chain-state", type=pathlib.Path, required=True)
    parser.add_argument("--github-commit", default="unknown")
    args = parser.parse_args()
    if args.mode == "initialize":
        if args.start_rst7 is None:
            parser.error("--start-rst7 is required for initialize")
        initialize(args.seed_index, args.start_rst7.resolve(), args.chain_state.resolve(), args.github_commit)
    elif args.mode == "prepare":
        if args.window_index is None or args.output is None or args.start_rst7 is None:
            parser.error("--window-index, --output, and --start-rst7 are required for prepare")
        prepare(args.seed_index, args.window_index, args.output.resolve(), args.start_rst7.resolve(),
                args.chain_state.resolve(), args.github_commit)
    else:
        if args.output is None or args.scratch is None:
            parser.error("--output and --scratch are required for audit")
        audit(args.output.resolve(), args.scratch.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
