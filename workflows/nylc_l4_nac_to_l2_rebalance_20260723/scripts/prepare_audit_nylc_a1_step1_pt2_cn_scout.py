#!/usr/bin/env python3
"""Prepare and audit one dual-seed A1 q_attack x q_PT2 x q_CN scout window."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import pathlib
import re
from typing import Any

TASK_ROOT = pathlib.Path(
    "/work/home/acshdt1dks/nylon_pa66_scnet_20260708/"
    "l4_nac_to_l2_rebalance_20260723"
)
FRAME_ROOT = (
    TASK_ROOT
    / "a1_activated_nac_20260726/pt2_preorganized_frame_extraction/"
    "attempt_62112503"
)
QMMM_ROOT = (
    TASK_ROOT
    / "a1_activated_nac_20260726/qmmm/"
    "a1_unified_core_dftb3_preflight/attempt_62011285"
)
PRMTOP = QMMM_ROOT / "prepared/system.prmtop"
QMMM_AUDIT = QMMM_ROOT / "prepared/qmmm_preflight_audit.json"
QMMM_PASS = QMMM_ROOT / "PASS.json"
EXPECTED_PRMTOP_SHA256 = (
    "a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0"
)
SOURCES = (
    {
        "seed": "seed26723",
        "candidate": "seed26723_t378_f189",
        "gro_sha256": "14477791ce14a35cef0adf9b802b562e091660526ca06de1132f6a74070faf10",
    },
    {
        "seed": "seed26737",
        "candidate": "seed26737_t676_f338",
        "gro_sha256": "a7924184ad3db4c13e0eab4929d46ee99621eacd523e899c0aca39c540350bc8",
    },
)

ATTACK_TARGETS_A = (2.35, 1.85)
PT2_TARGETS_A = ((1.10, 1.60), (1.35, 1.35))
CN_TARGETS_A = (1.40, 1.60)
GUIDE_STEPS = 200
TARGET_STEPS = 600
GUIDE_FORCE = 5.0
TARGET_FORCE = 10.0
TARGET_TOLERANCE_A = 0.15
PT2_Q_TOLERANCE_A = 0.20
RMS_GRADIENT_MAX = 0.5
ATTACK_ANGLE_RANGE_DEG = (90.0, 130.0)

EXPECTED_SYSTEM_ATOMS = 133589
EXPECTED_QM_ATOMS = 146
EXPECTED_ELECTRONS = 510
EXPECTED_LINK_ATOMS = 6
QMCHARGE = 0
SPIN = 1

THR267_N = 8949
THR267_OG1 = 8960
TRANSFERRED_HG1 = 8961
L2_C12 = 10287
L2_O2 = 10288
L2_N3 = 10289
C12_C11 = 10286
QM_RANGES = (
    "@7160-7174,7756-7771,8235-8242,8949-8963,"
    "9567-9573,9587-9592,10273-10351"
)
NON_QM_SOLUTE_HEAVY_MASK = f"!(:SOL,NA,CL)&!@H=&!({QM_RANGES})"
HARD_PATTERNS = {
    "sander_bomb": r"SANDER BOMB",
    "segmentation": r"segmentation",
    "forrtl": r"forrtl",
    "nan": r"\bnan\b",
    "fatal": r"FATAL",
}
PASS_STATUS = "PASS_TECHNICAL_A1_PT2_CN_SCOUT"
FAIL_STATUS = "NOT_EVALUATED_A1_PT2_CN_SCOUT"
SCOPE = "bounded_restrained_scout_not_ts_pmf_barrier_or_mechanism"


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: pathlib.Path, payload: dict[str, Any]) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def window_from_index(index: int) -> dict[str, Any]:
    if not 0 <= int(index) < 16:
        raise ValueError(f"array index {index} is outside 0..15")
    seed_index = int(index) // 8
    within = int(index) % 8
    attack_index = within // 4
    pt2_index = (within % 4) // 2
    cn_index = within % 2
    source = SOURCES[seed_index]
    nh_target, hn3_target = PT2_TARGETS_A[pt2_index]
    return {
        "array_index": int(index),
        "seed": source["seed"],
        "candidate": source["candidate"],
        "gro_sha256": source["gro_sha256"],
        "attack_index": attack_index,
        "pt2_index": pt2_index,
        "cn_index": cn_index,
        "attack_target_A": ATTACK_TARGETS_A[attack_index],
        "nh_target_A": nh_target,
        "hn3_target_A": hn3_target,
        "pt2_q_target_A": nh_target - hn3_target,
        "cn_target_A": CN_TARGETS_A[cn_index],
        "tag": (
            f"{source['seed']}_a{attack_index}_p{pt2_index}_c{cn_index}"
        ),
    }


def _distance(structure: Any, index1: int, index2: int) -> float:
    left, right = structure.atoms[index1 - 1], structure.atoms[index2 - 1]
    return math.sqrt(
        (left.xx - right.xx) ** 2
        + (left.xy - right.xy) ** 2
        + (left.xz - right.xz) ** 2
    )


def _angle(structure: Any, first: int, center: int, third: int) -> float:
    a, b, c = [structure.atoms[index - 1] for index in (first, center, third)]
    u = (a.xx - b.xx, a.xy - b.xy, a.xz - b.xz)
    v = (c.xx - b.xx, c.xy - b.xy, c.xz - b.xz)
    denominator = math.sqrt(sum(value * value for value in u)) * math.sqrt(
        sum(value * value for value in v)
    )
    cosine = sum(x * y for x, y in zip(u, v)) / denominator
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def geometry(structure: Any) -> dict[str, float]:
    n_h = _distance(structure, THR267_N, TRANSFERRED_HG1)
    h_n3 = _distance(structure, TRANSFERRED_HG1, L2_N3)
    carbonyl_angles = (
        _angle(structure, C12_C11, L2_C12, L2_O2),
        _angle(structure, L2_O2, L2_C12, L2_N3),
        _angle(structure, L2_N3, L2_C12, C12_C11),
    )
    return {
        "q_attack_A": _distance(structure, THR267_OG1, L2_C12),
        "nalpha_hg1_A": n_h,
        "hg1_n3_A": h_n3,
        "pt2_q_A": n_h - h_n3,
        "q_cn_A": _distance(structure, L2_C12, L2_N3),
        "c12_o2_A": _distance(structure, L2_C12, L2_O2),
        "attack_angle_deg": _angle(
            structure, L2_O2, L2_C12, THR267_OG1
        ),
        "carbonyl_pyramidalization_deg": 360.0 - sum(carbonyl_angles),
    }


def validate_authority(window: dict[str, Any]) -> tuple[str, pathlib.Path]:
    source_root = FRAME_ROOT / window["candidate"]
    source_pass = source_root / "PASS.json"
    source_manifest = source_root / "manifest.json"
    source_gro = source_root / "source.gro"
    for path in (
        source_pass, source_manifest, source_gro, QMMM_PASS, QMMM_AUDIT, PRMTOP
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"missing immutable source {path}")
    if sha256(source_gro) != window["gro_sha256"]:
        raise ValueError("source GRO SHA256 changed")
    source_status = json.loads(source_pass.read_text(encoding="utf-8"))
    if source_status.get("status") != "PASS_A1_PT2_FRAME_EXTRACTION_AUDIT":
        raise ValueError("source frame extraction is not PASS")
    if sha256(PRMTOP) != EXPECTED_PRMTOP_SHA256:
        raise ValueError("frozen system.prmtop SHA256 changed")
    authority = json.loads(QMMM_PASS.read_text(encoding="utf-8"))
    audit = json.loads(QMMM_AUDIT.read_text(encoding="utf-8"))
    if authority.get("status") != "PASS_A1_UNIFIED_CORE_DFTB3_NUMERICAL_PREFLIGHT":
        raise ValueError("unified QM numerical authority is not PASS")
    expected = {
        "qm_atom_count": EXPECTED_QM_ATOMS,
        "qmcharge": QMCHARGE,
        "electron_count_including_link_h": EXPECTED_ELECTRONS,
        "link_atom_count": EXPECTED_LINK_ATOMS,
    }
    for key, value in expected.items():
        if audit.get(key) != value:
            raise ValueError(f"frozen QM contract changed for {key}")
    scope = str(audit.get("qm_region_scope", ""))
    step1_qm_water_count = 0 if "no QM water" in scope else None
    if step1_qm_water_count != 0:
        raise ValueError("frozen Step1 contract no longer proves no QM water")
    qmmask = audit.get("qmmask")
    if not isinstance(qmmask, str) or len(qmmask.split(",")) != EXPECTED_QM_ATOMS:
        raise ValueError("frozen qmmask is not the 146-atom mask")
    return qmmask, source_gro


def qmmm_block(qmmask: str) -> str:
    return f"""&qmmm
  qmmask='{qmmask}',
  qmcharge={QMCHARGE},
  spin={SPIN},
  qm_theory='DFTB3',
  dftb_telec=200.0,
  qmshake=0,
/
"""


def amber_input(title: str, qmmask: str, steps: int, ncyc: int) -> str:
    return f"""{title}
&cntrl
  imin=1, ntmin=2, maxcyc={steps}, ncyc={ncyc}, dx0=0.005,
  ntb=1, cut=10.0, ntpr=5, ntxo=1,
  ifqnt=1, nmropt=1, ntr=1,
  restraint_wt=1.0,
  restraintmask='{NON_QM_SOLUTE_HEAVY_MASK}',
/
{qmmm_block(qmmask)}&wt type='END' /
DISANG=restraints.RST
DUMPAVE=restraint.dat
"""


def distance_restraint(
    first: int, second: int, target: float, force: float
) -> str:
    return (
        f"&rst iat={first},{second}, "
        f"r1={max(0.1, target - 0.35):.3f}, "
        f"r2={target - 0.05:.3f}, r3={target + 0.05:.3f}, r4=4.500, "
        f"rk2={force:.1f}, rk3={force:.1f}, /\n"
    )


def restraints(targets: dict[str, float], force: float) -> str:
    return "".join(
        (
            distance_restraint(
                THR267_OG1, L2_C12, targets["q_attack"], force
            ),
            distance_restraint(
                THR267_N, TRANSFERRED_HG1, targets["nalpha_hg1"], force
            ),
            distance_restraint(
                TRANSFERRED_HG1, L2_N3, targets["hg1_n3"], force
            ),
            distance_restraint(
                L2_C12, L2_N3, targets["q_cn"], force
            ),
        )
    )


def prepare(
    index: int,
    output: pathlib.Path,
    start_rst7: pathlib.Path,
    github_commit: str,
) -> None:
    window = window_from_index(index)
    qmmask, source_gro = validate_authority(window)
    output.mkdir(parents=True, exist_ok=False)
    if start_rst7.exists():
        raise FileExistsError(start_rst7)
    if not start_rst7.parent.is_dir():
        raise ValueError(f"scratch parent does not exist: {start_rst7.parent}")

    import parmed as pmd

    structure = pmd.load_file(str(PRMTOP), xyz=str(source_gro))
    if len(structure.atoms) != EXPECTED_SYSTEM_ATOMS or structure.box is None:
        raise ValueError("coordinate transplant changed atom count or box")
    expected_atoms = (
        (THR267_N, "THR", "N"),
        (THR267_OG1, "THR", "OG1"),
        (TRANSFERRED_HG1, "THR", "HG1"),
        (L2_C12, "L2", "C12"),
        (L2_O2, "L2", "O2"),
        (L2_N3, "L2", "N3"),
    )
    for index1, resname, name in expected_atoms:
        atom = structure.atoms[index1 - 1]
        if atom.name != name or atom.residue.name != resname:
            raise ValueError(f"atom-order identity changed at {index1}")
    start_geometry = geometry(structure)
    structure.save(str(start_rst7), overwrite=False)

    final_targets = {
        "q_attack": window["attack_target_A"],
        "q_pt2": window["pt2_q_target_A"],
        "nalpha_hg1": window["nh_target_A"],
        "hg1_n3": window["hn3_target_A"],
        "q_cn": window["cn_target_A"],
    }
    guide_targets = {
        "q_attack": (start_geometry["q_attack_A"] + final_targets["q_attack"]) / 2,
        "nalpha_hg1": (
            start_geometry["nalpha_hg1_A"] + final_targets["nalpha_hg1"]
        ) / 2,
        "hg1_n3": (start_geometry["hg1_n3_A"] + final_targets["hg1_n3"]) / 2,
        "q_cn": (start_geometry["q_cn_A"] + final_targets["q_cn"]) / 2,
    }
    guide = output / "guide"
    target = output / "target"
    guide.mkdir()
    target.mkdir()
    (guide / "stage.in").write_text(
        amber_input(
            f"NylC A1 PT2-CN guide {window['tag']}",
            qmmask,
            GUIDE_STEPS,
            50,
        ),
        encoding="utf-8",
    )
    (guide / "restraints.RST").write_text(
        restraints(guide_targets, GUIDE_FORCE), encoding="utf-8"
    )
    (target / "stage.in").write_text(
        amber_input(
            f"NylC A1 PT2-CN target {window['tag']}",
            qmmask,
            TARGET_STEPS,
            150,
        ),
        encoding="utf-8",
    )
    (target / "restraints.RST").write_text(
        restraints(final_targets, TARGET_FORCE), encoding="utf-8"
    )
    manifest = {
        "schema_version": 1,
        "status": "READY_A1_PT2_CN_SCOUT",
        "scientific_status": "NOT_EVALUATED_TS_PMF_BARRIER_MECHANISM",
        "scientific_scope": SCOPE,
        "github_commit": github_commit,
        "window": window,
        "source": {
            "gro": str(source_gro),
            "gro_sha256": sha256(source_gro),
            "system_prmtop": str(PRMTOP),
            "system_prmtop_sha256": sha256(PRMTOP),
            "coordinate_transplant": "ParMed prmtop plus same-order GRO to scratch rst7",
        },
        "qm_contract": {
            "qm_atom_count": EXPECTED_QM_ATOMS,
            "qmcharge": QMCHARGE,
            "electron_count_including_link_h": EXPECTED_ELECTRONS,
            "link_atom_count": EXPECTED_LINK_ATOMS,
            "step1_qm_water_count": 0,
            "qmmask": qmmask,
        },
        "start_geometry": start_geometry,
        "guide_targets_A": guide_targets,
        "final_targets_A": final_targets,
        "protocol": {
            "guide_steps": GUIDE_STEPS,
            "target_steps": TARGET_STEPS,
            "guide_force_kcal_mol_A2": GUIDE_FORCE,
            "target_force_kcal_mol_A2": TARGET_FORCE,
            "attack_angle_restrained": False,
            "pt2_native_difference_restrained": False,
            "pt2_implementation": (
                "paired Nalpha-HG1 and HG1-N3 distance restraints; "
                "observed q_PT2 is audited as their difference"
            ),
            "windows_are_independent": True,
        },
    }
    write_json(output / "SCOUT_MANIFEST.json", manifest)


def inspect_stage(prmtop: pathlib.Path, directory: pathlib.Path) -> dict[str, Any]:
    output = directory / "stage.out"
    restart = directory / "stage.rst7"
    text = (
        output.read_text(encoding="utf-8", errors="replace")
        if output.is_file()
        else ""
    )
    hard = {
        name: len(re.findall(pattern, text, re.I))
        for name, pattern in HARD_PATTERNS.items()
    }
    complete = "FINAL RESULTS" in text and bool(re.search(r"Run\s+done", text))
    scc = len(re.findall(r"Convergence could not be achieved", text, re.I))
    vlimit = len(re.findall(r"vlimit\s+exceeded", text, re.I))
    bond_overflow = len(re.findall(r"BOND\s*=\s*\*+", text, re.I))
    records = re.findall(
        r"^\s*(\d+)\s+(-?\d+\.\d+E[+-]\d+)\s+"
        r"(\d+\.\d+E[+-]\d+)\s+(\d+\.\d+E[+-]\d+)",
        text,
        re.M,
    )
    observed: dict[str, float] = {}
    if restart.is_file() and restart.stat().st_size:
        import parmed as pmd

        observed = geometry(pmd.load_file(str(prmtop), xyz=str(restart)))
    finite = observed and all(
        isinstance(value, (int, float)) and math.isfinite(value)
        for value in observed.values()
    )
    technical = bool(
        complete
        and scc == 0
        and vlimit == 0
        and bond_overflow == 0
        and sum(hard.values()) == 0
        and finite
    )
    return {
        "complete": complete,
        "scc_warnings": scc,
        "vlimit_warnings": vlimit,
        "bond_overflow": bond_overflow,
        "hard_error_hits": hard,
        "geometry": observed,
        "last_minimization_record": (
            {
                "nstep": int(records[-1][0]),
                "energy_kcal_mol": float(records[-1][1]),
                "rms": float(records[-1][2]),
                "gmax": float(records[-1][3]),
            }
            if records
            else None
        ),
        "technical_pass": technical,
    }


def audit(output: pathlib.Path) -> None:
    manifest_path = output / "SCOUT_MANIFEST.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "READY_A1_PT2_CN_SCOUT":
        raise ValueError("scout manifest is not READY")
    window = manifest["window"]
    stages = {
        name: inspect_stage(PRMTOP, output / name)
        for name in ("guide", "target")
    }
    technical = all(stage["technical_pass"] for stage in stages.values())
    final = stages["target"]["geometry"] if technical else {}
    target_record = stages["target"]["last_minimization_record"]
    residuals: dict[str, float] = {}
    if technical:
        residuals = {
            "q_attack_A": abs(final["q_attack_A"] - window["attack_target_A"]),
            "nalpha_hg1_A": abs(final["nalpha_hg1_A"] - window["nh_target_A"]),
            "hg1_n3_A": abs(final["hg1_n3_A"] - window["hn3_target_A"]),
            "pt2_q_A": abs(final["pt2_q_A"] - window["pt2_q_target_A"]),
            "q_cn_A": abs(final["q_cn_A"] - window["cn_target_A"]),
        }
    target_response = bool(
        technical
        and residuals["q_attack_A"] <= TARGET_TOLERANCE_A
        and residuals["nalpha_hg1_A"] <= TARGET_TOLERANCE_A
        and residuals["hg1_n3_A"] <= TARGET_TOLERANCE_A
        and residuals["pt2_q_A"] <= PT2_Q_TOLERANCE_A
        and residuals["q_cn_A"] <= TARGET_TOLERANCE_A
    )
    angle_pass = bool(
        technical
        and ATTACK_ANGLE_RANGE_DEG[0]
        <= final["attack_angle_deg"]
        <= ATTACK_ANGLE_RANGE_DEG[1]
    )
    rms_pass = bool(
        target_record
        and target_record["rms"] <= RMS_GRADIENT_MAX
    )
    candidate = bool(target_response and angle_pass and rms_pass)
    candidate_gate = (
        "PASS_A1_PT2_CN_SCOUT_WINDOW_CANDIDATE"
        if candidate
        else "NOT_SELECTED_A1_PT2_CN_SCOUT_WINDOW"
    )
    result = {
        "schema_version": 1,
        "status": PASS_STATUS if technical else FAIL_STATUS,
        "scientific_status": "NOT_EVALUATED_TS_PMF_BARRIER_MECHANISM",
        "scientific_scope": SCOPE,
        "window": window,
        "candidate_gate": candidate_gate,
        "selection_policy": (
            "Do not compare total energies across differently restrained windows; "
            "select zero to two candidates per seed after all eight terminate."
        ),
        "gates": {
            "technical": technical,
            "target_response": target_response,
            "attack_angle_90_130_deg": angle_pass,
            "target_rms_gradient_le_0p5": rms_pass,
        },
        "target_residuals_A": residuals,
        "stages": stages,
    }
    write_json(output / "RESULT.json", result)
    if technical:
        write_json(output / "PASS.json", result)
    else:
        write_json(output / "NOT_EVALUATED.json", result)
        raise RuntimeError("scout window did not technically PASS")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("prepare", "audit"), required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--start-rst7", type=pathlib.Path)
    parser.add_argument("--github-commit", default="unknown")
    args = parser.parse_args()
    if args.mode == "prepare":
        if args.start_rst7 is None:
            parser.error("--start-rst7 is required for prepare")
        prepare(
            args.index,
            args.output.resolve(),
            args.start_rst7.resolve(),
            args.github_commit,
        )
    else:
        audit(args.output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
