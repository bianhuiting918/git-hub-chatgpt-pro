#!/usr/bin/env python3
"""Prepare sequential local and full release tests from the pinned A1 attack seed."""
import argparse
import hashlib
import json
import pathlib

TASK_ROOT = pathlib.Path("/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723")
QMMM_SOURCE = TASK_ROOT / "a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285"
SEED_SOURCE = TASK_ROOT / "a1_activated_nac_20260726/qmmm/a1_unified_step1_qattack_extension/attempt_62021985"
START_RESTART = SEED_SOURCE / "scan/q06_1p45A/window.rst7"
EXPECTED_SEED_PASS_SHA256 = "3083663f502f522974f7eac56e1b3a4c7015933c7d949562364d359006c25224"
EXPECTED_RESTART_SHA256 = "548dca57cb623e69d5fefaa8e270267a5bd599bd5995113219cdf349c346e56c"
EXPECTED_PRMTOP_SHA256 = "a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0"
EXPECTED_QM_ATOMS = 146
EXPECTED_ELECTRONS = 510
QMCHARGE = 0
SPIN = 1
QM_RANGES = "@7160-7174,7756-7771,8235-8242,8949-8963,9567-9573,9587-9592,10273-10351"
NON_QM_SOLUTE_HEAVY_MASK = f"!(:SOL,NA,CL)&!@H=&!({QM_RANGES})"


def sha256(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def amber_input(title, qmmask, restrained):
    release_controls = (
        "  ntr=1, restraint_wt=1.0,\n"
        f"  restraintmask='{NON_QM_SOLUTE_HEAVY_MASK}',\n"
        if restrained
        else "  ntr=0,\n"
    )
    return f"""NylC A1 Step1 {title}
&cntrl
  imin=1, ntmin=1, maxcyc=500, ncyc=100, dx0=0.005,
  ntb=1, cut=10.0, ntpr=5, ntxo=1,
  ifqnt=1,
{release_controls}/
&qmmm
  qmmask='{qmmask}',
  qmcharge=0,
  spin=1,
  qm_theory='DFTB3',
  dftb_telec=200.0,
  qmshake=0,
/
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)

    source_pass = QMMM_SOURCE / "PASS.json"
    source_audit = QMMM_SOURCE / "prepared/qmmm_preflight_audit.json"
    source_prmtop = QMMM_SOURCE / "prepared/system.prmtop"
    seed_pass = SEED_SOURCE / "PASS.json"
    for path in (source_pass, source_audit, source_prmtop, seed_pass, START_RESTART):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"missing immutable release source: {path}")

    expected_hashes = {
        seed_pass: EXPECTED_SEED_PASS_SHA256,
        START_RESTART: EXPECTED_RESTART_SHA256,
        source_prmtop: EXPECTED_PRMTOP_SHA256,
    }
    for path, expected in expected_hashes.items():
        observed = sha256(path)
        if observed != expected:
            raise ValueError(f"source hash mismatch: {path}: {observed}")

    authority = json.loads(source_pass.read_text(encoding="utf-8"))
    audit = json.loads(source_audit.read_text(encoding="utf-8"))
    seed = json.loads(seed_pass.read_text(encoding="utf-8"))
    if authority.get("status") != "PASS_A1_UNIFIED_CORE_DFTB3_NUMERICAL_PREFLIGHT":
        raise ValueError("unified-core numerical source did not PASS")
    if seed.get("status") != "PASS_TECHNICAL_A1_UNIFIED_STEP1_QATTACK_EXTENSION":
        raise ValueError("attack extension did not technically PASS")
    if seed.get("attack_bracket_seed_gate") != "PASS_CONSTRAINED_ATTACK_BRACKET_SEED":
        raise ValueError("attack extension did not yield the required constrained seed")

    expected = {
        "qm_atom_count": EXPECTED_QM_ATOMS,
        "qmcharge": QMCHARGE,
        "spin": SPIN,
        "electron_count_including_link_h": EXPECTED_ELECTRONS,
        "link_atom_count": 6,
        "bond_count_gt_3A": 0,
    }
    for key, value in expected.items():
        if audit.get(key) != value:
            raise ValueError(f"unified-core contract mismatch for {key}")
    qmmask = audit.get("qmmask")
    if not isinstance(qmmask, str) or len(qmmask.split(",")) != EXPECTED_QM_ATOMS:
        raise ValueError("source qmmask is not the frozen 146-atom core")

    stages = [
        ("local_release", True),
        ("full_release", False),
    ]
    manifest_stages = []
    for index, (name, restrained) in enumerate(stages):
        directory = output / name
        directory.mkdir()
        (directory / "stage.in").write_text(
            amber_input(name.replace("_", " "), qmmask, restrained),
            encoding="utf-8",
        )
        manifest_stages.append(
            {
                "index": index,
                "name": name,
                "non_qm_solute_heavy_position_restraint_kcal_mol_A2": 1.0 if restrained else 0.0,
                "reactive_distance_restraint": False,
                "reactive_angle_restraint": False,
            }
        )

    manifest = {
        "schema_version": 1,
        "status": "READY_A1_STEP1_RELEASE_TEST",
        "scientific_status": "NOT_EVALUATED_TS_PMF_BARRIER",
        "source_seed": str(SEED_SOURCE),
        "source_seed_status": seed["status"],
        "source_seed_gate": seed["attack_bracket_seed_gate"],
        "source_restart": str(START_RESTART),
        "source_sha256": {
            "seed_PASS": sha256(seed_pass),
            "start_restart": sha256(START_RESTART),
            "system_prmtop": sha256(source_prmtop),
        },
        "qm_contract": {
            "explicit_qm_atoms": EXPECTED_QM_ATOMS,
            "qmcharge": QMCHARGE,
            "spin": SPIN,
            "electrons_including_links": EXPECTED_ELECTRONS,
            "link_h": 6,
            "step1_qm_water_count": 0,
        },
        "protocol": {
            "stages": ["local_release", "full_release"],
            "steps_per_stage": 500,
            "steepest_descent_steps": 100,
            "then_conjugate_gradient": True,
            "propagation": "sequential from constrained seed; no reactive distance or angle restraints",
        },
        "stages": manifest_stages,
        "interpretation": (
            "Release test for local attack-basin retention. It is not a TS, TI proof, "
            "committor, PMF, barrier or mechanism."
        ),
    }
    (output / "RELEASE_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"status": manifest["status"], "stages": manifest["protocol"]["stages"]}))


if __name__ == "__main__":
    main()
