#!/usr/bin/env python3
"""Prepare one pinned A1 Step1 q_attack extension window from the k=200 endpoint."""
import argparse
import hashlib
import json
import pathlib

TASK_ROOT = pathlib.Path("/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723")
QMMM_SOURCE = TASK_ROOT / "a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285"
RECOVERY_SOURCE = TASK_ROOT / "a1_activated_nac_20260726/qmmm/a1_unified_step1_qattack_recovery/attempt_62016429_2_k200"
START_RESTART = RECOVERY_SOURCE / "scan/q05_1p65A/window.rst7"
EXPECTED_RECOVERY_PASS_SHA256 = "80ae4121248f90d8488125ce69d51aa31198d19970f63bbe04fe65ee57503081"
EXPECTED_RESTART_SHA256 = "6406421d86ae35330cc6a480e5f5b234f9ff8591e97d55743e354095b42a5e74"
EXPECTED_PRMTOP_SHA256 = "a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0"
EXPECTED_QM_ATOMS = 146
EXPECTED_ELECTRONS = 510
QMCHARGE = 0
SPIN = 1
THR267_OG1 = 8960
L2_C12 = 10287
L2_O2 = 10288
L2_N3 = 10289
TARGET_A = 1.45
FORCE_CONSTANT = 200.0
QM_RANGES = "@7160-7174,7756-7771,8235-8242,8949-8963,9567-9573,9587-9592,10273-10351"
NON_QM_SOLUTE_HEAVY_MASK = f"!(:SOL,NA,CL)&!@H=&!({QM_RANGES})"


def sha256(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def amber_input(qmmask):
    return f"""NylC A1 Step1 q_attack extension k=200 target=1.45 A
&cntrl
  imin=1, ntmin=2, maxcyc=150, ncyc=150, dx0=0.005,
  ntb=1, cut=10.0, ntpr=1, ntxo=1,
  ifqnt=1, nmropt=1, ntr=1,
  restraint_wt=1.0,
  restraintmask='{NON_QM_SOLUTE_HEAVY_MASK}',
/
&qmmm
  qmmask='{qmmask}',
  qmcharge=0,
  spin=1,
  qm_theory='DFTB3',
  dftb_telec=200.0,
  qmshake=0,
/
&wt type='END' /
DISANG=window.RST
DUMPAVE=window_restraint.dat
"""


def restraint_text():
    return (
        f"&rst iat={THR267_OG1},{L2_C12}, "
        f"r1={TARGET_A-0.25:.3f}, r2={TARGET_A-0.025:.3f}, "
        f"r3={TARGET_A+0.025:.3f}, r4=4.500, "
        f"rk2={FORCE_CONSTANT:.1f}, rk3={FORCE_CONSTANT:.1f}, /\n"
        f"&rst iat={L2_O2},{L2_C12},{THR267_OG1}, "
        "r1=85.0, r2=95.0, r3=115.0, r4=125.0, "
        "rk2=20.0, rk3=20.0, /\n"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)

    source_pass = QMMM_SOURCE / "PASS.json"
    source_audit = QMMM_SOURCE / "prepared/qmmm_preflight_audit.json"
    source_prmtop = QMMM_SOURCE / "prepared/system.prmtop"
    recovery_pass = RECOVERY_SOURCE / "PASS.json"
    for path in (source_pass, source_audit, source_prmtop, recovery_pass, START_RESTART):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"missing immutable extension source: {path}")

    expected_hashes = {
        recovery_pass: EXPECTED_RECOVERY_PASS_SHA256,
        START_RESTART: EXPECTED_RESTART_SHA256,
        source_prmtop: EXPECTED_PRMTOP_SHA256,
    }
    for path, expected in expected_hashes.items():
        observed = sha256(path)
        if observed != expected:
            raise ValueError(f"source hash mismatch: {path}: {observed}")

    authority = json.loads(source_pass.read_text(encoding="utf-8"))
    audit = json.loads(source_audit.read_text(encoding="utf-8"))
    recovery = json.loads(recovery_pass.read_text(encoding="utf-8"))
    if authority.get("status") != "PASS_A1_UNIFIED_CORE_DFTB3_NUMERICAL_PREFLIGHT":
        raise ValueError("unified-core numerical source did not PASS")
    if recovery.get("status") != "PASS_TECHNICAL_A1_UNIFIED_STEP1_QATTACK_RECOVERY":
        raise ValueError("k=200 recovery source did not technically PASS")
    if recovery.get("attack_bracket_seed_gate") != "NOT_EVALUATED_ATTACK_BRACKET_NOT_REACHED":
        raise ValueError("extension source is not the expected near-bracket result")
    if recovery.get("force_constant_kcal_mol_A2") != FORCE_CONSTANT:
        raise ValueError("extension source is not k=200")

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

    name = "q06_1p45A"
    directory = output / name
    directory.mkdir()
    (directory / "window.in").write_text(amber_input(qmmask), encoding="utf-8")
    (directory / "window.RST").write_text(restraint_text(), encoding="utf-8")
    manifest = {
        "schema_version": 1,
        "status": "READY_A1_UNIFIED_STEP1_QATTACK_EXTENSION",
        "scientific_status": "NOT_EVALUATED_TS_PMF_BARRIER",
        "source_recovery": str(RECOVERY_SOURCE),
        "source_recovery_status": recovery["status"],
        "source_attack_gate": recovery["attack_bracket_seed_gate"],
        "source_restart": str(START_RESTART),
        "source_sha256": {
            "recovery_PASS": sha256(recovery_pass),
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
            "coordinate": "OG1-C12 distance",
            "target_A": TARGET_A,
            "distance_force_kcal_mol_A2": FORCE_CONSTANT,
            "distance_outer_r4_A": 4.5,
            "angle_flat_bottom_deg": [95.0, 115.0],
            "steps": 150,
            "minimizer": "steepest_descent_only",
            "non_qm_solute_heavy_position_restraint_kcal_mol_A2": 1.0,
            "propagation": "single extension from job-62016429 task-2 q05 restart",
        },
        "reactive_atoms_index1": {
            "thr267_og1": THR267_OG1,
            "l2_c12": L2_C12,
            "l2_o2": L2_O2,
            "l2_n3": L2_N3,
            "nalpha_h": [8950, 8951, 8961],
            "thr267_nalpha": 8949,
        },
        "windows": [{"index": 0, "name": name, "target_A": TARGET_A}],
        "interpretation": "Minimal constrained attack-bracket extension; not a TS, PMF, barrier or mechanism.",
    }
    (output / "SCAN_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"status": manifest["status"], "target_A": TARGET_A, "windows": 1}))


if __name__ == "__main__":
    main()
