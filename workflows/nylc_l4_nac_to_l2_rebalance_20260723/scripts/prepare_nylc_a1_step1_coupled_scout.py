#!/usr/bin/env python3
"""Prepare one A1 Step1 coupled attack/carbonyl scout and local release."""
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
THR267_OG1 = 8960
L2_C12 = 10287
L2_O2 = 10288
L2_N3 = 10289
ATTACK_TARGET_A = 1.45
ATTACK_FORCE = 200.0
CARBONYL_TARGETS_A = (1.30, 1.35, 1.40)
CARBONYL_FORCE = 500.0
QM_RANGES = "@7160-7174,7756-7771,8235-8242,8949-8963,9567-9573,9587-9592,10273-10351"
NON_QM_SOLUTE_HEAVY_MASK = f"!(:SOL,NA,CL)&!@H=&!({QM_RANGES})"


def sha256(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def qmmm_block(qmmask):
    return f"""&qmmm
  qmmask='{qmmask}',
  qmcharge=0,
  spin=1,
  qm_theory='DFTB3',
  dftb_telec=200.0,
  qmshake=0,
/
"""


def coupled_input(qmmask, carbonyl_target):
    return f"""NylC A1 coupled attack-carbonyl scout C-O={carbonyl_target:.2f} A
&cntrl
  imin=1, ntmin=2, maxcyc=300, ncyc=300, dx0=0.005,
  ntb=1, cut=10.0, ntpr=5, ntxo=1,
  ifqnt=1, nmropt=1, ntr=1,
  restraint_wt=1.0,
  restraintmask='{NON_QM_SOLUTE_HEAVY_MASK}',
/
{qmmm_block(qmmask)}&wt type='END' /
DISANG=coupled.RST
DUMPAVE=coupled_restraint.dat
"""


def local_release_input(qmmask):
    return f"""NylC A1 coupled-scout local release
&cntrl
  imin=1, ntmin=1, maxcyc=500, ncyc=100, dx0=0.005,
  ntb=1, cut=10.0, ntpr=5, ntxo=1,
  ifqnt=1, ntr=1,
  restraint_wt=1.0,
  restraintmask='{NON_QM_SOLUTE_HEAVY_MASK}',
/
{qmmm_block(qmmask)}"""


def restraint_text(carbonyl_target):
    return (
        f"&rst iat={THR267_OG1},{L2_C12}, "
        f"r1={ATTACK_TARGET_A-0.25:.3f}, r2={ATTACK_TARGET_A-0.025:.3f}, "
        f"r3={ATTACK_TARGET_A+0.025:.3f}, r4=4.500, "
        f"rk2={ATTACK_FORCE:.1f}, rk3={ATTACK_FORCE:.1f}, /\n"
        f"&rst iat={L2_O2},{L2_C12},{THR267_OG1}, "
        "r1=85.0, r2=95.0, r3=115.0, r4=125.0, "
        "rk2=20.0, rk3=20.0, /\n"
        f"&rst iat={L2_C12},{L2_O2}, "
        f"r1={carbonyl_target-0.20:.3f}, r2={carbonyl_target-0.015:.3f}, "
        f"r3={carbonyl_target+0.015:.3f}, r4={carbonyl_target+0.20:.3f}, "
        f"rk2={CARBONYL_FORCE:.1f}, rk3={CARBONYL_FORCE:.1f}, /\n"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--carbonyl-target", type=float, choices=CARBONYL_TARGETS_A, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)

    source_pass = QMMM_SOURCE / "PASS.json"
    source_audit = QMMM_SOURCE / "prepared/qmmm_preflight_audit.json"
    source_prmtop = QMMM_SOURCE / "prepared/system.prmtop"
    seed_pass = SEED_SOURCE / "PASS.json"
    for path in (source_pass, source_audit, source_prmtop, seed_pass, START_RESTART):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"missing immutable coupled-scout source: {path}")

    for path, expected in {
        seed_pass: EXPECTED_SEED_PASS_SHA256,
        START_RESTART: EXPECTED_RESTART_SHA256,
        source_prmtop: EXPECTED_PRMTOP_SHA256,
    }.items():
        observed = sha256(path)
        if observed != expected:
            raise ValueError(f"source hash mismatch: {path}: {observed}")

    authority = json.loads(source_pass.read_text(encoding="utf-8"))
    audit = json.loads(source_audit.read_text(encoding="utf-8"))
    seed = json.loads(seed_pass.read_text(encoding="utf-8"))
    if authority.get("status") != "PASS_A1_UNIFIED_CORE_DFTB3_NUMERICAL_PREFLIGHT":
        raise ValueError("unified-core numerical source did not PASS")
    if seed.get("status") != "PASS_TECHNICAL_A1_UNIFIED_STEP1_QATTACK_EXTENSION":
        raise ValueError("attack seed did not technically PASS")
    if seed.get("attack_bracket_seed_gate") != "PASS_CONSTRAINED_ATTACK_BRACKET_SEED":
        raise ValueError("source is not the frozen constrained attack seed")

    for key, value in {
        "qm_atom_count": EXPECTED_QM_ATOMS,
        "qmcharge": QMCHARGE,
        "spin": SPIN,
        "electron_count_including_link_h": EXPECTED_ELECTRONS,
        "link_atom_count": 6,
        "bond_count_gt_3A": 0,
    }.items():
        if audit.get(key) != value:
            raise ValueError(f"unified-core contract mismatch for {key}")
    qmmask = audit.get("qmmask")
    if not isinstance(qmmask, str) or len(qmmask.split(",")) != EXPECTED_QM_ATOMS:
        raise ValueError("source qmmask is not the frozen 146-atom core")

    coupled = output / "coupled"
    release = output / "local_release"
    coupled.mkdir()
    release.mkdir()
    (coupled / "stage.in").write_text(
        coupled_input(qmmask, args.carbonyl_target), encoding="utf-8"
    )
    (coupled / "coupled.RST").write_text(
        restraint_text(args.carbonyl_target), encoding="utf-8"
    )
    (release / "stage.in").write_text(local_release_input(qmmask), encoding="utf-8")

    manifest = {
        "schema_version": 1,
        "status": "READY_A1_STEP1_COUPLED_SCOUT",
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
            "attack_target_A": ATTACK_TARGET_A,
            "attack_force_kcal_mol_A2": ATTACK_FORCE,
            "carbonyl_target_A": args.carbonyl_target,
            "carbonyl_force_kcal_mol_A2": CARBONYL_FORCE,
            "angle_flat_bottom_deg": [95.0, 115.0],
            "stages": ["coupled", "local_release"],
            "coupled_steps": 300,
            "release_steps": 500,
            "release_reactive_restraints": False,
            "non_qm_solute_heavy_position_restraint_kcal_mol_A2": 1.0,
        },
        "interpretation": (
            "Exploratory coupled attack/carbonyl scout followed by local release. "
            "Not a TI, TS, committor, PMF, barrier or mechanism proof."
        ),
    }
    (output / "COUPLED_SCOUT_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": manifest["status"],
        "attack_target_A": ATTACK_TARGET_A,
        "carbonyl_target_A": args.carbonyl_target,
    }))


if __name__ == "__main__":
    main()
