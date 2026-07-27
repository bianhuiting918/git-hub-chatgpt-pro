#!/usr/bin/env python3
"""Prepare one force-constant replica of the unified A1 Step1 q_attack recovery scan."""
import argparse
import hashlib
import json
import pathlib

TASK_ROOT = pathlib.Path("/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723")
FAILED_BRACKET = TASK_ROOT / "a1_activated_nac_20260726/qmmm/a1_unified_step1_qattack/attempt_62012919"
QMMM_SOURCE = TASK_ROOT / "a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285"
START_RESTART = FAILED_BRACKET / "scan/q03_2p65A/window.rst7"
EXPECTED_SOURCE_STATUS = "PASS_A1_UNIFIED_CORE_DFTB3_NUMERICAL_PREFLIGHT"
EXPECTED_SCAN_STATUS = "PASS_TECHNICAL_A1_UNIFIED_STEP1_QATTACK_SCAN"
THR267_OG1 = 8960
L2_C12 = 10287
L2_O2 = 10288
L2_N3 = 10289
EXPECTED_QM_ATOMS = 146
QMCHARGE = 0
SPIN = 1
EXPECTED_ELECTRONS = 510
TARGETS_A = (2.65, 2.45, 2.25, 2.05, 1.85, 1.65)
FORCE_CONSTANTS = (50.0, 100.0, 200.0)
QM_RANGES = "@7160-7174,7756-7771,8235-8242,8949-8963,9567-9573,9587-9592,10273-10351"
NON_QM_SOLUTE_HEAVY_MASK = f"!(:SOL,NA,CL)&!@H=&!({QM_RANGES})"


def sha256(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def amber_input(title, qmmask):
    return f"""{title}
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


def restraint_text(target, force_constant):
    return (
        f"&rst iat={THR267_OG1},{L2_C12}, "
        f"r1={target-0.25:.3f}, r2={target-0.025:.3f}, "
        f"r3={target+0.025:.3f}, r4=4.500, "
        f"rk2={force_constant:.1f}, rk3={force_constant:.1f}, /\n"
        f"&rst iat={L2_O2},{L2_C12},{THR267_OG1}, "
        "r1=85.0, r2=95.0, r3=115.0, r4=125.0, "
        "rk2=20.0, rk3=20.0, /\n"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--force-constant", type=float, choices=FORCE_CONSTANTS, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)

    source_pass = QMMM_SOURCE / "PASS.json"
    source_audit = QMMM_SOURCE / "prepared/qmmm_preflight_audit.json"
    source_prmtop = QMMM_SOURCE / "prepared/system.prmtop"
    prior_pass = FAILED_BRACKET / "PASS.json"
    for path in (source_pass, source_audit, source_prmtop, prior_pass, START_RESTART):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"missing immutable recovery source: {path}")

    authority = json.loads(source_pass.read_text(encoding="utf-8"))
    audit = json.loads(source_audit.read_text(encoding="utf-8"))
    prior = json.loads(prior_pass.read_text(encoding="utf-8"))
    if authority.get("status") != EXPECTED_SOURCE_STATUS:
        raise ValueError("unified-core numerical source did not PASS")
    if prior.get("status") != EXPECTED_SCAN_STATUS:
        raise ValueError("job 62012919 is not the expected technical PASS")
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
            raise ValueError(f"unified-core source contract mismatch for {key}")
    qmmask = audit.get("qmmask")
    if not isinstance(qmmask, str) or len(qmmask.split(",")) != EXPECTED_QM_ATOMS:
        raise ValueError("source qmmask is not the frozen 146-atom core")

    windows = []
    for index, target in enumerate(TARGETS_A):
        name = f"q{index:02d}_{target:.2f}A".replace(".", "p")
        directory = output / name
        directory.mkdir()
        (directory / "window.in").write_text(
            amber_input(
                f"NylC A1 Step1 q_attack recovery k={args.force_constant:.1f} target={target:.2f} A",
                qmmask,
            ),
            encoding="utf-8",
        )
        (directory / "window.RST").write_text(
            restraint_text(target, args.force_constant), encoding="utf-8"
        )
        windows.append({"index": index, "name": name, "target_A": target})

    manifest = {
        "schema_version": 1,
        "status": "READY_A1_UNIFIED_STEP1_QATTACK_RECOVERY",
        "scientific_status": "NOT_EVALUATED_TS_PMF_BARRIER",
        "source_scan": str(FAILED_BRACKET),
        "source_scan_status": prior["status"],
        "source_restart": str(START_RESTART),
        "source_sha256": {
            "start_restart": sha256(START_RESTART),
            "system_prmtop": sha256(source_prmtop),
            "unified_preflight_PASS": sha256(source_pass),
            "prior_scan_PASS": sha256(prior_pass),
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
            "targets_A": list(TARGETS_A),
            "distance_force_kcal_mol_A2": args.force_constant,
            "distance_outer_r4_A": 4.5,
            "angle_flat_bottom_deg": [95.0, 115.0],
            "steps_per_window": 150,
            "minimizer": "steepest_descent_only",
            "non_qm_solute_heavy_position_restraint_kcal_mol_A2": 1.0,
            "propagation": "sequential from common job-62012919 q03 restart",
        },
        "reactive_atoms_index1": {
            "thr267_og1": THR267_OG1,
            "l2_c12": L2_C12,
            "l2_o2": L2_O2,
            "l2_n3": L2_N3,
            "nalpha_h": [8950, 8951, 8961],
            "thr267_nalpha": 8949,
        },
        "windows": windows,
        "interpretation": "Corrected constrained attack scouting; not a TS, PMF, barrier or mechanism.",
    }
    (output / "SCAN_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"status": manifest["status"], "force": args.force_constant, "windows": len(windows)}))


if __name__ == "__main__":
    main()
