#!/usr/bin/env python3
"""Prepare a gated NylC PA66-L2 DFTB3 numerical preflight from a free NAC frame."""
import argparse
import hashlib
import json
import os
import pathlib
import subprocess

import parmed as pmd

import extract_true_thr267_recapture_sources as geom


PROTEIN_ATOMS = 10272
THR267_OG1 = 8961
L2_FIRST = 10273
L2_LAST = 10351
L2_REACTIVE_C = 10287
L2_REACTIVE_O = 10288
QMCHARGE = 1
SPIN = 1
LINK_ATOMIC_NUMBER = 1
DFTB_TELEC_K = 200.0
ATOMIC_NUMBERS = {"H": 1, "C": 6, "N": 7, "O": 8}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(command, stdin=None):
    result = subprocess.run(command, input=stdin, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result


def element(atom):
    if atom.atomic_number in ATOMIC_NUMBERS.values():
        return next(name for name, number in ATOMIC_NUMBERS.items() if atom.atomic_number == number)
    letters = "".join(character for character in atom.name if character.isalpha()).upper()
    if not letters or letters[0] not in ATOMIC_NUMBERS:
        raise ValueError(f"unsupported QM element at atom {atom.idx + 1}: {atom.name}")
    return letters[0]


def qmmm_input(title, maxcyc, qmmask):
    return f"""{title}
&cntrl
  imin=1, maxcyc={maxcyc}, ncyc={max(1, maxcyc // 2)},
  ntb=1, cut=10.0, ntpr=1,
  ifqnt=1,
/
&qmmm
  qmmask='{qmmask}',
  qmcharge={QMCHARGE},
  spin=1,
  qm_theory='DFTB3',
  dftb_telec={DFTB_TELEC_K},
  qmshake=0,
/
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-root", type=pathlib.Path, required=True)
    parser.add_argument("--post-root", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--gmx", required=True)
    args = parser.parse_args()

    task = args.task_root.resolve()
    post = args.post_root.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    scientific = json.loads((post / "l2_free_1ns_audit.json").read_text())
    series = json.loads((post / "nac_energy_series.json").read_text())
    selected = series.get("lowest_potential_nac_frame")
    if scientific.get("scientific_status") != "PASS_UNRESTRAINED_L2_NAC_PRESENT" or not selected:
        result = {
            "schema_version": 1,
            "candidate_id": "nylc_C18_trueT267_freeGS",
            "status": "NOT_EVALUATED_NO_FREE_L2_NAC",
            "source_scientific_status": scientific.get("scientific_status"),
            "nac_frame_count": series.get("nac_frame_count", 0),
            "interpretation": "No fully unrestrained L2 NAC is available for Step1 DFTB3 preflight.",
        }
        (output / "NOT_EVALUATED.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        print(json.dumps(result, sort_keys=True))
        return 0

    free = task / "candidates/nylc_C18_trueT267_freeGS/runs/rebalance_after_em_double/npt300free"
    build = task / "candidates/nylc_C18_trueT267_freeGS/build"
    temporary = output / "source.tmp.gro"
    promoted = output / "selected_free_l2_nac.gro"
    run(
        [
            args.gmx, "trjconv", "-s", str(free / "run.tpr"), "-f", str(free / "run.xtc"),
            "-dump", f"{selected['time_ps']:.6f}", "-pbc", "mol", "-o", str(temporary),
        ],
        stdin="System\n",
    )

    atoms, box, _ = geom.read_gro(temporary)
    if len(atoms) != 133589:
        raise ValueError(f"unexpected atom count: {len(atoms)}")
    thr = geom.validate_identity(atoms, THR267_OG1, (267, "THR", "OG1"))
    carbon = atoms[L2_REACTIVE_C]
    oxygen = atoms[L2_REACTIVE_O]
    c_to_thr = geom.minimum_image(geom.subtract(thr["xyz_nm"], carbon["xyz_nm"]), box)
    c_to_o = geom.minimum_image(geom.subtract(oxygen["xyz_nm"], carbon["xyz_nm"]), box)
    distance = geom.norm(c_to_thr)
    angle = geom.angle_deg(c_to_thr, c_to_o)
    if abs(distance - selected["distance_nm"]) > 0.002 or abs(angle - selected["angle_deg"]) > 0.3:
        raise ValueError(f"extracted NAC geometry mismatch: {distance:.6f} nm/{angle:.3f} deg vs {selected}")
    if not (distance <= 0.35 and 95.0 <= angle <= 115.0):
        raise ValueError("selected frame is not an audited true-Thr267 NAC")
    os.replace(temporary, promoted)

    previous = pathlib.Path.cwd()
    os.chdir(build)
    try:
        structure = pmd.load_file(str(build / "topol.top"), xyz=str(promoted))
    finally:
        os.chdir(previous)
    if len(structure.atoms) != 133589:
        raise ValueError(f"ParmEd atom count mismatch: {len(structure.atoms)}")
    active = structure.atoms[THR267_OG1 - 1]
    if active.name != "OG1" or active.residue.name != "THR":
        raise ValueError(f"active identity mismatch: {active.residue.name}:{active.name}")
    ligand = list(structure.atoms[L2_FIRST - 1 : L2_LAST])
    if len(ligand) != 79 or {atom.residue.name for atom in ligand} != {"L2"}:
        raise ValueError("complete 79-atom L2 ligand ordering changed")
    if structure.atoms[L2_REACTIVE_C - 1] not in ligand or structure.atoms[L2_REACTIVE_O - 1] not in ligand:
        raise ValueError("reactive L2 C/O is outside the complete ligand QM region")

    active_topology_charge = sum(atom.charge for atom in active.residue.atoms)
    ligand_topology_charge = sum(atom.charge for atom in ligand)
    if abs(active_topology_charge - 1.0) > 1.0e-4:
        raise ValueError(f"processed N-terminal Thr267 charge is not +1: {active_topology_charge}")
    if abs(ligand_topology_charge) > 1.0e-4:
        raise ValueError(f"audited PA66-L2 charge is not neutral: {ligand_topology_charge}")
    qm_atoms = list(active.residue.atoms) + ligand
    qm_indices = {atom.idx for atom in qm_atoms}
    qm_elements = [element(atom) for atom in qm_atoms]
    electrons = sum(ATOMIC_NUMBERS[item] for item in qm_elements) + LINK_ATOMIC_NUMBER - QMCHARGE
    if electrons % 2:
        raise ValueError(f"odd electron count {electrons} is incompatible with spin=1")
    boundary = []
    for bond in structure.bonds:
        left, right = bond.atom1, bond.atom2
        if (left.idx in qm_indices) != (right.idx in qm_indices):
            boundary.append({
                "qm_atom": (left if left.idx in qm_indices else right).idx + 1,
                "mm_atom": (right if left.idx in qm_indices else left).idx + 1,
            })
    if len(boundary) != 1:
        raise ValueError(f"expected one QM/MM boundary bond, observed {boundary}")

    structure.save(str(output / "system.prmtop"), overwrite=True)
    structure.save(str(output / "selected_free_l2_nac.rst7"), overwrite=True)
    qmmask = "@" + ",".join(str(atom.idx + 1) for atom in qm_atoms)
    (output / "01_qmmm_one_step.in").write_text(
        qmmm_input("NylC PA66-L2 Step1 DFTB3 one-step numerical preflight", maxcyc=1, qmmask=qmmask)
    )
    (output / "02_qmmm_20_step.in").write_text(
        qmmm_input("NylC PA66-L2 Step1 DFTB3 20-step segmented minimization", maxcyc=20, qmmask=qmmask)
    )
    audit = {
        "schema_version": 1,
        "candidate_id": "nylc_C18_trueT267_freeGS",
        "status": "READY_FOR_DFTB3_SMOKE",
        "source_window": "fully_unrestrained_L2_NPT_1ns",
        "selected_time_ps": selected["time_ps"],
        "selected_potential_energy_kj_mol": selected["potential_energy_kj_mol"],
        "selected_distance_nm": distance,
        "selected_angle_deg": angle,
        "selected_gro_sha256": sha256(promoted),
        "qm_theory": "DFTB3",
        "slater_koster_set": "3ob-3-1",
        "scc_convergence_rescue": {
            "dftb_telec_K": DFTB_TELEC_K,
            "reason": "Amber18-supported low electronic temperature after zero-K job 61712026 failed at step 1 and 100 K job 61712561 retained one SCC warning in 20 steps",
            "interpretation": "numerical rescue only; warnings remain a hard failure",
        },
        "qmcharge": QMCHARGE,
        "spin": SPIN,
        "qm_atom_count": len(qm_atoms),
        "qm_elements": sorted(set(qm_elements)),
        "electron_count_including_link_h": electrons,
        "electron_parity": "even",
        "link_atom_count": 1,
        "active_residue_topology_charge": active_topology_charge,
        "ligand_topology_charge": ligand_topology_charge,
        "active_residue": f"{active.residue.name}:{active.residue.idx + 1}",
        "active_og1_atom": THR267_OG1,
        "ligand_atom_range": [L2_FIRST, L2_LAST],
        "reactive_atoms": {"thr267_og1": THR267_OG1, "l2_c": L2_REACTIVE_C, "l2_o": L2_REACTIVE_O},
        "boundary_bonds": boundary,
        "qmmask": qmmask,
        "qm_region_scope": "minimal_numerical_smoke",
        "production_qm_region_gate": (
            "Before Step1 TS/PMF, expand the mechanistic QM region to include Asp306 and Asp308 "
            "and test inclusion/sensitivity of Tyr146, Lys189 and Asn219."
        ),
        "interpretation": "Numerical preflight input only; not a TS, reaction coordinate, PMF, or barrier.",
    }
    (output / "qmmm_preflight_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    print(json.dumps(audit, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
