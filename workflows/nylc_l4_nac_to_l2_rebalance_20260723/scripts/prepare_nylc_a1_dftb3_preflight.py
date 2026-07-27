#!/usr/bin/env python3
"""Prepare the frozen NylC A1 representative frame for a minimal DFTB3 smoke."""
import argparse
import hashlib
import json
import os
import pathlib

import parmed as pmd

TASK_ROOT = pathlib.Path("/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723")
REPRESENTATIVE_ROOT = TASK_ROOT / "a1_activated_nac_20260726/representative_frame/attempt_61990814"
TOPOLOGY_ROOT = TASK_ROOT / "a1_activated_nac_20260726/em/attempt_61968026_1_61968026/nac_evt25_time1462ps/input"
EXPECTED_REPRESENTATIVE_STATUS = "PASS_A1_REPRESENTATIVE_FRAME_EXTRACTION"
EXPECTED_SCIENTIFIC_STATUS = "PASS_A1_REPRESENTATIVE_NAC_FRAME"
EXPECTED_GRO_SHA256 = "7b095b7d36a25327558b4bbab20a98a43edaab51d0cacaf7a888ff886bafe679"
EXPECTED_REPRESENTATIVE_AUDIT_SHA256 = "8a61a57be0537a5fb0aac8f2379452fdb6980468cb06659c0953224d3a0aefa3"
EXPECTED_ACTIVE_GLOBAL_RESID = 622
EXPECTED_ACTIVE_ORIGINAL_RESID = 267
EXPECTED_TOPOLOGY_HASHES = {
    "topol.top": "af98733e218a8f83d0a5c46120d9d230a16c222f76d230654d44b542288cc205",
    "topol_Protein_chain_H.itp": "8fd4398af1356b515720c1da3126d08b7795b24b3c112d98ef3235afa57c8179",
    "PA66_L2_GMX.itp": "b0e753c60fd4b71c282d21cc6106a15e73d91d12a20d80e92dd01516162eb301",
}
EXPECTED_SYSTEM_ATOMS = 133589
PROTEIN_ATOMS = 10272
THR267_OG1 = 8960
L2_FIRST = 10273
L2_LAST = 10351
L2_REACTIVE_C = 10287
L2_REACTIVE_O = 10288
L2_REACTIVE_N = 10289
QMCHARGE = 0
SPIN = 1
DFTB_TELEC_K = 200.0
EXPECTED_QM_ATOMS = 94
EXPECTED_ELECTRONS_WITH_LINKS = 314
MAX_ALLOWED_BOND_LENGTH_A = 2.0
ATOMIC_NUMBERS = {"H": 1, "C": 6, "N": 7, "O": 8}


def sha256(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def element(atom):
    number = int(getattr(atom, "atomic_number", 0) or 0)
    for symbol, expected in ATOMIC_NUMBERS.items():
        if number == expected:
            return symbol
    letters = "".join(character for character in atom.name if character.isalpha()).upper()
    if not letters or letters[0] not in ATOMIC_NUMBERS:
        raise ValueError(f"unsupported 3ob-3-1 QM element at atom {atom.idx + 1}: {atom.name}")
    return letters[0]


def qmmm_input(title, maxcyc, qmmask):
    return f"""{title}
&cntrl
  imin=1, maxcyc={maxcyc}, ncyc={max(1, maxcyc // 2)},
  ntb=1, cut=10.0, ntpr=1, ntxo=1,
  ifqnt=1,
/
&qmmm
  qmmask='{qmmask}',
  qmcharge={QMCHARGE},
  spin={SPIN},
  qm_theory='DFTB3',
  dftb_telec=200.0,
  qmshake=0,
/
"""


def atom_by_name(residue, name):
    matches = [atom for atom in residue.atoms if atom.name == name]
    if len(matches) != 1:
        raise ValueError(f"expected one {residue.name}:{name}, found {len(matches)}")
    return matches[0]


def derive_qm_contract(structure):
    if len(structure.atoms) != EXPECTED_SYSTEM_ATOMS:
        raise ValueError(f"system atom count {len(structure.atoms)} != {EXPECTED_SYSTEM_ATOMS}")
    if PROTEIN_ATOMS != L2_FIRST - 1:
        raise ValueError("frozen protein/L2 boundary is inconsistent")
    og1 = structure.atoms[THR267_OG1 - 1]
    if (
        og1.name != "OG1"
        or og1.residue.name != "THR"
        or og1.residue.idx + 1 != EXPECTED_ACTIVE_GLOBAL_RESID
    ):
        raise ValueError("Thr267 OG1 global identity changed")
    active = og1.residue
    if len(active.atoms) != 15:
        raise ValueError(f"A1 Thr267 atom count {len(active.atoms)} != 15")
    nalpha = atom_by_name(active, "N")
    n_hydrogens = sorted(atom.name for atom in nalpha.bond_partners if atom.name.startswith("H"))
    if n_hydrogens != ["H1", "H2", "HG1"]:
        raise ValueError(f"N bonded hydrogens {n_hydrogens} != H1/H2/HG1")
    og1_bonded = sorted(atom.name for atom in og1.bond_partners)
    if og1_bonded != ["CB"]:
        raise ValueError(f"OG1 bonded atoms {og1_bonded} != CB")
    ligand = list(structure.atoms[L2_FIRST - 1:L2_LAST])
    if len(ligand) != 79 or {atom.residue.name for atom in ligand} != {"L2"}:
        raise ValueError("complete neutral 79-atom L2 ordering changed")
    for index1, name in (
        (L2_REACTIVE_C, "C12"),
        (L2_REACTIVE_O, "O2"),
        (L2_REACTIVE_N, "N3"),
    ):
        atom = structure.atoms[index1 - 1]
        if atom.name != name or atom.residue.name != "L2":
            raise ValueError(f"reactive atom {index1} is not L2:{name}")
    active_charge = sum(atom.charge for atom in active.atoms)
    ligand_charge = sum(atom.charge for atom in ligand)
    if abs(active_charge) > 1.0e-4:
        raise ValueError(f"A1 Thr267 patch net charge is not 0: {active_charge}")
    if abs(ligand_charge) > 1.0e-4:
        raise ValueError(f"L2 net charge is not 0: {ligand_charge}")
    bond_lengths = []
    for bond in structure.bonds:
        left, right = bond.atom1, bond.atom2
        distance = sum(
            (float(left.xx[axis]) - float(right.xx[axis])) ** 2
            for axis in range(3)
        ) ** 0.5
        bond_lengths.append(distance)
    max_bond_length_A = max(bond_lengths)
    bond_count_gt_3A = sum(distance > 3.0 for distance in bond_lengths)
    if max_bond_length_A > MAX_ALLOWED_BOND_LENGTH_A or bond_count_gt_3A:
        raise ValueError(
            "PBC-split or otherwise invalid bonded geometry: "
            f"max_bond_length_A={max_bond_length_A}; "
            f"bond_count_gt_3A={bond_count_gt_3A}"
        )
    qm_atoms = list(active.atoms) + ligand
    if len(qm_atoms) != EXPECTED_QM_ATOMS or len({atom.idx for atom in qm_atoms}) != EXPECTED_QM_ATOMS:
        raise ValueError(f"QM atom count {len(qm_atoms)} != {EXPECTED_QM_ATOMS}")
    qm_indices = {atom.idx for atom in qm_atoms}
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
    elements = [element(atom) for atom in qm_atoms]
    link_atom_count = len(boundary)
    electrons = sum(ATOMIC_NUMBERS[symbol] for symbol in elements) + link_atom_count - QMCHARGE
    if electrons != EXPECTED_ELECTRONS_WITH_LINKS:
        raise ValueError(
            f"derived electron count {electrons} != {EXPECTED_ELECTRONS_WITH_LINKS}"
        )
    if electrons % 2 or SPIN != 1:
        raise ValueError("derived electron parity is incompatible with singlet spin")
    return {
        "qm_atoms": qm_atoms,
        "boundary_bonds": boundary,
        "elements": elements,
        "electron_count_including_link_h": electrons,
        "active_residue_topology_charge": active_charge,
        "ligand_topology_charge": ligand_charge,
        "n_bonded_hydrogens": n_hydrogens,
        "og1_bonded_atoms": og1_bonded,
        "max_bond_length_A": max_bond_length_A,
        "bond_count_gt_3A": bond_count_gt_3A,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--coordinate", type=pathlib.Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    coordinate = args.coordinate.resolve()
    if not coordinate.is_file():
        raise ValueError(f"whole-system coordinate is missing: {coordinate}")
    output.mkdir(parents=True, exist_ok=False)
    pass_path = REPRESENTATIVE_ROOT / "PASS.json"
    representative_audit_path = REPRESENTATIVE_ROOT / "A1_REPRESENTATIVE_FRAME_AUDIT.json"
    gro = REPRESENTATIVE_ROOT / "representative_354ps.gro"
    authority = json.loads(pass_path.read_text(encoding="utf-8"))
    if authority.get("status") != EXPECTED_REPRESENTATIVE_STATUS:
        raise ValueError("representative extraction status is not PASS")
    if authority.get("scientific_status") != EXPECTED_SCIENTIFIC_STATUS:
        raise ValueError("representative scientific status is not PASS")
    if sha256(gro) != EXPECTED_GRO_SHA256:
        raise ValueError("representative GRO hash differs from frozen authority")
    if sha256(representative_audit_path) != EXPECTED_REPRESENTATIVE_AUDIT_SHA256:
        raise ValueError("representative audit hash differs from frozen authority")
    representative_audit = json.loads(
        representative_audit_path.read_text(encoding="utf-8")
    )
    active_mapping = representative_audit.get("atom_mapping", {}).get("thr267_og1", {})
    expected_mapping = {
        "index1": THR267_OG1,
        "resname": "THR",
        "name": "OG1",
        "source_global_resid": EXPECTED_ACTIVE_GLOBAL_RESID,
        "extracted_original_resid": EXPECTED_ACTIVE_ORIGINAL_RESID,
    }
    for key, expected in expected_mapping.items():
        if active_mapping.get(key) != expected:
            raise ValueError(f"representative dual-namespace mapping mismatch for {key}")
    for name, expected in EXPECTED_TOPOLOGY_HASHES.items():
        observed = sha256(TOPOLOGY_ROOT / name)
        if observed != expected:
            raise ValueError(f"topology hash mismatch for {name}: {observed}")
    previous = pathlib.Path.cwd()
    os.chdir(TOPOLOGY_ROOT)
    try:
        structure = pmd.load_file(
            str(TOPOLOGY_ROOT / "topol.top"),
            xyz=str(coordinate),
        )
    finally:
        os.chdir(previous)
    contract = derive_qm_contract(structure)
    structure.save(str(output / "system.prmtop"), overwrite=True)
    structure.save(str(output / "representative_354ps.rst7"), overwrite=True)
    qmmask = "@" + ",".join(str(atom.idx + 1) for atom in contract["qm_atoms"])
    (output / "01_qmmm_one_step.in").write_text(
        qmmm_input("NylC A1 minimal DFTB3 one-step numerical preflight", maxcyc=1, qmmask=qmmask),
        encoding="utf-8",
    )
    (output / "02_qmmm_20_step.in").write_text(
        qmmm_input("NylC A1 minimal DFTB3 twenty-step numerical preflight", maxcyc=20, qmmask=qmmask),
        encoding="utf-8",
    )
    audit = {
        "schema_version": 1,
        "status": "READY_A1_DFTB3_NUMERICAL_PREFLIGHT",
        "source_root": str(REPRESENTATIVE_ROOT),
        "source_status": authority["status"],
        "source_scientific_status": authority["scientific_status"],
        "source_gro_sha256": EXPECTED_GRO_SHA256,
        "coordinate_sha256": sha256(coordinate),
        "coordinate_preprocessing": "gmx trjconv -pbc mol -ur compact",
        "representative_audit_sha256": EXPECTED_REPRESENTATIVE_AUDIT_SHA256,
        "active_global_resid": EXPECTED_ACTIVE_GLOBAL_RESID,
        "active_original_resid": EXPECTED_ACTIVE_ORIGINAL_RESID,
        "topology_sha256": EXPECTED_TOPOLOGY_HASHES,
        "qm_theory": "DFTB3",
        "slater_koster_set": "3ob-3-1",
        "dftb_telec_K": DFTB_TELEC_K,
        "qmcharge": QMCHARGE,
        "spin": SPIN,
        "qm_atom_count": len(contract["qm_atoms"]),
        "electron_count_including_link_h": contract["electron_count_including_link_h"],
        "link_atom_count": len(contract["boundary_bonds"]),
        "boundary_bonds": contract["boundary_bonds"],
        "active_residue_topology_charge": contract["active_residue_topology_charge"],
        "ligand_topology_charge": contract["ligand_topology_charge"],
        "n_bonded_hydrogens": contract["n_bonded_hydrogens"],
        "og1_bonded_atoms": contract["og1_bonded_atoms"],
        "max_bond_length_A": contract["max_bond_length_A"],
        "bond_count_gt_3A": contract["bond_count_gt_3A"],
        "qmmask": qmmask,
        "qm_region_scope": "minimal Thr267 plus complete L2 numerical preflight",
        "production_qm_region_gate": (
            "Before production Step1 work, expand at least through Asp306 and Asp308 "
            "and test network sensitivity including Tyr146, Lys189 and Asn219."
        ),
        "interpretation": (
            "Numerical preflight only; not a TS, reaction coordinate, PMF, barrier, "
            "or evidence for proton transfer."
        ),
    }
    (output / "qmmm_preflight_audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
