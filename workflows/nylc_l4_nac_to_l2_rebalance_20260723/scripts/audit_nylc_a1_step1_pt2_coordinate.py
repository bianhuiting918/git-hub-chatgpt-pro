#!/usr/bin/env python3
"""Audit whether the pinned A1 Step1 seed can support direct Nalpha-HG1 to leaving-N3 PT2."""
import argparse
import hashlib
import json
import math
import pathlib

import parmed as pmd

TASK_ROOT = pathlib.Path("/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723")
QMMM_SOURCE = TASK_ROOT / "a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285"
SEED_SOURCE = TASK_ROOT / "a1_activated_nac_20260726/qmmm/a1_unified_step1_qattack_extension/attempt_62021985"
START_RESTART = SEED_SOURCE / "scan/q06_1p45A/window.rst7"
PRMTOP = QMMM_SOURCE / "prepared/system.prmtop"
QMMM_AUDIT = QMMM_SOURCE / "prepared/qmmm_preflight_audit.json"
EXPECTED_RESTART_SHA256 = "548dca57cb623e69d5fefaa8e270267a5bd599bd5995113219cdf349c346e56c"
EXPECTED_PRMTOP_SHA256 = "a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0"
EXPECTED_QM_ATOMS = 146
N_ALPHA = 8949
THR267_OG1 = 8960
TRANSFERRED_HG1 = 8961
L2_C12 = 10287
L2_N3 = 10289
REACTION_ATOMS = (N_ALPHA, THR267_OG1, TRANSFERRED_HG1, L2_C12, L2_N3)


def sha256(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def xyz(atom):
    return (atom.xx, atom.xy, atom.xz)


def distance(structure, i, j):
    return math.dist(xyz(structure.atoms[i - 1]), xyz(structure.atoms[j - 1]))


def angle(structure, i, j, k):
    a, b, c = (structure.atoms[x - 1] for x in (i, j, k))
    u = tuple(x - y for x, y in zip(xyz(a), xyz(b)))
    v = tuple(x - y for x, y in zip(xyz(c), xyz(b)))
    denominator = math.sqrt(sum(x*x for x in u)) * math.sqrt(sum(x*x for x in v))
    cosine = sum(x*y for x, y in zip(u, v)) / denominator
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def qmmask_atoms(text):
    if not isinstance(text, str) or not text.startswith("@"):
        return set()
    return {int(item) for item in text[1:].split(",") if item.strip()}


def atom_record(structure, index1):
    atom = structure.atoms[index1 - 1]
    return {
        "index1": index1,
        "name": atom.name,
        "residue_name": atom.residue.name,
        "residue_number1": atom.residue.number + 1,
        "atomic_number": atom.atomic_number,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()

    for path in (PRMTOP, START_RESTART, QMMM_AUDIT):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"missing pinned PT2 source: {path}")
    hashes = {
        "restart": sha256(START_RESTART),
        "prmtop": sha256(PRMTOP),
    }
    hash_contract = (
        hashes["restart"] == EXPECTED_RESTART_SHA256
        and hashes["prmtop"] == EXPECTED_PRMTOP_SHA256
    )

    audit = json.loads(QMMM_AUDIT.read_text(encoding="utf-8"))
    qm_atoms = qmmask_atoms(audit.get("qmmask"))
    structure = pmd.load_file(str(PRMTOP), xyz=str(START_RESTART))
    bonds = {
        tuple(sorted((bond.atom1.idx + 1, bond.atom2.idx + 1)))
        for bond in structure.bonds
    }
    identities = {str(i): atom_record(structure, i) for i in REACTION_ATOMS}
    identity_contract = {
        "nalpha_is_N": identities[str(N_ALPHA)]["name"] == "N",
        "og1_is_OG1": identities[str(THR267_OG1)]["name"] == "OG1",
        "transferred_h_is_HG1": identities[str(TRANSFERRED_HG1)]["name"] == "HG1",
        "c12_is_C12": identities[str(L2_C12)]["name"] == "C12",
        "n3_is_N3": identities[str(L2_N3)]["name"] == "N3",
        "thr_atoms_same_residue": len({
            identities[str(i)]["residue_number1"]
            for i in (N_ALPHA, THR267_OG1, TRANSFERRED_HG1)
        }) == 1,
        "l2_atoms_same_residue": (
            identities[str(L2_C12)]["residue_number1"]
            == identities[str(L2_N3)]["residue_number1"]
        ),
    }
    bond_contract = {
        "nalpha_hg1_bond_present": tuple(sorted((N_ALPHA, TRANSFERRED_HG1))) in bonds,
        "og1_hg1_bond_absent": tuple(sorted((THR267_OG1, TRANSFERRED_HG1))) not in bonds,
        "c12_n3_bond_present": tuple(sorted((L2_C12, L2_N3))) in bonds,
    }
    qm_contract = {
        "qm_atom_count": len(qm_atoms),
        "expected_qm_atom_count": EXPECTED_QM_ATOMS,
        "all_reaction_atoms_in_qm": all(index in qm_atoms for index in REACTION_ATOMS),
        "qmcharge": audit.get("qmcharge"),
        "electron_count_including_link_h": audit.get("electron_count_including_link_h"),
        "link_atom_count": audit.get("link_atom_count"),
        "step1_qm_water_count": 0,
    }
    atom_or_qm_contract = (
        hash_contract
        and all(identity_contract.values())
        and all(bond_contract.values())
        and qm_contract["qm_atom_count"] == EXPECTED_QM_ATOMS
        and qm_contract["all_reaction_atoms_in_qm"]
        and qm_contract["qmcharge"] == 0
        and qm_contract["electron_count_including_link_h"] == 510
        and qm_contract["link_atom_count"] == 6
    )

    geometry = {
        "nalpha_n3_A": distance(structure, N_ALPHA, L2_N3),
        "hg1_n3_A": distance(structure, TRANSFERRED_HG1, L2_N3),
        "nalpha_hg1_A": distance(structure, N_ALPHA, TRANSFERRED_HG1),
        "nalpha_hg1_n3_deg": angle(structure, N_ALPHA, TRANSFERRED_HG1, L2_N3),
        "og1_c12_A": distance(structure, THR267_OG1, L2_C12),
        "c12_n3_A": distance(structure, L2_C12, L2_N3),
    }
    direct_preorganized = (
        atom_or_qm_contract
        and geometry["nalpha_n3_A"] <= 3.5
        and geometry["hg1_n3_A"] <= 2.5
        and geometry["nalpha_hg1_n3_deg"] >= 135.0
    )
    status = (
        "PASS_A1_PT2_ATOM_AND_QM_CONTRACT"
        if atom_or_qm_contract
        else "FAIL_A1_PT2_ATOM_OR_QM_CONTRACT"
    )
    relay_gate = (
        "PASS_A1_PT2_DIRECT_RELAY_PREORGANIZED"
        if direct_preorganized
        else "NOT_EVALUATED_A1_PT2_DIRECT_RELAY_NOT_PREORGANIZED"
    )
    result = {
        "schema_version": 1,
        "status": status,
        "scientific_status": "NOT_EVALUATED_TS_PMF_BARRIER_MECHANISM",
        "source": {
            "restart": str(START_RESTART),
            "prmtop": str(PRMTOP),
            "sha256": hashes,
            "hash_contract_pass": hash_contract,
        },
        "atoms": identities,
        "identity_contract": identity_contract,
        "bond_contract": bond_contract,
        "qm_contract": qm_contract,
        "geometry": geometry,
        "direct_relay_preorganization_definition": {
            "nalpha_n3_A_max": 3.5,
            "hg1_n3_A_max": 2.5,
            "nalpha_hg1_n3_deg_min": 135.0,
            "scope": "geometric preorganization only; not proton transfer or mechanism proof",
        },
        "direct_relay_gate": relay_gate,
        "decision": (
            "Eligible to design a bounded direct PT2 scout; do not call this proton transfer."
            if direct_preorganized
            else "Do not strongly pull HG1 directly to N3 from this seed; revisit relay geometry or candidate seed."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": status, "direct_relay_gate": relay_gate, "geometry": geometry}))
    raise SystemExit(0 if atom_or_qm_contract else 2)


if __name__ == "__main__":
    main()
