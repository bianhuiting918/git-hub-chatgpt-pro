#!/usr/bin/env python3
"""Prepare DFTB3 one-step smokes for the neutral ASH306/Asp308- microstate."""
import argparse
import hashlib
import json
import pathlib

import parmed as pmd

SYSTEM_PRMTOP_SHA256 = "bf6fa73b819da92c1cb0fdc6f97c17bd3db37e1a49c0127ecac39fd9a8d1b8ae"
SYSTEM_RST7_SHA256 = "80a8c9371cf59799c25927c5d1da95e6ecbb4662ba3615f547fa752cf825289d"
THR267_OG1 = 8961
ASH306_OD1 = 9572
ASP308_OD1 = 9592
L2_FIRST = 10274
L2_LAST = 10352
TYR146_N = 7156
LYS189_NZ = 7768
ASN219_ND2 = 8240
WATER_THR_ASN_1 = 50166
WATER_THR_ASN_2 = 51303
WATER_ASP_PAIR = 81219
EXPECTED_ATOMS = 133589
EXPECTED_REGIONS = {"core": {"qmcharge": 0, "atom_count": 111, "boundary_count": 3, "electrons": 388}, "network": {"qmcharge": 1, "atom_count": 162, "boundary_count": 7, "electrons": 570}}
DFTB_TELEC_K = 200.0

SIDECHAIN_NAMES = {
    "ash306": {"CB", "HB1", "HB2", "CG", "OD1", "OD2", "HD2"},
    "asp308": {"CB", "HB1", "HB2", "CG", "OD1", "OD2"},
    "lys189": {"CB", "HB1", "HB2", "CG", "HG1", "HG2", "CD", "HD1", "HD2", "CE", "HE1", "HE2", "NZ", "HZ1", "HZ2", "HZ3"},
    "asn219": {"CB", "HB1", "HB2", "CG", "OD1", "ND2", "HD21", "HD22"},
}


def sha256(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def residue(structure, atom_number, resname, atomname):
    atom = structure.atoms[atom_number - 1]
    if atom.residue.name != resname or atom.name != atomname:
        raise ValueError(f"identity mismatch at {atom_number}: {atom.residue.name}:{atom.name}")
    return atom.residue


def sidechain(residue_value, names):
    atoms = [atom for atom in residue_value.atoms if atom.name in names]
    missing = names - {atom.name for atom in atoms}
    if missing:
        raise ValueError(f"missing sidechain atoms in {residue_value.name}: {sorted(missing)}")
    return atoms


def boundaries(structure, qm_atoms):
    indices = {atom.idx for atom in qm_atoms}
    result = []
    for bond in structure.bonds:
        left, right = bond.atom1, bond.atom2
        if (left.idx in indices) != (right.idx in indices):
            qm_atom = left if left.idx in indices else right
            mm_atom = right if left.idx in indices else left
            result.append({
                "qm_atom": qm_atom.idx + 1,
                "qm_name": qm_atom.name,
                "mm_atom": mm_atom.idx + 1,
                "mm_name": mm_atom.name,
            })
    return result


def qmmm_input(title, qmmask, qmcharge):
    return f"""NylC {title}
&cntrl
  imin=1, maxcyc=1, ncyc=1,
  ntb=1, cut=10.0, ntpr=1,
  ifqnt=1,
/
&qmmm
  qmmask='{qmmask}',
  qmcharge={qmcharge},
  spin=1,
  qm_theory='DFTB3',
  dftb_telec=200.0,
  qmshake=0,
/
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=pathlib.Path, required=True)
    parser.add_argument("--output-dir", type=pathlib.Path, required=True)
    args = parser.parse_args()
    source = args.source_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    prmtop = source / "system.prmtop"
    rst7 = source / "system.rst7"
    if sha256(prmtop) != SYSTEM_PRMTOP_SHA256 or sha256(rst7) != SYSTEM_RST7_SHA256:
        raise ValueError("authoritative neutral ASH306 Amber input SHA mismatch")
    structure = pmd.load_file(str(prmtop), xyz=str(rst7))
    if len(structure.atoms) != EXPECTED_ATOMS:
        raise ValueError(f"atom count changed: {len(structure.atoms)}")
    total_charge = sum(atom.charge for atom in structure.atoms)
    if abs(total_charge) > 1.0e-4:
        raise ValueError(f"full-system charge is not neutral: {total_charge}")

    thr = residue(structure, THR267_OG1, "THR", "OG1")
    d306 = residue(structure, ASH306_OD1, "ASH", "OD1")
    d308 = residue(structure, ASP308_OD1, "ASP", "OD1")
    lys = residue(structure, LYS189_NZ, "LYS", "NZ")
    asn = residue(structure, ASN219_ND2, "ASN", "ND2")
    tyr = residue(structure, TYR146_N, "TYR", "N")
    ligand = list(structure.atoms[L2_FIRST - 1:L2_LAST])
    if len(ligand) != 79 or {atom.residue.name for atom in ligand} != {"L2"}:
        raise ValueError("complete 79-atom L2 region changed")
    water1 = list(residue(structure, WATER_THR_ASN_1, "SOL", "OW").atoms)
    water2 = list(residue(structure, WATER_THR_ASN_2, "SOL", "OW").atoms)
    water3 = list(residue(structure, WATER_ASP_PAIR, "SOL", "OW").atoms)

    charge_checks = {
        "thr267": sum(atom.charge for atom in thr.atoms),
        "ash306": sum(atom.charge for atom in d306.atoms),
        "asp308": sum(atom.charge for atom in d308.atoms),
        "l2": sum(atom.charge for atom in ligand),
    }
    expected_charges = {"thr267": 1.0, "ash306": 0.0, "asp308": -1.0, "l2": 0.0}
    for name, expected in expected_charges.items():
        if abs(charge_checks[name] - expected) > 1.0e-4:
            raise ValueError(f"{name} charge mismatch: {charge_checks[name]} vs {expected}")

    core = list(thr.atoms) + ligand + sidechain(d306, SIDECHAIN_NAMES["ash306"]) + sidechain(d308, SIDECHAIN_NAMES["asp308"]) + water3
    network = core + sidechain(lys, SIDECHAIN_NAMES["lys189"]) + sidechain(asn, SIDECHAIN_NAMES["asn219"]) + list(tyr.atoms) + water1 + water2
    regions = {"core": core, "network": network}
    audit_regions = {}
    for name, atoms in regions.items():
        expected = EXPECTED_REGIONS[name]
        boundary = boundaries(structure, atoms)
        electrons = sum(atom.atomic_number for atom in atoms) + len(boundary) - expected["qmcharge"]
        if len(atoms) != expected["atom_count"] or len(boundary) != expected["boundary_count"]:
            raise ValueError(f"{name} size/boundary changed: {len(atoms)}/{len(boundary)}")
        if electrons != expected["electrons"] or electrons % 2:
            raise ValueError(f"{name} electron count changed: {electrons}")
        qmmask = "@" + ",".join(str(atom.idx + 1) for atom in atoms)
        (output / f"{name}_one_step.in").write_text(qmmm_input(f"{name} ASH306 numerical smoke", qmmask, expected["qmcharge"]))
        (output / f"{name}.qmmask").write_text(qmmask + "\n")
        audit_regions[name] = {
            "qmcharge": expected["qmcharge"],
            "spin": 1,
            "explicit_qm_atom_count": len(atoms),
            "boundary_count": len(boundary),
            "boundary_bonds": boundary,
            "electron_count_including_link_h": electrons,
            "electron_parity": "even",
            "topology_partial_charge_sum": sum(atom.charge for atom in atoms),
            "qmmask": qmmask,
        }

    audit = {
        "schema_version": 1,
        "candidate_id": "nylc_C18_trueT267_freeGS",
        "status": "READY_FOR_ASH306_QM_REGION_SMOKE",
        "source_job_id": 61718715,
        "system_prmtop_sha256": SYSTEM_PRMTOP_SHA256,
        "system_rst7_sha256": SYSTEM_RST7_SHA256,
        "protonation_microstate": "Asp306H_Asp308-",
        "full_system_total_charge_e": total_charge,
        "residue_charge_checks": charge_checks,
        "qm_theory": "DFTB3",
        "slater_koster_set": "3ob-3-1",
        "dftb_telec_K": DFTB_TELEC_K,
        "regions": audit_regions,
        "interpretation": "Numerical one-step microstate smoke only; absolute energy is not compared across different proton/ion compositions and this is not a TS, RC, PMF, or barrier.",
    }
    (output / "ash306_qm_region_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    print(json.dumps(audit, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
