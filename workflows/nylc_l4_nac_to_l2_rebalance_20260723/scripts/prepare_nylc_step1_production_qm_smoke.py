#!/usr/bin/env python3
"""Prepare one-step DFTB3 smoke inputs for two paper-consistent production QM regions."""
import argparse
import hashlib
import json
import pathlib

import parmed as pmd

TYR146_N = 7156
LYS189_NZ = 7768
ASN219_ND2 = 8240
THR267_OG1 = 8961
ASP306_OD1 = 9572
ASP308_OD1 = 9591
L2_FIRST = 10273
L2_LAST = 10351
WATER_THR_ASN_1 = 50165
WATER_THR_ASN_2 = 51302
WATER_ASP_PAIR = 81218
EXPECTED_SOURCE_GRO_SHA256 = "182f68a41c72c490e0046323c9b8a0f8514ccf7af7dc4d578dcb87342c1c8e2d"
EXPECTED_REGIONS = {"core": {"qmcharge": -1, "atom_count": 110, "boundary_count": 3}, "network": {"qmcharge": 0, "atom_count": 161, "boundary_count": 7}}
SPIN = 1
DFTB_TELEC_K = 200.0

SIDECHAIN_NAMES = {
    "asp306": {"CB", "HB1", "HB2", "CG", "OD1", "OD2"},
    "asp308": {"CB", "HB1", "HB2", "CG", "OD1", "OD2"},
    "lys189": {"CB", "HB1", "HB2", "CG", "HG1", "HG2", "CD", "HD1", "HD2", "CE", "HE1", "HE2", "NZ", "HZ1", "HZ2", "HZ3"},
    "asn219": {"CB", "HB1", "HB2", "CG", "OD1", "ND2", "HD21", "HD22"},
}


def sha256(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


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
            result.append({"qm_atom": qm_atom.idx + 1, "qm_name": qm_atom.name, "mm_atom": mm_atom.idx + 1, "mm_name": mm_atom.name})
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
    gro = source / "selected_free_l2_nac.gro"
    if sha256(gro) != EXPECTED_SOURCE_GRO_SHA256:
        raise ValueError("434 ps selected GRO SHA mismatch")
    structure = pmd.load_file(str(source / "system.prmtop"), xyz=str(source / "selected_free_l2_nac.rst7"))
    if len(structure.atoms) != 133589:
        raise ValueError(f"atom count changed: {len(structure.atoms)}")

    thr = residue(structure, THR267_OG1, "THR", "OG1")
    d306 = residue(structure, ASP306_OD1, "ASP", "OD1")
    d308 = residue(structure, ASP308_OD1, "ASP", "OD1")
    lys = residue(structure, LYS189_NZ, "LYS", "NZ")
    asn = residue(structure, ASN219_ND2, "ASN", "ND2")
    tyr = residue(structure, TYR146_N, "TYR", "N")
    ligand = list(structure.atoms[L2_FIRST - 1:L2_LAST])
    if len(ligand) != 79 or {atom.residue.name for atom in ligand} != {"L2"}:
        raise ValueError("complete neutral 79-atom L2 region changed")
    water1 = list(residue(structure, WATER_THR_ASN_1, "SOL", "OW").atoms)
    water2 = list(residue(structure, WATER_THR_ASN_2, "SOL", "OW").atoms)
    water3 = list(residue(structure, WATER_ASP_PAIR, "SOL", "OW").atoms)
    core = list(thr.atoms) + ligand + sidechain(d306, SIDECHAIN_NAMES["asp306"]) + sidechain(d308, SIDECHAIN_NAMES["asp308"]) + water3
    network = core + sidechain(lys, SIDECHAIN_NAMES["lys189"]) + sidechain(asn, SIDECHAIN_NAMES["asn219"]) + list(tyr.atoms) + water1 + water2
    regions = {"core": core, "network": network}

    audit_regions = {}
    for name, atoms in regions.items():
        expected = EXPECTED_REGIONS[name]
        boundary = boundaries(structure, atoms)
        if len(atoms) != expected["atom_count"] or len(boundary) != expected["boundary_count"]:
            raise ValueError(f"{name} region size/boundary changed: {len(atoms)}/{len(boundary)}")
        qmmask = "@" + ",".join(str(atom.idx + 1) for atom in atoms)
        electrons = sum(atom.atomic_number for atom in atoms) + len(boundary) - expected["qmcharge"]
        if electrons % 2:
            raise ValueError(f"{name} has odd electron count {electrons}")
        (output / f"{name}_one_step.in").write_text(qmmm_input(f"{name} production-region numerical smoke", qmmask, expected["qmcharge"]))
        (output / f"{name}.qmmask").write_text(qmmask + "\n")
        audit_regions[name] = {
            "qmcharge": expected["qmcharge"],
            "spin": SPIN,
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
        "status": "READY_FOR_PRODUCTION_QM_REGION_SMOKE",
        "source_time_ps": 434,
        "source_gro_sha256": EXPECTED_SOURCE_GRO_SHA256,
        "qm_theory": "DFTB3",
        "slater_koster_set": "3ob-3-1",
        "dftb_telec_K": DFTB_TELEC_K,
        "protonation_microstate": "all_deprotonated_Asp306_Asp308",
        "regions": audit_regions,
        "water_oxygen_atoms": [WATER_THR_ASN_1, WATER_THR_ASN_2, WATER_ASP_PAIR],
        "direct_donor_gate": "SUPERSEDED_INVALID_DIRECT_DONOR_ASSUMPTION",
        "interpretation": "Numerical boundary/SCC smoke only; not a TS, RC, PMF, or barrier. Asp306 protonation remains to be built and compared.",
    }
    (output / "production_qm_region_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    print(json.dumps(audit, sort_keys=True))


if __name__ == "__main__":
    main()
