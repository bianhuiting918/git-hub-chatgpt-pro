#!/usr/bin/env python3
"""Prepare gradual DFTB3 restraints for the NylC Step1 tetrahedral-seed bracket."""
import argparse
import hashlib
import json
import pathlib

import parmed as pmd

SOURCE_NAME = "step1_dftb3_preflight_post61710861_job61712692"
SOURCE_RST7 = "selected_free_l2_nac.rst7"
EXPECTED_GRO_SHA256 = "182f68a41c72c490e0046323c9b8a0f8514ccf7af7dc4d578dcb87342c1c8e2d"
THR267_OG1 = 8961
THR267_HG1 = 8962
L2_C12 = 10287
L2_O2 = 10288
L2_N3 = 10289
ASP306_OD1 = 9572
ASP306_OD2 = 9573
WATER_ASP_PAIR = (81218, 81219, 81220)
EXPECTED_QM_ATOMS = 110
QMCHARGE = -1
WINDOWS = [
    {"name": "w00", "attack_A": 3.30, "proton_A": 3.40, "force": 20.0},
    {"name": "w01", "attack_A": 2.80, "proton_A": 2.70, "force": 25.0},
    {"name": "w02", "attack_A": 2.30, "proton_A": 1.90, "force": 30.0},
    {"name": "w03", "attack_A": 1.75, "proton_A": 1.10, "force": 35.0},
]
SIDECHAIN_NAMES = {"CB", "HB1", "HB2", "CG", "OD1", "OD2"}


def sha256(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def require_atom(structure, number, resname, atomname):
    atom = structure.atoms[number - 1]
    if atom.residue.name != resname or atom.name != atomname:
        raise ValueError(f"identity mismatch at {number}: {atom.residue.name}:{atom.name}")
    return atom


def sidechain(residue):
    atoms = [atom for atom in residue.atoms if atom.name in SIDECHAIN_NAMES]
    if {atom.name for atom in atoms} != SIDECHAIN_NAMES:
        raise ValueError(f"incomplete ASP sidechain at residue {residue.idx + 1}")
    return atoms


def boundaries(structure, atoms):
    selected = {atom.idx for atom in atoms}
    result = []
    for bond in structure.bonds:
        left, right = bond.atom1, bond.atom2
        if (left.idx in selected) != (right.idx in selected):
            q = left if left.idx in selected else right
            m = right if left.idx in selected else left
            result.append({"qm_atom": q.idx + 1, "mm_atom": m.idx + 1})
    return result


def amber_input(name, qmmask):
    return f"""NylC Step1 {name} constrained tetrahedral-seed bracket
&cntrl
  imin=1, maxcyc=100, ncyc=50,
  ntb=1, cut=10.0, ntpr=5,
  ifqnt=1, nmropt=1,
/
&qmmm
  qmmask='{qmmask}',
  qmcharge={QMCHARGE}, spin=1,
  qm_theory='DFTB3',
  dftb_telec=200.0,
  qmshake=0,
/
&wt type='END' /
DISANG=../input/{name}.RST
LISTOUT=POUT
"""


def distance_restraint(atom1, atom2, target, force):
    r1 = max(0.0, target - 0.30)
    r2 = max(0.0, target - 0.03)
    r3 = target + 0.03
    r4 = target + 0.30
    return f"&rst iat={atom1},{atom2}, r1={r1:.3f}, r2={r2:.3f}, r3={r3:.3f}, r4={r4:.3f}, rk2={force:.2f}, rk3={force:.2f}, /"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=pathlib.Path, required=True)
    parser.add_argument("--acceptor", choices=("OD1", "OD2"), required=True)
    parser.add_argument("--output-dir", type=pathlib.Path, required=True)
    args = parser.parse_args()
    source = args.source_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)

    source_pass = json.loads((source / "PASS.json").read_text())
    if source.name != SOURCE_NAME or source_pass.get("status") != "PASS_DFTB3_PREFLIGHT":
        raise ValueError("authoritative all-deprotonated DFTB3 preflight source is not locked")
    if sha256(source / "selected_free_l2_nac.gro") != EXPECTED_GRO_SHA256:
        raise ValueError("authoritative 434 ps explicit-water NAC SHA mismatch")
    structure = pmd.load_file(str(source / "system.prmtop"), xyz=str(source / SOURCE_RST7))
    if len(structure.atoms) != 133589:
        raise ValueError(f"full explicit-solvent atom count changed: {len(structure.atoms)}")

    thr = require_atom(structure, THR267_OG1, "THR", "OG1").residue
    require_atom(structure, THR267_HG1, "THR", "HG1")
    require_atom(structure, L2_C12, "L2", "C12")
    require_atom(structure, L2_O2, "L2", "O2")
    require_atom(structure, L2_N3, "L2", "N3")
    d306 = require_atom(structure, ASP306_OD1, "ASP", "OD1").residue
    require_atom(structure, ASP306_OD2, "ASP", "OD2")
    d308 = require_atom(structure, 9591, "ASP", "OD1").residue
    ligand = list(structure.atoms[10272:10351])
    water = list(require_atom(structure, WATER_ASP_PAIR[0], "SOL", "OW").residue.atoms)
    qm_atoms = list(thr.atoms) + ligand + sidechain(d306) + sidechain(d308) + water
    boundary = boundaries(structure, qm_atoms)
    if len(qm_atoms) != EXPECTED_QM_ATOMS or len(boundary) != 3:
        raise ValueError(f"production core changed: atoms={len(qm_atoms)} boundaries={len(boundary)}")
    electrons = sum(atom.atomic_number for atom in qm_atoms) + len(boundary) - QMCHARGE
    if electrons != 388 or electrons % 2:
        raise ValueError(f"unexpected electron count {electrons}")
    qmmask = "@" + ",".join(str(atom.idx + 1) for atom in qm_atoms)
    acceptor_atom = ASP306_OD1 if args.acceptor == "OD1" else ASP306_OD2

    for window in WINDOWS:
        name = window["name"]
        (output / f"{name}.in").write_text(amber_input(name, qmmask))
        restraints = [
            distance_restraint(THR267_OG1, L2_C12, window["attack_A"], window["force"]),
            distance_restraint(THR267_HG1, acceptor_atom, window["proton_A"], window["force"]),
        ]
        (output / f"{name}.RST").write_text("\n".join(restraints) + "\n")

    manifest = {
        "schema_version": 1,
        "candidate_id": "nylc_C18_trueT267_freeGS",
        "status": "READY_FOR_STEP1_CONSTRAINED_BRACKET",
        "source": SOURCE_NAME,
        "source_coordinate": SOURCE_RST7,
        "source_gro_sha256": EXPECTED_GRO_SHA256,
        "solvent_model": "full_explicit_TIP3P",
        "protonation_microstate": "Asp306-_Asp308-",
        "acceptor_hypothesis": args.acceptor,
        "acceptor_atom": acceptor_atom,
        "reactive_atoms": {"thr267_og1": THR267_OG1, "thr267_hg1": THR267_HG1, "l2_c12": L2_C12, "l2_o2": L2_O2, "l2_n3": L2_N3},
        "qm_region": {"atom_count": len(qm_atoms), "boundary_count": len(boundary), "qmcharge": QMCHARGE, "electron_count": electrons, "bridge_water_atoms": list(WATER_ASP_PAIR), "qmmask": qmmask},
        "windows": WINDOWS,
        "restraint_scope": "attack OG1-C12 and proton approach HG1-Asp306O only; OG1-HG1 bond is observed, not restrained",
        "interpretation": "Constrained minimization bracket seed only; not a TS, committor, PMF, or barrier.",
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, sort_keys=True))


if __name__ == "__main__":
    raise SystemExit(main())
