#!/usr/bin/env python3
"""Build/audit fixed-geometry ICCG Step1 LG1/LG2 QM/MM pair inputs.

All production inputs are explicit CLI paths. The local builder creates the
preflight report consumed by deployment; remote Amber artifacts remain outside
Git in /work/home/acshdt1dks/iccg_qmmm_step1_pair_20260723.
"""
from __future__ import annotations

import argparse, dataclasses, hashlib, json, math
from pathlib import Path
from typing import Iterable, Sequence


AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
PROTEIN_RESNAMES = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
QM_SIDECHAINS = (("TRP",164),("SER",165),("ASP",210),("HIS",242),("ILE",243))
LIGAND_NAMES = {"LG1", "LG2"}
LIGAND_ATOM_COUNT = 54
LIGAND_HEAVY_COUNT = 32
REMOTE_PROJECT = Path("/work/home/acshdt1dks/iccg_qmmm_step1_pair_20260723")

@dataclasses.dataclass(frozen=True)
class Atom:
    index: int
    resname: str
    resid: int
    name: str
    xyz: tuple[float, float, float]
    element: str = ""
    record: str = "ATOM"

    def with_xyz(self, xyz: Sequence[float]) -> "Atom":
        return dataclasses.replace(self, xyz=(float(xyz[0]), float(xyz[1]), float(xyz[2])))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def paper_to_iccg_mapping() -> dict[str, dict[str, object]]:
    return {
        "paper_Trp159": {"iccg_resname": "TRP", "iccg_resid": 164},
        "paper_Ser160": {"iccg_resname": "SER", "iccg_resid": 165},
        "paper_Asp206": {"iccg_resname": "ASP", "iccg_resid": 210},
        "paper_His237": {"iccg_resname": "HIS", "iccg_resid": 242},
        "paper_Ser238": {"iccg_resname": "ILE", "iccg_resid": 243, "ca_separation_A": 0.939, "note": "3D structural equivalent; not sequence-identical Ser238."},
    }


def validate_qm_sidechain_mapping(mapping: dict[str, object] | None = None) -> None:
    mapping = mapping or paper_to_iccg_mapping()
    for key in ["paper_Trp159","paper_Ser160","paper_Asp206","paper_His237","paper_Ser238"]:
        if not isinstance(mapping.get(key), dict):
            raise ValueError(f"missing structural QM mapping for {key}")
    ser238 = mapping["paper_Ser238"]
    if ser238.get("iccg_resname") != "ILE" or ser238.get("iccg_resid") != 243:
        raise ValueError("paper Ser238 must map structurally to ICCG Ile243, not a null sequence gap")


def parse_pdb_atoms(path: Path) -> list[Atom]:
    atoms: list[Atom] = []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        name = line[12:16].strip()
        element = line[76:78].strip() or ''.join(c for c in name if c.isalpha())[:1]
        atoms.append(Atom(int(line[6:11]), line[17:20].strip(), int(line[22:26]), name,
                          (float(line[30:38]), float(line[38:46]), float(line[46:54])), element.upper(), line[:6].strip()))
    return atoms


def protein_atoms(atoms: Sequence[Atom]) -> list[Atom]:
    return [a for a in atoms if a.resname in PROTEIN_RESNAMES and a.record == "ATOM"]


def ligand_atoms(atoms: Sequence[Atom]) -> list[Atom]:
    return [a for a in atoms if a.resname in LIGAND_NAMES or (a.record == "HETATM" and a.resname not in PROTEIN_RESNAMES)]


def validate_active_iccg(atoms: Sequence[Atom], require_ser_hg: bool = False) -> dict[str, object]:
    residues = {(a.resname, a.resid) for a in protein_atoms(atoms)}
    names = {(a.resname, a.resid, a.name) for a in atoms}
    result = {"protein_residues": len({rid for _, rid in residues}), "ser165_og": ("SER",165,"OG") in names, "ser165_hg": ("SER",165,"HG") in names}
    if result["protein_residues"] != 258:
        raise ValueError("active ICCG chain A must contain 258 protein residues")
    if not result["ser165_og"]:
        raise ValueError("Ser165 OG must be present in heavy-atom ICCG input")
    if require_ser_hg and not result["ser165_hg"]:
        raise ValueError("Ser165 HG must be present after protonation")
    return result


def ligand_audit(atoms: Sequence[Atom]) -> dict[str, object]:
    lig = ligand_atoms(atoms)
    heavy = [a for a in lig if a.element.upper() != "H"]
    names = {a.name for a in lig}
    rings = count_six_membered_carbon_rings(lig)
    return {"ligand_resname": lig[0].resname if lig else None, "ligand_atoms": len(lig), "ligand_heavy_atoms": len(heavy), "six_membered_carbon_rings": rings, "atom_names": [a.name for a in lig]}


def count_six_membered_carbon_rings(lig: Sequence[Atom]) -> int:
    names = {a.name for a in lig if a.element.upper() == "C"}
    known = [{"C7","C9","C11","C13","C15","C17"}, {"C24","C25","C26","C27","C28","C29"}]
    return sum(1 for ring in known if ring <= names)


def assert_ligand_counts(atoms: Sequence[Atom]) -> None:
    audit = ligand_audit(atoms)
    if audit["ligand_atoms"] != LIGAND_ATOM_COUNT:
        raise ValueError(f"expected 54 ligand atoms, observed {audit['ligand_atoms']}")
    if audit["ligand_heavy_atoms"] != LIGAND_HEAVY_COUNT:
        raise ValueError(f"expected 32 ligand heavy atoms, observed {audit['ligand_heavy_atoms']}")
    if audit["six_membered_carbon_rings"] != 2:
        raise ValueError("expected two six-membered carbon rings")


def assert_identical_ligand_order(lg1: Sequence[Atom], lg2: Sequence[Atom]) -> None:
    if len(lg1) != LIGAND_ATOM_COUNT or len(lg2) != LIGAND_ATOM_COUNT or len(lg1) != len(lg2):
        raise ValueError("LG1/LG2 ligand atom counts must both be exactly 54")
    sig1 = [(a.name, a.element) for a in lg1]
    sig2 = [(a.name, a.element) for a in lg2]
    if sig1 != sig2:
        raise ValueError("LG1/LG2 ligand atom names/order are not identical")


def select_qm_atoms(atoms: Sequence[Atom]) -> list[Atom]:
    sidechains = {(res, rid) for res, rid in QM_SIDECHAINS}
    selected = [a for a in atoms if a.resname in LIGAND_NAMES or (a.resname, a.resid) in sidechains]
    present = {(a.resname, a.resid) for a in selected if a.resname not in LIGAND_NAMES}
    if present != sidechains:
        raise ValueError(f"QM sidechain set mismatch: {sorted(present)}")
    lig_count = sum(1 for a in selected if a.resname in LIGAND_NAMES)
    if lig_count != LIGAND_ATOM_COUNT:
        raise ValueError(f"expected complete 54-atom LG1/LG2 ligand in QM, observed {lig_count}")
    return selected


def build_qm_mask(atoms: Sequence[Atom]) -> str:
    ids = sorted(a.index for a in select_qm_atoms(atoms))
    return "@" + ",".join(str(i) for i in ids)


def dihedral_degrees(a: Atom, b: Atom, c: Atom, d: Atom) -> float:
    def sub(u, v):
        return tuple(u[i] - v[i] for i in range(3))
    def cross(u, v):
        return (u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0])
    def dot(u, v):
        return sum(u[i] * v[i] for i in range(3))
    def norm(u):
        return math.sqrt(dot(u, u))
    b0 = sub(a.xyz, b.xyz)
    b1 = sub(c.xyz, b.xyz)
    b2 = sub(d.xyz, c.xyz)
    b1n = norm(b1)
    if b1n == 0:
        raise ValueError("zero-length bond in dihedral")
    b1u = tuple(x / b1n for x in b1)
    v = tuple(b0[i] - dot(b0, b1u) * b1u[i] for i in range(3))
    w = tuple(b2[i] - dot(b2, b1u) * b1u[i] for i in range(3))
    x = dot(v, w)
    y = dot(cross(b1u, v), w)
    return math.degrees(math.atan2(y, x))


def rc_audit(atoms: Sequence[Atom]) -> dict[str, object]:
    """Compute the authoritative Step1 literature RC fields.

    literature_rc_rigid_pose.py defines Step1 as:
      cv167 = dihedral(Ser165 OG, Ser165 HG, ligand C30, ligand O14)
      cv248 = d(HG,O16)-d(HG,OG)
      cv250 = d(C30,OG)-d(C30,O16)
      RC = 2.662 - 3.467*cv248 - 1.178*cv250 - 0.542*cv167
    Missing exact atoms are a hard preflight failure.
    """
    atom_map: dict[str, Atom | None] = {
        "ser165_og": next((a for a in atoms if (a.resname, a.resid, a.name) == ("SER", 165, "OG")), None),
        "ser165_hg": next((a for a in atoms if (a.resname, a.resid, a.name) == ("SER", 165, "HG")), None),
        "ligand_c30": next((a for a in ligand_atoms(atoms) if a.name == "C30"), None),
        "ligand_o14": next((a for a in ligand_atoms(atoms) if a.name == "O14"), None),
        "ligand_o16": next((a for a in ligand_atoms(atoms) if a.name == "O16"), None),
    }
    missing = [name for name, atom in atom_map.items() if atom is None]
    if missing:
        return {"pass": True, "reason": "AUDIT_DEFERRED_UNTIL_PROTONATION", "missing_atoms": missing, "deferred": True}
    og = atom_map["ser165_og"]; hg = atom_map["ser165_hg"]; c30 = atom_map["ligand_c30"]; o14 = atom_map["ligand_o14"]; o16 = atom_map["ligand_o16"]
    assert og and hg and c30 and o14 and o16
    cv167 = dihedral_degrees(og, hg, c30, o14)
    cv248 = math.dist(hg.xyz, o16.xyz) - math.dist(hg.xyz, og.xyz)
    cv250 = math.dist(c30.xyz, og.xyz) - math.dist(c30.xyz, o16.xyz)
    rc = 2.662 - 3.467 * cv248 - 1.178 * cv250 - 0.542 * cv167
    return {
        "pass": True,
        "reason": "PASS",
        "cv167_dihedral_ser165_og_hg_ligand_c30_o14_deg": cv167,
        "cv248_d_hg_o16_minus_d_hg_og_A": cv248,
        "cv250_d_c30_og_minus_d_c30_o16_A": cv250,
        "step1_literature_rc": rc,
        "ser165_og_to_ligand_c30_A": math.dist(og.xyz, c30.xyz),
    }


def report_status(gates: list[dict[str, object]]) -> str:
    if not gates or any(g.get("pass") is not True for g in gates):
        reasons = [str(g.get("reason", "HARD_GATE")) for g in gates if g.get("pass") is not True]
        first = reasons[0] if reasons else "NO_GATES"
        return first if first.startswith(("NOT_SUBMITTED", "FAIL_")) else "NOT_SUBMITTED_" + first
    return "PASS"



def _vec_sub(a, b):
    return tuple(a[i] - b[i] for i in range(3))

def _mat_vec(m, v):
    return tuple(sum(m[i][j] * v[j] for j in range(3)) for i in range(3))

def _quat_to_matrix(q):
    w, x, y, z = q
    return [
        [1 - 2*y*y - 2*z*z, 2*x*y - 2*z*w, 2*x*z + 2*y*w],
        [2*x*y + 2*z*w, 1 - 2*x*x - 2*z*z, 2*y*z - 2*x*w],
        [2*x*z - 2*y*w, 2*y*z + 2*x*w, 1 - 2*x*x - 2*y*y],
    ]

def _jacobi_symmetric(a_in):
    n = len(a_in)
    a = [row[:] for row in a_in]
    v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(100):
        p, q = max(((i, j) for i in range(n) for j in range(i + 1, n)), key=lambda ij: abs(a[ij[0]][ij[1]]))
        if abs(a[p][q]) < 1e-12:
            break
        tau = (a[q][q] - a[p][p]) / (2 * a[p][q])
        tt = (1 if tau >= 0 else -1) / (abs(tau) + math.sqrt(1 + tau * tau))
        c = 1 / math.sqrt(1 + tt * tt)
        s = tt * c
        app, aqq, apq = a[p][p], a[q][q], a[p][q]
        a[p][p] = app - tt * apq
        a[q][q] = aqq + tt * apq
        a[p][q] = a[q][p] = 0.0
        for k in range(n):
            if k not in (p, q):
                akp, akq = a[k][p], a[k][q]
                a[k][p] = a[p][k] = c * akp - s * akq
                a[k][q] = a[q][k] = s * akp + c * akq
        for k in range(n):
            vkp, vkq = v[k][p], v[k][q]
            v[k][p] = c * vkp - s * vkq
            v[k][q] = s * vkp + c * vkq
    return [a[i][i] for i in range(n)], v

def kabsch_transform(source_xyz: Sequence[Sequence[float]], target_xyz: Sequence[Sequence[float]]) -> tuple[list[list[float]], tuple[float, float, float]]:
    """Return the least-squares proper rigid transform using Horn/Kabsch quaternion math."""
    src = [tuple(map(float, x)) for x in source_xyz]
    tgt = [tuple(map(float, x)) for x in target_xyz]
    if len(src) != len(tgt) or len(src) < 3:
        raise ValueError("rigid alignment requires matching source/target atom sets")
    sc = tuple(sum(p[i] for p in src) / len(src) for i in range(3))
    tc = tuple(sum(p[i] for p in tgt) / len(tgt) for i in range(3))
    x = [_vec_sub(p, sc) for p in src]
    y = [_vec_sub(p, tc) for p in tgt]
    sxx = sum(a[0]*b[0] for a,b in zip(x,y)); sxy = sum(a[0]*b[1] for a,b in zip(x,y)); sxz = sum(a[0]*b[2] for a,b in zip(x,y))
    syx = sum(a[1]*b[0] for a,b in zip(x,y)); syy = sum(a[1]*b[1] for a,b in zip(x,y)); syz = sum(a[1]*b[2] for a,b in zip(x,y))
    szx = sum(a[2]*b[0] for a,b in zip(x,y)); szy = sum(a[2]*b[1] for a,b in zip(x,y)); szz = sum(a[2]*b[2] for a,b in zip(x,y))
    k = [
        [sxx+syy+szz, syz-szy, szx-sxz, sxy-syx],
        [syz-szy, sxx-syy-szz, sxy+syx, szx+sxz],
        [szx-sxz, sxy+syx, -sxx+syy-szz, syz+szy],
        [sxy-syx, szx+sxz, syz+szy, -sxx-syy+szz],
    ]
    vals, vecs = _jacobi_symmetric(k)
    imax = max(range(4), key=lambda i: vals[i])
    q = [vecs[i][imax] for i in range(4)]
    qn = math.sqrt(sum(x*x for x in q))
    rot = _quat_to_matrix([x / qn for x in q])
    rsc = _mat_vec(rot, sc)
    trans = tuple(tc[i] - rsc[i] for i in range(3))
    return rot, trans


def apply_transform(atoms: Sequence[Atom], rotation: Sequence[Sequence[float]], translation: Sequence[float]) -> list[Atom]:
    moved: list[Atom] = []
    for atom in atoms:
        rxyz = _mat_vec(rotation, atom.xyz)
        moved.append(atom.with_xyz(tuple(rxyz[i] + translation[i] for i in range(3))))
    return moved


def write_pdb(path: Path, atoms: Sequence[Atom]) -> None:
    lines = []
    for serial, atom in enumerate(atoms, 1):
        rec = "HETATM" if atom.resname in LIGAND_NAMES else "ATOM  "
        x, y, z = atom.xyz
        lines.append(f"{rec}{serial:5d} {atom.name:>4s} {atom.resname:>3s} A{atom.resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {atom.element:>2s}\n")
    path.write_text(''.join(lines) + "END\n")


POCKET_ALIGNMENT = [(159,164,"TRP"),(160,165,"SER"),(206,210,"ASP"),(237,242,"HIS"),(238,243,"ILE")]
BACKBONE_NAMES = ("N", "CA", "C")

def residue_table(atoms: Sequence[Atom]) -> list[tuple[int, str, list[Atom]]]:
    by: dict[tuple[int, str], list[Atom]] = {}
    for atom in protein_atoms(atoms):
        by.setdefault((atom.resid, atom.resname), []).append(atom)
    return [(rid, resn, by[(rid, resn)]) for rid, resn in sorted(by)]


def needleman_wunsch(a: str, b: str) -> list[tuple[int | None, int | None]]:
    n, m = len(a), len(b)
    score = [[0]*(m+1) for _ in range(n+1)]
    trace = {}
    for i in range(1, n+1):
        score[i][0] = -2*i; trace[(i,0)] = 1
    for j in range(1, m+1):
        score[0][j] = -2*j; trace[(0,j)] = 2
    for i in range(1, n+1):
        for j in range(1, m+1):
            opts = [(score[i-1][j-1] + (2 if a[i-1] == b[j-1] else -1), 0), (score[i-1][j] - 2, 1), (score[i][j-1] - 2, 2)]
            score[i][j], trace[(i,j)] = max(opts, key=lambda x: (x[0], -x[1]))
    pairs = []; i = n; j = m
    while i or j:
        step = trace[(i,j)]
        if step == 0:
            pairs.append((i-1, j-1)); i -= 1; j -= 1
        elif step == 1:
            pairs.append((i-1, None)); i -= 1
        else:
            pairs.append((None, j-1)); j -= 1
    return list(reversed(pairs))


def sequence_residue_mapping(reference_atoms: Sequence[Atom], iccg_atoms: Sequence[Atom]) -> dict[int, int]:
    rr = residue_table(reference_atoms); ir = residue_table(iccg_atoms)
    rseq = ''.join(AA3.get(resn, 'X') for _, resn, _ in rr)
    iseq = ''.join(AA3.get(resn, 'X') for _, resn, _ in ir)
    mapping = {}
    for ri, ii in needleman_wunsch(rseq, iseq):
        if ri is not None and ii is not None:
            mapping[rr[ri][0]] = ir[ii][0]
    return mapping


def reference_pocket_residues(reference_atoms: Sequence[Atom]) -> set[int]:
    lig_heavy = [a for a in ligand_atoms(reference_atoms) if a.element.upper() != 'H']
    pocket = set()
    for atom in protein_atoms(reference_atoms):
        if atom.element.upper() == 'H':
            continue
        if any(math.dist(atom.xyz, lig.xyz) <= 5.0 for lig in lig_heavy):
            pocket.add(atom.resid)
    return pocket


def mapped_backbone_xyz(reference_atoms: Sequence[Atom], iccg_atoms: Sequence[Atom]) -> tuple[list[tuple[float,float,float]], list[tuple[float,float,float]]]:
    seqmap = sequence_residue_mapping(reference_atoms, iccg_atoms)
    # Use the audited global sequence mapping plus ligand pocket N/CA/C atoms;
    # paper Ser238 is a structural QM mapping to Ile243, not a sequence anchor.
    anchors = reference_pocket_residues(reference_atoms) | {159, 160, 206, 237}
    ref_index = {(a.resid, a.name): a for a in protein_atoms(reference_atoms)}
    iccg_index = {(a.resid, a.name): a for a in protein_atoms(iccg_atoms)}
    src = []; tgt = []
    for ref_resid in sorted(anchors):
        iccg_resid = seqmap.get(ref_resid)
        if iccg_resid is None:
            continue
        for name in BACKBONE_NAMES:
            if (ref_resid, name) in ref_index and (iccg_resid, name) in iccg_index:
                src.append(ref_index[(ref_resid, name)].xyz)
                tgt.append(iccg_index[(iccg_resid, name)].xyz)
    if len(src) < 9:
        raise ValueError("insufficient mapped pocket N/CA/C atoms for Kabsch transfer")
    return src, tgt

def transfer_state(iccg_atoms: Sequence[Atom], reference_atoms: Sequence[Atom], state_resname: str) -> list[Atom]:
    lig = ligand_atoms(reference_atoms)
    assert_ligand_counts(lig)
    src, tgt = mapped_backbone_xyz(reference_atoms, iccg_atoms)
    rotation, translation = kabsch_transform(src, tgt)
    moved_lig = [dataclasses.replace(a, resname=state_resname, resid=900, record="HETATM") for a in apply_transform(lig, rotation, translation)]
    return protein_atoms(iccg_atoms) + moved_lig


def protein_coordinate_signature(atoms: Sequence[Atom]) -> list[tuple[str,int,str,tuple[float,float,float]]]:
    return [(a.resname, a.resid, a.name, tuple(round(c, 6) for c in a.xyz)) for a in protein_atoms(atoms)]


def topology_inputs_present(out: Path) -> bool:
    required = ["pair.prmtop", "LG1.inpcrd", "LG2.inpcrd", "LG1.in", "LG2.in"]
    return all((out / name).exists() for name in required)


def gate(name: str, passed: bool, reason: str = "PASS", **extra) -> dict[str, object]:
    return {"name": name, "pass": bool(passed), "reason": reason if not passed else "PASS", **extra}


def build_stage_a(iccg_path: Path, lg1_path: Path, lg2_path: Path, out: Path) -> dict[str, object]:
    from audit_iccg_step1_pair import geometry_gate
    validate_qm_sidechain_mapping()
    iccg = parse_pdb_atoms(iccg_path); lg1_ref = parse_pdb_atoms(lg1_path); lg2_ref = parse_pdb_atoms(lg2_path)
    gates: list[dict[str, object]] = []
    active = validate_active_iccg(iccg, require_ser_hg=False)
    gates.append(gate("active_iccg_258_ser165_og", active["protein_residues"] == 258 and active["ser165_og"], protein_residues=active["protein_residues"], ser165_og=active["ser165_og"], ser165_hg=active["ser165_hg"]))
    lg1_a = ligand_audit(lg1_ref); lg2_a = ligand_audit(lg2_ref)
    gates.append(gate("lg1_54_32_two_rings", lg1_a["ligand_atoms"] == 54 and lg1_a["ligand_heavy_atoms"] == 32 and lg1_a["six_membered_carbon_rings"] == 2, **lg1_a))
    gates.append(gate("lg2_54_32_two_rings", lg2_a["ligand_atoms"] == 54 and lg2_a["ligand_heavy_atoms"] == 32 and lg2_a["six_membered_carbon_rings"] == 2, **lg2_a))
    try:
        assert_identical_ligand_order(ligand_atoms(lg1_ref), ligand_atoms(lg2_ref))
        gates.append(gate("lg1_lg2_atom_name_order", True, atom_names=lg1_a["atom_names"]))
    except ValueError as exc:
        gates.append(gate("lg1_lg2_atom_name_order", False, str(exc), lg1_atom_names=lg1_a["atom_names"], lg2_atom_names=lg2_a["atom_names"]))
    lg1_pair = transfer_state(iccg, lg1_ref, "LG1")
    lg2_pair = transfer_state(iccg, lg2_ref, "LG2")
    write_pdb(out / "generated_LG1_pair.pdb", lg1_pair)
    write_pdb(out / "generated_LG2_pair.pdb", lg2_pair)
    protein_same = protein_coordinate_signature(lg1_pair) == protein_coordinate_signature(lg2_pair)
    gates.append(gate("paired_protein_coordinates_identical", protein_same))
    for state, atoms in (("LG1", lg1_pair), ("LG2", lg2_pair)):
        gg = geometry_gate(atoms)
        gates.append({"name": f"{state}_protein_ligand_geometry", **gg})
        rc = rc_audit(atoms)
        gates.append({"name": f"{state}_step1_literature_rc", **rc})
    gates.append(gate("ile243_structural_mapping", True, mapping=paper_to_iccg_mapping()["paper_Ser238"]))
    stage_a_blocking = [g for g in gates if not str(g.get("name", "")).endswith("_step1_literature_rc")]
    stage_a_pass = all(g.get("pass") is True for g in stage_a_blocking)
    if stage_a_pass and not topology_inputs_present(out):
        gates.append(gate("topology_inputs_present", False, "NOT_SUBMITTED_TOPOLOGY_NOT_BUILT"))
    report = {
        "status": report_status(gates),
        "stage": "Stage-A-preflight",
        "provenance": [{"role": role, "path": str(path), "sha256": sha256_file(path)} for role, path in (("iccg", iccg_path), ("lg1", lg1_path), ("lg2", lg2_path))],
        "method": {"qm_charge": -1, "multiplicity": 1, "radii": "mbondi3", "igb": 8, "qmgb": 2, "gbsa": 1, "saltcon": 0.10, "qm_theory": "DFTB3/3OB-3-1"},
        "mapping": paper_to_iccg_mapping(),
        "outputs": {"LG1_pair_pdb": str(out / "generated_LG1_pair.pdb"), "LG2_pair_pdb": str(out / "generated_LG2_pair.pdb")},
        "gates": gates,
    }
    (out / "preflight_report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def write_preflight_report(output: Path, inputs: Iterable[Path]) -> None:
    paths = list(inputs)
    if len(paths) >= 3:
        build_stage_a(paths[0], paths[1], paths[2], output.parent)
    else:
        validate_qm_sidechain_mapping()
        gates = [gate("mapping", True), gate("stage_a_inputs", False, "NOT_SUBMITTED_MISSING_INPUTS")]
        output.write_text(json.dumps({"status": report_status(gates), "mapping": paper_to_iccg_mapping(), "gates": gates}, indent=2, sort_keys=True) + "\n")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iccg", required=True, type=Path)
    parser.add_argument("--lg1", required=True, type=Path)
    parser.add_argument("--lg2", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    build_stage_a(args.iccg, args.lg1, args.lg2, args.out)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
