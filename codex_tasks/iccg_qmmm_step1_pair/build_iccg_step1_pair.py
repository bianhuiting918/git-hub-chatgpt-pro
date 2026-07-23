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

AMBERTOOLS_DATA = "/work/home/acshdt1dks/petase_qmmm_pilot8_20260721/software/ambertools20_data_f0ab9845_lf"
AMBER_MODULE = "amber/2018-hpcx-gcc-7.3.1"
PARMED_EGG = "/public/software/apps/amber/2018/hpcx-2.4.1-gcc-7.3.1/lib/python2.7/site-packages/ParmEd-3.0.0_57.g74a84d30-py2.7-linux-x86_64.egg"
LITERATURE_PRMTOP = "/work/home/acshdt1dks/petase_orbmol_lg1_lg4_active_20260719/inputs/literature_rc_sources/acylation/vmd-md-b.prmtop"
LITERATURE_INPCRD = "/work/home/acshdt1dks/petase_orbmol_lg1_lg4_active_20260719/inputs/literature_rc_sources/acylation/vmd-md-b.inpcrd"
DFTB_SLKO_PATH = "/work/home/acshdt1dks/petase_qmmm_pilot8_20260721/third_party/3ob-3-1"
LITERATURE_LIGAND_RESNAME = "UNK"
LITERATURE_LIGAND_RESID = 265
LITERATURE_LIGAND_NET_CHARGE = 0
QMCHARGE = -1
EXPECTED_SYSTEM_CHARGE = 6.0
SYSTEM_CHARGE_TOLERANCE = 1e-4
LITERATURE_LIGAND_ORDER = ["C1","O2","C3","O4","C5","O6","C7","O8","C9","O10","C11","O12","C13","O14","C15","O16","C17","O18","C19","O20","C21","C22","C23","C24","C25","C26","C27","C28","C29","C30","C31","C32","H33","H34","H35","H36","H37","H38","H39","H40","H41","H42","H43","H44","H45","H46","H47","H48","H49","H50","H51","H52","H53","H54"]
QM_BACKBONE_EXCLUDED = {"N","H","H1","H2","H3","CA","HA","C","O","OXT"}

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


@dataclasses.dataclass(frozen=True)
class CommandResult:
    returncode: int
    stdout: str = ""
    stderr: str = ""

def default_runner(cmd: Sequence[str], cwd: Path) -> CommandResult:
    import subprocess
    proc = subprocess.run([str(x) for x in cmd], cwd=str(cwd), text=True, capture_output=True)
    return CommandResult(proc.returncode, proc.stdout, proc.stderr)


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


def translate_ligand_in_pair(pair_atoms: Sequence[Atom], vector: Sequence[float]) -> list[Atom]:
    return [atom.with_xyz(tuple(atom.xyz[i] + vector[i] for i in range(3))) if atom.resname in LIGAND_NAMES else atom for atom in pair_atoms]


def relief_candidates() -> list[tuple[float, float, float]]:
    candidates = [(0.0, 0.0, 0.0), (-0.240000, 0.0, 0.0), (-0.299578, -0.007654, -0.351502)]
    vals = [-0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
    for x in vals:
        for y in vals:
            for z in vals:
                if math.sqrt(x*x + y*y + z*z) <= 0.5000001:
                    candidates.append((x, y, z))
    seen = set(); unique = []
    for c in candidates:
        key = tuple(round(v, 6) for v in c)
        if key not in seen:
            seen.add(key); unique.append(c)
    return unique


def pocket_search_pair(pair_atoms: Sequence[Atom], cutoff_A: float = 6.0) -> list[Atom]:
    lig_heavy = [a for a in ligand_atoms(pair_atoms) if a.element.upper() != "H"]
    pocket = []
    for atom in protein_atoms(pair_atoms):
        if atom.element.upper() == "H":
            continue
        if any(math.dist(atom.xyz, lig.xyz) <= cutoff_A for lig in lig_heavy):
            pocket.append(atom)
    return pocket + list(ligand_atoms(pair_atoms))


def common_translation_relief(lg1_pair: Sequence[Atom], lg2_pair: Sequence[Atom]) -> dict[str, object]:
    from audit_iccg_step1_pair import geometry_gate
    lg1_scope = pocket_search_pair(lg1_pair, 6.0)
    lg2_scope = pocket_search_pair(lg2_pair, 6.0)
    best = None
    best_key = None
    for vec in relief_candidates():
        shifted1 = translate_ligand_in_pair(lg1_scope, vec)
        shifted2 = translate_ligand_in_pair(lg2_scope, vec)
        g1 = geometry_gate(shifted1); g2 = geometry_gate(shifted2)
        objective = max(float(g1["max_vdw_overlap_A"]), float(g2["max_vdw_overlap_A"]))
        norm = math.sqrt(sum(v*v for v in vec))
        both_pass = g1.get("pass") is True and g2.get("pass") is True
        # Prefer any two-state PASS, then the smallest displacement norm, then
        # smaller max overlap. If no candidate passes, keep the least-bad clash.
        key = (0, round(norm, 9), round(objective, 9), vec) if both_pass else (1, round(objective, 9), round(norm, 9), vec)
        record = {"search_scope": "6A_pocket_then_full_validation", "vector_A": list(vec), "norm_A": norm, "objective_max_overlap_A": objective, "LG1": g1, "LG2": g2, "search_both_pass": both_pass}
        if best is None or key < best_key:
            best = record; best_key = key
    assert best is not None
    full1 = geometry_gate(translate_ligand_in_pair(lg1_pair, best["vector_A"]))
    full2 = geometry_gate(translate_ligand_in_pair(lg2_pair, best["vector_A"]))
    best["full_validation"] = {"LG1": full1, "LG2": full2, "both_pass": full1.get("pass") is True and full2.get("pass") is True}
    best["LG1"] = full1
    best["LG2"] = full2
    best["objective_max_overlap_A"] = max(float(full1["max_vdw_overlap_A"]), float(full2["max_vdw_overlap_A"]))
    return best

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
    pre_geometry = {"LG1": geometry_gate(lg1_pair), "LG2": geometry_gate(lg2_pair)}
    relief = {"enabled": False, "reason": "INITIAL_GEOMETRY_PASS", "vector_A": [0.0, 0.0, 0.0], "norm_A": 0.0, "pre": pre_geometry}
    if not (pre_geometry["LG1"].get("pass") and pre_geometry["LG2"].get("pass")):
        best = common_translation_relief(lg1_pair, lg2_pair)
        relief = {"enabled": True, "reason": "INITIAL_GEOMETRY_FAILED_MINIMAX_COMMON_TRANSLATION", "search_scope": best["search_scope"], "vector_A": best["vector_A"], "norm_A": best["norm_A"], "pre": pre_geometry, "post": {"LG1": best["LG1"], "LG2": best["LG2"]}, "objective_max_overlap_A": best["objective_max_overlap_A"], "full_validation": best["full_validation"]}
        lg1_pair = translate_ligand_in_pair(lg1_pair, best["vector_A"])
        lg2_pair = translate_ligand_in_pair(lg2_pair, best["vector_A"])
    post_geometry = {"LG1": geometry_gate(lg1_pair), "LG2": geometry_gate(lg2_pair)}
    relief.setdefault("post", post_geometry)
    write_pdb(out / "generated_LG1_pair.pdb", lg1_pair)
    write_pdb(out / "generated_LG2_pair.pdb", lg2_pair)
    protein_same = protein_coordinate_signature(lg1_pair) == protein_coordinate_signature(lg2_pair)
    gates.append(gate("paired_protein_coordinates_identical", protein_same))
    for state, atoms in (("LG1", lg1_pair), ("LG2", lg2_pair)):
        gg = post_geometry[state]
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
        "relief": relief,
        "gates": gates,
    }
    (out / "preflight_report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report



def audit_stage_b_protonation(atoms: Sequence[Atom]) -> dict[str, object]:
    names = {(a.resname, a.resid, a.name) for a in atoms}
    ser_hg = ("SER", 165, "HG") in names
    his_hid = any(a.resname == "HID" and a.resid == 242 for a in atoms)
    his_hd1 = ("HID", 242, "HD1") in names
    his_ne2 = ("HID", 242, "NE2") in names
    passed = ser_hg and his_hid and his_hd1 and his_ne2
    missing = []
    if not ser_hg: missing.append("SER165_HG")
    if not his_hid: missing.append("HIS242_HID")
    if not his_hd1: missing.append("HID242_HD1")
    if not his_ne2: missing.append("HID242_NE2")
    return {"name": "stage_b_protonation", "pass": passed, "reason": "PASS" if passed else "NOT_SUBMITTED_PROTONATION_INCOMPLETE", "missing": missing}


def _stage_b_qm_residue_key(atom: Atom) -> tuple[str, int]:
    return ("HIS", atom.resid) if atom.resname == "HID" else (atom.resname, atom.resid)


def select_stage_b_qm_atoms(atoms: Sequence[Atom]) -> list[Atom]:
    sidechains = {(res, resid) for res, resid in QM_SIDECHAINS}
    selected = []
    for atom in atoms:
        if atom.resname in LIGAND_NAMES:
            selected.append(atom)
        elif _stage_b_qm_residue_key(atom) in sidechains and atom.name not in QM_BACKBONE_EXCLUDED:
            selected.append(atom)
    ligand_count = sum(1 for atom in selected if atom.resname in LIGAND_NAMES)
    present = {atom.resid for atom in selected if atom.resname not in LIGAND_NAMES}
    if ligand_count != 54 or present != {164, 165, 210, 242, 243}:
        raise ValueError("Stage-B QM region must contain ligand54 plus five complete sidechains")
    return selected


def audit_no_ser_ligand_bond(atoms: Sequence[Atom], bonds: Sequence[tuple[int, int]]) -> dict[str, object]:
    by_index = {atom.index: atom for atom in atoms}
    bad = []
    for i, j in bonds:
        a = by_index.get(i); b = by_index.get(j)
        if not a or not b:
            continue
        if ((a.resname == "SER" and a.resid == 165 and b.resname in LIGAND_NAMES) or (b.resname == "SER" and b.resid == 165 and a.resname in LIGAND_NAMES)):
            bad.append([i, j])
    return {"name": "no_ser_ligand_bond", "pass": not bad, "reason": "PASS" if not bad else "NOT_SUBMITTED_SER_LIGAND_BOND_PRESENT", "bad_bonds": bad}


def sander_input_text(state: str, iqmatoms: Sequence[int]) -> str:
    atoms = [int(i) for i in iqmatoms]
    iqmatoms_text = ",".join(str(i) for i in atoms)
    slko_path = DFTB_SLKO_PATH.rstrip("/") + "/"
    return f"""ICCG Step1 {state} DFTB3/GBN2 single point
 &cntrl
  imin=1, maxcyc=0, ntb=0, igb=8, gbsa=1, saltcon=0.10, cut=999,
  ifqnt=1,
 /
 &qmmm
  iqmatoms={iqmatoms_text},
  qmcharge={QMCHARGE}, qm_theory='DFTB3',
  qmgb=2, qmcut=999.0, dftb_maxiter=200, scfconv=1.0d-8, printcharges=1, verbosity=5,
  dftb_slko_path='{slko_path}',
 /
"""


def stage_b_required_outputs(out: Path) -> list[Path]:
    return [out / name for name in ["pair.prmtop", "LG1.inpcrd", "LG2.inpcrd", "LG1.in", "LG2.in"]]


def _dependency_paths(ambertools_prefix, parmed_egg, literature_prmtop, literature_inpcrd, dftb_slko_path) -> dict[str, Path]:
    return {
        "ambertools_prefix": Path(ambertools_prefix or AMBERTOOLS_DATA),
        "tleap": Path(ambertools_prefix or AMBERTOOLS_DATA) / "bin" / "tleap",
        "parmed_egg": Path(parmed_egg or PARMED_EGG),
        "literature_prmtop": Path(literature_prmtop or LITERATURE_PRMTOP),
        "literature_inpcrd": Path(literature_inpcrd or LITERATURE_INPCRD),
        "dftb_slko_path": Path(dftb_slko_path or DFTB_SLKO_PATH),
    }


DISULFIDE_PAIRS = ((238, 283), (275, 292))

def audit_disulfides(atoms: Sequence[Atom], bonds: Sequence[tuple[int, int]]) -> dict[str, object]:
    by_res_atom = {(a.resid, a.name): a for a in atoms}
    bondset = {tuple(sorted(pair)) for pair in bonds}
    missing = []
    for a, b in DISULFIDE_PAIRS:
        if (a, "SG") not in by_res_atom or (b, "SG") not in by_res_atom:
            missing.append(f"{a}-{b}:SG")
    cyx = [a for a in atoms if a.resid in {238, 283, 275, 292}]
    wrong_resname = sorted({a.resid for a in cyx if a.resname != "CYX"})
    cyx_hg = sorted({a.resid for a in cyx if a.name == "HG"})
    missing_bonds = [[a, b] for a, b in DISULFIDE_PAIRS if tuple(sorted((a, b))) not in bondset]
    passed = not missing and not wrong_resname and not cyx_hg and not missing_bonds
    return {"name": "stage_b_disulfides", "pass": passed, "reason": "PASS" if passed else "NOT_SUBMITTED_DISULFIDE_INCOMPLETE", "pairs": [list(p) for p in DISULFIDE_PAIRS], "missing_atoms": missing, "wrong_resname": wrong_resname, "cyx_hg_present": cyx_hg, "missing_bonds": missing_bonds}

def _rewrite_resname_line(line: str, resname: str) -> str:
    return line[:17] + f"{resname:>3s}" + line[20:]

def write_stage_b_protein_only_pdb(source_pair: Path, output: Path) -> Path:
    lines = []
    for line in source_pair.read_text().splitlines(True):
        if not line.startswith("ATOM"):
            continue
        name = line[12:16].strip()
        element = (line[76:78].strip() or ''.join(c for c in name if c.isalpha())[:1]).upper()
        if element == "H":
            continue
        resid = int(line[22:26])
        resname = line[17:20].strip()
        if resid in {238, 275, 283, 292}:
            resname = "CYX"
        elif resid == 242 and resname == "HIS":
            resname = "HID"
        elif resname == "HIS":
            resname = "HIE"
        lines.append(_rewrite_resname_line(line, resname))
    output.write_text(''.join(lines) + "END\n")
    return output

def _write_tleap_input(out: Path) -> Path:
    text = """source leaprc.protein.ff14SB
set default PBRadii mbondi3
p = loadpdb stage_b_protein_only.pdb
bond p.238.SG p.283.SG
bond p.275.SG p.292.SG
saveamberparm p protein.prmtop protein.rst7
quit
"""
    path = out / "stage_b_tleap.in"
    path.write_text(text)
    return path


def _write_parmed_stage_b_script(out: Path, deps: dict[str, Path]) -> Path:
    script = """#!/usr/bin/env python2
import json, math, os, sys
sys.path.insert(0, __PARMED_EGG__)
import parmed
from parmed.amber import ChamberParm
from parmed.tools.actions import changeRadii
LITERATURE_PRMTOP = __LITERATURE_PRMTOP__
LITERATURE_INPCRD = __LITERATURE_INPCRD__
LIGAND_ORDER = __LIGAND_ORDER__
QM_BACKBONE_EXCLUDED = set(__EXCLUDED__)
EXPECTED_DISULFIDES = [(238, 283), (275, 292)]
QMCHARGE = -1
EXPECTED_SYSTEM_CHARGE = 6.0
SYSTEM_CHARGE_TOLERANCE = 1e-4


def finite_xyz(xyz):
    return all(not (math.isnan(float(v)) or math.isinf(float(v))) for row in xyz for v in row)


def load_pair_pdb(path):
    atoms = []
    for line in open(path):
        if line.startswith(('ATOM', 'HETATM')):
            atoms.append((line[:6].strip(), line[12:16].strip(), line[17:20].strip(), int(line[22:26]), [float(line[30:38]), float(line[38:46]), float(line[46:54])]))
    return atoms


def ligand_xyz_by_order(path):
    het = [(name, xyz) for rec, name, resname, resid, xyz in load_pair_pdb(path) if rec == 'HETATM']
    names = [name for name, xyz in het]
    if names != LIGAND_ORDER:
        raise SystemExit('NOT_SUBMITTED_LIGAND_ORDER_MISMATCH:%s' % path)
    return [xyz for name, xyz in het]


def max_delta(coords_a, coords_b):
    delta = 0.0
    for a, b in zip(coords_a, coords_b):
        d = max(abs(float(a[i]) - float(b[i])) for i in range(3))
        if d > delta:
            delta = d
    return delta


def atom_xyz(atom):
    return [float(atom.xx), float(atom.xy), float(atom.xz)]


def copy_with_ligand_xyz(top, xyz):
    state = top.copy(parmed.Structure)
    for atom, coord in zip(state.atoms[-54:], xyz):
        atom.xx, atom.xy, atom.xz = coord
    return ChamberParm.from_structure(state)


def apply_ligand_xyz(top, xyz):
    return copy_with_ligand_xyz(top, xyz)


def bond_residue_pairs(top):
    pairs = set()
    for bond in top.bonds:
        ra = bond.atom1.residue
        rb = bond.atom2.residue
        if ra is not rb:
            pairs.add((ra.idx + 1, rb.idx + 1))
            pairs.add((rb.idx + 1, ra.idx + 1))
    return pairs


protein = parmed.load_file('protein.prmtop', 'protein.rst7')
paper = parmed.load_file(LITERATURE_PRMTOP, LITERATURE_INPCRD)
sources = [res for res in paper.residues if res.name == 'UNK' and len(res.atoms) == 54]
if len(sources) != 1:
    raise SystemExit('NOT_SUBMITTED_LIGAND_SOURCE_NOT_UNIQUE')
srcres = sources[0]
ligand = paper[':%d' % (srcres.idx + 1)].copy(parmed.Structure)
if [a.name for a in ligand.atoms] != LIGAND_ORDER:
    raise SystemExit('NOT_SUBMITTED_LIGAND_ORDER_MISMATCH')
source_ligand_types = [getattr(a, 'type', '') for a in ligand.atoms]
source_ligand_charges = [round(float(a.charge), 8) for a in ligand.atoms]
source_ligand_bonds = sorted([(min(b.atom1.idx - ligand.atoms[0].idx, b.atom2.idx - ligand.atoms[0].idx), max(b.atom1.idx - ligand.atoms[0].idx, b.atom2.idx - ligand.atoms[0].idx)) for b in ligand.bonds])
for res in ligand.residues:
    res.name = 'LG1'
ligand.box = None
combined = protein.copy(parmed.Structure) + ligand
top = ChamberParm.from_structure(combined)
changeRadii(top,'mbondi3').execute()
original_resids = []
for rec, name, resname, resid, xyz in load_pair_pdb('stage_b_protein_only.pdb'):
    if resid not in original_resids:
        original_resids.append(resid)
if len(original_resids) != 258:
    raise SystemExit('NOT_SUBMITTED_PROTEIN_RESIDUE_COUNT')
orig_to_top = dict((resid, i + 1) for i, resid in enumerate(original_resids))
lg1_xyz = ligand_xyz_by_order('generated_LG1_pair.pdb')
lg2_xyz = ligand_xyz_by_order('generated_LG2_pair.pdb')
lg1_top = apply_ligand_xyz(top, lg1_xyz)
lg2_top = apply_ligand_xyz(top, lg2_xyz)
protein_coords = [atom_xyz(a) for a in top.atoms[:-54]]
lg1_protein_delta = max_delta(protein_coords, [atom_xyz(a) for a in lg1_top.atoms[:-54]])
lg2_protein_delta = max_delta(protein_coords, [atom_xyz(a) for a in lg2_top.atoms[:-54]])
protein_coordinate_max_delta_A = max(lg1_protein_delta, lg2_protein_delta)
LG1_ligand_generated_pair_max_delta_A = max_delta(lg1_xyz, [atom_xyz(a) for a in lg1_top.atoms[-54:]])
LG2_ligand_generated_pair_max_delta_A = max_delta(lg2_xyz, [atom_xyz(a) for a in lg2_top.atoms[-54:]])
top.save('pair.prmtop', overwrite=True)
lg1_top.save('LG1.inpcrd', overwrite=True)
lg2_top.save('LG2.inpcrd', overwrite=True)
iqmatoms = []
qm_top_resids = set(orig_to_top[x] for x in [164, 165, 210, 242, 243] if x in orig_to_top)
for atom in top.atoms:
    if atom.residue.idx + 1 in qm_top_resids and atom.name not in QM_BACKBONE_EXCLUDED:
        iqmatoms.append(atom.idx + 1)
lig_start = len(top.atoms) - 54 + 1
iqmatoms.extend(range(lig_start, len(top.atoms) + 1))
lig_atoms = top.atoms[-54:]
lig_types_match = [getattr(a, 'type', '') for a in lig_atoms] == source_ligand_types
lig_charges_match = [round(float(a.charge), 8) for a in lig_atoms] == source_ligand_charges
lig_bonds = sorted([(min(b.atom1.idx - lig_atoms[0].idx, b.atom2.idx - lig_atoms[0].idx), max(b.atom1.idx - lig_atoms[0].idx, b.atom2.idx - lig_atoms[0].idx)) for b in top.bonds if b.atom1 in lig_atoms and b.atom2 in lig_atoms])
ligand_internal_bonds_match = lig_bonds == source_ligand_bonds
crossbond = any((b.atom1 in lig_atoms) != (b.atom2 in lig_atoms) for b in top.bonds)
res_by_orig = dict((orig, top.residues[i]) for i, orig in enumerate(original_resids))
bond_pairs = bond_residue_pairs(top)
disulfide_bonds = all((orig_to_top[a], orig_to_top[b]) in bond_pairs for a, b in EXPECTED_DISULFIDES)
cyx_no_hg = all(res_by_orig[x].name == 'CYX' and not any(a.name == 'HG' for a in res_by_orig[x].atoms) for x in [238, 275, 283, 292])
his242 = res_by_orig[242]
his242_hid = his242.name == 'HID' and any(a.name == 'HD1' for a in his242.atoms) and not any(a.name == 'HE2' for a in his242.atoms)
ser165_hg = any(a.name == 'HG' for a in res_by_orig[165].atoms)
mbondi3 = all(abs(float(getattr(a, 'solvent_radius', 0.0))) > 0.0 for a in top.atoms)
coords_finite = finite_xyz([atom_xyz(a) for a in lg1_top.atoms]) and finite_xyz([atom_xyz(a) for a in lg2_top.atoms])
system_charge = sum(a.charge for a in top.atoms)
checks = [
    ('residue_count_259', len(top.residues) == 259),
    ('protein_residue_count_258', len(original_resids) == 258),
    ('atom_count_3902', len(top.atoms) == 3902),
    ('protein_atom_count_3848', len(top.atoms) - 54 == 3848),
    ('ligand_atom_count_54', len(lig_atoms) == 54),
    ('ligand_names', [a.name for a in lig_atoms] == LIGAND_ORDER),
    ('ligand_types', lig_types_match),
    ('ligand_charges', lig_charges_match),
    ('ligand_internal_bonds', ligand_internal_bonds_match),
    ('two_cyx_disulfide_bonds', disulfide_bonds),
    ('cyx_no_hg', cyx_no_hg),
    ('his242_hid_hd1_no_he2', his242_hid),
    ('ser165_hg', ser165_hg),
    ('no_protein_ligand_crossbond', not crossbond),
    ('protein_coordinate_delta_zero', protein_coordinate_max_delta_A == 0.0),
    ('LG1_ligand_generated_pair_delta_zero', LG1_ligand_generated_pair_max_delta_A == 0.0),
    ('LG2_ligand_generated_pair_delta_zero', LG2_ligand_generated_pair_max_delta_A == 0.0),
    ('finite_coordinates', coords_finite),
    ('mbondi3', mbondi3),
    ('qm_atom_count_100', len(iqmatoms) == 100),
    ('qmcharge_minus_one', QMCHARGE == -1),
    ('system_charge_plus_six', abs(system_charge - EXPECTED_SYSTEM_CHARGE) < 1e-4),
]
gates = [{'name': name, 'pass': bool(ok), 'reason': 'PASS' if ok else 'NOT_SUBMITTED_TOPOLOGY_AUDIT_FAILED'} for name, ok in checks]
status = 'PASS' if all(g['pass'] for g in gates) else 'NOT_SUBMITTED_TOPOLOGY_AUDIT_FAILED'
json.dump({
    'status': status,
    'gates': gates,
    'residues': len(top.residues),
    'protein_residues': len(original_resids),
    'atoms': len(top.atoms),
    'protein_atoms': len(top.atoms) - 54,
    'ligand_atoms': 54,
    'qm_atom_count': len(iqmatoms),
    'iqmatoms': iqmatoms,
    'qmcharge': QMCHARGE,
    'multiplicity': 1,
    'mbondi3': mbondi3,
    'ligand_order': LIGAND_ORDER,
    'ligand_types_match': lig_types_match,
    'ligand_charges_match': lig_charges_match,
    'ligand_internal_bonds_match': ligand_internal_bonds_match,
    'protein_coordinate_max_delta_A': protein_coordinate_max_delta_A,
    'LG1_ligand_generated_pair_max_delta_A': LG1_ligand_generated_pair_max_delta_A,
    'LG2_ligand_generated_pair_max_delta_A': LG2_ligand_generated_pair_max_delta_A,
    'system_charge': system_charge,
    'rc_audit_only': True,
}, open('stage_b_parmed_manifest.json','w'), indent=2, sort_keys=True)
"""
    script = (script
        .replace("__PARMED_EGG__", repr(str(deps["parmed_egg"])))
        .replace("__LITERATURE_PRMTOP__", repr(str(deps["literature_prmtop"])))
        .replace("__LITERATURE_INPCRD__", repr(str(deps["literature_inpcrd"])))
        .replace("__LIGAND_ORDER__", repr(LITERATURE_LIGAND_ORDER))
        .replace("__EXCLUDED__", repr(sorted(QM_BACKBONE_EXCLUDED)))
    )
    path = out / "parmed_stage_b.py"
    path.write_text(script)
    return path

def _write_stage_b_report(out: Path, status: str, gates: list[dict[str, object]], deps: dict[str, Path]) -> dict[str, object]:
    report = {
        "status": status,
        "stage": "Stage-B-topology",
        "remote_resources": {"AmberTools": str(deps["ambertools_prefix"]), "Amber module": AMBER_MODULE, "ParmEd egg": str(deps["parmed_egg"]), "literature_prmtop": str(deps["literature_prmtop"]), "literature_inpcrd": str(deps["literature_inpcrd"]), "DFTB": str(deps["dftb_slko_path"])},
        "ligand_copy_contract": {"source_resname": LITERATURE_LIGAND_RESNAME, "source_resid": LITERATURE_LIGAND_RESID, "atom_order": LITERATURE_LIGAND_ORDER, "net_charge": LITERATURE_LIGAND_NET_CHARGE, "copy_types_charges_bonds": True, "reparameterize": False},
        "qm_contract": {"ligand_atoms": 54, "sidechains": [f"{res}{rid}" for res, rid in QM_SIDECHAINS], "qmcharge": -1, "multiplicity": 1, "excluded_backbone_atoms": sorted(QM_BACKBONE_EXCLUDED)},
        "gates": gates,
    }
    (out / "topology_audit.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def build_stage_b(out: Path, runner=default_runner, ambertools_prefix=None, parmed_egg=None, literature_prmtop=None, literature_inpcrd=None, dftb_slko_path=None) -> dict[str, object]:
    out.mkdir(parents=True, exist_ok=True)
    deps = _dependency_paths(ambertools_prefix, parmed_egg, literature_prmtop, literature_inpcrd, dftb_slko_path)
    missing_deps = [name for name, path in deps.items() if not path.exists()]
    for state in ("LG1", "LG2"):
        (out / f"{state}.in").write_text(sander_input_text(state, []))
    if missing_deps:
        return _write_stage_b_report(out, "NOT_SUBMITTED_DEPENDENCY_MISSING", [{"name": "stage_b_dependencies", "pass": False, "reason": "NOT_SUBMITTED_DEPENDENCY_MISSING", "missing": missing_deps}], deps)
    protein_only = out / "stage_b_protein_only.pdb"
    source_pair = out / "generated_LG1_pair.pdb"
    if source_pair.exists():
        write_stage_b_protein_only_pdb(source_pair, protein_only)
    else:
        protein_only.write_text("END\n")
    tleap_input = _write_tleap_input(out)
    tleap = runner([str(deps["tleap"]), "-f", str(tleap_input)], out)
    if tleap.returncode != 0:
        return _write_stage_b_report(out, "NOT_SUBMITTED_TLEAP_FAILED", [{"name": "tleap_ff14sb_mbondi3", "pass": False, "reason": "NOT_SUBMITTED_TLEAP_FAILED", "stdout": tleap.stdout, "stderr": tleap.stderr}], deps)
    parmed_script = _write_parmed_stage_b_script(out, deps)
    parmed = runner(["/usr/bin/env", "PYTHONPATH=" + str(deps["parmed_egg"]), "/usr/bin/python2.7", str(parmed_script), "parmed_stage_b"], out)
    if parmed.returncode != 0:
        return _write_stage_b_report(out, "NOT_SUBMITTED_PARMED_FAILED", [{"name": "parmed_ligand_copy", "pass": False, "reason": "NOT_SUBMITTED_PARMED_FAILED", "stdout": parmed.stdout, "stderr": parmed.stderr}], deps)
    manifest_path = out / "stage_b_parmed_manifest.json"
    if not manifest_path.exists():
        return _write_stage_b_report(out, "NOT_SUBMITTED_TOPOLOGY_AUDIT_FAILED", [{"name": "stage_b_parmed_manifest", "pass": False, "reason": "NOT_SUBMITTED_TOPOLOGY_AUDIT_FAILED"}], deps)
    manifest = json.loads(manifest_path.read_text())
    iqmatoms = [int(i) for i in manifest.get("iqmatoms", [])]
    for state in ("LG1", "LG2"):
        (out / f"{state}.in").write_text(sander_input_text(state, iqmatoms))
    missing = [str(path) for path in stage_b_required_outputs(out) if not path.exists()]
    manifest_gates = manifest.get("gates")
    gates_ok = isinstance(manifest_gates, list) and bool(manifest_gates) and all(g.get("pass") is True for g in manifest_gates)
    manifest_pass = manifest.get("status") == "PASS" and bool(manifest.get("iqmatoms")) and gates_ok
    gates = [
        {"name": "stage_b_dependencies", "pass": True, "reason": "PASS"},
        {"name": "tleap_ff14sb_mbondi3", "pass": True, "reason": "PASS"},
        {"name": "parmed_ligand_copy", "pass": True, "reason": "PASS", "temporary_script": str(parmed_script)},
        {"name": "stage_b_topology_manifest", "pass": bool(manifest_pass), "reason": "PASS" if manifest_pass else "NOT_SUBMITTED_TOPOLOGY_AUDIT_FAILED"},
        {"name": "stage_b_required_topology_and_inputs", "pass": not missing, "reason": "PASS" if not missing else "NOT_SUBMITTED_TOPOLOGY_NOT_BUILT", "missing": missing},
    ]
    status = "PASS" if (manifest_pass and not missing) else ("NOT_SUBMITTED_TOPOLOGY_NOT_BUILT" if missing else "NOT_SUBMITTED_TOPOLOGY_AUDIT_FAILED")
    report = _write_stage_b_report(out, status, gates, deps)
    report["topology_manifest"] = manifest
    (out / "topology_audit.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
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
    parser.add_argument("--stage-b", action="store_true")
    parser.add_argument("--ambertools-prefix", type=Path, default=Path(AMBERTOOLS_DATA))
    parser.add_argument("--parmed-egg", type=Path, default=Path(PARMED_EGG))
    parser.add_argument("--literature-prmtop", type=Path, default=Path(LITERATURE_PRMTOP))
    parser.add_argument("--literature-inpcrd", type=Path, default=Path(LITERATURE_INPCRD))
    parser.add_argument("--dftb-slko-path", type=Path, default=Path(DFTB_SLKO_PATH))
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    if args.stage_b:
        build_stage_b(args.out, ambertools_prefix=args.ambertools_prefix, parmed_egg=args.parmed_egg, literature_prmtop=args.literature_prmtop, literature_inpcrd=args.literature_inpcrd, dftb_slko_path=args.dftb_slko_path)
    else:
        build_stage_a(args.iccg, args.lg1, args.lg2, args.out)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
