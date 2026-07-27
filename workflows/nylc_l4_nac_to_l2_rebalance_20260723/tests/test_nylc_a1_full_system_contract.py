#!/usr/bin/env python3
import importlib.util
import re
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
BUILDER = HERE / "scripts" / "build_nylc_a1_full_system.py"


def load_builder():
    spec = importlib.util.spec_from_file_location("a1_full", BUILDER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def atoms_from_itp(text):
    block = re.split(r"(?mi)^\s*\[\s*atoms\s*\]\s*$", text, maxsplit=1)[1]
    block = re.split(r"(?m)^\s*\[", block, maxsplit=1)[0]
    rows = {}
    for line in block.splitlines():
        data = line.split(";", 1)[0].split()
        if len(data) >= 8 and data[0].isdigit():
            rows[int(data[0])] = {
                "type": data[1], "resnr": int(data[2]), "name": data[4],
                "charge": float(data[6]),
            }
    return rows


def interaction_set(text, section, arity):
    result = set()
    parts = re.split(rf"(?mi)^\s*\[\s*{section}\s*\]\s*$", text)
    for block in parts[1:]:
        block = re.split(r"(?m)^\s*\[", block, maxsplit=1)[0]
        for line in block.splitlines():
            data = line.split(";", 1)[0].split()
            if len(data) >= arity + 1 and all(x.isdigit() for x in data[:arity]):
                result.add(tuple(int(x) for x in data[:arity]))
    return result


PARENT_ATOMS = [
("N3","N",-0.62540,14.01),("H","H1",0.1934,1.008),("H","H2",0.1934,1.008),
("CT","CA",0.0034,12.01),("HP","HA",0.1087,1.008),("CT","CB",0.4514,12.01),
("H1","HB",-0.0323,1.008),("CT","CG2",-0.2554,12.01),("HC","HG21",0.0627,1.008),
("HC","HG22",0.0627,1.008),("HC","HG23",0.0627,1.008),("OH","OG1",-0.6764,16.0),
("HO","HG1",0.407,1.008),("C","C",0.6163,12.01),("O","O",-0.5722,16.0),
("N","N",-0.4157,14.01),("H","H",0.2719,1.008),("CT","CA",-0.0389,12.01),
("H1","HA",0.1007,1.008),
]


def parent_itp():
    atom_lines = []
    for i, (typ, name, charge, mass) in enumerate(PARENT_ATOMS, 1):
        resnr = 267 if i <= 15 else 268
        atom_lines.append(f"{i} {typ} {resnr} THR {name} {i} {charge:.7f} {mass}")
    return """[ moleculetype ]
Protein_chain_H 3
[ atoms ]
""" + "\n".join(atom_lines) + """
[ bonds ]
1 2 1
1 3 1
1 4 1
4 5 1
4 6 1
4 14 1
6 7 1
6 8 1
6 12 1
12 13 1
14 15 1
14 16 1
16 17 1
16 18 1
18 19 1
[ pairs ]
1 7 1
1 8 1
1 12 1
1 15 1
1 16 1
2 5 1
2 6 1
2 14 1
3 5 1
3 6 1
3 14 1
4 13 1
7 13 1
8 13 1
12 14 1
[ angles ]
2 1 3 1
2 1 4 1
3 1 4 1
1 4 5 1
1 4 6 1
1 4 14 1
4 6 12 1
7 6 12 1
8 6 12 1
6 12 13 1
[ dihedrals ]
2 1 4 5 9
2 1 4 6 9
2 1 4 14 9
3 1 4 5 9
3 1 4 6 9
3 1 4 14 9
1 4 6 12 9
4 6 12 13 9
7 6 12 13 9
8 6 12 13 9
[ dihedrals ]
4 16 14 15 4
"""


MODEL_ATOMS = [
("nz","N",-1.106599,14.01),("hn","H1",0.3848,1.008),("hn","H2",0.3848,1.008),
("hn","HG1",0.3848,1.008),("c3","CA",0.0655,12.01),("hx","HA",0.0807,1.008),
("c3","CB",0.2468,12.01),("h1","HB",0.0497,1.008),("c3","CG2",-0.1011,12.01),
("hc","HG21",0.049033,1.008),("hc","HG22",0.049033,1.008),("hc","HG23",0.049033,1.008),
("o","OG1",-0.4895,16.0),("c","C",0.6131,12.01),("o","O",-0.6161,16.0),
("ns","NCAP",-0.5649,14.01),("hn","HCAP",0.3105,1.008),("c3","CCAP",0.0733,12.01),
("h1","HC1",0.0457,1.008),
]


def model_top():
    atom_lines = [
        f"{i} {typ} 1 TA1 {name} {i} {charge:.7f} {mass}"
        for i, (typ, name, charge, mass) in enumerate(MODEL_ATOMS, 1)
    ]
    return """[ atoms ]
""" + "\n".join(atom_lines) + """
[ bonds ]
1 5 1 0.14843 199836.208
7 13 1 0.13206 339782.640
1 2 1 0.10271 371815.344
1 3 1 0.10271 371815.344
1 4 1 0.10271 371815.344
[ pairs ]
1 8 1
1 9 1
1 13 1
1 15 1
1 16 1
2 6 1
2 7 1
2 14 1
3 6 1
3 7 1
3 14 1
4 6 1
4 7 1
4 14 1
6 13 1
8 13 1
9 13 1
10 13 1
11 13 1
12 13 1
13 14 1
[ angles ]
2 1 3 1 108.54 393.7144
2 1 4 1 108.54 393.7144
2 1 5 1 111.18 460.82576
3 1 4 1 108.54 393.7144
3 1 5 1 111.18 460.82576
4 1 5 1 111.18 460.82576
5 7 13 1 112.84 650.86304
8 7 13 1 116.75 485.59504
9 7 13 1 112.84 650.86304
[ dihedrals ]
1 5 7 13 1 0.0 0.6508444 3
2 1 5 6 1 0.0 0.4560560 3
2 1 5 7 1 0.0 0.6508444 3
2 1 5 14 1 0.0 0.6508444 3
3 1 5 6 1 0.0 0.4560560 3
3 1 5 7 1 0.0 0.6508444 3
3 1 5 14 1 0.0 0.6508444 3
4 1 5 6 1 0.0 0.4560560 3
4 1 5 7 1 0.0 0.6508444 3
4 1 5 14 1 0.0 0.6508444 3
6 5 7 13 1 0.0 0.6508444 3
13 7 5 14 1 0.0 0.6508444 3
"""


def gro_line(resid, residue, atom, atom_id, x, y, z):
    return f"{resid:5d}{residue:<5}{atom:>5}{atom_id:5d}{x:8.3f}{y:8.3f}{z:8.3f}"


def gro_fixture():
    rows = [gro_line(267, "THR", f"X{i}", i, i / 100.0, 0.0, 0.0) for i in range(1, 20)]
    return "fixture\n19\n" + "\n".join(rows) + "\n5.0 5.0 5.0\n"


MODEL_MOL2 = """@<TRIPOS>ATOM
1 N 0.000 0.000 0.000 nz 1 TA1 0.0
2 H1 0.000 1.000 0.000 hn 1 TA1 0.0
3 H2 0.000 0.000 1.000 hn 1 TA1 0.0
4 HG1 1.000 2.000 3.000 hn 1 TA1 0.0
@<TRIPOS>BOND
"""


@unittest.skipUnless(BUILDER.is_file(), "full-system A1 builder not implemented")
class A1FullSystemContractTests(unittest.TestCase):
    def test_patch_preserves_atom_count_and_integer_chain_charge(self):
        module = load_builder()
        patched, audit = module.patch_chain_itp(parent_itp(), model_top())
        self.assertEqual(audit["atom_count_before"], audit["atom_count_after"])
        self.assertAlmostEqual(audit["chain_charge_before_e"], audit["chain_charge_after_e"], places=6)
        self.assertAlmostEqual(audit["thr267_charge_after_e"], 0.0, places=6)
        self.assertAlmostEqual(audit["raw_a1_thr267_charge_e"], 0.044, places=6)
        self.assertEqual(audit["status"], "PASS_A1_CHAIN_PATCH")
        atoms = atoms_from_itp(patched)
        self.assertEqual((atoms[1]["type"], atoms[12]["type"], atoms[13]["type"]), ("nz", "o", "hn"))

    def test_patch_replaces_old_og_hg1_graph_with_n_hg1_graph(self):
        module = load_builder()
        patched, audit = module.patch_chain_itp(parent_itp(), model_top())
        bonds = {tuple(sorted(x)) for x in interaction_set(patched, "bonds", 2)}
        self.assertIn((1, 13), bonds)
        self.assertNotIn((12, 13), bonds)
        pairs = {tuple(sorted(x)) for x in interaction_set(patched, "pairs", 2)}
        for pair in [(5, 13), (6, 13), (13, 14)]:
            self.assertIn(pair, pairs)
        for pair in [(4, 13), (7, 13), (8, 13)]:
            self.assertNotIn(pair, pairs)
        angles = interaction_set(patched, "angles", 3)
        for angle in [(2, 1, 13), (3, 1, 13), (13, 1, 4)]:
            self.assertIn(angle, angles)
        self.assertNotIn((6, 12, 13), angles)
        dihedrals = interaction_set(patched, "dihedrals", 4)
        for dih in [(13, 1, 4, 5), (13, 1, 4, 6), (13, 1, 4, 14)]:
            self.assertIn(dih, dihedrals)
        self.assertNotIn((4, 6, 12, 13), dihedrals)
        self.assertEqual(audit["removed_bond"], [12, 13])
        self.assertEqual(audit["added_bond"], [1, 13])

    def test_coordinate_patch_changes_only_hg1(self):
        module = load_builder()
        patched, audit = module.patch_gro_coordinates(gro_fixture(), MODEL_MOL2, 1)
        before = module.parse_gro(gro_fixture())
        after = module.parse_gro(patched)
        changed = [
            i for i in before["atoms"]
            if before["atoms"][i]["xyz_nm"] != after["atoms"][i]["xyz_nm"]
        ]
        self.assertEqual(changed, [13])
        self.assertEqual(after["atoms"][13]["xyz_nm"], [0.1, 0.2, 0.3])
        self.assertEqual(audit["status"], "PASS_A1_COORDINATE_PATCH")


class A1FullSystemBuilderPresenceTest(unittest.TestCase):
    def test_full_system_builder_exists(self):
        self.assertTrue(BUILDER.is_file(), str(BUILDER))


if __name__ == "__main__":
    unittest.main()
