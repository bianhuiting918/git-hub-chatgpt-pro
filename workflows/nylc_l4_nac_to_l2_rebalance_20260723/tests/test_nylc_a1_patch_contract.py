#!/usr/bin/env python3
import importlib.util
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
SCRIPT = HERE / "scripts" / "prepare_nylc_a1_parameter_model.py"


def load_module():
    spec = importlib.util.spec_from_file_location("prepare_nylc_a1_parameter_model", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parent_fixture():
    coords = {
        "N": (0.000, 0.000, 0.000),
        "H1": (-0.080, 0.060, 0.000),
        "H2": (-0.080, -0.060, 0.000),
        "CA": (0.145, 0.000, 0.000),
        "HA": (0.175, 0.095, 0.000),
        "CB": (0.205, -0.075, 0.120),
        "HB": (0.165, -0.175, 0.120),
        "CG2": (0.365, -0.070, 0.120),
        "HG21": (0.400, 0.030, 0.120),
        "HG22": (0.400, -0.120, 0.205),
        "HG23": (0.400, -0.120, 0.035),
        "OG1": (0.155, -0.035, 0.250),
        "HG1": (0.080, -0.030, 0.285),
        "C": (0.205, 0.040, -0.130),
        "O": (0.160, 0.130, -0.190),
    }
    thr_atoms = {
        name: {"source_atom_id": i + 1, "coordinate_nm": xyz}
        for i, (name, xyz) in enumerate(coords.items())
    }
    next_atoms = {
        "N": {"source_atom_id": 16, "coordinate_nm": (0.315, -0.025, -0.170)},
        "H": {"source_atom_id": 17, "coordinate_nm": (0.350, -0.095, -0.120)},
        "CA": {"source_atom_id": 18, "coordinate_nm": (0.375, 0.005, -0.300)},
        "HA": {"source_atom_id": 19, "coordinate_nm": (0.330, 0.090, -0.350)},
    }
    bonds = {
        frozenset(pair)
        for pair in [
            ("N", "H1"), ("N", "H2"), ("N", "CA"),
            ("CA", "HA"), ("CA", "CB"), ("CA", "C"),
            ("CB", "HB"), ("CB", "CG2"), ("CB", "OG1"),
            ("CG2", "HG21"), ("CG2", "HG22"), ("CG2", "HG23"),
            ("OG1", "HG1"), ("C", "O"), ("C", "next:N"),
            ("next:N", "next:H"), ("next:N", "next:CA"), ("next:CA", "next:HA"),
        ]
    }
    return thr_atoms, next_atoms, bonds


def gromacs_text_fixture():
    thr, nxt, bonds = parent_fixture()
    atom_rows = []
    local_names = list(thr) + list(nxt)
    for local_id, name in enumerate(local_names, 1):
        is_next = local_id > len(thr)
        record = nxt[name] if is_next else thr[name]
        resnr = 268 if is_next else 267
        resname = "THR"
        atom_rows.append(
            f"{local_id:5d} CT {resnr:5d} {resname:4s} {name:5s} {local_id:5d} 0.0 1.0"
        )
    thr_name_to_local = {name: i + 1 for i, name in enumerate(thr)}
    next_name_to_local = {
        name: len(thr) + i + 1 for i, name in enumerate(nxt)
    }
    bond_rows = []
    for pair in bonds:
        left, right = tuple(pair)
        left_name = left.split(":", 1)[-1]
        right_name = right.split(":", 1)[-1]
        if left.startswith("next:"):
            left_id = next_name_to_local[left_name]
        else:
            left_id = thr_name_to_local[left_name]
        if right.startswith("next:"):
            right_id = next_name_to_local[right_name]
        else:
            right_id = thr_name_to_local[right_name]
        bond_rows.append(f"{left_id:5d} {right_id:5d} 1")
    itp = "[ atoms ]\n" + "\n".join(atom_rows) + "\n\n[ bonds ]\n" + "\n".join(bond_rows) + "\n"

    gro_rows = []
    all_records = [(267, "THR", name, rec) for name, rec in thr.items()]
    all_records += [(268, "THR", name, rec) for name, rec in nxt.items()]
    for atom_id, (resnr, resname, name, rec) in enumerate(all_records, 1):
        x, y, z = rec["coordinate_nm"]
        gro_rows.append(
            f"{resnr:5d}{resname:<5s}{name:>5s}{atom_id:5d}{x:8.3f}{y:8.3f}{z:8.3f}"
        )
    gro = "fixture\n" + str(len(gro_rows)) + "\n" + "\n".join(gro_rows) + "\n   4.0 4.0 4.0\n"
    return itp, gro


class A1ParameterModelContractTests(unittest.TestCase):
    def test_builder_transfers_existing_hg1_without_changing_thr_atom_count(self):
        module = load_module()
        thr, nxt, bonds = parent_fixture()
        model = module.build_a1_capped_model(thr, nxt, bonds)
        thr_model_atoms = [atom for atom in model.atoms if atom.source_role == "Thr267"]
        self.assertEqual(len(thr_model_atoms), len(thr))
        self.assertEqual(model.atom("HG1").source_atom_id, thr["HG1"]["source_atom_id"])
        self.assertNotEqual(model.atom("HG1").coordinate_nm, thr["HG1"]["coordinate_nm"])

    def test_builder_creates_nh3_alkoxide_bond_graph_and_neutral_formal_charge(self):
        module = load_module()
        thr, nxt, bonds = parent_fixture()
        model = module.build_a1_capped_model(thr, nxt, bonds)
        self.assertTrue(model.has_bond("N", "HG1"))
        self.assertFalse(model.has_bond("OG1", "HG1"))
        self.assertTrue(model.has_bond("C", "NCAP"))
        self.assertEqual(model.atom("N").formal_charge, 1)
        self.assertEqual(model.atom("OG1").formal_charge, -1)
        self.assertEqual(sum(atom.formal_charge for atom in model.atoms), 0)
        self.assertEqual(len(model.atoms), 21)
        self.assertEqual(len(model.bonds), 20)

    def test_builder_rejects_parent_without_og1_hg1_bond(self):
        module = load_module()
        thr, nxt, bonds = parent_fixture()
        bonds.remove(frozenset(("OG1", "HG1")))
        with self.assertRaisesRegex(module.ModelError, "OG1-HG1"):
            module.build_a1_capped_model(thr, nxt, bonds)

    def test_all_mol2_atom_names_are_element_safe_and_four_chars_or_less(self):
        module = load_module()
        thr, nxt, bonds = parent_fixture()
        model = module.build_a1_capped_model(thr, nxt, bonds)
        for atom in model.atoms:
            self.assertLessEqual(len(atom.name), 4)
            self.assertEqual(atom.name[0], atom.element[0])

    def test_mol2_records_formal_ion_pair_and_exact_atom_map(self):
        module = load_module()
        thr, nxt, bonds = parent_fixture()
        model = module.build_a1_capped_model(thr, nxt, bonds)
        mol2 = module.render_mol2(model)
        self.assertIn("@<TRIPOS>MOLECULE\nNTA1_CAP\n21 20", mol2)
        self.assertIn("N.4", mol2)
        self.assertIn("O.3", mol2)
        self.assertEqual(
            {row["model_atom"]: row["source_atom_id"] for row in model.atom_map if row["source_atom_id"] is not None}["HG1"],
            thr["HG1"]["source_atom_id"],
        )


    def test_gromacs_parser_uses_bonded_graph_and_global_offset(self):
        module = load_module()
        itp, gro = gromacs_text_fixture()
        model = module.model_from_gromacs_text(
            itp, gro, chain_first_global_atom=1, thr_resnr=267, next_resnr=268
        )
        self.assertEqual(model.atom("N").source_atom_id, 1)
        self.assertEqual(model.atom("HG1").source_atom_id, 13)
        self.assertTrue(model.has_bond("N", "HG1"))
        self.assertFalse(model.has_bond("OG1", "HG1"))

    def test_gromacs_parser_rejects_wrong_global_atom_identity(self):
        module = load_module()
        itp, gro = gromacs_text_fixture()
        lines = gro.splitlines()
        lines[2] = lines[2][:10] + f"{'XX':>5s}" + lines[2][15:]
        bad = "\n".join(lines) + "\n"
        with self.assertRaisesRegex(module.ModelError, "identity mismatch"):
            module.model_from_gromacs_text(
                itp, bad, chain_first_global_atom=1, thr_resnr=267, next_resnr=268
            )


if __name__ == "__main__":
    unittest.main()
