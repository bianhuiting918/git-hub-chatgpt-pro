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
        self.assertTrue(model.has_bond("C", "CAP_N"))
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


if __name__ == "__main__":
    unittest.main()
