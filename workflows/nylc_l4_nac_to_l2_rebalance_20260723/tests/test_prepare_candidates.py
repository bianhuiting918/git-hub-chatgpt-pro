#!/usr/bin/env python3
import sys
import tempfile
import unittest
from types import SimpleNamespace
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE / "scripts"))

from prepare_candidates import (  # noqa: E402
    make_posre,
    reaction_geometry_from_target_globals,
    rewrite_topology,
)


class PrepareCandidateTests(unittest.TestCase):
    def test_rewrite_nylc_topology_replaces_only_ligand_contract(self):
        source = """#include "amber99sb-ildn.ff/forcefield.itp"
#include "l4_nfree_GMX.itp"
#include "topol_Protein_chain_A.itp"
[ molecules ]
Protein_chain_A 1
l4_nfree 1
SOL 10
"""
        rewritten = rewrite_topology(source, "l4_nfree_GMX.itp", "l4_nfree")
        self.assertIn('#include "PA66_L2_GMX.itp"', rewritten)
        self.assertNotIn("l4_nfree_GMX.itp", rewritten)
        self.assertIn("PA66_L2 1", rewritten)
        self.assertNotIn("\nl4_nfree 1", rewritten)
        self.assertEqual(rewritten.count("POSRES_L2_1000"), 1)

    def test_rewrite_nyl12_removes_old_ligand_restraint_blocks(self):
        source = """#include "amber99sb-ildn.ff/forcefield.itp"
#include "pa66_l4_GMX.itp"
#ifdef POSRES_LIG
#include "posre_P66.itp"
#endif
#ifdef POSRES_LIG_WEAK
#include "posre_P66_weak.itp"
#endif
#include "protein_Protein_chain_A.itp"
[ molecules ]
Protein_chain_A 1
pa66_l4 1
"""
        rewritten = rewrite_topology(source, "pa66_l4_GMX.itp", "pa66_l4")
        self.assertNotIn("POSRES_LIG", rewritten)
        self.assertNotIn("posre_P66", rewritten)
        self.assertIn('#include "posre_l2_10.itp"', rewritten)
        self.assertIn("PA66_L2 1", rewritten)

    def test_posre_is_exactly_33_heavy_atoms(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "posre.itp"
            make_posre(path, 500)
            rows = [
                line for line in path.read_text(encoding="utf-8").splitlines()
                if line.strip() and not line.lstrip().startswith((";", "["))
            ]
            self.assertEqual(len(rows), 33)
            self.assertEqual([int(row.split()[0]) for row in rows], list(range(1, 34)))
            self.assertTrue(all(row.split()[2:] == ["500", "500", "500"] for row in rows))

    def test_nac_angle_uses_carbonyl_oxygen_not_amide_nitrogen(self):
        atoms = [
            SimpleNamespace(coordinate=(0.0, 0.0, 0.0)),
            SimpleNamespace(coordinate=(1.0, 0.0, 0.0)),
            SimpleNamespace(coordinate=(-1.0, 0.0, 0.0)),
            SimpleNamespace(coordinate=(0.0, 1.0, 0.0)),
        ]
        system = SimpleNamespace(atoms=atoms, box=(5.0, 5.0, 5.0))
        geometry = reaction_geometry_from_target_globals(
            system,
            {
                "carbonyl_c": 1,
                "carbonyl_o": 2,
                "amide_n": 3,
                "thr_og1": 4,
            },
        )
        self.assertAlmostEqual(geometry["distance_thr_og1_to_carbonyl_c_nm"], 1.0)
        self.assertAlmostEqual(
            geometry["angle_thr_og1_carbonyl_c_carbonyl_o_deg"], 90.0
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
