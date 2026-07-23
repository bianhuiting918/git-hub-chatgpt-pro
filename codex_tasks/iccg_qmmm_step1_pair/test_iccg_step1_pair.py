#!/usr/bin/env python3
"""Unit/integration tests for the ICCG Step1 LG1/LG2 pair gate."""
import json, tempfile, unittest
from pathlib import Path
import audit_iccg_step1_pair as audit
import build_iccg_step1_pair as build
import collect_iccg_step1_pair as collect

FIXTURES = Path(__file__).resolve().parents[1] / "iccg_lg1_dimer" / "inputs"
ICCG_HEAVY = FIXTURES / "ICCG_active_chainA_heavy.pdb"
LG1_REF = FIXTURES / "IsPETase_WT_LG1_reference.pdb"

class PairGateTests(unittest.TestCase):
    def test_real_sibling_fixtures_validate_258_54_32_two_rings_and_ile243_mapping(self):
        iccg = build.parse_pdb_atoms(ICCG_HEAVY)
        lg1 = build.parse_pdb_atoms(LG1_REF)
        self.assertEqual(build.validate_active_iccg(iccg)["protein_residues"], 258)
        self.assertFalse(build.validate_active_iccg(iccg)["ser165_hg"], "heavy-only fixture should not require HG before protonation")
        la = build.ligand_audit(lg1)
        self.assertEqual(la["ligand_resname"], "LG1")
        self.assertEqual(la["ligand_atoms"], 54)
        self.assertEqual(la["ligand_heavy_atoms"], 32)
        self.assertEqual(la["six_membered_carbon_rings"], 2)
        mapping = build.paper_to_iccg_mapping()["paper_Ser238"]
        self.assertEqual((mapping["iccg_resname"], mapping["iccg_resid"]), ("ILE", 243))
        self.assertAlmostEqual(mapping["ca_separation_A"], 0.939)

    def test_structural_mapping_rejects_sequence_gap_null_mapping(self):
        with self.assertRaises(ValueError):
            build.validate_qm_sidechain_mapping({"paper_Ser238": None})

    def test_strict_lg1_lg2_atom_order_and_topology_pairing(self):
        lg1 = [build.Atom(i, "LG1", 999, f"C{i}", (float(i), 0.0, 0.0), "C", "HETATM") for i in range(1, 55)]
        lg2 = [build.Atom(i, "LG2", 999, f"C{i}", (float(i), 1.0, 0.0), "C", "HETATM") for i in range(1, 55)]
        build.assert_identical_ligand_order(lg1, lg2)
        bad_count = lg2[:-1]
        with self.assertRaises(ValueError):
            build.assert_identical_ligand_order(lg1, bad_count)
        bad_order = list(reversed(lg2))
        with self.assertRaises(ValueError):
            build.assert_identical_ligand_order(lg1, bad_order)
        self.assertTrue(audit.topology_pairing_pass(["N", "CA", "L1"], ["N", "CA", "L1"]))
        self.assertFalse(audit.topology_pairing_pass([], []))
        self.assertFalse(audit.topology_pairing_pass(["N", "CA", "L1"], ["N", "CA", "L2"]))

    def test_lg1_lg2_five_residue_full_ligand_qm_selection(self):
        atoms = []
        idx = 1
        for resid, resname in [(164,"TRP"),(165,"SER"),(210,"ASP"),(242,"HIS"),(243,"ILE")]:
            for name in ["CA", "CB", "X"]:
                atoms.append(build.Atom(idx, resname, resid, name, (idx,0,0), "C", "ATOM")); idx += 1
        for n in range(54):
            atoms.append(build.Atom(idx, "LG2", 900, f"L{n+1}", (idx,0,0), "C", "HETATM")); idx += 1
        qm = build.select_qm_atoms(atoms)
        self.assertEqual(len([a for a in qm if a.resname == "LG2"]), 54)
        self.assertEqual({a.resid for a in qm if a.resname != "LG2"}, {164,165,210,242,243})

    def test_geometry_gate_only_protein_ligand_and_excludes_bonds(self):
        p = build.Atom(1, "ALA", 1, "C", (0,0,0), "C", "ATOM")
        p2 = build.Atom(2, "ALA", 1, "O", (9.0,0,0), "O", "ATOM")
        l = build.Atom(3, "LG1", 9, "C1", (0.9,0,0), "C", "HETATM")
        self.assertTrue(audit.geometry_gate([p, p2, l], bonds={(1,3)})["pass"])
        failed = audit.geometry_gate([p, p2, l], bonds=set())
        self.assertFalse(failed["pass"])
        self.assertEqual(failed["reason"], "FAIL_GEOMETRY_CLASH_NOT_LABEL")

    def test_preflight_report_name_and_fail_closed_submission(self):
        with tempfile.TemporaryDirectory() as td:
            out = Path(td)
            build.write_preflight_report(out / "preflight_report.json", [ICCG_HEAVY, LG1_REF])
            self.assertTrue((out / "preflight_report.json").exists())
            self.assertFalse(audit.can_submit(out / "preflight_report.json"))
            empty = out / "empty_gates.json"
            empty.write_text(json.dumps({"status": "PASS", "gates": []}))
            self.assertFalse(audit.can_submit(empty))
            missing = out / "missing_gates.json"
            missing.write_text(json.dumps({"status": "PASS"}))
            self.assertFalse(audit.can_submit(missing))
            ok = out / "ok.json"
            required = [
                "active_iccg_258_ser165_og",
                "lg1_54_32_two_rings",
                "lg2_54_32_two_rings",
                "lg1_lg2_atom_name_order",
                "paired_protein_coordinates_identical",
                "LG1_protein_ligand_geometry",
                "LG2_protein_ligand_geometry",
                "ile243_structural_mapping",
            ]
            ok.write_text(json.dumps({"status": "NOT_SUBMITTED_TOPOLOGY_NOT_BUILT", "gates": [{"name": name, "pass": True} for name in required] + [{"name": "topology_inputs_present", "pass": False}]}))
            self.assertTrue(audit.can_submit(ok))

    def test_rigid_transform_moves_ligand_not_protein(self):
        src = [(0.,0.,0.),(1.,0.,0.),(0.,1.,0.)]
        tgt = [(2.,3.,4.),(3.,3.,4.),(2.,4.,4.)]
        r, t = build.kabsch_transform(src, tgt)
        moved = build.apply_transform([build.Atom(1,"LG1",1,"C1",src[0],"C","HETATM")], r, t)[0]
        self.assertEqual(moved.xyz, tgt[0])

    def test_step1_rc_requires_exact_literature_atoms_and_hg(self):
        atoms = [
            build.Atom(1, "SER", 165, "OG", (0.0, 0.0, 0.0), "O", "ATOM"),
            build.Atom(2, "LG2", 900, "C30", (4.57, 0.0, 0.0), "C", "HETATM"),
            build.Atom(3, "LG2", 900, "O14", (4.57, 1.0, 0.0), "O", "HETATM"),
            build.Atom(4, "LG2", 900, "O16", (3.00, 0.0, 0.0), "O", "HETATM"),
        ]
        missing = build.rc_audit(atoms)
        self.assertTrue(missing["pass"])
        self.assertEqual(missing["reason"], "AUDIT_DEFERRED_UNTIL_PROTONATION")
        self.assertIn("ser165_hg", missing["missing_atoms"])
        atoms.append(build.Atom(5, "SER", 165, "HG", (0.96, 0.0, 0.0), "H", "ATOM"))
        rc = build.rc_audit(atoms)
        self.assertTrue(rc["pass"])
        self.assertIn("cv167_dihedral_ser165_og_hg_ligand_c30_o14_deg", rc)
        self.assertIn("cv248_d_hg_o16_minus_d_hg_og_A", rc)
        self.assertIn("cv250_d_c30_og_minus_d_c30_o16_A", rc)
        self.assertAlmostEqual(rc["ser165_og_to_ligand_c30_A"], 4.57)

    def test_lg2_expected_overlap_above_threshold_blocks_submission(self):
        protein = build.Atom(1, "SER", 165, "OG", (0.0, 0.0, 0.0), "O", "ATOM")
        ligand = build.Atom(2, "LG2", 900, "C30", (2.221, 0.0, 0.0), "C", "HETATM")
        gate = audit.geometry_gate([protein, ligand])
        self.assertFalse(gate["pass"])
        self.assertEqual(gate["reason"], "FAIL_GEOMETRY_CLASH_NOT_LABEL")
        self.assertGreater(gate["max_vdw_overlap_A"], 0.80)

    def test_rigid_transfer_preserves_literature_ligand_internal_distances(self):
        iccg = build.parse_pdb_atoms(ICCG_HEAVY)
        ref = build.parse_pdb_atoms(LG1_REF)
        before = build.ligand_atoms(ref)
        after = build.ligand_atoms(build.transfer_state(iccg, ref, "LG1"))
        self.assertEqual([a.name for a in before], [a.name for a in after])
        for i, atom_i in enumerate(before):
            for j, atom_j in enumerate(before[i + 1:], i + 1):
                d0 = ((atom_i.xyz[0]-atom_j.xyz[0])**2 + (atom_i.xyz[1]-atom_j.xyz[1])**2 + (atom_i.xyz[2]-atom_j.xyz[2])**2) ** 0.5
                d1 = ((after[i].xyz[0]-after[j].xyz[0])**2 + (after[i].xyz[1]-after[j].xyz[1])**2 + (after[i].xyz[2]-after[j].xyz[2])**2) ** 0.5
                self.assertAlmostEqual(d0, d1, places=6)

    def test_stage_a_main_transfers_lg1_pair_and_fails_closed_with_real_gates(self):
        with tempfile.TemporaryDirectory() as td:
            out = Path(td)
            rc = build.main(["--iccg", str(ICCG_HEAVY), "--lg1", str(LG1_REF), "--lg2", str(LG1_REF), "--out", str(out)])
            self.assertEqual(rc, 0)
            report = json.loads((out / "preflight_report.json").read_text())
            self.assertTrue((out / "generated_LG1_pair.pdb").stat().st_size > 0)
            self.assertTrue((out / "generated_LG2_pair.pdb").stat().st_size > 0)
            self.assertTrue(report["gates"], "Stage-A must produce non-empty gates")
            self.assertEqual(report["status"], "NOT_SUBMITTED_TOPOLOGY_NOT_BUILT")
            self.assertTrue(audit.can_submit(out / "preflight_report.json"), "Stage-A structural gates should pass despite deferred RC and missing topology")
            lg1_geo = next(g for g in report["gates"] if g["name"] == "LG1_protein_ligand_geometry")
            self.assertTrue(lg1_geo["pass"])
            self.assertAlmostEqual(lg1_geo["min_nonbonded_A"], 2.669, places=3)
            self.assertLessEqual(lg1_geo["max_vdw_overlap_A"], 0.80)
            lg1_rc = next(g for g in report["gates"] if g["name"] == "LG1_step1_literature_rc")
            self.assertTrue(lg1_rc["pass"])
            self.assertEqual(lg1_rc["reason"], "AUDIT_DEFERRED_UNTIL_PROTONATION")
            self.assertIn("ser165_hg", lg1_rc["missing_atoms"])

    def test_pass_json_not_written_from_exit_zero_only(self):
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "PASS.json"
            log = Path(td) / "sander.out"
            log.write_text("Amber exited with code 0\n")
            with self.assertRaises(ValueError):
                collect.write_pass_json(log, out, process_exit=0, geometry_pass=True)
            self.assertFalse(out.exists())

if __name__ == "__main__":
    unittest.main(verbosity=2)
