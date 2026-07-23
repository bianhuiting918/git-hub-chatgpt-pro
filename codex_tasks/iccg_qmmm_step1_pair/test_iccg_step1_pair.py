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

    def test_geometry_gate_scans_all_pairs_and_reports_true_worst_pair(self):
        atoms = [
            build.Atom(1, "ALA", 1, "CB", (0.0, 0.0, 0.0), "C", "ATOM"),
            build.Atom(2, "TYR", 95, "CB", (10.0, 0.0, 0.0), "C", "ATOM"),
            build.Atom(3, "LG2", 900, "O14", (2.50, 0.0, 0.0), "O", "HETATM"),
            build.Atom(4, "LG2", 900, "C30", (11.0, 0.0, 0.0), "C", "HETATM"),
        ]
        gate = audit.geometry_gate(atoms)
        self.assertFalse(gate["pass"])
        self.assertEqual(gate["worst_pair"]["protein_residue"], 95)
        self.assertEqual(gate["worst_pair"]["ligand_atom"], "C30")
        self.assertGreater(gate["max_vdw_overlap_A"], 2.0)

    def test_common_translation_relief_uses_one_vector_and_preserves_internal_distances(self):
        protein = [build.Atom(1, "ALA", 1, "CB", (0.0, 0.0, 0.0), "C", "ATOM")]
        lig1 = [build.Atom(2, "LG1", 900, "C1", (2.35, 0.0, 0.0), "C", "HETATM"), build.Atom(3, "LG1", 900, "C2", (3.85, 0.0, 0.0), "C", "HETATM")]
        lig2 = [build.Atom(4, "LG2", 900, "C1", (2.35, 0.0, 0.0), "C", "HETATM"), build.Atom(5, "LG2", 900, "C2", (3.85, 0.0, 0.0), "C", "HETATM")]
        best = build.common_translation_relief(protein + lig1, protein + lig2)
        self.assertLessEqual(best["norm_A"], 0.5000001)
        self.assertTrue(best["LG1"]["pass"])
        self.assertTrue(best["LG2"]["pass"])
        shifted = build.translate_ligand_in_pair(protein + lig1, best["vector_A"])
        before = ((lig1[0].xyz[0]-lig1[1].xyz[0])**2 + (lig1[0].xyz[1]-lig1[1].xyz[1])**2 + (lig1[0].xyz[2]-lig1[1].xyz[2])**2) ** 0.5
        after_lig = build.ligand_atoms(shifted)
        after = ((after_lig[0].xyz[0]-after_lig[1].xyz[0])**2 + (after_lig[0].xyz[1]-after_lig[1].xyz[1])**2 + (after_lig[0].xyz[2]-after_lig[1].xyz[2])**2) ** 0.5
        self.assertAlmostEqual(before, after, places=6)

    def test_common_translation_relief_prefers_minimum_pass_norm_over_lower_overlap(self):
        protein = [build.Atom(1, "ALA", 1, "CB", (0.0, 0.0, 0.0), "C", "ATOM")]
        lig1 = [build.Atom(2, "LG1", 900, "C1", (-2.361, 0.0, 0.0), "C", "HETATM")]
        lig2 = [build.Atom(3, "LG2", 900, "C1", (-2.361, 0.0, 0.0), "C", "HETATM")]
        best = build.common_translation_relief(protein + lig1, protein + lig2)
        self.assertEqual(best["search_scope"], "6A_pocket_then_full_validation")
        self.assertEqual([round(v, 6) for v in best["vector_A"]], [-0.24, 0.0, 0.0])
        self.assertAlmostEqual(best["norm_A"], 0.24, places=6)
        self.assertTrue(best["full_validation"]["both_pass"])
        self.assertLessEqual(best["objective_max_overlap_A"], 0.80)

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

class StageBTests(unittest.TestCase):
    def test_stage_b_literature_ligand_copy_contract(self):
        self.assertEqual(build.LITERATURE_LIGAND_RESNAME, "UNK")
        self.assertEqual(build.LITERATURE_LIGAND_RESID, 265)
        self.assertEqual(build.LITERATURE_LIGAND_NET_CHARGE, 0)
        self.assertEqual(build.LITERATURE_LIGAND_ORDER[:4], ["C1", "O2", "C3", "O4"])
        self.assertEqual(build.LITERATURE_LIGAND_ORDER[-1], "H54")
        self.assertEqual(len(build.LITERATURE_LIGAND_ORDER), 54)

    def test_stage_b_his242_hid_and_ser165_hg_required(self):
        atoms = [
            build.Atom(1, "SER", 165, "OG", (0,0,0), "O", "ATOM"),
            build.Atom(2, "HID", 242, "HD1", (0,0,0), "H", "ATOM"),
            build.Atom(3, "HID", 242, "NE2", (0,0,0), "N", "ATOM"),
        ]
        audit_b = build.audit_stage_b_protonation(atoms)
        self.assertFalse(audit_b["pass"])
        self.assertEqual(audit_b["reason"], "NOT_SUBMITTED_PROTONATION_INCOMPLETE")
        atoms.append(build.Atom(4, "SER", 165, "HG", (0,0,0), "H", "ATOM"))
        self.assertTrue(build.audit_stage_b_protonation(atoms)["pass"])

    def test_stage_b_qm_sidechains_exclude_backbone_and_include_ligand54(self):
        atoms = []
        idx = 1
        for resid, resname in [(164,"TRP"),(165,"SER"),(210,"ASP"),(242,"HID"),(243,"ILE")]:
            for name in ["N","H","CA","HA","C","O","CB","SC"]:
                atoms.append(build.Atom(idx, resname, resid, name, (0,0,0), "C", "ATOM")); idx += 1
        for name in build.LITERATURE_LIGAND_ORDER:
            atoms.append(build.Atom(idx, "LG1", 900, name, (0,0,0), "C" if not name.startswith("H") else "H", "HETATM")); idx += 1
        qm = build.select_stage_b_qm_atoms(atoms)
        self.assertEqual(sum(1 for a in qm if a.resname == "LG1"), 54)
        self.assertFalse(any(a.name in build.QM_BACKBONE_EXCLUDED for a in qm if a.resname != "LG1"))
        self.assertEqual({a.resid for a in qm if a.resname != "LG1"}, {164,165,210,242,243})

    def test_stage_b_no_ser_ligand_covalent_bond_gate(self):
        bonds = [(1, 2), (10, 11)]
        atoms = [build.Atom(1, "SER", 165, "OG", (0,0,0), "O", "ATOM"), build.Atom(2, "LG1", 900, "C30", (0,0,0), "C", "HETATM")]
        gate = build.audit_no_ser_ligand_bond(atoms, bonds)
        self.assertFalse(gate["pass"])
        self.assertEqual(gate["reason"], "NOT_SUBMITTED_SER_LIGAND_BOND_PRESENT")
        self.assertTrue(build.audit_no_ser_ligand_bond(atoms, [(10, 11)])["pass"])

    def test_stage_b_sander_inputs_and_missing_topology_fail_closed(self):
        with tempfile.TemporaryDirectory() as td:
            out = Path(td)
            deps = {}
            for key, name in {"ambertools_prefix":"ambertools", "parmed_egg":"ParmEd.egg", "literature_prmtop":"vmd-md-b.prmtop", "literature_inpcrd":"vmd-md-b.inpcrd", "dftb_slko_path":"3ob-3-1"}.items():
                path = out / name
                if key.endswith("path") or key == "ambertools_prefix": path.mkdir()
                else: path.write_text("fixture")
                deps[key] = path
            (deps["ambertools_prefix"] / "bin").mkdir(); (deps["ambertools_prefix"] / "bin" / "tleap").write_text("tleap")
            def fake_runner(cmd, cwd):
                joined = " ".join(str(x) for x in cmd)
                if "tleap" in joined:
                    (cwd / "protein.prmtop").write_text("protein topology"); (cwd / "protein.rst7").write_text("protein coords")
                if "parmed_stage_b" in joined:
                    (cwd / "stage_b_parmed_manifest.json").write_text(json.dumps({"status":"PASS", "iqmatoms":[1,2,3], "qm_atom_count":100, "residues":259, "atoms":3902, "mbondi3":True, "gates":[{"name":"dynamic","pass":True}]}))
                return build.CommandResult(0, "ok", "")
            report = build.build_stage_b(out, runner=fake_runner, **deps)
            self.assertEqual(report["status"], "NOT_SUBMITTED_TOPOLOGY_NOT_BUILT")
            self.assertFalse(report["gates"][-1]["pass"])
            text = build.sander_input_text("LG1", "@1,2,3")
            for needle in ["imin=1", "maxcyc=0", "igb=8", "ifqnt=1", "qmgb=2", "dftb_maxiter=200", "scfconv=1d-8", build.DFTB_SLKO_PATH]:
                self.assertIn(needle, text)

    def test_stage_b_invokes_tleap_and_parmed_runner_to_build_outputs(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deps = {}
            for key, name in {"ambertools_prefix":"ambertools", "parmed_egg":"ParmEd.egg", "literature_prmtop":"vmd-md-b.prmtop", "literature_inpcrd":"vmd-md-b.inpcrd", "dftb_slko_path":"3ob-3-1"}.items():
                path = root / name
                if key.endswith("path") or key == "ambertools_prefix":
                    path.mkdir()
                else:
                    path.write_text("fixture")
                deps[key] = path
            (deps["ambertools_prefix"] / "bin").mkdir()
            (deps["ambertools_prefix"] / "bin" / "tleap").write_text("#!/bin/sh\n")
            for state in ("LG1", "LG2"):
                (root / f"generated_{state}_pair.pdb").write_text("END\n")
            calls = []
            def fake_runner(cmd, cwd):
                calls.append(cmd)
                joined = " ".join(str(x) for x in cmd)
                if "tleap" in joined:
                    (cwd / "protein.prmtop").write_text("protein topology")
                    (cwd / "protein.rst7").write_text("protein coords")
                if "parmed_stage_b" in joined:
                    (cwd / "pair.prmtop").write_text("pair topology")
                    (cwd / "LG1.inpcrd").write_text("lg1 coords")
                    (cwd / "LG2.inpcrd").write_text("lg2 coords")
                    (cwd / "stage_b_parmed_manifest.json").write_text(json.dumps({"status":"PASS", "iqmatoms":[1,2,3], "qm_atom_count":100, "residues":259, "atoms":3902, "mbondi3":True, "gates":[{"name":"dynamic","pass":True}]}))
                return build.CommandResult(0, "ok", "")
            report = build.build_stage_b(root, runner=fake_runner, **deps)
            self.assertEqual(report["status"], "PASS")
            self.assertTrue(any("tleap" in " ".join(map(str, c)) for c in calls))
            self.assertTrue(any("parmed_stage_b" in " ".join(map(str, c)) for c in calls))
            self.assertTrue((root / "topology_audit.json").exists())
            self.assertTrue((root / "LG1.in").exists())
            self.assertTrue((root / "LG2.in").exists())
            self.assertTrue(report["gates"][0]["pass"])

    def test_stage_b_missing_dependency_fails_before_runner(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            calls = []
            report = build.build_stage_b(root, runner=lambda cmd, cwd: calls.append(cmd), ambertools_prefix=root / "missing")
            self.assertEqual(report["status"], "NOT_SUBMITTED_DEPENDENCY_MISSING")
            self.assertEqual(calls, [])

    def test_stage_b_tleap_marks_cyx_and_explicit_disulfide_bonds(self):
        with tempfile.TemporaryDirectory() as td:
            out = Path(td)
            tleap = build._write_tleap_input(out)
            text = tleap.read_text()
            self.assertNotIn("set p.238.name CYX", text)
            self.assertIn("bond p.238.SG p.283.SG", text)
            self.assertIn("bond p.275.SG p.292.SG", text)

    def test_stage_b_disulfide_audit_requires_cyx_no_hg_and_two_ss_bonds(self):
        atoms = [
            build.Atom(238, "CYX", 238, "SG", (0,0,0), "S", "ATOM"),
            build.Atom(283, "CYX", 283, "SG", (2.07,0,0), "S", "ATOM"),
            build.Atom(275, "CYX", 275, "SG", (0,4,0), "S", "ATOM"),
            build.Atom(292, "CYX", 292, "SG", (2.11,4,0), "S", "ATOM"),
        ]
        ok = build.audit_disulfides(atoms, [(238, 283), (275, 292)])
        self.assertTrue(ok["pass"])
        bad = atoms + [build.Atom(999, "CYX", 238, "HG", (0,0,1), "H", "ATOM")]
        self.assertFalse(build.audit_disulfides(bad, [(238, 283), (275, 292)])["pass"])

    def test_collector_parses_amber18_verbose_scc_and_final_results_energy(self):
        amber = """Header
QMMM SCC-DFTB: SCC-DFTB for step 1 converged in 2 cycles.
 DFTBESCF = -11256.4565 EGB = -2188.3072
 FINAL RESULTS
 NSTEP       ENERGY          RMS            GMAX
     1  -1.8090E+04   0.0E+00       4.3766E+02
 Maximum number of minimization cycles reached
"""
        with tempfile.TemporaryDirectory() as td:
            log = Path(td) / "lg2.out"
            log.write_text(amber)
            parsed = collect.parse_sander_output(log)
            self.assertTrue(parsed["scc_converged"])
            self.assertTrue(parsed["normal_termination"])
            self.assertAlmostEqual(parsed["TOTAL_ENERGY"], -1.8090e4)
            self.assertAlmostEqual(parsed["DFTBESCF"], -11256.4565)
            self.assertAlmostEqual(parsed["EGB"], -2188.3072)

    def test_postrun_geometry_cli_writes_json_before_collector_reads_it(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            pdb = root / "state.pdb"
            pdb.write_text("".join([
                "ATOM      1  CB  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n",
                "HETATM    2  C1  LG1 A 900       3.000   0.000   0.000  1.00  0.00           C\n",
                "END\n",
            ]))
            geom = root / "postrun_geometry.json"
            rc = audit.main([str(pdb), "--geometry-json", str(geom)])
            self.assertEqual(rc, 0)
            self.assertTrue(json.loads(geom.read_text())["pass"])

    def test_stage_b_protein_only_preprocesses_templates_before_tleap(self):
        lines = [
            "ATOM      1  SG  CYS A 238       0.000   0.000   0.000  1.00  0.00           S\n",
            "ATOM      2  HG  CYS A 238       0.000   0.000   1.000  1.00  0.00           H\n",
            "ATOM      3  ND1 HIS A 242       1.000   0.000   0.000  1.00  0.00           N\n",
            "ATOM      4  HD1 HIS A 242       1.000   0.000   1.000  1.00  0.00           H\n",
            "ATOM      5  NE2 HIS A 200       2.000   0.000   0.000  1.00  0.00           N\n",
            "ATOM      6  CB  ALA A 100       3.000   0.000   0.000  1.00  0.00           C\n",
            "HETATM    7  C1  LG1 A 900       4.000   0.000   0.000  1.00  0.00           C\n",
            "END\n",
        ]
        with tempfile.TemporaryDirectory() as td:
            out = Path(td)
            src = out / "generated_LG1_pair.pdb"
            src.write_text("".join(lines))
            protein = build.write_stage_b_protein_only_pdb(src, out / "stage_b_protein_only.pdb")
            text = protein.read_text()
            self.assertIn("CYX A 238", text)
            self.assertIn("HID A 242", text)
            self.assertIn("HIE A 200", text)
            self.assertNotIn(" HG ", text)
            self.assertNotIn("HETATM", text)
            tleap = build._write_tleap_input(out).read_text()
            self.assertNotIn("set p.238.name CYX", tleap)
            self.assertIn("bond p.238.SG p.283.SG", tleap)

    def test_stage_b_parmed_script_contains_real_api_and_save_calls(self):
        with tempfile.TemporaryDirectory() as td:
            out = Path(td)
            deps = {"parmed_egg": out / "ParmEd.egg", "literature_prmtop": out / "vmd-md-b.prmtop", "literature_inpcrd": out / "vmd-md-b.inpcrd"}
            for path in deps.values(): path.write_text("fixture")
            script = build._write_parmed_stage_b_script(out, {"parmed_egg": deps["parmed_egg"], "literature_prmtop": deps["literature_prmtop"], "literature_inpcrd": deps["literature_inpcrd"]})
            text = script.read_text()
            for needle in ["import parmed", "load_file", "ChamberParm.from_structure", "changeRadii", "save('pair.prmtop'", "save('LG1.inpcrd'", "save('LG2.inpcrd'", "iqmatoms"]:
                self.assertIn(needle, text)
            self.assertNotIn("Production intent", text)

    def test_stage_b_requires_auditable_manifest_and_writes_real_qm_mask(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deps = {}
            for key, name in {"ambertools_prefix":"ambertools", "parmed_egg":"ParmEd.egg", "literature_prmtop":"vmd-md-b.prmtop", "literature_inpcrd":"vmd-md-b.inpcrd", "dftb_slko_path":"3ob-3-1"}.items():
                path = root / name
                if key.endswith("path") or key == "ambertools_prefix": path.mkdir()
                else: path.write_text("fixture")
                deps[key] = path
            (deps["ambertools_prefix"] / "bin").mkdir(); (deps["ambertools_prefix"] / "bin" / "tleap").write_text("tleap")
            (root / "generated_LG1_pair.pdb").write_text("END\n"); (root / "generated_LG2_pair.pdb").write_text("END\n")
            def fake_runner(cmd, cwd):
                joined = " ".join(map(str, cmd))
                if "tleap" in joined:
                    (cwd / "protein.prmtop").write_text("protein"); (cwd / "protein.rst7").write_text("coords")
                if "parmed_stage_b" in joined:
                    for name in ["pair.prmtop","LG1.inpcrd","LG2.inpcrd"]: (cwd / name).write_text(name)
                    (cwd / "stage_b_parmed_manifest.json").write_text(json.dumps({"status":"PASS", "iqmatoms":[1,2,3], "qm_atom_count":100, "residues":259, "atoms":3902, "mbondi3":True, "gates":[{"name":"dynamic","pass":True}]}))
                return build.CommandResult(0, "ok", "")
            report = build.build_stage_b(root, runner=fake_runner, **deps)
            self.assertEqual(report["status"], "PASS")
            self.assertNotIn("STAGE_B_QM_MASK_PLACEHOLDER", (root / "LG1.in").read_text())
            self.assertIn("qmmask='@1,2,3'", (root / "LG1.in").read_text())
            self.assertEqual(report["topology_manifest"]["qm_atom_count"], 100)
            # Now prove touch-only outputs without manifest do not pass.
            root2 = root / "no_manifest"; root2.mkdir()
            for key, path in deps.items():
                pass
            (root2 / "generated_LG1_pair.pdb").write_text("END\n"); (root2 / "generated_LG2_pair.pdb").write_text("END\n")
            def touch_only(cmd, cwd):
                if "tleap" in " ".join(map(str, cmd)):
                    (cwd / "protein.prmtop").write_text("protein"); (cwd / "protein.rst7").write_text("coords")
                if "parmed_stage_b" in " ".join(map(str, cmd)):
                    for name in ["pair.prmtop","LG1.inpcrd","LG2.inpcrd"]: (cwd / name).write_text(name)
                return build.CommandResult(0, "ok", "")
            report2 = build.build_stage_b(root2, runner=touch_only, **deps)
            self.assertEqual(report2["status"], "NOT_SUBMITTED_TOPOLOGY_AUDIT_FAILED")

    def test_parmed_script_overrides_lg1_lg2_ligand_coordinates_separately(self):
        with tempfile.TemporaryDirectory() as td:
            out = Path(td)
            deps = {"parmed_egg": out / "ParmEd.egg", "literature_prmtop": out / "vmd-md-b.prmtop", "literature_inpcrd": out / "vmd-md-b.inpcrd"}
            for path in deps.values(): path.write_text("fixture")
            script = build._write_parmed_stage_b_script(out, deps).read_text()
            self.assertIn("generated_LG1_pair.pdb", script)
            self.assertIn("generated_LG2_pair.pdb", script)
            self.assertIn("apply_ligand_xyz", script)
            self.assertIn("LG1_ligand_generated_pair_max_delta_A", script)
            self.assertIn("LG2_ligand_generated_pair_max_delta_A", script)
            self.assertIn("protein_coordinate_max_delta_A", script)
            self.assertNotIn("Production coordinate override", script)

    def test_stage_b_uses_usr_bin_python27_with_pythonpath_env(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deps = {}
            for key, name in {"ambertools_prefix":"ambertools", "parmed_egg":"ParmEd.egg", "literature_prmtop":"vmd-md-b.prmtop", "literature_inpcrd":"vmd-md-b.inpcrd", "dftb_slko_path":"3ob-3-1"}.items():
                path = root / name
                if key.endswith("path") or key == "ambertools_prefix": path.mkdir()
                else: path.write_text("fixture")
                deps[key] = path
            (deps["ambertools_prefix"] / "bin").mkdir(); (deps["ambertools_prefix"] / "bin" / "tleap").write_text("tleap")
            calls = []
            def fake_runner(cmd, cwd):
                calls.append(cmd)
                joined = " ".join(map(str, cmd))
                if "tleap" in joined:
                    (cwd / "protein.prmtop").write_text("protein"); (cwd / "protein.rst7").write_text("coords")
                if "parmed_stage_b" in joined:
                    for name in ["pair.prmtop","LG1.inpcrd","LG2.inpcrd"]: (cwd / name).write_text(name)
                    (cwd / "stage_b_parmed_manifest.json").write_text(json.dumps({"status":"PASS", "iqmatoms":[1,2,3], "gates":[{"name":"dynamic","pass":True}], "qm_atom_count":100}))
                return build.CommandResult(0, "ok", "")
            build.build_stage_b(root, runner=fake_runner, **deps)
            parmed_cmd = next(c for c in calls if "parmed_stage_b" in " ".join(map(str, c)))
            self.assertEqual(parmed_cmd[0], "/usr/bin/env")
            self.assertIn(f"PYTHONPATH={deps['parmed_egg']}", parmed_cmd)
            self.assertIn("/usr/bin/python2.7", parmed_cmd)

    def test_stage_b_manifest_gates_must_all_pass_dynamically(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deps = {}
            for key, name in {"ambertools_prefix":"ambertools", "parmed_egg":"ParmEd.egg", "literature_prmtop":"vmd-md-b.prmtop", "literature_inpcrd":"vmd-md-b.inpcrd", "dftb_slko_path":"3ob-3-1"}.items():
                path = root / name
                if key.endswith("path") or key == "ambertools_prefix": path.mkdir()
                else: path.write_text("fixture")
                deps[key] = path
            (deps["ambertools_prefix"] / "bin").mkdir(); (deps["ambertools_prefix"] / "bin" / "tleap").write_text("tleap")
            def fake_runner(cmd, cwd):
                if "tleap" in " ".join(map(str, cmd)):
                    (cwd / "protein.prmtop").write_text("protein"); (cwd / "protein.rst7").write_text("coords")
                if "parmed_stage_b" in " ".join(map(str, cmd)):
                    for name in ["pair.prmtop","LG1.inpcrd","LG2.inpcrd"]: (cwd / name).write_text(name)
                    (cwd / "stage_b_parmed_manifest.json").write_text(json.dumps({"status":"PASS", "iqmatoms":[1,2,3], "gates":[{"name":"ligand_charges", "pass":False}], "qm_atom_count":100}))
                return build.CommandResult(0, "ok", "")
            report = build.build_stage_b(root, runner=fake_runner, **deps)
            self.assertEqual(report["status"], "NOT_SUBMITTED_TOPOLOGY_AUDIT_FAILED")

    def test_sbatch_uses_single_rank_for_fixed_geometry_single_point(self):
        text = Path("run_iccg_step1_pair.sbatch").read_text()
        self.assertIn("#SBATCH --ntasks=1", text)
        self.assertNotIn("#SBATCH --ntasks=16", text)


if __name__ == "__main__":
    unittest.main(verbosity=2)
