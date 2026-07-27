#!/usr/bin/env python3
import json
import unittest
from pathlib import Path

HERE=Path(__file__).resolve().parents[1]
MANIFEST=HERE/"manifests"/"nylc_a1_nac_audit_universe.json"
RUNNER=HERE/"scripts"/"run_nylc_a1_nac_audit.sh"
AUDITOR=HERE/"scripts"/"audit_nylc_a1_nac_replica.py"
SBATCH=HERE/"slurm"/"run_nylc_a1_nac_audit.sbatch"
MERGER=HERE/"scripts"/"merge_nylc_a1_nac_audits.py"
MERGE_SBATCH=HERE/"slurm"/"run_nylc_a1_nac_merge.sbatch"

class A1NACAuditContract(unittest.TestCase):
    def test_manifest_has_exact_nine_slot_universe(self):
        self.assertTrue(MANIFEST.is_file())
        payload=json.loads(MANIFEST.read_text())
        rows=payload["replicas"]
        self.assertEqual(len(rows),9)
        self.assertEqual([row["slot"] for row in rows],list(range(9)))
        self.assertEqual(
            sorted({row["candidate_id"] for row in rows}),
            ["nac_evt08_time1086ps","nac_evt18_time1206ps","nac_evt25_time1462ps"],
        )
        for candidate in {row["candidate_id"] for row in rows}:
            self.assertEqual(
                sorted(row["velocity_seed"] for row in rows if row["candidate_id"]==candidate),
                [26711,26723,26737],
            )
        self.assertEqual(len({(row["candidate_id"],row["velocity_seed"]) for row in rows}),9)

    def test_manifest_freezes_scientific_contract(self):
        payload=json.loads(MANIFEST.read_text())
        self.assertEqual(payload["microstate"],"A1_Thr267_Ogamma_minus_NalphaH3_plus")
        self.assertEqual(payload["gate_residues"],[261,262,263,264,265,266])
        self.assertFalse(payload["gate_includes_thr267"])
        self.assertEqual(payload["nac"]["distance_max_nm"],0.35)
        self.assertEqual(payload["nac"]["angle_min_deg"],95.0)
        self.assertEqual(payload["nac"]["angle_max_deg"],115.0)
        self.assertEqual(payload["analysis_window_ps"],[0.0,1000.0])
        for row in payload["replicas"]:
            self.assertTrue(row["fully_unrestrained"])
            self.assertIn("/npt300free",row["free_run_root"])
            self.assertIn("source_cycle.ndx",row["source_cycle_ndx"])
            self.assertIn("EQUILIBRATION_COMPLETE.json",row["completion_manifest"])


    def test_a1_audit_runtime_files_exist(self):
        self.assertTrue(RUNNER.is_file())
        self.assertTrue(AUDITOR.is_file())
        self.assertTrue(SBATCH.is_file())

    @unittest.skipUnless(AUDITOR.is_file(),"A1 auditor missing")
    def test_auditor_uses_joint_nac_and_excludes_proton_path_claims(self):
        text=AUDITOR.read_text()
        for token in ["A1_Thr267_Ogamma_minus_NalphaH3_plus","distance_max_nm",
                      "angle_min_deg","angle_max_deg","longest_continuous_nac",
                      "gate_opening_nm","thermodynamics","numerical_issue_counts"]:
            self.assertIn(token,text)
        self.assertNotIn("proton_preorganization",text)

    @unittest.skipUnless(RUNNER.is_file(),"runner missing")
    def test_runner_generates_primitives_and_preserves_failures(self):
        text=RUNNER.read_text()
        for token in ["generate_nylc_m1_ensemble_primitives.py",
                      "audit_nylc_a1_nac_replica.py","gmx energy",
                      "NOT_EVALUATED","refusing to overwrite",
                      "run_history.tsv","run_history.jsonl","sha256sum"]:
            self.assertIn(token,text)

    @unittest.skipUnless(SBATCH.is_file(),"sbatch missing")
    def test_sbatch_is_nine_way_immutable_array(self):
        text=SBATCH.read_text()
        for token in ["#SBATCH --array=0-8%9","#SBATCH -p xahcnormal",
                      "code_snapshots","SNAPSHOT_SHA256.tsv",
                      "A1_NAC_CODE_ROOT","SLURM_ARRAY_TASK_ID"]:
            self.assertIn(token,text)


    def test_merge_runtime_files_exist(self):
        self.assertTrue(MERGER.is_file())
        self.assertTrue(MERGE_SBATCH.is_file())

    @unittest.skipUnless(MERGER.is_file(),"merger missing")
    def test_merger_keeps_denominator_and_ranks_only_eligible_frames(self):
        text=MERGER.read_text()
        for token in ["expected_slots","missing_slots","NOT_EVALUATED_MISSING_AUDIT",
                      "replica_denominator","candidate_denominator",
                      "eligible_replica_count","nac_occupancy",
                      "potential_energy_kj_mol","selected_replica"]:
            self.assertIn(token,text)

    @unittest.skipUnless(MERGE_SBATCH.is_file(),"merge sbatch missing")
    def test_merge_job_is_afterany_and_records_history(self):
        text=MERGE_SBATCH.read_text()
        for token in ["AUDIT_ARRAY_JOB_ID","afterany","run_history.tsv",
                      "run_history.jsonl","merge_nylc_a1_nac_audits.py"]:
            self.assertIn(token,text)

if __name__=="__main__":
    unittest.main()
