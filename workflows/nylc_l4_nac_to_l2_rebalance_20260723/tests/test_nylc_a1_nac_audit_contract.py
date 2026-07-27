#!/usr/bin/env python3
import json
import unittest
from pathlib import Path

HERE=Path(__file__).resolve().parents[1]
MANIFEST=HERE/"manifests"/"nylc_a1_nac_audit_universe.json"

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

if __name__=="__main__":
    unittest.main()
