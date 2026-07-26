#!/usr/bin/env python3
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
SCRIPT = HERE / "scripts" / "audit_nylc_a1_parameter_source.py"
MANIFEST = HERE / "manifests" / "nylc_a1_activated_nac.authority.json"


def load_module():
    spec = importlib.util.spec_from_file_location("audit_nylc_a1_parameter_source", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class A1AuthorityContractTests(unittest.TestCase):
    def test_authority_freezes_atom_conserving_transfer(self):
        authority = json.loads(MANIFEST.read_text(encoding="utf-8"))
        self.assertEqual(authority["schema_version"], 1)
        self.assertEqual(authority["parent_audit"]["job_id"], 61841413)
        self.assertEqual(len(authority["parents"]), 3)
        self.assertEqual(
            [row["id"] for row in authority["parents"]],
            [
                "nac_evt18_time1206ps",
                "nac_evt25_time1462ps",
                "nac_evt08_time1086ps",
            ],
        )
        self.assertEqual(
            authority["target_microstate"],
            {
                "Thr267_Nalpha": "NH3+",
                "Thr267_Ogamma": "O-",
                "Asp306": "ASH",
                "Asp308": "ASP",
            },
        )
        self.assertEqual(
            authority["topology_delta"]["removed_bonds"],
            [["Thr267:OG1", "Thr267:HG1"]],
        )
        self.assertEqual(
            authority["topology_delta"]["added_bonds"],
            [["Thr267:N", "Thr267:HG1"]],
        )
        self.assertEqual(authority["invariants"]["atom_count_delta"], 0)
        self.assertEqual(authority["invariants"]["system_charge_delta_e"], 0)
        self.assertEqual(authority["force_field"]["name"], "amber99sb-ildn")
        self.assertEqual(
            authority["parameter_source"]["status"],
            "NOT_EVALUATED_PENDING_GENERATION",
        )

    def test_parent_coordinate_hashes_are_frozen(self):
        authority = json.loads(MANIFEST.read_text(encoding="utf-8"))
        expected = {
            "nac_evt18_time1206ps": "c8b9c8bf8988b7c7b1b470579894eb9b7501cc71bda6bb856e1c291d792e23dc",
            "nac_evt25_time1462ps": "e1c292b7da52d68bdc3c41213d934f8fabcdc5f3dd953282274d793598690fae",
            "nac_evt08_time1086ps": "6cddb322d02425c040afaa86dd9b2939e8dd933b9651b793b0a2ffa907611124",
        }
        self.assertEqual(
            {row["id"]: row["coordinate_sha256"] for row in authority["parents"]},
            expected,
        )

    def test_pending_parameter_source_refuses_mm_readiness(self):
        module = load_module()
        authority = json.loads(MANIFEST.read_text(encoding="utf-8"))
        with self.assertRaisesRegex(module.AuditError, "parameter source is not ready"):
            module.audit_authority(authority, require_ready=True)

    def test_ready_parameter_source_requires_complete_provenance(self):
        module = load_module()
        authority = json.loads(MANIFEST.read_text(encoding="utf-8"))
        authority["parameter_source"] = {
            "status": "PASS_A1_SCREENING_PATCH",
            "patch_path": "parameters/a1_patch.json",
        }
        with self.assertRaisesRegex(module.AuditError, "missing provenance"):
            module.audit_authority(authority, require_ready=True)

    def test_nonlocal_patch_modification_is_rejected(self):
        module = load_module()
        authority = json.loads(MANIFEST.read_text(encoding="utf-8"))
        patch = module.required_patch_fixture()
        patch["modified_atoms"].append("Asp306:OD2")
        authority["parameter_source"] = {
            "status": "PASS_A1_SCREENING_PATCH",
            "patch_path": "parameters/a1_patch.json",
            "provenance": patch["provenance"],
            "residue_charge_e": 0.0,
            "bonded_state": patch["bonded_state"],
            "modified_atoms": patch["modified_atoms"],
        }
        with self.assertRaisesRegex(module.AuditError, "nonlocal"):
            module.audit_authority(authority, require_ready=True)


if __name__ == "__main__":
    unittest.main()
