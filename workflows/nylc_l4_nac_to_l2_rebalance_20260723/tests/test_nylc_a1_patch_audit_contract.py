#!/usr/bin/env python3
import copy
import importlib.util
import math
import unittest
from pathlib import Path

import test_nylc_a1_patch_contract as model_contract

HERE = Path(__file__).resolve().parents[1]
AUDITOR = HERE / "scripts" / "audit_nylc_a1_patch.py"
BUILDER = HERE / "scripts" / "prepare_nylc_a1_parameter_model.py"
RUNNER = HERE / "scripts" / "run_nylc_a1_parameterization.sh"
SBATCH = HERE / "slurm" / "run_nylc_a1_parameterization.sbatch"


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def valid_inputs():
    builder = load(BUILDER, "prepare_nylc_a1_parameter_model_for_audit")
    thr, nxt, bonds = model_contract.parent_fixture()
    model = builder.build_a1_capped_model(thr, nxt, bonds)
    mol2 = builder.render_mol2(model)
    provenance = {
        "model_definition": "N-terminal Thr(A1)-N-methylamide capped model",
        "net_charge_e": 0,
        "charge_method": "AM1-BCC",
        "atom_types": "GAFF2",
        "bonded_parameters": "parmchk2",
        "tool_versions": {
            "antechamber": "AmberTools 18",
            "parmchk2": "AmberTools 18",
            "tleap": "AmberTools 18",
        },
        "input_sha256": "0" * 64,
        "output_sha256": "1" * 64,
        "validation_status": "PASS",
        "parmchk2_unresolved_count": 0,
        "tleap_load_status": "PASS",
        "gromacs_grompp_status": "PASS",
        "amber_energy_kcal_mol": -12.5,
        "gromacs_energy_kj_mol": -50.0,
    }
    return mol2, "MASS\n\nBOND\n\nANGLE\n\nDIHE\n\nIMPROPER\n\nNONBON\n", provenance


class A1PatchAuditContractTests(unittest.TestCase):
    def test_valid_neutral_a1_patch_passes(self):
        module = load(AUDITOR, "audit_nylc_a1_patch")
        mol2, frcmod, provenance = valid_inputs()
        result = module.audit_patch(mol2, frcmod, provenance)
        self.assertEqual(result["status"], "PASS_A1_SCREENING_PATCH")
        self.assertEqual(result["atom_count"], 21)
        self.assertAlmostEqual(result["charge_e"], 0.0, places=6)

    def test_unresolved_parmchk2_term_is_rejected(self):
        module = load(AUDITOR, "audit_nylc_a1_patch_unresolved")
        mol2, frcmod, provenance = valid_inputs()
        frcmod += "ATTN, need revision\n"
        with self.assertRaisesRegex(module.PatchAuditError, "unresolved"):
            module.audit_patch(mol2, frcmod, provenance)

    def test_nonintegral_charge_is_rejected(self):
        module = load(AUDITOR, "audit_nylc_a1_patch_charge")
        mol2, frcmod, provenance = valid_inputs()
        mol2 = mol2.replace("  1.000000", "  0.900000", 1)
        with self.assertRaisesRegex(module.PatchAuditError, "charge"):
            module.audit_patch(mol2, frcmod, provenance)

    def test_og1_hg1_bond_is_rejected(self):
        module = load(AUDITOR, "audit_nylc_a1_patch_bond")
        mol2, frcmod, provenance = valid_inputs()
        lines = mol2.splitlines()
        atom_ids = {}
        in_atoms = False
        for line in lines:
            if line == "@<TRIPOS>ATOM":
                in_atoms = True
                continue
            if line == "@<TRIPOS>BOND":
                break
            if in_atoms and line.strip():
                fields = line.split()
                atom_ids[fields[1]] = int(fields[0])
        bond_marker = lines.index("@<TRIPOS>BOND")
        lines.insert(bond_marker + 1, f"999 {atom_ids['OG1']} {atom_ids['HG1']} 1")
        bad = "\n".join(lines) + "\n"
        with self.assertRaisesRegex(module.PatchAuditError, "OG1-HG1"):
            module.audit_patch(bad, frcmod, provenance)

    def test_nonfinite_validation_energy_is_rejected(self):
        module = load(AUDITOR, "audit_nylc_a1_patch_energy")
        mol2, frcmod, provenance = valid_inputs()
        provenance = copy.deepcopy(provenance)
        provenance["amber_energy_kcal_mol"] = math.nan
        with self.assertRaisesRegex(module.PatchAuditError, "finite"):
            module.audit_patch(mol2, frcmod, provenance)

    def test_driver_names_required_ambertools_and_history(self):
        text = RUNNER.read_text(encoding="utf-8")
        for token in [
            "antechamber", "parmchk2", "tleap", "sqm", "msander",
            "run_history.tsv", "run_history.jsonl",
            "audit_nylc_a1_patch.py",
            "charge_normalization.json",
            "NTA1_CAP.am1bcc.raw.mol2",
            "editconf",
            "NTA1_CAP.gromacs.boxed.gro",
            "pbc                      = xyz",
            "continuation             = yes",
        ]:
            self.assertIn(token, text)

    def test_driver_requires_audited_task_local_amberclassic(self):
        text = RUNNER.read_text(encoding="utf-8")
        for token in [
            "ACTIVE_AMBERCLASSIC.json",
            "PASS_AMBERCLASSIC_INSTALL",
            "AmberClassic.sh",
            "AMBER_PREFIX",
            "source_archive_sha256",
            "set +u",
            "set -u",
            "compiler/gcc/11.4.0",
        ]:
            self.assertIn(token, text)
        self.assertNotIn("module load amber/2018", text)
        self.assertNotIn("\\nsander -O", text)

    def test_parameterization_runs_through_slurm(self):
        self.assertTrue(SBATCH.is_file())
        text = SBATCH.read_text(encoding="utf-8")
        for token in [
            "#SBATCH -p xahcnormal",
            "#SBATCH -c 4",
            "run_nylc_a1_parameterization.sh",
            "SLURM_JOB_ID",
        ]:
            self.assertIn(token, text)

    def test_exit_trap_disarms_before_returning_status(self):
        text = RUNNER.read_text(encoding="utf-8")
        finish_body = text.split("finish() {", 1)[1].split("}", 1)[0]
        self.assertIn("trap - EXIT", finish_body)

    def test_driver_records_phase_boundaries_before_ambertools(self):
        text = RUNNER.read_text(encoding="utf-8")
        self.assertIn("phase.tsv", text)
        for phase in [
            "before_module_purge",
            "after_module_purge",
            "after_active_manifest",
            "after_amberclassic_source",
            "after_tool_resolution",
        ]:
            self.assertIn(phase, text)


if __name__ == "__main__":
    unittest.main()
