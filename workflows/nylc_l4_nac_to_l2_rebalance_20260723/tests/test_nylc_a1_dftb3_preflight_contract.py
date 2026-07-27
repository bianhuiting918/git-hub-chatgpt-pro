#!/usr/bin/env python3
import importlib.util
import pathlib
import sys
import types
import unittest
from unittest import mock

HERE = pathlib.Path(__file__).resolve().parents[1]
PREPARE = HERE / "scripts" / "prepare_nylc_a1_dftb3_preflight.py"
AUDIT = HERE / "scripts" / "audit_nylc_a1_dftb3_smoke.py"
RUNNER = HERE / "scripts" / "run_nylc_a1_dftb3_preflight.sh"
SLURM = HERE / "slurm" / "run_nylc_a1_dftb3_preflight.sbatch"


class A1Dftb3PreflightContractTests(unittest.TestCase):
    def test_a1_specific_files_exist(self):
        for path in (PREPARE, AUDIT, RUNNER, SLURM):
            self.assertTrue(path.is_file(), path)

    def test_preparer_freezes_passed_representative_and_topology(self):
        text = PREPARE.read_text(encoding="utf-8")
        required = (
            "attempt_61990814",
            "PASS_A1_REPRESENTATIVE_FRAME_EXTRACTION",
            "PASS_A1_REPRESENTATIVE_NAC_FRAME",
            "7b095b7d36a25327558b4bbab20a98a43edaab51d0cacaf7a888ff886bafe679",
            "8a61a57be0537a5fb0aac8f2379452fdb6980468cb06659c0953224d3a0aefa3",
            "af98733e218a8f83d0a5c46120d9d230a16c222f76d230654d44b542288cc205",
            "8fd4398af1356b515720c1da3126d08b7795b24b3c112d98ef3235afa57c8179",
            "b0e753c60fd4b71c282d21cc6106a15e73d91d12a20d80e92dd01516162eb301",
        )
        for token in required:
            self.assertIn(token, text)
        self.assertNotIn("trjconv", text)
        self.assertNotIn("lowest_potential", text)

    def test_qm_contract_is_derived_and_fail_closed(self):
        text = PREPARE.read_text(encoding="utf-8")
        for token in (
            "PROTEIN_ATOMS = 10272",
            "THR267_OG1 = 8960",
            "EXPECTED_ACTIVE_GLOBAL_RESID = 622",
            "EXPECTED_ACTIVE_ORIGINAL_RESID = 267",
            "source_global_resid",
            "extracted_original_resid",
            "L2_FIRST = 10273",
            "L2_LAST = 10351",
            "L2_REACTIVE_C = 10287",
            "L2_REACTIVE_O = 10288",
            "L2_REACTIVE_N = 10289",
            "QMCHARGE = 0",
            "SPIN = 1",
            "EXPECTED_QM_ATOMS = 94",
            "EXPECTED_ELECTRONS_WITH_LINKS = 314",
            "expected one QM/MM boundary bond",
            "N bonded hydrogens",
            "OG1 bonded atoms",
            "active_residue_topology_charge",
            "ligand_topology_charge",
        ):
            self.assertIn(token, text)

    def test_amber_inputs_are_dftb3_200k_one_and_twenty_steps(self):
        text = PREPARE.read_text(encoding="utf-8")
        for token in (
            "qm_theory='DFTB3'",
            "dftb_telec=200.0",
            "qmcharge={QMCHARGE}",
            "spin={SPIN}",
            "maxcyc=1",
            "maxcyc=20",
            "3ob-3-1",
        ):
            self.assertIn(token, text)

    def test_auditor_fails_on_numerical_warnings(self):
        text = AUDIT.read_text(encoding="utf-8")
        for token in (
            "Convergence could not be achieved",
            "vlimit",
            "SANDER BOMB",
            "segmentation",
            "forrtl",
            r"\bnan\b",
            "FATAL",
            "PASS_A1_DFTB3_NUMERICAL_PREFLIGHT",
            "FAIL_A1_DFTB3_NUMERICAL_PREFLIGHT",
        ):
            self.assertIn(token, text)

    def test_runner_uses_amber18_and_never_runs_mm_or_selects_frames(self):
        text = RUNNER.read_text(encoding="utf-8")
        for token in (
            "Amber18",
            "module load amber/2018-hpcx-gcc-7.3.1",
            "GMXDATA=/public/software/apps/Gromacs-DCU2/2022.1/mpi/share/gromacs",
            'export AMBERHOME="$AMBER_RUNTIME"',
            "3ob-3-1",
            "sander",
            "01_qmmm_one_step.in",
            "02_qmmm_20_step.in",
            "run_history.tsv",
            "run_history.jsonl",
            "flock -x",
            "NOT_EVALUATED.json",
            "PASS.json",
            "code snapshot",
        ):
            self.assertIn(token, text)
        self.assertLess(text.index("export GMXDATA"), text.index("CURRENT=preparation"))
        for forbidden in ("gmx mdrun", "trjconv", "lowest_potential", 'source "$AMBER_RUNTIME/amber.sh"'):
            self.assertNotIn(forbidden, text)

    def test_slurm_is_cpu_only_bounded_and_immutable(self):
        text = SLURM.read_text(encoding="utf-8")
        for token in (
            "#SBATCH -p xahcnormal",
            "#SBATCH -c 4",
            "#SBATCH --mem-per-cpu=2500M",
            "#SBATCH -t 02:00:00",
            "A1_DFTB3_CODE_SOURCE",
            "SNAPSHOT_SHA256.tsv",
            "sha256sum -c",
        ):
            self.assertIn(token, text)
        self.assertNotIn("#SBATCH --gres", text)
        self.assertNotIn("sbatch ", text)

    def test_scope_is_numerical_only_and_production_qm_remains_expanded(self):
        combined = PREPARE.read_text(encoding="utf-8") + AUDIT.read_text(encoding="utf-8")
        for token in ("numerical preflight", "not a TS", "PMF", "barrier", "Asp306", "Asp308"):
            self.assertIn(token, combined)


    def test_qmmm_input_interpolates_complete_numeric_contract(self):
        fake_parmed = types.ModuleType("parmed")
        with mock.patch.dict(sys.modules, {"parmed": fake_parmed}):
            spec = importlib.util.spec_from_file_location("a1_dftb3_prepare", PREPARE)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
        one = module.qmmm_input("one", maxcyc=1, qmmask="@1,2,3")
        twenty = module.qmmm_input("twenty", maxcyc=20, qmmask="@1,2,3")
        self.assertIn("maxcyc=1", one)
        self.assertIn("maxcyc=20", twenty)
        for rendered in (one, twenty):
            self.assertIn("qmmask='@1,2,3'", rendered)
            self.assertIn("qmcharge=0", rendered)
            self.assertIn("spin=1", rendered)
            self.assertNotIn("{", rendered)
            self.assertNotIn("}", rendered)

    def test_runner_transactionally_records_terminal_history_and_demotes_pass(self):
        text = RUNNER.read_text(encoding="utf-8")
        for token in (
            "PROMOTED=0",
            "demote_promoted_outputs",
            "FAILED_NOT_PROMOTED_",
            "os.fsync",
            "truncate",
            "if ((PROMOTED))",
        ):
            self.assertIn(token, text)
        self.assertLess(text.rindex("append_history"), text.rindex("trap - EXIT"))


if __name__ == "__main__":
    unittest.main()
