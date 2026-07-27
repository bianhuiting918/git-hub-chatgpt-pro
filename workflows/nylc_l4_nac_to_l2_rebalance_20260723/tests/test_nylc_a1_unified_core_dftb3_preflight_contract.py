#!/usr/bin/env python3
import importlib.util
import pathlib
import sys
import types
import unittest
from unittest import mock

HERE = pathlib.Path(__file__).resolve().parents[1]
PREPARE = HERE / "scripts" / "prepare_nylc_a1_unified_core_dftb3_preflight.py"
AUDIT = HERE / "scripts" / "audit_nylc_a1_unified_core_dftb3_smoke.py"
RUNNER = HERE / "scripts" / "run_nylc_a1_unified_core_dftb3_preflight.sh"
SLURM = HERE / "slurm" / "run_nylc_a1_unified_core_dftb3_preflight.sbatch"


class A1UnifiedCoreDftb3ContractTests(unittest.TestCase):
    def test_files_exist(self):
        for path in (PREPARE, AUDIT, RUNNER, SLURM):
            self.assertTrue(path.is_file(), path)

    def test_protein_residue_set_and_formal_contract(self):
        text = PREPARE.read_text(encoding="utf-8")
        for token in (
            "EXPECTED_TYR146_GLOBAL_RESID = 501",
            "EXPECTED_LYS189_GLOBAL_RESID = 544",
            "EXPECTED_ASN219_GLOBAL_RESID = 574",
            "EXPECTED_ACTIVE_GLOBAL_RESID = 622",
            "EXPECTED_ASP306_GLOBAL_RESID = 661",
            "EXPECTED_ASP308_GLOBAL_RESID = 663",
            "QMCHARGE = 0",
            "EXPECTED_QM_ATOMS = 146",
            "EXPECTED_LINK_ATOMS = 6",
            "EXPECTED_ELECTRONS_WITH_LINKS = 510",
            "TYR146_SIDECHAIN_NAMES",
            "LYS189_SIDECHAIN_NAMES",
            "ASN219_SIDECHAIN_NAMES",
            "ASP306_SIDECHAIN_NAMES",
            "ASP308_SIDECHAIN_NAMES",
        ):
            self.assertIn(token, text)

    def test_exact_six_boundaries_are_frozen(self):
        text = PREPARE.read_text(encoding="utf-8")
        for token in (
            "(8962, 8964)",
            "(9567, 9565)",
            "(9587, 9585)",
            "(7160, 7158)",
            "(7756, 7754)",
            "(8235, 8233)",
        ):
            self.assertIn(token, text)

    def test_step1_step2_scope_is_explicit(self):
        combined = PREPARE.read_text(encoding="utf-8") + AUDIT.read_text(encoding="utf-8")
        for token in (
            "Tyr146",
            "Lys189",
            "Asn219",
            "Asp306",
            "Asp308",
            "complete L2",
            "Step1",
            "Step2",
            "no QM water",
            "Step2 water",
            "numerical preflight",
            "not a TS",
            "PMF",
        ):
            self.assertIn(token, combined)

    def test_input_interpolates_zero_charge_contract(self):
        fake_parmed = types.ModuleType("parmed")
        with mock.patch.dict(sys.modules, {"parmed": fake_parmed}):
            spec = importlib.util.spec_from_file_location("a1_unified_prepare", PREPARE)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
        rendered = module.qmmm_input("unified", maxcyc=20, qmmask="@1,2,3")
        self.assertIn("maxcyc=20", rendered)
        self.assertIn("qmmask='@1,2,3'", rendered)
        self.assertIn("qmcharge=0", rendered)
        self.assertIn("spin=1", rendered)
        self.assertNotIn("{", rendered)
        self.assertNotIn("}", rendered)

    def test_auditor_fails_closed(self):
        text = AUDIT.read_text(encoding="utf-8")
        for token in (
            "qmcharge",
            "qm_atom_count",
            "electron_count_including_link_h",
            "link_atom_count",
            "max_bond_length_A",
            "bond_count_gt_3A",
            "bond_overflow",
            "Convergence could not be achieved",
            "PASS_A1_UNIFIED_CORE_DFTB3_NUMERICAL_PREFLIGHT",
            "FAIL_A1_UNIFIED_CORE_DFTB3_NUMERICAL_PREFLIGHT",
        ):
            self.assertIn(token, text)

    def test_runner_is_immutable_whole_and_slurm_safe_mpi(self):
        text = RUNNER.read_text(encoding="utf-8")
        for token in (
            "gmx_mpi trjconv",
            "-pbc mol",
            "-ur compact",
            "mpirun --bind-to none",
            "SLURM_NTASKS",
            "run_history.tsv",
            "run_history.jsonl",
            "NOT_EVALUATED.json",
            "PASS.json",
        ):
            self.assertIn(token, text)
        for forbidden in ("gmx mdrun", "pmemd.cuda", "srun -n"):
            self.assertNotIn(forbidden, text)

    def test_slurm_uses_eight_cpu_ranks_without_dcu(self):
        text = SLURM.read_text(encoding="utf-8")
        for token in (
            "#SBATCH -p xahcnormal",
            "#SBATCH -n 8",
            "#SBATCH -c 1",
            "#SBATCH --mem-per-cpu=2500M",
            "A1_UNIFIED_CORE_DFTB3_CODE_SOURCE",
            "SNAPSHOT_SHA256.tsv",
        ):
            self.assertIn(token, text)
        self.assertNotIn("#SBATCH --gres", text)


if __name__ == "__main__":
    unittest.main()
