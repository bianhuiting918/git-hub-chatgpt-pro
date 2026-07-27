#!/usr/bin/env python3
import importlib.util
import pathlib
import sys
import types
import unittest
from unittest import mock

HERE = pathlib.Path(__file__).resolve().parents[1]
PREPARE = HERE / "scripts" / "prepare_nylc_a1_expanded_dftb3_preflight.py"
AUDIT = HERE / "scripts" / "audit_nylc_a1_expanded_dftb3_smoke.py"
RUNNER = HERE / "scripts" / "run_nylc_a1_expanded_dftb3_preflight.sh"
SLURM = HERE / "slurm" / "run_nylc_a1_expanded_dftb3_preflight.sbatch"


class A1ExpandedDftb3ContractTests(unittest.TestCase):
    def test_expanded_files_exist(self):
        for path in (PREPARE, AUDIT, RUNNER, SLURM):
            self.assertTrue(path.is_file(), path)

    def test_preparer_freezes_expanded_microstate_and_charge(self):
        text = PREPARE.read_text(encoding="utf-8")
        for token in (
            "attempt_61990814",
            "PASS_A1_REPRESENTATIVE_NAC_FRAME",
            "EXPECTED_ASP306_GLOBAL_RESID = 661",
            "EXPECTED_ASP308_GLOBAL_RESID = 663",
            'EXPECTED_ASP306_NAME = "ASH"',
            'EXPECTED_ASP308_NAME = "ASP"',
            "ASP306_CB = 9567",
            "ASP308_CB = 9587",
            "QMCHARGE = -1",
            "SPIN = 1",
            "EXPECTED_QM_ATOMS = 107",
            "EXPECTED_LINK_ATOMS = 3",
            "EXPECTED_ELECTRONS_WITH_LINKS = 378",
            "HD2",
            "max_bond_length_A",
            "bond_count_gt_3A",
            "MAX_ALLOWED_BOND_LENGTH_A = 2.0",
        ):
            self.assertIn(token, text)

    def test_preparer_selects_complete_ash_and_asp_sidechains(self):
        text = PREPARE.read_text(encoding="utf-8")
        for token in (
            '"CB", "HB1", "HB2", "CG", "OD1", "OD2", "HD2"',
            '"CB", "HB1", "HB2", "CG", "OD1", "OD2"',
            "expected three QM/MM boundary bonds",
            "asp306_sidechain_atoms",
            "asp308_sidechain_atoms",
            "active_residue_topology_charge",
            "asp306_residue_topology_charge",
            "asp308_residue_topology_charge",
        ):
            self.assertIn(token, text)

    def test_expanded_input_interpolates_numeric_contract(self):
        fake_parmed = types.ModuleType("parmed")
        with mock.patch.dict(sys.modules, {"parmed": fake_parmed}):
            spec = importlib.util.spec_from_file_location("a1_expanded_prepare", PREPARE)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
        rendered = module.qmmm_input("expanded", maxcyc=20, qmmask="@1,2,3")
        self.assertIn("maxcyc=20", rendered)
        self.assertIn("qmmask='@1,2,3'", rendered)
        self.assertIn("qmcharge=-1", rendered)
        self.assertIn("spin=1", rendered)
        self.assertNotIn("{", rendered)
        self.assertNotIn("}", rendered)

    def test_bond_length_uses_xyz_scalars(self):
        fake_parmed = types.ModuleType("parmed")
        with mock.patch.dict(sys.modules, {"parmed": fake_parmed}):
            spec = importlib.util.spec_from_file_location("a1_expanded_bond", PREPARE)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
        left = types.SimpleNamespace(xx=0.0, xy=0.0, xz=0.0)
        right = types.SimpleNamespace(xx=3.0, xy=4.0, xz=12.0)
        self.assertAlmostEqual(module.bond_length_A(left, right), 13.0)

    def test_auditor_fails_closed_for_expanded_contract(self):
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
            "SANDER BOMB",
            r"\bnan\b",
            "PASS_A1_EXPANDED_DFTB3_NUMERICAL_PREFLIGHT",
            "FAIL_A1_EXPANDED_DFTB3_NUMERICAL_PREFLIGHT",
        ):
            self.assertIn(token, text)

    def test_runner_is_whole_job_local_and_mpi(self):
        text = RUNNER.read_text(encoding="utf-8")
        for token in (
            "gmx_mpi trjconv",
            "-pbc mol",
            "-ur compact",
            "--coordinate",
            "sander.MPI",
            "mpirun",
            "SLURM_NTASKS",
            "-inf",
            "01_qmmm_one_step.mdinfo",
            "02_qmmm_20_step.mdinfo",
            "run_history.tsv",
            "run_history.jsonl",
            "NOT_EVALUATED.json",
            "PASS.json",
        ):
            self.assertIn(token, text)
        for forbidden in ("gmx mdrun", "lowest_potential", "pmemd.cuda", "srun -n"):
            self.assertNotIn(forbidden, text)

    def test_slurm_uses_eight_cpu_ranks_without_dcu_request(self):
        text = SLURM.read_text(encoding="utf-8")
        for token in (
            "#SBATCH -p xahcnormal",
            "#SBATCH -n 8",
            "#SBATCH -c 1",
            "#SBATCH --mem-per-cpu=2500M",
            "#SBATCH -t 02:00:00",
            "A1_EXPANDED_DFTB3_CODE_SOURCE",
            "SNAPSHOT_SHA256.tsv",
        ):
            self.assertIn(token, text)
        self.assertNotIn("#SBATCH --gres", text)

    def test_scope_remains_preflight_not_mechanism(self):
        combined = PREPARE.read_text(encoding="utf-8") + AUDIT.read_text(encoding="utf-8")
        for token in (
            "numerical preflight",
            "not a TS",
            "PMF",
            "barrier",
            "Asp306",
            "Asp308",
            "Tyr146",
            "Lys189",
            "Asn219",
        ):
            self.assertIn(token, combined)


if __name__ == "__main__":
    unittest.main()
