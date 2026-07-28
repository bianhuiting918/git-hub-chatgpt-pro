#!/usr/bin/env python3
import importlib.util
import pathlib
import unittest

FLOW = pathlib.Path(__file__).resolve().parents[1]
PREPARE = FLOW / "scripts" / "prepare_audit_nylc_a1_step1_attack_inherited.py"
RUNNER = FLOW / "scripts" / "run_nylc_a1_step1_attack_inherited.sh"
SLURM = FLOW / "slurm" / "run_nylc_a1_step1_attack_inherited.sbatch"


def load_module():
    if not PREPARE.is_file():
        raise unittest.SkipTest("production driver absent at TDD RED")
    spec = importlib.util.spec_from_file_location("_attack_inherited", PREPARE)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load attack-only driver")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class AttackInheritedContract(unittest.TestCase):
    def test_minimal_production_files_exist(self):
        for path in (PREPARE, RUNNER, SLURM):
            self.assertTrue(path.is_file(), path)

    def test_exact_dual_seed_sequential_targets_and_force_ramp(self):
        module = load_module()
        self.assertEqual(
            module.ATTACK_TARGETS_A,
            (3.0, 2.8, 2.6, 2.4, 2.2, 2.0, 1.8, 1.65, 1.50),
        )
        self.assertEqual(
            module.ATTACK_FORCE_KCAL_MOL_A2,
            (10.0, 10.0, 10.0, 10.0, 25.0, 25.0, 25.0, 50.0, 50.0),
        )
        self.assertEqual(module.seed_from_index(0)["seed"], "seed26723")
        self.assertEqual(module.seed_from_index(1)["seed"], "seed26737")
        with self.assertRaises(ValueError):
            module.seed_from_index(2)

    def test_frozen_sources_and_qm_contract_are_reused_not_rebuilt(self):
        module = load_module()
        text = PREPARE.read_text(encoding="utf-8")
        base = module.BASE
        self.assertIn("prepare_audit_nylc_a1_step1_pt2_cn_scout.py", text)
        self.assertEqual(base.FRAME_ROOT.name, "attempt_62112503")
        self.assertEqual(
            {source["candidate"] for source in base.SOURCES},
            {"seed26723_t378_f189", "seed26737_t676_f338"},
        )
        self.assertEqual(
            base.EXPECTED_PRMTOP_SHA256,
            "a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0",
        )
        self.assertEqual(base.EXPECTED_QM_ATOMS, 146)
        self.assertEqual(base.QMCHARGE, 0)
        self.assertEqual(base.EXPECTED_ELECTRONS, 510)
        self.assertEqual(base.EXPECTED_LINK_ATOMS, 6)
        self.assertNotIn("62021985", text)
        self.assertNotIn("tleap", text.lower())

    def test_only_attack_is_biased_and_outer_boundary_is_safe(self):
        module = load_module()
        restraint = module.attack_restraint(1.65, 50.0)
        self.assertIn("iat=8960,10287", restraint)
        self.assertIn("r4=4.500", restraint)
        for forbidden in ("8949,8961", "8961,10289", "10287,10289"):
            self.assertNotIn(forbidden, restraint)
        text = PREPARE.read_text(encoding="utf-8")
        self.assertIn("PT2_AND_CN_ARE_AUDIT_ONLY = True", text)
        self.assertIn("attack_angle_restrained", text)
        self.assertNotIn("target + 0.25", text)

    def test_runner_is_inherited_but_seeds_are_parallel(self):
        if not RUNNER.is_file() or not SLURM.is_file():
            self.skipTest("production runner absent at TDD RED")
        runner = RUNNER.read_text(encoding="utf-8")
        slurm = SLURM.read_text(encoding="utf-8")
        for token in (
            'CURRENT_RST7="$START_RST7"',
            'CURRENT_RST7="$NEXT_RST7"',
            "run_history.tsv",
            "run_history.jsonl",
            "CHAIN_MANIFEST.json",
            "SHA256.tsv",
        ):
            self.assertIn(token, runner)
        for token in (
            "#SBATCH --array=0-1%2",
            "#SBATCH -n 8",
            "#SBATCH --mem-per-cpu=2500M",
            "A1_ATTACK_CODE_SOURCE",
            "A1_ATTACK_GITHUB_COMMIT",
        ):
            self.assertIn(token, slurm)
        self.assertNotIn("#SBATCH --gres", slurm)

    def test_audit_separates_technical_restrained_and_scientific_gates(self):
        module = load_module()
        text = PREPARE.read_text(encoding="utf-8")
        for token in (
            "PASS_TECHNICAL_A1_ATTACK_INHERITED_CHAIN",
            "TETRAHEDRAL_LIKE_RESTRAINED",
            "FORCED_CLOSE_CONTACT",
            "c12_out_of_plane_A",
            "carbonyl_elongation_from_baseline_A",
            "reactant_guard",
            "NOT_EVALUATED_TS_PMF_BARRIER_MECHANISM",
            "DO_NOT_START_PMF",
        ):
            self.assertIn(token, text)


if __name__ == "__main__":
    unittest.main()
