#!/usr/bin/env python3
import importlib.util
import pathlib
import unittest

FLOW = pathlib.Path(__file__).resolve().parents[1]
PREPARE = FLOW / "scripts" / "prepare_audit_nylc_a1_step1_pt2_cn_scout.py"
RUNNER = FLOW / "scripts" / "run_nylc_a1_step1_pt2_cn_scout.sh"
SLURM = FLOW / "slurm" / "run_nylc_a1_step1_pt2_cn_scout.sbatch"


class PT2CNScoutContract(unittest.TestCase):
    def test_minimal_files_exist(self):
        for path in (PREPARE, RUNNER, SLURM):
            self.assertTrue(path.is_file(), path)

    def test_frozen_sources_are_new_dual_frames_and_existing_qm(self):
        text = PREPARE.read_text(encoding="utf-8")
        for token in (
            "attempt_62112503",
            "seed26723_t378_f189",
            "seed26737_t676_f338",
            "14477791ce14a35cef0adf9b802b562e091660526ca06de1132f6a74070faf10",
            "a7924184ad3db4c13e0eab4929d46ee99621eacd523e899c0aca39c540350bc8",
            "attempt_62011285",
            "a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0",
            "qm_atom_count", "146", "qmcharge", "510", "link_atom_count", "6",
            "step1_qm_water_count",
        ):
            self.assertIn(token, text)
        self.assertNotIn("62021985", text)

    def test_grid_is_exact_seed_major_cartesian_product(self):
        spec = importlib.util.spec_from_file_location("_pt2_grid", PREPARE)
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        windows = [module.window_from_index(index) for index in range(16)]
        keys = {
            (
                item["seed"],
                item["attack_target_A"],
                item["pt2_q_target_A"],
                item["cn_target_A"],
            )
            for item in windows
        }
        self.assertEqual(len(keys), 16)
        self.assertEqual(sum(item["seed"] == "seed26723" for item in windows), 8)
        self.assertEqual(sum(item["seed"] == "seed26737" for item in windows), 8)
        self.assertEqual(module.ATTACK_TARGETS_A, (2.35, 1.85))
        self.assertEqual(module.PT2_TARGETS_A, ((1.10, 1.60), (1.35, 1.35)))
        self.assertEqual(module.CN_TARGETS_A, (1.40, 1.60))

    def test_restraint_and_audit_contract_is_bounded(self):
        text = PREPARE.read_text(encoding="utf-8")
        for token in (
            "THR267_N = 8949", "THR267_OG1 = 8960",
            "TRANSFERRED_HG1 = 8961", "L2_C12 = 10287",
            "L2_O2 = 10288", "L2_N3 = 10289",
            "GUIDE_STEPS = 200", "TARGET_STEPS = 600",
            "GUIDE_FORCE = 5.0", "TARGET_FORCE = 10.0",
            "q_attack", "q_pt2", "q_cn", "pt2_q_A",
            "FINAL RESULTS", "Run", "SANDER BOMB",
            "PASS_TECHNICAL_A1_PT2_CN_SCOUT",
            "NOT_EVALUATED_A1_PT2_CN_SCOUT",
        ):
            self.assertIn(token, text)

    def test_runner_and_slurm_are_reproducible_and_small(self):
        runner = RUNNER.read_text(encoding="utf-8")
        slurm = SLURM.read_text(encoding="utf-8")
        for token in ("run_history.tsv", "run_history.jsonl", "SHA256.tsv"):
            self.assertIn(token, runner)
        for token in (
            "#SBATCH --array=0-15%4", "#SBATCH -n 8",
            "#SBATCH --mem-per-cpu=2500M", "#SBATCH -t 01:00:00",
            "A1_PT2_CN_CODE_SOURCE", "A1_PT2_CN_GITHUB_COMMIT",
        ):
            self.assertIn(token, slurm)
        self.assertNotIn("#SBATCH --gres", slurm)


if __name__ == "__main__":
    unittest.main()
