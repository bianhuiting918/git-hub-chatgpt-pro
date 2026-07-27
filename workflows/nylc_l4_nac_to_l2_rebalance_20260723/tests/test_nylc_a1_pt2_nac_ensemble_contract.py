#!/usr/bin/env python3
import pathlib
import unittest

HERE = pathlib.Path(__file__).resolve().parents[1]
ANALYZER = HERE / "scripts" / "analyze_nylc_a1_pt2_nac_ensemble.py"
RUNNER = HERE / "scripts" / "run_nylc_a1_pt2_nac_ensemble.sh"
SLURM = HERE / "slurm" / "run_nylc_a1_pt2_nac_ensemble.sbatch"


class A1PT2NACEnsembleContract(unittest.TestCase):
    def test_files_exist(self):
        for path in (ANALYZER, RUNNER, SLURM):
            self.assertTrue(path.is_file(), path)

    def test_analyzer_freezes_denominator_and_atoms(self):
        text = ANALYZER.read_text(encoding="utf-8")
        for token in (
            "EXPECTED_REPLICA_DENOMINATOR = 9",
            "N_ALPHA = 8949",
            "THR267_OG1 = 8960",
            "TRANSFERRED_HG1 = 8961",
            "L2_C12 = 10287",
            "L2_O2 = 10288",
            "L2_N3 = 10289",
            "61976600",
            "generate.stdout",
            "replica_audit.json",
        ):
            self.assertIn(token, text)

    def test_joint_nac_and_pt2_gates_are_explicit(self):
        text = ANALYZER.read_text(encoding="utf-8")
        for token in (
            "NAC_DISTANCE_MAX_NM = 0.35",
            "NAC_ANGLE_MIN_DEG = 95.0",
            "NAC_ANGLE_MAX_DEG = 115.0",
            "PT2_DONOR_ACCEPTOR_MAX_NM = 0.35",
            "PT2_H_ACCEPTOR_MAX_NM = 0.25",
            "PT2_ANGLE_MIN_DEG = 135.0",
            "PASS_A1_PT2_PREORGANIZED_NAC_REPRODUCED",
            "NOT_SUPPORTED_A1_PT2_PREORGANIZED_NAC_NOT_REPRODUCED",
            "NOT_EVALUATED_PROTON_TRANSFER_TS_PMF_BARRIER_MECHANISM",
        ):
            self.assertIn(token, text)

    def test_source_hashes_and_all_replicas_are_fail_closed(self):
        text = ANALYZER.read_text(encoding="utf-8")
        for token in (
            "sha256",
            "tpr_sha256",
            "xtc_sha256",
            "actual_replica_denominator",
            "independent_positive_replica_count",
            "candidate_frames.tsv",
        ):
            self.assertIn(token, text)

    def test_runner_and_slurm_preserve_audit_boundaries(self):
        runner = RUNNER.read_text(encoding="utf-8")
        slurm = SLURM.read_text(encoding="utf-8")
        for token in (
            "run_history.tsv",
            "run_history.jsonl",
            "NOT_EVALUATED.json",
            "PASS.json",
            "scientific_status",
        ):
            self.assertIn(token, runner)
        for token in (
            "#SBATCH -n 4",
            "#SBATCH --mem-per-cpu=2500M",
            "#SBATCH -t 01:00:00",
            "A1_PT2_NAC_CODE_SOURCE",
            "A1_PT2_NAC_GITHUB_COMMIT",
        ):
            self.assertIn(token, slurm)
        self.assertNotIn("#SBATCH --gres", slurm)


if __name__ == "__main__":
    unittest.main()
