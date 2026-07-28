#!/usr/bin/env python3
import pathlib
import unittest

HERE = pathlib.Path(__file__).resolve().parents[1]
ANALYZER = HERE / "scripts" / "extract_audit_nylc_a1_pt2_frames.py"
RUNNER = HERE / "scripts" / "run_extract_audit_nylc_a1_pt2_frames.sh"
SLURM = HERE / "slurm" / "run_extract_audit_nylc_a1_pt2_frames.sbatch"


class A1PT2FrameExtractionContract(unittest.TestCase):
    def test_files_exist(self):
        for path in (ANALYZER, RUNNER, SLURM):
            self.assertTrue(path.is_file(), path)

    def test_two_sources_are_frozen_independently(self):
        text = ANALYZER.read_text(encoding="utf-8")
        for token in (
            "seed26723", "378.0", "189",
            "c60078a92c2ace51facde4ef64e453f690177fc88b4d6363427f935944fa2e43",
            "1a54f1b5b9f139b746c22d9e0f7e9a4a94eb8154bf2b881888986eedca933d89",
            "seed26737", "676.0", "338",
            "dbd19a399547319d10630430ed494d33f6271bab0f30cb0de5a466c6af13ba20",
            "fcba14da98b331368061dcd990f2467628ad77b9b7a9c4ce88090f92e0831b05",
        ):
            self.assertIn(token, text)

    def test_extraction_is_atomic_and_non_overwriting(self):
        runner = RUNNER.read_text(encoding="utf-8")
        analyzer = ANALYZER.read_text(encoding="utf-8")
        for token in ("source.tmp.gro", "source.gro", "-pbc mol", "-ur compact", "-dump"):
            self.assertIn(token, runner)
        for token in ("exist_ok=False", "source.tmp.gro", "source.gro", "manifest.json"):
            self.assertIn(token, analyzer)

    def test_atom_bond_time_geometry_and_box_gates_are_fail_closed(self):
        text = ANALYZER.read_text(encoding="utf-8")
        for token in (
            "N_ALPHA = 8949", "THR267_OG1 = 8960", "TRANSFERRED_HG1 = 8961",
            "L2_C12 = 10287", "L2_O2 = 10288", "L2_N3 = 10289",
            "NAC_DISTANCE_MAX_NM = 0.35", "NAC_ANGLE_MIN_DEG = 95.0",
            "NAC_ANGLE_MAX_DEG = 115.0", "PT2_DONOR_ACCEPTOR_MAX_NM = 0.35",
            "PT2_H_ACCEPTOR_MAX_NM = 0.25", "PT2_ANGLE_MIN_DEG = 135.0",
            "time_tolerance_ps", "atom_count", "box_nm", "minimum_contact_nm",
            "nalpha_hg1", "og1_hg1_absent", "c12_n3", "sha256",
            "NOT_EVALUATED.json", "PASS.json",
        ):
            self.assertIn(token, text)

    def test_runner_records_history_and_slurm_is_small_cpu_only(self):
        runner = RUNNER.read_text(encoding="utf-8")
        slurm = SLURM.read_text(encoding="utf-8")
        for token in ("run_history.tsv", "run_history.jsonl", "SHA256.tsv"):
            self.assertIn(token, runner)
        for token in (
            "#SBATCH -n 4", "#SBATCH --mem-per-cpu=2500M", "#SBATCH -t 01:00:00",
            "A1_PT2_FRAME_CODE_SOURCE", "A1_PT2_FRAME_GITHUB_COMMIT",
        ):
            self.assertIn(token, slurm)
        self.assertNotIn("#SBATCH --gres", slurm)


if __name__ == "__main__":
    unittest.main()
