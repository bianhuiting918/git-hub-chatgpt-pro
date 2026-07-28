#!/usr/bin/env python3
"""Contract for the NylC A1 Step1 acyl-endpoint stability fallback.

This is intentionally RED until the three production files named below exist.
It checks externally observable configuration, release gating, and provenance
rather than prescribing the implementation's internal control flow.
"""

import importlib.util
import pathlib
import unittest


FLOW = pathlib.Path(__file__).resolve().parents[1]
PREPARE = FLOW / "scripts" / "prepare_audit_nylc_a1_step1_acyl_endpoint_stability.py"
RUNNER = FLOW / "scripts" / "run_nylc_a1_step1_acyl_endpoint_stability.sh"
SLURM = FLOW / "slurm" / "run_nylc_a1_step1_acyl_endpoint_stability.sbatch"


def load_module():
    if not PREPARE.is_file():
        raise unittest.SkipTest("production driver absent at TDD RED")
    spec = importlib.util.spec_from_file_location("_a1_acyl_endpoint", PREPARE)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load acyl-endpoint preparation driver")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class A1AcylEndpointStabilityContract(unittest.TestCase):
    def test_minimal_production_files_exist(self):
        for path in (PREPARE, RUNNER, SLURM):
            self.assertTrue(path.is_file(), path)

    def test_two_frozen_seeds_and_qmmm_identity_are_reused(self):
        module = load_module()
        text = PREPARE.read_text(encoding="utf-8")
        base = module.BASE
        self.assertEqual(module.seed_from_index(0)["seed"], "seed26723")
        self.assertEqual(module.seed_from_index(1)["seed"], "seed26737")
        with self.assertRaises(ValueError):
            module.seed_from_index(2)
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
        self.assertTrue(getattr(module, "NO_QM_WATER", False))
        self.assertNotIn("62021985", text)
        self.assertNotIn("tleap", text.lower())

    def test_reactant_n3_graph_and_reactive_atom_identity_are_guarded(self):
        module = load_module()
        self.assertEqual(
            module.REACTIVE_ATOMS,
            {
                "nalpha": 8949,
                "og1": 8960,
                "hg1": 8961,
                "c11": 10286,
                "c12": 10287,
                "o2": 10288,
                "n3": 10289,
            },
        )
        guard = module.validate_n3_reactant_graph
        self.assertTrue(guard({10287: "C", 10286: "C", 8949: "H"}))
        for invalid in (
            {10287: "C", 10286: "C"},
            {10287: "C", 10286: "C", 8949: "N"},
            {10286: "C", 8949: "H", 9000: "C"},
        ):
            with self.subTest(neighbors=invalid):
                self.assertFalse(guard(invalid))

    def test_construction_targets_and_release_lengths_are_exact(self):
        module = load_module()
        self.assertEqual(module.INTERMEDIATE_TARGETS_A,
                         {"attack": 1.60, "hg1_n3": 1.40, "cn": 1.75})
        self.assertEqual(module.PRODUCT_TARGETS_A,
                         {"attack": 1.50, "hg1_n3": 1.05, "cn": 2.20})
        self.assertEqual(module.BUILD_FORCE_KCAL_MOL_A2, (25.0, 50.0))
        self.assertEqual(module.INTERMEDIATE_MINIMIZATION, (400, 100))
        self.assertEqual(module.PRODUCT_MINIMIZATION, (800, 200))
        self.assertEqual(module.LOCAL_RELEASE_MINIMIZATION, (800, 200))
        self.assertEqual(module.FULL_RELEASE_MINIMIZATION, (800, 200))
        self.assertEqual(module.LOCAL_RELEASE_STEPS, 800)
        self.assertEqual(module.FULL_RELEASE_STEPS, 800)
        self.assertEqual(module.RELEASE_MD_STEPS, 500)
        self.assertEqual(module.RELEASE_DT_PS, 0.0005)
        self.assertEqual(module.RELEASE_TEMP_K, 300.0)
        self.assertEqual(module.RELEASE_NTWX, 10)

    def test_only_three_reactive_distances_are_biased_and_product_is_unforced(self):
        module = load_module()
        for target, force in (
            (module.INTERMEDIATE_TARGETS_A, 25.0),
            (module.PRODUCT_TARGETS_A, 50.0),
        ):
            restraints = module.reactive_restraints(target, force)
            self.assertEqual(len(restraints), 3)
            joined = "\n".join(restraints)
            for atoms in ("8960,10287", "8961,10289", "10287,10289"):
                self.assertIn(atoms, joined)
            for forbidden in ("8949,8961", "10287,10288", "iat=10287,10288"):
                self.assertNotIn(forbidden, joined)
        text = PREPARE.read_text(encoding="utf-8").lower()
        self.assertNotIn("improper", text)
        self.assertNotIn("reaction angle", text)

    def test_product_gate_requires_release_stable_acyl_chemistry(self):
        module = load_module()
        self.assertEqual(module.PRODUCT_GATE["attack_A"], (1.40, 1.65))
        self.assertEqual(module.PRODUCT_GATE["hg1_n3_A"], (0.95, 1.20))
        self.assertEqual(module.PRODUCT_GATE["nalpha_hg1_A_min"], 1.55)
        self.assertEqual(module.PRODUCT_GATE["qpt_A_min"], 0.45)
        self.assertEqual(module.PRODUCT_GATE["c12_n3_A_min"], 2.05)
        self.assertEqual(module.PRODUCT_GATE["c12_o2_A"], (1.18, 1.30))
        self.assertEqual(module.PRODUCT_GATE["product_out_of_plane_A_max"], 0.12)
        self.assertEqual(module.PRODUCT_GATE["product_angle_sum_deg_min"], 350.0)
        self.assertEqual(module.PRODUCT_GATE["last_fraction"], 0.40)
        self.assertEqual(module.PRODUCT_GATE["occupancy_min"], 0.80)
        product_like = {
            "attack_A": 1.50, "hg1_n3_A": 1.05, "nalpha_hg1_A": 1.60,
            "qPT_A": 0.50, "c12_n3_A": 2.10, "c12_o2_A": 1.24,
            "product_out_of_plane_A": 0.10, "product_angle_sum_deg": 355.0,
            "hg1_nearest_qm_heavy_atom": 10289,
        }
        self.assertTrue(module.is_product_like(product_like))
        product_like["hg1_nearest_qm_heavy_atom"] = 8960
        self.assertFalse(module.is_product_like(product_like))

    def test_runner_stages_scratch_trajectory_and_failure_audit(self):
        if not RUNNER.is_file():
            self.skipTest("production runner absent at TDD RED")
        runner = RUNNER.read_text(encoding="utf-8")
        expected_stages = ("intermediate", "product", "local_release", "full_release", "release_md")
        self.assertLess(
            max(runner.index(stage) for stage in expected_stages[:-1]),
            runner.index(expected_stages[-1]),
        )
        for stage in expected_stages:
            self.assertIn(stage, runner)
        for token in (
            "SLURM_TMPDIR",
            "run_history.tsv",
            "run_history.jsonl",
            "ENDPOINT_MANIFEST.json",
            "RESULT.json",
            "NOT_EVALUATED_A1_ACYL_ENDPOINT_STABILITY",
            "trap",
            "exit 1",
            "flock",
        ):
            self.assertIn(token, runner)
        self.assertNotIn("cp \"$RELEASE_TRAJ\"", runner)

    def test_seed_and_cross_seed_classifications_preserve_pmf_boundary(self):
        module = load_module()
        text = PREPARE.read_text(encoding="utf-8")
        expected_classes = {
            "PERSISTS_ACYL_PRODUCT",
            "RETURNS_TETRAHEDRAL",
            "RETURNS_REACTANT",
            "ZWITTERIONIC_CLEAVAGE",
            "MISROUTED_PROTON",
            "RESTRAINT_DEPENDENT_PRODUCT",
            "NOT_EVALUATED_TECHNICAL_FAILURE",
        }
        self.assertEqual(set(module.SEED_CLASSIFICATIONS), expected_classes)
        for token in (
            "PASS_ACYL_PRODUCT_ENDPOINT_REPRODUCED",
            "NOT_REPRODUCED_A1_ACYL_PRODUCT_ENDPOINT",
            "FAIL_NO_RELEASE_STABLE_A1_ACYL_PRODUCT_ENDPOINT",
            "NOT_EVALUATED_TECHNICAL_A1_ACYL_PRODUCT_ENDPOINT",
            "NOT_EVALUATED_TS_PMF_BARRIER_MECHANISM",
            "DO_NOT_START_PMF",
        ):
            self.assertIn(token, text)
        self.assertEqual(module.cross_seed_status(
            ["PERSISTS_ACYL_PRODUCT", "PERSISTS_ACYL_PRODUCT"]), 
            "PASS_ACYL_PRODUCT_ENDPOINT_REPRODUCED")
        self.assertEqual(module.cross_seed_status(
            ["PERSISTS_ACYL_PRODUCT", "RETURNS_REACTANT"]),
            "NOT_REPRODUCED_A1_ACYL_PRODUCT_ENDPOINT")
        self.assertEqual(module.cross_seed_status(
            ["RETURNS_REACTANT", "MISROUTED_PROTON"]),
            "FAIL_NO_RELEASE_STABLE_A1_ACYL_PRODUCT_ENDPOINT")
        self.assertEqual(module.cross_seed_status(
            ["NOT_EVALUATED_TECHNICAL_FAILURE", "PERSISTS_ACYL_PRODUCT"]),
            "NOT_EVALUATED_TECHNICAL_A1_ACYL_PRODUCT_ENDPOINT")

    def test_slurm_two_seed_resources_and_immutable_snapshot_are_required(self):
        if not SLURM.is_file():
            self.skipTest("production Slurm wrapper absent at TDD RED")
        slurm = SLURM.read_text(encoding="utf-8")
        for token in (
            "#SBATCH --array=0-1%2",
            "#SBATCH -n 8",
            "#SBATCH --time=08:00:00",
            "#SBATCH --mem-per-cpu=2500M",
            "A1_ACYL_ENDPOINT_CODE_SOURCE",
            "A1_ACYL_ENDPOINT_GITHUB_COMMIT",
            "GITHUB_COMMIT",
            "SNAPSHOT_SHA256.tsv",
            "sha256sum -c",
        ):
            self.assertIn(token, slurm)
        self.assertNotIn("#SBATCH --gres", slurm)
        self.assertNotIn("#SBATCH --gpus", slurm)

    def test_candidate_attempts_are_not_overwritten_and_minimum_audit_persists(self):
        if not PREPARE.is_file() or not RUNNER.is_file():
            self.skipTest("production files absent at TDD RED")
        text = PREPARE.read_text(encoding="utf-8") + RUNNER.read_text(encoding="utf-8")
        for token in (
            "ENDPOINT_MANIFEST.json",
            "SHA256.tsv",
            "run_history.tsv",
            "run_history.jsonl",
            "constructed_product.rst7",
            "released_endpoint.rst7",
        ):
            self.assertIn(token, text)
        self.assertIn("already exists", text.lower())
        self.assertNotIn("rm -rf", text)


if __name__ == "__main__":
    unittest.main()
