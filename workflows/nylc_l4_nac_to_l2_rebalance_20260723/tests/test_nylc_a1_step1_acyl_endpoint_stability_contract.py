#!/usr/bin/env python3
"""Behavioral contract for the NylC A1 Step1 acyl-endpoint fallback.

The module stays RED until all three production files exist.  Once present,
the test exercises release, product persistence, provenance, and scheduler
behavior through the driver's small observable helper interfaces.
"""

import importlib.util
import pathlib
import tempfile
import unittest
from unittest import mock


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
    def test_exact_minimal_production_surface_exists(self):
        missing = [str(path) for path in (PREPARE, RUNNER, SLURM) if not path.is_file()]
        self.assertFalse(missing, "missing required production files: " + ", ".join(missing))

    def test_frozen_authority_is_delegated_to_base_without_invented_source_values(self):
        module = load_module()
        base = module.BASE
        self.assertEqual(
            base.EXPECTED_PRMTOP_SHA256,
            "a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0",
        )
        self.assertEqual(
            {source["candidate"] for source in base.SOURCES},
            {"seed26723_t378_f189", "seed26737_t676_f338"},
        )
        self.assertEqual(base.EXPECTED_QM_ATOMS, 146)
        self.assertEqual(base.QMCHARGE, 0)
        self.assertEqual(base.EXPECTED_ELECTRONS, 510)
        self.assertEqual(base.EXPECTED_LINK_ATOMS, 6)
        self.assertTrue(module.NO_QM_WATER)
        source = base.SOURCES[0]
        with mock.patch.object(base, "validate_authority", return_value={"validated": True}) as validate:
            self.assertEqual(module.validate_source_authority(source), {"validated": True})
        validate.assert_called_once_with(source)

    def test_reactant_n3_graph_and_atoms_are_hard_preconditions(self):
        module = load_module()
        self.assertEqual(
            module.REACTIVE_ATOMS,
            {"nalpha": 8949, "og1": 8960, "hg1": 8961, "c11": 10286,
             "c12": 10287, "o2": 10288, "n3": 10289},
        )
        guard = module.validate_n3_reactant_graph
        self.assertTrue(guard({10287: "C", 10286: "C", 8949: "H"}))
        for neighbors in (
            {10287: "C", 10286: "C"},
            {10287: "C", 10286: "C", 8949: "N"},
            {10286: "C", 8949: "H", 9000: "C"},
        ):
            with self.subTest(neighbors=neighbors):
                self.assertFalse(guard(neighbors))

    def test_construction_and_release_stage_specs_remove_reactive_bias(self):
        module = load_module()
        self.assertEqual(module.INTERMEDIATE_TARGETS_A,
                         {"attack": 1.60, "hg1_n3": 1.40, "cn": 1.75})
        self.assertEqual(module.PRODUCT_TARGETS_A,
                         {"attack": 1.50, "hg1_n3": 1.05, "cn": 2.20})
        self.assertEqual(module.BUILD_FORCE_KCAL_MOL_A2, (25.0, 50.0))
        for stage, expected_cycles in (
            ("intermediate", (400, 100)),
            ("product", (800, 200)),
            ("local_release", (800, 200)),
            ("full_release", (800, 200)),
        ):
            with self.subTest(stage=stage):
                spec = module.stage_spec(stage, "seed26723")
                self.assertEqual((spec["maxcyc"], spec["ncyc"]), expected_cycles)
        for stage in ("intermediate", "product"):
            with self.subTest(stage=stage):
                self.assertEqual(
                    set(module.stage_spec(stage, "seed26723")["reactive_restraints"]),
                    {"attack", "hg1_n3", "cn"},
                )
        local = module.stage_spec("local_release", "seed26723")
        self.assertEqual(local["reactive_restraints"], ())
        self.assertTrue(local["environment_restraint"])
        full = module.stage_spec("full_release", "seed26723")
        self.assertEqual(full["reactive_restraints"], ())
        self.assertEqual(full["ntr"], 0)
        self.assertEqual(full["nmropt"], 0)
        self.assertIsNone(full["disang"])

    def test_release_md_accepts_only_valid_full_release_and_locks_seed_controls(self):
        module = load_module()
        valid_full_release = {
            "stage": "full_release", "technical_pass": True,
            "restart_path": "full_release.rst7", "restart_sha256": "a" * 64,
        }
        for seed, expected_ig in (("seed26723", 26723), ("seed26737", 26737)):
            with self.subTest(seed=seed):
                md = module.release_md_spec(seed, valid_full_release)
                self.assertEqual(md["input_restart_sha256"], valid_full_release["restart_sha256"])
                self.assertEqual(md["ig"], expected_ig)
                self.assertEqual(md["nstlim"], 500)
                self.assertEqual(md["dt"], 0.0005)
                self.assertEqual(md["temp0"], 300.0)
                self.assertEqual(md["ntwx"], 10)
                self.assertEqual((md["ntt"], md["gamma_ln"], md["ntc"], md["ntf"], md["ntpr"]),
                                 (3, 2.0, 1, 1, 10))
        for invalid in (
            {"stage": "product", "technical_pass": True, "restart_path": "x", "restart_sha256": "a" * 64},
            {"stage": "full_release", "technical_pass": False, "restart_path": "x", "restart_sha256": "a" * 64},
            {"stage": "full_release", "technical_pass": True, "restart_path": "", "restart_sha256": ""},
        ):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    module.release_md_spec("seed26723", invalid)

    def test_product_geometry_uses_o2_og1_c11_plane_and_all_gate_boundaries(self):
        module = load_module()
        coordinates = {
            10288: (0.0, 0.0, 0.0), 8960: (1.0, 0.0, 0.0),
            10286: (0.0, 1.0, 0.0), 10287: (0.2, 0.3, 0.10),
        }
        measured = module.measure_product_geometry(coordinates)
        self.assertAlmostEqual(measured["product_out_of_plane_A"], 0.10, places=7)
        good = {
            "attack_A": 1.50, "hg1_n3_A": 1.05, "nalpha_hg1_A": 1.55,
            "qPT_A": 0.45, "c12_n3_A": 2.05, "c12_o2_A": 1.24,
            "product_out_of_plane_A": 0.12, "product_angle_sum_deg": 350.0,
            "hg1_nearest_qm_heavy_atom": 10289,
        }
        self.assertTrue(module.is_product_like(good))
        negatives = {
            "attack_low": {"attack_A": 1.399},
            "attack_high": {"attack_A": 1.651},
            "hg1_n3_low": {"hg1_n3_A": 0.949},
            "hg1_n3_high": {"hg1_n3_A": 1.201},
            "nalpha_hg1": {"nalpha_hg1_A": 1.549},
            "qpt": {"qPT_A": 0.449},
            "cn": {"c12_n3_A": 2.049},
            "carbonyl_low": {"c12_o2_A": 1.179},
            "carbonyl_high": {"c12_o2_A": 1.301},
            "oop": {"product_out_of_plane_A": 0.121},
            "angle_sum": {"product_angle_sum_deg": 349.999},
            "nearest_qm_heavy": {"hg1_nearest_qm_heavy_atom": 8960},
        }
        for name, update in negatives.items():
            with self.subTest(gate=name):
                frame = dict(good)
                frame.update(update)
                self.assertFalse(module.is_product_like(frame))

    def test_persistence_uses_only_final_40_percent_and_requires_80_percent_product(self):
        module = load_module()
        good = {
            "attack_A": 1.50, "hg1_n3_A": 1.05, "nalpha_hg1_A": 1.55,
            "qPT_A": 0.45, "c12_n3_A": 2.05, "c12_o2_A": 1.24,
            "product_out_of_plane_A": 0.12, "product_angle_sum_deg": 350.0,
            "hg1_nearest_qm_heavy_atom": 10289,
        }
        bad = dict(good, qPT_A=0.44)
        self.assertTrue(module.release_md_persists([bad] * 6 + [good] * 4))
        self.assertFalse(module.release_md_persists([good] * 6 + [good] * 3 + [bad]))
        self.assertFalse(module.release_md_persists([good] * 4 + [bad] * 6))

    def test_classification_and_pmf_boundary_remain_nonadvancing(self):
        module = load_module()
        self.assertEqual(
            set(module.SEED_CLASSIFICATIONS),
            {"PERSISTS_ACYL_PRODUCT", "RETURNS_TETRAHEDRAL", "RETURNS_REACTANT",
             "ZWITTERIONIC_CLEAVAGE", "MISROUTED_PROTON",
             "RESTRAINT_DEPENDENT_PRODUCT", "NOT_EVALUATED_TECHNICAL_FAILURE"},
        )
        for outcomes, expected in (
            (["PERSISTS_ACYL_PRODUCT"] * 2, "PASS_ACYL_PRODUCT_ENDPOINT_REPRODUCED"),
            (["PERSISTS_ACYL_PRODUCT", "RETURNS_REACTANT"],
             "NOT_REPRODUCED_A1_ACYL_PRODUCT_ENDPOINT"),
            (["RETURNS_REACTANT", "MISROUTED_PROTON"],
             "FAIL_NO_RELEASE_STABLE_A1_ACYL_PRODUCT_ENDPOINT"),
            (["NOT_EVALUATED_TECHNICAL_FAILURE", "PERSISTS_ACYL_PRODUCT"],
             "NOT_EVALUATED_TECHNICAL_A1_ACYL_PRODUCT_ENDPOINT"),
        ):
            with self.subTest(outcomes=outcomes):
                self.assertEqual(module.cross_seed_status(outcomes), expected)
        self.assertFalse(module.AUTO_START_PATH_SAMPLING)
        self.assertIn("NOT_EVALUATED_TS_PMF_BARRIER_MECHANISM", module.TERMINAL_BOUNDARIES)
        self.assertIn("DO_NOT_START_PMF", module.TERMINAL_BOUNDARIES)

    def test_persistence_policy_prevents_trajectory_leaks_and_candidate_overwrite(self):
        module = load_module()
        allowed = {
            "ENDPOINT_MANIFEST.json", "RESULT.json", "PASS.json", "NOT_EVALUATED.json",
            "SHA256.tsv", "run_history.tsv", "run_history.jsonl",
            "constructed_product.rst7", "released_endpoint.rst7",
        }
        for name in allowed:
            with self.subTest(allowed=name):
                self.assertTrue(module.is_persistable_artifact(name))
        for name in ("release.nc", "release.mdcrd", "release.dcd", "checkpoint.nc",
                     "system.prmtop", "frame.xtc"):
            with self.subTest(rejected=name):
                self.assertFalse(module.is_persistable_artifact(name))
        with tempfile.TemporaryDirectory() as tmp:
            candidate = pathlib.Path(tmp) / "candidate"
            candidate.mkdir()
            with self.assertRaises(FileExistsError):
                module.require_fresh_candidate_dir(candidate)

    def test_runner_and_slurm_enforce_operational_safety_and_snapshot_coverage(self):
        if not RUNNER.is_file() or not SLURM.is_file():
            self.skipTest("production runner or Slurm wrapper absent at TDD RED")
        runner = RUNNER.read_text(encoding="utf-8")
        slurm = SLURM.read_text(encoding="utf-8")
        self.assertIn("SLURM_TMPDIR", runner)
        for artifact in ("ENDPOINT_MANIFEST.json", "RESULT.json", "SHA256.tsv",
                         "run_history.tsv", "run_history.jsonl"):
            self.assertIn(artifact, runner)
        self.assertIn("flock", runner)
        self.assertIn("trap", runner)
        self.assertNotIn("rm -rf", runner)
        for token in (
            "#SBATCH --array=0-1%2", "#SBATCH -N 1", "#SBATCH -n 8",
            "#SBATCH --time=08:00:00", "#SBATCH --mem-per-cpu=2500M",
            "A1_ACYL_ENDPOINT_CODE_SOURCE", "A1_ACYL_ENDPOINT_GITHUB_COMMIT",
            "GITHUB_COMMIT", "SNAPSHOT_SHA256.tsv", "sha256sum -c",
            PREPARE.name, RUNNER.name, SLURM.name,
            "prepare_audit_nylc_a1_step1_pt2_cn_scout.py",
        ):
            self.assertIn(token, slurm)
        self.assertNotIn("#SBATCH --gres", slurm)
        self.assertNotIn("#SBATCH --gpus", slurm)
        invocation_text = (PREPARE.read_text(encoding="utf-8") + runner + slurm).lower()
        self.assertNotRegex(invocation_text, r"\b(?:neb|string|umbrella)\b.*(?:sander|sbatch)")
        self.assertFalse(module.AUTO_START_PATH_SAMPLING)


if __name__ == "__main__":
    unittest.main()
