#!/usr/bin/env python3
"""Contract for the approved NylC A1 bidirectional Step1 follow-up."""
from __future__ import annotations

import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"


def load(name: str):
    path = SCRIPTS / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class BidirectionalFollowupContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reverse = load("prepare_audit_nylc_a1_step1_reverse_inherited_followup.py")
        cls.forward = load("prepare_audit_nylc_a1_step1_forward_force_calibration.py")

    def test_reverse_sources_and_serial_inheritance_are_frozen(self):
        desc = self.reverse.describe()
        self.assertEqual(desc["array_task_count"], 2)
        self.assertEqual(desc["max_windows_per_chain"], 8)
        self.assertTrue(desc["strict_serial_restart_inheritance"])
        self.assertFalse(desc["failed_restart_inheritance"])
        self.assertEqual(
            [self.reverse.task_spec(i)["source_task_index"] for i in range(2)],
            [6, 8],
        )
        self.assertEqual(
            [self.reverse.task_spec(i)["force_scale"] for i in range(2)],
            [1.0, 4.0],
        )
        self.assertEqual(
            [self.reverse.task_spec(i)["route_label"] for i in range(2)],
            ["LOW_BIAS", "TARGET_CLOSE"],
        )
        self.assertEqual(
            [self.reverse.task_spec(i)["source_restart_sha256"] for i in range(2)],
            [
                "54ecabec3cbb993a81b36ce86630355e2143797481d71726ca83edd4d65862df",
                "2f30ffd0746cc7c53f36ac102b446256ea575e24a2d5bfa9da15a1455db96a7f",
            ],
        )
        self.assertEqual(desc["active_groups"], ["C12-N3", "C12-O2", "qPT"])
        self.assertEqual(desc["monitored_unrestrained_groups"], ["OG1-C12"])
        self.assertEqual(
            desc["per_window_deltas_A"],
            {"cn": -0.08, "carbonyl": 0.03, "nalpha_hg1": -0.05, "hg1_n3": 0.05},
        )

    def test_forward_matrix_and_mechanisms_are_frozen(self):
        desc = self.forward.describe()
        self.assertEqual(desc["array_task_count"], 12)
        self.assertEqual(desc["force_scales"], [1.0, 2.0, 4.0])
        self.assertEqual(
            desc["mechanisms"],
            ["ADDITION_FIRST_FORWARD", "FULLY_CONCERTED_FORWARD"],
        )
        self.assertEqual(
            desc["active_groups"]["ADDITION_FIRST_FORWARD"],
            ["OG1-C12", "C12-O2"],
        )
        self.assertEqual(
            desc["monitored_unrestrained_groups"]["ADDITION_FIRST_FORWARD"],
            ["C12-N3", "qPT"],
        )
        self.assertEqual(
            desc["active_groups"]["FULLY_CONCERTED_FORWARD"],
            ["OG1-C12", "C12-N3", "C12-O2", "qPT"],
        )
        self.assertEqual(
            desc["first_window_deltas_A"],
            {
                "attack": -0.04,
                "cn": 0.08,
                "carbonyl": 0.03,
                "nalpha_hg1": 0.05,
                "hg1_n3": -0.05,
            },
        )
        specs = [self.forward.task_spec(i) for i in range(12)]
        self.assertEqual(
            {(item["seed"], item["mechanism"], item["force_scale"]) for item in specs},
            {
                (seed, mechanism, scale)
                for seed in ("seed26723", "seed26737")
                for mechanism in (
                    "ADDITION_FIRST_FORWARD",
                    "FULLY_CONCERTED_FORWARD",
                )
                for scale in (1.0, 2.0, 4.0)
            },
        )

    def test_frozen_qmmm_contract_and_no_automatic_downstream(self):
        for module in (self.reverse, self.forward):
            desc = module.describe()
            contract = desc["frozen_step1_contract"]
            self.assertEqual(contract["qm_atom_count"], 146)
            self.assertEqual(contract["qmcharge"], 0)
            self.assertEqual(contract["electron_count_including_link_h"], 510)
            self.assertEqual(contract["link_atom_count"], 6)
            self.assertEqual(contract["step1_qm_water_count"], 0)
            self.assertEqual(desc["automatic_downstream_action"], "NONE")
            self.assertTrue(desc["restrained_structures_are_not_ts"])


if __name__ == "__main__":
    unittest.main()
