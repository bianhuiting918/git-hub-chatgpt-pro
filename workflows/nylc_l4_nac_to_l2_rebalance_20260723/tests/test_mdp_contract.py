#!/usr/bin/env python3
import unittest
from pathlib import Path

FLOW = Path(__file__).resolve().parents[1]
MDP = FLOW / "mdp"


def parse_mdp(name):
    values = {}
    for raw in (MDP / name).read_text(encoding="utf-8").splitlines():
        body = raw.split(";", 1)[0].strip()
        if not body or "=" not in body:
            continue
        key, value = body.split("=", 1)
        values[key.strip().lower()] = value.strip()
    return values


class MdpContractTests(unittest.TestCase):
    def test_stage_schedule_and_restraint_contract(self):
        expected = {
            "em.mdp": ("steep", 50000, "POSRES_L2_1000"),
            "nvt50.mdp": ("md", 50000, "POSRES_L2_1000"),
            "nvt150.mdp": ("md", 50000, "POSRES_L2_500"),
            "nvt300.mdp": ("md", 100000, "POSRES_L2_100"),
            "npt300r.mdp": ("md", 100000, "POSRES_L2_100"),
            "npt300rel.mdp": ("md", 100000, "POSRES_L2_10"),
        }
        for name, (integrator, nsteps, restraint) in expected.items():
            values = parse_mdp(name)
            self.assertEqual(values["integrator"].lower(), integrator)
            self.assertEqual(int(values["nsteps"]), nsteps)
            self.assertIn(restraint, values.get("define", ""))

    def test_temperature_ramp_contract(self):
        for name, temperature in (
            ("nvt50.mdp", 50),
            ("nvt150.mdp", 150),
            ("nvt300.mdp", 300),
            ("npt300r.mdp", 300),
            ("npt300rel.mdp", 300),
        ):
            values = parse_mdp(name)
            self.assertEqual(float(values["ref-t"]), temperature)
        self.assertEqual(parse_mdp("nvt50.mdp")["gen-vel"].lower(), "yes")
        self.assertEqual(parse_mdp("nvt150.mdp")["gen-vel"].lower(), "no")

    def test_final_window_is_one_ns_fully_unrestrained_npt(self):
        values = parse_mdp("npt300free.mdp")
        self.assertEqual(float(values["dt"]) * int(values["nsteps"]), 1000.0)
        self.assertNotIn("define", values)
        self.assertEqual(values["pcoupl"].lower(), "parrinello-rahman")
        self.assertEqual(float(values["ref-t"]), 300.0)
        self.assertEqual(float(values["ref-p"]), 1.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
