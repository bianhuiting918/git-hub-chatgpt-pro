#!/usr/bin/env python3
import importlib.util
import itertools
import math
import pathlib
import unittest

SCRIPT = (
    pathlib.Path(__file__).resolve().parents[1]
    / "scripts/prepare_audit_nylc_a1_step2_qmwater_endpoint.py"
)
SPEC = importlib.util.spec_from_file_location("nylc_a1_step2_qmwater", SCRIPT)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load {SCRIPT}")
DRIVER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DRIVER)


def explicit_minimum_image(left, right, box):
    a, b, c, alpha_deg, beta_deg, gamma_deg = map(float, box)
    alpha, beta, gamma = map(math.radians, (alpha_deg, beta_deg, gamma_deg))
    va = (a, 0.0, 0.0)
    vb = (b * math.cos(gamma), b * math.sin(gamma), 0.0)
    cx = c * math.cos(beta)
    cy = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
    vc = (cx, cy, math.sqrt(c * c - cx * cx - cy * cy))
    displacement = tuple(x - y for x, y in zip(left, right))
    candidates = []
    for shift in itertools.product(range(-2, 3), repeat=3):
        vector = tuple(
            displacement[axis]
            - shift[0] * va[axis]
            - shift[1] * vb[axis]
            - shift[2] * vc[axis]
            for axis in range(3)
        )
        candidates.append((sum(value * value for value in vector), shift, vector))
    return min(candidates)[2]


class TriclinicMinimumImageContract(unittest.TestCase):
    def test_orthorhombic_behavior_is_preserved(self):
        cell = DRIVER._periodic_cell((10.0, 20.0, 30.0, 90.0, 90.0, 90.0))
        observed = DRIVER._delta((9.8, 19.7, 0.2), (0.2, 0.3, 29.7), cell)
        self.assertEqual(len(observed), 3)
        for actual, expected in zip(observed, (-0.4, -0.6, 0.5)):
            self.assertAlmostEqual(actual, expected, places=10)

    def test_scnet_60_60_90_boxes_match_explicit_lattice_enumeration(self):
        for box in (
            (123.7404, 123.7404, 123.7403877, 60.0000502, 60.0000502, 90.0),
            (123.8014, 123.8014, 123.801464, 60.0000171, 60.0000171, 90.0),
        ):
            with self.subTest(box=box):
                left = (119.0, 112.0, 109.0)
                right = (4.0, 8.0, 7.0)
                observed = DRIVER._delta(left, right, DRIVER._periodic_cell(box))
                expected = explicit_minimum_image(left, right, box)
                for actual, reference in zip(observed, expected):
                    self.assertAlmostEqual(actual, reference, places=9)

    def test_invalid_or_degenerate_boxes_are_rejected(self):
        invalid = (
            (0.0, 10.0, 10.0, 90.0, 90.0, 90.0),
            (10.0, 10.0, 10.0, 0.0, 90.0, 90.0),
            (10.0, 10.0, 10.0, 90.0, 90.0, 180.0),
            (10.0, 10.0, float("nan"), 90.0, 90.0, 90.0),
        )
        for box in invalid:
            with self.subTest(box=box):
                with self.assertRaises(ValueError):
                    DRIVER._periodic_cell(box)


if __name__ == "__main__":
    unittest.main()
