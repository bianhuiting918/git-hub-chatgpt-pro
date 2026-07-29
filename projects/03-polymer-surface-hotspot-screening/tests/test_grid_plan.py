import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
PLANNER = PROJECT / "scripts" / "build_grid_plan.py"


def load_planner():
    spec = importlib.util.spec_from_file_location("build_grid_plan", PLANNER)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class GridPlanTests(unittest.TestCase):
    def test_small_receptor_uses_one_identical_geometry_grid(self):
        planner = load_planner()
        plan = planner.build_plan(
            bbox_min=[0.0, 0.0, 0.0],
            bbox_max=[30.0, 24.0, 18.0],
            spacing=0.75,
            margin=5.0,
            overlap=8.0,
            max_npts=126,
        )
        self.assertEqual(plan["status"], "GRID_PLAN_PASS")
        self.assertEqual(plan["channels"], ["A", "C", "OA", "HD"])
        self.assertEqual(len(plan["tiles"]), 1)
        tile = plan["tiles"][0]
        self.assertTrue(all(value % 2 == 0 for value in tile["npts"]))
        self.assertTrue(all(value <= 126 for value in tile["npts"]))
        self.assertEqual(
            {channel: tile["geometry_sha256"] for channel in plan["channels"]},
            tile["channel_geometry_sha256"],
        )

    def test_large_receptor_tiles_cover_margin_with_required_overlap(self):
        planner = load_planner()
        plan = planner.build_plan(
            bbox_min=[0.0, 0.0, 0.0],
            bbox_max=[190.0, 20.0, 20.0],
            spacing=0.75,
            margin=5.0,
            overlap=8.0,
            max_npts=126,
        )
        self.assertGreater(len(plan["tiles"]), 1)
        x_segments = sorted({(t["lower"][0], t["upper"][0]) for t in plan["tiles"]})
        self.assertLessEqual(x_segments[0][0], -5.0)
        self.assertGreaterEqual(x_segments[-1][1], 195.0)
        for left, right in zip(x_segments, x_segments[1:]):
            self.assertGreaterEqual(left[1] - right[0], 8.0 - 1e-6)
        for tile in plan["tiles"]:
            self.assertTrue(all(n <= 126 for n in tile["npts"]))

    def test_rejects_channel_geometry_mismatch(self):
        planner = load_planner()
        geometries = {
            "A": {"center": [0, 0, 0], "npts": [20, 20, 20], "spacing": 0.75},
            "C": {"center": [0, 0, 0], "npts": [20, 20, 20], "spacing": 0.75},
            "OA": {"center": [0, 0, 0], "npts": [22, 20, 20], "spacing": 0.75},
            "HD": {"center": [0, 0, 0], "npts": [20, 20, 20], "spacing": 0.75},
        }
        with self.assertRaisesRegex(ValueError, "nonidentical channel geometry"):
            planner.validate_channel_geometries(geometries)

    def test_cli_writes_deterministic_plan_from_receptor_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            metadata = root / "metadata.json"
            metadata.write_text(json.dumps({
                "status": "RECEPTOR_PASS",
                "prepared_bbox": {"min": [-10, -12, -14], "max": [35, 31, 29]},
                "prepared_pdbqt_sha256": "a" * 64,
            }), encoding="utf-8")
            outputs = []
            for suffix in ("a", "b"):
                output = root / f"plan_{suffix}.json"
                completed = subprocess.run([
                    sys.executable, str(PLANNER),
                    "--receptor-metadata", str(metadata),
                    "--output", str(output),
                    "--spacing", "0.75",
                    "--margin", "5.0",
                    "--overlap", "8.0",
                    "--max-npts", "126",
                ], capture_output=True, text=True)
                self.assertEqual(completed.returncode, 0, completed.stderr)
                outputs.append(output.read_bytes())
            self.assertEqual(outputs[0], outputs[1])


if __name__ == "__main__":
    unittest.main()
