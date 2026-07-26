import importlib.util
import tempfile
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
RUNNER = PROJECT / "scripts" / "run_grid_maps.py"


def load_runner():
    spec = importlib.util.spec_from_file_location("run_grid_maps", RUNNER)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def write_map(path, spacing=0.75, npts=(20, 22, 24), center=(1.0, 2.0, 3.0)):
    path.write_text(
        "GRID_PARAMETER_FILE tile.gpf\n"
        "GRID_DATA_FILE tile.maps.fld\n"
        "MACROMOLECULE receptor.pdbqt\n"
        f"SPACING {spacing:.3f}\n"
        f"NELEMENTS {npts[0]} {npts[1]} {npts[2]}\n"
        f"CENTER {center[0]:.3f} {center[1]:.3f} {center[2]:.3f}\n"
        "0\n", encoding="ascii"
    )


class GridMapTests(unittest.TestCase):
    def test_gpf_requests_all_four_channels_with_one_geometry(self):
        runner = load_runner()
        text = runner.render_gpf(
            receptor_name="receptor.pdbqt",
            parameter_name="AD4.1_bound.dat",
            prefix="tile",
            receptor_types=["A", "C", "HD", "N", "OA", "SA"],
            center=[1.0, 2.0, 3.0],
            npts=[20, 22, 24],
            spacing=0.75,
        )
        self.assertIn("ligand_types A C OA HD", text)
        for channel in ("A", "C", "OA", "HD"):
            self.assertIn(f"map tile.{channel}.map", text)
        self.assertIn("npts 20 22 24", text)
        self.assertIn("spacing 0.750000", text)
        self.assertIn("gridcenter 1.000000 2.000000 3.000000", text)

    def test_map_header_validation_accepts_identical_geometry(self):
        runner = load_runner()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = {}
            for channel in ("A", "C", "OA", "HD"):
                path = root / f"tile.{channel}.map"
                write_map(path)
                paths[channel] = path
            geometry = runner.validate_map_set(paths)
            self.assertEqual(geometry["npts"], [20, 22, 24])
            self.assertEqual(geometry["spacing"], 0.75)
            self.assertEqual(geometry["center"], [1.0, 2.0, 3.0])

    def test_map_header_validation_rejects_mismatch_and_missing(self):
        runner = load_runner()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = {}
            for channel in ("A", "C", "OA", "HD"):
                path = root / f"tile.{channel}.map"
                write_map(path, npts=(22, 22, 24) if channel == "OA" else (20, 22, 24))
                paths[channel] = path
            with self.assertRaisesRegex(ValueError, "nonidentical map geometry"):
                runner.validate_map_set(paths)
            paths["OA"].unlink()
            with self.assertRaisesRegex(ValueError, "missing map"):
                runner.validate_map_set(paths)


if __name__ == "__main__":
    unittest.main()
