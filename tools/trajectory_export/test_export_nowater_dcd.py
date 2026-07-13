import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).parent))


class SampleFramesTest(unittest.TestCase):
    def test_returns_requested_number_of_evenly_spaced_indices(self):
        from export_nowater_dcd import sample_indices

        indices = sample_indices(73, 50)

        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[-1], 72)
        self.assertEqual(len(indices), 50)
        self.assertEqual(len(set(indices)), 50)


if __name__ == "__main__":
    unittest.main()

