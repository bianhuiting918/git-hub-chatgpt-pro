import re
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
INSTALL = PROJECT / "scripts" / "install_phase1.sh"
LOCK = PROJECT / "config" / "analysis-requirements.lock"


class InstallPhase1ContractTests(unittest.TestCase):
    def test_installer_is_server_local_and_guarded(self):
        text = INSTALL.read_text(encoding="utf-8")
        self.assertIn("/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725", text)
        self.assertIn("/work/home/acshdt1dks/", text)
        self.assertIn("MAX_JOBS", text)
        self.assertRegex(text, r"MAX_JOBS[^\n]*64")
        self.assertNotRegex(text, r"(^|\s)sudo(\s|$)")
        self.assertIn("envs/surface-screen-py311", text)
        self.assertIn("cache/pip", text)

    def test_analysis_python_is_created_with_ssl_capable_conda(self):
        text = INSTALL.read_text(encoding="utf-8")
        self.assertIn("/work/home/acshdt1dks/anaconda3/bin/conda", text)
        self.assertIn('create --prefix "$ENV_DIR" --yes', text)
        self.assertRegex(text, r"python=3\.11")
        ssl_check = '"$ENV_DIR/bin/python" -c "import ssl"'
        pip_install = '"$ENV_DIR/bin/python" -m pip install'
        self.assertIn(ssl_check, text)
        self.assertIn(pip_install, text)
        self.assertLess(text.index(ssl_check), text.index(pip_install))
        self.assertNotIn("-m venv", text)
        self.assertNotIn("rm -rf", text)

    def test_bootstrap_pip_is_not_pinned_to_unavailable_defaults_version(self):
        text = INSTALL.read_text(encoding="utf-8")
        self.assertIn('"python=3.11" "pip"', text)
        self.assertNotIn('"pip=24.3.1"', text)

    def test_conda_uses_copy_mode_for_hpc_filesystem(self):
        text = INSTALL.read_text(encoding="utf-8")
        self.assertIn('CONDA_ALWAYS_COPY=true CONDA_PKGS_DIRS="$CONDA_CACHE"', text)

    def test_rdkit_pin_is_available_on_sugon_pypi_mirror(self):
        lock_text = LOCK.read_text(encoding="utf-8")
        self.assertIn("rdkit==2024.3.2", lock_text)
        self.assertNotIn("rdkit==2024.9.6", lock_text)

    def test_each_run_uses_a_unique_log_directory(self):
        text = INSTALL.read_text(encoding="utf-8")
        stamp = "run_stamp=$(date -u +%Y%m%dT%H%M%SZ)"
        log_dir = 'LOG_DIR="$ROOT/logs/install_phase1/$run_stamp"'
        self.assertIn(stamp, text)
        self.assertIn(log_dir, text)
        self.assertLess(text.index(stamp), text.index(log_dir))

    def test_slurm_spool_copy_uses_explicit_repo_project_path(self):
        text = INSTALL.read_text(encoding="utf-8")
        self.assertIn(
            "/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725/"
            "repo/projects/03-polymer-surface-hotspot-screening",
            text,
        )
        self.assertNotIn('dirname "${BASH_SOURCE[0]}"', text)

    def test_gate_is_written_only_after_version_and_import_checks(self):
        text = INSTALL.read_text(encoding="utf-8")
        checks = [
            '"$ADFR_BIN/autosite"',
            '"$ADFR_BIN/autogrid4"',
            '"$ADFR_BIN/prepare_receptor"',
            "import numpy",
            "import rdkit",
            "import freesasa",
        ]
        gate = "INSTALL_PASS.json"
        self.assertIn(gate, text)
        gate_position = text.rindex(gate)
        for token in checks:
            self.assertIn(token, text)
            self.assertLess(text.index(token), gate_position)

    def test_analysis_dependencies_are_exactly_pinned(self):
        lines = [
            line.strip()
            for line in LOCK.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
        required = {
            "numpy",
            "scipy",
            "pandas",
            "biopython",
            "pyyaml",
            "rdkit",
            "freesasa",
            "pytest",
        }
        observed = set()
        for line in lines:
            self.assertRegex(line, r"^[A-Za-z0-9_.-]+==[^=\s]+$")
            observed.add(re.split(r"==", line, maxsplit=1)[0].lower())
        self.assertEqual(observed, required)


if __name__ == "__main__":
    unittest.main()
