import json
import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "run_surface_field_parallel_chunk.py"


def write_manifest(path: Path, rows: int) -> None:
    header = (
        "candidate_id\tmaterial_family\tsequence_md5\tinput_status\t"
        "selected_chain\treceptor_path\treceptor_sha256\n"
    )
    body = []
    for index in range(rows):
        body.append(
            f"candidate_{index}\tPET\t{index:032x}\tREADY_FOR_STRUCTURE_EVALUATION\t"
            "A\t/unused/receptor.pdb\t" + "0" * 64 + "\n"
        )
    path.write_text(header + "".join(body), encoding="utf-8")


def write_fake_worker(path: Path) -> None:
    path.write_text(
        """#!/usr/bin/env python3
import argparse
import json
import time
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("--row-index", type=int, required=True)
p.add_argument("--manifest")
p.add_argument("--family")
p.add_argument("--output-root")
p.add_argument("--project-root")
p.add_argument("--tool-root")
p.add_argument("--adfr-root")
a = p.parse_args()
root = Path(a.output_root)
root.mkdir(parents=True, exist_ok=True)
(root / f"start_{a.row_index}.json").write_text(json.dumps({"t": time.time()}))
time.sleep(0.35)
(root / f"end_{a.row_index}.json").write_text(json.dumps({"t": time.time()}))
raise SystemExit(7 if a.row_index == 2 else 0)
""",
        encoding="utf-8",
    )


def test_parallel_chunk_skips_primary_pass_and_isolates_row_failure(tmp_path):
    manifest = tmp_path / "manifest.tsv"
    write_manifest(manifest, 5)
    worker = tmp_path / "fake_worker.py"
    write_fake_worker(worker)

    primary = tmp_path / "primary"
    passed = primary / f"{0:032x}"
    passed.mkdir(parents=True)
    (passed / "FIELD_PASS.json").write_text(
        json.dumps({"status": "FIELD_PASS"}), encoding="utf-8"
    )
    recovery = tmp_path / "recovery"

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--manifest",
            str(manifest),
            "--family",
            "PET",
            "--output-root",
            str(recovery),
            "--project-root",
            str(tmp_path),
            "--tool-root",
            str(tmp_path),
            "--adfr-root",
            str(tmp_path),
            "--row-worker",
            str(worker),
            "--start-row",
            "0",
            "--end-row",
            "5",
            "--workers",
            "4",
            "--existing-pass-root",
            str(primary),
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 10
    assert "skipped_primary_pass=1" in result.stdout
    assert "failed=1" in result.stdout
    assert not (recovery / "start_0.json").exists()
    assert (recovery / "start_1.json").exists()
    assert (recovery / "start_2.json").exists()
    assert (recovery / "start_3.json").exists()
    assert (recovery / "start_4.json").exists()

    starts = [
        json.loads((recovery / f"start_{i}.json").read_text())["t"]
        for i in (1, 2, 3, 4)
    ]
    ends = [
        json.loads((recovery / f"end_{i}.json").read_text())["t"]
        for i in (1, 2, 3, 4)
    ]
    assert max(starts) < min(ends)


def test_slurm_contract_requests_four_cpus():
    text = (
        Path(__file__).parents[1]
        / "slurm"
        / "surface_field_parallel_chunk_array.sbatch"
    ).read_text(encoding="utf-8")
    assert "#SBATCH --cpus-per-task=4" in text
    assert "--workers \"$SLURM_CPUS_PER_TASK\"" in text
    assert 'EXISTING_PASS_ROOT=${EXISTING_PASS_ROOT:?EXISTING_PASS_ROOT is required}' in text
