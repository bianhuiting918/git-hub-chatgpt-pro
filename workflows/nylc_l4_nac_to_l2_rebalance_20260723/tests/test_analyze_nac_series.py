import json
import subprocess
import sys
from pathlib import Path

FLOW = Path(__file__).resolve().parents[1]
SCRIPT = FLOW / "scripts" / "analyze_nac_series.py"


def _write_xvg(path: Path, values):
    path.write_text(
        "\n".join(["# synthetic", '@ title "test"'] + [f"{time} {value}" for time, value in values]) + "\n"
    )


def test_reports_occupancy_and_longest_continuous_residence(tmp_path):
    distance = tmp_path / "distance.xvg"
    angle = tmp_path / "angle.xvg"
    output = tmp_path / "audit.json"
    times = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0]
    _write_xvg(distance, list(zip(times, [0.34, 0.33, 0.36, 0.32, 0.31, 0.34, 0.37])))
    _write_xvg(angle, list(zip(times, [100, 110, 105, 94, 100, 115, 108])))

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--distance-xvg",
            str(distance),
            "--angle-xvg",
            str(angle),
            "--output",
            str(output),
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0, result.stderr
    audit = json.loads(output.read_text())
    assert audit["frame_count"] == 7
    assert audit["nac_frame_count"] == 4
    assert audit["nac_occupancy"] == 4 / 7
    assert audit["longest_continuous_nac"]["start_ps"] == 8.0
    assert audit["longest_continuous_nac"]["end_ps"] == 10.0
    assert audit["longest_continuous_nac"]["duration_ps"] == 2.0
    assert audit["longest_continuous_nac"]["frame_count"] == 2


def test_rejects_misaligned_time_series(tmp_path):
    distance = tmp_path / "distance.xvg"
    angle = tmp_path / "angle.xvg"
    output = tmp_path / "audit.json"
    _write_xvg(distance, [(0.0, 0.3), (2.0, 0.3)])
    _write_xvg(angle, [(0.0, 100), (3.0, 100)])

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--distance-xvg",
            str(distance),
            "--angle-xvg",
            str(angle),
            "--output",
            str(output),
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0
    assert "time series do not align" in result.stderr
