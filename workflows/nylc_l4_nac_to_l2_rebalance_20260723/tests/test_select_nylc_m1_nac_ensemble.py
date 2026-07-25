import sys
from pathlib import Path

import numpy as np

FLOW = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(FLOW / "scripts"))

from select_nylc_m1_nac_ensemble import (
    group_nac_events,
    kabsch_rmsd,
    read_gro,
    select_twelve,
    _write_frame_index,
)


def _synthetic_85_nac_rows():
    frame_counts = [3, 3, 3, 3, 2, 2, 2, 2, 2] + [1] * 63
    rows = []
    time_ps = 1000.0
    for event_index, frame_count in enumerate(frame_counts):
        for member_index in range(frame_count):
            rows.append(
                {
                    "time_ps": time_ps + 2.0 * member_index,
                    "distance_nm": 0.30 + 0.0001 * event_index,
                    "angle_deg": 105.0,
                    "potential_energy_kj_mol": -1000.0 - event_index,
                    "nac": True,
                }
            )
        time_ps += 2.0 * frame_count + 2.0
    return rows


def test_groups_85_frames_into_72_time_events():
    events = group_nac_events(_synthetic_85_nac_rows(), sample_interval_ps=2.0)

    assert len(events) == 72
    assert sum(event["frame_count"] for event in events) == 85
    assert [event["frame_count"] for event in events[:9]] == [3, 3, 3, 3, 2, 2, 2, 2, 2]


def _selection_fixture():
    rows = _synthetic_85_nac_rows()
    events = group_nac_events(rows, sample_interval_ps=2.0)
    for event in events:
        for member in event["members"]:
            member["distance_nm"] = 0.34
            member["angle_deg"] = 112.0
            member["potential_energy_kj_mol"] = -100.0 - member["frame_index"]

    # Persistence events are deliberately also central and low energy. They may
    # occupy only their persistence category, never central_low_energy.
    for event in events[:9]:
        for member in event["members"]:
            member["distance_nm"] = 0.30
            member["angle_deg"] = 105.0
            member["potential_energy_kj_mol"] = -10000.0 - member["frame_index"]

    # Make six single-frame independent central/energy candidates.
    for event in events[9:15]:
        member = event["members"][0]
        member["distance_nm"] = 0.30
        member["angle_deg"] = 105.0
        member["potential_energy_kj_mol"] = -5000.0 - member["frame_index"]

    # Preserve the historical 1462 ps source as an explicit control event.
    legacy = events[-1]
    legacy["members"][0]["time_ps"] = 1462.0
    legacy["start_ps"] = 1462.0
    legacy["end_ps"] = 1462.0

    frame_count = sum(event["frame_count"] for event in events)
    rmsd_nm = np.full((frame_count, frame_count), 0.25)
    np.fill_diagonal(rmsd_nm, 0.0)
    return events, rmsd_nm


def test_selects_deterministic_4_plus_4_plus_3_plus_1_without_event_reuse():
    events, rmsd_nm = _selection_fixture()

    selected = select_twelve(events, rmsd_nm)

    assert len(selected) == 12
    assert len({item["event_id"] for item in selected}) == 12
    assert sum(item["category"] == "three_frame_medoid" for item in selected) == 4
    assert sum(item["category"] == "two_frame_medoid" for item in selected) == 4
    assert sum(item["category"] == "central_low_energy" for item in selected) == 3
    assert [item for item in selected if item["time_ps"] == 1462.0][0]["category"] == "legacy_control"
    assert select_twelve(events, rmsd_nm) == selected


def test_persistence_event_cannot_also_fill_central_energy_category():
    events, rmsd_nm = _selection_fixture()

    selected = select_twelve(events, rmsd_nm)
    persistence_ids = {
        item["event_id"]
        for item in selected
        if item["category"] in {"three_frame_medoid", "two_frame_medoid"}
    }
    central_ids = {
        item["event_id"] for item in selected if item["category"] == "central_low_energy"
    }

    assert persistence_ids.isdisjoint(central_ids)


def test_reads_gro_and_kabsch_is_translation_rotation_invariant(tmp_path):
    gro = tmp_path / "tiny.gro"
    gro.write_text(
        "tiny\n"
        "3\n"
        "    1ALA      N    1   0.000   0.000   0.000\n"
        "    1ALA     CA    2   0.100   0.000   0.000\n"
        "    1ALA      C    3   0.000   0.100   0.000\n"
        "   1.00000   1.00000   1.00000\n"
    )
    atoms, box = read_gro(gro)

    reference = np.array([atom["coord_nm"] for atom in atoms])
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    mobile = reference @ rotation + np.array([2.0, -3.0, 0.5])

    assert len(atoms) == 3
    assert box == (1.0, 1.0, 1.0)
    assert kabsch_rmsd(reference, mobile) < 1e-12


def test_selection_batch_does_not_launch_mpi_gromacs_inside_an_srun_step():
    batch = FLOW / "slurm" / "run_nylc_m1_ensemble_select.sbatch"
    text = batch.read_text()

    assert 'srun -n 1 -c "${SLURM_CPUS_PER_TASK:-8}" "$PYTHON"' not in text
    assert '"$PYTHON" \\\n    "$FLOW/scripts/select_nylc_m1_nac_ensemble.py"' in text


def test_gromacs_frame_index_file_converts_zero_based_source_indices(tmp_path):
    frame_index = tmp_path / "frames.ndx"

    _write_frame_index(frame_index, [0, 14, 5830])

    assert frame_index.read_text().splitlines() == ["[ frames ]", "1 15 5831"]
