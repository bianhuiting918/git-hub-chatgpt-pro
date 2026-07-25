import pathlib
import sys

import numpy as np
import pytest

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from assemble_nylc_m1_final_audit_input import (
    cluster_candidate,
    load_raw_frames,
)


def write_xvg(path, rows):
    path.write_text(
        "\n".join(" ".join(str(value) for value in row) for row in rows) + "\n"
    )


def test_load_raw_frames_joins_only_identical_time_axes(tmp_path):
    write_xvg(tmp_path / "nac_distance.xvg", [(100, 0.32), (102, 0.50)])
    write_xvg(tmp_path / "nac_angle.xvg", [(100, 105), (102, 90)])
    write_xvg(tmp_path / "gate_opening.xvg", [(100, 1.1), (102, 1.2)])
    write_xvg(tmp_path / "pocket_state.xvg", [(100, 6, 0.8), (102, 2, 1.3)])

    frames = load_raw_frames(tmp_path)

    assert frames == [
        {
            "time_ps": 100.0,
            "distance_nm": 0.32,
            "angle_deg": 105.0,
            "gate_opening_nm": 1.1,
            "pocket_contact_count": 6,
            "ligand_pocket_com_nm": 0.8,
        },
        {
            "time_ps": 102.0,
            "distance_nm": 0.5,
            "angle_deg": 90.0,
            "gate_opening_nm": 1.2,
            "pocket_contact_count": 2,
            "ligand_pocket_com_nm": 1.3,
        },
    ]

    write_xvg(tmp_path / "nac_angle.xvg", [(100, 105), (103, 90)])
    with pytest.raises(ValueError, match="time axes"):
        load_raw_frames(tmp_path)


def fingerprint(path, times, local_offset=0.0):
    alignment = np.asarray(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
            for _ in times
        ],
        dtype=float,
    )
    local = np.asarray(
        [
            [[0.0, 0.0, 0.0], [1.0 + local_offset, 0.0, 0.0]]
            for _ in times
        ],
        dtype=float,
    )
    np.savez_compressed(
        path,
        time_ps=np.asarray(times, dtype=float),
        alignment_A=alignment,
        local_A=local,
    )


def primitive_frames(times):
    return {
        round(float(time), 6): {
            "time_ps": float(time),
            "distance_nm": 0.32,
            "angle_deg": 105.0,
            "gate_opening_nm": 1.1,
            "pocket_contact_count": 6,
            "ligand_pocket_com_nm": 0.8,
            "proton_route_annotation": "OgH_to_Nalpha_preorg",
        }
        for time in times
    }


def test_cluster_candidate_reports_reproduced_medoid_without_importing_ranker(tmp_path):
    paths = {}
    frame_maps = {}
    for seed, offset in ((26711, 0.00), (26723, 0.02), (26737, 8.0)):
        path = tmp_path / f"{seed}.npz"
        fingerprint(path, [100.0, 102.0], local_offset=offset)
        paths[seed] = path
        frame_maps[seed] = primitive_frames([100.0, 102.0])

    def coordinate_resolver(seed, time_ps):
        return {
            "source_coordinate_sha256": "a" * 64,
            "source_coordinate_path": f"/medoids/{seed}_{time_ps:.0f}.gro",
        }

    summary = cluster_candidate(
        "candidate_a",
        paths,
        frame_maps,
        coordinate_resolver=coordinate_resolver,
        cutoff_nm=0.20,
    )

    assert summary["reproduced"] is True
    assert summary["cluster_count"] == 1
    assert summary["medoid"]["membership_by_replica"] == {
        "26711": 2,
        "26723": 2,
    }
    assert summary["medoid"]["reaction_geometry"] == {
        "distance_nm": 0.32,
        "angle_deg": 105.0,
    }
    assert summary["medoid"]["source_coordinate_sha256"] == "a" * 64


def test_cluster_candidate_rejects_single_replica_recurrence(tmp_path):
    path = tmp_path / "26711.npz"
    fingerprint(path, [100.0, 102.0, 104.0, 106.0])

    summary = cluster_candidate(
        "candidate_a",
        {26711: path},
        {26711: primitive_frames([100.0, 102.0, 104.0, 106.0])},
        coordinate_resolver=lambda seed, time: {},
    )

    assert summary == {
        "reproduced": False,
        "cluster_count": 0,
        "medoid": None,
    }
