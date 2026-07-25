import pathlib
import sys

import numpy as np
import pytest

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from rank_nylc_m1_ensemble import count_reproduced_clusters, gate_stage_b, rank_stage_a


SEEDS = (26711, 26723, 26737)


def replica(candidate_id, seed, nac_frames, occupancy, longest, cluster_count=0, status="PASS"):
    scientific = "PASS_REPLICA_NAC_PRESENT" if nac_frames else "FAIL_REPLICA_NO_NAC"
    if status == "FAIL_TECHNICAL":
        scientific = "NOT_EVALUATED_NUMERICAL_FAIL"
    if status == "FAIL_UNBOUND":
        scientific = "FAIL_UNBOUND"
    return {
        "candidate_id": candidate_id,
        "velocity_seed": seed,
        "technical_status": "FAIL" if status == "FAIL_TECHNICAL" else "PASS",
        "scientific_status": scientific,
        "nac": {
            "nac_frame_count_after_20ps": nac_frames,
            "analysis_frame_count_after_20ps": 100,
            "nac_occupancy_after_20ps": occupancy,
            "longest_event_after_20ps": {"duration_ps": longest},
        },
        "cross_replica_cluster_count": cluster_count,
    }


def candidate(candidate_id, values, *, cluster_count=0, status="PASS"):
    return [
        replica(candidate_id, seed, n, occ, longest, cluster_count, status)
        for seed, (n, occ, longest) in zip(SEEDS, values)
    ]


def test_stage_a_ranks_conformations_not_replicas_and_excludes_all_zero_nac():
    audits = []
    audits += candidate("c_high_repro", [(8, 0.08, 6.0), (7, 0.07, 5.0), (6, 0.06, 4.0)], cluster_count=2)
    audits += candidate("c_two_repro", [(8, 0.08, 6.0), (7, 0.07, 5.0), (0, 0.00, 0.0)], cluster_count=1)
    audits += candidate("c_one_repro", [(9, 0.09, 8.0), (0, 0.00, 0.0), (0, 0.00, 0.0)])
    audits += candidate("c_low_occ", [(2, 0.02, 2.0), (2, 0.02, 2.0), (1, 0.01, 0.0)])
    audits += candidate("c_extra_a", [(3, 0.03, 2.0), (1, 0.01, 0.0), (0, 0.00, 0.0)])
    audits += candidate("c_extra_b", [(1, 0.01, 0.0), (1, 0.01, 0.0), (1, 0.01, 0.0)])
    audits += candidate("c_unbound", [(5, 0.05, 4.0), (4, 0.04, 2.0), (3, 0.03, 2.0)], status="FAIL_UNBOUND")
    audits += candidate("c_technical", [(5, 0.05, 4.0), (4, 0.04, 2.0), (3, 0.03, 2.0)], status="FAIL_TECHNICAL")
    audits += candidate("c_zero", [(0, 0.00, 0.0), (0, 0.00, 0.0), (0, 0.00, 0.0)])

    result = rank_stage_a(audits, top_n=6)

    selected = result["selected_candidate_ids"]
    assert selected[0] == "c_high_repro"
    assert len(selected) == 6
    assert len(set(selected)) == 6
    assert "c_zero" not in selected
    assert "c_unbound" not in selected
    assert "c_technical" not in selected
    assert result["excluded"]["c_zero"] == "FAIL_UNRESTRAINED_NO_NAC_EARLY"
    assert result["ranking"][0]["replica_count_with_nac_after_20ps"] == 3
    assert result["ranking"][0]["pooled_nac_occupancy_after_20ps"] == 21 / 300


def test_fewer_than_six_evaluable_candidates_is_explicit():
    audits = candidate("c_pass", [(2, 0.02, 2.0), (0, 0.0, 0.0), (0, 0.0, 0.0)])
    audits += candidate("c_zero", [(0, 0.0, 0.0), (0, 0.0, 0.0), (0, 0.0, 0.0)])

    result = rank_stage_a(audits, top_n=6)

    assert result["selected_candidate_ids"] == ["c_pass"]
    assert result["selection_status"] == "NOT_EVALUATED_STAGEB_FEWER_THAN_SIX"


def write_fingerprints(path, alignment_frames, local_frames):
    np.savez_compressed(
        path,
        time_ps=np.arange(len(alignment_frames), dtype=float),
        alignment_A=np.asarray(alignment_frames, dtype=float),
        local_A=np.asarray(local_frames, dtype=float),
    )


def test_reproduced_cluster_requires_four_frames_and_two_replicas(tmp_path):
    alignment = np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    local = np.asarray([[0.0, 0.0, 1.0], [1.0, 1.0, 1.0]])
    rotation = np.asarray([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    shifted_alignment = alignment @ rotation + np.asarray([5.0, -2.0, 3.0])
    shifted_local = local @ rotation + np.asarray([5.0, -2.0, 3.0])

    paths = []
    for seed, a, loc in (
        (26711, alignment, local),
        (26723, shifted_alignment, shifted_local),
    ):
        path = tmp_path / f"{seed}.npz"
        write_fingerprints(path, [a, a], [loc, loc])
        paths.append((seed, path))

    assert count_reproduced_clusters(paths, cutoff_nm=0.20) == 1


def test_single_replica_cluster_is_not_reproduced(tmp_path):
    alignment = np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    local = np.asarray([[0.0, 0.0, 1.0], [1.0, 1.0, 1.0]])
    path = tmp_path / "single.npz"
    write_fingerprints(path, [alignment] * 4, [local] * 4)

    assert count_reproduced_clusters([(26711, path)], cutoff_nm=0.20) == 0


def stage_b_replica(
    seed,
    *,
    nac_frames=5,
    denominator=100,
    longest=4.0,
    bound=True,
    technical=True,
    reproduced_clusters=1,
):
    return {
        "candidate_id": "stage_b_candidate",
        "velocity_seed": seed,
        "technical_status": "PASS" if technical else "FAIL",
        "scientific_status": "PASS_REPLICA_NAC_PRESENT",
        "analysis": {"window_ps": [100.0, 1000.0]},
        "nac": {
            "nac_frame_count": nac_frames,
            "frame_count": denominator,
            "nac_occupancy": nac_frames / denominator,
            "longest_event": {"duration_ps": longest},
        },
        "bound_state": {
            "final_frame_retained": bound,
            "severe_clash": False,
        },
        "thermodynamics": {"stable": True},
        "cross_replica_cluster_count": reproduced_clusters,
    }


def passing_stage_b_audits():
    return [stage_b_replica(seed) for seed in SEEDS]


def test_stage_b_exact_pass_gate():
    decision = gate_stage_b(passing_stage_b_audits())

    assert decision["replicas_with_nac_after_100ps"] >= 2
    assert decision["pooled_nac_occupancy_100_1000ps"] >= 0.01
    assert decision["longest_continuous_nac_ps"] >= 4.0
    assert decision["bound_replica_count"] == 3
    assert decision["technical_failure_count"] == 0
    assert decision["cross_replica_cluster_reproduced"] is True
    assert decision["scientific_status"] == "PASS_UNRESTRAINED_M1_ENSEMBLE_NAC"
    assert decision["failed_gates"] == []


@pytest.mark.parametrize(
    "audits,reason,expected_status",
    [
        (
            [
                stage_b_replica(SEEDS[0], nac_frames=30),
                stage_b_replica(SEEDS[1], nac_frames=0),
                stage_b_replica(SEEDS[2], nac_frames=0),
            ],
            "INSUFFICIENT_NAC_REPLICA_REPRODUCTION",
            "FAIL_UNRESTRAINED_M1_ENSEMBLE_NAC",
        ),
        (
            [stage_b_replica(seed, nac_frames=1, denominator=200) for seed in SEEDS],
            "POOLED_NAC_OCCUPANCY_BELOW_0.01",
            "FAIL_UNRESTRAINED_M1_ENSEMBLE_NAC",
        ),
        (
            [stage_b_replica(seed, longest=2.0) for seed in SEEDS],
            "LONGEST_NAC_EVENT_BELOW_4PS",
            "FAIL_UNRESTRAINED_M1_ENSEMBLE_NAC",
        ),
        (
            [
                stage_b_replica(SEEDS[0], bound=False),
                stage_b_replica(SEEDS[1]),
                stage_b_replica(SEEDS[2]),
            ],
            "UNBOUND_REPLICA",
            "FAIL_UNRESTRAINED_M1_ENSEMBLE_NAC",
        ),
        (
            [
                stage_b_replica(SEEDS[0], technical=False),
                stage_b_replica(SEEDS[1]),
                stage_b_replica(SEEDS[2]),
            ],
            "TECHNICAL_FAILURE",
            "NOT_EVALUATED_TECHNICAL_FAILURE",
        ),
        (
            [stage_b_replica(seed, reproduced_clusters=0) for seed in SEEDS],
            "NO_REPRODUCED_CROSS_REPLICA_CLUSTER",
            "FAIL_UNRESTRAINED_M1_ENSEMBLE_NAC",
        ),
    ],
)
def test_stage_b_preserves_each_failed_gate_reason(audits, reason, expected_status):
    decision = gate_stage_b(audits)

    assert decision["scientific_status"] == expected_status
    assert reason in decision["failed_gates"]
