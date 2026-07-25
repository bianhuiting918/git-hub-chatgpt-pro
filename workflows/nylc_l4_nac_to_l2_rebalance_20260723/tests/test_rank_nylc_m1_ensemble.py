import pathlib
import sys

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from rank_nylc_m1_ensemble import rank_stage_a


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
    audits += candidate("c_unbound", [(5, 0.05, 4.0), (4, 0.04, 2.0), (3, 0.03, 2.0)], status="FAIL_UNBOUND")
    audits += candidate("c_technical", [(5, 0.05, 4.0), (4, 0.04, 2.0), (3, 0.03, 2.0)], status="FAIL_TECHNICAL")
    audits += candidate("c_zero", [(0, 0.00, 0.0), (0, 0.00, 0.0), (0, 0.00, 0.0)])

    result = rank_stage_a(audits, top_n=6)

    selected = result["selected_candidate_ids"]
    assert selected[0] == "c_high_repro"
    assert len(selected) == 6
    assert len(set(selected)) == 6
    assert "c_zero" not in selected
    assert result["excluded"]["c_zero"] == "FAIL_UNRESTRAINED_NO_NAC_EARLY"
    assert result["ranking"][0]["replica_count_with_nac_after_20ps"] == 3
    assert result["ranking"][0]["pooled_nac_occupancy_after_20ps"] == 21 / 300


def test_fewer_than_six_evaluable_candidates_is_explicit():
    audits = candidate("c_pass", [(2, 0.02, 2.0), (0, 0.0, 0.0), (0, 0.0, 0.0)])
    audits += candidate("c_zero", [(0, 0.0, 0.0), (0, 0.0, 0.0), (0, 0.0, 0.0)])

    result = rank_stage_a(audits, top_n=6)

    assert result["selected_candidate_ids"] == ["c_pass"]
    assert result["selection_status"] == "NOT_EVALUATED_STAGEB_FEWER_THAN_SIX"
