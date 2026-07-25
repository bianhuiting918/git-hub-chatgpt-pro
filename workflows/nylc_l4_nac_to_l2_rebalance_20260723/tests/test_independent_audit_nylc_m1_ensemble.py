import copy
import pathlib
import sys

import pytest

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from independent_audit_nylc_m1_ensemble import audit_final


SEEDS = (26711, 26723, 26737)


def stage_b_replica(candidate_id, seed):
    frames = [
        {"time_ps": time, "nac": True, "bound": True}
        for time in (100.0, 102.0, 104.0, 106.0)
    ]
    return {
        "candidate_id": candidate_id,
        "velocity_seed": seed,
        "path": f"/task/{candidate_id}/seed{seed}/stageB_attempt_new",
        "complete_json": {"status": "PASS_TECHNICAL_STAGE_B"},
        "free_tpr_contract": {
            "position_restraints": 0,
            "distance_restraints": 0,
            "processed_mdp_nonempty_define_values": [],
        },
        "numerical": {
            "finished_mdrun": True,
            "fatal": 0,
            "lincs_warning": 0,
            "settle_problem": 0,
            "nan": 0,
        },
        "thermodynamics": {"stable": True},
        "primitive_frames": frames,
        "advertised_scientific_status": "PASS_UNRESTRAINED_M1_ENSEMBLE_NAC",
    }


def valid_payload():
    candidates = [f"c{index:02d}" for index in range(12)]
    stage_a = []
    slot = 0
    for candidate in candidates:
        for seed in SEEDS:
            stage_a.append(
                {
                    "stageA_slot": slot,
                    "candidate_id": candidate,
                    "velocity_seed": seed,
                    "technical_status": "PASS",
                    "scientific_status": "PASS_REPLICA_NAC_PRESENT",
                }
            )
            slot += 1
    selected = [candidates[0]]
    stage_b = [stage_b_replica(selected[0], seed) for seed in SEEDS]
    return {
        "schema_version": 1,
        "authority": {
            "selected_candidate_ids": candidates,
            "gate_definition": "NylC residues 261-266; Thr267 excluded",
        },
        "stageA": {
            "expected_replica_count": 36,
            "replica_records": stage_a,
            "selected_candidate_ids": selected,
        },
        "stageB": {
            "replica_records": stage_b,
            "cluster_summaries": {
                selected[0]: {
                    "reproduced": True,
                    "cluster_count": 1,
                    "medoid": {
                        "source_replica_seed": 26711,
                        "time_ps": 102.0,
                        "source_coordinate_sha256": "a" * 64,
                        "reaction_geometry": {"distance_nm": 0.32, "angle_deg": 105.0},
                        "gate_opening_nm": 1.1,
                        "proton_route_annotation": "OgH_to_Nalpha_preorg",
                        "membership_by_replica": {"26711": 2, "26723": 2},
                    },
                }
            },
        },
        "forbidden_job_ids": ["61801874", "61803121"],
    }


def test_independent_audit_passes_complete_recurrent_candidate():
    result = audit_final(valid_payload())

    assert result["technical_status"] == "PASS"
    assert result["candidate_universe"] == 12
    assert result["stageA_expected_replicas"] == 36
    assert result["stageB_expected_replicas"] == 3
    assert result["scientific_pass_count"] == 1
    assert result["qmmm_eligible_count"] == 1
    assert result["qmmm_eligibility"][0]["medoid_is_qmmm_optimized_gs"] is False


@pytest.mark.parametrize(
    "mutator,match",
    [
        (
            lambda payload: payload["authority"]["selected_candidate_ids"].__setitem__(1, "c00"),
            "duplicate candidate",
        ),
        (
            lambda payload: payload["stageA"]["replica_records"][1].__setitem__(
                "velocity_seed", 26711
            ),
            "duplicate seed",
        ),
        (
            lambda payload: payload["stageA"]["replica_records"].pop(),
            "36 Stage A",
        ),
        (
            lambda payload: payload["stageB"]["replica_records"][0].__setitem__(
                "candidate_id", "c01"
            ),
            "absent from Stage A top six",
        ),
        (
            lambda payload: payload["authority"].__setitem__(
                "gate_definition", "NylC residues 261-267 including Thr267"
            ),
            "Thr267",
        ),
        (
            lambda payload: payload["stageB"]["replica_records"][0][
                "free_tpr_contract"
            ].__setitem__("position_restraints", 1),
            "restraint",
        ),
        (
            lambda payload: payload["stageB"]["replica_records"][0].__setitem__(
                "path", "/task/legacy/61801874/output"
            ),
            "forbidden",
        ),
        (
            lambda payload: payload["stageB"]["replica_records"][0].pop(
                "primitive_frames"
            ),
            "primitive",
        ),
    ],
)
def test_independent_audit_rejects_corrupt_or_superseded_inputs(mutator, match):
    payload = valid_payload()
    mutator(payload)

    with pytest.raises(ValueError, match=match):
        audit_final(payload)


def test_advertised_pass_is_rejected_when_primitive_occupancy_is_below_one_percent():
    payload = valid_payload()
    frames = [
        {"time_ps": 100.0 + 2.0 * index, "nac": index == 0, "bound": True}
        for index in range(200)
    ]
    for record in payload["stageB"]["replica_records"]:
        record["primitive_frames"] = copy.deepcopy(frames)

    with pytest.raises(ValueError, match="advertised PASS"):
        audit_final(payload)
