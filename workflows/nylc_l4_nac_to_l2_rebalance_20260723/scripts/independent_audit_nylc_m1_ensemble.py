#!/usr/bin/env python3
"""Independent final audit of the NylC M1 real-NAC ensemble."""

from __future__ import annotations

import argparse
import csv
import json
import math
import pathlib
import re
import statistics
from collections import defaultdict

EXPECTED_SEEDS = (26711, 26723, 26737)
PASS_STATUS = "PASS_UNRESTRAINED_M1_ENSEMBLE_NAC"


def _longest_event(frames: list[dict]) -> float:
    if not frames:
        return 0.0
    times = [float(frame["time_ps"]) for frame in frames]
    if len(times) > 1:
        spacing = statistics.median(
            right - left for left, right in zip(times, times[1:])
        )
    else:
        spacing = 2.0
    best = 0.0
    start = None
    previous = None
    for frame in frames:
        time = float(frame["time_ps"])
        active = bool(frame["nac"])
        if active and (start is None or (previous is not None and time - previous > spacing + 1.0e-6)):
            if start is not None and previous is not None:
                best = max(best, previous - start)
            start = time
        elif not active and start is not None:
            if previous is not None:
                best = max(best, previous - start)
            start = None
        previous = time
    if start is not None and previous is not None:
        best = max(best, previous - start)
    return best


def _validate_authority(payload: dict) -> tuple[list[str], list[str]]:
    authority = payload["authority"]
    candidates = [str(value) for value in authority["selected_candidate_ids"]]
    if len(candidates) != 12:
        raise ValueError("candidate universe must contain exactly 12 candidates")
    if len(candidates) != len(set(candidates)):
        raise ValueError("duplicate candidate IDs in authority")
    gate = str(authority["gate_definition"])
    if "Thr267 excluded" not in gate or "261-266" not in gate:
        raise ValueError("gate definition must use residues 261-266 with Thr267 excluded")
    forbidden = [str(value) for value in payload.get("forbidden_job_ids", [])]
    return candidates, forbidden


def _validate_stage_a(payload: dict, candidates: list[str]) -> list[str]:
    stage_a = payload["stageA"]
    expected = int(stage_a.get("expected_replica_count", 36))
    records = stage_a["replica_records"]
    if expected != 36 or len(records) != 36:
        raise ValueError("expected exactly 36 Stage A replica records")
    slots = [int(record["stageA_slot"]) for record in records]
    if sorted(slots) != list(range(36)):
        raise ValueError("36 Stage A slots are incomplete or duplicated")
    pairs = [(str(record["candidate_id"]), int(record["velocity_seed"])) for record in records]
    if len(pairs) != len(set(pairs)):
        raise ValueError("duplicate seed within a Stage A candidate")
    grouped = defaultdict(set)
    for candidate, seed in pairs:
        if candidate not in candidates:
            raise ValueError("Stage A candidate is absent from authority")
        grouped[candidate].add(seed)
    for candidate in candidates:
        if grouped[candidate] != set(EXPECTED_SEEDS):
            raise ValueError(f"duplicate seed or incomplete seed triplet for {candidate}")
    selected = [str(value) for value in stage_a["selected_candidate_ids"]]
    if len(selected) > 6 or len(selected) != len(set(selected)):
        raise ValueError("Stage A top six contains duplicate or excessive candidates")
    if not set(selected).issubset(candidates):
        raise ValueError("Stage A top six contains a candidate absent from authority")
    return selected


def _validate_path(path: str, forbidden: list[str]) -> None:
    for job_id in forbidden:
        if job_id and job_id in path:
            raise ValueError(f"forbidden legacy or cancelled output path contains {job_id}")


def _replica_metrics(record: dict, forbidden: list[str]) -> dict:
    _validate_path(str(record["path"]), forbidden)
    if "primitive_frames" not in record:
        raise ValueError("COMPLETE.json without primitive frame evidence is not auditable")
    contract = record["free_tpr_contract"]
    if (
        int(contract.get("position_restraints", -1)) != 0
        or int(contract.get("distance_restraints", -1)) != 0
        or contract.get("processed_mdp_nonempty_define_values")
    ):
        raise ValueError("fully unrestrained TPR contract declares a restraint")
    numerical = record["numerical"]
    technical = bool(numerical.get("finished_mdrun")) and not any(
        int(numerical.get(key, 0))
        for key in ("fatal", "lincs_warning", "settle_problem", "nan")
    )
    frames = [
        frame
        for frame in record["primitive_frames"]
        if 100.0 - 1.0e-6 <= float(frame["time_ps"]) <= 1000.0 + 1.0e-6
    ]
    if not frames:
        raise ValueError("primitive Stage B window has no frames")
    times = [float(frame["time_ps"]) for frame in frames]
    if any(not math.isfinite(value) for value in times) or any(
        right <= left for left, right in zip(times, times[1:])
    ):
        raise ValueError("primitive frame times are non-finite, duplicated or unordered")
    nac_count = sum(bool(frame["nac"]) for frame in frames)
    bound = all(bool(frame["bound"]) for frame in frames)
    return {
        "candidate_id": str(record["candidate_id"]),
        "velocity_seed": int(record["velocity_seed"]),
        "technical_pass": technical,
        "thermodynamic_pass": bool(record["thermodynamics"]["stable"]),
        "frame_count": len(frames),
        "nac_frame_count": nac_count,
        "nac_occupancy": nac_count / len(frames),
        "longest_nac_event_ps": _longest_event(frames),
        "bound": bound,
        "advertised_scientific_status": record.get("advertised_scientific_status"),
        "path": str(record["path"]),
    }


def _validate_medoid(candidate_id: str, cluster: dict) -> dict:
    if not cluster.get("reproduced"):
        return {}
    medoid = dict(cluster.get("medoid") or {})
    required = {
        "source_replica_seed",
        "time_ps",
        "source_coordinate_sha256",
        "reaction_geometry",
        "gate_opening_nm",
        "proton_route_annotation",
        "membership_by_replica",
    }
    missing = sorted(required - set(medoid))
    if missing:
        raise ValueError(f"{candidate_id}: recurrent cluster medoid lacks {missing}")
    if not re.fullmatch(r"[0-9a-f]{64}", str(medoid["source_coordinate_sha256"])):
        raise ValueError(f"{candidate_id}: invalid medoid coordinate checksum")
    membership = {str(key): int(value) for key, value in medoid["membership_by_replica"].items()}
    if sum(count >= 2 for count in membership.values()) < 2:
        raise ValueError(f"{candidate_id}: cluster is not reproduced across two replicas")
    medoid["membership_by_replica"] = membership
    medoid["medoid_is_qmmm_optimized_gs"] = False
    medoid["caveat"] = (
        "Cross-replica classical-MM NAC medoid; not yet a QM/MM-optimized ground state."
    )
    return medoid


def audit_final(payload: dict) -> dict:
    candidates, forbidden = _validate_authority(payload)
    selected = _validate_stage_a(payload, candidates)
    records = payload["stageB"]["replica_records"]
    expected_stage_b = 3 * len(selected)
    if len(records) != expected_stage_b:
        raise ValueError(
            f"Stage B expected {expected_stage_b} replica records, got {len(records)}"
        )

    grouped = defaultdict(list)
    seen = set()
    for record in records:
        candidate = str(record["candidate_id"])
        seed = int(record["velocity_seed"])
        if candidate not in selected:
            raise ValueError("Stage B candidate is absent from Stage A top six")
        pair = (candidate, seed)
        if pair in seen:
            raise ValueError("duplicate seed in Stage B candidate")
        seen.add(pair)
        grouped[candidate].append(_replica_metrics(record, forbidden))
    for candidate in selected:
        if {item["velocity_seed"] for item in grouped[candidate]} != set(EXPECTED_SEEDS):
            raise ValueError(f"{candidate}: incomplete Stage B seed triplet")

    clusters = payload["stageB"].get("cluster_summaries", {})
    decisions = []
    failures = []
    eligible = []
    for candidate in selected:
        replicas = grouped[candidate]
        technical_failures = sum(not item["technical_pass"] for item in replicas)
        thermo_failures = sum(not item["thermodynamic_pass"] for item in replicas)
        valid = [item for item in replicas if item["technical_pass"]]
        replicas_with_nac = sum(item["nac_frame_count"] > 0 for item in valid)
        pooled_count = sum(item["nac_frame_count"] for item in valid)
        pooled_denominator = sum(item["frame_count"] for item in valid)
        occupancy = pooled_count / pooled_denominator if pooled_denominator else 0.0
        longest = max((item["longest_nac_event_ps"] for item in valid), default=0.0)
        bound_count = sum(item["bound"] for item in replicas)
        cluster = clusters.get(candidate, {})
        cluster_reproduced = bool(cluster.get("reproduced"))

        failed_gates = []
        if technical_failures:
            failed_gates.append("TECHNICAL_FAILURE")
        if replicas_with_nac < 2:
            failed_gates.append("INSUFFICIENT_NAC_REPLICA_REPRODUCTION")
        if occupancy < 0.01:
            failed_gates.append("POOLED_NAC_OCCUPANCY_BELOW_0.01")
        if longest < 4.0:
            failed_gates.append("LONGEST_NAC_EVENT_BELOW_4PS")
        if bound_count < 3:
            failed_gates.append("UNBOUND_REPLICA")
        if thermo_failures:
            failed_gates.append("THERMODYNAMIC_INSTABILITY")
        if not cluster_reproduced:
            failed_gates.append("NO_REPRODUCED_CROSS_REPLICA_CLUSTER")
        if technical_failures:
            status = "NOT_EVALUATED_TECHNICAL_FAILURE"
        elif failed_gates:
            status = "FAIL_UNRESTRAINED_M1_ENSEMBLE_NAC"
        else:
            status = PASS_STATUS

        advertised_pass = any(
            item["advertised_scientific_status"] == PASS_STATUS for item in replicas
        )
        if advertised_pass and status != PASS_STATUS:
            raise ValueError(
                f"{candidate}: advertised PASS disagrees with independently recomputed primitive occupancy or gates"
            )

        decision = {
            "candidate_id": candidate,
            "scientific_status": status,
            "failed_gates": failed_gates,
            "replicas_with_nac_after_100ps": replicas_with_nac,
            "pooled_nac_frame_count_100_1000ps": pooled_count,
            "pooled_denominator_100_1000ps": pooled_denominator,
            "pooled_nac_occupancy_100_1000ps": occupancy,
            "longest_continuous_nac_ps": longest,
            "bound_replica_count": bound_count,
            "technical_failure_count": technical_failures,
            "thermodynamic_failure_count": thermo_failures,
            "cross_replica_cluster_reproduced": cluster_reproduced,
        }
        decisions.append(decision)
        if status == PASS_STATUS:
            medoid = _validate_medoid(candidate, cluster)
            eligible.append(
                {
                    **decision,
                    "recurrent_nac_cluster_count": int(cluster.get("cluster_count", 0)),
                    "medoid_is_qmmm_optimized_gs": False,
                    "medoid": medoid,
                }
            )
        else:
            failures.append(
                {
                    "candidate_id": candidate,
                    "scientific_status": status,
                    "reasons": failed_gates,
                }
            )

    eligible.sort(
        key=lambda item: (
            -item["replicas_with_nac_after_100ps"],
            -item["pooled_nac_occupancy_100_1000ps"],
            -item["longest_continuous_nac_ps"],
            item["candidate_id"],
        )
    )
    eligible = eligible[:3]
    return {
        "schema_version": 1,
        "technical_status": "PASS",
        "scientific_status": (
            PASS_STATUS if eligible else "FAIL_NO_QMMM_ELIGIBLE_RECURRENT_NAC"
        ),
        "candidate_universe": len(candidates),
        "stageA_expected_replicas": 36,
        "stageB_expected_replicas": expected_stage_b,
        "stageA_selected_conformation_count": len(selected),
        "scientific_pass_count": sum(
            decision["scientific_status"] == PASS_STATUS for decision in decisions
        ),
        "qmmm_eligible_count": len(eligible),
        "qmmm_eligibility": eligible,
        "candidate_decisions": decisions,
        "failure_ledger": failures,
        "scientific_boundary": (
            "Eligibility nominates recurrent classical-MM NAC clusters for QM/MM preflight; "
            "it is not a TS, PMF, barrier, or QM/MM-optimized GS."
        ),
    }


def _write_outputs(result: dict, output_dir: pathlib.Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=False)
    (output_dir / "FINAL_AUDIT.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )
    (output_dir / "qmmm_eligible_top3.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "count": result["qmmm_eligible_count"],
                "candidates": result["qmmm_eligibility"],
                "scientific_boundary": result["scientific_boundary"],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    (output_dir / "failure_ledger.json").write_text(
        json.dumps(result["failure_ledger"], indent=2, sort_keys=True) + "\n"
    )
    with (output_dir / "candidate_gate.tsv").open("w", newline="") as handle:
        fields = [
            "candidate_id",
            "scientific_status",
            "replicas_with_nac_after_100ps",
            "pooled_nac_occupancy_100_1000ps",
            "longest_continuous_nac_ps",
            "bound_replica_count",
            "technical_failure_count",
            "cross_replica_cluster_reproduced",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(result["candidate_decisions"])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=pathlib.Path, required=True)
    parser.add_argument("--output-dir", type=pathlib.Path, required=True)
    args = parser.parse_args()
    result = audit_final(json.loads(args.input.read_text()))
    _write_outputs(result, args.output_dir)
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
