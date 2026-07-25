#!/usr/bin/env python3
"""Aggregate three Stage A replicas per conformation and select at most six."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
from collections import defaultdict

EXPECTED_SEEDS = (26711, 26723, 26737)
UNBOUND_STATUSES = {"FAIL_UNBOUND", "FAIL_SEVERE_CLASH"}


def _candidate_summary(candidate_id: str, audits: list[dict]) -> tuple[dict | None, str | None]:
    seeds = [int(item["velocity_seed"]) for item in audits]
    if len(seeds) != len(set(seeds)):
        raise ValueError(f"{candidate_id}: duplicate velocity seed")
    if set(seeds) != set(EXPECTED_SEEDS):
        raise ValueError(f"{candidate_id}: expected seeds {EXPECTED_SEEDS}, got {sorted(seeds)}")

    technical_failures = sum(item["technical_status"] != "PASS" for item in audits)
    unbound = sum(item["scientific_status"] in UNBOUND_STATUSES for item in audits)
    evaluable = [
        item
        for item in audits
        if item["technical_status"] == "PASS"
        and item["scientific_status"] not in UNBOUND_STATUSES
        and not item["scientific_status"].startswith("NOT_EVALUATED")
    ]
    if not evaluable:
        if unbound == len(audits):
            return None, "FAIL_UNBOUND"
        return None, "NOT_EVALUATED_ALL_REPLICAS_FAILED"

    nac_frames = [
        int(item["nac"]["nac_frame_count_after_20ps"])
        for item in evaluable
    ]
    denominators = [
        int(item["nac"]["analysis_frame_count_after_20ps"])
        for item in evaluable
    ]
    if sum(nac_frames) == 0:
        return None, "FAIL_UNRESTRAINED_NO_NAC_EARLY"

    replica_count = sum(value > 0 for value in nac_frames)
    pooled_denominator = sum(denominators)
    summary = {
        "candidate_id": candidate_id,
        "replica_count_with_nac_after_20ps": replica_count,
        "pooled_nac_frame_count_after_20ps": sum(nac_frames),
        "pooled_denominator_after_20ps": pooled_denominator,
        "pooled_nac_occupancy_after_20ps": (
            sum(nac_frames) / pooled_denominator if pooled_denominator else 0.0
        ),
        "longest_continuous_nac_ps": max(
            float(item["nac"]["longest_event_after_20ps"]["duration_ps"])
            for item in evaluable
        ),
        "reproduced_cluster_count": max(
            int(item.get("cross_replica_cluster_count", 0)) for item in audits
        ),
        "unbound_replica_count": unbound,
        "technical_failure_count": technical_failures,
        "evaluable_replica_count": len(evaluable),
        "replica_seeds": sorted(seeds),
        "replica_scientific_statuses": {
            str(item["velocity_seed"]): item["scientific_status"] for item in audits
        },
    }
    return summary, None


def rank_stage_a(candidate_audits: list[dict], top_n: int = 6) -> dict:
    if not 0 <= top_n <= 6:
        raise ValueError("top_n must be between zero and six")
    grouped = defaultdict(list)
    for audit in candidate_audits:
        grouped[str(audit["candidate_id"])].append(audit)
    if not grouped:
        raise ValueError("no Stage A replica audits")

    ranking = []
    excluded = {}
    for candidate_id in sorted(grouped):
        summary, reason = _candidate_summary(candidate_id, grouped[candidate_id])
        if reason:
            excluded[candidate_id] = reason
        else:
            ranking.append(summary)

    ranking.sort(
        key=lambda item: (
            -item["replica_count_with_nac_after_20ps"],
            -item["pooled_nac_occupancy_after_20ps"],
            -item["longest_continuous_nac_ps"],
            -item["reproduced_cluster_count"],
            item["unbound_replica_count"],
            item["technical_failure_count"],
            item["candidate_id"],
        )
    )
    for index, item in enumerate(ranking, 1):
        item["rank"] = index
    selected = ranking[:top_n]
    return {
        "schema_version": 1,
        "stage": "M1_StageA_100ps_fully_unrestrained",
        "candidate_count": len(grouped),
        "replica_audit_count": len(candidate_audits),
        "ranking": ranking,
        "selected_candidate_ids": [item["candidate_id"] for item in selected],
        "selected_conformations": selected,
        "selection_status": (
            "PASS_STAGEA_TOP6_SELECTED"
            if len(selected) == top_n
            else "NOT_EVALUATED_STAGEB_FEWER_THAN_SIX"
        ),
        "excluded": excluded,
        "scientific_boundary": (
            "Stage A ranking is triage only; Stage B 100-1000 ps gates QM/MM eligibility."
        ),
    }


def _write_outputs(result: dict, audits: list[dict], output_dir: pathlib.Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=False)
    (output_dir / "stageA_replica_audits.json").write_text(
        json.dumps(audits, indent=2, sort_keys=True) + "\n"
    )
    (output_dir / "stageA_top6_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "status": result["selection_status"],
                "selected_candidate_ids": result["selected_candidate_ids"],
                "selected_conformations": result["selected_conformations"],
                "velocity_seeds": list(EXPECTED_SEEDS),
                "stageB_extension_ps": 900.0,
                "stageB_primary_window_ps": [100.0, 1000.0],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    (output_dir / "stageA_failures.json").write_text(
        json.dumps(result["excluded"], indent=2, sort_keys=True) + "\n"
    )
    with (output_dir / "stageA_conformation_ranking.tsv").open("w", newline="") as handle:
        fields = [
            "rank",
            "candidate_id",
            "replica_count_with_nac_after_20ps",
            "pooled_nac_frame_count_after_20ps",
            "pooled_denominator_after_20ps",
            "pooled_nac_occupancy_after_20ps",
            "longest_continuous_nac_ps",
            "reproduced_cluster_count",
            "unbound_replica_count",
            "technical_failure_count",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(result["ranking"])
    (output_dir / "stageA_ranking_summary.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audits", type=pathlib.Path, required=True)
    parser.add_argument("--output-dir", type=pathlib.Path, required=True)
    parser.add_argument("--top-n", type=int, default=6)
    args = parser.parse_args()
    payload = json.loads(args.audits.read_text())
    audits = payload["replica_audits"] if isinstance(payload, dict) else payload
    result = rank_stage_a(audits, top_n=args.top_n)
    _write_outputs(result, audits, args.output_dir)
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
