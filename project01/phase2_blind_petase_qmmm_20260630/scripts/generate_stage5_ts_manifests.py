#!/usr/bin/env python3
"""Generate blind PETase Stage 5/6 TS refinement and ensemble manifests."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


REFINEMENT_FIELDS = [
    "ts_guess_id",
    "path_id",
    "reaction_stage",
    "source_pose_id",
    "source_structure_path",
    "refined_structure_path",
    "qm_region_version",
    "refinement_method",
    "imaginary_mode_status",
    "endpoint_check_status",
    "refinement_status",
    "rejection_reason",
    "source",
]

ENSEMBLE_FIELDS = [
    "ensemble_candidate_id",
    "ts_guess_id",
    "path_id",
    "reaction_stage",
    "source_pose_id",
    "refined_structure_path",
    "cluster_id",
    "ensemble_status",
    "diversity_notes",
    "next_step",
]

COMMITTOR_FIELDS = [
    "committor_job_id",
    "ensemble_candidate_id",
    "path_id",
    "reaction_stage",
    "structure_path",
    "shooting_velocity_seed_plan",
    "target_trial_count",
    "committor_status",
    "notes",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def reaction_stage_from_path(path_id: str) -> str:
    if path_id.startswith("AC_"):
        return "acylation"
    if path_id.startswith("DE_"):
        return "deacylation"
    return "unknown"


def accepted_guesses(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    accepted = []
    for row in rows:
        if row.get("ts_guess_id") in {"", "pending"}:
            continue
        if row.get("structure_path") in {"", "pending"}:
            continue
        if row.get("validation_status") != "endpoint_checked":
            continue
        accepted.append(row)
    return accepted


def build_refinement_rows(guesses: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = []
    for guess in guesses:
        stage = reaction_stage_from_path(guess["path_id"])
        rows.append(
            {
                "ts_guess_id": guess["ts_guess_id"],
                "path_id": guess["path_id"],
                "reaction_stage": stage,
                "source_pose_id": guess["source_pose_id"],
                "source_structure_path": guess["structure_path"],
                "refined_structure_path": "pending",
                "qm_region_version": "pending_from_stage4_qm_region",
                "refinement_method": "qmmm_ts_optimization_or_constrained_saddle_refinement",
                "imaginary_mode_status": "not_checked",
                "endpoint_check_status": guess["validation_status"],
                "refinement_status": "not_started",
                "rejection_reason": "none",
                "source": "blind_stage4_ts_like_guess",
            }
        )
    return rows


def build_ensemble_rows(refinement_rows: list[dict[str, str]], reaction_stage: str) -> list[dict[str, str]]:
    rows = []
    prefix = "AC_TS" if reaction_stage == "acylation" else "DE_TS"
    stage_rows = [row for row in refinement_rows if row["reaction_stage"] == reaction_stage]
    for index, row in enumerate(stage_rows, start=1):
        rows.append(
            {
                "ensemble_candidate_id": f"{prefix}_{index:04d}",
                "ts_guess_id": row["ts_guess_id"],
                "path_id": row["path_id"],
                "reaction_stage": reaction_stage,
                "source_pose_id": row["source_pose_id"],
                "refined_structure_path": row["refined_structure_path"],
                "cluster_id": "pending",
                "ensemble_status": "candidate_pending_refinement",
                "diversity_notes": "preserve_until_refinement_and_clustering",
                "next_step": "run_ts_refinement_then_imaginary_mode_and_endpoint_checks",
            }
        )
    if not rows:
        rows.append(
            {
                "ensemble_candidate_id": "pending",
                "ts_guess_id": "pending",
                "path_id": "pending",
                "reaction_stage": reaction_stage,
                "source_pose_id": "pending",
                "refined_structure_path": "pending",
                "cluster_id": "pending",
                "ensemble_status": "not_generated",
                "diversity_notes": "no_endpoint_checked_stage4_guess_available",
                "next_step": "generate_stage4_ts_like_guess",
            }
        )
    return rows


def build_committor_queue(ac_rows: list[dict[str, str]], de_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    queue = []
    for row in ac_rows + de_rows:
        if row["ensemble_candidate_id"] == "pending":
            continue
        queue.append(
            {
                "committor_job_id": f"COM_{row['ensemble_candidate_id']}",
                "ensemble_candidate_id": row["ensemble_candidate_id"],
                "path_id": row["path_id"],
                "reaction_stage": row["reaction_stage"],
                "structure_path": row["refined_structure_path"],
                "shooting_velocity_seed_plan": "pending_after_refined_ts_structure_exists",
                "target_trial_count": "pending",
                "committor_status": "not_ready_refinement_required",
                "notes": "queue entry created from blind TS-like guess; do not run committor before refinement checks",
            }
        )
    if not queue:
        queue.append(
            {
                "committor_job_id": "pending",
                "ensemble_candidate_id": "pending",
                "path_id": "pending",
                "reaction_stage": "pending",
                "structure_path": "pending",
                "shooting_velocity_seed_plan": "pending",
                "target_trial_count": "pending",
                "committor_status": "not_generated",
                "notes": "no TS ensemble candidates available",
            }
        )
    return queue


def generate(ts_like_manifest: Path, out_root: Path) -> list[Path]:
    guesses = accepted_guesses(read_tsv(ts_like_manifest))
    refinement_rows = build_refinement_rows(guesses)
    ac_rows = build_ensemble_rows(refinement_rows, "acylation")
    de_rows = build_ensemble_rows(refinement_rows, "deacylation")
    committor_rows = build_committor_queue(ac_rows, de_rows)

    refinement_path = out_root / "05_ts_refinement" / "ts_refinement_manifest.tsv"
    ac_path = out_root / "06_ts_ensemble" / "acylation_ts_ensemble.tsv"
    de_path = out_root / "06_ts_ensemble" / "deacylation_ts_ensemble.tsv"
    committor_path = out_root / "07_committor" / "committor_queue.tsv"

    write_tsv(refinement_path, REFINEMENT_FIELDS, refinement_rows)
    write_tsv(ac_path, ENSEMBLE_FIELDS, ac_rows)
    write_tsv(de_path, ENSEMBLE_FIELDS, de_rows)
    write_tsv(committor_path, COMMITTOR_FIELDS, committor_rows)
    return [refinement_path, ac_path, de_path, committor_path]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ts-like-manifest", default="work/blind_work/04_qmmm_exploration/ts_like_guess_manifest.tsv")
    parser.add_argument("--out-root", default="work/blind_work")
    args = parser.parse_args()
    generated = generate(Path(args.ts_like_manifest), Path(args.out_root))
    for path in generated:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
