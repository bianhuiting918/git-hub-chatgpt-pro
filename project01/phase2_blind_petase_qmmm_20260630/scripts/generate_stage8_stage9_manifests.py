#!/usr/bin/env python3
"""Generate blind PETase Stage 8 free-energy and Stage 9 validation manifests."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


PMF_FIELDS = [
    "pmf_job_id",
    "accepted_ts_id",
    "path_id",
    "reaction_stage",
    "structure_path",
    "reaction_coordinate_source",
    "sampling_method",
    "window_or_string_plan",
    "analysis_method",
    "sampling_status",
    "barrier_kcal_mol",
    "source",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def accepted_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    accepted = []
    for row in rows:
        if row.get("acceptance_status") != "accepted":
            continue
        if row.get("structure_path") in {"", "pending"}:
            continue
        accepted.append(row)
    return accepted


def pmf_rows(rows: list[dict[str, str]], reaction_stage: str) -> list[dict[str, str]]:
    stage_rows = [row for row in rows if row.get("reaction_stage") == reaction_stage]
    if not stage_rows:
        return [
            {
                "pmf_job_id": "pending",
                "accepted_ts_id": "pending",
                "path_id": "pending",
                "reaction_stage": reaction_stage,
                "structure_path": "pending",
                "reaction_coordinate_source": "not_ready_no_accepted_ts",
                "sampling_method": "pending",
                "window_or_string_plan": "pending",
                "analysis_method": "pending",
                "sampling_status": "not_ready",
                "barrier_kcal_mol": "pending",
                "source": "blind_workflow_pending",
            }
        ]
    output = []
    prefix = "PMF_AC" if reaction_stage == "acylation" else "PMF_DE"
    for index, row in enumerate(stage_rows, start=1):
        output.append(
            {
                "pmf_job_id": f"{prefix}_{index:04d}",
                "accepted_ts_id": row["accepted_ts_id"],
                "path_id": row["path_id"],
                "reaction_stage": reaction_stage,
                "structure_path": row["structure_path"],
                "reaction_coordinate_source": "blind_successful_trajectories_and_candidate_cv_sets",
                "sampling_method": "umbrella_or_string_or_metadynamics_to_be_selected",
                "window_or_string_plan": "pending_from_stage7_endpoint_trajectories",
                "analysis_method": "MBAR_or_WHAM",
                "sampling_status": "not_started",
                "barrier_kcal_mol": "pending",
                "source": "blind_accepted_ts",
            }
        )
    return output


def write_barrier_summary(path: Path, ac_ready: bool, de_ready: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    status = "not_started" if ac_ready or de_ready else "not_ready"
    path.write_text(
        "\n".join(
            [
                "# Blind PETase Free-Energy Barrier Summary",
                "",
                f"Status: {status}",
                "",
                "Boundary: barrier estimates must come from blind PMF/free-energy calculations only.",
                "Do not use article-derived values, rate constants, rate-limiting conclusions, or reaction coordinates.",
                "",
                "## Required Evidence",
                "",
                "- Accepted TS ensemble candidate and committor evidence.",
                "- Reaction coordinate derived from blind trajectories/CV screening.",
                "- Sampling window or string definition with input hashes.",
                "- MBAR/WHAM or equivalent uncertainty estimate.",
                "- Separate acylation and deacylation barrier estimates before any final paper comparison.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def write_final_validation_docs(compare: Path, audit: Path, ready: bool) -> None:
    compare.parent.mkdir(parents=True, exist_ok=True)
    status = "ready_for_blind_vs_paper_after_pmfs" if ready else "not_ready"
    compare.write_text(
        "\n".join(
            [
                "# Blind Vs Paper Comparison",
                "",
                f"Status: {status}",
                "",
                "Do not open paper results until blind TS ensembles and preliminary PMFs are recorded.",
                "",
                "## Compare After Unlock",
                "",
                "- Mechanism sequence and proton relay.",
                "- TS ensemble geometry and diversity.",
                "- Acylation and deacylation barriers with uncertainties.",
                "- Rate-limiting assignment from blind barriers only.",
                "- Disagreements traceable to structure, substrate model, protonation, QM region, sampling, or RC choices.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    audit.write_text(
        "\n".join(
            [
                "# Discrepancy Audit",
                "",
                f"Status: {status}",
                "",
                "not_ready until blind PMF/barrier outputs exist.",
                "",
                "## Audit Axes",
                "",
                "- Structure template and missing-residue handling.",
                "- Substrate model and pose source.",
                "- Protonation state and tautomer branch.",
                "- QM region and boundary treatment.",
                "- Sampling method and convergence.",
                "- Reaction coordinate construction.",
                "- TS validation and committor evidence.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def generate(accepted_manifest: Path, out_root: Path) -> list[Path]:
    rows = accepted_rows(read_tsv(accepted_manifest))
    ac_rows = pmf_rows(rows, "acylation")
    de_rows = pmf_rows(rows, "deacylation")
    ac_ready = any(row["source"] == "blind_accepted_ts" for row in ac_rows)
    de_ready = any(row["source"] == "blind_accepted_ts" for row in de_rows)

    ac_path = out_root / "08_free_energy" / "acylation_pmf.tsv"
    de_path = out_root / "08_free_energy" / "deacylation_pmf.tsv"
    barrier_path = out_root / "08_free_energy" / "barrier_summary.md"
    compare_path = out_root / "09_paper_validation" / "blind_vs_paper_comparison.md"
    audit_path = out_root / "09_paper_validation" / "discrepancy_audit.md"

    write_tsv(ac_path, PMF_FIELDS, ac_rows)
    write_tsv(de_path, PMF_FIELDS, de_rows)
    write_barrier_summary(barrier_path, ac_ready, de_ready)
    write_final_validation_docs(compare_path, audit_path, ac_ready and de_ready)
    return [ac_path, de_path, barrier_path, compare_path, audit_path]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accepted-ts-manifest", default="work/blind_work/07_committor/accepted_ts_manifest.tsv")
    parser.add_argument("--out-root", default="work/blind_work")
    args = parser.parse_args()
    generated = generate(Path(args.accepted_ts_manifest), Path(args.out_root))
    for path in generated:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
