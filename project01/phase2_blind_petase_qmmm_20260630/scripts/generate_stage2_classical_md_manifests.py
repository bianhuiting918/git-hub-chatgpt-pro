#!/usr/bin/env python3
"""Generate blind PETase Stage 2 classical MD manifests from accepted GS poses."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


QUEUE_FIELDS = [
    "md_job_id",
    "source_pose_id",
    "stage",
    "template_pdb",
    "substrate_model",
    "structure_path",
    "replicate_index",
    "ensemble",
    "temperature_K",
    "pressure_bar",
    "planned_equilibration",
    "planned_production",
    "cluster_input",
    "status",
    "source",
]

PRODUCTIVE_FIELDS = [
    "conformer_id",
    "source_pose_id",
    "md_job_id",
    "stage",
    "cluster_id",
    "representative_structure_path",
    "frame_count",
    "geometry_pass_fraction",
    "selection_status",
    "source",
]

REJECTED_FIELDS = [
    "rejected_item_id",
    "source_pose_id",
    "md_job_id",
    "stage",
    "rejection_reason",
    "status",
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


def reaction_stage(pose_id: str) -> str:
    if pose_id.startswith("AC_"):
        return "acylation"
    if pose_id.startswith("DE_"):
        return "deacylation"
    if pose_id.startswith("REACTIVE_"):
        return "preacylation_michaelis_complex"
    return "unknown"


def accepted_gs_pose_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    accepted = []
    for row in rows:
        structure_path = row.get("relaxed_structure_path", "")
        if row.get("pass_fail") != "pass":
            continue
        if structure_path in {"", "pending"}:
            continue
        accepted.append(row)
    return accepted


def queue_rows(accepted_rows: list[dict[str, str]], replicates: int) -> list[dict[str, str]]:
    if not accepted_rows:
        return [
            {
                "md_job_id": "pending",
                "source_pose_id": "pending",
                "stage": "pending",
                "template_pdb": "pending",
                "substrate_model": "pending",
                "structure_path": "pending",
                "replicate_index": "pending",
                "ensemble": "pending",
                "temperature_K": "pending",
                "pressure_bar": "pending",
                "planned_equilibration": "pending",
                "planned_production": "pending",
                "cluster_input": "not_ready_no_accepted_gs_pose",
                "status": "not_ready",
                "source": "blind_workflow_pending",
            }
        ]
    output = []
    for row in accepted_rows:
        pose_id = row["pose_id"]
        stage = reaction_stage(pose_id)
        for replicate in range(1, replicates + 1):
            output.append(
                {
                    "md_job_id": f"MD_{pose_id}_R{replicate:02d}",
                    "source_pose_id": pose_id,
                    "stage": stage,
                    "template_pdb": row.get("template_pdb", ""),
                    "substrate_model": row.get("substrate_model", ""),
                    "structure_path": row.get("relaxed_structure_path", ""),
                    "replicate_index": str(replicate),
                    "ensemble": "NPT_after_NVT_equilibration",
                    "temperature_K": "300",
                    "pressure_bar": "1",
                    "planned_equilibration": "minimize_heat_density_equilibrate",
                    "planned_production": "short_replicate_for_active_site_clustering",
                    "cluster_input": "md_trajectory_after_run",
                    "status": "not_started",
                    "source": "blind_accepted_gs_pose",
                }
            )
    return output


def productive_rows(queue: list[dict[str, str]]) -> list[dict[str, str]]:
    ready_jobs = [row for row in queue if row["source"] == "blind_accepted_gs_pose"]
    if not ready_jobs:
        return [
            {
                "conformer_id": "pending",
                "source_pose_id": "pending",
                "md_job_id": "pending",
                "stage": "pending",
                "cluster_id": "pending",
                "representative_structure_path": "pending",
                "frame_count": "pending",
                "geometry_pass_fraction": "pending",
                "selection_status": "not_ready",
                "source": "awaiting_accepted_gs_pose",
            }
        ]
    return [
        {
            "conformer_id": "pending",
            "source_pose_id": "pending",
            "md_job_id": "pending",
            "stage": "pending",
            "cluster_id": "pending",
            "representative_structure_path": "pending",
            "frame_count": "pending",
            "geometry_pass_fraction": "pending",
            "selection_status": "not_ready",
            "source": "awaiting_md_replicates",
        }
    ]


def rejected_rows() -> list[dict[str, str]]:
    return [
        {
            "rejected_item_id": "pending",
            "source_pose_id": "pending",
            "md_job_id": "pending",
            "stage": "pending",
            "rejection_reason": "pending_md_geometry_clustering",
            "status": "not_ready",
            "source": "blind_workflow_pending",
        }
    ]


def write_protocol(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "# Stage 2 Classical MD Protocol",
                "",
                "Purpose: generate a blind ensemble of PETase active-site conformers from accepted Stage 1 ground-state poses.",
                "",
                "Boundary: do not use article trajectories, article transition-state frames, article reaction coordinates, or article barrier information.",
                "",
                "Gate inputs:",
                "",
                "- Stage 1 GS pose row with `pass_fail=pass`.",
                "- Non-pending relaxed structure path with recorded input provenance.",
                "- Ligand atom labels and protonation branch already reviewed.",
                "",
                "Required workflow:",
                "",
                "1. Minimize, heat, equilibrate, and run independent short MD replicates for each accepted pose.",
                "2. Cluster active-site frames using blind catalytic-geometry descriptors.",
                "3. Score clusters by nucleophile distance, attack angle, His acceptor geometry, and oxyanion stabilization.",
                "4. Write accepted cluster representatives to `productive_conformer_manifest.tsv`.",
                "5. Write rejected clusters or unstable poses to `rejected_pose_manifest.tsv` with explicit reasons.",
                "6. Only productive conformers may feed Stage 4 low-cost QM/MM scans.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def generate(gs_pose_manifest: Path, out_root: Path, replicates: int) -> list[Path]:
    if replicates < 1:
        raise ValueError("--replicates must be at least 1")
    accepted = accepted_gs_pose_rows(read_tsv(gs_pose_manifest))
    queue = queue_rows(accepted, replicates)

    stage_dir = out_root / "02_classical_md"
    queue_path = stage_dir / "md_replicate_queue.tsv"
    productive_path = stage_dir / "productive_conformer_manifest.tsv"
    rejected_path = stage_dir / "rejected_pose_manifest.tsv"
    protocol_path = stage_dir / "stage2_classical_md_protocol.md"

    write_tsv(queue_path, QUEUE_FIELDS, queue)
    write_tsv(productive_path, PRODUCTIVE_FIELDS, productive_rows(queue))
    write_tsv(rejected_path, REJECTED_FIELDS, rejected_rows())
    write_protocol(protocol_path)
    return [queue_path, productive_path, rejected_path, protocol_path]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gs-pose-manifest", default="work/blind_work/01_system_setup/gs_pose_manifest.tsv")
    parser.add_argument("--out-root", default="work/blind_work")
    parser.add_argument("--replicates", type=int, default=3)
    args = parser.parse_args()
    generated = generate(Path(args.gs_pose_manifest), Path(args.out_root), args.replicates)
    for path in generated:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
