#!/usr/bin/env python3
"""Run or audit the blind PETase Stage 1/2 gates on a compute host."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path


STATUS_FIELDS = ["gate", "status", "required_before", "command", "output_path", "note"]


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def append_status(rows: list[dict[str, str]], gate: str, status: str, required_before: str, command: str, output_path: Path | str, note: str) -> None:
    rows.append(
        {
            "gate": gate,
            "status": status,
            "required_before": required_before,
            "command": command,
            "output_path": str(output_path),
            "note": note,
        }
    )


def run_command(command: list[str], cwd: Path) -> tuple[int, str]:
    result = subprocess.run(command, cwd=str(cwd), text=True, capture_output=True)
    combined = "\n".join(part for part in [result.stdout.strip(), result.stderr.strip()] if part)
    return result.returncode, combined.replace("\t", " ").replace("\r", " ").replace("\n", " | ")[:1000]


def accepted_gs_pose_count(path: Path) -> int:
    count = 0
    for row in read_tsv(path):
        if row.get("pass_fail") == "pass" and row.get("relaxed_structure_path") not in {"", "pending", None}:
            count += 1
    return count


def ready_ligand_count(path: Path) -> int:
    count = 0
    for row in read_tsv(path):
        if row.get("status") != "built_needs_scissile_candidate_selection":
            continue
        if row.get("pdb_path") in {"", "pending", None}:
            continue
        if row.get("atom_label_path") in {"", "pending", None}:
            continue
        count += 1
    return count


def write_next_actions(path: Path, rows: list[dict[str, str]]) -> None:
    blocked = [row for row in rows if row["status"] == "blocked"]
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Blind Stage 1/2 Gate Next Actions",
        "",
        "This runner records the current blind workflow state from structure-derived inputs.",
        "",
        "Boundary: do not use article trajectories, article transition-state coordinates, article reaction coordinates, article selected CVs, article barriers, or article mechanism conclusions.",
        "",
        "## Current Blocking Items",
        "",
    ]
    if blocked:
        for row in blocked:
            lines.append(f"- {row['gate']}: {row['note']}")
    else:
        lines.append("- None recorded by this runner.")
    lines.extend(
        [
            "",
            "## Required Order",
            "",
            "1. Confirm compute-environment tools and versions.",
            "2. Build ligand conformers and atom-label tables from blind substrate definitions.",
            "3. Run protonation review on cleaned PETase coordinates.",
            "4. Generate blind pose-generation queues from prepared structures and ligand atom labels.",
            "5. Generate and score Stage 1 ground-state poses.",
            "6. Queue Stage 2 classical MD only from accepted ground-state poses.",
            "7. Do not activate Stage 4 QM/MM scans until productive Stage 2 conformers exist.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate(args: argparse.Namespace) -> int:
    project_root = Path(args.project_root)
    phase_root = project_root / args.phase_subdir
    out_root = Path(args.out_root)
    status_dir = out_root / "00_run_status"
    rows: list[dict[str, str]] = []

    scripts_dir = phase_root / "scripts"
    local_work = project_root / "work"
    if not scripts_dir.exists() and local_work.exists():
        scripts_dir = local_work

    if not args.skip_shell_probes:
        probe = scripts_dir / "probe_stage1_compute_environment.sh"
        if not probe.exists():
            append_status(rows, "environment_probe", "blocked", "ligand_build", str(probe), status_dir, "environment probe script missing")
        elif shutil.which("bash") is None:
            append_status(rows, "environment_probe", "blocked", "ligand_build", f"bash {probe}", status_dir, "bash_missing_on_host")
        else:
            probe_out = out_root / "01_system_setup" / "environment_probe"
            code, note = run_command(["bash", str(probe), str(probe_out)], project_root)
            append_status(rows, "environment_probe", "completed" if code == 0 else "failed", "ligand_build", f"bash {probe} {probe_out}", probe_out, note)

    smiles = phase_root / "blind_work" / "01_system_setup" / "ligand_smiles.tsv"
    ligand_out = out_root / "01_system_setup" / "ligand_build"
    ligand_script = scripts_dir / "build_stage1_ligands_rdkit.py"
    if not smiles.exists():
        append_status(rows, "ligand_smiles", "blocked", "ligand_build", str(smiles), smiles, "required blind ligand SMILES table missing")
    elif not ligand_script.exists():
        append_status(rows, "ligand_build", "blocked", "protonation_gate", str(ligand_script), ligand_out, "ligand build script missing")
    else:
        code, note = run_command(
            [
                sys.executable,
                str(ligand_script),
                "--smiles-table",
                str(smiles),
                "--out-dir",
                str(ligand_out),
                "--conformers",
                str(args.conformers),
            ],
            project_root,
        )
        status = "completed" if code == 0 else "blocked"
        append_status(rows, "ligand_build", status, "protonation_gate", f"{sys.executable} {ligand_script}", ligand_out / "ligand_build_manifest.tsv", note)

    pose_script = scripts_dir / "generate_stage1_pose_generation_queue.py"
    prepared_manifest = phase_root / "blind_work" / "01_system_setup" / "prepared_structure_manifest.tsv"
    triad_manifest = phase_root / "blind_work" / "01_system_setup" / "ser_his_asp_triad_candidates.tsv"
    ligand_manifest = ligand_out / "ligand_build_manifest.tsv"
    fallback_ligand_manifest = phase_root / "blind_work" / "01_system_setup" / "ligand_build" / "ligand_build_manifest.tsv"
    if not ligand_manifest.exists() and fallback_ligand_manifest.exists():
        ligand_manifest = fallback_ligand_manifest
    pose_queue_path = out_root / "01_system_setup" / "pose_generation_queue.tsv"
    ready_ligands = ready_ligand_count(ligand_manifest)
    if ready_ligands == 0:
        append_status(rows, "pose_generation_queue", "blocked", "gs_pose_generation", str(pose_script), pose_queue_path, "requires ready_ligand_build_manifest before pose queue generation")
    elif not pose_script.exists():
        append_status(rows, "pose_generation_queue", "blocked", "gs_pose_generation", str(pose_script), pose_queue_path, "pose generation queue script missing")
    else:
        code, note = run_command(
            [
                sys.executable,
                str(pose_script),
                "--prepared-manifest",
                str(prepared_manifest),
                "--triad-manifest",
                str(triad_manifest),
                "--ligand-build-manifest",
                str(ligand_manifest),
                "--out-root",
                str(out_root),
            ],
            project_root,
        )
        append_status(rows, "pose_generation_queue", "completed" if code == 0 else "failed", "gs_pose_generation", f"{sys.executable} {pose_script}", pose_queue_path, note)

    input_pdb = phase_root / "blind_work" / "01_system_setup" / "prepared_initial_pdbs" / args.primary_pdb
    protonation_out = out_root / "01_system_setup" / "protonation_gate"
    protonation_script = scripts_dir / "run_stage1_protonation_gate.sh"
    if args.skip_shell_probes:
        pass
    elif not input_pdb.exists():
        append_status(rows, "protonation_input", "blocked", "protonation_gate", str(input_pdb), input_pdb, "cleaned primary PDB missing")
    elif not protonation_script.exists():
        append_status(rows, "protonation_gate", "blocked", "gs_pose_generation", str(protonation_script), protonation_out, "protonation gate script missing")
    elif shutil.which("bash") is None:
        append_status(rows, "protonation_gate", "blocked", "gs_pose_generation", f"bash {protonation_script} {input_pdb}", protonation_out, "bash_missing_on_host")
    else:
        code, note = run_command(["bash", str(protonation_script), str(input_pdb), str(protonation_out)], project_root)
        append_status(rows, "protonation_gate", "completed" if code == 0 else "failed", "gs_pose_generation", f"bash {protonation_script} {input_pdb}", protonation_out, note)

    gs_manifest = out_root / "01_system_setup" / "gs_pose_manifest.tsv"
    if not gs_manifest.exists():
        gs_manifest = phase_root / "blind_work" / "01_system_setup" / "gs_pose_manifest.tsv"
    accepted_count = accepted_gs_pose_count(gs_manifest)
    if accepted_count == 0:
        append_status(rows, "gs_pose_acceptance", "blocked", "stage2_md_queue", str(gs_manifest), gs_manifest, "no accepted_gs_pose with pass_fail=pass and real relaxed_structure_path")
    else:
        append_status(rows, "gs_pose_acceptance", "completed", "stage2_md_queue", str(gs_manifest), gs_manifest, f"accepted_gs_pose_count={accepted_count}")

    stage2_script = scripts_dir / "generate_stage2_classical_md_manifests.py"
    if accepted_count == 0:
        append_status(rows, "stage2_md_queue", "blocked", "classical_md", str(stage2_script), out_root / "02_classical_md", "requires accepted_gs_pose before MD queue generation")
    elif not stage2_script.exists():
        append_status(rows, "stage2_md_queue", "blocked", "classical_md", str(stage2_script), out_root / "02_classical_md", "Stage 2 generator missing")
    else:
        code, note = run_command(
            [
                sys.executable,
                str(stage2_script),
                "--gs-pose-manifest",
                str(gs_manifest),
                "--out-root",
                str(out_root),
                "--replicates",
                str(args.replicates),
            ],
            project_root,
        )
        append_status(rows, "stage2_md_queue", "completed" if code == 0 else "failed", "classical_md", f"{sys.executable} {stage2_script}", out_root / "02_classical_md" / "md_replicate_queue.tsv", note)

    status_path = status_dir / "stage1_stage2_gate_status.tsv"
    next_actions = status_dir / "stage1_stage2_next_actions.md"
    write_tsv(status_path, rows)
    write_next_actions(next_actions, rows)
    print(status_path)
    print(next_actions)
    return 2 if any(row["status"] == "blocked" for row in rows) else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--phase-subdir", default="project01/phase2_blind_petase_qmmm_20260630")
    parser.add_argument("--out-root", default="project01/phase2_blind_petase_qmmm_20260630/blind_work")
    parser.add_argument("--primary-pdb", default="6EQE_chainA_initial_clean_v2.pdb")
    parser.add_argument("--conformers", type=int, default=20)
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--skip-shell-probes", action="store_true")
    return generate(parser.parse_args())


if __name__ == "__main__":
    raise SystemExit(main())