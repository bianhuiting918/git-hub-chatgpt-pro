#!/usr/bin/env python3
"""Independent, non-destructive audit of the completed NylC L2 rebalance."""
import argparse
import json
import os
import pathlib
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone


STAGES = ("nvt50", "nvt150", "nvt300", "npt300r", "npt300rel", "npt300free")
FINISHED_MARKER = "Finished " + "md" + "run"
EXPECTED_BACKUP = "run_history.tsv.pre_repair_61708900"
HARD_PATTERNS = {
    "fatal": r"(?i)fatal(?:\s+error)?",
    "lincs_warning": r"(?i)lincs\s+warning",
    "settle_problem": r"(?i)(?:cannot\s+be\s+settled|settle[^\n]*(?:warning|error|failed))",
    "nan": r"(?i)\bnan\b",
}


def run(command, stdout_path, stderr_path, input_text=None, cwd=None):
    with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
        return subprocess.run(
            command,
            input=input_text,
            text=True,
            stdout=stdout,
            stderr=stderr,
            cwd=cwd,
            check=True,
        )


def scan_log(path):
    text = path.read_text(errors="replace")
    counts = {name: len(re.findall(pattern, text)) for name, pattern in HARD_PATTERNS.items()}
    return {
        "finished_mdrun": FINISHED_MARKER in text,
        "numerical_issue_counts": counts,
        "technical_pass": FINISHED_MARKER in text and not any(counts.values()),
    }


def repair_history(task_root, main_job):
    path = task_root / "run_history.tsv"
    text = path.read_text(errors="strict")
    marker = "\trebalance_after_em_double\t"
    index = text.find(marker)
    if index < 0:
        return {"status": "NO_LITERAL_ESCAPE_RECORDS", "changed": False}
    start = text.rfind("\n", 0, index) + 1
    prefix, affected = text[:start], text[start:]
    if f"job={main_job}" not in affected:
        raise RuntimeError("literal escape block does not belong to expected main job")
    if "\t" in prefix or "\n" in prefix:
        raise RuntimeError("literal escapes exist outside the expected repair block")
    backup = task_root / (EXPECTED_BACKUP if str(main_job) == "61708900" else f"run_history.tsv.pre_repair_{main_job}")
    if backup.exists():
        raise RuntimeError(f"refusing to overwrite history backup: {backup}")
    shutil.copy2(path, backup)
    repaired = affected.replace("\t", "	").replace("\n", "
")
    temporary = path.with_suffix(".tsv.repairing")
    temporary.write_text(prefix + repaired)
    temporary.replace(path)
    return {
        "status": "REPAIRED_LITERAL_ESCAPES",
        "changed": True,
        "backup": str(backup),
        "literal_tab_count": affected.count("\t"),
        "literal_newline_count": affected.count("\n"),
    }


def record(task_root, event, state, detail, job):
    stamp = datetime.now(timezone.utc).astimezone().isoformat()
    with (task_root / "run_history.tsv").open("a") as handle:
        handle.write(f"{stamp}	{event}	{state}	candidate=nylc_C18_trueT267_freeGS;job={job};{detail}
")
    payload = {
        "timestamp": stamp,
        "event": event,
        "state": state,
        "detail": detail,
        "candidate_id": "nylc_C18_trueT267_freeGS",
        "slurm_job_id": str(job),
    }
    with (task_root / "run_history.jsonl").open("a") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "
")


def link(source, destination):
    destination.symlink_to(source)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-root", type=pathlib.Path, required=True)
    parser.add_argument("--main-job", default="61708900")
    parser.add_argument("--post-job", default=os.environ.get("SLURM_JOB_ID", "manual"))
    parser.add_argument("--gmx", required=True)
    args = parser.parse_args()

    task = args.task_root.resolve()
    flow = task / "repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
    candidate = task / "candidates/nylc_C18_trueT267_freeGS"
    build = candidate / "build"
    runs = candidate / "runs/rebalance_after_em_double"
    free = runs / "npt300free"
    post = runs / f"postprocess_job_{args.post_job}"
    post.mkdir(parents=True, exist_ok=False)
    record(task, "rebalance_after_em_double_postprocess", "START", f"main_job={args.main_job}", args.post_job)

    history_repair = repair_history(task, args.main_job)
    stage_results = {}
    for stage in STAGES:
        log = runs / stage / "run.log"
        stage_results[stage] = scan_log(log) if log.is_file() else {
            "finished_mdrun": False,
            "numerical_issue_counts": {},
            "technical_pass": False,
            "missing_log": str(log),
        }
    (post / "stage_numerical_audit.json").write_text(
        json.dumps(stage_results, indent=2, sort_keys=True) + "
"
    )

    required = [free / name for name in ("run.xtc", "run.tpr", "run.edr", "run.gro", "run.log")]
    missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    if missing or not all(item["technical_pass"] for item in stage_results.values()):
        summary = {
            "schema_version": 1,
            "candidate_id": "nylc_C18_trueT267_freeGS",
            "main_job_id": args.main_job,
            "postprocess_job_id": str(args.post_job),
            "status": "NOT_EVALUATED_INCOMPLETE_OR_NUMERICAL_FAIL",
            "missing_final_artifacts": missing,
            "history_repair": history_repair,
            "stages": stage_results,
        }
        (post / "postprocess_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "
")
        record(task, "rebalance_after_em_double_postprocess", "NOT_EVALUATED", f"summary={post}/postprocess_summary.json", args.post_job)
        raise SystemExit(2)

    for name in ("run.gro", "run.log"):
        link(free / name, post / name)

    run([args.gmx, "distance", "-f", str(free / "run.xtc"), "-s", str(free / "run.tpr"),
         "-select", "atomnr 8961 plus atomnr 10287", "-oall", str(post / "nac_distance.xvg")],
        post / "nac_distance.stdout", post / "nac_distance.stderr")
    (post / "nac_angle.ndx").write_text("[ nac_angle ]
8961 10287 10288
")
    run([args.gmx, "angle", "-f", str(free / "run.xtc"), "-n", str(post / "nac_angle.ndx"),
         "-ov", str(post / "nac_angle.xvg"), "-type", "angle"],
        post / "nac_angle.stdout", post / "nac_angle.stderr", input_text="0
")
    run([args.gmx, "energy", "-f", str(free / "run.edr"), "-o", str(post / "potential_energy.xvg")],
        post / "potential.stdout", post / "potential.stderr", input_text="Potential
0
")
    run([args.gmx, "energy", "-f", str(free / "run.edr"), "-o", str(post / "thermo.xvg")],
        post / "thermo.stdout", post / "thermo.stderr", input_text="Temperature
Pressure
Volume
0
")
    run([args.gmx, "distance", "-f", str(free / "run.xtc"), "-s", str(free / "run.tpr"),
         "-n", str(build / "source_cycle.ndx"), "-select", 'com of group "Core" plus com of group "Gate"',
         "-oxyz", str(post / "gate_core_vector.xvg")],
        post / "gate.stdout", post / "gate.stderr")

    run([sys.executable, str(flow / "scripts/analyze_nac_series.py"),
         "--distance-xvg", str(post / "nac_distance.xvg"),
         "--angle-xvg", str(post / "nac_angle.xvg"),
         "--potential-xvg", str(post / "potential_energy.xvg"),
         "--distance-max-nm", "0.35", "--angle-min-deg", "95", "--angle-max-deg", "115",
         "--output", str(post / "nac_energy_series.json")],
        post / "series.stdout", post / "series.stderr")
    run([sys.executable, str(flow / "scripts/audit_nylc_true_thr267_l2_free.py"),
         "--run-root", str(post), "--job-id", args.main_job],
        post / "audit.stdout", post / "audit.stderr")

    scientific = json.loads((post / "l2_free_1ns_audit.json").read_text())
    summary = {
        "schema_version": 1,
        "candidate_id": "nylc_C18_trueT267_freeGS",
        "main_job_id": args.main_job,
        "postprocess_job_id": str(args.post_job),
        "status": "PASS_POSTPROCESS_TECHNICAL",
        "history_repair": history_repair,
        "gate_definition": "NylC residues 261-266; Thr267 excluded",
        "scientific_audit": scientific,
        "stages": stage_results,
    }
    (post / "postprocess_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "
")
    record(task, "rebalance_after_em_double_postprocess", "PASS_TECHNICAL",
           f"main_job={args.main_job};scientific_audit={post}/l2_free_1ns_audit.json;gate=residues 261-266 excluding Thr267",
           args.post_job)
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
