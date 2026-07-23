#!/usr/bin/env python3
"""Strict per-candidate grompp preflight with failure isolation and audits."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def grompp_command(gmx: str, mdp: Path, build: Path, out: Path) -> List[str]:
    return [
        gmx,
        "grompp",
        "-f",
        str(mdp),
        "-c",
        str(build / "rebuilt.gro"),
        "-r",
        str(build / "rebuilt.gro"),
        "-p",
        str(build / "topol.top"),
        "-o",
        str(out / "em.tpr"),
        "-po",
        str(out / "em.processed.mdp"),
        "-pp",
        str(out / "processed.top"),
        "-maxwarn",
        "0",
    ]


def scan_diagnostics(text: str) -> Dict[str, int]:
    patterns = {
        "fatal": r"(?i)fatal(?:\s+error)?",
        "lincs_warning": r"(?i)lincs\s+warning",
        "settle_problem": r"(?i)(?:cannot\s+be\s+settled|settle[^\n]*(?:warning|error|failed))",
        "nan": r"(?i)\bnan\b",
    }
    return {name: len(re.findall(pattern, text)) for name, pattern in patterns.items()}


def _load_prepare_rows(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _record(task_root: Path, event: str, state: str, detail: Dict[str, Any]) -> None:
    stamp = datetime.now(timezone.utc).astimezone().isoformat()
    tsv = task_root / "run_history.tsv"
    with tsv.open("a", encoding="utf-8") as handle:
        handle.write(
            "\t".join(
                [
                    stamp,
                    event,
                    state,
                    json.dumps(detail, sort_keys=True, separators=(",", ":")),
                ]
            )
            + "\n"
        )
    jsonl = task_root / "run_history.jsonl"
    payload = {
        "timestamp": stamp,
        "event": event,
        "state": state,
        "detail": detail,
    }
    with jsonl.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")


def preflight_all(task_root: Path, flow_root: Path, gmx: str) -> Dict[str, Any]:
    prepare_rows = _load_prepare_rows(task_root / "audit" / "prepare_summary.tsv")
    mdp = flow_root / "mdp" / "em.mdp"
    rows = []
    for prepare in prepare_rows:
        candidate_id = prepare["candidate_id"]
        if prepare["state"] != "BUILD_PASS":
            rows.append(
                {
                    "candidate_id": candidate_id,
                    "state": "NOT_EVALUATED_BUILD_FAIL",
                    "returncode": "",
                    "reason": prepare["reason"],
                }
            )
            continue
        candidate_root = task_root / "candidates" / candidate_id
        build = candidate_root / "build"
        out = candidate_root / "preflight"
        out.mkdir(parents=True, exist_ok=True)
        stdout_path = out / "grompp.stdout"
        stderr_path = out / "grompp.stderr"
        command = grompp_command(gmx, mdp, build, out)
        _record(task_root, "grompp_preflight", "START", {
            "candidate_id": candidate_id,
            "command": command,
        })
        result = subprocess.run(
            command,
            cwd=build,
            text=True,
            capture_output=True,
            check=False,
        )
        stdout_path.write_text(result.stdout, encoding="utf-8")
        stderr_path.write_text(result.stderr, encoding="utf-8")
        diagnostics = scan_diagnostics(result.stdout + "\n" + result.stderr)
        outputs_exist = (out / "em.tpr").is_file() and (out / "em.tpr").stat().st_size > 0
        state = (
            "GROMPP_PASS"
            if result.returncode == 0
            and outputs_exist
            and diagnostics["fatal"] == 0
            and diagnostics["nan"] == 0
            else "GROMPP_FAIL"
        )
        audit = {
            "schema_version": 1,
            "candidate_id": candidate_id,
            "state": state,
            "returncode": result.returncode,
            "command": command,
            "diagnostics": diagnostics,
            "outputs_exist": outputs_exist,
            "input_hashes": {
                "gro_sha256": sha256(build / "rebuilt.gro"),
                "topology_sha256": sha256(build / "topol.top"),
                "l2_itp_sha256": sha256(build / "PA66_L2_GMX.itp"),
                "em_mdp_sha256": sha256(mdp),
            },
            "output_hashes": {
                "em_tpr_sha256": sha256(out / "em.tpr") if outputs_exist else None,
            },
            "stderr_tail": result.stderr.splitlines()[-40:],
        }
        (out / "grompp_preflight.json").write_text(
            json.dumps(audit, indent=2) + "\n", encoding="utf-8"
        )
        reason = "" if state == "GROMPP_PASS" else (
            f"grompp_rc={result.returncode}; diagnostics={diagnostics}"
        )
        rows.append(
            {
                "candidate_id": candidate_id,
                "state": state,
                "returncode": str(result.returncode),
                "reason": reason,
            }
        )
        _record(task_root, "grompp_preflight", state, {
            "candidate_id": candidate_id,
            "returncode": result.returncode,
            "diagnostics": diagnostics,
        })

    summary = task_root / "audit" / "grompp_summary.tsv"
    lines = ["candidate_id\tstate\treturncode\treason"]
    lines.extend(
        "\t".join(
            row[key].replace("\t", " ").replace("\n", " ")
            for key in ("candidate_id", "state", "returncode", "reason")
        )
        for row in rows
    )
    summary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {
        "candidate_count": len(rows),
        "grompp_pass_count": sum(row["state"] == "GROMPP_PASS" for row in rows),
        "rows": rows,
        "summary": str(summary),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-root", type=Path, required=True)
    parser.add_argument("--flow-root", type=Path, required=True)
    parser.add_argument("--gmx", required=True)
    args = parser.parse_args()
    print(json.dumps(preflight_all(args.task_root, args.flow_root, args.gmx), indent=2))


if __name__ == "__main__":
    main()
