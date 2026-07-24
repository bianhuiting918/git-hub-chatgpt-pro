#!/usr/bin/env python3
"""Parse fixed-coordinate Amber QM/MM single-point outputs and summarize TS-GS gaps."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
from pathlib import Path


FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"


def read_rst7_coordinates(path: Path) -> list[float]:
    lines = path.read_text(errors="replace").splitlines()
    if len(lines) < 2:
        raise ValueError(f"invalid rst7: {path}")
    natom = int(lines[1].split()[0])
    values: list[float] = []
    for line in lines[2:]:
        for token in line.split():
            try:
                values.append(float(token))
            except ValueError:
                pass
            if len(values) >= 3 * natom:
                return values
    raise ValueError(f"not enough coordinates in {path}")


def parse_final_energy(text: str) -> float:
    marker = text.rfind("FINAL RESULTS")
    region = text[marker:] if marker >= 0 else text
    match = re.search(
        rf"(?m)^\s*NSTEP\s+ENERGY\s+RMS\s+GMAX.*?^\s*\d+\s+({FLOAT})\s+{FLOAT}\s+{FLOAT}",
        region,
        re.S,
    )
    if match:
        return float(match.group(1))
    candidates = re.findall(rf"(?m)^\s*ENERGY\s*=\s*({FLOAT})", region)
    if candidates:
        return float(candidates[-1])
    raise ValueError("final minimization ENERGY not found")


def parse_cases(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def analyze(task: Path) -> dict:
    cases = parse_cases(task / "input" / "cases.tsv")
    rows = []
    for case in cases:
        work = task / "work" / case["case_id"]
        out = work / "sp.out"
        stderr = work / "launcher.stderr"
        text = out.read_text(errors="replace") if out.exists() else ""
        err_text = stderr.read_text(errors="replace") if stderr.exists() else ""
        hard_patterns = {
            "amber_error": len(re.findall(r"(?i)\bERROR\b|FATAL", text + "\n" + err_text)),
            "nan": len(re.findall(r"(?i)\bNaN\b", text + "\n" + err_text)),
            "scc_failure": len(re.findall(r"(?i)SCC.*(?:fail|not converg)", text + "\n" + err_text)),
            "vlimit": len(re.findall(r"(?i)vlimit", text + "\n" + err_text)),
        }
        try:
            energy = parse_final_energy(text)
            energy_status = "PASS"
        except Exception as exc:
            energy = math.nan
            energy_status = f"FAIL:{type(exc).__name__}:{exc}"

        output_rst = work / "sp.rst7"
        coord_delta = math.nan
        if output_rst.exists():
            before = read_rst7_coordinates(Path(case["coordinate"]))
            after = read_rst7_coordinates(output_rst)
            if len(before) == len(after):
                coord_delta = max(abs(a - b) for a, b in zip(before, after))
        row = {
            **case,
            "energy_kcal_mol": energy,
            "energy_status": energy_status,
            "max_coordinate_delta_angstrom": coord_delta,
            **hard_patterns,
            "output": str(out),
        }
        rows.append(row)

    out_tsv = task / "audit" / "singlepoint_energies.tsv"
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with out_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    steps = {}
    for step in ("step1", "step2"):
        subset = [row for row in rows if row["step"] == step]
        ts = [row for row in subset if row["state"] == "TS" and math.isfinite(row["energy_kcal_mol"])]
        gs = [row for row in subset if row["state"] == "GS" and math.isfinite(row["energy_kcal_mol"])]
        if len(ts) != 1 or not gs:
            steps[step] = {"status": "NOT_EVALUATED", "reason": "missing finite TS or GS energies"}
            continue
        ets = ts[0]["energy_kcal_mol"]
        gaps = [ets - row["energy_kcal_mol"] for row in gs]
        steps[step] = {
            "status": "DIAGNOSTIC_ONLY_FIXED_COORDINATE_DELTA_E",
            "ts_case": ts[0]["case_id"],
            "ts_energy_kcal_mol": ets,
            "gs_case_count": len(gs),
            "gs_energy_mean_kcal_mol": statistics.mean(row["energy_kcal_mol"] for row in gs),
            "gs_energy_sd_kcal_mol": statistics.stdev(row["energy_kcal_mol"] for row in gs) if len(gs) > 1 else 0.0,
            "delta_e_vs_each_gs_kcal_mol": {
                row["case_id"]: ets - row["energy_kcal_mol"] for row in gs
            },
            "delta_e_mean_kcal_mol": statistics.mean(gaps),
            "delta_e_sd_kcal_mol": statistics.stdev(gaps) if len(gaps) > 1 else 0.0,
            "delta_e_min_kcal_mol": min(gaps),
            "delta_e_max_kcal_mol": max(gaps),
        }

    result = {
        "schema_version": 1,
        "method": "Amber sander QM/MM fixed-coordinate energy; DFTB3 using the installed Amber DFTB parameterization",
        "interpretation": "Potential-energy diagnostic only; not a free-energy barrier and not a stationary-point frequency result.",
        "all_cases_finite": all(math.isfinite(row["energy_kcal_mol"]) for row in rows),
        "all_coordinates_unchanged": all(
            math.isfinite(row["max_coordinate_delta_angstrom"])
            and row["max_coordinate_delta_angstrom"] <= 1e-8
            for row in rows
        ),
        "hard_error_total": sum(
            row["amber_error"] + row["nan"] + row["scc_failure"] for row in rows
        ),
        "steps": steps,
    }
    (task / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.task)


if __name__ == "__main__":
    main()
