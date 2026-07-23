#!/usr/bin/env python3
"""Collect Amber single-point outputs without confusing exit zero for PASS."""
from __future__ import annotations
import argparse, json, math, re
from pathlib import Path

def parse_sander_output(path: Path) -> dict[str, object]:
    text = path.read_text(errors="replace")
    def finite_field(name: str) -> float | None:
        matches = re.findall(rf"{name}\s*=\s*([-+0-9.Ee]+)", text)
        if not matches:
            return None
        value = float(matches[-1])
        return value if math.isfinite(value) else None
    def final_results_energy() -> float | None:
        idx = text.rfind("FINAL RESULTS")
        if idx < 0:
            return None
        block = text[idx:].splitlines()
        for i, line in enumerate(block):
            if "NSTEP" in line and "ENERGY" in line:
                for data in block[i+1:]:
                    nums = re.findall(r"[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?", data)
                    if len(nums) >= 2:
                        value = float(nums[1])
                        return value if math.isfinite(value) else None
        return None
    scc = bool(re.search(r"SCC-DFTB for step\s+\d+\s+converged in\s+\d+\s+cycles", text)) or "SCC CONVERGED" in text or "DFTB SCC converged" in text
    total = finite_field("TOTAL_ENERGY")
    if total is None:
        total = final_results_energy()
    return {
        "normal_termination": "FINAL RESULTS" in text or "Normal termination" in text,
        "scc_converged": scc,
        "DFTBESCF": finite_field("DFTBESCF"),
        "EGB": finite_field("EGB"),
        "TOTAL_ENERGY": total,
    }

def geometry_report_pass(path: Path) -> bool:
    report = json.loads(path.read_text())
    return report.get("pass") is True or report.get("status") == "PASS"

def write_pass_json(log: Path, out: Path, process_exit: int, geometry_pass: bool) -> None:
    parsed = parse_sander_output(log)
    required = [parsed["normal_termination"], parsed["scc_converged"], parsed["DFTBESCF"] is not None, parsed["EGB"] is not None, parsed["TOTAL_ENERGY"] is not None, geometry_pass]
    if process_exit != 0 or not all(required):
        raise ValueError("PASS.json requires normal Amber termination, converged SCC, finite energies, and geometry PASS; exit 0 alone is insufficient")
    out.write_text(json.dumps({"status": "PASS", **parsed}, indent=2, sort_keys=True) + "\n")

def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--log", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--process-exit", required=True, type=int)
    p.add_argument("--geometry-report", required=True, type=Path)
    args = p.parse_args(argv)
    write_pass_json(args.log, args.out, args.process_exit, geometry_report_pass(args.geometry_report))
    return 0
if __name__ == "__main__":
    raise SystemExit(main())
