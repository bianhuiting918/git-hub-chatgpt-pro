#!/usr/bin/env python3
"""Audit and salvage the immutable completed double-precision CG EM result."""
import argparse
import hashlib
import json
import pathlib
import re


HARD_PATTERNS = {
    "fatal": r"(?i)fatal(?:\s+error)?",
    "lincs_warning": r"(?i)lincs\s+warning",
    "settle_problem": r"(?i)(?:cannot\s+be\s+settled|settle[^\n]*(?:warning|error|failed))",
    "nan": r"(?i)(?:^|[^A-Za-z0-9_])nan(?:[^A-Za-z0-9_]|$)",
}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def last_float(pattern, text):
    values = re.findall(pattern, text, flags=re.MULTILINE)
    if not values:
        raise ValueError("required log value not found: " + pattern)
    return float(values[-1])


def gro_atom_count(path):
    with path.open(errors="replace") as handle:
        handle.readline()
        return int(handle.readline().strip())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=pathlib.Path, required=True)
    parser.add_argument("--job-id", default="61708647")
    parser.add_argument("--output", type=pathlib.Path)
    args = parser.parse_args()

    root = args.run_root.resolve()
    output = args.output or root / "PASS.json"
    required = [root / name for name in ("run.log", "run.gro", "run.tpr", "run.processed.mdp")]
    missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise SystemExit("missing or empty immutable EM artifacts: " + ", ".join(missing))

    log = (root / "run.log").read_text(errors="replace")
    mdp = (root / "run.processed.mdp").read_text(errors="replace")
    fmax = last_float(r"Maximum force\s*=\s*([+0-9.eE-]+)", log)
    potential = last_float(r"Potential Energy\s*=\s*([+0-9.eE-]+)", log)
    issue_counts = {name: len(re.findall(pattern, log)) for name, pattern in HARD_PATTERNS.items()}
    checks = {
        "finished_mdrun": bool(re.search(r"Finished mdrun", log)),
        "cg_converged": bool(re.search(r"Polak-Ribiere Conjugate Gradients converged to Fmax < 500", log)),
        "fmax_within_threshold": fmax <= 500.0,
        "double_precision": bool(re.search(r"Precision:\s+double", log)),
        "gromacs_2023_1": bool(re.search(r"GROMACS version:\s+2023\.1", log)),
        "cg_integrator": bool(re.search(r"(?m)^integrator\s*=\s*cg\s*$", mdp)),
        "flexible_water": "-DFLEXIBLE" in mdp,
        "no_hard_numerical_issues": not any(issue_counts.values()),
        "atom_count_133589": gro_atom_count(root / "run.gro") == 133589,
    }
    passed = all(checks.values())
    result = {
        "schema_version": 1,
        "candidate_id": "nylc_C18_trueT267_freeGS",
        "stage": "em_cg_flexible_double",
        "status": "PASS_TECHNICAL_EM" if passed else "NOT_EVALUATED_EM_AUDIT_FAIL",
        "source_slurm_job_id": str(args.job_id),
        "scheduler_disposition": "FAILED_POSTPROCESS_ONLY",
        "immutable_run_root": str(root),
        "integrator": "cg",
        "flexible_water": True,
        "precision": "double",
        "gromacs_version": "2023.1",
        "fmax_kj_mol_nm": fmax,
        "threshold_kj_mol_nm": 500.0,
        "potential_energy_kj_mol": potential,
        "atom_count": gro_atom_count(root / "run.gro"),
        "numerical_issue_counts": issue_counts,
        "checks": checks,
        "sha256": {path.name: sha256(path) for path in required},
        "scientific_nac_evidence": False,
        "note": "Technical EM gate only; no restrained or minimized structure is scientific NAC evidence.",
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))
    raise SystemExit(0 if passed else 2)


if __name__ == "__main__":
    main()
