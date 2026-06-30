#!/usr/bin/env python3
"""Run post-sequence entrance gates for a target manifest.

Rows with missing predicted PDB files are marked NOT_EVALUATED, not FAIL.
Rows with evaluated PDB files receive PASS or FAIL from postseq_entrance_gate.py.
The output manifest is the authoritative input for deciding which samples can
enter the current sequence panel. PLACER and QMMM are deferred for this stage.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path


def read_tsv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--fixed-residues", required=True)
    parser.add_argument("--gate-script", required=True)
    parser.add_argument("--out-manifest", required=True)
    parser.add_argument("--json-dir", required=True)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--require-ligand", action="store_true")
    args = parser.parse_args()

    rows = read_tsv(Path(args.targets))
    out_rows = []
    for row in rows:
        sample_id = row["sample_id"]
        predicted_pdb = Path(row["postseq_predicted_pdb"])
        gate_json = Path(args.json_dir) / f"{sample_id}.postseq_entrance_gate.json"
        out = dict(row)
        out["postseq_entrance_gate_json"] = str(gate_json)
        if not predicted_pdb.exists() or predicted_pdb.stat().st_size == 0:
            out.update(
                {
                    "postseq_structure_status": "NOT_EVALUATED_MISSING_PDB",
                    "entrance_gate_status": "NOT_EVALUATED",
                    "entrance_gate_fail_reasons": "missing_predicted_pdb",
                    "global_backbone_rmsd_A": "",
                    "fixed_backbone_rmsd_A": "",
                    "catalytic_heavy_rmsd_A": "",
                    "candidate_has_ligand": "",
                }
            )
            out_rows.append(out)
            continue

        cmd = [
            args.python,
            args.gate_script,
            "--reference",
            args.reference,
            "--candidate",
            str(predicted_pdb),
            "--fixed-residues",
            args.fixed_residues,
            "--out",
            str(gate_json),
        ]
        if args.require_ligand:
            cmd.append("--require-ligand")
        proc = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if proc.returncode != 0:
            out.update(
                {
                    "postseq_structure_status": "GATE_SCRIPT_ERROR",
                    "entrance_gate_status": "NOT_EVALUATED",
                    "entrance_gate_fail_reasons": f"gate_script_rc_{proc.returncode}",
                    "gate_script_stderr": proc.stderr[-500:].replace("\n", " "),
                    "global_backbone_rmsd_A": "",
                    "fixed_backbone_rmsd_A": "",
                    "catalytic_heavy_rmsd_A": "",
                    "candidate_has_ligand": "",
                }
            )
            out_rows.append(out)
            continue

        result = json.loads(gate_json.read_text())
        metrics = result.get("metrics", {})
        out.update(
            {
                "postseq_structure_status": result.get("evaluation_status", "EVALUATED_WITH_OUTPUT"),
                "entrance_gate_status": result.get("status", "NOT_EVALUATED"),
                "entrance_gate_fail_reasons": ";".join(result.get("fail_reasons", [])),
                "global_backbone_rmsd_A": metrics.get("global_backbone_rmsd_A"),
                "fixed_backbone_rmsd_A": metrics.get("fixed_backbone_rmsd_A"),
                "catalytic_heavy_rmsd_A": metrics.get("catalytic_heavy_rmsd_A"),
                "max_abs_protein_key_distance_delta_A": metrics.get("max_abs_protein_key_distance_delta_A"),
                "candidate_has_ligand": metrics.get("candidate_has_ligand"),
                "fixed_backbone_mean_bfactor_or_plddt": metrics.get("fixed_backbone_mean_bfactor_or_plddt"),
                "catalytic_mean_bfactor_or_plddt": metrics.get("catalytic_mean_bfactor_or_plddt"),
            }
        )
        out_rows.append(out)

    extra_fields = [
        "postseq_entrance_gate_json",
        "entrance_gate_fail_reasons",
        "global_backbone_rmsd_A",
        "fixed_backbone_rmsd_A",
        "catalytic_heavy_rmsd_A",
        "max_abs_protein_key_distance_delta_A",
        "candidate_has_ligand",
        "fixed_backbone_mean_bfactor_or_plddt",
        "catalytic_mean_bfactor_or_plddt",
        "gate_script_stderr",
    ]
    fields = list(rows[0].keys()) if rows else []
    for field in extra_fields:
        if field not in fields:
            fields.append(field)
    write_tsv(Path(args.out_manifest), out_rows, fields)

    counts = {}
    for row in out_rows:
        key = row.get("entrance_gate_status", "NA")
        counts[key] = counts.get(key, 0) + 1
    print(json.dumps({"rows": len(out_rows), "entrance_gate_status_counts": counts}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
