#!/usr/bin/env python3
"""Merge all nine A1 NAC audits without dropping failed or missing slots."""

from __future__ import annotations
import argparse
import csv
import json
import math
import pathlib
import statistics
from collections import defaultdict

def load_records(universe:dict,audit_root:pathlib.Path,audit_array_job_id:str)->tuple[list[dict],list[int]]:
    records=[]; missing_slots=[]
    expected_slots=list(range(9))
    rows={int(row["slot"]):row for row in universe["replicas"]}
    if sorted(rows)!=expected_slots:
        raise ValueError(f"expected_slots mismatch: {sorted(rows)}")
    for slot in expected_slots:
        row=rows[slot]
        base_pattern=f"attempt_{audit_array_job_id}_{slot}_*/{row['candidate_id']}/seed{row['velocity_seed']}"
        roots=sorted(audit_root.glob(base_pattern))
        if len(roots)>1:
            raise ValueError(f"slot {slot} has multiple audit roots: {roots}")
        if not roots:
            missing_slots.append(slot)
            record={"schema_version":1,"slot":slot,"candidate_id":row["candidate_id"],"velocity_seed":row["velocity_seed"],"technical_status":"FAIL","scientific_status":"NOT_EVALUATED_MISSING_AUDIT","reason":"no audit output directory"}
        else:
            root=roots[0]; audit=root/"replica_audit.json"; failed=root/"NOT_EVALUATED.json"
            if audit.is_file():
                record=json.loads(audit.read_text())
            elif failed.is_file():
                record=json.loads(failed.read_text())
            else:
                missing_slots.append(slot)
                record={"schema_version":1,"slot":slot,"candidate_id":row["candidate_id"],"velocity_seed":row["velocity_seed"],"technical_status":"FAIL","scientific_status":"NOT_EVALUATED_MISSING_AUDIT","reason":"audit root lacks result record"}
            record["audit_root"]=str(root)
        record["slot"]=slot
        records.append(record)
    if len(records)!=9 or [r["slot"] for r in records]!=expected_slots:
        raise ValueError("replica_denominator is not exactly nine ordered slots")
    return records,missing_slots

def candidate_summaries(records:list[dict])->list[dict]:
    grouped=defaultdict(list)
    for record in records: grouped[record["candidate_id"]].append(record)
    output=[]
    for candidate in sorted(grouped):
        rows=sorted(grouped[candidate],key=lambda x:x["velocity_seed"])
        eligible=[r for r in rows if r.get("scientific_status")=="PASS_REPLICA_NAC_PRESENT"]
        occupancies=[float(r["nac"]["nac_occupancy"]) for r in eligible]
        longest=[float(r["nac"]["longest_event"]["duration_ps"]) for r in eligible]
        potentials=[float(r["potential_energy_kj_mol"]["mean"]) for r in eligible if r.get("potential_energy_kj_mol")]
        output.append({
            "candidate_id":candidate,
            "candidate_denominator":3,
            "technical_pass_count":sum(r.get("technical_status")=="PASS" for r in rows),
            "eligible_replica_count":len(eligible),
            "nac_positive_seed_count":sum(r.get("nac",{}).get("nac_frame_count",0)>0 for r in rows),
            "mean_eligible_nac_occupancy":statistics.fmean(occupancies) if occupancies else None,
            "max_longest_continuous_nac_ps":max(longest) if longest else 0.0,
            "mean_eligible_potential_energy_kj_mol":statistics.fmean(potentials) if potentials else None,
            "replica_statuses":[{"slot":r["slot"],"velocity_seed":r["velocity_seed"],"technical_status":r.get("technical_status"),"scientific_status":r.get("scientific_status")} for r in rows],
        })
    if len(output)!=3:
        raise ValueError("candidate denominator is not exactly three")
    return output

def select_replica(records:list[dict],candidate_rows:list[dict]):
    pass_count={r["candidate_id"]:r["eligible_replica_count"] for r in candidate_rows}
    eligible=[r for r in records if r.get("scientific_status")=="PASS_REPLICA_NAC_PRESENT"]
    def key(r):
        potential=r.get("potential_energy_kj_mol") or {}
        energy=float(potential.get("mean",math.inf))
        return (
            -pass_count[r["candidate_id"]],
            -float(r["nac"]["nac_occupancy"]),
            -float(r["nac"]["longest_event"]["duration_ps"]),
            energy,
            int(r["slot"]),
        )
    return min(eligible,key=key) if eligible else None

def write_tsv(path:pathlib.Path,records:list[dict])->None:
    columns=["slot","candidate_id","velocity_seed","technical_status","scientific_status","frame_count","nac_frame_count","nac_occupancy","longest_nac_ps","potential_mean_kj_mol","gate_mean_nm","temperature_mean_K","pressure_mean_bar"]
    with path.open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=columns,delimiter="\t"); writer.writeheader()
        for r in records:
            nac=r.get("nac",{}); potential=r.get("potential_energy_kj_mol") or {}; gate=r.get("gate_opening_nm") or {}; thermo=r.get("thermodynamics") or {}
            writer.writerow({
                "slot":r["slot"],"candidate_id":r["candidate_id"],"velocity_seed":r["velocity_seed"],
                "technical_status":r.get("technical_status"),"scientific_status":r.get("scientific_status"),
                "frame_count":(r.get("analysis") or {}).get("frame_count"),
                "nac_frame_count":nac.get("nac_frame_count"),"nac_occupancy":nac.get("nac_occupancy"),
                "longest_nac_ps":(nac.get("longest_event") or {}).get("duration_ps"),
                "potential_mean_kj_mol":potential.get("mean"),"gate_mean_nm":gate.get("mean"),
                "temperature_mean_K":(thermo.get("temperature_K") or {}).get("mean"),
                "pressure_mean_bar":(thermo.get("pressure_bar") or {}).get("mean"),
            })

def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--universe",type=pathlib.Path,required=True)
    ap.add_argument("--audit-root",type=pathlib.Path,required=True)
    ap.add_argument("--audit-array-job-id",required=True)
    ap.add_argument("--output-dir",type=pathlib.Path,required=True)
    args=ap.parse_args()
    if args.output_dir.exists(): raise SystemExit(f"refusing to overwrite {args.output_dir}")
    args.output_dir.mkdir(parents=True)
    universe=json.loads(args.universe.read_text())
    records,missing_slots=load_records(universe,args.audit_root,args.audit_array_job_id)
    candidates=candidate_summaries(records)
    selected=select_replica(records,candidates)
    eligible_replica_count=sum(r.get("scientific_status")=="PASS_REPLICA_NAC_PRESENT" for r in records)
    technical_pass_count=sum(r.get("technical_status")=="PASS" for r in records)
    if selected:
        scientific_status="PASS_A1_NAC_CANDIDATE_AVAILABLE"
    elif technical_pass_count==0:
        scientific_status="NOT_EVALUATED_NO_TECHNICAL_AUDITS"
    else:
        scientific_status="FAIL_NO_ELIGIBLE_A1_NAC"
    summary={
        "schema_version":1,
        "status":"PASS_TECHNICAL_A1_NAC_MERGE" if not missing_slots else "PASS_WITH_NOT_EVALUATED_MISSING_AUDIT",
        "scientific_status":scientific_status,
        "audit_array_job_id":str(args.audit_array_job_id),
        "expected_slots":list(range(9)),
        "missing_slots":missing_slots,
        "replica_denominator":9,
        "candidate_denominator":3,
        "eligible_replica_count":eligible_replica_count,
        "technical_pass_count":technical_pass_count,
        "selection_order":["candidate eligible seed count descending","replica nac_occupancy descending","longest continuous NAC descending","potential_energy_kj_mol mean ascending"],
        "selected_replica":selected,
        "candidate_summaries":candidates,
        "scientific_scope":"classical-MM A1 NAC screening; potential energy is not a QM/MM barrier",
    }
    (args.output_dir/"replica_audits.json").write_text(json.dumps(records,indent=2,sort_keys=True)+"\n")
    (args.output_dir/"candidate_summaries.json").write_text(json.dumps(candidates,indent=2,sort_keys=True)+"\n")
    (args.output_dir/"A1_NAC_AUDIT_SUMMARY.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    write_tsv(args.output_dir/"replica_summary.tsv",records)
    print(json.dumps(summary,sort_keys=True))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
