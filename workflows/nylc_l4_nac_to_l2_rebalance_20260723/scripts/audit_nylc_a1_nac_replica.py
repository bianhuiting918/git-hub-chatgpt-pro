#!/usr/bin/env python3
"""Audit one fully unrestrained NylC A1 replica without proton-path inference."""

from __future__ import annotations
import argparse
import json
import math
import pathlib
import statistics

from analyze_nac_series import audit_series, read_xvg
from audit_nylc_m1_ensemble_replica import (
    _median_spacing,
    _read_table,
    _scan_log,
    _stats,
    _thermodynamics,
    _validate_times,
    _window,
)

MICROSTATE="A1_Thr267_Ogamma_minus_NalphaH3_plus"

def audit_replica(run_root:pathlib.Path,manifest:dict,window:tuple[float,float])->dict:
    root=pathlib.Path(run_root)
    start,end=map(float,window)
    if not manifest.get("fully_unrestrained"):
        raise ValueError("replica manifest is not fully unrestrained")
    if manifest.get("microstate")!=MICROSTATE:
        raise ValueError("replica manifest is not the A1 microstate")
    contract=manifest["nac_contract"]
    distance_max_nm=float(contract["distance_max_nm"])
    angle_min_deg=float(contract["angle_min_deg"])
    angle_max_deg=float(contract["angle_max_deg"])
    distance_all=read_xvg(root/"nac_distance.xvg")
    angle_all=read_xvg(root/"nac_angle.xvg")
    if len(distance_all)<2:
        raise ValueError("at least two NAC frames are required")
    sample_interval=statistics.median(r[0]-l[0] for l,r in zip(distance_all,distance_all[1:]))
    nac=audit_series(
        distance_all,angle_all,distance_max_nm,angle_min_deg,angle_max_deg,
        analysis_start_ps=start,analysis_end_ps=end,sample_interval_ps=sample_interval,
    )
    after_start=min(end,start+20.0)
    after=audit_series(
        distance_all,angle_all,distance_max_nm,angle_min_deg,angle_max_deg,
        analysis_start_ps=after_start,analysis_end_ps=end,sample_interval_ps=sample_interval,
    )
    nac["analysis_frame_count_after_20ps"]=after["frame_count"]
    nac["nac_frame_count_after_20ps"]=after["nac_frame_count"]
    nac["nac_occupancy_after_20ps"]=after["nac_occupancy"]
    nac["longest_event"]=nac.pop("longest_continuous_nac") or {"duration_ps":0.0,"frame_count":0}
    nac["longest_event_after_20ps"]=after["longest_continuous_nac"] or {"duration_ps":0.0,"frame_count":0}

    distance=_window(_read_table(root/"nac_distance.xvg",2),start,end)
    angle=_window(_read_table(root/"nac_angle.xvg",2),start,end)
    gate=_window(_read_table(root/"gate_opening.xvg",2),start,end)
    thermo=_window(_read_table(root/"thermo.xvg",4),start,end)
    pocket=_window(_read_table(root/"pocket_state.xvg",3),start,end)
    for label,rows in (("angle",angle),("gate",gate),("thermo",thermo),("pocket",pocket)):
        _validate_times(distance,rows,label)

    contact=json.loads((root/"minimum_contact.json").read_text())
    minimum_protein=float(contact["minimum_ligand_protein_heavy_nm"])
    minimum_water=float(contact["minimum_ligand_water_heavy_nm"])
    if not all(math.isfinite(v) for v in (minimum_protein,minimum_water)):
        raise ValueError("minimum contact contains non-finite value")
    severe_clash=minimum_protein<0.18
    retained=[row[1]>=3.0 and row[2]<=1.2 for row in pocket]
    numerical=_scan_log(root/"run.log")
    numerical_pass=numerical.pop("finished_mdrun") and not any(numerical.values())
    thermo_audit=_thermodynamics(thermo)

    if not numerical_pass:
        scientific="NOT_EVALUATED_NUMERICAL_FAIL"
    elif severe_clash:
        scientific="FAIL_SEVERE_CLASH"
    elif not retained[-1]:
        scientific="FAIL_UNBOUND"
    elif not thermo_audit["stable"]:
        scientific="NOT_EVALUATED_THERMODYNAMIC_FAIL"
    elif nac["nac_frame_count"]==0:
        scientific="FAIL_REPLICA_NO_NAC"
    else:
        scientific="PASS_REPLICA_NAC_PRESENT"

    potential=None
    if (root/"potential_energy.xvg").is_file():
        rows=_window(_read_table(root/"potential_energy.xvg",2),start,end)
        _validate_times(distance,rows,"potential")
        potential=_stats(row[1] for row in rows)

    return {
        "schema_version":1,
        "candidate_id":manifest["candidate_id"],
        "velocity_seed":int(manifest["velocity_seed"]),
        "microstate":MICROSTATE,
        "gate_definition":"NylC residues 261-266; Thr267 excluded",
        "nac_contract":{
            "joint_required":True,
            "distance_max_nm":distance_max_nm,
            "angle_min_deg":angle_min_deg,
            "angle_max_deg":angle_max_deg,
        },
        "fully_unrestrained":True,
        "analysis":{
            "window_ps":[start,end],
            "frame_count":len(distance),
            "sampling_interval_ps":_median_spacing(distance),
            "actual_denominator":len(distance),
        },
        "technical_status":"PASS" if numerical_pass else "FAIL",
        "scientific_status":scientific,
        "nac":nac,
        "distance_nm":_stats(row[1] for row in distance),
        "angle_deg":_stats(row[1] for row in angle),
        "gate_opening_nm":_stats(row[1] for row in gate),
        "thermodynamics":thermo_audit,
        "potential_energy_kj_mol":potential,
        "bound_state":{
            "definition":"contacts_le_0.45nm>=3 and ligand_pocket_com<=1.2nm",
            "retained_frame_count":sum(retained),
            "retained_occupancy":sum(retained)/len(retained),
            "final_frame_retained":retained[-1],
            "minimum_ligand_protein_heavy_nm":minimum_protein,
            "minimum_ligand_water_heavy_nm":minimum_water,
            "severe_clash":severe_clash,
        },
        "numerical_issue_counts":numerical,
        "scientific_scope":"classical-MM A1 NAC stability only; no proton-path or QM/MM barrier inference",
    }

def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--run-root",type=pathlib.Path,required=True)
    ap.add_argument("--manifest",type=pathlib.Path,required=True)
    ap.add_argument("--window-start-ps",type=float,required=True)
    ap.add_argument("--window-end-ps",type=float,required=True)
    ap.add_argument("--output",type=pathlib.Path,required=True)
    args=ap.parse_args()
    result=audit_replica(args.run_root,json.loads(args.manifest.read_text()),(args.window_start_ps,args.window_end_ps))
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,sort_keys=True))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
