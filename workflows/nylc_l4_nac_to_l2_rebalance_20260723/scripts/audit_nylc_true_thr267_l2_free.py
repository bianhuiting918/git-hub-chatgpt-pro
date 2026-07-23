#!/usr/bin/env python3
"""Audit the final fully unrestrained NylC PA66-L2 NPT window."""
import argparse
import json
import pathlib
import re
import statistics
import sys

SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import audit_true_thr267_recapture_pilot as recapture
import extract_true_thr267_recapture_sources as geom

THR267_OG1 = 8961
L2_REACTIVE_C = 10287
L2_REACTIVE_O = 10288
GATE_AXIS = (0.4904295935403325, 0.813080787402625, -0.3136533866174433)
GATE_BASELINE = -1.51109


def read_xvg(path):
    rows=[]
    for line in pathlib.Path(path).read_text(errors="replace").splitlines():
        if line.strip() and line[0] not in "#@":
            rows.append([float(value) for value in line.split()])
    return rows


def end_geometry(path):
    atoms,box,_=geom.read_gro(path)
    thr=geom.validate_identity(atoms,THR267_OG1,(267,"THR","OG1"))
    carbon=atoms[L2_REACTIVE_C]
    oxygen=atoms[L2_REACTIVE_O]
    c_to_thr=geom.minimum_image(geom.subtract(thr["xyz_nm"],carbon["xyz_nm"]),box)
    c_to_o=geom.minimum_image(geom.subtract(oxygen["xyz_nm"],carbon["xyz_nm"]),box)
    return geom.norm(c_to_thr),geom.angle_deg(c_to_thr,c_to_o)


def summary(rows, labels):
    result={}
    for column,label in enumerate(labels,1):
        values=[row[column] for row in rows]
        result[label]={"mean":statistics.fmean(values),"stdev":statistics.pstdev(values),"min":min(values),"max":max(values)}
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--run-root",type=pathlib.Path,required=True)
    parser.add_argument("--job-id",required=True)
    args=parser.parse_args()
    root=args.run_root
    series=json.loads((root/"nac_energy_series.json").read_text())
    log=(root/"run.log").read_text(errors="replace")
    counts={
      "fatal":len(re.findall(r"(?i)fatal(?:\s+error)?",log)),
      "lincs_warning":len(re.findall(r"(?i)lincs\s+warning",log)),
      "settle_problem":len(re.findall(r"(?i)(?:cannot\s+be\s+settled|settle[^\n]*(?:warning|error|failed))",log)),
      "nan":len(re.findall(r"(?i)\bnan\b",log)),
    }
    finished=bool(re.search(r"Finished mdrun",log))
    atoms,box,_=geom.read_gro(root/"run.gro")
    contact=recapture.minimum_ligand_protein_distance(atoms,box)
    distance,angle=end_geometry(root/"run.gro")
    gate=[]
    for row in read_xvg(root/"gate_core_vector.xvg"):
        absolute=sum(row[i+1]*GATE_AXIS[i] for i in range(3))
        gate.append(absolute-GATE_BASELINE)
    numerical_pass=finished and not any(counts.values())
    contact_pass=contact["distance_nm"]>=0.10 and contact["minimum_heavy_atom_contact"]["distance_nm"]>=0.18
    if not numerical_pass:
        scientific="NOT_EVALUATED_NUMERICAL_FAIL"
    elif not contact_pass:
        scientific="FAIL_UNRESTRAINED_L2_SEVERE_CONTACT"
    elif series["nac_frame_count"]>0:
        scientific="PASS_UNRESTRAINED_L2_NAC_PRESENT"
    else:
        scientific="FAIL_UNRESTRAINED_L2_NO_NAC"
    result={
      "schema_version":1,
      "candidate_id":"nylc_C18_trueT267_freeGS",
      "slurm_job_id":args.job_id,
      "window_type":"fully_unrestrained_NPT_1ns",
      "technical_status":"PASS" if numerical_pass else "FAIL",
      "scientific_status":scientific,
      "nac":series,
      "end_geometry":{"distance_nm":round(distance,6),"angle_deg":round(angle,3),"joint_nac":bool(distance<=0.35 and 95.0<=angle<=115.0)},
      "minimum_ligand_protein_contact":contact,
      "gate_definition":"NylC residues 261-266; Thr267 excluded",
      "gate_opening_nm":{"frame_count":len(gate),"mean":statistics.fmean(gate),"stdev":statistics.pstdev(gate),"min":min(gate),"max":max(gate)},
      "thermodynamics":summary(read_xvg(root/"thermo.xvg"),("temperature_K","pressure_bar","volume_nm3")),
      "numerical_issue_counts":counts,
      "step1_qmmm_gate":"ELIGIBLE_FOR_DFTB3_PREFLIGHT" if scientific=="PASS_UNRESTRAINED_L2_NAC_PRESENT" else "NOT_ELIGIBLE_NO_FREE_L2_NAC",
    }
    (root/"l2_free_1ns_audit.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,sort_keys=True))
    raise SystemExit(0 if numerical_pass else 2)


if __name__=="__main__":
    main()
