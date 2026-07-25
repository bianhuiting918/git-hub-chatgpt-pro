#!/bin/bash
#SBATCH -J nylc_m1_audit20
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem-per-cpu=4000M
#SBATCH -t 03:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/nylc_m1_audit20_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/nylc_m1_audit20_%j.err

set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
PYMDA="$TASK/tools/proton_geometry/venv/bin/python"
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
M1_REBAL_JOB="${M1_REBAL_JOB:?}"
M1_EXT_JOB="${M1_EXT_JOB:?}"
candidate=nylc_C18_trueT267_freeGS
build="$TASK/candidates/$candidate/microstates/M1_nalphaH2_Asp306H/build_job_61801089/build"
base="$TASK/candidates/$candidate/microstates/M1_nalphaH2_Asp306H/runs/rebalance_job_$M1_REBAL_JOB/npt300free"
ext="$TASK/candidates/$candidate/microstates/M1_nalphaH2_Asp306H/runs/free20ns_job_$M1_EXT_JOB"
job="${SLURM_JOB_ID:?}"
out="$TASK/candidates/$candidate/microstates/M1_nalphaH2_Asp306H/analysis/full20ns_job_$job"
mkdir -p "$TASK/logs"
test ! -e "$out"
mkdir -p "$out"

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  ( flock -x 9; printf '%s\tnylc_m1_full20ns_audit\t%s\tcandidate=%s;job=%s;%s\n' "$stamp" "$state" "$candidate" "$job" "$detail" >&9 ) 9>>"$TASK/run_history.tsv"
  ( flock -x 8; "$PY" - "$stamp" "$state" "$detail" "$job" >&8 <<'PY'
import json,sys
stamp,state,detail,job=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":"nylc_m1_full20ns_audit","state":state,"detail":detail,"candidate_id":"nylc_C18_trueT267_freeGS","microstate":"M1_NalphaH2_Asp306H","slurm_job_id":job},sort_keys=True))
PY
  ) 8>>"$TASK/run_history.jsonl"
}
fail() { local rc=$?; trap - ERR; set +e; record FAIL_TECHNICAL "exit_code=$rc;out=$out"; exit "$rc"; }
trap fail ERR
record START "rebalance_job=$M1_REBAL_JOB;extension_job=$M1_EXT_JOB;window=0_to_20ns;fully_unrestrained=yes"

test -s "$base/run.tpr"
test -s "$base/run.xtc"
test -s "$base/run.edr"
test -s "$base/run.log"
mapfile -t ext_xtc < <(find "$ext" -maxdepth 1 -type f -name 'run*.xtc' -printf '%p\n' | sort)
mapfile -t ext_edr < <(find "$ext" -maxdepth 1 -type f -name 'run*.edr' -printf '%p\n' | sort)
mapfile -t ext_log < <(find "$ext" -maxdepth 1 -type f -name 'run*.log' -printf '%p\n' | sort)
mapfile -t ext_gro < <(find "$ext" -maxdepth 1 -type f -name 'run*.gro' -printf '%p\n' | sort)
test "${#ext_xtc[@]}" -eq 1
test "${#ext_edr[@]}" -eq 1
test "${#ext_log[@]}" -eq 1
test "${#ext_gro[@]}" -eq 1
test -s "$ext/run.tpr"

set +u
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
set -u
export GMX_MAXBACKUP=-1

printf '[ nac_angle ]\n8960 10287 10288\n' >"$out/nac_angle.ndx"
analyze_segment() {
  local tag=$1 tpr=$2 xtc=$3 edr=$4
  "$GMX" distance -f "$xtc" -s "$tpr" -select 'atomnr 8960 plus atomnr 10287' -oall "$out/${tag}_distance.xvg" >"$out/${tag}_distance.stdout" 2>"$out/${tag}_distance.stderr"
  printf '0\n' | "$GMX" angle -f "$xtc" -n "$out/nac_angle.ndx" -ov "$out/${tag}_angle.xvg" -type angle >"$out/${tag}_angle.stdout" 2>"$out/${tag}_angle.stderr"
  printf 'Potential\n0\n' | "$GMX" energy -f "$edr" -o "$out/${tag}_potential.xvg" >"$out/${tag}_potential.stdout" 2>"$out/${tag}_potential.stderr"
  printf 'Temperature\nPressure\nVolume\n0\n' | "$GMX" energy -f "$edr" -o "$out/${tag}_thermo.xvg" >"$out/${tag}_thermo.stdout" 2>"$out/${tag}_thermo.stderr"
  "$GMX" distance -f "$xtc" -s "$tpr" -n "$build/source_cycle.ndx" -select 'com of group "Core" plus com of group "Gate"' -oxyz "$out/${tag}_gate.xvg" >"$out/${tag}_gate.stdout" 2>"$out/${tag}_gate.stderr"
}
analyze_segment base "$base/run.tpr" "$base/run.xtc" "$base/run.edr"
analyze_segment ext "$ext/run.tpr" "${ext_xtc[0]}" "${ext_edr[0]}"

"$PY" - "$out" <<'PY'
import pathlib,sys
root=pathlib.Path(sys.argv[1])
def read(path):
    rows=[]
    for line in path.read_text(errors="replace").splitlines():
        s=line.strip()
        if s and s[0] not in "#@":
            rows.append([float(x) for x in s.split()])
    return rows
def merge(name):
    by_time={}
    for tag in ("base","ext"):
        for row in read(root/f"{tag}_{name}.xvg"):
            by_time[round(row[0],6)]=row
    rows=[by_time[k] for k in sorted(by_time)]
    (root/f"combined_{name}.xvg").write_text("\n".join(" ".join(f"{x:.10g}" for x in row) for row in rows)+"\n")
for name in ("distance","angle","potential","thermo","gate"):
    merge(name)
PY

"$PY" "$FLOW/scripts/analyze_nac_series.py" --distance-xvg "$out/combined_distance.xvg" --angle-xvg "$out/combined_angle.xvg" --potential-xvg "$out/combined_potential.xvg" --distance-max-nm 0.35 --angle-min-deg 95 --angle-max-deg 115 --output "$out/nac_energy_series.json"
"$PY" - "$out/case.json" "$ext/run.tpr" "$base/run.xtc" "${ext_xtc[0]}" "$M1_REBAL_JOB" "$M1_EXT_JOB" <<'PY'
import json,sys
path,tpr,base_xtc,ext_xtc,reb,extjob=sys.argv[1:]
cfg={"schema_version":1,"analysis_scope":"fully_unrestrained_explicit_water_fixed_topology_M1_full_0_to_20ns","thresholds_angstrom_degree":{"favorable":{"donor_acceptor_max":3.5,"hydrogen_acceptor_max":2.5,"dha_min_deg":135.0},"strong":{"donor_acceptor_max":3.0,"hydrogen_acceptor_max":2.2,"dha_min_deg":150.0},"nac":{"attack_distance_max":3.5,"attack_angle_min_deg":95.0,"attack_angle_max_deg":115.0}},"cases":[{"id":"nylc_C18_M1_NalphaH2_Asp306H_full20ns","enzyme":"NylC-GYAQ","substrate":"PA66-L2","window":"fully_unrestrained_NPT_0_to_20ns","scientific_role":"fixed_topology_MM_preorganization_not_proton_transfer","source_job_ids":[reb,extjob],"atom_mapping_authority":"M1_build_job_61801089_PASS.json","thr267_og1_index1":8960,"carbonyl_c_index1":10287,"carbonyl_o_index1":10288,"tpr":tpr,"xtc":[base_xtc,ext_xtc]}]}
open(path,"w").write(json.dumps(cfg,indent=2,sort_keys=True)+"\n")
PY
"$PYMDA" "$FLOW/scripts/analyze_nylc_m1_proton_geometry.py" --cases "$out/case.json" --outdir "$out/geometry" >"$out/geometry.stdout" 2>"$out/geometry.stderr"

"$PY" - "$FLOW/scripts" "$out/geometry/combined_summary.json" "$out/nac_energy_series.json" "$out/combined_thermo.xvg" "$out/combined_gate.xvg" "$base/run.log" "${ext_log[0]}" "${ext_gro[0]}" "$out/AUDIT.json" "$job" "$M1_REBAL_JOB" "$M1_EXT_JOB" <<'PY'
import json,math,pathlib,re,statistics,sys
scripts,geom_path,nac_path,thermo_path,gate_path,base_log,ext_log,final_gro,out_path,job,reb,extjob=sys.argv[1:]
sys.path.insert(0,scripts)
import audit_true_thr267_recapture_pilot as recapture
import extract_true_thr267_recapture_sources as grogeom
geom=json.load(open(geom_path))
assert geom["technical_status"]=="PASS" and not geom["failures"] and len(geom["results"])==1
r=geom["results"][0]
assert r["topology_state"]["microstate_label"]=="M1_NalphaH2_Asp306H"
nac=json.load(open(nac_path))
def rows(path):
    out=[]
    for line in open(path,errors="replace"):
        s=line.strip()
        if s and s[0] not in "#@": out.append([float(x) for x in s.split()])
    return out
thermo_rows=rows(thermo_path)
labels=("temperature_K","pressure_bar","volume_nm3")
thermo={}
for i,label in enumerate(labels,1):
    v=[x[i] for x in thermo_rows]
    thermo[label]={"mean":statistics.fmean(v),"stdev":statistics.pstdev(v),"min":min(v),"max":max(v)}
vol=[x[3] for x in thermo_rows]; half=max(1,len(vol)//2)
volume_half_drift=abs(statistics.fmean(vol[:half])-statistics.fmean(vol[-half:]))/statistics.fmean(vol)
thermo_pass=295.0<=thermo["temperature_K"]["mean"]<=305.0 and abs(thermo["pressure_bar"]["mean"])<=100.0 and volume_half_drift<=0.02
axis=(0.4904295935403325,0.813080787402625,-0.3136533866174433); baseline=-1.51109
gate=[sum(row[i+1]*axis[i] for i in range(3))-baseline for row in rows(gate_path)]
atoms,box,_=grogeom.read_gro(pathlib.Path(final_gro))
contact=recapture.minimum_ligand_protein_distance(atoms,box,ligand_resname="L2")
contact_pass=contact["distance_nm"]>=0.10 and contact["minimum_heavy_atom_contact"]["distance_nm"]>=0.18
text="\n".join(open(p,errors="replace").read() for p in (base_log,ext_log))
issues={"fatal":len(re.findall(r"(?i)fatal(?:\s+error)?",text)),"lincs_warning":len(re.findall(r"(?i)lincs\s+warning",text)),"settle_problem":len(re.findall(r"(?i)(?:cannot\s+be\s+settled|settle[^\n]*(?:warning|error|failed))",text)),"nan":len(re.findall(r"(?i)\bnan\b",text)),"finished_mdrun_count":len(re.findall(r"Finished mdrun",text))}
numerical_pass=issues["finished_mdrun_count"]==2 and not any(issues[k] for k in ("fatal","lincs_warning","settle_problem","nan"))
trajectory_pass=r["trajectory"]["frame_count"]>=9990 and r["trajectory"]["time_start_ps"]<=2.0 and r["trajectory"]["time_end_ps"]>=19998.0
route=r["route_statistics"]["OgH_to_Nalpha_preorg"]
nalpha_nac_favorable=route["nac_frames"]["favorable_frames"]
nac_present=nac["nac_frame_count"]>0
technical_pass=numerical_pass and thermo_pass and contact_pass and trajectory_pass
postrelay_eligible=technical_pass and nac_present and nalpha_nac_favorable>0
if not numerical_pass:
    scientific="NOT_EVALUATED_NUMERICAL_FAIL"
elif not thermo_pass:
    scientific="NOT_EVALUATED_THERMODYNAMIC_INSTABILITY"
elif not contact_pass:
    scientific="FAIL_M1_20NS_SEVERE_CONTACT"
elif not trajectory_pass:
    scientific="NOT_EVALUATED_INCOMPLETE_20NS"
elif not nac_present:
    scientific="FAIL_M1_20NS_NO_NAC"
elif nalpha_nac_favorable==0:
    scientific="FAIL_M1_NAC_WITHOUT_NALPHA_PROTON_ACCEPTOR_PREORGANIZATION"
else:
    scientific="PASS_M1_NAC_AND_NALPHA_PROTON_ACCEPTOR_PREORGANIZATION"
payload={"schema_version":1,"status":"PASS_TECHNICAL_FULL20NS_AUDIT" if technical_pass else "NOT_EVALUATED_TECHNICAL_GATE","scientific_status":scientific,"qmmm_gate":"ELIGIBLE_POST_RELAY_QMMM_CONSTRUCTION" if postrelay_eligible else "NOT_ELIGIBLE_POST_RELAY_QMMM","candidate_id":"nylc_C18_trueT267_freeGS","microstate":"M1_NalphaH2_Asp306H","slurm_job_id":job,"source_job_ids":{"rebalance":reb,"extension":extjob},"trajectory":r["trajectory"],"nac":nac,"proton_geometry":r["route_statistics"],"topology_state":r["topology_state"],"thermodynamics":thermo,"volume_half_drift_fraction":volume_half_drift,"thermodynamic_gate_pass":thermo_pass,"gate_definition":"NylC residues 261-266; Thr267 excluded","gate_opening_nm":{"frame_count":len(gate),"mean":statistics.fmean(gate),"stdev":statistics.pstdev(gate),"min":min(gate),"max":max(gate)},"minimum_ligand_protein_contact":contact,"contact_gate_pass":contact_pass,"numerical_issue_counts":issues,"trajectory_gate_pass":trajectory_pass,"postrelay_gate_definition":"technical+thermo+contact+complete20ns, at least one NAC frame, and at least one favorable Thr-OH to Nalpha geometry within NAC frames","fixed_topology_warning":"Classical MM geometry can nominate preorganization but cannot execute proton transfer or establish a barrier."}
open(out_path,"w").write(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(json.dumps({"status":payload["status"],"scientific_status":scientific,"qmmm_gate":payload["qmmm_gate"],"nac_frames":nac["nac_frame_count"],"nalpha_nac_favorable":nalpha_nac_favorable},sort_keys=True))
PY

technical="$("$PY" - "$out/AUDIT.json" <<'PY'
import json,sys
d=json.load(open(sys.argv[1])); print(d["status"])
PY
)"
scientific="$("$PY" - "$out/AUDIT.json" <<'PY'
import json,sys
d=json.load(open(sys.argv[1])); print(d["scientific_status"])
PY
)"
qmgate="$("$PY" - "$out/AUDIT.json" <<'PY'
import json,sys
d=json.load(open(sys.argv[1])); print(d["qmmm_gate"])
PY
)"
if [[ "$technical" == PASS_TECHNICAL_FULL20NS_AUDIT ]]; then
  cp "$out/AUDIT.json" "$out/PASS.json"
else
  cp "$out/AUDIT.json" "$out/NOT_EVALUATED.json"
fi

nac_count="$("$PY" - "$out/nac_energy_series.json" <<'PY'
import json,sys
print(json.load(open(sys.argv[1]))["nac_frame_count"])
PY
)"
if [[ "$nac_count" -gt 0 ]]; then
  selected_time="$("$PY" - "$out/nac_energy_series.json" <<'PY'
import json,sys
print(json.load(open(sys.argv[1]))["lowest_potential_nac_frame"]["time_ps"])
PY
)"
  if awk -v t="$selected_time" 'BEGIN{exit !(t<=1000.0)}'; then
    source_tpr="$base/run.tpr"; source_xtc="$base/run.xtc"
  else
    source_tpr="$ext/run.tpr"; source_xtc="${ext_xtc[0]}"
  fi
  printf 'System\n' | "$GMX" trjconv -f "$source_xtc" -s "$source_tpr" -dump "$selected_time" -o "$out/lowest_potential_NAC.gro" >"$out/select.stdout" 2>"$out/select.stderr"
  sha256sum "$out/lowest_potential_NAC.gro" >"$out/lowest_potential_NAC.sha256"
fi
record COMPLETE "technical=$technical;scientific=$scientific;qmmm_gate=$qmgate;out=$out/AUDIT.json"
trap - ERR
