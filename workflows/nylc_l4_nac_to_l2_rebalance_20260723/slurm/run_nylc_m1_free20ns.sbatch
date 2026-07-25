#!/bin/bash
#SBATCH -J nylc_m1_free20
#SBATCH -p xahdnormal
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -c 8
#SBATCH --gres=dcu:2
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 12:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/nylc_m1_free20_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/nylc_m1_free20_%j.err
set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
M1_REBAL_JOB="${M1_REBAL_JOB:?submit with --export=ALL,M1_REBAL_JOB=<rebalance job>}"
BASE="$TASK/candidates/nylc_C18_trueT267_freeGS/microstates/M1_nalphaH2_Asp306H/runs/rebalance_job_$M1_REBAL_JOB/npt300free"
ANALYZER="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/analyze_nylc_m1_proton_geometry.py"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
PYMDA="$TASK/tools/proton_geometry/venv/bin/python"
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
job="${SLURM_JOB_ID:?}"
OUT="$TASK/candidates/nylc_C18_trueT267_freeGS/microstates/M1_nalphaH2_Asp306H/runs/free20ns_job_${job}"
mkdir -p "$TASK/logs"
record() {
 local state=$1 detail=$2 stamp
 stamp="$(date -Is)"
 ( flock -x 9; printf '%s\tnylc_m1_free20ns\t%s\tcandidate=nylc_C18_trueT267_freeGS;job=%s;%s\n' "$stamp" "$state" "$job" "$detail" >&9 ) 9>>"$TASK/run_history.tsv"
 ( flock -x 8; "$PY" - "$stamp" "$state" "$detail" "$job" >&8 <<'PY'
import json,sys
stamp,state,detail,job=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":"nylc_m1_free20ns","state":state,"detail":detail,"candidate_id":"nylc_C18_trueT267_freeGS","slurm_job_id":job},sort_keys=True))
PY
 ) 8>>"$TASK/run_history.jsonl"
}
fail() { local rc=$?; trap - ERR; set +e; record FAIL_TECHNICAL "exit_code=$rc;out=$OUT"; exit "$rc"; }
trap fail ERR
record START "source_rebalance_job=$M1_REBAL_JOB;extension_ps=19000;scope=fully_unrestrained_fixed_topology_MM_geometry_screen_only"
"$PY" - "$BASE/PASS.json" <<'PY'
import json,sys
raw=open(sys.argv[1]).read()
d,end=json.JSONDecoder().raw_decode(raw)
assert raw[end:].strip() in ("", r"\n")
assert d["status"]=="PASS_TECHNICAL_STAGE" and d["stage"]=="npt300free" and d["restrained"] is False
assert d["microstate"]=="M1_NalphaH2_Asp306H"
PY
test ! -e "$OUT"
mkdir "$OUT"
set +u
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
set -u
export GMX_MAXBACKUP=-1
"$GMX" convert-tpr -s "$BASE/run.tpr" -extend 19000 -o "$OUT/run.tpr" >"$OUT/convert_tpr.stdout" 2>"$OUT/convert_tpr.stderr"
mpirun -np 2 "$GMX" mdrun -s "$OUT/run.tpr" -cpi "$BASE/run.cpt" -deffnm "$OUT/run" -noappend -ntomp "${SLURM_CPUS_PER_TASK:-8}" -nb gpu -pme gpu -npme 1 -bonded gpu -pin on >"$OUT/mdrun.stdout" 2>"$OUT/mdrun.stderr"
mapfile -t logs < <(find "$OUT" -maxdepth 1 -type f -name 'run*.log' -printf '%p\n' | sort)
mapfile -t xtcs < <(find "$OUT" -maxdepth 1 -type f -name 'run*.xtc' -printf '%p\n' | sort)
test "${#logs[@]}" -eq 1
test "${#xtcs[@]}" -eq 1
"$PY" - "${logs[0]}" "$OUT/mdrun.stderr" "$OUT/numerical_audit.json" <<'PY'
import json,re,sys
text="\n".join(open(p,errors="replace").read() for p in sys.argv[1:3])
out={"fatal":len(re.findall(r"(?i)fatal(?:\s+error)?",text)),"lincs_warning":len(re.findall(r"(?i)lincs\s+warning",text)),"settle_problem":len(re.findall(r"(?i)(?:cannot\s+be\s+settled|settle[^\n]*(?:warning|error|failed))",text)),"nan":len(re.findall(r"(?i)\bnan\b",text)),"finished_mdrun":bool(re.search(r"Finished mdrun",text))}
open(sys.argv[3],"w").write(json.dumps(out,indent=2,sort_keys=True)+"\n")
assert out["finished_mdrun"] and not any(out[k] for k in ("fatal","lincs_warning","settle_problem","nan"))
PY
"$PY" - "$OUT/case.json" "$OUT/run.tpr" "${xtcs[0]}" "$job" <<'PY'
import json,sys
out,tpr,xtc,job=sys.argv[1:]
cfg={"schema_version":1,"analysis_scope":"fully_unrestrained_explicit_water_fixed_topology_MM_geometry_only","excluded_superseded":[],"thresholds_angstrom_degree":{"favorable":{"donor_acceptor_max":3.5,"hydrogen_acceptor_max":2.5,"dha_min_deg":135.0},"strong":{"donor_acceptor_max":3.0,"hydrogen_acceptor_max":2.2,"dha_min_deg":150.0},"nac":{"attack_distance_max":3.5,"attack_angle_min_deg":95.0,"attack_angle_max_deg":115.0}},"cases":[{"id":"nylc_C18_M1_NalphaH2_Asp306H_free20ns_job_"+job,"enzyme":"NylC-GYAQ","substrate":"PA66-L2","window":"fully_unrestrained_NPT_extension_1_to_20ns","scientific_role":"M1_fixed_topology_MM_geometry_extension_not_proton_transfer","source_job_id":job,"atom_mapping_authority":"M1_build_job_61801089_PASS.json","thr267_og1_index1":8960,"carbonyl_c_index1":10287,"carbonyl_o_index1":10288,"tpr":tpr,"xtc":xtc}]}
open(out,"w").write(json.dumps(cfg,indent=2,sort_keys=True)+"\n")
PY
"$PYMDA" "$ANALYZER" --cases "$OUT/case.json" --outdir "$OUT/analysis" >"$OUT/analysis.stdout" 2>"$OUT/analysis.stderr"
"$PY" - "$OUT/analysis/combined_summary.json" "$OUT/PASS.json" "$job" "$M1_REBAL_JOB" <<'PY'
import json,sys
src,out,job,parent=sys.argv[1:]
d=json.load(open(src)); assert d["technical_status"]=="PASS" and not d["failures"] and len(d["results"])==1
r=d["results"][0]; t=r["trajectory"]
assert t["frame_count"]>=9000 and t["time_start_ps"]>=1000.0 and t["time_end_ps"]>=19998.0
assert r["topology_state"]["microstate_label"]=="M1_NalphaH2_Asp306H"
routes={}
for name in ("OgH_to_Asp306","OgH_to_Asp308","OgH_to_water","OgH_to_Nalpha_preorg"):
 routes[name]={"all_frames":r["route_statistics"][name]["all_frames"]["favorable_frames"],"nac_frames":r["route_statistics"][name]["nac_frames"]["favorable_frames"]}
payload={"status":"PASS_TECHNICAL_M1_20NS_GEOMETRY","scientific_status":"M1_GEOMETRY_SCREEN_COMPLETE_NOT_PROTON_TRANSFER","candidate_id":"nylc_C18_trueT267_freeGS","slurm_job_id":job,"source_rebalance_job_id":parent,"extension_ps":19000,"fully_unrestrained":True,"fixed_topology_microstate":"M1_NalphaH2_Asp306H","microstate_parameterization":"classical_neutral_Nalpha_proxy_plus_ASH306","trajectory":t,"nac":r["nac"],"favorable_route_counts":routes,"interpretation":"Distances and angles rank preorganization only; classical fixed-topology MD cannot transfer protons or establish a barrier."}
open(out,"w").write(json.dumps(payload,indent=2,sort_keys=True)+"\n")
PY
record PASS_TECHNICAL "status=PASS_TECHNICAL_M1_20NS_GEOMETRY;scientific=M1_GEOMETRY_SCREEN_COMPLETE_NOT_PROTON_TRANSFER;out=$OUT/PASS.json"
trap - ERR