#!/bin/bash
#SBATCH -J nylc_m1_p400
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 01:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/nylc_m1_p400_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/nylc_m1_p400_%j.err
set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
PYMDA="$TASK/tools/proton_geometry/venv/bin/python"
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
source_job=61801789
run="$TASK/candidates/nylc_C18_trueT267_freeGS/microstates/M1_nalphaH2_Asp306H/runs/rebalance_job_$source_job/npt300free"
job="${SLURM_JOB_ID:?}"
out="$TASK/candidates/nylc_C18_trueT267_freeGS/microstates/M1_nalphaH2_Asp306H/analysis/partial_0_400ps_job_$job"
mkdir -p "$TASK/logs"
test ! -e "$out"; mkdir -p "$out"
record() {
 local state=$1 detail=$2 stamp; stamp="$(date -Is)"
 ( flock -x 9; printf '%s\tnylc_m1_partial400_audit\t%s\tcandidate=nylc_C18_trueT267_freeGS;job=%s;%s\n' "$stamp" "$state" "$job" "$detail" >&9 ) 9>>"$TASK/run_history.tsv"
 ( flock -x 8; "$PY" - "$stamp" "$state" "$detail" "$job" >&8 <<'PY'
import json,sys
stamp,state,detail,job=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":"nylc_m1_partial400_audit","state":state,"detail":detail,"candidate_id":"nylc_C18_trueT267_freeGS","microstate":"M1_NalphaH2_Asp306H","slurm_job_id":job},sort_keys=True))
PY
 ) 8>>"$TASK/run_history.jsonl"
}
fail(){ local rc=$?; trap - ERR; set +e; record FAIL_TECHNICAL "rc=$rc;out=$out"; exit "$rc"; }
trap fail ERR
record START "source_job=$source_job;window=0_400ps;diagnostic_only=yes"
test -s "$run/run.tpr"; test -s "$run/run.xtc"; test -s "$run/run.edr"
set +u; source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh; set -u
export GMX_MAXBACKUP=-1
printf 'System\n' | "$GMX" trjconv -f "$run/run.xtc" -s "$run/run.tpr" -b 0 -e 400 -o "$out/snapshot_0_400ps.xtc" >"$out/trjconv.stdout" 2>"$out/trjconv.stderr"
"$GMX" distance -f "$out/snapshot_0_400ps.xtc" -s "$run/run.tpr" -select 'atomnr 8960 plus atomnr 10287' -oall "$out/distance.xvg" >"$out/distance.stdout" 2>"$out/distance.stderr"
printf '[ nac_angle ]\n8960 10287 10288\n' >"$out/nac_angle.ndx"
printf '0\n' | "$GMX" angle -f "$out/snapshot_0_400ps.xtc" -n "$out/nac_angle.ndx" -ov "$out/angle.xvg" -type angle >"$out/angle.stdout" 2>"$out/angle.stderr"
printf 'Potential\n0\n' | "$GMX" energy -f "$run/run.edr" -b 0 -e 400 -o "$out/potential.xvg" >"$out/potential.stdout" 2>"$out/potential.stderr"
"$PY" "$FLOW/scripts/analyze_nac_series.py" --distance-xvg "$out/distance.xvg" --angle-xvg "$out/angle.xvg" --potential-xvg "$out/potential.xvg" --distance-max-nm 0.35 --angle-min-deg 95 --angle-max-deg 115 --output "$out/nac.json"
"$PY" - "$out/case.json" "$run/run.tpr" "$out/snapshot_0_400ps.xtc" <<'PY'
import json,sys
path,tpr,xtc=sys.argv[1:]
cfg={"schema_version":1,"analysis_scope":"partial_0_400ps_diagnostic_only","thresholds_angstrom_degree":{"favorable":{"donor_acceptor_max":3.5,"hydrogen_acceptor_max":2.5,"dha_min_deg":135.0},"strong":{"donor_acceptor_max":3.0,"hydrogen_acceptor_max":2.2,"dha_min_deg":150.0},"nac":{"attack_distance_max":3.5,"attack_angle_min_deg":95.0,"attack_angle_max_deg":115.0}},"cases":[{"id":"nylc_C18_M1_partial400","enzyme":"NylC-GYAQ","substrate":"PA66-L2","window":"fully_unrestrained_0_400ps_partial","scientific_role":"diagnostic_only_not_final_gate","source_job_id":"61801789","atom_mapping_authority":"M1_build_job_61801089_PASS.json","thr267_og1_index1":8960,"carbonyl_c_index1":10287,"carbonyl_o_index1":10288,"tpr":tpr,"xtc":xtc}]}
open(path,"w").write(json.dumps(cfg,indent=2,sort_keys=True)+"\n")
PY
"$PYMDA" "$FLOW/scripts/analyze_nylc_m1_proton_geometry.py" --cases "$out/case.json" --outdir "$out/geometry" >"$out/geometry.stdout" 2>"$out/geometry.stderr"
"$PY" - "$out/nac.json" "$out/geometry/combined_summary.json" "$out/PARTIAL_DIAGNOSTIC.json" "$job" <<'PY'
import json,sys
nac=json.load(open(sys.argv[1])); r=json.load(open(sys.argv[2]))["results"][0]
assert r["topology_state"]["microstate_label"]=="M1_NalphaH2_Asp306H"
routes={}
for name in ("OgH_to_Nalpha_preorg","OgH_to_water","OgH_to_Asp306","OgH_to_Asp308"):
    routes[name]={"all_frames":r["route_statistics"][name]["all_frames"],"nac_frames":r["route_statistics"][name]["nac_frames"]}
p={"schema_version":1,"status":"PASS_TECHNICAL_PARTIAL_DIAGNOSTIC","scientific_status":"NOT_EVALUATED_PARTIAL_400PS","qmmm_gate":"NOT_EVALUATED_REQUIRES_COMPLETE_20NS","candidate_id":"nylc_C18_trueT267_freeGS","microstate":"M1_NalphaH2_Asp306H","slurm_job_id":sys.argv[4],"source_job_id":"61801789","trajectory":r["trajectory"],"nac":nac,"proton_geometry":routes,"fixed_topology_warning":"Geometry only; no proton transfer or barrier inference."}
open(sys.argv[3],"w").write(json.dumps(p,indent=2,sort_keys=True)+"\n")
print(json.dumps({"status":p["status"],"nac_frames":nac["nac_frame_count"],"nac_occupancy":nac["nac_occupancy"],"nalpha_all_favorable":routes["OgH_to_Nalpha_preorg"]["all_frames"]["favorable_frames"],"nalpha_nac_favorable":routes["OgH_to_Nalpha_preorg"]["nac_frames"]["favorable_frames"]},sort_keys=True))
PY
record COMPLETE "status=PASS_TECHNICAL_PARTIAL_DIAGNOSTIC;scientific=NOT_EVALUATED_PARTIAL_400PS;out=$out/PARTIAL_DIAGNOSTIC.json"
trap - ERR
