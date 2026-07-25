#!/bin/bash
#SBATCH -J nylc_m1_rebal
#SBATCH -p xahdnormal
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -c 8
#SBATCH --gres=dcu:2
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 06:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/nylc_m1_rebalance_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/nylc_m1_rebalance_%j.err

set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
candidate=nylc_C18_trueT267_freeGS
M1_EM_JOB="${M1_EM_JOB:?submit with --export=ALL,M1_EM_JOB=<passed EM job>}"
build="$TASK/candidates/$candidate/microstates/M1_nalphaH2_Asp306H/build_job_61801089/build"
em="$TASK/candidates/$candidate/microstates/M1_nalphaH2_Asp306H/runs/em_job_$M1_EM_JOB"
em_pass="$em/PASS.json"
job="${SLURM_JOB_ID:?}"
runs="$TASK/candidates/$candidate/microstates/M1_nalphaH2_Asp306H/runs/rebalance_job_$job"
mkdir -p "$TASK/logs"
test ! -e "$runs"
mkdir -p "$runs"
set +u
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
set -u
export GMX_MAXBACKUP=-1

record() {
  local event=$1 state=$2 detail=$3 stamp
  stamp="$(date -Is)"
  ( flock -x 9; printf '%s\t%s\t%s\tcandidate=%s;job=%s;%s\n' "$stamp" "$event" "$state" "$candidate" "$job" "$detail" >&9 ) 9>>"$TASK/run_history.tsv"
  ( flock -x 8; "$PY" - "$stamp" "$event" "$state" "$detail" "$job" >&8 <<'PY'
import json,sys
stamp,event,state,detail,job=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":event,"state":state,"detail":detail,"candidate_id":"nylc_C18_trueT267_freeGS","microstate":"M1_NalphaH2_Asp306H","slurm_job_id":job},sort_keys=True))
PY
  ) 8>>"$TASK/run_history.jsonl"
}
scan_log() {
  "$PY" - "$1" "$2" <<'PY'
import json,re,sys
text=open(sys.argv[1],errors="replace").read()
out={"fatal":len(re.findall(r"(?i)fatal(?:\s+error)?",text)),"lincs_warning":len(re.findall(r"(?i)lincs\s+warning",text)),"settle_problem":len(re.findall(r"(?i)(?:cannot\s+be\s+settled|settle[^\n]*(?:warning|error|failed))",text)),"nan":len(re.findall(r"(?i)\bnan\b",text)),"finished_mdrun":bool(re.search(r"Finished mdrun",text))}
open(sys.argv[2],"w").write(json.dumps(out,indent=2,sort_keys=True)+"\n")
raise SystemExit(0 if out["finished_mdrun"] and not any(out[k] for k in ("fatal","lincs_warning","settle_problem","nan")) else 1)
PY
}
run_stage() {
  local name=$1 coord=$2 checkpoint=${3:-} restrained=${4:-yes}
  local stage="$runs/$name" mdp="$FLOW/mdp/${name}_m1.mdp"
  mkdir -p "$stage"
  local args=(grompp -f "$mdp" -c "$coord" -p "$build/topol.top" -o "$stage/run.tpr" -po "$stage/run.processed.mdp" -maxwarn 0)
  if [[ "$restrained" == yes ]]; then args+=(-r "$build/system_M1.gro"); fi
  if [[ -n "$checkpoint" ]]; then args+=(-t "$checkpoint"); fi
  record "$name" START "restrained=$restrained;parent_em=$M1_EM_JOB;microstate=M1_NalphaH2_Asp306H"
  (cd "$build" && "$GMX" "${args[@]}" >"$stage/grompp.stdout" 2>"$stage/grompp.stderr")
  mpirun -np 2 "$GMX" mdrun -s "$stage/run.tpr" -deffnm "$stage/run" -ntomp "${SLURM_CPUS_PER_TASK:-8}" -nb gpu -pme gpu -npme 1 -bonded gpu -pin on >"$stage/mdrun.stdout" 2>"$stage/mdrun.stderr"
  test -s "$stage/run.gro"
  scan_log "$stage/run.log" "$stage/numerical_audit.json"
  "$PY" - "$stage/PASS.json" "$name" "$restrained" "$M1_EM_JOB" <<'PY'
import json,sys
path,stage,restrained,parent=sys.argv[1:]
open(path,"w").write(json.dumps({"status":"PASS_TECHNICAL_STAGE","candidate_id":"nylc_C18_trueT267_freeGS","microstate":"M1_NalphaH2_Asp306H","stage":stage,"restrained":restrained=="yes","parent_em_job_id":parent,"scientific_nac_evidence":stage=="npt300free" and restrained=="no"},indent=2,sort_keys=True)+"\n")
PY
  record "$name" PASS_TECHNICAL "restrained=$restrained;parent_em=$M1_EM_JOB"
}

failure=1
current=initialization
trap 'rc=$?; if [[ $failure -ne 0 ]]; then record nylc_m1_rebalance FAIL_STAGE "stage=$current;rc=$rc"; fi; exit $rc' EXIT
"$PY" - "$em_pass" "$M1_EM_JOB" <<'PY'
import json,sys
d=json.load(open(sys.argv[1]))
assert d["status"]=="PASS_TECHNICAL_EM"
assert d["slurm_job_id"]==sys.argv[2]
assert d["microstate"]=="M1_NalphaH2_Asp306H"
assert d["fmax_kj_mol_nm"]<=d["threshold_kj_mol_nm"]
PY
test -s "$em/run.gro"
record nylc_m1_rebalance START "parent_em=$M1_EM_JOB;stages=heavy_restrained_lowT_then_gradual_release_then_1ns_free"

current=nvt50
run_stage nvt50 "$em/run.gro" "" yes
current=nvt150
run_stage nvt150 "$runs/nvt50/run.gro" "$runs/nvt50/run.cpt" yes
current=nvt300
run_stage nvt300 "$runs/nvt150/run.gro" "$runs/nvt150/run.cpt" yes
current=npt300r
run_stage npt300r "$runs/nvt300/run.gro" "$runs/nvt300/run.cpt" yes
current=npt300rel
run_stage npt300rel "$runs/npt300r/run.gro" "$runs/npt300r/run.cpt" yes
current=npt300free
run_stage npt300free "$runs/npt300rel/run.gro" "$runs/npt300rel/run.cpt" no

current=scientific_audit
free="$runs/npt300free"
"$GMX" distance -f "$free/run.xtc" -s "$free/run.tpr" -select 'atomnr 8960 plus atomnr 10287' -oall "$free/nac_distance.xvg" >"$free/nac_distance.stdout" 2>"$free/nac_distance.stderr"
printf '[ nac_angle ]\n8960 10287 10288\n' >"$free/nac_angle.ndx"
printf '0\n' | "$GMX" angle -f "$free/run.xtc" -n "$free/nac_angle.ndx" -ov "$free/nac_angle.xvg" -type angle >"$free/nac_angle.stdout" 2>"$free/nac_angle.stderr"
printf 'Potential\n0\n' | "$GMX" energy -f "$free/run.edr" -o "$free/potential_energy.xvg" >"$free/potential.stdout" 2>"$free/potential.stderr"
printf 'Temperature\nPressure\nVolume\n0\n' | "$GMX" energy -f "$free/run.edr" -o "$free/thermo.xvg" >"$free/thermo.stdout" 2>"$free/thermo.stderr"
"$GMX" distance -f "$free/run.xtc" -s "$free/run.tpr" -n "$build/source_cycle.ndx" -select 'com of group "Core" plus com of group "Gate"' -oxyz "$free/gate_core_vector.xvg" >"$free/gate.stdout" 2>"$free/gate.stderr"
"$PY" "$FLOW/scripts/analyze_nac_series.py" --distance-xvg "$free/nac_distance.xvg" --angle-xvg "$free/nac_angle.xvg" --potential-xvg "$free/potential_energy.xvg" --distance-max-nm 0.35 --angle-min-deg 95 --angle-max-deg 115 --output "$free/nac_energy_series.json"
"$PY" "$FLOW/scripts/audit_nylc_m1_free.py" --run-root "$free" --job-id "$job" >"$free/audit.stdout" 2>"$free/audit.stderr"
scientific="$("$PY" - "$free/m1_free_1ns_audit.json" <<'PY'
import json,sys
print(json.load(open(sys.argv[1]))["scientific_status"])
PY
)"
record nylc_m1_rebalance TECHNICAL_COMPLETE "scientific_status=$scientific;audit=$free/m1_free_1ns_audit.json;gate=residues_261_266_excluding_Thr267"
failure=0
trap - EXIT
