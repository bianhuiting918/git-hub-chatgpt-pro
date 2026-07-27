#!/usr/bin/env bash
set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
GMX=/public/software/apps/gromacs/2022.2/hpcx-gcc7.3.1/bin/gmx_mpi
PARAM="$TASK_ROOT/a1_activated_nac_20260726/parameterization/job_61902961"
INDEX="${1:?usage: run_nylc_a1_full_system_preflight.sh ARRAY_INDEX}"
GITHUB_COMMIT="${A1_GITHUB_COMMIT:-unknown}"

CANDIDATES=(nac_evt18_time1206ps nac_evt25_time1462ps nac_evt08_time1086ps)
SOURCES=(
"$TASK_ROOT/ensemble/final_audit_job_61841413/assembled/medoids/nac_evt18_time1206ps/seed26711_time258.000ps.gro"
"$TASK_ROOT/ensemble/final_audit_job_61841413/assembled/medoids/nac_evt25_time1462ps/seed26711_time606.000ps.gro"
"$TASK_ROOT/ensemble/final_audit_job_61841413/assembled/medoids/nac_evt08_time1086ps/seed26711_time100.000ps.gro"
)
BUILDS=(
"$TASK_ROOT/ensemble/candidates/nac_evt18_time1206ps/build_job_61813799_6/build"
"$TASK_ROOT/ensemble/candidates/nac_evt25_time1462ps/build_job_61813799_11/build"
"$TASK_ROOT/ensemble/candidates/nac_evt08_time1086ps/build_job_61813799_5/build"
)
if [[ ! "$INDEX" =~ ^[0-2]$ ]]; then printf 'invalid candidate index: %s\n' "$INDEX" >&2; exit 2; fi
CANDIDATE="${CANDIDATES[$INDEX]}"
SOURCE="${SOURCES[$INDEX]}"
PARENT_BUILD="${BUILDS[$INDEX]}"
ATTEMPT="${SLURM_ARRAY_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-$INDEX}_${SLURM_JOB_ID:-manual}"
OUT="$TASK_ROOT/a1_activated_nac_20260726/full_system_preflight/attempt_$ATTEMPT/$CANDIDATE"
if [[ -e "$OUT" ]]; then printf 'refusing to overwrite %s\n' "$OUT" >&2; exit 2; fi
mkdir -p "$OUT"

EVENT=nylc_a1_full_system_preflight
COMMAND="run_nylc_a1_full_system_preflight.sh $INDEX"
STATE=FAIL_TECHNICAL
DETAIL="candidate=$CANDIDATE; initializing"
append_history() {
  local now
  now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
  (
    flock -x 9
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$now" "$EVENT" "$COMMAND" "$STATE" "$GITHUB_COMMIT" "$DETAIL" >>"$TASK_ROOT/run_history.tsv"
    "$PY" - "$now" "$EVENT" "$COMMAND" "$STATE" "$GITHUB_COMMIT" "$DETAIL" "$OUT" >>"$TASK_ROOT/run_history.jsonl" <<'PY'
import json,sys
now,event,command,state,commit,detail,output=sys.argv[1:]
print(json.dumps({"time":now,"event":event,"command":command,"state":state,
"git_commit":commit,"detail":detail,"output":output},sort_keys=True))
PY
  ) 9>"$TASK_ROOT/.run_history.lock"
}
finish() {
  code=$?
  trap - EXIT
  if [[ $code -ne 0 ]]; then DETAIL="candidate=$CANDIDATE; exit_code=$code; output=$OUT"; fi
  append_history
  exit "$code"
}
trap finish EXIT

"$PY" - "$PARAM/PATCH_AUDIT.json" "$PARAM/parameter_provenance.json" <<'PY'
import json,sys
patch=json.load(open(sys.argv[1])); provenance=json.load(open(sys.argv[2]))
if patch.get("status")!="PASS_A1_SCREENING_PATCH": raise SystemExit("parameter patch is not PASS_A1_SCREENING_PATCH")
if provenance.get("validation_status")!="PASS": raise SystemExit("parameter provenance is not PASS")
if provenance.get("parmchk2_unresolved_count")!=0: raise SystemExit("unresolved parmchk2 terms")
PY

"$PY" "$FLOW/scripts/build_nylc_a1_full_system.py" \
 --source-gro "$SOURCE" --parent-chain-itp "$PARENT_BUILD/topol_Protein_chain_H.itp" \
 --model-top "$PARAM/NTA1_CAP.gromacs.top" --model-mol2 "$PARAM/model/NTA1_CAP.input.mol2" \
 --chain-first-global-atom 8949 --candidate-id "$CANDIDATE" --output-dir "$OUT/patch"

mkdir "$OUT/system"
cp "$PARENT_BUILD"/*.itp "$OUT/system/"
cp "$PARENT_BUILD/topol.top" "$OUT/system/"
cp "$OUT/patch/topol_Protein_chain_H.itp" "$OUT/system/topol_Protein_chain_H.itp"
cp "$OUT/patch/system_A1.gro" "$OUT/system/system_A1.gro"
cat >"$OUT/system/preflight.mdp" <<'EOF'
integrator = md
dt = 0.001
nsteps = 0
cutoff-scheme = Verlet
nstlist = 10
rlist = 1.0
rcoulomb = 1.0
rvdw = 1.0
coulombtype = PME
vdwtype = Cut-off
pbc = xyz
continuation = yes
gen-vel = no
comm-mode = None
constraints = h-bonds
constraint-algorithm = lincs
nstenergy = 1
EOF

module purge
module load gromacs/2022.2-hpcx-gcc-7.3.1
export GMX_MAXBACKUP=-1
cd "$OUT/system"
"$GMX" grompp -f preflight.mdp -c system_A1.gro -p topol.top -o preflight.tpr \
 -pp processed.top -po preflight.expanded.mdp -maxwarn 0 >grompp.stdout 2>grompp.stderr
OMP_NUM_THREADS=1 "$GMX" mdrun -s preflight.tpr -deffnm gmx_single -ntomp 1 >mdrun.stdout 2>mdrun.stderr
printf 'Potential\n0\n' | "$GMX" energy -f gmx_single.edr -o potential.xvg >energy.stdout 2>energy.stderr
POTENTIAL="$(awk '!/^[@#]/{value=$2} END{print value}' potential.xvg)"
[[ "$POTENTIAL" =~ ^[-+0-9.eE]+$ ]]

"$PY" "$FLOW/scripts/audit_nylc_a1_full_system.py" \
 --gro system_A1.gro --chain-itp topol_Protein_chain_H.itp \
 --build-audit "$OUT/patch/A1_FULL_SYSTEM_BUILD.json" \
 --potential-energy-kj-mol "$POTENTIAL" --grompp-maxwarn-zero \
 --output A1_FULL_SYSTEM_PREFLIGHT.json
sha256sum topol_Protein_chain_H.itp system_A1.gro preflight.mdp preflight.tpr \
 gmx_single.edr potential.xvg A1_FULL_SYSTEM_PREFLIGHT.json >sha256.tsv
STATE=PASS_TECHNICAL
DETAIL="candidate=$CANDIDATE; status=PASS_A1_FULL_SYSTEM_PREFLIGHT; output=$OUT"
