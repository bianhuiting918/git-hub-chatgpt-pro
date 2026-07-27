#!/usr/bin/env bash
set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
GMX=/public/software/apps/gromacs/2022.2/hpcx-gcc7.3.1/bin/gmx_mpi
INDEX="${1:?usage: run_nylc_a1_em.sh INDEX}"; GITHUB_COMMIT="${A1_GITHUB_COMMIT:-unknown}"
CANDIDATES=(nac_evt18_time1206ps nac_evt25_time1462ps nac_evt08_time1086ps)
INPUTS=(
"$TASK_ROOT/a1_activated_nac_20260726/full_system_preflight/attempt_61966928_0_61966930/nac_evt18_time1206ps"
"$TASK_ROOT/a1_activated_nac_20260726/full_system_preflight/attempt_61966928_1_61966931/nac_evt25_time1462ps"
"$TASK_ROOT/a1_activated_nac_20260726/full_system_preflight/attempt_61966928_2_61966928/nac_evt08_time1086ps"
)
[[ "$INDEX" =~ ^[0-2]$ ]] || exit 2
CANDIDATE="${CANDIDATES[$INDEX]}"; INPUT="${INPUTS[$INDEX]}"
ATTEMPT="${SLURM_ARRAY_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-$INDEX}_${SLURM_JOB_ID:-manual}"
OUT="$TASK_ROOT/a1_activated_nac_20260726/em/attempt_$ATTEMPT/$CANDIDATE"
[[ ! -e "$OUT" ]] || { printf 'refusing to overwrite %s\n' "$OUT" >&2; exit 2; }
mkdir -p "$OUT/input" "$OUT/em_hrelax" "$OUT/em_free"
EVENT=nylc_a1_em; COMMAND="run_nylc_a1_em.sh $INDEX"; STATE=FAIL_TECHNICAL; DETAIL="candidate=$CANDIDATE; initializing"
append_history(){ local now; now="$(date '+%Y-%m-%dT%H:%M:%S%z')"; (
flock -x 9
printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$now" "$EVENT" "$COMMAND" "$STATE" "$GITHUB_COMMIT" "$DETAIL" >>"$TASK_ROOT/run_history.tsv"
"$PY" - "$now" "$EVENT" "$COMMAND" "$STATE" "$GITHUB_COMMIT" "$DETAIL" "$OUT" >>"$TASK_ROOT/run_history.jsonl" <<'PY'
import json,sys
now,event,command,state,commit,detail,output=sys.argv[1:]
print(json.dumps({"time":now,"event":event,"command":command,"state":state,"git_commit":commit,"detail":detail,"output":output},sort_keys=True))
PY
) 9>"$TASK_ROOT/.run_history.lock"; }
finish(){ code=$?; trap - EXIT; if [[ $code -ne 0 ]]; then DETAIL="candidate=$CANDIDATE; exit_code=$code; output=$OUT"; fi; append_history; exit "$code"; }
trap finish EXIT
"$PY" - "$INPUT/system/A1_FULL_SYSTEM_PREFLIGHT.json" <<'PY'
import json,sys
if json.load(open(sys.argv[1])).get("status")!="PASS_A1_FULL_SYSTEM_PREFLIGHT": raise SystemExit("input preflight is not PASS")
PY
cp "$INPUT/system"/*.itp "$OUT/input/"
cp "$INPUT/system/topol.top" "$INPUT/system/system_A1.gro" "$OUT/input/"
cp "$INPUT/patch/A1_FULL_SYSTEM_BUILD.json" "$OUT/input/"

cat >"$OUT/em_hrelax/em.mdp" <<'EOF'
define = -DPOSRES -DPOSRES_L2_1000
integrator = steep
emtol = 1000
emstep = 0.0002
nsteps = 10000
cutoff-scheme = Verlet
nstlist = 20
rlist = 1.0
coulombtype = PME
rcoulomb = 1.0
vdwtype = Cut-off
rvdw = 1.0
DispCorr = EnerPres
pbc = xyz
constraints = none
nstenergy = 100
nstlog = 100
EOF
cat >"$OUT/em_free/em.mdp" <<'EOF'
integrator = steep
emtol = 500
emstep = 0.0001
nsteps = 50000
cutoff-scheme = Verlet
nstlist = 20
rlist = 1.0
coulombtype = PME
rcoulomb = 1.0
vdwtype = Cut-off
rvdw = 1.0
DispCorr = EnerPres
pbc = xyz
constraints = none
nstenergy = 100
nstlog = 100
EOF
module purge
module load gromacs/2022.2-hpcx-gcc-7.3.1
export GMX_MAXBACKUP=-1 OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
cd "$OUT/em_hrelax"
"$GMX" grompp -f em.mdp -c ../input/system_A1.gro -r ../input/system_A1.gro -p ../input/topol.top -o run.tpr -maxwarn 0 >grompp.stdout 2>grompp.stderr
"$GMX" mdrun -s run.tpr -deffnm run -ntomp "$OMP_NUM_THREADS" >mdrun.stdout 2>mdrun.stderr
cd "$OUT/em_free"
"$GMX" grompp -f em.mdp -c ../em_hrelax/run.gro -p ../input/topol.top -o run.tpr -maxwarn 0 >grompp.stdout 2>grompp.stderr
"$GMX" mdrun -s run.tpr -deffnm run -ntomp "$OMP_NUM_THREADS" >mdrun.stdout 2>mdrun.stderr
"$PY" "$FLOW/scripts/audit_nylc_a1_em.py" --gro run.gro --chain-itp ../input/topol_Protein_chain_H.itp --build-audit ../input/A1_FULL_SYSTEM_BUILD.json --log run.log --output A1_EM_AUDIT.json
sha256sum run.gro run.log run.edr A1_EM_AUDIT.json >sha256.tsv
STATE=PASS_TECHNICAL; DETAIL="candidate=$CANDIDATE; status=PASS_A1_EM; output=$OUT"
