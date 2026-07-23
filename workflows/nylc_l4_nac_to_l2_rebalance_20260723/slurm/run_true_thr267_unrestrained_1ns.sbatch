#!/bin/bash
#SBATCH -J nylc_t267_1ns
#SBATCH -p xahdnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=dcu:1
#SBATCH --mem-per-cpu=2000M
#SBATCH -t 02:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/true_thr267_1ns_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/true_thr267_1ns_%j.err

set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PY=/work/home/acshdt1dks/anaconda3/bin/python3.9
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
BRANCH=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/apo_gate_l4_three_carbonyl_20260715/cyclic_gate_nac_20260719/branches/nylc_gyaq/C18
PARENT="$TASK/candidates/nylc_C18_trueT267_recapture/unrestrained_pilots/job_61705307/seed_26703"
PARENT_GRO_SHA=ab6fb4083f02e9bf1c2af3677500ceccc139ba08882e32750c63136bf3320b6f
PARENT_CPT_SHA=0e8a93774b8c0836bee7a36b041676112c014cd1eec9e96055a2ade971b6d459
job="${SLURM_JOB_ID:?}"
root="$TASK/candidates/nylc_C18_trueT267_recapture/unrestrained_1ns/job_$job"
mkdir -p "$root" "$TASK/logs"
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
export GMX_MAXBACKUP=-1

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  (
    flock -x 9
    printf '%s\ttrue_thr267_unrestrained_1ns\t%s\tcandidate=nylc_C18_trueT267_recapture;job=%s;%s\n' \
      "$stamp" "$state" "$job" "$detail" >&9
  ) 9>>"$TASK/run_history.tsv"
  (
    flock -x 8
    "$PY" - "$stamp" "$state" "$detail" "$job" >&8 <<'PY'
import json,sys
stamp,state,detail,job=sys.argv[1:]
print(json.dumps({
  "timestamp":stamp,
  "event":"true_thr267_unrestrained_1ns",
  "state":state,
  "detail":detail,
  "candidate_id":"nylc_C18_trueT267_recapture",
  "slurm_job_id":job,
  "parent_array_job_id":"61705307",
  "parent_seed":26703,
},sort_keys=True))
PY
  ) 8>>"$TASK/run_history.jsonl"
}

failure=1
stage=initialization
trap 'rc=$?; if [[ $failure -ne 0 ]]; then record FAIL "stage=$stage;rc=$rc"; fi; exit $rc' EXIT
test "$(sha256sum "$PARENT/free_npt100.gro" | awk '{print $1}')" = "$PARENT_GRO_SHA"
test "$(sha256sum "$PARENT/free_npt100.cpt" | awk '{print $1}')" = "$PARENT_CPT_SHA"
test "$("$PY" -c "import json;print(json.load(open('$PARENT/unrestrained_pilot_audit.json'))['scientific_status'])")" = PASS_UNRESTRAINED_PILOT_NAC_PRESENT
for f in run.tpr run.gro run.xtc unrestrained_1ns_audit.json; do test ! -e "$root/$f"; done

cat >"$root/run.mdp" <<'EOF'
integrator               = md
dt                       = 0.002
nsteps                   = 500000
continuation             = yes
constraints              = h-bonds
constraint-algorithm     = lincs
lincs-iter               = 1
lincs-order              = 4
cutoff-scheme            = Verlet
nstlist                  = 20
rlist                    = 1.0
coulombtype              = PME
rcoulomb                 = 1.0
vdwtype                  = Cut-off
rvdw                     = 1.0
DispCorr                 = EnerPres
tcoupl                   = V-rescale
tc-grps                  = System
tau-t                    = 1.0
ref-t                    = 300
pcoupl                   = Parrinello-Rahman
pcoupltype               = isotropic
tau-p                    = 5.0
ref-p                    = 1.0
compressibility          = 4.5e-5
pbc                      = xyz
gen-vel                  = no
nstxout-compressed       = 1000
nstenergy                = 1000
nstlog                   = 5000
comm-mode                = Linear
nstcomm                  = 100
pull                     = no
EOF

record START "parent=61705307_seed26703;window=fully_unrestrained_NPT_1ns"
stage=grompp
(
  cd "$BRANCH"
  "$GMX" grompp -f "$root/run.mdp" -c "$PARENT/free_npt100.gro" -t "$PARENT/free_npt100.cpt" \
    -p "$BRANCH/topol.top" -o "$root/run.tpr" -po "$root/run.processed.mdp" \
    -maxwarn 0 >"$root/grompp.stdout" 2>"$root/grompp.stderr"
)
stage=mdrun
mpirun -np 1 "$GMX" mdrun -s "$root/run.tpr" -deffnm "$root/run" \
  -ntomp "${SLURM_CPUS_PER_TASK:-8}" -nb gpu -pme gpu -bonded gpu -pin on \
  >"$root/mdrun.stdout" 2>"$root/mdrun.stderr"

stage=audit
"$GMX" distance -f "$root/run.xtc" -s "$root/run.tpr" \
  -select 'atomnr 8961 plus atomnr 10297' -oall "$root/nac_distance.xvg" \
  >"$root/nac_distance.stdout" 2>"$root/nac_distance.stderr"
printf '[ nac_angle ]\n8961 10297 10298\n' >"$root/nac_angle.ndx"
printf '0\n' | "$GMX" angle -f "$root/run.xtc" -n "$root/nac_angle.ndx" \
  -ov "$root/nac_angle.xvg" -type angle >"$root/nac_angle.stdout" 2>"$root/nac_angle.stderr"
printf 'Potential\n0\n' | "$GMX" energy -f "$root/run.edr" \
  -o "$root/potential_energy.xvg" >"$root/potential.stdout" 2>"$root/potential.stderr"
printf 'Temperature\nPressure\nVolume\n0\n' | "$GMX" energy -f "$root/run.edr" \
  -o "$root/thermo.xvg" >"$root/thermo.stdout" 2>"$root/thermo.stderr"
"$GMX" distance -f "$root/run.xtc" -s "$root/run.tpr" -n "$BRANCH/cycle.ndx" \
  -select 'com of group "Core" plus com of group "Gate"' -oxyz "$root/gate_core_vector.xvg" \
  >"$root/gate.stdout" 2>"$root/gate.stderr"
"$PY" "$FLOW/scripts/analyze_nac_series.py" \
  --distance-xvg "$root/nac_distance.xvg" --angle-xvg "$root/nac_angle.xvg" \
  --potential-xvg "$root/potential_energy.xvg" --distance-max-nm 0.35 \
  --angle-min-deg 95 --angle-max-deg 115 --output "$root/nac_energy_series.json"
"$PY" "$FLOW/scripts/audit_true_thr267_unrestrained_pilot.py" \
  --run-root "$root" --seed 26703 --array-job-id "$job" --stem run \
  --window-type fully_unrestrained_NPT_1ns --output-name unrestrained_1ns_audit.json \
  >"$root/audit.stdout" 2>"$root/audit.stderr"
record PASS_TECHNICAL "audit=$root/unrestrained_1ns_audit.json"
failure=0
trap - EXIT
