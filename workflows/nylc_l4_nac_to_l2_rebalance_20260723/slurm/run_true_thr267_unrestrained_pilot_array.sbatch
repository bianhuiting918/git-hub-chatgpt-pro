#!/bin/bash
#SBATCH -J nylc_t267_free
#SBATCH -p xahdnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=dcu:1
#SBATCH --mem-per-cpu=2000M
#SBATCH -t 01:00:00
#SBATCH --array=0-2%3
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/true_thr267_free_%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/true_thr267_free_%A_%a.err

set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PY=/work/home/acshdt1dks/anaconda3/bin/python3.9
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
BRANCH=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/apo_gate_l4_three_carbonyl_20260715/cyclic_gate_nac_20260719/branches/nylc_gyaq/C18
SELECTED="$TASK/candidates/nylc_C18_trueT267_recapture/selected_low_energy_nac/release2_61705002/selected_lowest_potential_nac.gro"
EXPECTED_SHA=d0f157575193c4662f5463c14f5d9075e817c18142aa74f0f7456b8bc732b4f7
seed=$((26701 + SLURM_ARRAY_TASK_ID))
root="$TASK/candidates/nylc_C18_trueT267_recapture/unrestrained_pilots/job_${SLURM_ARRAY_JOB_ID}/seed_${seed}"
mkdir -p "$root" "$TASK/logs"
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
export GMX_MAXBACKUP=-1

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  (
    flock -x 9
    printf '%s\ttrue_thr267_unrestrained_pilot\t%s\tcandidate=nylc_C18_trueT267_recapture;job=%s;task=%s;seed=%s;%s\n' \
      "$stamp" "$state" "${SLURM_ARRAY_JOB_ID}" "${SLURM_ARRAY_TASK_ID}" "$seed" "$detail" >&9
  ) 9>>"$TASK/run_history.tsv"
  (
    flock -x 8
    "$PY" - "$stamp" "$state" "$detail" "${SLURM_ARRAY_JOB_ID}" "${SLURM_ARRAY_TASK_ID}" "$seed" >&8 <<'PY'
import json,sys
stamp,state,detail,job,task,seed=sys.argv[1:]
print(json.dumps({
  "timestamp":stamp,
  "event":"true_thr267_unrestrained_pilot",
  "state":state,
  "detail":detail,
  "candidate_id":"nylc_C18_trueT267_recapture",
  "slurm_array_job_id":job,
  "slurm_array_task_id":task,
  "velocity_seed":int(seed),
},sort_keys=True))
PY
  ) 8>>"$TASK/run_history.jsonl"
}

failure=1
stage=initialization
trap 'rc=$?; if [[ $failure -ne 0 ]]; then record FAIL "stage=$stage;rc=$rc"; fi; exit $rc' EXIT
test "$(sha256sum "$SELECTED" | awk '{print $1}')" = "$EXPECTED_SHA"
for f in free_nvt20.tpr free_nvt20.gro free_npt100.tpr free_npt100.gro unrestrained_pilot_audit.json; do
  test ! -e "$root/$f"
done

cat >"$root/free_nvt20.mdp" <<EOF
integrator               = md
dt                       = 0.002
nsteps                   = 10000
continuation             = no
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
pcoupl                   = no
pbc                      = xyz
gen-vel                  = yes
gen-temp                 = 300
gen-seed                 = $seed
nstxout-compressed       = 500
nstenergy                = 500
nstlog                   = 500
comm-mode                = Linear
nstcomm                  = 100
pull                     = no
EOF

cat >"$root/free_npt100.mdp" <<'EOF'
integrator               = md
dt                       = 0.002
nsteps                   = 50000
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
nstxout-compressed       = 500
nstenergy                = 500
nstlog                   = 500
comm-mode                = Linear
nstcomm                  = 100
pull                     = no
EOF

record START "source=restrained_lowest_potential_NAC;window=fully_unrestrained"
stage=grompp_nvt20
(
  cd "$BRANCH"
  "$GMX" grompp -f "$root/free_nvt20.mdp" -c "$SELECTED" -p "$BRANCH/topol.top" \
    -o "$root/free_nvt20.tpr" -po "$root/free_nvt20.processed.mdp" -maxwarn 0 \
    >"$root/grompp_nvt20.stdout" 2>"$root/grompp_nvt20.stderr"
)
stage=mdrun_nvt20
mpirun -np 1 "$GMX" mdrun -s "$root/free_nvt20.tpr" -deffnm "$root/free_nvt20" \
  -ntomp "${SLURM_CPUS_PER_TASK:-8}" -nb gpu -pme gpu -bonded gpu -pin on \
  >"$root/mdrun_nvt20.stdout" 2>"$root/mdrun_nvt20.stderr"

stage=grompp_npt100
(
  cd "$BRANCH"
  "$GMX" grompp -f "$root/free_npt100.mdp" -c "$root/free_nvt20.gro" -t "$root/free_nvt20.cpt" \
    -p "$BRANCH/topol.top" -o "$root/free_npt100.tpr" -po "$root/free_npt100.processed.mdp" \
    -maxwarn 0 >"$root/grompp_npt100.stdout" 2>"$root/grompp_npt100.stderr"
)
stage=mdrun_npt100
mpirun -np 1 "$GMX" mdrun -s "$root/free_npt100.tpr" -deffnm "$root/free_npt100" \
  -ntomp "${SLURM_CPUS_PER_TASK:-8}" -nb gpu -pme gpu -bonded gpu -pin on \
  >"$root/mdrun_npt100.stdout" 2>"$root/mdrun_npt100.stderr"

stage=audit
"$GMX" distance -f "$root/free_npt100.xtc" -s "$root/free_npt100.tpr" \
  -select 'atomnr 8961 plus atomnr 10297' -oall "$root/nac_distance.xvg" \
  >"$root/nac_distance.stdout" 2>"$root/nac_distance.stderr"
printf '[ nac_angle ]\n8961 10297 10298\n' >"$root/nac_angle.ndx"
printf '0\n' | "$GMX" angle -f "$root/free_npt100.xtc" -n "$root/nac_angle.ndx" \
  -ov "$root/nac_angle.xvg" -type angle >"$root/nac_angle.stdout" 2>"$root/nac_angle.stderr"
printf 'Potential\n0\n' | "$GMX" energy -f "$root/free_npt100.edr" \
  -o "$root/potential_energy.xvg" >"$root/potential.stdout" 2>"$root/potential.stderr"
printf 'Temperature\nPressure\nVolume\n0\n' | "$GMX" energy -f "$root/free_npt100.edr" \
  -o "$root/thermo.xvg" >"$root/thermo.stdout" 2>"$root/thermo.stderr"
"$GMX" distance -f "$root/free_npt100.xtc" -s "$root/free_npt100.tpr" -n "$BRANCH/cycle.ndx" \
  -select 'com of group "Core" plus com of group "Gate"' -oxyz "$root/gate_core_vector.xvg" \
  >"$root/gate.stdout" 2>"$root/gate.stderr"

"$PY" "$FLOW/scripts/analyze_nac_series.py" \
  --distance-xvg "$root/nac_distance.xvg" --angle-xvg "$root/nac_angle.xvg" \
  --potential-xvg "$root/potential_energy.xvg" --distance-max-nm 0.35 \
  --angle-min-deg 95 --angle-max-deg 115 --output "$root/nac_energy_series.json"
"$PY" "$FLOW/scripts/audit_true_thr267_unrestrained_pilot.py" \
  --run-root "$root" --seed "$seed" --array-job-id "${SLURM_ARRAY_JOB_ID}" \
  >"$root/audit.stdout" 2>"$root/audit.stderr"
record PASS_TECHNICAL "audit=$root/unrestrained_pilot_audit.json"
failure=0
trap - EXIT
