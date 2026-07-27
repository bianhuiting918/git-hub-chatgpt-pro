#!/bin/bash
#SBATCH -J nylc_a1_cpl
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -c 1
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 02:00:00
#SBATCH --array=0-2
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.err
set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
: "${A1_COUPLED_SCOUT_GITHUB_COMMIT:?set immutable GitHub commit}"
: "${A1_COUPLED_SCOUT_CODE_SOURCE:?set deployed immutable code root}"
CODE_SOURCE="$(realpath "$A1_COUPLED_SCOUT_CODE_SOURCE")"
[[ -s "$CODE_SOURCE/GITHUB_COMMIT" && -s "$CODE_SOURCE/SNAPSHOT_SHA256.tsv" ]]
[[ "$(tr -d '[:space:]' <"$CODE_SOURCE/GITHUB_COMMIT")" == "$A1_COUPLED_SCOUT_GITHUB_COMMIT" ]]
(
    cd "$CODE_SOURCE"
    sha256sum -c SNAPSHOT_SHA256.tsv
)
TARGETS=(1.30 1.35 1.40)
CARBONYL_TARGET="${TARGETS[${SLURM_ARRAY_TASK_ID}]}"
SNAP="$TASK_ROOT/a1_activated_nac_20260726/code_snapshots/a1_step1_coupled_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
[[ ! -e "$SNAP" ]] || { printf 'refusing to overwrite %s\n' "$SNAP" >&2; exit 2; }
mkdir -p "$SNAP/scripts"
for name in prepare_nylc_a1_step1_coupled_scout.py audit_nylc_a1_step1_coupled_scout.py run_nylc_a1_step1_coupled_scout.sh; do
    cp "$CODE_SOURCE/scripts/$name" "$SNAP/scripts/$name"
done
cp "$CODE_SOURCE/GITHUB_COMMIT" "$SNAP/GITHUB_COMMIT"
(
    cd "$SNAP"
    sha256sum scripts/* GITHUB_COMMIT >SNAPSHOT_SHA256.tsv
)
chmod -R a-w "$SNAP"
export A1_COUPLED_SCOUT_CODE_ROOT="$SNAP"
export A1_COUPLED_SCOUT_ATTEMPT="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
export A1_COUPLED_SCOUT_CARBONYL_TARGET="$CARBONYL_TARGET"
bash "$SNAP/scripts/run_nylc_a1_step1_coupled_scout.sh"
