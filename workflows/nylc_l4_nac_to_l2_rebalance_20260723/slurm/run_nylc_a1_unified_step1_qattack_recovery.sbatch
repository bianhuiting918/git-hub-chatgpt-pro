#!/bin/bash
#SBATCH -J nylc_a1_qrec
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
: "${A1_QATTACK_RECOVERY_GITHUB_COMMIT:?set immutable GitHub commit}"
: "${A1_QATTACK_RECOVERY_CODE_SOURCE:?set deployed immutable code root}"
CODE_SOURCE="$(realpath "$A1_QATTACK_RECOVERY_CODE_SOURCE")"
[[ -s "$CODE_SOURCE/GITHUB_COMMIT" && -s "$CODE_SOURCE/SNAPSHOT_SHA256.tsv" ]]
[[ "$(tr -d '[:space:]' <"$CODE_SOURCE/GITHUB_COMMIT")" == "$A1_QATTACK_RECOVERY_GITHUB_COMMIT" ]]
(
    cd "$CODE_SOURCE"
    sha256sum -c SNAPSHOT_SHA256.tsv
)
FORCES=(50 100 200)
FORCE="${FORCES[${SLURM_ARRAY_TASK_ID}]}"
SNAP="$TASK_ROOT/a1_activated_nac_20260726/code_snapshots/a1_qattack_recovery_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
[[ ! -e "$SNAP" ]] || { printf 'refusing to overwrite %s\n' "$SNAP" >&2; exit 2; }
mkdir -p "$SNAP/scripts"
for name in prepare_nylc_a1_unified_step1_qattack_recovery.py audit_nylc_a1_unified_step1_qattack_recovery.py run_nylc_a1_unified_step1_qattack_recovery.sh; do
    cp "$CODE_SOURCE/scripts/$name" "$SNAP/scripts/$name"
done
cp "$CODE_SOURCE/GITHUB_COMMIT" "$SNAP/GITHUB_COMMIT"
(
    cd "$SNAP"
    sha256sum scripts/* GITHUB_COMMIT >SNAPSHOT_SHA256.tsv
)
chmod -R a-w "$SNAP"
export A1_QATTACK_RECOVERY_CODE_ROOT="$SNAP"
export A1_QATTACK_RECOVERY_ATTEMPT="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
export A1_QATTACK_RECOVERY_FORCE="$FORCE"
bash "$SNAP/scripts/run_nylc_a1_unified_step1_qattack_recovery.sh"
