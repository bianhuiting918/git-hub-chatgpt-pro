#!/bin/bash
#SBATCH -J nylc_a1_em
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=2000M
#SBATCH -t 02:00:00
#SBATCH --array=0-2%3
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.err
set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
SNAP="$TASK_ROOT/a1_activated_nac_20260726/code_snapshots/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${SLURM_JOB_ID}"
mkdir -p "$SNAP"
for name in run_nylc_a1_em.sh audit_nylc_a1_em.py audit_nylc_a1_full_system.py build_nylc_a1_full_system.py; do
    cp "$FLOW/scripts/$name" "$SNAP/$name"
done
sha256sum "$SNAP"/run_nylc_a1_em.sh "$SNAP"/*.py >"$SNAP/SNAPSHOT_SHA256.tsv"
export A1_SCRIPT_ROOT="$SNAP"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
bash "$SNAP/run_nylc_a1_em.sh" "$SLURM_ARRAY_TASK_ID"
