#!/bin/bash
#SBATCH -J nylc_a1_param
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem-per-cpu=2000M
#SBATCH -t 02:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.err

set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
OUT_DIR="$TASK_ROOT/a1_activated_nac_20260726/parameterization/job_${SLURM_JOB_ID:?SLURM_JOB_ID is required}"
export GITHUB_COMMIT="${A1_GITHUB_COMMIT:-2c0ad36333c041cc8ced2dc520504f0e6d49d145}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
bash "$FLOW/scripts/run_nylc_a1_parameterization.sh" "$OUT_DIR"
