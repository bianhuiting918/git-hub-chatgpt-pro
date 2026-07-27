#!/bin/bash
#SBATCH -J nylc_a1_equil
#SBATCH -p xahdnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=dcu:1
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 06:00:00
#SBATCH --array=0-8%9
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.err
set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
SNAP="$TASK_ROOT/a1_activated_nac_20260726/code_snapshots/equil_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${SLURM_JOB_ID}"
mkdir -p "$SNAP/mdp"
cp "$FLOW/scripts/run_nylc_a1_equilibration.sh" "$SNAP/"
for name in nvt50_m1.mdp nvt150_m1.mdp nvt300_m1.mdp npt300r_m1.mdp npt300rel_m1.mdp npt300free_m1.mdp; do
    cp "$FLOW/mdp/$name" "$SNAP/mdp/$name"
done
sha256sum "$SNAP/run_nylc_a1_equilibration.sh" "$SNAP"/mdp/*.mdp >"$SNAP/SNAPSHOT_SHA256.tsv"
export A1_EQ_CODE_ROOT="$SNAP"
bash "$SNAP/run_nylc_a1_equilibration.sh" "$SLURM_ARRAY_TASK_ID"
