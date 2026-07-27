#!/bin/bash
#SBATCH -J nylc_a1_nac
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 02:00:00
#SBATCH --array=0-8%9
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.err
set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
SNAP="$TASK_ROOT/a1_activated_nac_20260726/code_snapshots/nac_audit_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${SLURM_JOB_ID}"
mkdir -p "$SNAP/scripts" "$SNAP/manifests"
for name in run_nylc_a1_nac_audit.sh audit_nylc_a1_nac_replica.py generate_nylc_m1_ensemble_primitives.py analyze_nylc_m1_proton_geometry.py audit_nylc_m1_ensemble_replica.py analyze_nac_series.py; do
    cp "$FLOW/scripts/$name" "$SNAP/scripts/$name"
done
cp "$FLOW/manifests/nylc_a1_nac_audit_universe.json" "$SNAP/manifests/"
sha256sum "$SNAP"/scripts/* "$SNAP"/manifests/* >"$SNAP/SNAPSHOT_SHA256.tsv"
export A1_NAC_CODE_ROOT="$SNAP"
bash "$SNAP/scripts/run_nylc_a1_nac_audit.sh" "$SLURM_ARRAY_TASK_ID"
