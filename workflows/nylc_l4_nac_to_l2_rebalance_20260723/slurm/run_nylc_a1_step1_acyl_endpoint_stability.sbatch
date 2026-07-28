#!/bin/bash
#SBATCH -J nylc_a1acyl
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem-per-cpu=2500M
#SBATCH --time=08:00:00
#SBATCH --array=0-1%2
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.err

set -euo pipefail
: "${A1_ACYL_ENDPOINT_CODE_SOURCE:?set A1_ACYL_ENDPOINT_CODE_SOURCE}"
: "${A1_ACYL_ENDPOINT_GITHUB_COMMIT:?set A1_ACYL_ENDPOINT_GITHUB_COMMIT}"

CODE_SOURCE="$(realpath "$A1_ACYL_ENDPOINT_CODE_SOURCE")"
test -s "$CODE_SOURCE/GITHUB_COMMIT"
test -s "$CODE_SOURCE/SNAPSHOT_SHA256.tsv"
test "$(tr -d '[:space:]' <"$CODE_SOURCE/GITHUB_COMMIT")" = "$A1_ACYL_ENDPOINT_GITHUB_COMMIT"
for file in     scripts/prepare_audit_nylc_a1_step1_acyl_endpoint_stability.py     scripts/run_nylc_a1_step1_acyl_endpoint_stability.sh     slurm/run_nylc_a1_step1_acyl_endpoint_stability.sbatch     scripts/prepare_audit_nylc_a1_step1_pt2_cn_scout.py
do
    grep -Fq "$file" "$CODE_SOURCE/SNAPSHOT_SHA256.tsv"
done
(
    cd "$CODE_SOURCE"
    sha256sum -c SNAPSHOT_SHA256.tsv
)

export A1_ACYL_ENDPOINT_CODE_ROOT="$CODE_SOURCE"
bash "$CODE_SOURCE/scripts/run_nylc_a1_step1_acyl_endpoint_stability.sh"
