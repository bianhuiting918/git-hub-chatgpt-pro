#!/bin/bash
#SBATCH -J nylc_a2qmwat
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -c 1
#SBATCH --mem-per-cpu=2500M
#SBATCH --time=08:00:00
#SBATCH --array=0-1
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.err

set -euo pipefail
: "${A1_STEP2_QMWATER_CODE_SOURCE:?set immutable workflow code source}"
: "${A1_STEP2_QMWATER_GITHUB_COMMIT:?set full runtime GitHub commit}"

CODE_SOURCE="$(realpath "$A1_STEP2_QMWATER_CODE_SOURCE")"
[[ "$A1_STEP2_QMWATER_GITHUB_COMMIT" =~ ^[0-9a-f]{40}$ ]]
test -s "$CODE_SOURCE/GITHUB_COMMIT"
test -s "$CODE_SOURCE/SNAPSHOT_SHA256.tsv"
test "$(tr -d '[:space:]' <"$CODE_SOURCE/GITHUB_COMMIT")" =     "$A1_STEP2_QMWATER_GITHUB_COMMIT"

for file in     scripts/prepare_audit_nylc_a1_step2_qmwater_endpoint.py     scripts/run_nylc_a1_step2_qmwater_endpoint.sh     slurm/run_nylc_a1_step2_qmwater_endpoint.sbatch     scripts/prepare_audit_nylc_a1_step1_acyl_endpoint_stability.py     scripts/prepare_audit_nylc_a1_step1_pt2_cn_scout.py     audit/nylc_a1_acyl_endpoint_connectivity_sensitivity_62216380_v1.json
do
    grep -Fq "$file" "$CODE_SOURCE/SNAPSHOT_SHA256.tsv"
done
(
    cd "$CODE_SOURCE"
    sha256sum -c SNAPSHOT_SHA256.tsv
)

export A1_STEP2_QMWATER_CODE_ROOT="$CODE_SOURCE"
bash "$CODE_SOURCE/scripts/run_nylc_a1_step2_qmwater_endpoint.sh"
