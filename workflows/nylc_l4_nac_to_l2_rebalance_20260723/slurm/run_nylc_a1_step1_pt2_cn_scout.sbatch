#!/bin/bash
#SBATCH -J nylc_a1_pt2cn
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 01:00:00
#SBATCH --array=0-15%4
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.err

set -euo pipefail
: "${A1_PT2_CN_CODE_SOURCE:?set A1_PT2_CN_CODE_SOURCE}"
: "${A1_PT2_CN_GITHUB_COMMIT:?set A1_PT2_CN_GITHUB_COMMIT}"

CODE_SOURCE="$(realpath "$A1_PT2_CN_CODE_SOURCE")"
test -s "$CODE_SOURCE/GITHUB_COMMIT"
test -s "$CODE_SOURCE/SNAPSHOT_SHA256.tsv"
test "$(tr -d '[:space:]' <"$CODE_SOURCE/GITHUB_COMMIT")" = \
    "$A1_PT2_CN_GITHUB_COMMIT"
grep -Fq \
    "scripts/prepare_audit_nylc_a1_step1_pt2_cn_scout.py" \
    "$CODE_SOURCE/SNAPSHOT_SHA256.tsv"
grep -Fq \
    "scripts/run_nylc_a1_step1_pt2_cn_scout.sh" \
    "$CODE_SOURCE/SNAPSHOT_SHA256.tsv"
(
    cd "$CODE_SOURCE"
    sha256sum -c SNAPSHOT_SHA256.tsv
)

export A1_PT2_CN_CODE_ROOT="$CODE_SOURCE"
bash "$CODE_SOURCE/scripts/run_nylc_a1_step1_pt2_cn_scout.sh"
