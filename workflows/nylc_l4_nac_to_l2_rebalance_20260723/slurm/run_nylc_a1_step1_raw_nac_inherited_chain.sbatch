#!/usr/bin/env bash
#SBATCH --job-name=nylc_a1nacchain
#SBATCH --array=0-3
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2500M
#SBATCH --time=48:00:00
#SBATCH --output=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-nylc_a1nacchain-%A_%a.out
#SBATCH --error=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-nylc_a1nacchain-%A_%a.err

set -euo pipefail
: "${A1_RAW_NAC_CHAIN_CODE_ROOT:?export immutable workflow root}"
: "${A1_RAW_NAC_CHAIN_GITHUB_COMMIT:?export full GitHub commit}"
test -s "$A1_RAW_NAC_CHAIN_CODE_ROOT/GITHUB_COMMIT"
test "$(tr -d '\r\n' < "$A1_RAW_NAC_CHAIN_CODE_ROOT/GITHUB_COMMIT")" = "$A1_RAW_NAC_CHAIN_GITHUB_COMMIT"
test -s "$A1_RAW_NAC_CHAIN_CODE_ROOT/SNAPSHOT_SHA256.tsv"
exec bash "$A1_RAW_NAC_CHAIN_CODE_ROOT/scripts/run_nylc_a1_step1_raw_nac_inherited_chain.sh"
