#!/bin/bash
#SBATCH -J nylc_a1_pt2
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 1
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 01:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.err
set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
: "${A1_PT2_NAC_GITHUB_COMMIT:?set immutable GitHub commit}"
: "${A1_PT2_NAC_CODE_SOURCE:?set deployed immutable code root}"
CODE_SOURCE="$(realpath "$A1_PT2_NAC_CODE_SOURCE")"
[[ -s "$CODE_SOURCE/GITHUB_COMMIT" && -s "$CODE_SOURCE/SNAPSHOT_SHA256.tsv" ]]
[[ "$(tr -d '[:space:]' <"$CODE_SOURCE/GITHUB_COMMIT")" == "$A1_PT2_NAC_GITHUB_COMMIT" ]]
(
    cd "$CODE_SOURCE"
    sha256sum -c SNAPSHOT_SHA256.tsv
)
SNAP="$TASK_ROOT/a1_activated_nac_20260726/code_snapshots/a1_pt2_nac_${SLURM_JOB_ID}"
[[ ! -e "$SNAP" ]] || { printf 'refusing to overwrite %s\n' "$SNAP" >&2; exit 2; }
mkdir -p "$SNAP/scripts"
cp "$CODE_SOURCE/scripts/analyze_nylc_a1_pt2_nac_ensemble.py" "$SNAP/scripts/"
cp "$CODE_SOURCE/scripts/run_nylc_a1_pt2_nac_ensemble.sh" "$SNAP/scripts/"
cp "$CODE_SOURCE/GITHUB_COMMIT" "$SNAP/GITHUB_COMMIT"
(
    cd "$SNAP"
    sha256sum scripts/* GITHUB_COMMIT >SNAPSHOT_SHA256.tsv
)
chmod -R a-w "$SNAP"
export A1_PT2_NAC_CODE_ROOT="$SNAP"
export A1_PT2_NAC_ATTEMPT="${SLURM_JOB_ID}"
bash "$SNAP/scripts/run_nylc_a1_pt2_nac_ensemble.sh"
