#!/bin/bash
#SBATCH -J nylc_a1u_dftb3
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -c 1
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 02:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.err
set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
: "${A1_UNIFIED_CORE_GITHUB_COMMIT:?set immutable GitHub commit}"
: "${A1_UNIFIED_CORE_DFTB3_CODE_SOURCE:?set deployed immutable code root}"
printf 'expanded_preflight commit=%s code_source=%s\n' "$A1_UNIFIED_CORE_GITHUB_COMMIT" "$A1_UNIFIED_CORE_DFTB3_CODE_SOURCE"
if [[ ! -d "$A1_UNIFIED_CORE_DFTB3_CODE_SOURCE" ]]; then
    printf 'code source is not visible on compute node: %s\n' "$A1_UNIFIED_CORE_DFTB3_CODE_SOURCE" >&2
    exit 3
fi
CODE_SOURCE="$(realpath "$A1_UNIFIED_CORE_DFTB3_CODE_SOURCE")"
if [[ ! -s "$CODE_SOURCE/GITHUB_COMMIT" || ! -s "$CODE_SOURCE/SNAPSHOT_SHA256.tsv" ]]; then
    printf 'snapshot metadata is missing under %s\n' "$CODE_SOURCE" >&2
    exit 3
fi
if [[ "$(tr -d '[:space:]' <"$CODE_SOURCE/GITHUB_COMMIT")" != "$A1_UNIFIED_CORE_GITHUB_COMMIT" ]]; then
    printf 'snapshot commit mismatch under %s\n' "$CODE_SOURCE" >&2
    exit 3
fi
(
    cd "$CODE_SOURCE"
    sha256sum -c SNAPSHOT_SHA256.tsv
)
SNAP="$TASK_ROOT/a1_activated_nac_20260726/code_snapshots/a1_unified_core_dftb3_${SLURM_JOB_ID}"
[[ ! -e "$SNAP" ]] || { printf 'refusing to overwrite code snapshot %s\n' "$SNAP" >&2; exit 2; }
mkdir -p "$SNAP/scripts"
for name in prepare_nylc_a1_unified_core_dftb3_preflight.py audit_nylc_a1_unified_core_dftb3_smoke.py run_nylc_a1_unified_core_dftb3_preflight.sh; do
    grep -Fq "scripts/$name" "$CODE_SOURCE/SNAPSHOT_SHA256.tsv"
    cp "$CODE_SOURCE/scripts/$name" "$SNAP/scripts/$name"
done
cp "$CODE_SOURCE/GITHUB_COMMIT" "$SNAP/GITHUB_COMMIT"
(
    cd "$SNAP"
    sha256sum scripts/* GITHUB_COMMIT >SNAPSHOT_SHA256.tsv
)
chmod -R a-w "$SNAP"
export A1_UNIFIED_CORE_DFTB3_CODE_ROOT="$SNAP"
export A1_UNIFIED_CORE_DFTB3_ATTEMPT="$SLURM_JOB_ID"
bash "$SNAP/scripts/run_nylc_a1_unified_core_dftb3_preflight.sh"
