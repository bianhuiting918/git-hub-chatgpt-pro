#!/bin/bash
#SBATCH -J nylc_a1_pt2_frames
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 01:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.err

set -euo pipefail

: "${A1_PT2_FRAME_CODE_SOURCE:?submit with A1_PT2_FRAME_CODE_SOURCE}"
: "${A1_PT2_FRAME_GITHUB_COMMIT:?submit with A1_PT2_FRAME_GITHUB_COMMIT}"

CODE_SOURCE="$(realpath "$A1_PT2_FRAME_CODE_SOURCE")"
[[ -d "$CODE_SOURCE/scripts" ]]
[[ -s "$CODE_SOURCE/GITHUB_COMMIT" ]]
[[ -s "$CODE_SOURCE/SNAPSHOT_SHA256.tsv" ]]
[[ "$(tr -d '[:space:]' <"$CODE_SOURCE/GITHUB_COMMIT")" == \
    "$A1_PT2_FRAME_GITHUB_COMMIT" ]] || {
    printf 'deployed commit marker differs from A1_PT2_FRAME_GITHUB_COMMIT\n' >&2
    exit 2
}
grep -Fq \
    'scripts/extract_audit_nylc_a1_pt2_frames.py' \
    "$CODE_SOURCE/SNAPSHOT_SHA256.tsv"
grep -Fq \
    'scripts/run_extract_audit_nylc_a1_pt2_frames.sh' \
    "$CODE_SOURCE/SNAPSHOT_SHA256.tsv"
(
    cd "$CODE_SOURCE"
    sha256sum -c SNAPSHOT_SHA256.tsv
)

export A1_PT2_FRAME_CODE_ROOT="$CODE_SOURCE"
bash "$CODE_SOURCE/scripts/run_extract_audit_nylc_a1_pt2_frames.sh"
