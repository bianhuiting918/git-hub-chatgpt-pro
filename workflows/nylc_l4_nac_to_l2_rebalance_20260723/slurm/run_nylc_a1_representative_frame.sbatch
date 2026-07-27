#!/bin/bash
#SBATCH -J nylc_a1_rep354
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 00:30:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.err

set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
: "${A1_GITHUB_COMMIT:?submit with A1_GITHUB_COMMIT set to the immutable GitHub commit SHA}"
: "${A1_REP_CODE_SOURCE:?submit with A1_REP_CODE_SOURCE set to a deployed immutable code root}"
CODE_SOURCE="$(realpath "$A1_REP_CODE_SOURCE")"
[[ -d "$CODE_SOURCE/scripts" && -d "$CODE_SOURCE/mdp" ]]
[[ -s "$CODE_SOURCE/GITHUB_COMMIT" ]]
[[ -s "$CODE_SOURCE/SNAPSHOT_SHA256.tsv" ]]
[[ "$(tr -d '[:space:]' <"$CODE_SOURCE/GITHUB_COMMIT")" == "$A1_GITHUB_COMMIT" ]] || {
    printf 'deployed code commit marker does not match A1_GITHUB_COMMIT\n' >&2
    exit 2
}
(
    cd "$CODE_SOURCE"
    sha256sum -c SNAPSHOT_SHA256.tsv
)
SNAP="$TASK_ROOT/a1_activated_nac_20260726/code_snapshots/representative_frame_${SLURM_JOB_ID}"

[[ ! -e "$SNAP" ]] || {
    printf 'refusing to overwrite code snapshot %s\n' "$SNAP" >&2
    exit 2
}
mkdir -p "$SNAP/scripts" "$SNAP/mdp"
for name in \
    run_nylc_a1_representative_frame.sh \
    audit_nylc_a1_representative_frame.py \
    generate_nylc_m1_ensemble_primitives.py; do
    cp "$CODE_SOURCE/scripts/$name" "$SNAP/scripts/$name"
done
cp "$CODE_SOURCE/mdp/em_cg_flexible_m1.mdp" "$SNAP/mdp/em_cg_flexible_m1.mdp"
cp "$CODE_SOURCE/GITHUB_COMMIT" "$SNAP/GITHUB_COMMIT"
(
    cd "$SNAP"
    sha256sum scripts/* mdp/* GITHUB_COMMIT >SNAPSHOT_SHA256.tsv
)
chmod -R a-w "$SNAP"

export A1_REP_CODE_ROOT="$SNAP"
export A1_REP_ATTEMPT="$SLURM_JOB_ID"
bash "$SNAP/scripts/run_nylc_a1_representative_frame.sh"
