#!/bin/bash
#SBATCH -J nylc_a1_rep354
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.err

set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
: "\${A1_GITHUB_COMMIT:?submit with A1_GITHUB_COMMIT set to the immutable GitHub commit SHA}"
SNAP="$TASK_ROOT/a1_activated_nac_20260726/code_snapshots/representative_frame_\${SLURM_JOB_ID}"

[[ ! -e "$SNAP" ]] || {
    printf 'refusing to overwrite code snapshot %s\n' "$SNAP" >&2
    exit 2
}
mkdir -p "$SNAP/scripts" "$SNAP/mdp"
for name in \
    run_nylc_a1_representative_frame.sh \
    audit_nylc_a1_representative_frame.py \
    generate_nylc_m1_ensemble_primitives.py; do
    cp "$FLOW/scripts/$name" "$SNAP/scripts/$name"
done
cp "$FLOW/mdp/em_cg_flexible_m1.mdp" "$SNAP/mdp/em_cg_flexible_m1.mdp"
sha256sum "$SNAP"/scripts/* "$SNAP"/mdp/* >"$SNAP/SNAPSHOT_SHA256.tsv"
chmod -R a-w "$SNAP"

export A1_REP_CODE_ROOT="$SNAP"
export A1_REP_ATTEMPT="$SLURM_JOB_ID"
bash "$SNAP/scripts/run_nylc_a1_representative_frame.sh"
