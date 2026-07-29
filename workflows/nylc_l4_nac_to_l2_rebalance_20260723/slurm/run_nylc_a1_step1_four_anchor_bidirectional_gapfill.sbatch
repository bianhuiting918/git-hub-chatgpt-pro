#!/bin/bash
#SBATCH -J nylc_a1fa8
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem-per-cpu=2500M
#SBATCH --time=24:00:00
#SBATCH --array=0-15
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.err
set -euo pipefail
: "${A1_FOUR_ANCHOR_CODE_SOURCE:?set immutable code source}"
: "${A1_FOUR_ANCHOR_GITHUB_COMMIT:?set immutable GitHub commit}"

CODE_SOURCE="$(realpath "$A1_FOUR_ANCHOR_CODE_SOURCE")"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
test -s "$CODE_SOURCE/GITHUB_COMMIT"
test -s "$CODE_SOURCE/SNAPSHOT_SHA256.tsv"
test "$(tr -d '[:space:]' < "$CODE_SOURCE/GITHUB_COMMIT")" = "$A1_FOUR_ANCHOR_GITHUB_COMMIT"

for file in  scripts/prepare_audit_nylc_a1_step1_four_anchor_bidirectional_gapfill.py  scripts/run_nylc_a1_step1_four_anchor_bidirectional_gapfill.sh  slurm/run_nylc_a1_step1_four_anchor_bidirectional_gapfill.sbatch  scripts/prepare_audit_nylc_a1_step1_product_reverse_boundary.py  scripts/prepare_audit_nylc_a1_step1_acyl_endpoint_stability.py  scripts/prepare_audit_nylc_a1_step1_tetra_parallel_scout.py; do
    grep -Fq "$file" "$CODE_SOURCE/SNAPSHOT_SHA256.tsv"
    test -s "$CODE_SOURCE/$file"
done
(
    cd "$CODE_SOURCE"
    sha256sum -c SNAPSHOT_SHA256.tsv
)

DESCRIPTION=$("$PY" "$CODE_SOURCE/scripts/prepare_audit_nylc_a1_step1_four_anchor_bidirectional_gapfill.py" --mode describe)
"$PY" - "$DESCRIPTION" <<'PY'
import json,sys
x=json.loads(sys.argv[1])
assert x["anchors"] == ["R","I1","I2","P"]
assert x["directions"] == ["REACTANT","PRODUCT"]
assert x["seed_count"] == 2
assert x["array_task_count"] == 16
assert x["stage_restart_persisted"] is True
assert x["restrained_structures_are_not_ts"] is True
PY
bash -n "$CODE_SOURCE/scripts/run_nylc_a1_step1_four_anchor_bidirectional_gapfill.sh"
bash -n "$CODE_SOURCE/slurm/run_nylc_a1_step1_four_anchor_bidirectional_gapfill.sbatch"

export A1_FOUR_ANCHOR_CODE_ROOT="$CODE_SOURCE"
bash "$CODE_SOURCE/scripts/run_nylc_a1_step1_four_anchor_bidirectional_gapfill.sh"
