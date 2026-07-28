#!/bin/bash
#SBATCH -J nylc_a1revb
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem-per-cpu=2500M
#SBATCH --time=12:00:00
#SBATCH --array=0-1
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.err

set -euo pipefail

: "\${A1_PRODUCT_REVERSE_CODE_SOURCE:?set A1_PRODUCT_REVERSE_CODE_SOURCE}"
: "\${A1_PRODUCT_REVERSE_GITHUB_COMMIT:?set A1_PRODUCT_REVERSE_GITHUB_COMMIT}"

CODE_SOURCE="$(realpath "$A1_PRODUCT_REVERSE_CODE_SOURCE")"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python

test -d "$CODE_SOURCE/scripts"
test -d "$CODE_SOURCE/slurm"
test -s "$CODE_SOURCE/GITHUB_COMMIT"
test -s "$CODE_SOURCE/SNAPSHOT_SHA256.tsv"
test "$(tr -d '[:space:]' <"$CODE_SOURCE/GITHUB_COMMIT")" = \
    "$A1_PRODUCT_REVERSE_GITHUB_COMMIT"

for file in \
    scripts/prepare_audit_nylc_a1_step1_product_reverse_boundary.py \
    scripts/run_nylc_a1_step1_product_reverse_boundary.sh \
    slurm/run_nylc_a1_step1_product_reverse_boundary.sbatch \
    scripts/prepare_audit_nylc_a1_step1_acyl_endpoint_stability.py \
    scripts/prepare_audit_nylc_a1_step1_tetra_parallel_scout.py \
    audit/nylc_a1_acyl_endpoint_connectivity_sensitivity_62216380_v1.json
do
    grep -Fq "$file" "$CODE_SOURCE/SNAPSHOT_SHA256.tsv"
    test -s "$CODE_SOURCE/$file"
done

(
    cd "$CODE_SOURCE"
    sha256sum -c SNAPSHOT_SHA256.tsv
)

export PYTHONDONTWRITEBYTECODE=1
DESCRIPTION="$(
    "$PY" "$CODE_SOURCE/scripts/prepare_audit_nylc_a1_step1_product_reverse_boundary.py" \
        --mode describe
)"
"$PY" - "$DESCRIPTION" <<'PY'
import json
import sys

description = json.loads(sys.argv[1])
assert description["source_authority_commit"] == (
    "19f327b422fd5081b119e4ff4c1883edefd85ba4"
)
assert description["source_job"] == "62216380"
assert description["seed_count"] == 2
assert description["window_count_per_seed"] == 4
assert description["candidate_limit_per_seed"] == 2
assert description["release_trajectory_policy"] == "scratch_only"
assert description["restrained_endpoint_shooting_forbidden"] is True
PY

bash -n "$CODE_SOURCE/scripts/run_nylc_a1_step1_product_reverse_boundary.sh"
bash -n "$CODE_SOURCE/slurm/run_nylc_a1_step1_product_reverse_boundary.sbatch"

export A1_PRODUCT_REVERSE_CODE_ROOT="$CODE_SOURCE"
bash "$CODE_SOURCE/scripts/run_nylc_a1_step1_product_reverse_boundary.sh"
