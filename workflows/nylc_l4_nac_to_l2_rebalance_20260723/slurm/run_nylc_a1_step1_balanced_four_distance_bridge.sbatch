#!/bin/bash
#SBATCH -J nylc_a1b4d
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem-per-cpu=2500M
#SBATCH --time=12:00:00
#SBATCH --array=0-5
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.err
set -euo pipefail
: "${A1_BALANCED_BRIDGE_CODE_SOURCE:?set A1_BALANCED_BRIDGE_CODE_SOURCE}"
: "${A1_BALANCED_BRIDGE_GITHUB_COMMIT:?set A1_BALANCED_BRIDGE_GITHUB_COMMIT}"
CODE_SOURCE="$(realpath "$A1_BALANCED_BRIDGE_CODE_SOURCE")"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
test -d "$CODE_SOURCE/scripts"; test -d "$CODE_SOURCE/slurm"
test -s "$CODE_SOURCE/GITHUB_COMMIT"; test -s "$CODE_SOURCE/SNAPSHOT_SHA256.tsv"
test "$(tr -d '[:space:]' <"$CODE_SOURCE/GITHUB_COMMIT")" = "$A1_BALANCED_BRIDGE_GITHUB_COMMIT"
for file in scripts/prepare_audit_nylc_a1_step1_balanced_four_distance_bridge.py scripts/run_nylc_a1_step1_balanced_four_distance_bridge.sh slurm/run_nylc_a1_step1_balanced_four_distance_bridge.sbatch scripts/prepare_audit_nylc_a1_step1_product_reverse_boundary.py scripts/prepare_audit_nylc_a1_step1_acyl_endpoint_stability.py scripts/prepare_audit_nylc_a1_step1_tetra_parallel_scout.py audit/nylc_a1_acyl_endpoint_connectivity_sensitivity_62216380_v1.json; do
 grep -Fq "$file" "$CODE_SOURCE/SNAPSHOT_SHA256.tsv"; test -s "$CODE_SOURCE/$file"
done
( cd "$CODE_SOURCE"; sha256sum -c SNAPSHOT_SHA256.tsv )
export PYTHONDONTWRITEBYTECODE=1
DESCRIPTION="$("$PY" "$CODE_SOURCE/scripts/prepare_audit_nylc_a1_step1_balanced_four_distance_bridge.py" --mode describe)"
"$PY" - "$DESCRIPTION" <<'PY'
import json,sys
x=json.loads(sys.argv[1])
assert x["seed_count"] == 2 and x["branch_count"] == 3 and x["array_task_count"] == 6
assert x["candidate_limit_per_seed"] == 1 and x["continuous_frame_requirement"] == 3
assert x["bias_representation"] == "four_independent_distance_restraints_not_a_linear_combination_cv"
assert x["restrained_endpoint_shooting_forbidden"] is True
PY
bash -n "$CODE_SOURCE/scripts/run_nylc_a1_step1_balanced_four_distance_bridge.sh"
bash -n "$CODE_SOURCE/slurm/run_nylc_a1_step1_balanced_four_distance_bridge.sbatch"
export A1_BALANCED_BRIDGE_CODE_ROOT="$CODE_SOURCE"
bash "$CODE_SOURCE/scripts/run_nylc_a1_step1_balanced_four_distance_bridge.sh"
