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
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
: "${A1_GITHUB_COMMIT:?submit with A1_GITHUB_COMMIT set to the immutable GitHub commit SHA}"
: "${A1_REP_CODE_SOURCE:?submit with A1_REP_CODE_SOURCE set to a deployed immutable code root}"
CODE_SOURCE="$(realpath "$A1_REP_CODE_SOURCE")"
[[ -d "$CODE_SOURCE/scripts" && -d "$CODE_SOURCE/mdp" ]]
[[ -s "$CODE_SOURCE/GITHUB_COMMIT" ]]
[[ -s "$CODE_SOURCE/SNAPSHOT_SHA256.tsv" ]]
dependency=analyze_nylc_m1_proton_geometry.py
grep -Fq "scripts/$dependency" "$CODE_SOURCE/SNAPSHOT_SHA256.tsv" || {
    printf 'deployed source manifest omits %s\n' "$dependency" >&2
    exit 2
}
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
    generate_nylc_m1_ensemble_primitives.py \
    analyze_nylc_m1_proton_geometry.py; do
    cp "$CODE_SOURCE/scripts/$name" "$SNAP/scripts/$name"
done
cp "$CODE_SOURCE/mdp/em_cg_flexible_m1.mdp" "$SNAP/mdp/em_cg_flexible_m1.mdp"
cp "$CODE_SOURCE/GITHUB_COMMIT" "$SNAP/GITHUB_COMMIT"
(
    cd "$SNAP"
    sha256sum scripts/* mdp/* GITHUB_COMMIT >SNAPSHOT_SHA256.tsv
)
PYTHONPATH="$SNAP/scripts" "$PY" - \
    "$SNAP/scripts/audit_nylc_a1_representative_frame.py" <<'PY'
import importlib.util
import pathlib
import sys
source = pathlib.Path(sys.argv[1])
spec = importlib.util.spec_from_file_location("_a1_rep_snapshot_smoke", source)
if spec is None or spec.loader is None:
    raise SystemExit("cannot load snapshotted representative-frame auditor")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
if not callable(module._load_existing_gate_opening()):
    raise SystemExit("snapshotted gate-opening primitive is not callable")
PY
chmod -R a-w "$SNAP"

export A1_REP_CODE_ROOT="$SNAP"
export A1_REP_ATTEMPT="$SLURM_JOB_ID"
bash "$SNAP/scripts/run_nylc_a1_representative_frame.sh"
