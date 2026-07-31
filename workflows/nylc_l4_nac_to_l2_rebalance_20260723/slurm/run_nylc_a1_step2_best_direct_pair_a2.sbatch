#!/bin/bash
#SBATCH -J nylc_a2pair
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -c 1
#SBATCH --mem-per-cpu=2500M
#SBATCH --time=08:00:00
#SBATCH --array=0-1
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.err

set -euo pipefail
: "${A1_STEP2_DIRECT_A2_CODE_SOURCE:?set immutable workflow code source}"
: "${A1_STEP2_DIRECT_A2_GITHUB_COMMIT:?set full runtime GitHub commit}"

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT="$(realpath "$A1_STEP2_DIRECT_A2_CODE_SOURCE")"
COMMIT="$A1_STEP2_DIRECT_A2_GITHUB_COMMIT"
JOB_ID="${SLURM_JOB_ID:?run through sbatch}"
ARRAY_JOB_ID="${SLURM_ARRAY_JOB_ID:?run as array}"
INDEX="${SLURM_ARRAY_TASK_ID:?run as array task}"
SOURCE_TASKS=(4 6)
SOURCE_SHAS=(
  cbe9fbddb0ea8c3a6f12c9ae03d19864b6f4d6e2582901fe35ffa828d0f0c1ca
  5c4227dc8752f1a181e5b17d2890e6bce5a8a790f3ed5f91a181445e72006efe
)
SOURCE_TASK="${SOURCE_TASKS[$INDEX]}"
SOURCE_SHA="${SOURCE_SHAS[$INDEX]}"
OUT_ROOT="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step2_best_direct_pair_a2_preflight"
OUT="$OUT_ROOT/attempt_${ARRAY_JOB_ID}_${INDEX}"
SCRATCH="${SLURM_TMPDIR:-/tmp}/nylc_a2_best_direct_${ARRAY_JOB_ID}_${INDEX}"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
DRIVER="$CODE_ROOT/scripts/prepare_audit_nylc_a1_step2_best_direct_pair_a2.py"
BASE_DRIVER="$CODE_ROOT/scripts/prepare_audit_nylc_a1_step2_qmwater_endpoint.py"
PRMTOP="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285/prepared/system.prmtop"
SOURCE="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step2_water_network_nearmiss_continuation/attempt_62507122_${SOURCE_TASK}/direct_event_0.rst7"
TERMINAL=0

write_hashes() {
    [[ -d "$OUT" ]] || return 0
    (
        cd "$OUT"
        for name in A2_MANIFEST.json READY.json WATER_SELECTION.json \
            WATER_CANDIDATES.tsv A2_LEG_0.json A2_LEG_1.json RESULT.json \
            PASS.json FAIL.json NOT_EVALUATED.json \
            a2_leg0_endpoint.rst7 a2_leg1_endpoint.rst7 \
            A2_LEG_0.stage.out A2_LEG_0.stage.mdinfo \
            A2_LEG_0.stage.rst7 A2_LEG_0.a2.mdcrd A2_LEG_0.engine.rc \
            A2_LEG_1.stage.out A2_LEG_1.stage.mdinfo \
            A2_LEG_1.stage.rst7 A2_LEG_1.a2.mdcrd A2_LEG_1.engine.rc
        do
            [[ -f "$name" ]] && sha256sum "$name"
        done | sort -k2
    ) >"$OUT/SHA256.tsv"
}

append_history() {
    local state="$1" detail="$2" now
    now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$now" "nylc_a1_step2_best_direct_pair_a2" \
            "run_nylc_a1_step2_best_direct_pair_a2.sbatch" "$state" \
            "$COMMIT" "$detail" >>"$TASK_ROOT/run_history.tsv"
        "$PY" - "$TASK_ROOT/run_history.jsonl" "$now" "$state" "$detail" \
            "$COMMIT" "$OUT" <<'PY'
import json,sys
path,now,state,detail,commit,out=sys.argv[1:]
with open(path,"a",encoding="utf-8") as handle:
    handle.write(json.dumps({
        "time":now,
        "event":"nylc_a1_step2_best_direct_pair_a2",
        "command":"run_nylc_a1_step2_best_direct_pair_a2.sbatch",
        "state":state,
        "git_commit":commit,
        "detail":detail,
        "output":out,
    },sort_keys=True)+"\n")
PY
    ) 9>>"$TASK_ROOT/.run_history.lock"
}

finish() {
    local rc=$?
    trap - EXIT
    if (( rc != 0 && ! TERMINAL )); then
        mkdir -p "$OUT"
        "$PY" - "$OUT" "$rc" "$SOURCE_TASK" <<'PY' || true
import json,pathlib,sys
out=pathlib.Path(sys.argv[1]); rc=int(sys.argv[2]); source_task=int(sys.argv[3])
payload={
 "schema_version":1,
 "status":"NOT_EVALUATED_STEP2_A2_RUNTIME_FAILURE",
 "classification":"NOT_EVALUATED_STEP2_A2_RUNTIME_FAILURE",
 "seed":"seed26737",
 "source_task":source_task,
 "exit_code":rc,
 "NEXT":"STEP2_PRODUCT_ENDPOINT_BLOCKED_PENDING_A2_PASS",
 "product_endpoint_implemented":False,
}
for name in ("RESULT.json","NOT_EVALUATED.json"):
    (out/name).write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
PY
        write_hashes || true
        append_history NOT_EVALUATED \
            "job=$JOB_ID;source_sha=$SOURCE_SHA;exit_code=$rc;output=$OUT" || true
    fi
    exit "$rc"
}
trap finish EXIT

[[ "$COMMIT" =~ ^[0-9a-f]{40}$ ]]
test ! -e "$OUT"
test -s "$CODE_ROOT/GITHUB_COMMIT"
test -s "$CODE_ROOT/SNAPSHOT_SHA256.tsv"
test "$(tr -d '[:space:]' <"$CODE_ROOT/GITHUB_COMMIT")" = "$COMMIT"
(
    cd "$CODE_ROOT"
    sha256sum -c SNAPSHOT_SHA256.tsv
)
for file in "$DRIVER" "$BASE_DRIVER" "$PRMTOP" "$SOURCE"; do test -s "$file"; done
test "$(sha256sum "$SOURCE" | awk '{print $1}')" = "$SOURCE_SHA"
mkdir -p "$OUT_ROOT"

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$TASK_ROOT/../nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime"
export OMP_NUM_THREADS=1
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"

append_history STARTED \
    "job=$JOB_ID;seed=seed26737;source_task=$SOURCE_TASK;source_sha=$SOURCE_SHA;water=13046,13047,13048;qm=149;q=0;e=518;link=6;output=$OUT"

"$PY" "$DRIVER" --mode prepare --seed-index 0 --output "$OUT" \
    --scratch "$SCRATCH" --code-root "$CODE_ROOT" --github-commit "$COMMIT"

if [[ ! -s "$OUT/READY.json" ]]; then
    test -s "$OUT/RESULT.json"
    write_hashes
    append_history NOT_EVALUATED \
        "job=$JOB_ID;seed=seed26737;source_task=$SOURCE_TASK;prepare_terminal=1;output=$OUT"
    TERMINAL=1
    trap - EXIT
    exit 0
fi

declare -a PIDS=()
for LEG in 0 1; do
    STAGE="$SCRATCH/a2_leg$LEG"
    test -s "$STAGE/stage.in"
    (
        cd "$STAGE"
        set +e
        mpirun --bind-to none -np 8 sander.MPI -O \
            -i stage.in -o stage.out -p "$PRMTOP" -c "$SOURCE" \
            -r stage.rst7 -x a2.mdcrd -inf stage.mdinfo
        printf '%s\n' "$?" >engine.rc
        exit 0
    ) &
    PIDS[$LEG]=$!
done
for LEG in 0 1; do wait "${PIDS[$LEG]}"; done

# Persist raw evidence before auditing so future parser changes are reauditable.
for LEG in 0 1; do
    STAGE="$SCRATCH/a2_leg$LEG"
    for NAME in stage.out stage.mdinfo stage.rst7 a2.mdcrd engine.rc; do
        test -s "$STAGE/$NAME"
        cp "$STAGE/$NAME" "$OUT/A2_LEG_${LEG}.$NAME"
    done
done

for LEG in 0 1; do
    "$PY" "$DRIVER" --mode audit-leg --output "$OUT" \
        --scratch "$SCRATCH" --leg "$LEG"
    if "$PY" - "$OUT/A2_LEG_$LEG.json" <<'PY'
import json,sys
raise SystemExit(0 if json.load(open(sys.argv[1]))["technical_pass"] else 1)
PY
    then
        cp "$SCRATCH/a2_leg$LEG/stage.rst7" \
            "$OUT/a2_leg${LEG}_endpoint.rst7"
    fi
done

"$PY" "$DRIVER" --mode finalize --output "$OUT"
test -s "$OUT/RESULT.json"
write_hashes
STATUS="$("$PY" - "$OUT/RESULT.json" <<'PY'
import json,sys
print(json.load(open(sys.argv[1]))["status"])
PY
)"
append_history COMPLETED \
    "job=$JOB_ID;seed=seed26737;source_task=$SOURCE_TASK;status=$STATUS;legs=2x8;output=$OUT"
TERMINAL=1
trap - EXIT
printf '%s output=%s\n' "$STATUS" "$OUT"
