#!/usr/bin/env bash
set -uo pipefail

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT="${A1_PT2_FRAME_CODE_ROOT:?set A1_PT2_FRAME_CODE_ROOT}"
GITHUB_COMMIT="${A1_PT2_FRAME_GITHUB_COMMIT:?set A1_PT2_FRAME_GITHUB_COMMIT}"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
GMX=/public/software/apps/gromacs/2022.2/hpcx-gcc7.3.1/bin/gmx_mpi
ANALYZER="$CODE_ROOT/scripts/extract_audit_nylc_a1_pt2_frames.py"
ATTEMPT="${SLURM_JOB_ID:-manual_$(date -u '+%Y%m%dT%H%M%SZ')}"
OUT="$TASK_ROOT/a1_activated_nac_20260726/pt2_preorganized_frame_extraction/attempt_$ATTEMPT"

[[ -x "$PY" && -x "$GMX" && -s "$ANALYZER" ]] || {
    printf 'missing Python, GROMACS, or analyzer\n' >&2
    exit 2
}
[[ ! -e "$OUT" ]] || {
    printf 'refusing to overwrite %s\n' "$OUT" >&2
    exit 2
}
mkdir -p "$OUT"

EVENT=nylc_a1_pt2_frame_extraction
COMMAND=run_extract_audit_nylc_a1_pt2_frames.sh
STATE=STARTED
DETAIL="attempt=$ATTEMPT;output=$OUT;scientific_status=NOT_EVALUATED"

append_history() {
    local now
    now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        "$PY" - "$TASK_ROOT/run_history.tsv" "$TASK_ROOT/run_history.jsonl" \
            "$now" "$EVENT" "$COMMAND" "$STATE" "$GITHUB_COMMIT" \
            "$DETAIL" "$OUT" <<'PY'
import json
import os
import sys

tsv_path, jsonl_path, now, event, command, state, commit, detail, output = sys.argv[1:]
tsv_row = "\t".join((now, event, command, state, commit, detail)) + "\n"
jsonl_row = json.dumps({
    "time": now,
    "event": event,
    "command": command,
    "state": state,
    "git_commit": commit,
    "detail": detail,
    "output": output,
}, sort_keys=True) + "\n"
with open(tsv_path, "a+", encoding="utf-8") as tsv, open(
    jsonl_path, "a+", encoding="utf-8"
) as jsonl:
    tsv.seek(0, os.SEEK_END)
    jsonl.seek(0, os.SEEK_END)
    tsv_offset = tsv.tell()
    jsonl_offset = jsonl.tell()
    try:
        tsv.write(tsv_row)
        tsv.flush()
        os.fsync(tsv.fileno())
        jsonl.write(jsonl_row)
        jsonl.flush()
        os.fsync(jsonl.fileno())
    except BaseException:
        tsv.seek(tsv_offset)
        tsv.truncate()
        tsv.flush()
        os.fsync(tsv.fileno())
        jsonl.seek(jsonl_offset)
        jsonl.truncate()
        jsonl.flush()
        os.fsync(jsonl.fileno())
        raise
PY
    ) 9>"$TASK_ROOT/.run_history.lock"
}

append_history

module purge
module load gromacs/2022.2-hpcx-gcc-7.3.1
export GMX_MAXBACKUP=-1 OMP_NUM_THREADS=1

run_candidate() {
    local candidate=$1
    local time_ps=$2
    local tpr=$3
    local xtc=$4
    local tpr_sha=$5
    local xtc_sha=$6
    local candidate_dir="$OUT/$candidate"

    "$PY" "$ANALYZER" \
        --candidate "$candidate" \
        --candidate-dir "$candidate_dir" \
        --prepare || return 1

    if ! {
        printf '%s  %s\n' "$tpr_sha" "$tpr"
        printf '%s  %s\n' "$xtc_sha" "$xtc"
    } | sha256sum -c -; then
        "$PY" "$ANALYZER" \
            --candidate "$candidate" \
            --candidate-dir "$candidate_dir" \
            --fail --stage source_hashes \
            --detail "runner source SHA256 verification failed" || true
        return 1
    fi

    if ! printf 'System\n' | "$GMX" trjconv \
        -s "$tpr" \
        -f "$xtc" \
        -o "$candidate_dir/source.tmp.gro" \
        -dump "$time_ps" \
        -pbc mol \
        -ur compact; then
        "$PY" "$ANALYZER" \
            --candidate "$candidate" \
            --candidate-dir "$candidate_dir" \
            --fail --stage extraction \
            --detail "gmx trjconv -dump failed" || true
        return 1
    fi

    if ! "$PY" "$ANALYZER" \
        --candidate "$candidate" \
        --candidate-dir "$candidate_dir" \
        --github-commit "$GITHUB_COMMIT" \
        --audit; then
        return 1
    fi

    [[ -s "$candidate_dir/source.gro" ]]
    [[ -s "$candidate_dir/manifest.json" ]]
    [[ -s "$candidate_dir/PASS.json" ]]
    [[ -s "$candidate_dir/SHA256.tsv" ]]
}

A_TPR="$TASK_ROOT/a1_activated_nac_20260726/equilibration/nac_evt25_time1462ps/seed26723/attempt_61970146_4_61970151/npt300free/run.tpr"
A_XTC="$TASK_ROOT/a1_activated_nac_20260726/equilibration/nac_evt25_time1462ps/seed26723/attempt_61970146_4_61970151/npt300free/run.xtc"
B_TPR="$TASK_ROOT/a1_activated_nac_20260726/equilibration/nac_evt25_time1462ps/seed26737/attempt_61970146_5_61970146/npt300free/run.tpr"
B_XTC="$TASK_ROOT/a1_activated_nac_20260726/equilibration/nac_evt25_time1462ps/seed26737/attempt_61970146_5_61970146/npt300free/run.xtc"
A_TPR_SHA=c60078a92c2ace51facde4ef64e453f690177fc88b4d6363427f935944fa2e43
A_XTC_SHA=1a54f1b5b9f139b746c22d9e0f7e9a4a94eb8154bf2b881888986eedca933d89
B_TPR_SHA=dbd19a399547319d10630430ed494d33f6271bab0f30cb0de5a466c6af13ba20
B_XTC_SHA=fcba14da98b331368061dcd990f2467628ad77b9b7a9c4ce88090f92e0831b05

passed=0
failed=0
if run_candidate \
    seed26723_t378_f189 378.0 \
    "$A_TPR" "$A_XTC" "$A_TPR_SHA" "$A_XTC_SHA"; then
    passed=$((passed + 1))
else
    failed=$((failed + 1))
fi
if run_candidate \
    seed26737_t676_f338 676.0 \
    "$B_TPR" "$B_XTC" "$B_TPR_SHA" "$B_XTC_SHA"; then
    passed=$((passed + 1))
else
    failed=$((failed + 1))
fi

if ((failed == 0 && passed == 2)); then
    STATE=PASS
    DETAIL="attempt=$ATTEMPT;passed=2;failed=0;gate=PASS_A1_PT2_DUAL_FRAME_EXTRACTION_AUDIT;scope=preorganization_only;output=$OUT"
    append_history
    printf 'PASS_A1_PT2_DUAL_FRAME_EXTRACTION_AUDIT %s\n' "$OUT"
    exit 0
fi

STATE=NOT_EVALUATED
DETAIL="attempt=$ATTEMPT;passed=$passed;failed=$failed;gate=NOT_EVALUATED_A1_PT2_DUAL_FRAME_EXTRACTION;output=$OUT"
append_history
printf 'NOT_EVALUATED_A1_PT2_DUAL_FRAME_EXTRACTION passed=%d failed=%d\n' \
    "$passed" "$failed" >&2
exit 1
