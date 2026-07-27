#!/usr/bin/env bash
set -euo pipefail

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
CODE_ROOT="\${A1_REP_CODE_ROOT:-$FLOW}"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
GMX=/public/software/apps/gromacs/2022.2/hpcx-gcc7.3.1/bin/gmx_mpi
GITHUB_COMMIT="\${A1_GITHUB_COMMIT:-unknown}"

SOURCE_ROOT="$TASK_ROOT/a1_activated_nac_20260726/equilibration/nac_evt25_time1462ps/seed26723/attempt_61970146_4_61970151/npt300free"
SOURCE_TPR="$SOURCE_ROOT/run.tpr"
SOURCE_XTC="$SOURCE_ROOT/run.xtc"
TOPOLOGY_ROOT="$TASK_ROOT/a1_activated_nac_20260726/em/attempt_61968026_1_61968026/nac_evt25_time1462ps/input"
SOURCE_NDX="$TASK_ROOT/ensemble/candidates/nac_evt25_time1462ps/build_job_61813799_11/build/source_cycle.ndx"
SOURCE_MDP="$CODE_ROOT/mdp/em_cg_flexible_m1.mdp"
EXPECTED_TPR_SHA256=c60078a92c2ace51facde4ef64e453f690177fc88b4d6363427f935944fa2e43
EXPECTED_XTC_SHA256=1a54f1b5b9f139b746c22d9e0f7e9a4a94eb8154bf2b881888986eedca933d89
ATTEMPT="\${A1_REP_ATTEMPT:-\${SLURM_JOB_ID:-manual_$(date -u '+%Y%m%dT%H%M%SZ')}}"
OUT="$TASK_ROOT/a1_activated_nac_20260726/representative_frame/attempt_$ATTEMPT"

[[ ! -e "$OUT" ]] || {
    printf 'refusing to overwrite %s\n' "$OUT" >&2
    exit 2
}
mkdir -p "$OUT"

EVENT=nylc_a1_representative_frame
COMMAND="run_nylc_a1_representative_frame.sh"
STATE=STARTED
CURRENT=initialization
DETAIL="source=evt25;seed=26723;time_ps=354;scientific_status=NOT_EVALUATED;output=$OUT"

append_history() {
    local now
    now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$now" "$EVENT" "$COMMAND" "$STATE" "$GITHUB_COMMIT" "$DETAIL" \
            >>"$TASK_ROOT/run_history.tsv"
        "$PY" - "$now" "$EVENT" "$COMMAND" "$STATE" "$GITHUB_COMMIT" \
            "$DETAIL" "$OUT" >>"$TASK_ROOT/run_history.jsonl" <<'PY'
import json
import sys
now, event, command, state, commit, detail, output = sys.argv[1:]
print(json.dumps({
    "time": now,
    "event": event,
    "command": command,
    "state": state,
    "git_commit": commit,
    "detail": detail,
    "output": output,
}, sort_keys=True))
PY
    ) 9>"$TASK_ROOT/.run_history.lock"
}

finish() {
    local code=$?
    trap - EXIT
    STATE=NOT_EVALUATED
    DETAIL="source=evt25;seed=26723;stage=$CURRENT;exit_code=$code;scientific_status=NOT_EVALUATED;output=$OUT"
    "$PY" - "$OUT/NOT_EVALUATED.json" "$CURRENT" "$code" \
        "$OUT/A1_REPRESENTATIVE_FRAME_AUDIT.json" <<'PY'
import json
import pathlib
import sys
output, stage, code, audit_path = sys.argv[1:]
audit_status = None
path = pathlib.Path(audit_path)
if path.is_file():
    try:
        audit_status = json.loads(path.read_text(encoding="utf-8")).get(
            "scientific_status"
        )
    except (OSError, ValueError):
        audit_status = None
payload = {
    "schema_version": 1,
    "status": "NOT_EVALUATED_A1_REPRESENTATIVE_FRAME",
    "technical_status": "FAIL",
    "scientific_status": "NOT_EVALUATED",
    "audit_scientific_status": audit_status,
    "failed_stage": stage,
    "exit_code": int(code),
    "promoted": False,
}
pathlib.Path(output).write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
PY
    append_history
    exit "$code"
}
trap finish EXIT
append_history

CURRENT=input_validation
for path in \
    "$SOURCE_TPR" \
    "$SOURCE_XTC" \
    "$SOURCE_NDX" \
    "$SOURCE_MDP" \
    "$TOPOLOGY_ROOT/topol.top"; do
    [[ -s "$path" ]] || {
        printf 'missing required input %s\n' "$path" >&2
        exit 2
    }
done
mapfile -d '' ITP_FILES < <(
    find "$TOPOLOGY_ROOT" -maxdepth 1 -type f -name '*.itp' -print0 | sort -z
)
((\${#ITP_FILES[@]} > 0)) || {
    printf 'no topology ITP files under %s\n' "$TOPOLOGY_ROOT" >&2
    exit 2
}

CURRENT=source_hashes
printf '%s  %s\n' "$EXPECTED_TPR_SHA256" "$SOURCE_TPR" | sha256sum -c -
printf '%s  %s\n' "$EXPECTED_XTC_SHA256" "$SOURCE_XTC" | sha256sum -c -
sha256sum \
    "$SOURCE_TPR" \
    "$SOURCE_XTC" \
    "$SOURCE_NDX" \
    "$SOURCE_MDP" \
    "$TOPOLOGY_ROOT/topol.top" \
    "\${ITP_FILES[@]}" >"$OUT/input_sha256.tsv"

module purge
module load gromacs/2022.2-hpcx-gcc-7.3.1
export GMX_MAXBACKUP=-1 OMP_NUM_THREADS="\${SLURM_CPUS_PER_TASK:-4}"

CURRENT=extract_system_354ps
printf 'System\n' | "$GMX" trjconv \
    -s "$SOURCE_TPR" \
    -f "$SOURCE_XTC" \
    -o "$OUT/source.tmp.gro" \
    -dump 354 \
    >"$OUT/extract.stdout" 2>"$OUT/extract.stderr"
[[ -s "$OUT/source.tmp.gro" ]] || {
    printf 'System extraction did not create source.tmp.gro\n' >&2
    exit 2
}

CURRENT=representative_frame_audit
"$PY" "$CODE_ROOT/scripts/audit_nylc_a1_representative_frame.py" \
    --tpr "$SOURCE_TPR" \
    --xtc "$SOURCE_XTC" \
    --gro "$OUT/source.tmp.gro" \
    --ndx "$SOURCE_NDX" \
    --time-ps 354 \
    --tpr-sha256 "$EXPECTED_TPR_SHA256" \
    --xtc-sha256 "$EXPECTED_XTC_SHA256" \
    --output "$OUT/A1_REPRESENTATIVE_FRAME_AUDIT.json" \
    >"$OUT/audit.stdout" 2>"$OUT/audit.stderr"
"$PY" - "$OUT/A1_REPRESENTATIVE_FRAME_AUDIT.json" <<'PY'
import json
import sys
payload = json.load(open(sys.argv[1], encoding="utf-8"))
if payload.get("scientific_status") != "PASS_A1_REPRESENTATIVE_NAC_FRAME":
    raise SystemExit("representative-frame auditor did not PASS")
if not payload.get("gates") or not all(payload["gates"].values()):
    raise SystemExit("representative-frame audit contains a failed gate")
PY

CURRENT=render_unrestrained_preflight_mdp
"$PY" - "$SOURCE_MDP" "$OUT/preflight.unrestrained.mdp" <<'PY'
import pathlib
import sys
source = pathlib.Path(sys.argv[1])
output = pathlib.Path(sys.argv[2])
kept = []
removed_defines = []
for line in source.read_text(encoding="utf-8").splitlines():
    body = line.split(";", 1)[0].strip()
    if "=" in body:
        key, value = (part.strip() for part in body.split("=", 1))
        if key == "define":
            removed_defines.append(value)
            continue
    kept.append(line)
expected = "-DPOSRES -DPOSRES_L2_1000 -DFLEXIBLE"
if removed_defines != [expected]:
    raise SystemExit(
        "source EM MDP restrained define differs from the frozen expected value"
    )
rendered = "\n".join(kept).rstrip() + "\n"
for line in rendered.splitlines():
    body = line.split(";", 1)[0].strip()
    if "=" in body:
        key, value = (part.strip() for part in body.split("=", 1))
        if key.lower() == "define" and value:
            raise SystemExit("rendered preflight MDP retained a nonempty define")
output.write_text(rendered, encoding="utf-8")
PY

CURRENT=unrestrained_grompp
(
    cd "$TOPOLOGY_ROOT"
    "$GMX" grompp \
        -f "$OUT/preflight.unrestrained.mdp" \
        -c "$OUT/source.tmp.gro" \
        -p "$TOPOLOGY_ROOT/topol.top" \
        -o "$OUT/preflight.tpr" \
        -po "$OUT/preflight.processed.mdp" \
        -maxwarn 0 \
        >"$OUT/grompp.stdout" 2>"$OUT/grompp.stderr"
)
"$GMX" dump -s "$OUT/preflight.tpr" \
    >"$OUT/preflight.tpr.dump" 2>"$OUT/preflight.tpr_dump.stderr"

CURRENT=unrestrained_preflight_gate
"$PY" - \
    "$OUT/preflight.tpr.dump" \
    "$OUT/preflight.processed.mdp" \
    "$OUT/PREFLIGHT_GATE.json" <<'PY'
import json
import pathlib
import re
import sys
dump = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8", errors="replace")
mdp = pathlib.Path(sys.argv[2]).read_text(encoding="utf-8", errors="replace")
position_restraints = sum(
    map(int, re.findall(r"#posres_xA\s*=\s*(\d+)", dump))
)
distance_restraints = len(re.findall(r"(?i)\bDISRES\b", dump))
nonempty_defines = []
for line in mdp.splitlines():
    body = line.split(";", 1)[0].strip()
    if "=" not in body:
        continue
    key, value = (part.strip() for part in body.split("=", 1))
    if key.lower() == "define" and value:
        nonempty_defines.append(value)
passed = not position_restraints and not distance_restraints and not nonempty_defines
payload = {
    "schema_version": 1,
    "status": (
        "PASS_TECHNICAL_UNRESTRAINED_GROMPP"
        if passed
        else "FAIL_TECHNICAL_RESTRAINTS_PRESENT"
    ),
    "grompp_maxwarn": 0,
    "position_restraints": position_restraints,
    "distance_restraints": distance_restraints,
    "nonempty_defines": nonempty_defines,
}
pathlib.Path(sys.argv[3]).write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
if not passed:
    raise SystemExit("preflight TPR/processed MDP contains restraints or defines")
PY

CURRENT=promote_artifacts
mv "$OUT/source.tmp.gro" "$OUT/representative_354ps.gro"
"$GMX" editconf \
    -f "$OUT/representative_354ps.gro" \
    -o "$OUT/representative_354ps.pdb" \
    >"$OUT/pdb.stdout" 2>"$OUT/pdb.stderr"
[[ -s "$OUT/representative_354ps.gro" ]]
[[ -s "$OUT/representative_354ps.pdb" ]]
chmod 0444 "$OUT/representative_354ps.gro" "$OUT/representative_354ps.pdb"

CURRENT=final_hashes
{
    sha256sum \
        "$SOURCE_TPR" \
        "$SOURCE_XTC" \
        "$SOURCE_NDX" \
        "$SOURCE_MDP" \
        "$TOPOLOGY_ROOT/topol.top" \
        "\${ITP_FILES[@]}"
    sha256sum \
        "$OUT/A1_REPRESENTATIVE_FRAME_AUDIT.json" \
        "$OUT/preflight.unrestrained.mdp" \
        "$OUT/preflight.processed.mdp" \
        "$OUT/preflight.tpr" \
        "$OUT/PREFLIGHT_GATE.json" \
        "$OUT/representative_354ps.gro" \
        "$OUT/representative_354ps.pdb"
} >"$OUT/SHA256.tsv"

"$PY" - \
    "$OUT/PASS.json" \
    "$OUT/A1_REPRESENTATIVE_FRAME_AUDIT.json" \
    "$OUT/PREFLIGHT_GATE.json" \
    "$OUT/representative_354ps.gro" \
    "$OUT/representative_354ps.pdb" \
    "$GITHUB_COMMIT" <<'PY'
import hashlib
import json
import pathlib
import sys
output, audit_path, preflight_path, gro_path, pdb_path, commit = sys.argv[1:]
audit = json.load(open(audit_path, encoding="utf-8"))
preflight = json.load(open(preflight_path, encoding="utf-8"))
if audit.get("scientific_status") != "PASS_A1_REPRESENTATIVE_NAC_FRAME":
    raise SystemExit("cannot write PASS without representative-frame PASS")
if preflight.get("status") != "PASS_TECHNICAL_UNRESTRAINED_GROMPP":
    raise SystemExit("cannot write PASS without unrestrained grompp PASS")
def sha256(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()
payload = {
    "schema_version": 1,
    "status": "PASS_A1_REPRESENTATIVE_FRAME_EXTRACTION",
    "scientific_status": "PASS_A1_REPRESENTATIVE_NAC_FRAME",
    "scientific_scope": (
        "classical_fixed_topology_preorganization_not_proton_transfer"
    ),
    "source_candidate": "nac_evt25_time1462ps",
    "velocity_seed": 26723,
    "selected_time_ps": 354.0,
    "grompp_maxwarn": 0,
    "preflight_status": preflight["status"],
    "promoted": True,
    "git_commit": commit,
    "outputs": {
        "gro": {"path": gro_path, "sha256": sha256(gro_path)},
        "pdb": {"path": pdb_path, "sha256": sha256(pdb_path)},
        "audit": {"path": audit_path, "sha256": sha256(audit_path)},
    },
}
pathlib.Path(output).write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
PY

STATE=PASS_TECHNICAL
DETAIL="source=evt25;seed=26723;time_ps=354;scientific_status=PASS_A1_REPRESENTATIVE_NAC_FRAME;output=$OUT"
trap - EXIT
append_history
