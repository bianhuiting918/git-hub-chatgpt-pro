#!/usr/bin/env bash
set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
TOOLS_ROOT="$TASK_ROOT/tools"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
SOURCE_JOB=61898124
SOURCE_PREFIX="$TOOLS_ROOT/amberclassic_job_61898124_0b35bfeb9602"
JOB_ID="${SLURM_JOB_ID:?SLURM_JOB_ID is required}"
PREFIX="$TOOLS_ROOT/amberclassic_activation_job_${JOB_ID}_0b35bfeb9602"
ACTIVE="$TOOLS_ROOT/ACTIVE_AMBERCLASSIC.json"
STATE=FAIL_TECHNICAL
DETAIL=initializing
EVENT=nylc_recover_amberclassic_install
COMMAND="recover_nylc_amberclassic_install.sh source_job=$SOURCE_JOB prefix=$PREFIX"
if [[ -e "$PREFIX" ]]; then printf 'Refusing to overwrite activation attempt: %s\n' "$PREFIX" >&2; exit 2; fi
mkdir -p "$PREFIX"
append_history() {
  local now
  now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
  (
    flock -x 9
    printf '%s\t%s\t%s\t%s\t%s\n' "$now" "$EVENT" "$COMMAND" "$STATE" "$DETAIL" >>"$TASK_ROOT/run_history.tsv"
    "$PY" - "$now" "$STATE" "$DETAIL" "$JOB_ID" "$PREFIX" "$SOURCE_PREFIX" >>"$TASK_ROOT/run_history.jsonl" <<'PY'
import json, sys
now, state, detail, job, prefix, source_prefix = sys.argv[1:]
print(json.dumps({"time":now,"event":"nylc_recover_amberclassic_install","command":"recover_nylc_amberclassic_install.sh","state":state,"detail":detail,"slurm_job_id":job,"prefix":prefix,"source_prefix":source_prefix}, sort_keys=True))
PY
  ) 9>"$TASK_ROOT/.run_history.lock"
}
finish() {
  code=$?
  trap - EXIT
  if [[ $code -ne 0 ]]; then DETAIL="exit_code=$code; prefix=$PREFIX; source_prefix=$SOURCE_PREFIX"; fi
  append_history
  exit "$code"
}
trap finish EXIT
"$PY" - "$SOURCE_PREFIX/dependency_relocation.json" <<'PY'
import json, pathlib, sys
record=json.loads(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"))
if record.get("status")!="PASS_DEPENDENCY_PATH_RELOCATION" or record.get("remaining_stale_occurrences")!=0:
    raise SystemExit("source dependency relocation did not pass")
PY
for rel in bin dat lib include AmberClassic.sh dependency_relocation.json compiler_versions.txt source.sha256; do
  test -e "$SOURCE_PREFIX/$rel"
  ln -s "$SOURCE_PREFIX/$rel" "$PREFIX/$rel"
done
source "$PREFIX/AmberClassic.sh"
for exe in antechamber parmchk2 resp respgen tleap sqm msander; do
  resolved="$(command -v "$exe")"
  case "$resolved" in "$PREFIX"/bin/*) ;; *) printf 'Resolved %s outside recovered prefix: %s\n' "$exe" "$resolved" >&2; exit 4 ;; esac
  printf '%s\n' "$resolved" >"$PREFIX/$exe.path"
done
antechamber -h >"$PREFIX/antechamber_help.stdout" 2>"$PREFIX/antechamber_help.stderr"
"$PY" - "$PREFIX" "$ACTIVE" "$SOURCE_PREFIX" "$SOURCE_JOB" "$JOB_ID" <<'PY'
import hashlib, json, pathlib, sys
prefix=pathlib.Path(sys.argv[1]); active=pathlib.Path(sys.argv[2]); source_prefix=pathlib.Path(sys.argv[3]); source_job,job=sys.argv[4:]
required=[prefix/"dependency_relocation.json",prefix/"dat"/"antechamber"/"CONNECT.TPL",prefix/"bin"/"antechamber",prefix/"bin"/"parmchk2",prefix/"bin"/"resp",prefix/"bin"/"respgen",prefix/"bin"/"tleap",prefix/"bin"/"sqm",prefix/"bin"/"msander"]
source_sha=(source_prefix/"source.sha256").read_text(encoding="utf-8").split()[0]
record={"schema_version":1,"status":"PASS_AMBERCLASSIC_INSTALL","source_repository":"Amber-MD/AmberClassic","source_commit":"0b35bfeb96026ffa4e5876391a0828f39b3cfc8d","source_archive_sha256":source_sha,"slurm_job_id":int(job),"recovery_of_slurm_job":int(source_job),"source_prefix":str(source_prefix),"prefix":str(prefix),"required_files":{str(path.relative_to(prefix)):hashlib.sha256(path.read_bytes()).hexdigest() for path in required}}
(prefix/"PASS.json").write_text(json.dumps(record,indent=2,sort_keys=True)+"\n",encoding="utf-8")
active.write_text(json.dumps(record,indent=2,sort_keys=True)+"\n",encoding="utf-8")
PY
STATE=PASS_TECHNICAL
DETAIL="status=PASS_AMBERCLASSIC_INSTALL; recovery_of_slurm_job=$SOURCE_JOB; prefix=$PREFIX"
