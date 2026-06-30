#!/usr/bin/env bash
set -u

input_pdb="${1:-}"
out_dir="${2:-work/blind_work/01_system_setup/protonation_gate}"

if [[ -z "$input_pdb" ]]; then
  echo "Usage: bash run_stage1_protonation_gate.sh INPUT_PDB [OUT_DIR]" >&2
  exit 2
fi

mkdir -p "$out_dir"
report="$out_dir/protonation_gate_report.tsv"
md="$out_dir/protonation_gate_report.md"
: > "$report"
printf "category\ttool\tstatus\tpath_or_version\tcommand\toutput_path\toutput_sha256\tnote\n" >> "$report"

hash_file() {
  local file="$1"
  if [[ ! -s "$file" ]]; then
    printf ""
  elif command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$file" | awk '{print $1}'
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$file" | awk '{print $1}'
  else
    printf "sha256_tool_missing"
  fi
}

record_missing() {
  local category="$1"
  local tool="$2"
  printf "%s\t%s\tmissing\t\tcommand_not_found:%s\t\t\t\n" "$category" "$tool" "$tool" >> "$report"
}

record_file() {
  local category="$1"
  local tool="$2"
  local status="$3"
  local version="$4"
  local command="$5"
  local output="$6"
  local note="$7"
  local digest
  digest="$(hash_file "$output")"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$category" "$tool" "$status" "$version" "$command" "$output" "$digest" "$note" >> "$report"
}

if [[ ! -s "$input_pdb" ]]; then
  printf "input\tpdb\tfailed\t\tinput_missing_or_empty:%s\t\t\t\n" "$input_pdb" >> "$report"
  exit 2
fi

input_sha="$(hash_file "$input_pdb")"
printf "input\tpdb\tavailable\t%s\tprovided_input\t%s\t%s\tblind_stage1_clean_structure\n" "$input_pdb" "$input_pdb" "$input_sha" >> "$report"

propka_bin=""
if command -v propka >/dev/null 2>&1; then
  propka_bin="propka"
elif command -v propka3 >/dev/null 2>&1; then
  propka_bin="propka3"
fi

if [[ -n "$propka_bin" ]]; then
  propka_version="$($propka_bin --version 2>&1 | head -n 1 | tr '\t' ' ')"
  propka_log="$out_dir/propka_stdout_stderr.log"
  (cd "$out_dir" && "$propka_bin" "$(cd "$(dirname "$input_pdb")" && pwd)/$(basename "$input_pdb")") > "$propka_log" 2>&1
  record_file "pka" "$propka_bin" "ran_check_output" "$propka_version" "$propka_bin $input_pdb" "$propka_log" "inspect pka file for catalytic Asp/His states"
else
  record_missing "pka" "propka"
fi

if command -v pdb2pqr >/dev/null 2>&1; then
  pdb2pqr_version="$(pdb2pqr --version 2>&1 | head -n 1 | tr '\t' ' ')"
  pqr_out="$out_dir/protonated_pdb2pqr_ph7.pqr"
  pqr_log="$out_dir/pdb2pqr_stdout_stderr.log"
  pdb2pqr --ff=AMBER --with-ph=7.0 "$input_pdb" "$pqr_out" > "$pqr_log" 2>&1
  if [[ -s "$pqr_out" ]]; then
    record_file "protonation" "pdb2pqr" "produced_output" "$pdb2pqr_version" "pdb2pqr --ff=AMBER --with-ph=7.0 $input_pdb $pqr_out" "$pqr_out" "review residue states before topology generation"
  else
    record_file "protonation" "pdb2pqr" "ran_no_output" "$pdb2pqr_version" "pdb2pqr --ff=AMBER --with-ph=7.0 $input_pdb $pqr_out" "$pqr_log" "see stderr log"
  fi
else
  record_missing "protonation" "pdb2pqr"
fi

if command -v reduce >/dev/null 2>&1; then
  reduce_version="$(reduce -version 2>&1 | head -n 1 | tr '\t' ' ')"
  reduce_out="$out_dir/reduce_build_h.pdb"
  reduce -BUILD "$input_pdb" > "$reduce_out" 2> "$out_dir/reduce_stderr.log"
  if [[ -s "$reduce_out" ]]; then
    record_file "hydrogen" "reduce" "produced_output" "$reduce_version" "reduce -BUILD $input_pdb > $reduce_out" "$reduce_out" "compare His tautomer choices with pKa output"
  else
    record_file "hydrogen" "reduce" "ran_no_output" "$reduce_version" "reduce -BUILD $input_pdb > $reduce_out" "$out_dir/reduce_stderr.log" "see stderr log"
  fi
else
  record_missing "hydrogen" "reduce"
fi

{
  echo "# Stage 1 Protonation Gate Report"
  echo
  echo "Generated: $(date -Is)"
  echo
  echo "Input PDB: $input_pdb"
  echo "Input SHA256: $input_sha"
  echo
  echo "Boundary: this gate uses the cleaned PETase structure and pH/protonation tools only. It does not use paper TS structures, paper reaction coordinates, selected CVs, trajectories, barriers, or mechanism conclusions."
  echo
  echo "## Required Manual Review"
  echo
  echo "- Catalytic Asp state."
  echo "- Catalytic His tautomer/protonation."
  echo "- Nearby Cys/disulfide consistency."
  echo "- Remote His tautomer naming for topology."
  echo "- Any branch that differs from the primary hypothesis manifest."
  echo
  echo "## Probe Table"
  echo
  echo '```tsv'
  cat "$report"
  echo '```'
} > "$md"

echo "Wrote $report"
echo "Wrote $md"
