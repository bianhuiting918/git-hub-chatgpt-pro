#!/usr/bin/env bash
set -u

out_dir="${1:-blind_work/01_system_setup/environment_probe}"
mkdir -p "$out_dir"

report="$out_dir/stage1_compute_environment_probe.tsv"
: > "$report"
printf "category\tname\tstatus\tpath_or_version\tnote\n" >> "$report"

probe_cmd() {
  local category="$1"
  local name="$2"
  local cmd="$3"
  local version_cmd="${4:-}"
  local path
  path="$(command -v "$cmd" 2>/dev/null || true)"
  if [[ -n "$path" ]]; then
    local version="available"
    if [[ -n "$version_cmd" ]]; then
      version="$(eval "$version_cmd" 2>&1 | head -n 3 | tr '\t' ' ' | tr '\n' ';')"
    fi
    printf "%s\t%s\tavailable\t%s\t%s\n" "$category" "$name" "$path" "$version" >> "$report"
  else
    printf "%s\t%s\tmissing\t\tcommand_not_found:%s\n" "$category" "$name" "$cmd" >> "$report"
  fi
}

probe_python_module() {
  local category="$1"
  local module="$2"
  python3 - "$module" >> "$report" 2>/dev/null <<'PY'
import importlib.util
import sys
module = sys.argv[1]
status = "available" if importlib.util.find_spec(module) else "missing"
print(f"python_module\t{module}\t{status}\t{sys.executable}\tpython={sys.version.split()[0]}")
PY
  if [[ "${PIPESTATUS[0]}" -ne 0 ]]; then
    printf "%s\t%s\tmissing\t\tpython3_import_check_failed\n" "$category" "$module" >> "$report"
  fi
}

probe_cmd "shell" "bash" "bash" "bash --version"
probe_cmd "python" "python3" "python3" "python3 --version"
probe_cmd "chemistry" "obabel" "obabel" "obabel -V"
probe_cmd "chemistry" "babel" "babel" "babel -V"
probe_cmd "protonation" "propka" "propka" "propka --version"
probe_cmd "protonation" "propka3" "propka3" "propka3 --version"
probe_cmd "protonation" "pdb2pqr" "pdb2pqr" "pdb2pqr --version"
probe_cmd "protonation" "reduce" "reduce" "reduce -version"
probe_cmd "ambertools" "tleap" "tleap" "tleap -h"
probe_cmd "ambertools" "antechamber" "antechamber" "antechamber -h"
probe_cmd "ambertools" "parmchk2" "parmchk2" "parmchk2 -h"
probe_cmd "docking" "vina" "vina" "vina --version"
probe_cmd "docking" "smina" "smina" "smina --version"
probe_cmd "qmmm" "gmx" "gmx" "gmx --version"
probe_cmd "qmmm" "cp2k" "cp2k" "cp2k --version"

if command -v python3 >/dev/null 2>&1; then
  probe_python_module "python_module" "rdkit"
  probe_python_module "python_module" "openbabel"
  probe_python_module "python_module" "Bio"
fi

{
  echo "# Stage 1 Compute Environment Probe"
  echo
  echo "Generated: $(date -Is)"
  echo
  echo "Host: $(hostname 2>/dev/null || echo unknown)"
  echo "User: $(id -un 2>/dev/null || echo unknown)"
  echo "Working directory: $(pwd)"
  echo
  echo "## Probe Table"
  echo
  echo '```tsv'
  cat "$report"
  echo '```'
} > "$out_dir/stage1_compute_environment_probe.md"

echo "Wrote $report"
echo "Wrote $out_dir/stage1_compute_environment_probe.md"
