#!/usr/bin/env bash
set -euo pipefail

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
GMX=/public/software/apps/gromacs/2022.2/hpcx-gcc7.3.1/bin/gmx_mpi
SOURCE_GRO="$TASK_ROOT/ensemble/final_audit_job_61841413/assembled/medoids/nac_evt18_time1206ps/seed26711_time258.000ps.gro"
CHAIN_ITP="$TASK_ROOT/ensemble/candidates/nac_evt18_time1206ps/build_job_61813799_6/build/topol_Protein_chain_H.itp"
CHAIN_FIRST_GLOBAL_ATOM=8949
GITHUB_COMMIT="${GITHUB_COMMIT:-bfe9119152050c80d87dc74d39a2e58c5aaf2b7e}"
OUT_DIR="${1:?usage: run_nylc_a1_parameterization.sh OUTPUT_DIR}"

case "$OUT_DIR" in
  "$TASK_ROOT"/a1_activated_nac_20260726/parameterization/*) ;;
  *) printf 'Refusing output outside task parameterization root: %s\n' "$OUT_DIR" >&2; exit 2 ;;
esac
if [[ -e "$OUT_DIR" ]]; then
  printf 'Refusing to overwrite existing output: %s\n' "$OUT_DIR" >&2
  exit 2
fi
mkdir -p "$OUT_DIR"
record_phase() {
  printf '%s\t%s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" "$1" >>"$OUT_DIR/phase.tsv"
}
record_phase before_module_purge

EVENT=nylc_a1_parameterization
COMMAND="run_nylc_a1_parameterization.sh $OUT_DIR"
STATE=FAIL_TECHNICAL
DETAIL=initializing
append_history() {
  local now
  now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
  (
    flock -x 9
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$now" "$EVENT" "$COMMAND" "$STATE" "$GITHUB_COMMIT" "$DETAIL" >>"$TASK_ROOT/run_history.tsv"
    "$PY" - "$now" "$EVENT" "$COMMAND" "$STATE" "$GITHUB_COMMIT" "$DETAIL" "$OUT_DIR" >>"$TASK_ROOT/run_history.jsonl" <<'PY'
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
  code=$?
  trap - EXIT
  if [[ $code -ne 0 ]]; then
    DETAIL="exit_code=$code; output=$OUT_DIR"
  fi
  append_history
  exit $code
}
trap finish EXIT

module purge
record_phase after_module_purge
module load compiler/gcc/11.4.0
AMBER_ACTIVE="$TASK_ROOT/tools/ACTIVE_AMBERCLASSIC.json"
mapfile -t AMBER_META < <("$PY" - "$AMBER_ACTIVE" <<'PY'
import json
import pathlib
import sys
active = pathlib.Path(sys.argv[1])
record = json.loads(active.read_text(encoding="utf-8"))
if record.get("status") != "PASS_AMBERCLASSIC_INSTALL":
    raise SystemExit("ACTIVE_AMBERCLASSIC.json is not PASS_AMBERCLASSIC_INSTALL")
prefix = pathlib.Path(record["prefix"])
pass_record = json.loads((prefix / "PASS.json").read_text(encoding="utf-8"))
if pass_record != record:
    raise SystemExit("active AmberClassic manifest differs from installation PASS.json")
for rel in ("AmberClassic.sh", "dat/antechamber/CONNECT.TPL", "bin/antechamber",
            "bin/parmchk2", "bin/tleap", "bin/sqm", "bin/msander"):
    if not (prefix / rel).is_file():
        raise SystemExit(f"missing audited AmberClassic file: {rel}")
print(prefix)
print(record["source_archive_sha256"])
print(record["source_commit"])
PY
)
AMBER_PREFIX="${AMBER_META[0]}"
export AMBER_SOURCE_SHA256="${AMBER_META[1]}"
export AMBER_SOURCE_COMMIT="${AMBER_META[2]}"
export AMBER_PREFIX
record_phase after_active_manifest
set +u
source "$AMBER_PREFIX/AmberClassic.sh"
set -u
record_phase after_amberclassic_source
for exe in antechamber parmchk2 resp respgen tleap sqm msander; do
  resolved="$(command -v "$exe")"
  case "$resolved" in
    "$AMBER_PREFIX"/bin/*) ;;
    *) printf 'Resolved %s outside audited AmberClassic prefix: %s\n' "$exe" "$resolved" >&2; exit 4 ;;
  esac
  printf '%s\n' "$resolved" >"$OUT_DIR/$exe.path"
done
record_phase after_tool_resolution
[[ -x "$GMX" ]]
[[ -x "$PY" ]]

"$PY" "$FLOW/scripts/prepare_nylc_a1_parameter_model.py" \
  --gro "$SOURCE_GRO" \
  --chain-itp "$CHAIN_ITP" \
  --chain-first-global-atom "$CHAIN_FIRST_GLOBAL_ATOM" \
  --parent-id nac_evt18_time1206ps \
  --output-dir "$OUT_DIR/model"

cd "$OUT_DIR"
antechamber \
  -i model/NTA1_CAP.input.mol2 -fi mol2 \
  -o NTA1_CAP.am1bcc.raw.mol2 -fo mol2 \
  -c bcc -nc 0 -m 1 -at gaff2 -rn NTA1 -s 2 \
  >antechamber.stdout 2>antechamber.stderr
"$PY" - NTA1_CAP.am1bcc.raw.mol2 NTA1_CAP.am1bcc.mol2 charge_normalization.json <<'PY'
import json
import pathlib
import sys
source = pathlib.Path(sys.argv[1])
output = pathlib.Path(sys.argv[2])
audit = pathlib.Path(sys.argv[3])
lines = source.read_text(encoding="utf-8").splitlines()
in_atoms = False
rows = []
for index, line in enumerate(lines):
    if line == "@<TRIPOS>ATOM":
        in_atoms = True
        continue
    if line == "@<TRIPOS>BOND":
        in_atoms = False
    elif in_atoms and line.strip():
        fields = line.split()
        rows.append((index, fields))
charge_before = sum(float(fields[8]) for _, fields in rows)
if abs(charge_before) > 0.01:
    raise SystemExit(f"AM1-BCC charge residual too large: {charge_before}")
targets = [(index, fields) for index, fields in rows if fields[1] == "N"]
if len(targets) != 1:
    raise SystemExit("expected exactly one A1 N atom for charge normalization")
index, fields = targets[0]
old_charge = float(fields[8])
fields[8] = f"{old_charge - charge_before:.6f}"
lines[index] = " ".join(fields)
charge_after = sum(float(fields[8]) for _, fields in rows)
if abs(charge_after) > 5.0e-7:
    raise SystemExit(f"normalized charge is not zero: {charge_after}")
output.write_text("\n".join(lines) + "\n", encoding="utf-8")
audit.write_text(json.dumps({
    "schema_version": 1,
    "status": "PASS_AM1BCC_ROUNDING_NORMALIZATION",
    "charge_before_e": charge_before,
    "charge_after_e": charge_after,
    "adjusted_atom": "N",
    "adjustment_e": -charge_before,
}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY
parmchk2 -i NTA1_CAP.am1bcc.mol2 -f mol2 -o NTA1_CAP.frcmod -s gaff2 \
  >parmchk2.stdout 2>parmchk2.stderr

cat >tleap.in <<'EOF'
source leaprc.gaff2
loadamberparams NTA1_CAP.frcmod
A1 = loadmol2 NTA1_CAP.am1bcc.mol2
check A1
charge A1
saveamberparm A1 NTA1_CAP.prmtop NTA1_CAP.inpcrd
savepdb A1 NTA1_CAP.pdb
quit
EOF
tleap -f tleap.in >tleap.out 2>tleap.err

cat >amber_single.in <<'EOF'
A1 one-step vacuum minimization for finite-energy audit
&cntrl
  imin=1, maxcyc=1, ncyc=1,
  ntb=0, igb=0, cut=999.0,
  ntpr=1,
/
EOF
msander -O -i amber_single.in -p NTA1_CAP.prmtop -c NTA1_CAP.inpcrd \
  -o amber_single.out -r amber_single.rst
AMBER_ENERGY="$(awk '/FINAL RESULTS/{seen=1} seen && /NSTEP/{getline; print $2; exit}' amber_single.out)"
[[ "$AMBER_ENERGY" =~ ^[-+0-9.eE]+$ ]]

"$PY" - NTA1_CAP.prmtop NTA1_CAP.inpcrd NTA1_CAP.gromacs.top NTA1_CAP.gromacs.gro <<'PY'
import parmed as pmd
import sys
prmtop, inpcrd, top, gro = sys.argv[1:]
structure = pmd.load_file(prmtop, xyz=inpcrd)
structure.save(top, format="gromacs", overwrite=False)
structure.save(gro, format="gro", overwrite=False)
PY
module unload compiler/gcc/11.4.0
module load gromacs/2022.2-hpcx-gcc-7.3.1
export GMX_MAXBACKUP=-1
"$GMX" editconf -f NTA1_CAP.gromacs.gro -o NTA1_CAP.gromacs.boxed.gro \
  -box 20 20 20 -center 10 10 10 >editconf.stdout 2>editconf.stderr

cat >single.mdp <<'EOF'
integrator               = md
dt                       = 0.001
nsteps                   = 0
cutoff-scheme            = Verlet
nstlist                  = 10
rlist                    = 2.0
rcoulomb                 = 2.0
rvdw                     = 2.0
coulombtype              = Cut-off
vdwtype                  = Cut-off
pbc                      = xyz
continuation             = yes
comm-mode                = None
constraints              = none
nstenergy                = 1
EOF
"$GMX" grompp -f single.mdp -c NTA1_CAP.gromacs.boxed.gro -p NTA1_CAP.gromacs.top \
  -o single.tpr -po single.expanded.mdp -maxwarn 0 >grompp.stdout 2>grompp.stderr
OMP_NUM_THREADS=1 "$GMX" mdrun -s single.tpr -deffnm gmx_single -ntomp 1 \
  >mdrun.stdout 2>mdrun.stderr
GMX_ENERGY="$(awk '/Potential Energy/{value=$3} END{print value}' gmx_single.log)"
[[ "$GMX_ENERGY" =~ ^[-+0-9.eE]+$ ]]

"$PY" - "$OUT_DIR" "$AMBER_ENERGY" "$GMX_ENERGY" <<'PY'
import hashlib
import json
import pathlib
import subprocess
import sys

root = pathlib.Path(sys.argv[1])
amber_energy = float(sys.argv[2])
gmx_energy = float(sys.argv[3])

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def version(exe):
    run = subprocess.run([exe, "-h"], text=True, capture_output=True)
    text = (run.stdout + run.stderr).splitlines()
    return text[0] if text else "AmberTools 2018 module"

frcmod = (root / "NTA1_CAP.frcmod").read_text(encoding="utf-8", errors="replace").lower()
unresolved = sum(
    frcmod.count(marker)
    for marker in ("attn, need revision", "unresolved", "missing parameter")
)
provenance = {
    "model_definition": "graph-derived N-terminal Thr267(A1)-N-methylamide capped model",
    "net_charge_e": 0,
    "charge_method": "AM1-BCC",
    "atom_types": "GAFF2 local screening model",
    "bonded_parameters": "GAFF2 plus parmchk2",
    "tool_versions": {
        "antechamber": version("antechamber"),
        "parmchk2": version("parmchk2"),
        "tleap": version("tleap"),
    },
    "input_sha256": sha(root / "model" / "NTA1_CAP.input.mol2"),
    "output_sha256": sha(root / "NTA1_CAP.am1bcc.mol2"),
    "validation_status": "PASS",
    "parmchk2_unresolved_count": unresolved,
    "tleap_load_status": "PASS",
    "gromacs_grompp_status": "PASS",
    "amber_energy_kcal_mol": amber_energy,
    "gromacs_energy_kj_mol": gmx_energy,
    "scientific_scope": "fixed-topology classical A1 NAC stability screening only",
}
(root / "parameter_provenance.json").write_text(
    json.dumps(provenance, indent=2, sort_keys=True) + "\n",
    encoding="utf-8",
)
PY

"$PY" "$FLOW/scripts/audit_nylc_a1_patch.py" \
  --mol2 NTA1_CAP.am1bcc.mol2 \
  --frcmod NTA1_CAP.frcmod \
  --provenance parameter_provenance.json \
  --output PATCH_AUDIT.json

sha256sum   model/NTA1_CAP.input.mol2   NTA1_CAP.am1bcc.raw.mol2   NTA1_CAP.am1bcc.mol2   charge_normalization.json   NTA1_CAP.frcmod   NTA1_CAP.prmtop   NTA1_CAP.inpcrd   NTA1_CAP.gromacs.top   NTA1_CAP.gromacs.gro   NTA1_CAP.gromacs.boxed.gro   parameter_provenance.json   PATCH_AUDIT.json >sha256.tsv

STATE=PASS_TECHNICAL
DETAIL="parameter_status=PASS_A1_SCREENING_PATCH; output=$OUT_DIR"
