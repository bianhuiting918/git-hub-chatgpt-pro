#!/bin/bash
#SBATCH -J nylc_a1_acpdb
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem-per-cpu=1000M
#SBATCH -t 00:10:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.err

set -uo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
INPUT="$TASK_ROOT/a1_activated_nac_20260726/parameterization/job_61887449/model/NTA1_CAP.input.mol2"
OUT="$TASK_ROOT/a1_activated_nac_20260726/diagnostics/antechamber_pdb_job_${SLURM_JOB_ID:?}"
mkdir -p "$OUT"
"$PY" - "$INPUT" "$OUT/NTA1_CAP.input.pdb" <<'PY'
import pathlib
import sys
source = pathlib.Path(sys.argv[1]).read_text().splitlines()
out = pathlib.Path(sys.argv[2])
section = None
atoms = {}
bonds = []
for line in source:
    if line.startswith("@<TRIPOS>"):
        section = line[9:].lower()
        continue
    fields = line.split()
    if section == "atom" and len(fields) >= 6:
        idx = int(fields[0])
        name = fields[1]
        short = {
            "CAP_N": "NT", "CAP_H": "HT", "CAP_C": "CM",
            "CAP_H1": "HM1", "CAP_H2": "HM2", "CAP_H3": "HM3",
        }.get(name, name)
        element = "".join(ch for ch in fields[5] if ch.isalpha()).split(".")[0]
        element = element[:2].title()
        atoms[idx] = (short, float(fields[2]), float(fields[3]), float(fields[4]), element)
    elif section == "bond" and len(fields) >= 4:
        bonds.append((int(fields[1]), int(fields[2])))
lines = []
for idx in sorted(atoms):
    name, x, y, z, element = atoms[idx]
    lines.append(
        f"ATOM  {idx:5d} {name:>4s} NTA A   1    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}"
    )
neighbors = {idx: [] for idx in atoms}
for left, right in bonds:
    neighbors[left].append(right)
    neighbors[right].append(left)
for idx in sorted(neighbors):
    if neighbors[idx]:
        lines.append("CONECT" + f"{idx:5d}" + "".join(f"{j:5d}" for j in sorted(neighbors[idx])))
lines.extend(["TER", "END"])
out.write_text("\n".join(lines) + "\n")
PY
module purge
module load amber/2018-hpcx-gcc-7.3.1
antechamber -i "$OUT/NTA1_CAP.input.pdb" -fi pdb \
  -o "$OUT/read_charge.mol2" -fo mol2 \
  -c rc -nc 0 -m 1 -at gaff2 -rn NTA1 -s 2 \
  >"$OUT/read_charge.stdout" 2>"$OUT/read_charge.stderr"
parse_rc=$?
printf 'parse_rc\t%s\n' "$parse_rc" >"$OUT/result.tsv"
now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
state=PASS_DIAGNOSTIC_PDB_PARSE
if [[ $parse_rc -ne 0 ]]; then state=FAIL_DIAGNOSTIC_PDB_PARSE; fi
printf '%s\t%s\t%s\t%s\t%s\n' "$now" nylc_a1_antechamber_pdb "antechamber -fi pdb -c rc" "$state" "job=${SLURM_JOB_ID};parse_rc=$parse_rc;output=$OUT" >>"$TASK_ROOT/run_history.tsv"
"$PY" - "$now" "$state" "$parse_rc" "$OUT" "${SLURM_JOB_ID}" >>"$TASK_ROOT/run_history.jsonl" <<'PY'
import json
import sys
now, state, rc, output, job = sys.argv[1:]
print(json.dumps({
    "time": now,
    "event": "nylc_a1_antechamber_pdb",
    "command": "antechamber -fi pdb -c rc",
    "state": state,
    "parse_rc": int(rc),
    "output": output,
    "slurm_job_id": job,
}, sort_keys=True))
PY
exit "$parse_rc"
