#!/bin/bash
set -euo pipefail
echo "SUPERSEDED_WRONG_NUCLEOPHILE_IDENTITY: atom 8896 is Thr262 OG1, not catalytic Thr267 OG1" >&2
exit 42
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PYTHON=/work/home/acshdt1dks/anaconda3/bin/python3.9
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
job_id="${1:?usage: audit_nylc_c18_free_1ns.sh SLURM_JOB_ID}"
candidate=nylc_c18_11854ps
candidate_root="$TASK_ROOT/candidates/$candidate"
root="$candidate_root/runs_step1_gs/free_1ns_$job_id"
index="$candidate_root/build/source_cycle.ndx"
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
export GMX_MAXBACKUP=-1

test -s "$root/UNRESTRAINED_1NS_COMPLETE.json"
test -s "$root/run.xtc"
test -s "$root/run.tpr"
test -s "$root/run.edr"
test -s "$index"

# NylC NAC atoms: Thr267 OG1 = 8896, L2 reactive carbonyl C/O = 10287/10288.
"$GMX" distance -f "$root/run.xtc" -s "$root/run.tpr" -select 'atomnr 8896 plus atomnr 10287' -oall "$root/nac_distance.xvg" >"$root/nac_distance.stdout" 2>"$root/nac_distance.stderr"
printf '[ nac_angle ]\n8896 10287 10288\n' >"$root/nac_angle.ndx"
printf '0\n' | "$GMX" angle -f "$root/run.xtc" -n "$root/nac_angle.ndx" -ov "$root/nac_angle.xvg" -type angle >"$root/nac_angle.stdout" 2>"$root/nac_angle.stderr"

printf 'Potential\n0\n' | "$GMX" energy -f "$root/run.edr" -o "$root/potential_energy.xvg" >"$root/potential_energy.stdout" 2>"$root/potential_energy.stderr"
printf 'Temperature\nPressure\nVolume\n0\n' | "$GMX" energy -f "$root/run.edr" -o "$root/thermo.xvg" >"$root/thermo.stdout" 2>"$root/thermo.stderr"

# Gate is residues 261-266 only; source_cycle.ndx Gate excludes Thr267.
"$GMX" distance -f "$root/run.xtc" -s "$root/run.tpr" -n "$index" -select 'com of group "Core" plus com of group "Gate"' -oxyz "$root/gate_core_vector.xvg" >"$root/gate_core_vector.stdout" 2>"$root/gate_core_vector.stderr"

"$PYTHON" "$FLOW/scripts/analyze_nac_series.py" \
  --distance-xvg "$root/nac_distance.xvg" \
  --angle-xvg "$root/nac_angle.xvg" \
  --potential-xvg "$root/potential_energy.xvg" \
  --distance-max-nm 0.35 \
  --angle-min-deg 95 \
  --angle-max-deg 115 \
  --output "$root/nac_series_energy_audit.json"

"$PYTHON" - "$root" "$job_id" <<'PY'
import json,math,statistics,sys
from pathlib import Path
root=Path(sys.argv[1]); job=sys.argv[2]

def read_xvg(name):
    rows=[]
    for line in (root/name).read_text().splitlines():
        if not line.strip() or line[0] in "#@":
            continue
        rows.append([float(x) for x in line.split()])
    return rows

gate=read_xvg("gate_core_vector.xvg")
v=(0.4904295935403325,0.813080787402625,-0.3136533866174433)
baseline=-1.51109
grows=[]
for row in gate:
    t,x,y,z=row[:4]
    absolute=x*v[0]+y*v[1]+z*v[2]
    grows.append((t,absolute-baseline))
thermo=read_xvg("thermo.xvg")
labels=("temperature_K","pressure_bar","volume_nm3")
thermo_summary={}
for col,label in enumerate(labels,1):
    vals=[r[col] for r in thermo]
    thermo_summary[label]={
        "mean":statistics.fmean(vals),
        "stdev":statistics.pstdev(vals),
        "min":min(vals),
        "max":max(vals),
    }
nac=json.loads((root/"nac_series_energy_audit.json").read_text())
numerical=json.loads((root/"numerical_audit.json").read_text())
out={
 "schema_version":1,
 "candidate_id":"nylc_c18_11854ps",
 "slurm_job_id":job,
 "window_type":"fully_unrestrained_NPT",
 "nac":nac,
 "gate_definition":"NylC residues 261-266; Thr267 excluded",
 "gate_opening_nm":{
   "frame_count":len(grows),
   "mean":statistics.fmean(x[1] for x in grows),
   "stdev":statistics.pstdev(x[1] for x in grows),
   "min":min(x[1] for x in grows),
   "max":max(x[1] for x in grows),
 },
 "thermodynamics":thermo_summary,
 "numerical_audit":numerical,
 "scientific_status":"NOT_EVALUATED_PENDING_SCIENTIFIC_INTERPRETATION"
}
(root/"step1_free_1ns_audit.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
PY

selected_time="$("$PYTHON" - "$root/nac_series_energy_audit.json" <<'PY'
import json,sys
x=json.load(open(sys.argv[1])).get("lowest_potential_nac_frame")
print("" if x is None else x["time_ps"])
PY
)"
if [[ -n "$selected_time" ]]; then
  selected="$candidate_root/selected_step1_gs/c18_free_1ns_${job_id}_lowest_potential_nac.gro"
  printf '0\n' | "$GMX" trjconv -f "$root/run.xtc" -s "$root/run.tpr" -dump "$selected_time" -o "$selected" >"$root/extract_selected.stdout" 2>"$root/extract_selected.stderr"
  selected_sha="$(sha256sum "$selected" | awk '{print $1}')"
  "$PYTHON" - "$root/step1_free_1ns_audit.json" "$selected_time" "$selected" "$selected_sha" <<'PY'
import json,sys
path,time,gro,sha=sys.argv[1:]
x=json.load(open(path))
x["selected_lowest_potential_nac_frame"]={"time_ps":float(time),"remote_gro":gro,"gro_sha256":sha}
open(path,"w").write(json.dumps(x,indent=2,sort_keys=True)+"\n")
PY
fi

cat "$root/step1_free_1ns_audit.json"
