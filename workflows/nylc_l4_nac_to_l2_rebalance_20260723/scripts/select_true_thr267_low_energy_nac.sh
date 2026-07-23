#!/bin/bash
set -euo pipefail

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PYTHON=/work/home/acshdt1dks/anaconda3/bin/python3.9
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
run_root="${1:?usage: select_true_thr267_low_energy_nac.sh RUN_ROOT OUTPUT_ROOT}"
output_root="${2:?usage: select_true_thr267_low_energy_nac.sh RUN_ROOT OUTPUT_ROOT}"

source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
export GMX_MAXBACKUP=-1
mkdir -p "$output_root"

test -s "$run_root/release2.xtc"
test -s "$run_root/release2.tpr"
test -s "$run_root/release2.edr"
test -s "$run_root/release2_pullx.xvg"
test ! -e "$output_root/selected_lowest_potential_nac.gro"
test ! -e "$output_root/selected_lowest_potential_nac.json"

# pullx columns are time, true Thr267 OG1--reactive carbonyl distance, and
# Thr267 OG1--C--O angle. Preserve comments while splitting the two observables.
awk '/^[#@]/{print;next}{print $1,$2}' "$run_root/release2_pullx.xvg" >"$output_root/nac_distance.xvg"
awk '/^[#@]/{print;next}{print $1,$3}' "$run_root/release2_pullx.xvg" >"$output_root/nac_angle.xvg"
printf 'Potential\n0\n' | "$GMX" energy \
  -f "$run_root/release2.edr" \
  -o "$output_root/potential_energy.xvg" \
  >"$output_root/energy.stdout" 2>"$output_root/energy.stderr"

"$PYTHON" "$FLOW/scripts/analyze_nac_series.py" \
  --distance-xvg "$output_root/nac_distance.xvg" \
  --angle-xvg "$output_root/nac_angle.xvg" \
  --potential-xvg "$output_root/potential_energy.xvg" \
  --distance-max-nm 0.35 \
  --angle-min-deg 95 \
  --angle-max-deg 115 \
  --output "$output_root/nac_energy_series.json"

selected_time="$("$PYTHON" - "$output_root/nac_energy_series.json" <<'PY'
import json,sys
record=json.load(open(sys.argv[1]))["lowest_potential_nac_frame"]
if record is None:
    raise SystemExit("no true-Thr267 NAC frame exists in the restrained source trajectory")
print(record["time_ps"])
PY
)"

printf '0\n' | "$GMX" trjconv \
  -f "$run_root/release2.xtc" \
  -s "$run_root/release2.tpr" \
  -dump "$selected_time" \
  -o "$output_root/selected_lowest_potential_nac.gro" \
  >"$output_root/trjconv.stdout" 2>"$output_root/trjconv.stderr"

selected_sha="$(sha256sum "$output_root/selected_lowest_potential_nac.gro" | awk '{print $1}')"
"$PYTHON" - "$output_root/nac_energy_series.json" "$output_root/selected_lowest_potential_nac.json" "$selected_time" "$selected_sha" "$run_root" <<'PY'
import json,sys
series_path,out_path,time,sha,source_root=sys.argv[1:]
series=json.load(open(series_path))
chosen=series["lowest_potential_nac_frame"]
out={
    "schema_version":1,
    "candidate_id":"nylc_C18_trueT267_recapture",
    "selection_rule":"minimum instantaneous potential energy among frames satisfying true Thr267 NAC distance <= 0.35 nm and angle 95-115 deg",
    "source_run_root":source_root,
    "source_window":"restrained_release2",
    "scientific_status":"NOT_A_SCIENTIFIC_GS_RESTRAINED_SOURCE",
    "selected_time_ps":float(time),
    "selected_potential_energy_kj_mol":chosen["potential_energy_kj_mol"],
    "selected_distance_nm":chosen["distance_nm"],
    "selected_angle_deg":chosen["angle_deg"],
    "selected_gro":"selected_lowest_potential_nac.gro",
    "selected_gro_sha256":sha,
    "trajectory_nac":{
        "frame_count":series["frame_count"],
        "nac_frame_count":series["nac_frame_count"],
        "occupancy":series["occupancy"],
        "longest_continuous_nac_ps":series["longest_continuous_nac_ps"],
    },
    "next_gate":"fully_unrestrained multi-seed NPT pilot; this restrained selection alone cannot establish GS stability",
}
open(out_path,"w").write(json.dumps(out,indent=2,sort_keys=True)+"\n")
PY

cat "$output_root/selected_lowest_potential_nac.json"
