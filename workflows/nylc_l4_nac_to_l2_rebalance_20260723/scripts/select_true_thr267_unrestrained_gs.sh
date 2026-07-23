#!/bin/bash
set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PY=/work/home/acshdt1dks/anaconda3/bin/python3.9
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
pilot="$TASK/candidates/nylc_C18_trueT267_recapture/unrestrained_pilots/job_61705307/seed_26703"
followup="$TASK/candidates/nylc_C18_trueT267_recapture/unrestrained_1ns/job_61705692/unrestrained_1ns_audit.json"
output_root="$TASK/candidates/nylc_C18_trueT267_freeGS"
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
export GMX_MAXBACKUP=-1
mkdir -p "$output_root"

test ! -e "$output_root/source.gro"
test ! -e "$output_root/source_manifest.json"
test "$("$PY" -c "import json;print(json.load(open('$pilot/unrestrained_pilot_audit.json'))['scientific_status'])")" = PASS_UNRESTRAINED_PILOT_NAC_PRESENT
test -s "$followup"

selected_time="$("$PY" - "$pilot/nac_energy_series.json" <<'PY'
import json,sys
x=json.load(open(sys.argv[1]))["lowest_potential_nac_frame"]
if x is None:
    raise SystemExit("no unrestrained true-Thr267 NAC frame")
print(x["time_ps"])
PY
)"
printf '0\n' | "$GMX" trjconv -f "$pilot/free_npt100.xtc" -s "$pilot/free_npt100.tpr" \
  -dump "$selected_time" -o "$output_root/source.gro" \
  >"$output_root/trjconv.stdout" 2>"$output_root/trjconv.stderr"
sha="$(sha256sum "$output_root/source.gro" | awk '{print $1}')"

"$PY" - "$FLOW/scripts" "$pilot/nac_energy_series.json" "$pilot/unrestrained_pilot_audit.json" "$followup" "$output_root/source.gro" "$output_root/source_manifest.json" "$sha" <<'PY'
import json,pathlib,sys
script_dir,series_path,pilot_path,followup_path,gro_path,out_path,sha=sys.argv[1:]
sys.path.insert(0,script_dir)
import audit_true_thr267_recapture_pilot as audit
import extract_true_thr267_recapture_sources as geom
series=json.load(open(series_path))
pilot=json.load(open(pilot_path))
followup=json.load(open(followup_path))
chosen=series["lowest_potential_nac_frame"]
gro=pathlib.Path(gro_path)
distance,angle=audit.geometry(gro,geom.CANDIDATES["nylc_C18_trueT267_recapture"])
atoms,box,_=geom.read_gro(gro)
contact=audit.minimum_ligand_protein_distance(atoms,box)
if not (distance <= 0.35 and 95.0 <= angle <= 115.0):
    raise SystemExit("extracted frame does not satisfy true-Thr267 NAC")
if abs(distance-chosen["distance_nm"]) > 0.002 or abs(angle-chosen["angle_deg"]) > 0.2:
    raise SystemExit("extracted GRO geometry does not reproduce audited trajectory frame")
out={
  "schema_version":1,
  "candidate_id":"nylc_C18_trueT267_freeGS",
  "status":"VERIFIED_RARE_UNRESTRAINED_NAC_GS_CANDIDATE",
  "scientific_caveat":"NOT_ENSEMBLE_STABLE",
  "selection_rule":"lowest instantaneous potential energy among true-Thr267 NAC frames in the 100 ps fully unrestrained NPT seed 26703",
  "source":{
    "array_job_id":61705307,
    "velocity_seed":26703,
    "time_ps":chosen["time_ps"],
    "window_type":"fully_unrestrained_NPT_100ps",
    "gro":str(gro),
    "gro_sha256":sha,
  },
  "geometry":{
    "true_thr267_og1_to_carbonyl_c_nm":round(distance,6),
    "true_thr267_attack_angle_deg":round(angle,3),
    "joint_nac":True,
  },
  "potential_energy_kj_mol":chosen["potential_energy_kj_mol"],
  "minimum_ligand_protein_contact":contact,
  "gate_definition":"NylC residues 261-266; Thr267 excluded",
  "pilot_nac":pilot["nac"],
  "one_ns_followup":{
    "job_id":61705692,
    "frame_count":followup["nac"]["frame_count"],
    "nac_frame_count":followup["nac"]["nac_frame_count"],
    "occupancy":followup["nac"]["nac_occupancy"],
    "scientific_status":followup["scientific_status"],
  },
  "next_gate":"L4-to-L2 graph truncation and topology/charge/valence/contact/grompp audit; later fully unrestrained L2 validation remains required",
}
open(out_path,"w").write(json.dumps(out,indent=2,sort_keys=True)+"\n")
PY
cat "$output_root/source_manifest.json"
