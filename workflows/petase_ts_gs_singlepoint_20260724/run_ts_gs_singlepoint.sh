#!/usr/bin/env bash
set -euo pipefail

ROOT=/Dell/Dell14/bianht/petase_blind_qmmm/repo/project01/phase2_blind_petase_qmmm_20260630
BLIND="$ROOT/blind_work"
TASK="$BLIND/10_energy_diagnostics/attempt155_petase_step1_step2_ts_gs_dftb3_singlepoint"
SANDER=/Dell/Dell14/bianht/micromamba/envs/petase_stage1/bin/sander
MAX_PARALLEL="${MAX_PARALLEL:-4}"

mkdir -p "$TASK"/{input,structures,work,audit}
HISTORY="$TASK/run_history.tsv"
if [[ ! -s "$HISTORY" ]]; then
  printf 'timestamp\tevent\tcase_id\tstatus\tdetail\n' > "$HISTORY"
fi

STEP1_TOPO="$BLIND/05_atesa_acylation/attempt058_qm65_gap_boundary_from057_8x4_1500_corrected_atoms/complex.prmtop"
STEP2_TOPO="$BLIND/07_focused_committor/deacylation_seed01_seed04_longcommit_attempt007_2x16/complex.prmtop"

cat > "$TASK/input/step1_sp.in" <<'EOF'
PETase step1 fixed-coordinate QM/MM single point
&cntrl
 imin=1, maxcyc=0, ntmin=1,
 ntx=1, irest=0, ntb=1, cut=8.0,
 ifqnt=1, igb=0, nmropt=0, ntpr=1,
/
&qmmm
 qmmask='@1884,1885,1886,1887,1888,1889,1890,1891,1892,1893,1894,1895,1896,1897,1898,1899,1900,1901,1908,1909,1910,1911,1912,2533,2534,2535,2536,2537,2538,2984,2985,2986,2987,2988,2989,2990,2991,2992,2993,2994,3001,3002,3003,3004,3005,3837,3838,3839,3840,3841,3842,3843,3844,3845,3846,3847,3848,3849,3850,3851,3852,3853,3854,3855,3856',
 qmcharge=-1, spin=1, qm_theory='DFTB3',
 dftb_telec=100, dftb_maxiter=200, itrmax=3000,
 qm_ewald=0, qmshake=0, qmcut=8.0, writepdb=0,
/
EOF

cat > "$TASK/input/step2_sp.in" <<'EOF'
PETase step2 fixed-coordinate QM/MM single point
&cntrl
 imin=1, maxcyc=0, ntmin=1,
 ntx=1, irest=0, ntb=1, cut=8.0,
 ifqnt=1, igb=0, nmropt=0, ntpr=1,
/
&qmmm
 qmmask='@1908,1909,1910,1911,1912,2533,2534,2535,2536,2537,2538,2984,2985,2986,2987,2988,2989,2990,2991,2992,2993,2994,3837,3838,3839,3840,3841,3842,3843,3844,3845,3846,3847,3848,3849,3850,3851,3852,3853,3854,3855,3856,3863,3864,3865',
 qmcharge=-1, qm_theory='DFTB3', qmshake=0,
 writepdb=0, qm_ewald=1,
/
EOF

cat > "$TASK/input/cases.tsv" <<EOF
case_id	step	state	topology	mdin	coordinate
s1_ts_g05	step1	TS	$STEP1_TOPO	$TASK/input/step1_sp.in	$BLIND/05_atesa_acylation/attempt060_qm65_g05_ts_core_rc_lmax_prep/g05_bf03_cand02_rep46_f001_rep12_f001.rst7
s1_gs_michaelis_low	step1	GS	$STEP1_TOPO	$TASK/input/step1_sp.in	$BLIND/05_atesa_acylation/attempt117_hybrid_rc_anchor_refinement_plan/anchors/michaelis_lowrc_w02.rst7
s1_gs_michaelis_mid	step1	GS	$STEP1_TOPO	$TASK/input/step1_sp.in	$BLIND/05_atesa_acylation/attempt117_hybrid_rc_anchor_refinement_plan/anchors/michaelis_midrc_w10r.rst7
s1_gs_michaelis_high	step1	GS	$STEP1_TOPO	$TASK/input/step1_sp.in	$BLIND/05_atesa_acylation/attempt117_hybrid_rc_anchor_refinement_plan/anchors/michaelis_highrc_b175.rst7
s2_ts_seed1	step2	TS	$STEP2_TOPO	$TASK/input/step2_sp.in	$BLIND/07_focused_committor/deacylation_seed01_seed04_longcommit_attempt007_2x16/seeds/s01.rst7
s2_gs_reactant_q08	step2	GS	$STEP2_TOPO	$TASK/input/step2_sp.in	$BLIND/07_focused_committor/deacylation_seed01_seed04_longcommit_attempt007_2x16/work/s01_q08/md.rst7
s2_gs_reactant_q09	step2	GS	$STEP2_TOPO	$TASK/input/step2_sp.in	$BLIND/07_focused_committor/deacylation_seed01_seed04_longcommit_attempt007_2x16/work/s01_q09/md.rst7
s2_gs_reactant_q12	step2	GS	$STEP2_TOPO	$TASK/input/step2_sp.in	$BLIND/07_focused_committor/deacylation_seed01_seed04_longcommit_attempt007_2x16/work/s01_q12/md.rst7
s2_gs_reactant_q15	step2	GS	$STEP2_TOPO	$TASK/input/step2_sp.in	$BLIND/07_focused_committor/deacylation_seed01_seed04_longcommit_attempt007_2x16/work/s01_q15/md.rst7
EOF

sha256sum "$SANDER" "$STEP1_TOPO" "$STEP2_TOPO" "$TASK/input/step1_sp.in" "$TASK/input/step2_sp.in" > "$TASK/audit/input_sha256.tsv"

run_case() {
  local case_id="$1" step="$2" state="$3" topology="$4" mdin="$5" coordinate="$6"
  local dir="$TASK/work/$case_id"
  mkdir -p "$dir"
  printf '%s\tSTART\t%s\tRUNNING\t%s %s\n' "$(date -Is)" "$case_id" "$step" "$state" >> "$HISTORY"
  set +e
  "$SANDER" -O -i "$mdin" -o "$dir/sp.out" -p "$topology" -c "$coordinate" -r "$dir/sp.rst7" -inf "$dir/sp.mdinfo" > "$dir/launcher.stdout" 2> "$dir/launcher.stderr"
  local rc=$?
  set -e
  printf '%s\tFINISH\t%s\trc=%s\t%s %s\n' "$(date -Is)" "$case_id" "$rc" "$step" "$state" >> "$HISTORY"
  return "$rc"
}
export -f run_case
export TASK SANDER HISTORY

if [[ -n "${CASE_FILTER:-}" ]]; then
  awk -F '\t' -v id="$CASE_FILTER" 'NR > 1 && $1 == id' "$TASK/input/cases.tsv" |
    xargs -P 1 -n 6 bash -c 'run_case "$@"' _
  grep -E 'FINAL RESULTS|NSTEP|ENERGY|ERROR|NaN|SCC' "$TASK/work/$CASE_FILTER/sp.out" |
    tail -20 > "$TASK/audit/${CASE_FILTER}_smoke_excerpt.txt" || true
else
  tail -n +2 "$TASK/input/cases.tsv" |
    xargs -P "$MAX_PARALLEL" -n 6 bash -c 'run_case "$@"' _
  python3 "$TASK/analyze_ts_gs_singlepoint.py" --task "$TASK"
fi
