#!/usr/bin/env bash
#SBATCH --job-name=autosite_smoke
#SBATCH --partition=xahcnormal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=01:00:00
#SBATCH --output=/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725/logs/%x_%j.out
#SBATCH --error=/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725/logs/%x_%j.err

set -euo pipefail
ROOT=/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725
ADFR="$ROOT/software/ADFRsuite-1.0/bin"
RUN="$ROOT/results/autosite_technical_smoke_20260726"
SUMMARY="$RUN/summary.tsv"

if [ -e "$RUN" ]; then
  printf '%s FAIL_EXISTING_OUTPUT path=%s\n' "$(date -Is)" "$RUN"
  exit 2
fi
mkdir -p "$RUN"
printf 'candidate_id\tinput_pdb\tinput_sha256\tprepare_rc\tautosite_rc\toutput_files\ttechnical_status\n' > "$SUMMARY"

run_one() {
  candidate_id=$1
  input_pdb=$2
  candidate_dir="$RUN/$candidate_id"
  mkdir -p "$candidate_dir"
  input_sha=$(sha256sum "$input_pdb" | awk '{print $1}')
  receptor="$candidate_dir/receptor.pdbqt"

  set +e
  "$ADFR/prepare_receptor" -r "$input_pdb" -o "$receptor" \
    -A checkhydrogens -U nphs_lps_waters_nonstdres \
    >"$candidate_dir/prepare.stdout.log" 2>"$candidate_dir/prepare.stderr.log"
  prepare_rc=$?
  set -e

  autosite_rc=NA
  status=NOT_EVALUATED_RECEPTOR_PREP
  if [ "$prepare_rc" -eq 0 ] && [ -s "$receptor" ]; then
    set +e
    "$ADFR/autosite" -r "$receptor" -o "$candidate_dir/autosite" -n 10 \
      >"$candidate_dir/autosite.stdout.log" 2>"$candidate_dir/autosite.stderr.log"
    autosite_rc=$?
    set -e
    if [ "$autosite_rc" -eq 0 ]; then
      status=PASS_AUTOSITE_TECHNICAL_SMOKE
    else
      status=NOT_EVALUATED_AUTOSITE
    fi
  fi
  output_files=$(find "$candidate_dir" -type f | wc -l)
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$candidate_id" "$input_pdb" "$input_sha" "$prepare_rc" "$autosite_rc" "$output_files" "$status" \
    >> "$SUMMARY"
}

run_one PET_6ILW /work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/inputs/receptors/SEQMD5_8ddcd9d7d7bbf5872959258e293cb55b.pdb
run_one NYLC_3AXG /work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/structures/saprot_nylonase_structure_recovery_20260725_v4/priority813_public/SEQMD5_b48a5e7e77b3b427ca89ac65da361050.pdb

sha256sum "$SUMMARY" > "$RUN/SHA256SUMS"
if grep -q 'NOT_EVALUATED' "$SUMMARY"; then
  printf '%s COMPLETE_WITH_TECHNICAL_FAILURES output=%s\n' "$(date -Is)" "$RUN"
  exit 1
fi
printf '%s PASS_TECHNICAL_ONLY output=%s\n' "$(date -Is)" "$RUN"
