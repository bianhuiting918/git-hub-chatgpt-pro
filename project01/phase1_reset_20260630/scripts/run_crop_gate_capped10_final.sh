#!/usr/bin/env bash
set -euo pipefail

ROOT=/data/bht/project01_phase1_reset_gpu
PYTHON=/data/bht/design_tools/envs/placer_env/bin/python
RUN_DIR="$ROOT/placer_runs/capped10_n50_bonded"
OUT_DIR="$ROOT/placer_crop_gate_capped10_final"
expected=50

csv_count=$(find "$RUN_DIR" -maxdepth 1 -name '*.csv' | wc -l)
pdb_count=$(find "$RUN_DIR" -maxdepth 1 -name '*_model.pdb' | wc -l)

if [[ "$csv_count" -lt "$expected" || "$pdb_count" -lt "$expected" ]]; then
  echo "NOT_READY csv=$csv_count pdb=$pdb_count expected=$expected" >&2
  exit 2
fi

rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

"$PYTHON" "$ROOT/scripts/two_layer_placer_crop_gate.py" \
  --reference "$ROOT/manifests/denovo_SER_hydrolase_full_input.pdb" \
  --run-dir capped10_bonded "$RUN_DIR" \
  --out-dir "$OUT_DIR" \
  --expected-conformers-per-sequence 50 \
  --minimum-crop-pass-conformers-per-sequence 10 | tee "$OUT_DIR/run.log"
