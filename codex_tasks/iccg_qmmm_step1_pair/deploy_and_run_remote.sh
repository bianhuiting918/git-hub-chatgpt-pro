#!/bin/bash
# Deploy scripts to the remote project directory and submit only after PASS preflight.
set -euo pipefail
REMOTE_PROJECT=/work/home/acshdt1dks/iccg_qmmm_step1_pair_20260723
mkdir -p "$REMOTE_PROJECT"
cp build_iccg_step1_pair.py audit_iccg_step1_pair.py collect_iccg_step1_pair.py run_iccg_step1_pair.sbatch "$REMOTE_PROJECT"/
cd "$REMOTE_PROJECT"
python build_iccg_step1_pair.py \
  --iccg /work/home/acshdt1dks/iccg_sequence_design_20260723/inputs/ICCG_active_chainA.pdb \
  --lg1 /work/home/acshdt1dks/petase_qmmm_pilot8_20260721/inputs/reference/IsPETase_WT_LG1.pdb \
  --lg2 /work/home/acshdt1dks/petase_qmmm_pilot8_20260721/inputs/reference/IsPETase_WT_LG2.pdb \
  --out "$REMOTE_PROJECT"
status=$(python - <<'PY'
import json
from pathlib import Path
r=json.loads(Path('preflight_report.json').read_text())
print(r.get('status','NOT_SUBMITTED_NO_STATUS'))
PY
)
printf '%s\n' "$status" > SUBMISSION_STATUS.txt
printf '2026-07-23T00:00:00Z\tREMOTE_STAGE_A\tbuild_iccg_step1_pair.py\t%s\tRemote Stage-A preflight status.\n' "$status" >> RUN_HISTORY.tsv
if python audit_iccg_step1_pair.py preflight_report.json --submit-check; then
  if [[ -f pair.prmtop && -f LG1.inpcrd && -f LG2.inpcrd && -f LG1.in && -f LG2.in ]]; then
    sbatch run_iccg_step1_pair.sbatch
  else
    printf 'NOT_SUBMITTED_TOPOLOGY_NOT_BUILT\n' > SUBMISSION_STATUS.txt
    printf '2026-07-23T00:00:00Z\tREMOTE_SUBMIT\tsbatch run_iccg_step1_pair.sbatch\tNOT_SUBMITTED_TOPOLOGY_NOT_BUILT\tStage-A passed but Amber topology/input files are absent.\n' >> RUN_HISTORY.tsv
  fi
else
  exit 0
fi
