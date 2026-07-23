#!/bin/bash
# Deploy scripts to the remote project directory and submit only after PASS preflight and topology audit.
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
stage_a_status=$(python - <<'PY'
import json
from pathlib import Path
print(json.loads(Path('preflight_report.json').read_text()).get('status','NOT_SUBMITTED_NO_STAGE_A_STATUS'))
PY
)
printf '%s\n' "$stage_a_status" > SUBMISSION_STATUS.txt
printf '2026-07-23T00:00:00Z\tREMOTE_STAGE_A\tbuild_iccg_step1_pair.py\t%s\tRemote Stage-A preflight status.\n' "$stage_a_status" >> RUN_HISTORY.tsv
if ! python audit_iccg_step1_pair.py preflight_report.json --submit-check; then
  exit 0
fi
python build_iccg_step1_pair.py \
  --iccg /work/home/acshdt1dks/iccg_sequence_design_20260723/inputs/ICCG_active_chainA.pdb \
  --lg1 /work/home/acshdt1dks/petase_qmmm_pilot8_20260721/inputs/reference/IsPETase_WT_LG1.pdb \
  --lg2 /work/home/acshdt1dks/petase_qmmm_pilot8_20260721/inputs/reference/IsPETase_WT_LG2.pdb \
  --out "$REMOTE_PROJECT" --stage-b \
  --ambertools-prefix /work/home/acshdt1dks/petase_qmmm_pilot8_20260721/software/ambertools20_data_f0ab9845_lf \
  --parmed-egg /public/software/apps/amber/2018/hpcx-2.4.1-gcc-7.3.1/lib/python2.7/site-packages/ParmEd-3.0.0_57.g74a84d30-py2.7-linux-x86_64.egg \
  --literature-prmtop /work/home/acshdt1dks/petase_orbmol_lg1_lg4_active_20260719/inputs/literature_rc_sources/acylation/vmd-md-b.prmtop \
  --literature-inpcrd /work/home/acshdt1dks/petase_orbmol_lg1_lg4_active_20260719/inputs/literature_rc_sources/acylation/vmd-md-b.inpcrd \
  --dftb-slko-path /work/home/acshdt1dks/petase_qmmm_pilot8_20260721/third_party/3ob-3-1
topology_status=$(python - <<'PY'
import json
from pathlib import Path
print(json.loads(Path('topology_audit.json').read_text()).get('status','NOT_SUBMITTED_NO_TOPOLOGY_STATUS'))
PY
)
printf '%s\n' "$topology_status" > SUBMISSION_STATUS.txt
printf '2026-07-23T00:00:00Z\tREMOTE_STAGE_B\tbuild_iccg_step1_pair.py --stage-b\t%s\tRemote topology audit status.\n' "$topology_status" >> RUN_HISTORY.tsv
if [[ "$topology_status" != "PASS" ]]; then
  exit 0
fi
sbatch run_iccg_step1_pair.sbatch
