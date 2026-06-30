# Project 01 Phase1 Sequence-Generation Plan - 2026-06-30

> **For agentic workers:** keep this file as the controlling runbook for the current Phase1 sequence-generation stage until all 90/80/70/60/50 bins have 10 distinct postseq entrance-pass sequences. Do not start QMMM in this stage.

**Goal:** Generate sequence-diverse enzyme backgrounds at 90/80/70/60/50 identity until each similarity bin has 10 distinct sequences that pass the postseq protein/pocket entrance gate.

**Architecture:** Run sequence/backbone generation, predict full protein structures, and run the postseq entrance gate immediately after each batch. Treat sequence generation as a closed loop: if a bin lacks enough postseq entrance-pass samples, regenerate that bin with stronger active-pocket constraints and repeat prediction/gating. Use two separated generation tracks: natural-scaffold conserved-site redesign and Baker-style active-site constrained de novo scaffold generation.

**Tech Stack:** ESMFold/fair-esm + OpenFold on GPU for post-sequence structures; project Python gate scripts for Kabsch/RMSD/key-distance evaluation; GitHub stores only markdown, scripts, TSV/JSON manifests, and summaries.

## Global Constraints

- Do not upload PDB ensembles, model weights, trajectories, large logs, or QMMM raw outputs to GitHub.
- Missing files, failed downloads, or not-yet-generated batches are `NOT_EVALUATED`, not `FAIL`.
- Current stage ends at sequence-panel generation; do not launch QMMM in this stage.
- PLACER is optional downstream diagnosis in this stage, not a hard filter for accepting sequences.
- A sequence is accepted for the current stage when its predicted full-protein structure passes the postseq entrance gate.
- 80/70/60 bins must be regenerated if they do not reach the target postseq entrance-pass count.
- Natural-scaffold regeneration must use MSA/conservation evidence to fix catalytic, pocket, shell, and conserved-core residues before generating low-similarity variants.
- De novo scaffold regeneration must start from fixed active-site/reactive geometry and generate new backbones around that motif; do not treat this as mutation of one fixed natural backbone.

## Initial Audited State Before Round02

Evaluated universe:

- 50 post-sequence structure targets: 10 each for 90/80/70/60/50 identity bins.
- All 50 now have GPU ESMFold predicted structures and postseq entrance-gate outputs.

Postseq entrance gate results:

| Bin | Evaluated | PASS | FAIL | Still needed to reach 10 PASS |
| --- | ---: | ---: | ---: | ---: |
| 90 | 10 | 9 | 1 | 1 |
| 80 | 10 | 8 | 2 | 2 |
| 70 | 10 | 4 | 6 | 6 |
| 60 | 10 | 0 | 10 | 10 |
| 50 | 10 | 0 | 10 | 10 |

Interpretation:

- 90/80/70 already contain entrance-pass samples but still need enough rows to complete the 10-per-bin sequence panel.
- In the initial all50 batch, 60/50 were evaluated FAIL under the entrance gate, not merely missing.
- The completed PLACER pilot is diagnostic only. `CROP_STRICT_PASS = 0/300` is not a reason to reject postseq entrance-pass sequences in this stage, because TS-like conformers are not expected to behave like ligand low-energy poses.

## Current Audited State After Round02/Round03/Round04/Round05

Updated GPU status:

| Batch | Evaluated | ESMFold OK | Entrance PASS | Interpretation |
| --- | ---: | ---: | ---: | --- |
| round02 controlled mutation-count candidates | 164 | 164 | 7 | PASS only in 90%; random low-identity mutations drifted the pocket. |
| round02d actual-bin LigandMPNN refilter | 66 | 66 | 18 | Cooperative LigandMPNN designs restored 70/60/50 passes. |
| round03 actual-bin production | 152 | 152 | 33 | Completed 80/70 and added one 60. |
| round04 empirical 60 subsets | 144 | 144 | 5 | Used empirical 60-pass mutation position sets. |
| round05 template2 60 subset | 72 | 72 | 8 | Completed 60 bin. |

Current combined accepted distinct sequence-panel counts:

| Bin | Accepted distinct PASS | Target | Still needed |
| --- | ---: | ---: | ---: |
| 90 | 16 | 10 | 0 |
| 80 | 11 | 10 | 0 |
| 70 | 38 | 10 | 0 |
| 60 | 16 | 10 | 0 |
| 50 | 11 | 10 | 0 |

The current sequence-generation target is complete. Do not generate more variants for the denovo_SER_hydrolase reference unless replacement rows are needed after manual review.

Use `/data/bht/project01_phase1_reset_gpu/manifests/postseq_entrance_pass_sequence_panel_capped10.tsv` as the 10-per-bin current-stage panel.

## Current Stage Acceptance Targets

Use these targets before the current sequence-generation stage can be considered complete:

- `postseq_entrance_pass_sequences_per_similarity_bin = 10`.
- `accepted_sequence_for_current_stage = postseq_entrance_gate == PASS`.
- `defer_placer_filtering = true`.
- `defer_qmmm_calculation = true`.

## Threshold Provenance

Current strict gate values are project screening thresholds, not all direct literature thresholds. The active gate for this stage is the postseq protein/pocket entrance gate. PLACER ligand/reaction-pose filters are retained only as downstream diagnostic records and must not be used to reject sequences in the current sequence-generation stage.

Current threshold labels:

```text
global_backbone_rmsd_max_A = 2.50                  # PROJECT_ENTRANCE_GATE
fixed_backbone_rmsd_max_A = 1.00                   # PROJECT_ENTRANCE_GATE
catalytic_heavy_rmsd_max_A = 0.75                  # PROJECT_ENTRANCE_GATE
protein_key_distance_abs_delta_max_A = 0.75        # PROJECT_ENTRANCE_GATE
ligand_heavy_rmsd_max_A = DEFERRED                 # not used in current sequence-generation stage
ligand_key_distance_abs_delta_max_A = DEFERRED     # not used in current sequence-generation stage
minimum_placer_crop_pass_conformers_per_sequence = DEFERRED
```

Do not use PLACER failure as a sequence-rejection rule in this stage. The immediate question is whether the protein background preserves the intended active pocket strongly enough to later test explicit TS conformer stability.

If a bin misses any target, return to the nearest upstream generation stage:

- Entrance gate shortage: regenerate sequences/backbones for that bin.
- PLACER crop shortage: record as downstream diagnostic only; do not regenerate solely because of PLACER failure in this stage.
- QMMM shortage: out of scope for this stage.

Per-sequence acceptance rule for the current stage:

```text
ACCEPT_SEQUENCE_FOR_SEQUENCE_PANEL =
  postseq_entrance_gate == PASS
```

Per-bin completion rule:

```text
COMPLETE_BIN =
  count(distinct ACCEPT_SEQUENCE_FOR_SEQUENCE_PANEL sequences in bin) >= 10
```

## Task 1: Regenerate Low-Pass Bins With Two Tracks

### Track A: natural-scaffold conserved-site redesign

For natural scaffold work, first build or import an MSA for the chosen enzyme family/scaffold. Use it to classify residues:

```text
FIXED = catalytic residues + ligand/direct-contact residues + MSA-conserved core + pocket-shell positions that preserve key geometry
MUTABLE = nonconserved background residues, expanded cautiously after gate evidence
```

This track should now be used for new natural-scaffold work, not to fill the completed denovo_SER_hydrolase panel. If a natural scaffold is selected, first build the MSA/conservation map and then generate a separate manifest.

```text
current serine-hydrolase existing/reference-scaffold sequence panel: COMPLETE`r`ntrue active-site constrained de novo/new-backbone panel: NOT_STARTED`r`nnew natural-scaffold panel: NOT_STARTED
```

The current seed/queue files are:

```text
project01/phase1_reset_20260630/phase1_dual_track_seed_manifest.tsv
project01/phase1_reset_20260630/phase1_next_generation_queue_20260630.md
project01/phase1_reset_20260630/natural_scaffold_msa_generation_protocol_20260630.md
```

### Track B: active-site constrained de novo scaffold generation

For Baker-style serine-hydrolase diversity, fixed active-site/reactive geometry comes before backbone generation. This true new-backbone Route B is not yet complete; the completed capped10 panel only covers the existing/reference-scaffold sequence route. Generate new backbones around the motif, then design sequences for those backbones. Keep these candidates in a separate manifest until the downstream comparison between natural and de novo backgrounds is explicitly defined.

**Files:**

- Read: `/data/bht/project01_phase1_reset_gpu/postseq_structure_gate/tables/all50_entrance_gate.tsv`
- Read: `/data/bht/project01_phase1_reset_gpu/manifests/postseq_structure_targets.gpu_local.tsv`
- Create: `/data/bht/project01_phase1_reset_gpu/manifests/regeneration_queue_round02.tsv`
- Create: `project01/phase1_reset_20260630/regeneration_round02_summary.md`

**Interfaces:**

- Consumes: bin-level pass counts from `all50_entrance_gate.tsv`.
- Produces: a queue of bins and requested additional sequences for the next generation round.

- [ ] **Step 1: Build the round-02 regeneration queue**

Run on GPU host:

```bash
cd /data/bht/project01_phase1_reset_gpu
/data/bht/design_tools/envs/rfaa_venv/bin/python - <<'PY'
import csv
from pathlib import Path

gate = Path("postseq_structure_gate/tables/all50_entrance_gate.tsv")
target = 10
counts = {str(b): {"PASS": 0, "FAIL": 0, "NOT_EVALUATED": 0} for b in [90, 80, 70, 60, 50]}
with gate.open() as f:
    for row in csv.DictReader(f, delimiter="\t"):
        status = row["status"]
        bucket = status if status in {"PASS", "FAIL"} else "NOT_EVALUATED"
        counts[row["bin"]][bucket] += 1

out = Path("manifests/regeneration_queue_round02.tsv")
with out.open("w", newline="") as f:
    fields = ["bin", "current_pass", "current_fail", "current_not_evaluated", "target_pass", "new_sequences_requested", "policy"]
    writer = csv.DictWriter(f, fields, delimiter="\t")
    writer.writeheader()
    for bin_id in ["90", "80", "70", "60", "50"]:
        current_pass = counts[bin_id]["PASS"]
        needed = max(0, target - current_pass)
        if needed == 0:
            policy = "NO_REGENERATION_NEEDED"
        elif bin_id in {"60", "50"}:
            policy = "REGENERATE_WITH_STRONGER_ACTIVE_POCKET_CONSTRAINTS"
        else:
            policy = "REGENERATE_UNTIL_TARGET_PASS_COUNT"
        writer.writerow({
            "bin": bin_id,
            "current_pass": current_pass,
            "current_fail": counts[bin_id]["FAIL"],
            "current_not_evaluated": counts[bin_id]["NOT_EVALUATED"],
            "target_pass": target,
            "new_sequences_requested": needed,
            "policy": policy,
        })
print(out)
PY
```

Expected output:

```text
manifests/regeneration_queue_round02.tsv
```

- [ ] **Step 2: Generate additional sequences only for bins with shortage**

Use the existing fixed-active-pocket generation workflow, but set bin-specific sequence identity targets:

```text
90: complete; do not generate unless replacement is needed
80: complete; do not generate unless replacement is needed
70: complete; do not generate unless replacement is needed
60: complete; do not generate unless replacement is needed
50: complete; do not generate unless replacement is needed
```

If replacement rows are ever needed, use empirical tolerated-position redesign rather than relaxing the gate:

```text
do_not_relax_postseq_entrance_gate = true
prefer_empirical_tolerated_position_redesign_for_60 = true
```

- [ ] **Step 3: Predict and gate each new batch immediately on GPU**

For each new batch, write new predicted PDB paths under:

```text
/data/bht/project01_phase1_reset_gpu/postseq_structure_models_round02/
```

Then run:

```bash
cd /data/bht/project01_phase1_reset_gpu
/data/bht/design_tools/envs/rfaa_venv/bin/python scripts/run_postseq_gate_manifest.py \
  --manifest manifests/postseq_structure_targets_round02.gpu_local.tsv \
  --reference manifests/denovo_SER_hydrolase_full_input.pdb \
  --fixed-residues manifests/fixed_active_or_direct_contact_residues.txt \
  --out-manifest manifests/postseq_entrance_gate_manifest_round02.tsv \
  --queue-out manifests/postseq_entrance_pass_sequence_panel_round02.gpu.tsv \
  --gate-json-dir postseq_structure_gate/json_round02
```

Expected result:

```text
New rows with status PASS, FAIL, or NOT_EVALUATED_*.
No missing prediction is counted as FAIL.
```

## Task 2: Maintain the Sequence Panel Manifest

**Files:**

- Read: `/data/bht/project01_phase1_reset_gpu/postseq_structure_gate/tables/all50_entrance_gate.tsv`
- Read: `/data/bht/project01_phase1_reset_gpu/manifests/postseq_entrance_gate_manifest_round02.tsv`
- Create: `/data/bht/project01_phase1_reset_gpu/manifests/postseq_entrance_pass_sequence_panel.tsv`
- Create: `project01/phase1_reset_20260630/sequence_panel_status_round02.md`

**Interfaces:**

- Consumes: postseq entrance-gate rows from all completed batches.
- Produces: one lightweight manifest of accepted sequences for the current stage.

- [ ] **Step 1: Merge entrance-pass rows across batches**

For every evaluated batch, append rows where `status == PASS` to:

```text
/data/bht/project01_phase1_reset_gpu/manifests/postseq_entrance_pass_sequence_panel.tsv
```

The manifest must include at least:

```text
sample_id
bin
sequence
sequence_identity_target
predicted_pdb
postseq_gate_json
global_backbone_rmsd_A
fixed_backbone_rmsd_A
catalytic_heavy_rmsd_A
max_abs_protein_key_delta_A
status
```

- [ ] **Step 2: Stop when each bin has 10 sequence-panel rows**

Completion criterion:

```text
for each bin in 90/80/70/60/50:
    count(distinct sample_id where status == PASS) >= 10
```

If any bin remains below 10, return to Task 1 Step 2 and regenerate more candidates for that bin.

## Task 3: Deferred PLACER/TS/QMMM Work

**Files:**

- Read later: `/data/bht/project01_phase1_reset_gpu/manifests/postseq_entrance_pass_sequence_panel.tsv`
- Create later: TS conformer and QMMM manifests outside the current stage.

**Interfaces:**

- Consumes: accepted sequence panel after the current stage is complete.
- Produces: future TS-stability/QMMM inputs.

- [ ] **Step 1: Do not launch QMMM during current sequence generation**

Current decision:

```text
QMMM_OUT_OF_SCOPE_FOR_CURRENT_STAGE = true
PLACER_REQUIRED_FOR_SEQUENCE_ACCEPTANCE = false
```

- [ ] **Step 2: Keep previous PLACER pilot as diagnostic only**

The completed PLACER pilot showed `CROP_STRICT_PASS = 0/300`, but this is not a sequence-generation blocker. Interpretation:

```text
TS conformers are not expected to be ligand low-energy poses.
The next scientific question is whether explicit TS-like conformer ensembles can be stabilized by accepted protein backgrounds.
Do not regenerate a sequence solely because PLACER failed to reproduce the ligand reference pose.
```

Future downstream work should start from the accepted sequence panel and define an explicit TS conformer ensemble or constrained TS-like geometry before QMMM.

## Task 4: GitHub Synchronization



**Files:**

- Modify: `project01/phase1_reset_20260630/phase1_reset_status_20260630.md`
- Modify or create: `project01/phase1_reset_20260630/*_summary.md`
- Modify or create: lightweight TSV/JSON manifests only when they are small.

- [ ] **Step 1: Sync only lightweight evidence**

Commit:

```text
markdown status files
gate summary TSV/JSON files
small manifest files
scripts
```

Do not commit:

```text
PDB files
PLACER model outputs
model checkpoints
trajectory files
large logs
QM/MM raw outputs
```

- [ ] **Step 2: Keep denominators explicit in every status update**

Every update must include:

```text
evaluated universe
PASS count
FAIL count
NOT_EVALUATED count
thresholds used
next upstream action if target pass count is not met
```

