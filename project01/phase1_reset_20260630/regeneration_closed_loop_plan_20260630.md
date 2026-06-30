# Project 01 Phase1 Regeneration Closed-Loop Plan - 2026-06-30

> **For agentic workers:** keep this file as the controlling runbook for Phase1 until all 90/80/70/60/50 bins have 10 distinct sequences that pass every screening stage and are ready for QMMM. Do not stop a bin after one failed batch.

**Goal:** Generate sequence-diverse enzyme backgrounds at 90/80/70/60/50 identity until each similarity bin has 10 distinct sequences that pass all filters and can enter QMMM calculation.

**Architecture:** Run prediction on the GPU host and run all lightweight gates on the same host immediately after each output is written. Treat every stage as a closed loop: if a bin lacks enough passing samples, regenerate that bin with tighter active-pocket constraints, then repeat prediction, gate, PLACER, and RMSD evaluation.

**Tech Stack:** ESMFold/fair-esm + OpenFold on GPU for post-sequence structures; PLACER on GPU for ligand/crop conformers; project Python gate scripts for Kabsch/RMSD/key-distance evaluation; GitHub stores only markdown, scripts, TSV/JSON manifests, and summaries.

## Global Constraints

- Do not upload PDB ensembles, model weights, trajectories, large logs, or QMMM raw outputs to GitHub.
- Missing files, failed downloads, or not-yet-generated batches are `NOT_EVALUATED`, not `FAIL`.
- A sequence enters PLACER only after post-sequence entrance gate `PASS`.
- A PLACER conformer enters full-ligand completion only after crop ligand/key-distance gate `PASS`.
- A completed structure enters QMMM manifest only after strict all-atom ligand RMSD and reaction geometry gates pass.
- 50/60/70 bins must be regenerated if they do not reach the target pass count.
- A sequence is not accepted unless at least 10 of its 50 PLACER conformers pass the PLACER crop/RMSD gate.
- If a sequence fails the per-sequence PLACER threshold, discard that sequence for final accounting and generate a replacement sequence from the upstream sequence/backbone generation stage.

## Current Audited State

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

- 90/80/70 have at least one entrance-pass sample and can feed pilot PLACER runs.
- 70 is not complete; it needs regeneration until at least 10 entrance-pass sequences are available.
- 60/50 are evaluated FAIL under the current entrance gate, not merely missing. They must be regenerated with stronger active-pocket preservation.

## Final Acceptance Targets

Use these targets before Phase1 can be considered complete:

- `final_accepted_sequences_per_similarity_bin = 10`.
- `placer_samples_per_postseq_pass_sequence = 50`.
- `minimum_placer_crop_pass_conformers_per_sequence = 10`.
- `minimum_qmmm_ready_strict_pass_conformers_per_accepted_sequence = 1`.

## Threshold Provenance

Current strict gate values are project screening thresholds, not all direct literature thresholds. PLACER screening must inherit the same protein/pocket geometry gate used immediately before PLACER; ligand/reaction-pose checks are added on top of that inherited gate.

- `ligand_heavy_rmsd_max_A = 0.75` is a project-level conservative ligand-pose preservation threshold.
- It is motivated by the Baker serine-hydrolase work reporting sub-angstrom backbone/design agreement and active-site all-atom RMSD values around the sub-angstrom range for successful designs, but this exact `0.75 A` ligand heavy-atom RMSD cutoff has not been confirmed as a direct Baker paper threshold.
- The public serine-hydrolase design repository uses PLACER/ChemNet-style geometry analysis with RMSD, hydrogen-bond, angle, dihedral, and uncertainty features; it does not by itself establish `0.75 A` as the final ligand-RMSD hard cutoff for this project.
- Until the threshold is recalibrated against the paper's exact filter code or an internal validation set, report it as `PROJECT_STRICT_GATE`, not `BAKER_LITERATURE_GATE`.

Current threshold labels:

```text
global_backbone_rmsd_max_A = 2.50                  # PROJECT_ENTRANCE_GATE
fixed_backbone_rmsd_max_A = 1.00                   # PROJECT_ENTRANCE_GATE
catalytic_heavy_rmsd_max_A = 0.75                  # PROJECT_ENTRANCE_GATE
protein_key_distance_abs_delta_max_A = 0.75        # PROJECT_ENTRANCE_GATE
ligand_heavy_rmsd_max_A = 0.75                     # PROJECT_STRICT_GATE, not confirmed Baker hard cutoff
ligand_key_distance_abs_delta_max_A = 0.50         # PROJECT_STRICT_GATE
minimum_placer_crop_pass_conformers_per_sequence = 10 / 50  # PROJECT_ROBUSTNESS_GATE
```

PLACER gate consistency rule:

```text
PLACER_CROP_STRICT_PASS =
  inherited_postseq_protein_gate == PASS
  and ligand_reaction_pose_gate == PASS
```

where `inherited_postseq_protein_gate` uses the same global backbone, fixed-pocket backbone, catalytic heavy-atom, and protein key-distance thresholds as the postseq entrance gate.

If a bin misses any target, return to the nearest upstream generation stage:

- Entrance gate shortage: regenerate sequences/backbones for that bin.
- Per-sequence PLACER crop shortage: discard that sequence from final accounting and regenerate a replacement sequence. Do not count a sequence with fewer than 10 crop-pass conformers as accepted, even if one conformer looks promising.
- Full-ligand strict gate shortage: revise ligand completion/grafting, then rerun strict gate before QMMM.

Per-sequence acceptance rule:

```text
ACCEPT_SEQUENCE_FOR_QMMM_POOL =
  postseq_entrance_gate == PASS
  and PLACER_conformers_generated == 50
  and PLACER_crop_RMSD_gate_PASS_count >= 10
  and full_ligand_strict_QMMM_ready_count >= 1
```

Per-bin completion rule:

```text
COMPLETE_BIN =
  count(distinct ACCEPT_SEQUENCE_FOR_QMMM_POOL sequences in bin) >= 10
```

## Task 1: Regenerate Low-Pass Bins

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
90: request 1 additional entrance-pass sequence
80: request 2 additional entrance-pass sequences
70: request 6 additional entrance-pass sequences
60: request 10 additional entrance-pass sequences with stronger fixed-pocket constraints
50: request 10 additional entrance-pass sequences with stronger fixed-pocket constraints
```

For 60/50, increase candidate oversampling rather than relaxing the gate. Start with at least 5x oversampling per needed pass:

```text
60: generate at least 50 new candidates before postseq gate
50: generate at least 50 new candidates before postseq gate
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
  --queue-out manifests/placer_n50_entrance_pass_queue_round02.gpu.tsv \
  --gate-json-dir postseq_structure_gate/json_round02
```

Expected result:

```text
New rows with status PASS, FAIL, or NOT_EVALUATED_*.
No missing prediction is counted as FAIL.
```

## Task 2: PLACER for Entrance-Pass Samples

**Files:**

- Read: `/data/bht/project01_phase1_reset_gpu/manifests/placer_n50_entrance_pass_queue.gpu.tsv`
- Read: `/data/bht/project01_phase1_reset_gpu/manifests/placer_n50_entrance_pass_queue_round02.gpu.tsv`
- Create: `/data/bht/project01_phase1_reset_gpu/manifests/placer_n50_holo_inputs.round02.tsv`
- Create: `/data/bht/project01_phase1_reset_gpu/placer_runs/`

**Interfaces:**

- Consumes: entrance-pass samples only.
- Produces: PLACER model PDB/CSV files outside GitHub and lightweight manifest/status files for GitHub.

- [ ] **Step 1: Build holo inputs for every entrance-pass sample**

Append reference `bu2` HETATM/CONECT records to each entrance-pass protein prediction:

```text
ligand_policy = append_reference_bu2_refpose_before_PLACER
```

This is only an initial pose for PLACER. It is not a final ligand conformation.

- [ ] **Step 2: Run PLACER n=50 only for entrance-pass samples**

Use the reference Ser128-OG to bu2-C1 bond length:

```text
A-128-SER-OG:X-1-bu2-C1:1.533
```

Representative command:

```bash
cd /data/bht/project01_phase1_reset_gpu
CUDA_VISIBLE_DEVICES=0 /data/bht/design_tools/envs/placer_env/bin/python \
  /data/bht/design_tools/src/PLACER/run_PLACER.py \
  -f manifests/placer_n50_selected.ifile \
  -o placer_runs/selected_n50_bonded \
  -n 50 \
  --predict_ligand bu2 \
  --bonds A-128-SER-OG:X-1-bu2-C1:1.533 \
  --suffix bonded_bu2_n50 \
  --rerank prmsd
```

Expected result:

```text
One *_model.pdb and one *.csv per entrance-pass sample.
```

## Task 3: Crop Ligand RMSD and Reaction Geometry Gate

**Files:**

- Read: PLACER `*_model.pdb`
- Read: reference `/data/bht/project01_phase1_reset_gpu/manifests/denovo_SER_hydrolase_full_input.pdb`
- Create: `/data/bht/project01_phase1_reset_gpu/placer_crop_gate/placer_crop_ligand_gate.tsv`
- Create: `/data/bht/project01_phase1_reset_gpu/manifests/full_ligand_completion_queue.tsv`

**Interfaces:**

- Consumes: PLACER crop outputs.
- Produces: only crop-pass conformers for full-ligand completion.

- [ ] **Step 1: Evaluate ligand and reaction geometry before completion**

For every PLACER conformer, compute the same protein/pocket metrics used in the postseq entrance gate:

```text
global backbone RMSD after Kabsch alignment
fixed-pocket backbone RMSD
catalytic heavy-atom RMSD
protein key-distance deltas
```

Then compute the PLACER-added ligand/reaction metrics:

```text
bu2 heavy-atom RMSD vs reference pose after the same protein/pocket alignment
Ser128-OG to bu2-C1 distance delta
ligand-contact/oxyanion distances if present in the crop
```

Gate:

```text
inherited_postseq_protein_gate:
  global_backbone_rmsd_max_A = 2.50
  fixed_backbone_rmsd_max_A = 1.00
  catalytic_heavy_rmsd_max_A = 0.75
  protein_key_distance_abs_delta_max_A = 0.75

placer_ligand_addon_gate:
  ligand_heavy_rmsd_max_A = 0.75
  ser128_og_to_bu2_c1_abs_delta_max_A = 0.50
```

Expected labels:

```text
CROP_STRICT_PASS
FAIL_INHERITED_POSTSEQ_PROTEIN_GATE
FAIL_PLACER_LIGAND_ADDON_GATE
FAIL_BOTH_PROTEIN_AND_LIGAND_GATES
NOT_EVALUATED_PARSE_FAIL
```

- [ ] **Step 2: Apply the per-sequence PLACER acceptance threshold**

For each entrance-pass sequence, count crop-pass conformers across exactly 50 PLACER outputs:

```text
sequence_placer_crop_pass_count = count(CROP_STRICT_PASS conformers among 50)
```

Gate:

```text
if sequence_placer_crop_pass_count >= 10:
    sequence_status = PLACER_SEQUENCE_ACCEPTED_FOR_FULL_LIGAND_COMPLETION
else:
    sequence_status = DISCARD_SEQUENCE_REGENERATE_FROM_UPSTREAM
```

Do not rescue a sequence with 1-9 crop-pass conformers by sending those few conformers to QMMM. The final dataset requires 10 robust conformers per accepted sequence before full-ligand completion and QMMM preparation.

- [ ] **Step 3: If crop pass is insufficient, return to upstream generation**

If a sequence has fewer than `minimum_placer_crop_pass_conformers_per_sequence = 10`, do not complete ligand or launch QMMM for that sequence. Instead:

```text
mark sequence as DISCARD_SEQUENCE_REGENERATE_FROM_UPSTREAM
generate a replacement sequence for the same similarity bin
run postseq structure prediction and entrance gate for the replacement
only then run PLACER n=50 again
```

Additional PLACER resampling may be used only as a diagnostic before discarding the sequence. It does not change final accounting unless the sequence has a clean, documented 50-conformer PLACER run with at least 10 crop-pass conformers.

## Task 4: Full-Ligand Completion and QMMM Queue

**Files:**

- Read: `/data/bht/project01_phase1_reset_gpu/manifests/full_ligand_completion_queue.tsv`
- Create: `/data/bht/project01_phase1_reset_gpu/full_ligand_completed/`
- Create: `/data/bht/project01_phase1_reset_gpu/manifests/qmmm_ready_strict_pass.tsv`
- Create: `project01/phase1_reset_20260630/qmmm_ready_manifest_summary.md`

**Interfaces:**

- Consumes: crop-pass conformers from sequences with at least 10/50 PLACER crop-pass conformers.
- Produces: strict-pass full protein/full ligand QMMM manifest rows grouped by accepted sequence.

- [ ] **Step 1: Complete only crop-pass conformers**

Preserve PLACER crop heavy-atom coordinates during ligand completion. Record:

```text
crop_to_completed_preserved_atom_rmsd_A
full_ligand_heavy_rmsd_A
reaction_key_distance_deltas_A
```

- [ ] **Step 2: Generate QMMM rows only for strict pass**

QM region:

```text
complete bu2 ligand
Ser128 side-chain reaction fragment
His95 imidazole fragment
Ser126 side-chain fragment when included by gate/proton relay
```

MM region:

```text
remaining protein environment and fixed background residues
```

Status labels:

```text
READY_STRICT_PASS
PLACER_SEQUENCE_ACCEPTED_FOR_FULL_LIGAND_COMPLETION
DISCARD_SEQUENCE_REGENERATE_FROM_UPSTREAM
BLOCKED_CROP_GEOMETRY_FAIL
BLOCKED_FULL_LIGAND_RMSD_FAIL
NOT_EVALUATED
```

## Task 5: GitHub Synchronization

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
