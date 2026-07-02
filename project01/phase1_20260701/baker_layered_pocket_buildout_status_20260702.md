# Project 01 Baker Layered Pocket Build-Out Reset - 2026-07-02

## Goal

Regenerate the Baker serine-hydrolase sequence campaign with a layered pocket-first strategy instead of relying on one full scaffold (`sample_1000_refined_0`) or one-shot full-scaffold generation only.

The new route is:

1. L0 fixed catalytic unit: Baker `theozyme.pdb` with `bn1`, catalytic Ser/His/Asp motif, oxyanion-hole motif, and motif sequence supplied.
2. L1 pocket-support scaffold sweep: generate multiple scaffold families around the fixed catalytic/theozyme unit using shorter or varied contigs.
3. L1 gate: count only present PDBs as evaluated; absent future outputs are NOT_EVALUATED. PASS requires motif CA RMSD <= 1.0 A, motif pair max delta <= 1.0 A, and ligand records present.
4. L2 extension/refinement: only L1 PASS structures are used for further build-out/refinement. This stage will be defined after L1 PASS structures exist.
5. Sequence regeneration: LigandMPNN generates larger multi-round sequence panels from multiple L1/L2 PASS scaffolds, not from one scaffold only.
6. Sequence bins: keep the project bins 90/80/70/60/50, but record per-parent scaffold, per-round seed, candidate count, selected count, and downstream gate status separately.

## First Layer Contig Families

The first L1 sweep intentionally tries several contig scales. It does not assume compact contigs will work; it measures them.

- `compact`: `6,A56-60,18,A83-85,8,A113-115,32,B145-147,6`
- `medium`: `8,A56-60,24,A83-85,10,A113-115,48,B145-147,8`
- `near_original`: `10,A56-60,30,A83-85,12,A113-115,60,B145-147,8`

The original one-shot contig was `12,A56-60,36,A83-85,15,A113-115,73,B145-147,10`. The new sweep asks whether shorter/support-focused scaffolds can preserve the catalytic geometry before sequence design.

## Current Precheck Evidence

The old batch50 run is still useful as a diagnostic, but it is not the final layered route. A precheck gate on the 7 PDB files already present in `outputs/ca_rfd_baker_theozyme_formal_constraints_batch50_20260702` gave `PASS=7/7` under the motif CA gate. This means the public checkpoint can preserve the minimal theozyme/motif geometry in generated backbones; the earlier failure was downstream sequence/structure validation, not immediate motif-generation failure.

Remote precheck files:

- `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_theozyme_formal_constraints_batch50_20260702_motif_gate_layered_precheck.tsv`
- `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_theozyme_formal_constraints_batch50_20260702_motif_gate_layered_precheck_summary.json`

## GPU Launch Policy

Do not start a second CA_RFDiffusion job while another compute process is saturating the GPU unless explicitly requested. The L1 launcher writes a BLOCKED status if GPU compute processes exist or utilization is above threshold.

## Success Criteria For This Reset

This reset is successful only when evidence files show:

- L1 evaluated denominator and PASS/FAIL/NOT_EVALUATED counts per contig family.
- At least one contig family has enough L1 PASS structures to justify L2 or sequence generation.
- Sequence panels are regenerated from multiple L1/L2 PASS scaffolds with per-bin counts, not from one scaffold only.
- Downstream reports distinguish sequence pool count from structure-gated accepted count.

## 2026-07-02 Current Execution Update

The active target is now explicitly: regenerate sequences only after motif-gated layered scaffolds exist, using larger multi-round candidate pools and repeated attempts until each 90/80/70/60/50 bin has enough accepted sequences for downstream structural evaluation.

Sequence-generation policy after L1/L2 PASS scaffolds:

- Parent universe: PASS rows from an L1/L2 motif gate TSV; non-existent future outputs are NOT_EVALUATED, not FAIL.
- Multi-scaffold input: use up to 5 parent scaffolds per launch by default, not one scaffold.
- Multi-round attempts: run at least 2 seed rounds by default; increase rounds if a bin lacks enough accepted sequences.
- Candidate pool size per parent and round: select up to 100 records for 90% and 200 records for each of 80/70/60/50 before downstream structure checks.
- Identity targets are computed from each parent scaffold's actual sequence length. If the scaffold has fewer mutable residues than requested, the target mutation count is capped at the actual redesignable positions and recorded as `target_mutations`; the original request is recorded as `requested_mutations`.
- The practical acceptance target remains 10 final accepted sequences per bin, but the raw generated sequence pool is intentionally larger so failed structure gates can be replaced by more attempts.

Remote evidence from `bht@192.168.10.38:/data/bht/project01_baker_serhyd_routeB_20260701` at `2026-07-02T17:15:29+08:00`:

- GPU state: utilization 100%, memory 30403 MiB / 81920 MiB, compute process count 10.
- Existing old batch50 process: PID 413631 still running, elapsed 05:10:13 at check time.
- Existing old batch50 PDB count: 8 present PDB files.
- L1 medium launch attempt: `CONTIG_SET=medium NUM_DESIGNS=20 DESIGN_STARTNUM=2000 scripts/launch_baker_layered_l1_contig_sweep.sh` wrote `BLOCKED_GPU_COMPUTE_PROCESS_PRESENT`; no new L1 job was started.
- L1 status file: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_layered_l1_medium_publicckpt_20260702_status.json`.
- Dynamic original-contig precheck on present old batch50 outputs: evaluated 8 PDB files, counts `PASS=8`, status `DONE`.
- Dynamic precheck summary: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_theozyme_formal_constraints_batch50_20260702_motif_gate_layered_dynamic_precheck_summary.json`.

Remote script verification:

- `scripts/gate_ca_rfdiffusion_theozyme_motif.py` supports dynamic motif maps through `--contig`.
- `scripts/generate_baker_layered_multiscaffold_ligandmpnn_bins.py` is installed on the GPU host, executable, UTF-8 without BOM, and passed `py_compile` together with the gate script.
- The generator is prepared but should not be run against the old original-contig precheck unless explicitly choosing that as a fallback. The intended next run is from L1/L2 PASS scaffold gate TSVs.

## 2026-07-02 L1 Queue Driver Update

A lightweight L1 queue driver was added so repeated heartbeats can progress through multiple pocket-support contig attempts instead of repeatedly trying only `medium`.

Queue policy:

- Default contig order: `compact medium near_original`.
- Default designs per contig: `NUM_DESIGNS=20`.
- Default start numbers: `compact=3000`, `medium=4000`, `near_original=5000`.
- Fixed run IDs are used for resumability: `ca_rfd_baker_layered_l1_<contig_set>_publicckpt_20260702`.
- If a run is already `LAUNCHED` and its PID is alive, the queue writes `MONITORING` and does not start a duplicate.
- If a run has PDB outputs, the queue runs the dynamic motif gate and writes `<run_id>_motif_gate.tsv` plus `<run_id>_motif_gate_summary.json`.
- If the GPU is busy, the queue records `BLOCKED_*` and exits cleanly.

Remote verification at `2026-07-02T17:22:32+08:00`:

- Remote script: `/data/bht/project01_baker_serhyd_routeB_20260701/scripts/run_baker_layered_l1_queue.sh`.
- Syntax check: `bash -n scripts/run_baker_layered_l1_queue.sh` passed.
- Queue execution: `CONTIG_SETS='compact medium near_original' NUM_DESIGNS=20 START_BASE=3000 scripts/run_baker_layered_l1_queue.sh`.
- Queue status: `BLOCKED_GPU_COMPUTE_PROCESS_PRESENT`.
- Current attempted set: `compact`.
- Queue status file: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_layered_l1_queue_status.json`.
- Compact launcher status file: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_layered_l1_compact_publicckpt_20260702_status.json`.

The active heartbeat now calls the queue driver first. It will move from `compact` to `medium` to `near_original` as each run either produces gateable PDBs or remains unevaluated due to GPU blocking.

## 2026-07-02 New-Route-Only Switch

The user clarified that previous Project 01 routes can stop and only the new layered route should run.

Remote process action at `2026-07-02T17:26:59+08:00`:

- Matched old Project 01 process: PID `413631`, old original-contig `ca_rfd_baker_theozyme_formal_constraints_batch50_20260702` batch50 run.
- Action: sent `TERM` with `kill 413631`.
- Verification: `ps -p 413631` returned no process; PID `413631` disappeared from `nvidia-smi --query-compute-apps`.
- Other GPU Python processes `209606` and `523038` were inspected and belong to other users/directories, so they were not touched.

New layered route launch:

- Command: `FORCE=1 ALLOW_SHARED_GPU=1 CONTIG_SET=compact NUM_DESIGNS=20 DESIGN_STARTNUM=3000 scripts/launch_baker_layered_l1_contig_sweep.sh`.
- Result: `LAUNCHED`.
- Run ID: `ca_rfd_baker_layered_l1_compact_publicckpt_20260702`.
- PID: `555939`.
- Log: `/data/bht/project01_baker_serhyd_routeB_20260701/logs/ca_rfd_baker_layered_l1_compact_publicckpt_20260702.log`.
- Output prefix: `/data/bht/project01_baker_serhyd_routeB_20260701/outputs/ca_rfd_baker_layered_l1_compact_publicckpt_20260702/sample`.
- Status file: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_layered_l1_compact_publicckpt_20260702_status.json`.
- Queue status after launch: `MONITORING`, current set `compact`, PID `555939`.
- PDB count immediately after launch: 0, so no L1 motif gate has been evaluated yet.

Heartbeat policy was updated accordingly: do not restart old original-contig batch50; monitor/advance only the new L1 queue with `FORCE=1 ALLOW_SHARED_GPU=1`.
