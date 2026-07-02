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
## 2026-07-02 Compact L1 Running Check

Remote check at `2026-07-02T17:36:51+08:00`:

- Active new-route PID: `555939`.
- Run ID: `ca_rfd_baker_layered_l1_compact_publicckpt_20260702`.
- State: running (`STAT=Rl`), elapsed `08:58`, CPU about `95.9%`, RSS about `763808 KB`.
- GPU evidence: PID `555939` is present in `nvidia-smi --query-compute-apps` and uses about `1056 MiB`.
- Output directory: `/data/bht/project01_baker_serhyd_routeB_20260701/outputs/ca_rfd_baker_layered_l1_compact_publicckpt_20260702`.
- Current compact L1 PDB count: `0`.
- Current log still shows `Making design ... sample_3000` and checkpoint load, with no traceback observed.

Screening/evaluation status:

- Evaluated universe for compact L1 motif gate: `0` PDB files, because no compact L1 PDB has been written yet.
- PASS: not evaluated.
- FAIL: not evaluated.
- NOT_EVALUATED: future compact L1 outputs that do not exist yet; these are not failures.

Next action remains: keep monitoring PID `555939`; when at least one compact L1 PDB appears, run the dynamic motif gate and then let the queue advance to `medium` / `near_original` or to sequence generation if enough PASS scaffolds exist.
## 2026-07-02 Compact L1 Denoising Progress

Remote check at `2026-07-02T17:42:50+08:00`:

- Active new-route PID: `555939`.
- Run ID: `ca_rfd_baker_layered_l1_compact_publicckpt_20260702`.
- State: running (`STAT=Rl`), elapsed `14:58`, CPU about `114%`, RSS about `1894260 KB`.
- Old original-contig batch50 process remains absent.
- Current compact L1 PDB count: `0`.
- Log progress: model loading completed and denoising is active for `sample_3000`; observed progress from `t=50` to `t=45`, about `5/50` denoising steps complete.
- Approximate observed step time: about 43-45 seconds per denoising step under shared GPU load.

Screening/evaluation status:

- Evaluated universe for compact L1 motif gate remains `0` PDB files.
- PASS: not evaluated.
- FAIL: not evaluated.
- NOT_EVALUATED: compact L1 outputs that have not been written yet; these are not failures.

Next action: keep monitoring PID `555939`; run the compact L1 dynamic motif gate immediately after any PDB appears, even before all 20 designs complete.
## 2026-07-02 Compact L1 Denoising t40 Check

Remote check at `2026-07-02T17:46:43+08:00`:

- Active new-route PID: `555939`.
- Run ID: `ca_rfd_baker_layered_l1_compact_publicckpt_20260702`.
- State: running (`STAT=Rl`), elapsed `18:36`, CPU about `115%`, RSS about `1894952 KB`.
- Old original-contig batch50 process remains absent.
- Current compact L1 PDB count: `0`.
- Log progress: denoising advanced to `t=40`, about `10/50` denoising steps complete for `sample_3000`.
- Queue status refreshed to `MONITORING` at `2026-07-02T17:46:43`.

Screening/evaluation status remains strict:

- Evaluated universe for compact L1 motif gate: `0` PDB files.
- PASS: not evaluated.
- FAIL: not evaluated.
- NOT_EVALUATED: compact L1 outputs not yet written.

Next action: continue monitoring until the first compact L1 PDB appears, then run dynamic motif gate on the present PDB universe immediately.
## 2026-07-02 Queue Early-Gate Fix

A queue-driver logic issue was fixed: the previous queue checked whether a CA_RFDiffusion PID was alive before checking whether any PDB had already been written. That could delay motif gating until all 20 compact designs finished.

New queue behavior:

- For each contig run, count present PDB files first.
- If `pdb_count > 0`, run the dynamic motif gate immediately on the present PDB universe.
- If the CA_RFDiffusion PID is still alive after gating, write queue status `GATED_RUNNING` and keep monitoring/regating as more PDBs appear.
- If `pdb_count = 0` and the PID is alive, write `MONITORING` with `pdb_count=0`.

Remote verification at `2026-07-02T17:56:59+08:00`:

- Remote script: `/data/bht/project01_baker_serhyd_routeB_20260701/scripts/run_baker_layered_l1_queue.sh`.
- `bash -n scripts/run_baker_layered_l1_queue.sh` passed.
- Queue execution: `FORCE=1 ALLOW_SHARED_GPU=1 CONTIG_SETS='compact medium near_original' NUM_DESIGNS=20 START_BASE=3000 scripts/run_baker_layered_l1_queue.sh`.
- Current queue status: `MONITORING`.
- Current set: `compact`.
- Current PID: `555939`.
- Current `pdb_count`: `0`.
- Next action in queue status: `wait_for_current_l1_to_write_pdb_or_finish`.

This preserves the strict evaluated/not-evaluated rule: no PDB means compact L1 motif gate evaluated universe remains `0`, not FAIL.
## 2026-07-02 New-Route-Only Confirmation

Remote check at `2026-07-02T18:00:30+08:00`:

- Old original-contig batch50 route remains stopped; no matching `formal_constraints_batch50` Project 01 process was present.
- Active new-route PID remains `555939`.
- Run ID: `ca_rfd_baker_layered_l1_compact_publicckpt_20260702`.
- State: running (`STAT=Rl`), elapsed `32:37`, CPU about `115%`, RSS about `1897496 KB`.
- Current compact L1 PDB count: `0`.
- Log progress: denoising advanced to `t=21`, about `29/50` denoising steps complete for `sample_3000`.
- Queue command rerun with the new-route-only settings:
  `FORCE=1 ALLOW_SHARED_GPU=1 CONTIG_SETS='compact medium near_original' NUM_DESIGNS=20 START_BASE=3000 scripts/run_baker_layered_l1_queue.sh`.
- Queue status remains `MONITORING`, current set `compact`, current PID `555939`, `pdb_count=0`, next action `wait_for_current_l1_to_write_pdb_or_finish`.

Current interpretation:

- Previous Project 01 routes are no longer being advanced.
- Only the layered pocket-first compact L1 route is active.
- Evaluated universe for compact L1 motif gate is still `0` PDB files, so there is no PASS/FAIL result yet.

Next action: continue monitoring PID `555939`; as soon as the first compact L1 PDB is written, run the dynamic motif gate on that present output universe and record PASS/FAIL counts.
## 2026-07-02 First Compact L1 Gate PASS

Remote check at `2026-07-02T18:18:09+08:00`:

- Active new-route PID remains `555939`; the same compact L1 run is still generating additional designs.
- First compact L1 output appeared:
  `/data/bht/project01_baker_serhyd_routeB_20260701/outputs/ca_rfd_baker_layered_l1_compact_publicckpt_20260702/sample_3000.pdb`.
- Output sizes: `sample_3000.pdb` about `31230` bytes; `sample_3000.trb` about `23508` bytes.
- Queue script ran the dynamic motif gate immediately on the present PDB universe.
- Queue status changed to `GATED_RUNNING`, with `pdb_count=1` and next action `continue_monitoring_current_l1_and_regate_present_pdbs`.

Compact L1 motif gate result:

- Evaluated universe: `1` present compact L1 PDB file; absent future samples remain `NOT_EVALUATED`.
- PASS: `1` (`sample_3000`).
- FAIL: `0` among evaluated outputs.
- Gate definition: 14 motif CA Kabsch RMSD `<= 1.0 A`, max pair-distance delta `<= 1.0 A`, and ligand `bn1` present.
- `sample_3000` metrics:
  - motif CA RMSD: `0.0776 A`.
  - max pair-distance delta: `0.2107 A`.
  - mean pair-distance delta: `0.0547 A`.
  - ligand atom records: `22`.

Evidence files on GPU host:

- Summary JSON: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_layered_l1_compact_publicckpt_20260702_motif_gate_summary.json`.
- Gate TSV: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_layered_l1_compact_publicckpt_20260702_motif_gate.tsv`.
- Queue status: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_layered_l1_queue_status.json`.

Next action: keep the compact L1 run active until more PDBs accumulate; re-run the dynamic gate on the present PDB universe and then let the queue advance to the next layered contig set according to `scripts/run_baker_layered_l1_queue.sh`.
## 2026-07-02 First PASS Scaffold Sequence Panel

Remote check at `2026-07-02T18:26:18+08:00`:

- Sequence generation was started from the first motif-gated compact L1 PASS scaffold, `sample_3000`.
- Parent gate evidence: compact L1 gate TSV `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_layered_l1_compact_publicckpt_20260702_motif_gate.tsv`, where `sample_3000` is `PASS`.
- Sequence panel run ID: `baker_layered_multiscaffold_ligandmpnn_from_compact_sample3000_20260702_182250`.
- Sequence generation route: `generate_baker_layered_multiscaffold_ligandmpnn_bins.py` using parent-limit `1`, rounds `2`, bins `90/80/70/60/50`.
- Fixed motif positions in generated scaffold numbering: `14` residues from contig `6,A56-60,18,A83-85,8,A113-115,32,B145-147,6`.
- Parent scaffold sequence length: `84` residues.
- Status: `DONE`, started `2026-07-02T18:22:50`, finished `2026-07-02T18:25:17`.

Selected sequence counts by target identity bin:

- `90`: `23` selected sequences.
- `80`: `212` selected sequences.
- `70`: `350` selected sequences.
- `60`: `399` selected sequences.
- `50`: `400` selected sequences.
- Total selected: `1384` sequences.

Per-round selected counts:

- Round 01: 90=`12`, 80=`109`, 70=`177`, 60=`199`, 50=`200`.
- Round 02: 90=`11`, 80=`103`, 70=`173`, 60=`200`, 50=`200`.

Audit interpretation:

- Evaluated universe for this sequence panel is one motif-gated parent scaffold, `sample_3000`.
- These are LigandMPNN sequence-layer candidates that pass the script filter for fixed motif/no-native constraints and target mutation binning.
- These are not yet structure-gated final candidates. They still need structure prediction and pocket/motif geometry screening before being counted as final per-bin accepted designs.
- The 90% bin has fewer than the nominal `100` selected target because the filter found only `23` passing records across two rounds for this compact parent.

Remote evidence files:

- Summary JSON: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_layered_multiscaffold_ligandmpnn_from_compact_sample3000_20260702_182250_summary.json`.
- Selected sequence TSV: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_layered_multiscaffold_ligandmpnn_from_compact_sample3000_20260702_182250_selected.tsv`.
- Sequence design directory: `/data/bht/project01_baker_serhyd_routeB_20260701/sequence_design/baker_layered_multiscaffold_ligandmpnn_from_compact_sample3000_20260702_182250`.

Current scaffold generation status remains active:

- Compact L1 CA_RFDiffusion PID `555939` is still running.
- Current compact L1 evaluated PDB universe remains `1` PDB; PASS=`1`, FAIL=`0`, future compact outputs are `NOT_EVALUATED`.
- Queue status remains `GATED_RUNNING` and will re-gate as additional compact L1 PDBs are written.

Next action: keep monitoring compact L1; when more compact PASS scaffolds appear, rerun/extend multiscaffold LigandMPNN panels and then advance to structure prediction plus motif/pocket gate for selected sequences.
## 2026-07-02 Compact L1 Second PASS and Two-Parent Sequence Panel

Remote check at `2026-07-02T18:53:06+08:00` and sequence panel check at `2026-07-02T19:01:29+08:00`:

Compact L1 scaffold gate:

- Active compact L1 CA_RFDiffusion PID remains `555939` and is continuing to generate later samples.
- Current compact L1 output PDB universe: `2` present PDB files, `sample_3000` and `sample_3001`.
- Evaluated universe: those `2` present PDB files only; future compact outputs remain `NOT_EVALUATED`.
- PASS: `2`.
- FAIL: `0` among evaluated outputs.
- Gate definition: 14 motif CA Kabsch RMSD `<= 1.0 A`, max pair-distance delta `<= 1.0 A`, ligand `bn1` present.
- `sample_3000`: motif CA RMSD `0.0776 A`, max pair-distance delta `0.2107 A`, mean pair-distance delta `0.0547 A`, ligand atom records `22`.
- `sample_3001`: motif CA RMSD `0.0762 A`, max pair-distance delta `0.1946 A`, mean pair-distance delta `0.0479 A`, ligand atom records `22`.

Evidence files on GPU host:

- Gate summary JSON: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_layered_l1_compact_publicckpt_20260702_motif_gate_summary.json`.
- Gate TSV: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_layered_l1_compact_publicckpt_20260702_motif_gate.tsv`.
- Queue status: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_layered_l1_queue_status.json`.

Two-parent sequence panel:

- Sequence panel run ID: `baker_layered_multiscaffold_ligandmpnn_from_compact_pass2_20260702_185358`.
- Parent scaffold universe for this panel: two motif-gated compact L1 PASS parents, `sample_3000` and `sample_3001`.
- Route: `generate_baker_layered_multiscaffold_ligandmpnn_bins.py`, parent-limit `2`, rounds `2`, bins `90/80/70/60/50`.
- Fixed motif positions: `14` scaffold-numbered motif residues from contig `6,A56-60,18,A83-85,8,A113-115,32,B145-147,6`.
- Parent scaffold length: `84` residues.
- Status: `DONE`, started `2026-07-02T18:53:59`, finished `2026-07-02T18:58:59`.

Selected sequence counts by target identity bin:

- `90`: `59` selected sequences.
- `80`: `441` selected sequences.
- `70`: `750` selected sequences.
- `60`: `799` selected sequences.
- `50`: `800` selected sequences.
- Total selected: `2849` sequences.

Selected sequence counts by parent and bin:

- `sample_3000`: 90=`23`, 80=`212`, 70=`350`, 60=`399`, 50=`400`.
- `sample_3001`: 90=`36`, 80=`229`, 70=`400`, 60=`400`, 50=`400`.

Audit interpretation:

- The scaffold gate PASS/FAIL denominator is the current present compact L1 PDB universe: `2` structures.
- The sequence panel denominator is the two PASS parent scaffolds used by the LigandMPNN panel run.
- These `2849` records are sequence-layer candidates passing the script's fixed motif/no-native and target mutation bin filters. They are not final structure-gated candidates.
- The 90% bin remains below the nominal high-diversity target because the current two compact parents yielded `59` filtered 90% records across two rounds. More PASS parents or additional rounds are needed for a larger 90% pool.

Remote sequence evidence files:

- Summary JSON: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_layered_multiscaffold_ligandmpnn_from_compact_pass2_20260702_185358_summary.json`.
- Selected sequence TSV: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_layered_multiscaffold_ligandmpnn_from_compact_pass2_20260702_185358_selected.tsv`.
- Sequence design directory: `/data/bht/project01_baker_serhyd_routeB_20260701/sequence_design/baker_layered_multiscaffold_ligandmpnn_from_compact_pass2_20260702_185358`.

Next action: continue monitoring compact L1 while `sample_3002` and later outputs are generated; each new PDB must be motif-gated before being added as a parent scaffold for additional sequence panel rounds.
