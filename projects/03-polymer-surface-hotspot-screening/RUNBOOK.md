# RUNBOOK: PET and nylon surface hotspot screening

## Purpose

Run the Phase 1 AutoSite/AutoGrid4 coarse-grained screen on Sugon while preserving exact input denominators, structure provenance, technical failures, and scientific interpretation boundaries.

This runbook does not authorize recomputation of existing PASS artifacts.

## Fixed locations

Repository:

    https://github.com/bianhuiting918/git-hub-chatgpt-pro

Branch during development:

    codex/polymer-surface-hotspot-screening-phase1

Remote CPU project:

    /work/home/acshdt1dks

Remote Phase 1 root:

    /work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725

Authoritative read-only structure manifests:

    PET (8,343 structure-linked rows):
    /work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/outputs/saprot_structure_conditional_likelihood_20260723/run_manifest.tsv

    Nylon primary authority (3,497 routes; use only files readable on Sugon):
    /work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/manifests/saprot_nylonase_4556_20260725_confscale_v2/run_manifest_all_3497.tsv

    Nylon recovered and preflight-ready (735 rows):
    /work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/manifests/saprot_nylonase_recovered_20260725_v7/run_manifest_ready_735.tsv

CPU limit:

    64

## Connection

Use the configured SSH route:

    ssh eshell111.hpccube.com

Do not bypass the SSH configuration. If the configured BindAddress is unavailable, stop and restore the approved network route.

The Sugon server does not need outbound GitHub access. Deploy exact files from a GitHub commit through the connected workstation into:

    /work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725/deployments/COMMIT_SHA/

For every deployed file, compare the server-side `git hash-object` value with the GitHub blob SHA before execution. Never modify the server repository checkout to work around an outbound network failure.

## Required remote directory tree

After confirming the resolved root, create:

    software/downloads
    software/ADFRsuite
    envs
    cache
    inputs
    manifests
    work
    results/gates
    logs
    scripts
    audits

No directory is created before the root check succeeds.

## Phase 0: read-only inventory

From a commit-pinned source snapshot on the server:

    cd /work/home/acshdt1dks
    bash polymer_surface_hotspot_screen_20260725/scripts/inventory_remote.sh

Record:

- hostname and filesystem;
- scheduler visibility;
- existing tools;
- exact manifest locations;
- free space;
- whether the phase root already exists.

Expected: no mutation.

## Phase 1: software installation

The installer is unprivileged and idempotent.

    cd /work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725
    bash scripts/install_phase1.sh

Required components:

- ADFRsuite containing AutoSite, AutoGrid4, prepare_receptor, and AGFR;
- isolated Python 3.11 analysis environment;
- NumPy, SciPy, pandas, Biopython, PyYAML, RDKit, FreeSASA, pytest.

The installer must write download URL, timestamp, byte count, and SHA256 for every installer/archive.

Pass artifact:

    results/gates/INSTALL_PASS.json

If absent, do not proceed.

## Phase 2: probe manifest

The primary map scan has no explicit ligand file.

The explicit probe manifest is for validation stages and chemistry controls. From a blob-verified commit snapshot, build and validate it with the pinned RDKit environment:

    /work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725/envs/surface-screen-py311/bin/python \
      deployments/COMMIT_SHA/projects/03-polymer-surface-hotspot-screening/scripts/build_probe_library.py \
      --manifest deployments/COMMIT_SHA/projects/03-polymer-surface-hotspot-screening/config/probes.tsv \
      --output-dir results/probes_v1

The builder retains an existing `PROBES_PASS.json` without overwrite and refuses partial pre-existing outputs. Stage 3 capped oligomers remain deferred.

Required checks:

- every active SMILES parses;
- formal charge matches;
- heavy-atom and total-atom counts match;
- canonical SMILES and SDF checksums are recorded;
- deferred oligomers are not silently generated from guesses.

## Phase 3: structure manifest

Generate manifests only from frozen authoritative inputs:

    python scripts/build_structure_manifest.py \
      --pet-manifest /work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/outputs/saprot_structure_conditional_likelihood_20260723/run_manifest.tsv \
      --nylon-manifest /work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/manifests/saprot_nylonase_4556_20260725_confscale_v2/run_manifest_all_3497.tsv \
      --nylon-recovered-manifest /work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/manifests/saprot_nylonase_recovered_20260725_v7/run_manifest_ready_735.tsv \
      --output-dir manifests \
      --seed 20260725

Expected outputs:

    manifests/phase1a_pdb_experimental.tsv
    manifests/phase1b_public_exact_sequence_predicted.tsv
    manifests/phase1c_project_esmfold.tsv
    manifests/phase1d_benchmark_shard.tsv
    manifests/manifest_summary.json
    manifests/SHA256SUMS

Before use, inspect:

- PET versus nylon counts;
- PDB experimental versus public predicted versus project ESMFold counts;
- legacy EXACT_OR_PUBLIC_STRUCTURE rows reclassified by row-level provenance;
- activity evidence labels;
- sequence MD5 counts;
- structure checksums;
- catalytic mapping source;
- biological assembly context;
- exclusions and NOT_EVALUATED categories.

## Phase 3.5: exterior-connected surface shell

Only run after a checksum-verified `MAPS_PASS.json`. Use the commit-pinned shell extractor and keep PET and nylon result roots separate:

    python deployments/COMMIT_SHA/projects/03-polymer-surface-hotspot-screening/scripts/extract_surface_shell.py \
      --maps-root results/grid_maps_smoke_RUN_ID/FAMILY_RECORD \
      --receptor-pdbqt results/receptor_determinism_RUN_ID/FAMILY_RECORD/run1/receptor.pdbqt \
      --receptor-pdb results/receptor_determinism_RUN_ID/FAMILY_RECORD/run1/selected.pdb \
      --output-dir results/surface_shell_smoke_RUN_ID/FAMILY_RECORD \
      --record-id FAMILY_RECORD

Frozen shell definition:

- boundary-connected solvent only;
- 26-connectivity;
- receptor occupancy from atomic van der Waals radii plus a 1.40 A solvent probe;
- primary near-surface band 0.00-2.25 A;
- sensitivity bands 0.00-1.50 A and 0.00-3.00 A;
- compact coordinates, A/C/OA/HD values, CSR neighbor graph, nearest receptor residue, and tile provenance.

The primary shell gate is `SHELL_PASS.json`. `SHELL_SASA_DISAGREEMENT` is a technical/scientific review flag, not evidence of inactivity. Verify `SHA256SUMS` before patch scoring. Never overwrite a partial or PASS shell directory.

## Phase 4: smoke test

Minimum smoke controls:

- PET 6ILW;
- at least one additional PET PDB-experimental structure;
- NylC 3AXG standardized-chain context;
- NylC 3AXG biological-assembly context;
- NylB 1WYB;
- PA66 hydrolase Nyl50 9DYS;
- at least one PDB-experimental inactive or screen-negative control if present in the frozen manifest.

Run with at most 8 CPUs:

    bash scripts/run_smoke.sh \
      --manifest manifests/phase1a_smoke.tsv \
      --cpus 8 \
      --run-id smoke_20260725_a

Repeat under a second run ID for determinism.

Required outputs:

    results/smoke_20260725_a/summary.tsv
    results/smoke_20260725_a/patches.tsv
    results/smoke_20260725_a/failures.tsv
    results/smoke_20260725_a/runtime.tsv
    results/gates/SMOKE_PASS.json

Do not interpret results if SMOKE_PASS.json is absent.

## Phase 5: benchmark shard

Only after smoke PASS:

    bash scripts/submit_benchmark.sh \
      --manifest manifests/phase1d_benchmark_shard.tsv \
      --max-total-cpus 64 \
      --run-id benchmark_20260725_a

Monitor read-only:

    squeue -u "$USER"
    sacct -j JOB_ID --format=JobID,State,ExitCode,Elapsed,MaxRSS,AllocCPUS

Do not cancel, resubmit, or overwrite completed records without explicit authorization.

## Phase 6: merge and audit

After all array tasks finish or have terminal failure records:

    python scripts/audit_results.py \
      --manifest manifests/phase1d_benchmark_shard.tsv \
      --results results/benchmark_20260725_a \
      --audit-dir audits/benchmark_20260725_a

Required gate:

    results/gates/BENCHMARK_PASS.json

COMPLETE, scheduler exit 0, or presence of all shard files is not sufficient.

## Output interpretation

Primary outputs are channel-specific:

- A aromatic-carbon catalytic enrichment;
- C aliphatic-carbon catalytic enrichment;
- OA acceptor-oxygen catalytic enrichment;
- HD donor-hydrogen catalytic enrichment;
- maximum equal-area off-target patch;
- within-protein catalytic percentile;
- global sticky-surface fraction.

Secondary composites:

- PET = mean(z_A, z_C, z_OA);
- nylon = mean(z_C, z_OA, z_HD).

These are screening proxies, not binding energies.

## Failure handling

Examples:

- NOT_EVALUATED_MISSING_STRUCTURE
- NOT_EVALUATED_STRUCTURE_PARSE
- NOT_EVALUATED_CATALYTIC_MAPPING
- NOT_EVALUATED_RECEPTOR_PREP
- NOT_EVALUATED_GRID
- NOT_EVALUATED_SURFACE_SHELL
- NOT_EVALUATED_PATCH_GENERATION
- CATALYTIC_REGION_NO_AUTOSITE_CLUSTER
- SHELL_SASA_DISAGREEMENT

Later records continue after a single-record failure. Retry policy is decided only after failure counts are reviewed.

## Run history

Append one TSV row per actual run to:

    logs/run_history.tsv

Columns:

    timestamp_local
    timestamp_utc
    git_commit
    run_id
    script
    parameters
    input_manifest
    input_sha256
    output_path
    scheduler_job_id
    exit_status
    n_input
    n_evaluated
    n_failed
    n_not_evaluated
    summary

Never write secrets into this file.

## Recovery

- Re-run inventory before resuming after a connection outage.
- Use checksum-validated record completion markers.
- Never infer completion from a nonempty directory.
- Never overwrite existing run IDs.
- Create a new run ID for parameter changes.
- Preserve failed logs until the failure ledger and audit are complete.
