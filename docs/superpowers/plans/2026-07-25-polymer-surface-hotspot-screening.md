# Polymer Surface Hotspot Screening Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan.

**Goal:** Install and validate a reproducible AutoSite/AutoGrid4 surface-affinity workflow on Sugon, freeze PET and nylon probe/structure manifests, benchmark exact experimental controls, and prepare an audited high-throughput launch.

**Architecture:** ADFRsuite prepares receptors and produces AutoSite plus whole-receptor C/OA/HD maps. A Python 3 analysis environment extracts an exterior-connected surface shell, constructs connected equal-area patches, maps independently annotated catalytic residues to the shell, and reports catalytic enrichment versus the strongest off-target patch. Manifests, logs, compact tables, and independent audit artifacts remain on the server; source and plans live in GitHub.

**Tech Stack:** Bash, Slurm, ADFRsuite 1.0/1.1rc1 as verified at installation, AutoGrid4.2, AutoSite 1.0/1.1, Python 3.11, NumPy, SciPy, pandas, Biopython, RDKit, FreeSASA, pytest, PyYAML.

## Global Constraints

- Remote CPU root: /Dell/Dell14/bianht/enzyme_scaffold_search_v2.
- Phase 1 remote root: /Dell/Dell14/bianht/enzyme_scaffold_search_v2/polymer_surface_hotspot_screen_20260725.
- Maximum CPU allocation: 64.
- Do not run ESMFold on the CPU host.
- Do not modify or recompute existing PASS artifacts.
- Do not mix PET and nylon canonical denominators.
- Preserve experimental/PDB, exact-sequence predicted, and other structure sources as separate strata.
- Use NOT_EVALUATED_* for missing, unparsable, or unmapped records.
- Do not report AutoGrid or composite scores as binding free energy or enzymatic activity.
- All downloads, environments, caches, intermediate maps, and outputs stay on Sugon.
- Every actual server run appends a run-history record.
- Use GitHub commits for source changes; copy a commit-pinned snapshot to the server for execution.

## Task 1: Create the project skeleton

**Files:**

- Create: projects/03-polymer-surface-hotspot-screening/README.md
- Create: projects/03-polymer-surface-hotspot-screening/RUNBOOK.md
- Create: projects/03-polymer-surface-hotspot-screening/config/defaults.yaml
- Create: projects/03-polymer-surface-hotspot-screening/config/probes.tsv
- Create: projects/03-polymer-surface-hotspot-screening/scripts/install_phase1.sh
- Create: projects/03-polymer-surface-hotspot-screening/scripts/inventory_remote.sh
- Create: projects/03-polymer-surface-hotspot-screening/scripts/build_structure_manifest.py
- Create: projects/03-polymer-surface-hotspot-screening/scripts/prepare_receptor.sh
- Create: projects/03-polymer-surface-hotspot-screening/scripts/run_autosite_maps.sh
- Create: projects/03-polymer-surface-hotspot-screening/scripts/extract_surface_shell.py
- Create: projects/03-polymer-surface-hotspot-screening/scripts/score_surface_patches.py
- Create: projects/03-polymer-surface-hotspot-screening/scripts/run_smoke.sh
- Create: projects/03-polymer-surface-hotspot-screening/scripts/submit_benchmark.sh
- Create: projects/03-polymer-surface-hotspot-screening/scripts/audit_results.py
- Create: projects/03-polymer-surface-hotspot-screening/tests/test_surface_shell.py
- Create: projects/03-polymer-surface-hotspot-screening/tests/test_patch_scoring.py
- Create: projects/03-polymer-surface-hotspot-screening/tests/test_manifest.py

**Step 1: Verify repository policy**

Read:

    README.md
    AGENTS.md
    docs/superpowers/specs/2026-07-25-polymer-surface-hotspot-screening-design.md

Expected: remote-compute, local-read-only, scientific-gate, and provenance constraints agree.

**Step 2: Add the skeleton files**

Create files with executable server paths and no local-data assumptions.

**Step 3: Review tree**

Verify all expected paths exist in the branch and no installer, raw structure, map, cache, or secret is committed.

**Step 4: Commit**

Commit message:

    feat: scaffold polymer surface hotspot screening

## Task 2: Implement an idempotent remote inventory

**Files:**

- Modify: projects/03-polymer-surface-hotspot-screening/scripts/inventory_remote.sh
- Modify: projects/03-polymer-surface-hotspot-screening/RUNBOOK.md

**Step 1: Define the failing test**

The inventory script must fail if the resolved project root is outside:

    /Dell/Dell14/bianht/enzyme_scaffold_search_v2

It must not write until the root check passes.

**Step 2: Implement read-only inventory**

Record:

- hostname, kernel, CPU, memory, filesystem free space;
- module system availability;
- conda, mamba, micromamba, Python, gcc, make;
- agfr, autosite, autogrid4, prepare_receptor;
- freesasa, obabel;
- Slurm commands and account/partition visibility;
- existing phase-root status;
- existing exact structure manifests and final handoff directories.

Output to stdout before any installation.

**Step 3: Verify**

Run:

    bash scripts/inventory_remote.sh

Expected: exit 0, exact remote root printed, no files changed.

**Step 4: Commit**

Commit message:

    feat: add read-only Sugon inventory

## Task 3: Implement pinned server-local installation

**Files:**

- Modify: projects/03-polymer-surface-hotspot-screening/scripts/install_phase1.sh
- Modify: projects/03-polymer-surface-hotspot-screening/RUNBOOK.md
- Create: projects/03-polymer-surface-hotspot-screening/config/environment-lock.yml

**Step 1: Write installation assertions**

The installer must reject:

- empty project root;
- root outside the authorized CPU project;
- destination symlinks escaping the phase root;
- use of sudo;
- more than 64 build threads.

**Step 2: Install ADFRsuite without system modification**

Download the official Linux x86_64 ADFRsuite package from the CCSB/Scripps source into:

    software/downloads/

Record:

- source URL;
- retrieval timestamp;
- byte size;
- SHA256.

Install into:

    software/ADFRsuite/

Do not rely on system Python. Export the suite bin directory only inside wrapper scripts.

Because the suite contains an insulated Python 2.7 runtime, test it independently from the Python 3 analysis environment.

**Step 3: Install Python 3 analysis environment**

Create:

    envs/surface-screen-py311/

Pin at least:

- python 3.11;
- numpy;
- scipy;
- pandas;
- biopython;
- pyyaml;
- rdkit;
- freesasa-c;
- pytest.

Write an explicit environment lock after solving.

**Step 4: Run version checks**

Required commands:

    AutoSite --help
    autogrid4 -h
    prepare_receptor -h
    freesasa --version
    python -c "import numpy, scipy, pandas, Bio, yaml, rdkit"

Record complete stdout/stderr and exit status.

**Step 5: Create installation gate**

Write:

    results/gates/INSTALL_PASS.json

only if every required check passes. Otherwise write INSTALL_FAIL.json with failure categories. Do not claim a pass from partial installation.

**Step 6: Commit**

Commit message:

    feat: add pinned server-local installer

## Task 4: Freeze the explicit-probe manifest

**Files:**

- Modify: projects/03-polymer-surface-hotspot-screening/config/probes.tsv
- Create: projects/03-polymer-surface-hotspot-screening/scripts/build_probe_library.py
- Create: projects/03-polymer-surface-hotspot-screening/tests/test_probe_library.py

**Step 1: Define manifest columns**

Required columns:

    probe_id
    material_family
    role
    preferred_name
    canonical_smiles
    formal_charge
    heavy_atom_count
    total_atom_count
    primary_or_control
    stage
    rationale
    source
    status

**Step 2: Add PET probes**

Include:

- DMT as the primary neutral repeat-chemistry proxy;
- BHET as a secondary larger proxy;
- toluene as aromatic diagnostic;
- ethyl acetate as ester diagnostic;
- terephthalic acid and ethylene glycol as raw-monomer controls only.

**Step 3: Add nylon probes**

Include:

- N-methylacetamide as amide diagnostic;
- n-hexane as aliphatic diagnostic;
- N,N'-dimethyladipamide as neutral PA66 adipic-side proxy;
- N,N'-hexamethylenebisacetamide as neutral PA66 diamine-side proxy;
- caprolactam, adipic acid, and hexamethylenediamine as raw-monomer controls only;
- capped PA6 and PA66 oligomers as deferred Stage 3 entries.

**Step 4: Test chemistry**

Use RDKit to:

- parse every non-deferred SMILES;
- calculate formula, formal charge, heavy-atom count, total-atom count, rotatable bonds;
- generate canonical SMILES;
- reject mismatches between declared and computed properties;
- generate a checksumed SDF library on the server.

**Step 5: Commit**

Commit message:

    data: freeze PET and nylon probe definitions

## Task 5: Freeze the structure manifest

**Files:**

- Modify: projects/03-polymer-surface-hotspot-screening/scripts/build_structure_manifest.py
- Create: projects/03-polymer-surface-hotspot-screening/config/structure_manifest_schema.tsv
- Modify: projects/03-polymer-surface-hotspot-screening/tests/test_manifest.py
- Modify: projects/03-polymer-surface-hotspot-screening/RUNBOOK.md

**Step 1: Read authoritative inputs**

Read but do not modify:

    /Dell/Dell14/bianht/enzyme_scaffold_search_v2/results/FINAL_LAYER_SCORE_DATA_HANDOFF_20260715

Locate the existing PET and nylon assayed panels, structure manifests, activity labels, structure links, and catalytic annotations. Treat the legacy EXACT_OR_PUBLIC_STRUCTURE field as mixed-source until each row is reclassified from PDB/wwPDB metadata or the named prediction provider.

**Step 2: Define manifest columns**

Required columns:

    record_id
    material_family
    candidate_universe
    sequence_md5
    activity_class
    activity_evidence_source
    structure_path
    structure_sha256
    structure_source
    pdb_id
    chain_id
    receptor_context
    biological_assembly_id
    catalytic_residues
    oxyanion_residues
    catalytic_annotation_source
    residue_mapping_status
    include_phase
    exclusion_reason

**Step 3: Enforce provenance**

Tests must reject:

- duplicate record_id;
- structure path without checksum;
- experimental activity without a source field;
- missing structure treated as inactive;
- PET and nylon universe labels merged;
- predicted structure labeled as experimental;
- catalytic residue mappings outside the selected chain;
- NylC forced into the cross-family D-D-T annotation.

**Step 4: Write frozen manifests**

Generate:

    manifests/phase1a_pdb_experimental.tsv
    manifests/phase1b_public_exact_sequence_predicted.tsv
    manifests/phase1c_project_esmfold.tsv
    manifests/phase1d_benchmark_shard.tsv
    manifests/manifest_summary.json

The benchmark shard selection uses a fixed seed and stores the selection rationale.

**Step 5: Commit**

Commit message:

    data: add auditable structure manifest builder

## Task 6: Prepare receptors deterministically

**Files:**

- Modify: projects/03-polymer-surface-hotspot-screening/scripts/prepare_receptor.sh
- Create: projects/03-polymer-surface-hotspot-screening/scripts/check_receptor.py
- Create: projects/03-polymer-surface-hotspot-screening/tests/test_receptor_checks.py

**Step 1: Write failing fixture tests**

Fixtures cover:

- alternate locations;
- multiple models;
- missing catalytic residue;
- a retained catalytic metal;
- an unwanted crystallization ligand;
- a chain mapping mismatch.

**Step 2: Implement standardized receptor preparation**

Produce one prepared PDBQT per manifest row plus metadata JSON containing:

- input/output checksums;
- selected model and chain;
- removed hetero groups;
- retained cofactors;
- atom counts;
- bounding box;
- catalytic mapping;
- warnings;
- status.

**Step 3: Verify determinism**

Prepare each smoke receptor twice and compare PDBQT checksums and metadata, excluding timestamps.

**Step 4: Commit**

Commit message:

    feat: add deterministic receptor preparation

## Task 7: Generate AutoSite and whole-receptor maps

**Files:**

- Modify: projects/03-polymer-surface-hotspot-screening/scripts/run_autosite_maps.sh
- Create: projects/03-polymer-surface-hotspot-screening/scripts/build_grid_plan.py
- Create: projects/03-polymer-surface-hotspot-screening/tests/test_grid_plan.py

**Step 1: Test grid planning**

Tests cover:

- small receptor requiring one grid;
- large receptor requiring tiles;
- exact overlap;
- nonidentical C/OA/HD geometry rejection;
- tile escaping receptor margin.

**Step 2: Implement grid plan**

Default:

    spacing = 0.75 angstrom
    margin = 5.0 angstrom
    overlap = 8.0 angstrom

The script probes the installed AutoGrid4 build for dimension limits and stores the measured limit in installation metadata.

**Step 3: Generate C, OA, and HD maps**

Use identical geometry. Run AutoSite. Compress logs. Keep raw maps for smoke and benchmark until audit passes.

**Step 4: Validate**

Reject a receptor if:

- required map missing;
- dimensions differ;
- values are nonfinite;
- receptor does not fit inside the union of tiles;
- AutoSite exits nonzero.

A missing catalytic-overlap cluster is a descriptive flag, not a map-generation failure.

**Step 5: Commit**

Commit message:

    feat: generate tiled AutoSite and AutoGrid maps

## Task 8: Extract the exterior surface shell

**Files:**

- Modify: projects/03-polymer-surface-hotspot-screening/scripts/extract_surface_shell.py
- Modify: projects/03-polymer-surface-hotspot-screening/tests/test_surface_shell.py

**Step 1: Add synthetic tests**

Use synthetic sphere, two-lobed object, and object with an enclosed cavity.

Assert:

- boundary flood-fill selects exterior solvent;
- enclosed cavity is excluded from the primary shell;
- shell thickness is correct within one grid cell;
- rigid translation does not alter summary values;
- tiling and untiled extraction agree in overlap.

**Step 2: Implement occupancy and flood-fill**

Use atomic van der Waals radii, map geometry, boundary-connected solvent flood-fill, and a frozen near-surface distance band.

Store compact shell arrays only:

    coordinates
    C/OA/HD values
    graph-neighbor index
    nearest receptor atom/residue
    tile provenance

**Step 3: FreeSASA cross-check**

Compare total exposed receptor area and exposed-residue set with FreeSASA. Large disagreement triggers SHELL_SASA_DISAGREEMENT for review.

**Step 4: Commit**

Commit message:

    feat: extract exterior-connected affinity shell

## Task 9: Score catalytic and off-target patches

**Files:**

- Modify: projects/03-polymer-surface-hotspot-screening/scripts/score_surface_patches.py
- Modify: projects/03-polymer-surface-hotspot-screening/tests/test_patch_scoring.py

**Step 1: Add tests**

Cover:

- catalytic hotspot strongest;
- off-target hotspot strongest;
- globally sticky surface;
- equal means but different top-20-percent means;
- disconnected surface;
- missing catalytic mapping;
- deterministic farthest-point seeds.

**Step 2: Implement nested catalytic patches**

Generate catalytic_core, catalytic_neighborhood, and catalytic_extended from surface-graph distance to independently annotated catalytic anchors.

**Step 3: Implement equal-area off-target patches**

Use deterministic farthest-point sampling and connected graph growth. Exclude the catalytic region plus buffer.

**Step 4: Preserve primary channel metrics**

Report C, OA, and HD separately. Add PET and nylon secondary composites only after calibration statistics are frozen.

**Step 5: Add interpretation flags**

Generate POCKET_ENRICHED, GLOBAL_STICKY, OFF_POCKET_HOTSPOT, CATALYTIC_REGION_NO_AUTOSITE_CLUSTER, and technical NOT_EVALUATED flags.

**Step 6: Commit**

Commit message:

    feat: compare catalytic and off-target surface patches

## Task 10: Run PDB-experimental smoke tests

**Files:**

- Modify: projects/03-polymer-surface-hotspot-screening/scripts/run_smoke.sh
- Modify: projects/03-polymer-surface-hotspot-screening/RUNBOOK.md

**Step 1: Freeze smoke manifest**

Minimum:

- PET 6ILW;
- one additional PET PDB-experimental family control;
- NylC 3AXG standardized-chain receptor;
- NylC 3AXG biological assembly context;
- NylB 1WYB;
- PA66 hydrolase Nyl50 9DYS;
- at least one experimentally inactive or screen-negative PDB-experimental structure if available.

**Step 2: Run in foreground or a small Slurm job**

Use no more than 8 CPUs for smoke.

**Step 3: Rerun determinism test**

Run the same manifest twice into separate directories and compare compact result hashes and metrics.

**Step 4: Audit smoke**

Required outputs:

    results/smoke/summary.tsv
    results/smoke/failures.tsv
    results/smoke/runtime.tsv
    results/gates/SMOKE_PASS.json or SMOKE_FAIL.json

**Step 5: Stop condition**

Do not submit benchmark if smoke gate fails.

**Step 6: Commit**

Commit message:

    test: add exact-structure smoke workflow

## Task 11: Run the stratified benchmark

**Files:**

- Modify: projects/03-polymer-surface-hotspot-screening/scripts/submit_benchmark.sh
- Create: projects/03-polymer-surface-hotspot-screening/slurm/benchmark_array.sbatch
- Modify: projects/03-polymer-surface-hotspot-screening/RUNBOOK.md

**Step 1: Freeze benchmark manifest**

Select deterministic PET and nylon records separately across:

- family;
- activity class;
- structure source;
- length bin;
- existing project score bin.

Store selected and unselected denominator summaries.

**Step 2: Submit bounded Slurm arrays**

- CPU total never exceeds 64;
- one receptor per array task or measured safe chunk;
- task-local scratch only under the phase root;
- resume skips only checksum-validated completed records;
- each failure writes a categorized row and does not block later records.

**Step 3: Monitor read-only**

Inspect scheduler state and artifacts without canceling or resubmitting unless explicitly authorized.

**Step 4: Merge by frozen precedence**

A record is rankable only when all required stages pass. Missing data remain NOT_EVALUATED.

**Step 5: Commit**

Commit message:

    feat: add bounded benchmark array workflow

## Task 12: Independent audit and production decision

**Files:**

- Modify: projects/03-polymer-surface-hotspot-screening/scripts/audit_results.py
- Create: projects/03-polymer-surface-hotspot-screening/tests/test_audit.py
- Modify: projects/03-polymer-surface-hotspot-screening/RUNBOOK.md

**Step 1: Recompute denominators**

Audit independently recomputes:

- manifest row counts;
- unique sequence MD5 counts;
- unique structure checksums;
- result row counts;
- evaluated/failure/NOT_EVALUATED counts;
- PET/nylon separation;
- exact/predicted separation;
- best off-target patch selection;
- catalytic percentiles;
- summary-table provenance.

**Step 2: Validate scientific gates**

Check stability to:

- repeated run;
- grid translation;
- selected spacing changes;
- tile boundary;
- monomer versus assembly context;
- structure source.

**Step 3: Emit gate artifact**

Write exactly one:

    results/gates/BENCHMARK_PASS.json
    results/gates/BENCHMARK_FAIL.json

A scheduler exit 0 or complete shard set is insufficient.

**Step 4: Production recommendation**

If PASS, freeze a production manifest and estimate CPU-hours/storage before submission. If FAIL, list which method assumption failed and the next bounded experiment.

**Step 5: Commit**

Commit message:

    audit: add independent benchmark gate

## Task 13: Open and maintain a Draft PR

**Step 1: Open Draft PR**

Title:

    research: add PET and nylon surface hotspot screening phase 1

The body must state:

- what is implemented;
- what remains remote and unexecuted;
- current SSH/connectivity status;
- scientific non-claims;
- exact files and gates expected from the server run.

**Step 2: Review**

Check that no local or server-only artifacts are committed and all planned commands are reproducible.

**Step 3: Keep Draft status**

Do not mark ready or merge until PDB-experimental smoke and benchmark audit artifacts are reviewed.
