# Polymer-surface hotspot screening

This project develops an auditable, ligand-light first-pass screen for asking whether PET- or nylon-related chemical affinity is enriched near an enzyme's catalytic region rather than at equally sized off-target surface patches.

The first-stage method is:

1. prepare a standardized receptor structure;
2. compute whole-receptor AutoGrid4 C, OA, and HD affinity maps;
3. use AutoSite only to identify and characterize high-affinity clusters;
4. extract an exterior-connected grid shell around the entire receptor;
5. score PET-like and nylon-like chemical compatibility on that shell;
6. define a catalytic surface region from independently annotated catalytic residues, not from the highest-scoring material hotspot;
7. compare the catalytic region with equal-area connected off-target patches;
8. preserve channel-level scores, uncertainty labels, structure provenance, and failure reasons.

## Scientific question

Primary question:

> Is material-related chemical compatibility enriched around the independently annotated catalytic region compared with the strongest equally sized non-catalytic surface patch?

This is a ranking proxy. It is not a prediction of polymer adsorption free energy, productive binding geometry, catalytic activity, kcat, KM, kcat/KM, or degradation rate.

## Phase 1 scope

Phase 1 is limited to PET hydrolases and nylon-related hydrolases with traceable experimental activity labels and available experimental or explicitly identified predicted structures.

It starts with a small benchmark, not the full candidate universe:

- PET exact experimental structures and active/inactive controls from the frozen project manifest;
- nylon PDB experimental structures, initially NylC 3AXG, NylB 1WYB, and PA66 hydrolase Nyl50 9DYS, plus additional controls only after row-level provenance verification;
- a structure-source-stratified extension using public exact-sequence predicted structures and project ESMFold structures only after the PDB-experimental benchmark passes;
- candidate-scale production only after runtime, failure-rate, and positive-control gates are met.

The existing PET and nylon candidate universes remain separate. Candidate IDs, exact experimental sequences, structure-linked records, and final plotted points must never be reported as the same denominator. The legacy handoff label EXACT_OR_PUBLIC_STRUCTURE is also not treated as synonymous with PDB experiment because it includes public predicted models.

## Material representations

### First-stage map templates

No ligand file is needed for the primary surface scan.

AutoGrid4 map energies are converted to favorable affinity values by multiplying by minus one and then kept as separate channels.

- PET map template: C and OA channels.
- Nylon map template: C, OA, and HD channels.

Primary selection is Pareto-based across channel-specific catalytic enrichment metrics. A fixed equal-weight composite is reported only as a secondary convenience score:

- PET secondary composite = mean(z_C, z_OA)
- Nylon secondary composite = mean(z_C, z_OA, z_HD)

All z values are robustly standardized using a frozen calibration set. Within-protein percentile enrichment is also reported so that localization and cross-protein magnitude remain distinguishable.

### Explicit validation probes

These are not used to define the primary Phase 1 result. They are reserved for later AutoLigand or docking validation.

PET:

- dimethyl terephthalate (DMT): primary neutral repeat-chemistry proxy;
- bis(2-hydroxyethyl) terephthalate (BHET): secondary larger proxy;
- toluene or benzene: aromatic diagnostic;
- ethyl acetate: ester diagnostic.

Nylon:

- N-methylacetamide: neutral amide diagnostic;
- n-hexane: aliphatic diagnostic;
- N,N'-dimethyladipamide: neutral adipic-side PA66 fragment;
- N,N'-hexamethylenebisacetamide: neutral diamine-side PA66 fragment;
- linear capped PA6 and PA66 oligomers: later-stage probes only.

Adipic acid, hexamethylenediamine, caprolactam, terephthalic acid, and ethylene glycol may be retained as raw-monomer chemistry controls, but they are not primary polymer-surface proxies because their charge state or connectivity differs from the polymer repeat.

## Core outputs

Each receptor receives:

- structure source and exact-sequence provenance;
- receptor preparation status;
- biological assembly context label;
- catalytic residue annotation source;
- AutoSite cluster table;
- C, OA, and HD catalytic-region metrics;
- equal-area maximum off-target patch metrics;
- catalytic-minus-off-target differences;
- catalytic percentiles among all connected surface patches;
- off-target hotspot count and area;
- global sticky-surface fraction;
- runtime and peak-memory estimates;
- explicit NOT_EVALUATED reason where needed.

A result is scientifically rankable only if receptor preparation, map generation, catalytic annotation, surface-shell extraction, and equal-area comparison all pass.

## Repository layout

- docs/superpowers/specs: scientific design and interpretation contract.
- docs/superpowers/plans: executable implementation plan.
- projects/03-polymer-surface-hotspot-screening/RUNBOOK.md: server commands and recovery instructions.
- projects/03-polymer-surface-hotspot-screening/config: frozen configurations and probe metadata.
- projects/03-polymer-surface-hotspot-screening/scripts: server-side installation, manifest, smoke-test, and production scripts.
- projects/03-polymer-surface-hotspot-screening/tests: small synthetic and known-structure tests.

Large structures, installers, affinity maps, caches, and result shards stay on the Sugon server and are never committed.

## Remote project location

Server project root:

    /Dell/Dell14/bianht/enzyme_scaffold_search_v2

Phase 1 working directory:

    /Dell/Dell14/bianht/enzyme_scaffold_search_v2/polymer_surface_hotspot_screen_20260725

All environments, installers, caches, scripts copied for execution, manifests, logs, and results must remain under that remote directory. CPU usage is capped at 64 threads. No ESMFold inference is started on the CPU host.

## Status

- Research design: approved for implementation.
- GitHub branch: codex/polymer-surface-hotspot-screening-phase1.
- Remote installation: blocked until the configured SSH route is reachable.
- Production calculation: not started.
