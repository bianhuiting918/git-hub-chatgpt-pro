# Polymer Surface Hotspot Screening Design

Date: 2026-07-25  
Status: Approved for Phase 1 implementation  
Owner: enzyme-scaffold screening project

## 1. Objective

Build a reproducible, high-throughput, coarse-grained calculation that compares PET-like or nylon-like chemical affinity near an independently annotated catalytic region with affinity at equally sized regions elsewhere on the same protein surface.

The design intentionally does not assume a complete polymer binding pose. It asks a narrower question:

> Does the catalytic neighborhood exhibit material-related chemical-field enrichment relative to the strongest competing surface patch?

This is appropriate for initial prioritization when polymer binding modes are unknown and full oligomer docking or molecular dynamics is too expensive.

## 2. Claims and non-claims

### Supported outputs

- relative catalytic-region enrichment within a protein;
- comparison with equal-area off-target surface patches;
- channel-specific PET-like and nylon-like chemical compatibility;
- ranking within a frozen, consistently prepared cohort;
- identification of proteins that are globally sticky rather than catalytically localized;
- identification of off-pocket hotspots requiring later interpretation.

### Unsupported outputs

The calculation does not directly provide:

- polymer adsorption free energy;
- Kd, kon, or koff;
- a productive catalytic pose;
- solid-state crystallinity or chain-extraction cost;
- catalytic rate, substrate turnover, or degradation yield;
- proof that an off-pocket hotspot is harmful;
- equivalence between predicted structures and experimental structures.

Every result table and report must preserve these boundaries.

## 3. Cohort and denominator contract

The existing enzyme-scaffold project contains multiple non-equivalent universes. They must remain separate.

### PET-related universes

- canonical candidate IDs;
- exact experimental sequences;
- PDB experimental structures;
- public exact-sequence predicted structures;
- project exact-sequence ESMFold structures;
- structure-linked candidates;
- final plotted or filtered points.

### Nylon-related universes

- canonical candidate IDs;
- exact experimental sequences;
- exact public structures;
- exact-sequence predicted structures;
- structure-linked candidates;
- final plotted or filtered points.

Every manifest and report must state:

- material family;
- universe name;
- number of input records;
- number of unique sequences by MD5;
- number of unique structure files by checksum;
- experimental activity class;
- structure source;
- biological assembly handling;
- evaluated, failed, and NOT_EVALUATED counts.

PET and nylon candidate universes must not be merged before family-specific metrics are complete.

## 4. Phased cohort expansion

### Phase 1A: PDB experimental-structure controls

Begin with PDB experimental structures and experimentally characterized positive and negative controls. The legacy handoff field EXACT_OR_PUBLIC_STRUCTURE is mixed-source and must be reclassified row by row before inclusion.

Minimum PET positive-control seeds already present in the authoritative handoff include:

- 6ILW;
- 4EB0;
- 4CG1;
- 4WFI;
- 3VIS.

The exact manifest, rather than this short seed list, is authoritative for the full PET control panel.

Minimum nylon PDB experimental controls:

- NylC 3AXG, analyzed both as a standardized chain-level receptor and in its documented oligomeric context;
- NylB 1WYB, an X-ray structure of 6-aminohexanoate-dimer hydrolase;
- Nyl50 9DYS, an X-ray structure of a PA66 hydrolase bound to tetraethylene glycol.

The frozen project manifests supply candidate and activity provenance, while RCSB/wwPDB metadata determine whether a coordinate set is genuinely experimental.

### Phase 1B: public exact-sequence predicted structures

Only after Phase 1A gates pass, add experimentally tested sequences having public exact-sequence AFDB or other computed models. Label provider, model version, and exact-sequence mapping explicitly.

### Phase 1C: project exact-sequence ESMFold structures

Add project ESMFold structures as a separate stratum after Phase 1B. Do not pool PDB, public predicted, and project ESMFold structures silently.

### Phase 1D: candidate-scale benchmark shard

Select a deterministic, stratified shard from each canonical candidate universe. Stratify by family, structure source, sequence-length bin, and existing project score bin. The shard is used to estimate runtime, map-size distribution, failure rate, and score stability.

### Phase 1E: production

Launch the full structure-linked universe only after installation, smoke, benchmark, and independent audit gates pass. Missing or unparsable structures remain NOT_EVALUATED and do not become biological failures.

## 5. Structure representation

### Standardized comparison receptor

The primary high-throughput ranking uses one standardized protein chain per record, selected by the frozen structure manifest. This matches the existing first-chain structure-linked workflow and avoids mixing different oligomer sizes in one ranking.

### Biological-assembly sensitivity receptor

Known experimental assemblies are evaluated separately where biological assembly is mechanistically important. NylC is a required assembly-aware control because autoproteolysis and oligomerization affect its exposed surfaces.

The two receptor contexts generate separate result namespaces:

- monomer_standardized;
- biological_assembly_context.

They are never pooled into a single rank.

## 6. Receptor preparation

For each structure:

1. verify file checksum and provenance;
2. select the intended model and chain or biological assembly;
3. resolve alternate locations deterministically by highest occupancy, then lexical altloc;
4. remove crystallization additives and non-mechanistic ligands for the apo surface scan;
5. retain catalytic or structural metals only when a frozen manifest explicitly requires them;
6. retain missing-residue annotations;
7. add hydrogens and AutoDock atom types with ADFRsuite prepare_receptor;
8. record warnings, atom count, net charge proxy, bounding box, and prepared PDBQT checksum.

A preparation error produces a specific NOT_EVALUATED_PREP_* reason.

## 7. Catalytic annotation

Catalytic regions are defined independently of predicted material affinity.

Each structure manifest row contains:

- catalytic nucleophile;
- catalytic acid/base residues;
- oxyanion-hole residues when known;
- annotation source;
- residue-number mapping status;
- confidence class.

PET seed residue annotations come from the authoritative project handoff and must be aligned to the actual processed chain before use.

For nylon cross-family candidate comparisons, the current versioned D-D-T definition may be retained as a motif anchor, but it is not automatically treated as a complete substrate-binding pocket. NylC uses its known N-terminal Thr catalytic mechanism and must not be forced into the D-D-T convention.

If catalytic residue mapping fails, the record is NOT_EVALUATED_CATALYTIC_MAPPING rather than inactive.

## 8. Whole-receptor affinity maps

AutoGrid4 maps are computed for:

- C: hydrophobic carbon probe;
- OA: ligand hydrogen-bond acceptor probe;
- HD: ligand hydrogen-bond donor probe.

AutoSite is run from the same receptor preparation to identify high-affinity clusters and feature types. AutoSite clusters are descriptive; they do not define the catalytic anchor and are not the final score.

### Grid geometry

- default spacing: 0.75 angstrom;
- receptor bounding box plus 5.0 angstrom margin;
- identical channel geometry for C, OA, and HD;
- deterministic tiling when a receptor exceeds the verified AutoGrid4 dimension limit;
- tile overlap: at least 8.0 angstrom;
- retain only exterior surface-shell values after merging;
- in overlap regions, retain the most favorable map energy per channel and record the overlap count.

The installed AutoGrid4 build is tested for supported grid dimensions before production. The configuration is updated from measured capabilities, not assumed.

## 9. Exterior-connected surface shell

The comparison uses affinity-map grid points near the exterior solvent-accessible receptor surface.

Algorithm:

1. rasterize receptor occupancy from atomic van der Waals radii on the map grid;
2. flood-fill solvent voxels from grid boundaries;
3. exclude enclosed internal cavities from the primary exterior surface analysis;
4. select exterior solvent voxels within a fixed distance band of receptor occupancy;
5. retain the same surface-shell definition across all proteins;
6. construct a 26-neighbor graph on retained voxels;
7. assign each voxel C, OA, and HD affinity values.

FreeSASA is used as an independent per-residue SASA and total-SASA quality check, not as the sole source of patch geometry.

Sensitivity analyses test at least two shell thicknesses and two grid spacings on the PDB-experimental controls.

## 10. Material channel definitions

AutoGrid energies are lower when more favorable. Define affinity values:

    A_C  = -E_C
    A_OA = -E_OA
    A_HD = -E_HD

Preserve raw energies and transformed affinities.

### PET-like compatibility

Primary channels:

- A_C;
- A_OA.

Secondary composite after robust calibration:

    PET_composite = mean(z_C, z_OA)

PET contains ester acceptors and aromatic/hydrophobic surfaces but no repeat-unit hydrogen-bond donor, so HD is not included in the primary PET composite.

### Nylon-like compatibility

Primary channels:

- A_C;
- A_OA;
- A_HD.

Secondary composite after robust calibration:

    NYLON_composite = mean(z_C, z_OA, z_HD)

The composite is an intentionally simple proxy. Primary selection remains multi-objective so one favorable channel cannot fully mask an unfavorable channel.

## 11. Normalization

Two distinct normalizations are reported.

### Within-protein localization

Each surface voxel and connected patch receives a percentile relative to patches from the same receptor. This answers whether the catalytic region is unusually favorable within that protein.

### Cross-protein magnitude

Raw map-derived patch summaries are standardized against a frozen PDB-experimental calibration cohort using median and median absolute deviation. The calibration statistics are versioned and are not recalculated after inspecting candidate results.

This separation prevents per-protein normalization from erasing absolute differences while preventing global magnitude from being confused with localization.

## 12. Catalytic surface region

For each annotated catalytic anchor:

1. find the nearest exterior surface-shell voxel;
2. grow a connected surface patch by graph distance to a frozen target radius;
3. optionally add an AutoSite cluster only if it overlaps the catalytic neighborhood under a predeclared distance criterion;
4. cap or expand the final region to a frozen target voxel count or approximate area for same-area comparisons.

Three nested regions are retained:

- catalytic_core;
- catalytic_neighborhood;
- catalytic_extended.

If AutoSite has no overlapping cluster, report CATALYTIC_REGION_NO_AUTOSITE_CLUSTER and continue with the anchor-defined surface region.

## 13. Equal-area off-target comparison

Generate connected off-target patches with the same voxel count and shell definition as each catalytic region.

Patch seeds are chosen by deterministic farthest-point sampling over the surface graph. Patches overlapping the catalytic exclusion zone are removed.

For every material/channel combination report:

- catalytic patch mean;
- catalytic top-20-percent mean;
- catalytic patch maximum;
- maximum equally sized off-target patch;
- median off-target patch;
- catalytic minus maximum off-target;
- catalytic percentile among all patches;
- number and total approximate area of off-target patches exceeding a frozen threshold;
- fraction of the entire surface exceeding the threshold.

The maximum off-target comparator prevents a large background surface from diluting a single strong competing hotspot.

## 14. Interpretation classes

Results receive descriptive, non-mechanistic flags.

- POCKET_ENRICHED: catalytic region exceeds the frozen enrichment gate.
- GLOBAL_STICKY: large fraction of the whole surface is favorable.
- OFF_POCKET_HOTSPOT: an equally sized off-target patch exceeds the catalytic region.
- CATALYTIC_REGION_NO_AUTOSITE_CLUSTER: anchor-defined region lacks an AutoSite cluster.
- STRUCTURE_CONTEXT_SENSITIVE: monomer and assembly-aware results differ beyond the frozen stability tolerance.
- NOT_EVALUATED_*: technical or provenance failure.

OFF_POCKET_HOTSPOT is not automatically penalized biologically. It may be nonproductive adsorption, a polymer-capture region, an oligomerization interface, or a structure artifact.

## 15. Explicit ligand validation

Only a small top-ranked and control set proceeds to explicit ligand placement.

### PET probe set

- DMT;
- BHET;
- aromatic diagnostic;
- ester diagnostic.

### Nylon probe set

- N-methylacetamide;
- n-hexane;
- N,N'-dimethyladipamide;
- N,N'-hexamethylenebisacetamide.

Raw synthesis monomers are chemistry controls. Large oligomers are deferred until the surface-map method demonstrates useful enrichment on known controls.

AutoLigand fixed-volume pseudo-ligands may be used before conventional docking. Explicit ligand scores are validation columns and do not overwrite primary map results.

## 16. Controls

### Positive controls

- known PET hydrolase catalytic grooves;
- NylC PDB experimental structure;
- NylB 1WYB and Nyl50 9DYS;
- additional PDB experimental nylon hydrolase or substrate-complex structures after row-level provenance verification.

### Negative controls

- experimentally inactive or screen-negative exact-sequence proteins;
- catalytic-site mutants when exact structures exist;
- same protein with randomized catalytic anchor positions as a localization null;
- channel-label permutations as a scoring null.

### Technical controls

- identical structure processed twice for determinism;
- rigid rotation and translation of a receptor;
- grid-spacing sensitivity;
- tile-boundary sensitivity;
- monomer versus biological assembly;
- structure-source stratification.

## 17. Acceptance gates

### Installation gate

Pass only if:

- versions and checksums are recorded;
- AutoSite and AutoGrid4 execute on a bundled test receptor;
- prepare_receptor generates deterministic atom counts and checksum;
- FreeSASA runs;
- the environment does not modify system Python;
- all artifacts reside under the remote project directory.

### Smoke gate

Pass only if:

- at least one PET and one nylon PDB experimental structure complete;
- C, OA, and HD maps share geometry;
- exterior shell is nonempty and exterior-connected;
- catalytic residues map successfully;
- equal-area off-target patches are generated;
- rerun metrics are stable within the frozen tolerance;
- no silent failure is converted to zero.

### Benchmark gate

Pass only if:

- input denominator and structure source are frozen;
- technical success rate meets the predeclared threshold;
- failures are categorized;
- positive controls show interpretable catalytic localization in at least one material-relevant channel;
- negative controls do not all collapse into the same score range;
- runtime and storage permit production within the 64-CPU constraint;
- an independent audit validates row counts, checksums, best-patch selection, and summary denominators.

Failure of biological separation does not mean the software failed; it means the scientific proxy is not yet validated.

## 18. Server and reproducibility policy

Remote root:

    /Dell/Dell14/bianht/enzyme_scaffold_search_v2/polymer_surface_hotspot_screen_20260725

Required subdirectories:

    software/
    envs/
    cache/
    inputs/
    manifests/
    work/
    results/
    logs/
    scripts/
    audits/

Each actual run appends one tab-separated record with:

- UTC and local timestamp;
- git commit;
- script;
- parameters;
- input manifest and checksum;
- output path;
- scheduler job ID;
- exit status;
- evaluated/pass/fail/NOT_EVALUATED counts;
- short result summary.

No credentials are stored. Raw maps may be deleted only by an explicit, separately audited retention step after compact feature tables and checksums pass.

## 19. Decision rule for production expansion

Do not launch the entire structure-linked collection merely because smoke jobs exit zero.

Expand only when:

1. the PDB-experimental control panel has been audited;
2. results are stable to grid translation and reasonable spacing changes;
3. catalytic-region enrichment is not explained solely by pocket depth;
4. monomer and assembly-context differences are characterized;
5. runtime and storage are measured;
6. a deterministic candidate shard passes;
7. the production manifest is frozen and checksummed.

## 20. Phase 1 deliverables

- pinned installation script and checksums;
- server-local environment and software inventory;
- frozen explicit-probe manifest;
- frozen structure manifest with provenance and catalytic annotations;
- PDB-experimental-control smoke results;
- stratified benchmark results;
- compact per-protein and per-patch tables;
- failure ledger;
- independent audit JSON;
- PASS or FAIL gate artifact;
- updated RUNBOOK and append-only run history;
- production launch recommendation.
