# Layered Folddisco Score-Only Plan

Date: 2026-07-03
Project: enzyme scaffold search v2

## Goal

Run a multi-enzyme, layered Folddisco scoring stage before Pythia / GRASE-style pocket scoring.

The purpose is not to discard candidates at this stage. Sequence search and Foldseek are used for recall; Folddisco is used to measure whether each candidate structure contains local residue geometry similar to the seed pocket at several motif depths. The output is a scored table for downstream Pythia-pocket / ligand-pocket prioritization.

## Enzyme Scope

This stage advances three enzyme families together:

- Nylonase / nylon hydrolase candidates
- PETase / cutinase-like polyester hydrolase candidates
- Polyurethane / carbamate hydrolase candidates

Laccase is not part of this L0-L3 pass unless a separate laccase-specific motif design is added later.

## Core Decisions

1. Keep all candidates for now.
   - `filter_decision = NOT_FILTERED_FOR_PYTHIA` for scored rows.
   - `NO_HIT`, `PARTIAL_OR_LOOSE_HIT`, and `NOT_EVALUATED` must remain separate.

2. Folddisco native hits and project scores are different labels.
   - `FOLDDISCO_NATIVE_HIT`: Folddisco returned a target under the query settings.
   - `LAYER_SCORE`: project-defined score using matched residue count and motif geometry.
   - A missing download, missing structure, parse failure, or query error is `NOT_EVALUATED`, not a biological fail.

3. Sequence identity alone is not enough.
   - A candidate can be retained because it is recalled by sequence/Foldseek.
   - Folddisco score must use geometry evidence: matched count, query motif size, and RMSD / RMAD-like output from Folddisco.

4. Use loose per-layer cutoffs during score-only phase.
   - These cutoffs are support-score denominators, not hard filters.
   - L0: 2.0 A
   - L1: 2.5 A
   - L2: 4.0 A
   - L3: 5.0 A
   - Multi-layer variants use the loosest cutoff among included layers.

5. Respect CPU usage on the CPU server.
   - Total CPU consumption must not exceed 64 cores.
   - Default runner should use sequential Folddisco jobs with `--threads 16` or at most four concurrent 16-thread jobs.
   - Any future parallelization must explicitly enforce `max_parallel_jobs * threads_per_job <= 64`.

## Layer Semantics

- L0: minimal catalytic or experimentally central functional residues.
- L1: first-shell residues that shape the catalytic environment.
- L2: broader substrate-binding / pocket-lining residues.
- L3: context residues, processing/interface residues, loops, or secondary support positions.

L0 is not always available. For example, Nyl50 / Nyl10 pocket definitions are useful substrate-pocket evidence, but they are not treated as native catalytic L0 evidence.

## Search Variants

For de novo broad online / database search, broad motifs can saturate top-N limits. Therefore use multiple targeted queries:

- L0
- L0+L1
- L0+L2
- L0+L1+L2
- FULL where appropriate

For after-recall local scoring of sequence/Foldseek/Folddisco candidates, score more variants because the candidate universe is already bounded:

- L0
- L1
- L2
- L3
- L0+L1
- L0+L2
- L0+L1+L2
- FULL where available

Do not use L2-only or L3-only as the primary broad de novo search unless the family lacks a reliable L0 motif and the result is clearly labeled as auxiliary pocket evidence.

## Scoring Formula

For each query result row:

```text
coverage_score = matched_count / query_count
geometry_score = max(0, min(1, 1 - min_rmsd / loose_cutoff_for_variant))
layer_score = coverage_score * geometry_score
```

Fields to preserve:

- enzyme_class
- seed_id
- query_variant
- included_layers
- motif
- query_count
- matched_count
- coverage_score
- min_rmsd
- rmsd_loose_cutoff
- geometry_score
- layer_score
- folddisco_native_hit
- status
- filter_decision

Status vocabulary:

- `SCORED`: Folddisco row parsed and scored.
- `NO_HIT`: structure/query evaluated but no target row returned.
- `NOT_EVALUATED`: structure unavailable, query failed, parse failed, or input missing.
- `PARTIAL_OR_LOOSE_HIT`: optional downstream label for a native hit with incomplete motif coverage or loose geometry.
- `STRICT_PASS`: reserved for later strict filtering, not used as the main output in this phase.

## Motif Definitions v0.1

### PETase / Cutinase

L0 catalytic triad seeds:

| Seed | Structure | L0 |
|---|---|---|
| PET_01 | 6ILW | A160 Ser, A206 Asp, A237 His |
| PET_02 | 4EB0 | A165 Ser, A210 Asp, A242 His |
| PET_03 | 4CG1 | A130 Ser, A176 Asp, A208 His |
| PET_04 | 4WFI / 4WFJ | A176 Ser, A222 Asp, A254 His |
| PET_05 | 3VIS | A169 Ser, A215 Asp, A247 His |

L1 first shell:

| Seed | L1 |
|---|---|
| 6ILW | A87 Tyr, A161 Met, A236 Ser |
| 4EB0 | A95 Tyr, A166 Met, A241 Ser |
| 4CG1 | A60 Tyr, A131 Met, A207 Thr |
| 4WFI / 4WFJ | A106 Phe, A177 Met, A253 Thr |
| 3VIS | A99 Tyr, A170 Met, A246 Ser |

L2 pocket residues:

| Seed | L2 |
|---|---|
| 6ILW | A159 Trp, A185 Trp, A238 Ser |
| 4EB0 | A190 Trp, A243 Phe |
| 4CG1 | A155 Trp, A209 Phe |
| 4WFI / 4WFJ | A201 Trp, A255 Phe |
| 3VIS | A194 Trp, A248 Phe; PE4-near residues A99, A100, A170, A194, A217, A248 are optional context |

PETase L3 remains optional loop/stability context and should not be hard-filtered until exact positions are finalized.

### Nylonase

NylC WT / GYAQ / HP native catalytic motif:

| Seed | Structure | L0 |
|---|---|---|
| NylC_WT | 5XYO | A189 Lys, A219 Asn, A267 Thr |
| NylC_GYAQ | 5Y0M | A189 Lys, A219 Asn, A267 Thr |
| NylC_HP | ColabFold local UniRef30 rank001 | A189 Lys, A219 Asn, A267 Thr |

NylC support layers:

| Layer | Residues |
|---|---|
| L1 | A191 Asp, A192 Trp, A195 Thr, A235 Asn, A306 Asp, A308 Asp |
| L2 | A144 Val, A145 Ile, A146 Tyr, A192 Trp, A250 Tyr, A301 Phe, A305 Met |
| L3 | A134, A137, A267 processing/interface motif; T227-like positions only as context where applicable |

Nyl50 / Nyl10 auxiliary pocket motifs:

| Seed | Evidence class | Layer | Motif |
|---|---|---|---|
| Nyl50 | native pocket core, no T227 | L2 | A58 Ser, A89 Tyr, A96 Ile, A98 Leu |
| Nyl50 | auxiliary/interface context | L3 | A227 Thr, B41 Asn, B68 Ser, B103 Val, B104 Ile |
| Nyl10 | transferred from Nyl50 | L2 | A59, A90, A97, A99 |

Old NylB catalytic motif can be carried as a separate auxiliary nylonase seed:

| Layer | Motif |
|---|---|
| L0 | A112 Ser, A115 Lys, A215 Tyr |
| L1 | A170 Tyr |
| L2 | A186 Trp, A331 Trp, A181 Asp, A266 Asn |

### Polyurethane / Carbamate Hydrolase

UMG-SP1 / 8S7Z:

| Layer | Motif |
|---|---|
| L0 | A171 Ser, A176 Ser, A77 Lys |
| L1 | A168 Glu, A172 Asp, A119 Lys, A120 Thr, A121 Asn, A201 His |
| L2 | A80 Phe, A156 Ser, A178 Arg, A179 Ile, A215 Val, A216 Ile, A217 Gly, A298 Tyr, A306 Phe |
| L3 | Representative context residues from A152-A181 and A198-A219 loops, e.g. A152, A181, A198, A219 |

UMG-SP2 / 9FZW:

| Layer | Motif |
|---|---|
| L0 | A166 Ser, A185 Ser, A91 Lys |
| L1 | A182 Glu, A186 Asp, A134 Ser, A135 Asn, A142 Asp, A215 His |
| L2 | A103 Trp, A130 Phe, A143 Trp, A190 Ser, A228 Ser, A382 Trp |
| L3 | Representative context residues from A117-A146, A182-A194, and A212-A233 loops, e.g. A117, A146, A212, A233 |

Aes72 and AbPURase can be added after structure mapping. Until then, they should not be used as hard-numbered Folddisco motifs.

## Execution Plan

1. Generate `layered_motif_queries.tsv` from the motif definition table.
2. Run Folddisco local scoring over the existing after-recall candidate structures from sequence/Foldseek/Folddisco recall.
3. Use CPU-capped execution: default `threads=16`, sequential jobs, hard cap `max_total_cpu=64`.
4. Parse Folddisco native result rows into score tables.
5. Produce summaries by enzyme class, seed, layer variant, status, and top candidates.
6. Feed scored candidates into Pythia-pocket / GRASE-style scoring without throwing away low Folddisco-score candidates at this stage.

## Expected Outputs

Remote project path:

```text
/Dell/Dell14/bianht/enzyme_scaffold_search_v2
```

Output directory:

```text
results/layered_folddisco/
```

Expected files:

- `inputs/layered_motif_queries.tsv`
- `runs/after_recall/query_result_manifest.tsv`
- `runs/after_recall/scores/layered_folddisco_hit_scores.tsv`
- `runs/after_recall/scores/layered_folddisco_score_summary.tsv`
- `runs/after_recall/logs/nohup_layered_after_recall.log`

## Reporting Requirements

Every report must state:

- evaluated universe and denominator
- which enzyme classes were evaluated
- Folddisco native-hit counts separately from project scores
- `SCORED`, `NO_HIT`, and `NOT_EVALUATED` counts separately
- whether any row was actually filtered before Pythia; current answer should be no
- CPU cap used for the run

## Immediate TODO

- [x] Commit this plan to GitHub.
- [ ] Create CPU-server motif definition TSV.
- [ ] Create query builder.
- [ ] Create after-recall Folddisco runner with CPU cap enforcement.
- [ ] Create score parser and summary generator.
- [ ] Start one CPU-capped run for nylonase, PETase, and polyurethane/carbamate hydrolase.
- [ ] Integrate scores with Pythia-pocket candidate tables after scoring is available.
