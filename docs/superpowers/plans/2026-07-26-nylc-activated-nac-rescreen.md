# NylC activated Thr267 NAC rescreen implementation plan

> **For Codex:** Execute this plan in order. Keep topology/parameter provenance, technical completion, and scientific gates separate. All compute runs on SCNet; GitHub receives only code, documentation, manifests, and compact text audits.

**Goal:** Build the atom-count-conserving NylC A1 microstate `Thr267 NalphaH3+ / Ogamma- / Asp306H / Asp308-` from three independently selected M1 NAC medoids, run audited fixed-topology MM screening, and hand only reproducible stable NAC medoids to DFTB3/3OB-3-1 QM/MM.

**Architecture:** Freeze the parent coordinates and microstate change in a machine-readable authority manifest. Generate and independently audit a local A1 residue patch before modifying any full system. A parsed-topology transformer transfers the existing HG1 from Oγ to Nα, preserves every other atom, applies the audited local patch, and emits structural/topological audits. A candidate-isolated SCNet workflow runs build -> preflight -> EM -> restrained preparation -> 3 x 1 ns fully unrestrained NPT. An independent final audit ranks recurrent NAC clusters; a separate QM/MM preflight consumes only passing medoids.

**Toolchain:** Python 3, pytest, GROMACS 2022.1/DCU, AmberTools/Amber 2018 (`antechamber`, `parmchk2`, `resp`, `respgen`, `tleap`, `sqm`, `sander`), Slurm, JSON/TSV/SHA256.

---

## Task 1: Freeze the A1 authority and parameter-source gate

**Files**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/manifests/nylc_a1_activated_nac.authority.json`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/audit_nylc_a1_parameter_source.py`
- Test: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_a1_parameter_source_contract.py`

**Step 1: Write the failing contract test**

The test must assert that the authority manifest freezes:

- final audit job `61841413`;
- exact three parent labels, source paths, times, and SHA256 values;
- parent M1 and target A1 protonation strings;
- one removed bond `Thr267:OG1-HG1`;
- one added bond `Thr267:N-HG1`;
- unchanged atom count and unchanged total system charge;
- unchanged `Asp306H` and `Asp308-`;
- force field `amber99sb-ildn`;
- explicit refusal when neither a native residue patch nor a generated-and-audited patch is supplied.

Example contract:

```python
def test_authority_freezes_atom_conserving_transfer(authority):
    assert authority["target_microstate"] == {
        "Thr267_Nalpha": "NH3+",
        "Thr267_Ogamma": "O-",
        "Asp306": "ASH",
        "Asp308": "ASP",
    }
    assert authority["topology_delta"]["removed_bonds"] == [["Thr267:OG1", "Thr267:HG1"]]
    assert authority["topology_delta"]["added_bonds"] == [["Thr267:N", "Thr267:HG1"]]
    assert authority["invariants"]["atom_count_delta"] == 0
    assert authority["invariants"]["system_charge_delta_e"] == 0
```

**Step 2: Run the test and confirm it fails**

SCNet command:

```bash
cd /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/repo
python3 -m pytest -q workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_a1_parameter_source_contract.py
```

Expected: FAIL because the manifest and audit module do not yet exist.

**Step 3: Implement the minimum authority audit**

`audit_nylc_a1_parameter_source.py` must parse JSON and return a nonzero exit code unless:

- every parent coordinate SHA matches;
- the requested native patch exists, or a generated patch contains provenance for molecule definition, net charge, charge method, atom types, bonded parameters, tool versions, input/output hashes, and validation status;
- the patch residue charge is exactly zero within (10^{-5}) e;
- the patch exposes N-HG1 and lacks OG1-HG1;
- no nonlocal atom is modified.

The first checked-in manifest must set `parameter_status: NOT_EVALUATED_PENDING_GENERATION`; it must not claim MM readiness.

**Step 4: Re-run tests**

Expected: PASS for schema/invariant tests and an intentional refusal fixture for missing parameter provenance.

**Step 5: Commit**

```text
docs(nylc): freeze activated Thr NAC authority
```

## Task 2: Generate and audit the A1 local residue patch on SCNet

**Files**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/prepare_nylc_a1_parameter_model.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/run_nylc_a1_parameterization.sh`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/audit_nylc_a1_patch.py`
- Test: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_a1_patch_contract.py`
- Update after remote audit: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/manifests/nylc_a1_activated_nac.authority.json`

**Step 1: Add failing tests**

Tests must verify the capped model contains the N-terminal Thr fragment plus the next-peptide cap, has formal charge zero for the A1 Thr residue, and preserves the atom correspondence `N/H1/H2/HG1/CA/CB/OG1/CG2/C/O`.

Reject all of these:

- relabeling OG1 from hydroxyl to alkoxide without a parameter source;
- deleting HG1 or adding a new atom;
- changing Asp306/Asp308;
- importing GAFF atom types into the protein without an explicit crosswalk and bonded-term audit;
- charge rounding that changes residue or system integral charge.

**Step 2: Build the isolated capped model**

Use a graph-derived model from the authoritative evt18 medoid. Place HG1 tetrahedrally at Nα before charge generation. Keep a machine-readable atom map back to the full-system atom IDs.

**Step 3: Run reproducible AmberTools parameter generation**

SCNet script header:

```bash
#!/usr/bin/env bash
set -euo pipefail
module purge
module load amber/2018-hpcx-gcc-7.3.1
for exe in antechamber parmchk2 resp respgen tleap sqm; do
  command -v "$exe" >/dev/null
done
```

Primary classical-screening route:

- generate AM1-BCC charges with AmberTools `sqm/antechamber` for the closed-shell, net-zero capped A1 model;
- derive missing bonded terms with `parmchk2`;
- export exact command lines, AmberTools version, input SHA256, output SHA256, formal charges, fitted charges, atom types, and `parmchk2` warnings;
- map only the A1 local residue terms into a named patch; do not replace the rest of Amber99SB-ILDN.

This AM1-BCC patch is a screening Hamiltonian, not a claim that A1 is a native Amber99SB-ILDN residue. The later DFTB3 QM/MM optimization is the scientific electronic-structure gate.

**Step 4: Audit before full-system use**

`audit_nylc_a1_patch.py` must check:

- residue charge 0.00000 e within tolerance;
- Nα has three H bonds; OG1 has no H bond;
- expected valences and no duplicate bonded terms;
- every required bond/angle/dihedral/LJ term resolves;
- no `ATTN, need revision`, untyped atom, or unresolved `parmchk2` placeholder remains;
- `tleap` can load the capped model and report the expected charge;
- independent GROMACS single-residue `grompp` succeeds after conversion;
- vacuum single-point energy is finite in both Amber and GROMACS; the audit records the values but does not require equality across engines.

If any check fails, set `parameter_status: FAIL_A1_PARAMETER_AUDIT` and route the three parents to direct QM/MM preflight; do not submit classical A1 dynamics.

**Step 5: Update authority manifest with hashes and audit status**

Only after the remote audit passes may the manifest say `PASS_A1_SCREENING_PATCH`.

**Step 6: Commit compact artifacts**

Commit scripts, manifest, parameter text small enough for GitHub, and compact audit JSON. Do not commit full-system topology, trajectories, checkpoints, or credentials.

```text
feat(nylc): add audited activated Thr screening patch
```

## Task 3: Implement the atom-count-conserving A1 full-system builder

**Files**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/build_nylc_a1_activated_nac.py`
- Test: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_build_nylc_a1_activated_nac.py`
- Reuse: `scripts/prepare_nylc_nalpha_h2_ash306.py` parsing helpers where safe

**Step 1: Write failing graph-transformation tests**

Construct a minimal ITP/GRO fixture and assert:

```python
assert before.atom_count == after.atom_count
assert after.bonded("Thr267:N", "Thr267:HG1")
assert not after.bonded("Thr267:OG1", "Thr267:HG1")
assert after.residue_charge("Thr267") == pytest.approx(0.0, abs=1e-5)
assert after.system_charge == pytest.approx(before.system_charge, abs=1e-5)
assert unchanged_atoms(before, after, exclude={"Thr267:HG1"})
```

Also test refusal for wrong residue number, duplicate Thr candidates, wrong parent microstate, mismatched SHA, and absent parameter audit PASS.

**Step 2: Implement parsed topology transformation**

The builder must:

- locate the unique Thr267 N/OG1/HG1 by residue and bonded graph;
- verify the input is M1: N-H1/H2, OG1-HG1, Asp306 HD2, Asp308 deprotonated;
- calculate a tetrahedral HG1 coordinate from the N-CA/H1/H2 frame;
- update only HG1 coordinates;
- apply the audited A1 atom-type/charge/bonded patch;
- preserve molecule ordering, periodic box, L2 substrate, water, and ions;
- emit `build_audit.json`, `atom_mapping.tsv`, `bond_diff.tsv`, and SHA256 records.

**Step 3: Add chemistry and geometry audits**

Record:

- atom and residue counts;
- total and per-residue charge;
- N-HG1 and nearest nonbonded distances;
- Thr Oγ–substrate C distance and O-C-Oγ angle;
- gate 261–266 opening only; Thr267 must not be included in the gate group;
- Asp306/Asp308 proton identities;
- minimum nonbonded distance with exclusions based on bonded graph.

**Step 4: Run unit tests and a no-submit remote dry run on evt18**

Expected: tests PASS; dry run writes only into a fresh candidate-specific build directory and performs no MD.

**Step 5: Commit**

```text
feat(nylc): build atom-conserving A1 full systems
```

## Task 4: Add candidate-isolated preflight and MM workflow

**Files**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/prepare_nylc_a1_ensemble.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_a1_build_array.sbatch`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_a1_mm_array.sbatch`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/submit_nylc_a1_pipeline.sh`
- Test: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_a1_ensemble_contract.py`
- Update: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/RUNBOOK.md`

**Step 1: Write failing workflow-contract tests**

Assert exactly 3 parents x 3 seeds, fresh candidate/seed directories, dependency-preserving Slurm submission, and skip/lock guards. Assert no stage can overwrite M1 or superseded outputs.

**Step 2: Generate stages**

For each parent:

1. A1 build and topology audit;
2. GROMACS `grompp` preflight;
3. steepest-descent EM;
4. low-temperature heavy-atom restrained relaxation;
5. gradual heating;
6. staged restraint release;
7. 1 ns fully unrestrained NPT for three deterministic seeds.

The final 1 ns window is a fresh continuation with zero positional restraints. Restraint stages are preparation only.

**Step 3: Add run history and failure isolation**

Every job appends:

- UTC time;
- git commit;
- script and parameters;
- input/output paths and hashes;
- Slurm ID and array index;
- exit code;
- technical status;
- scientific status or `NOT_EVALUATED_*`.

A failed parent/replica must not block unrelated array elements.

**Step 4: Read-only queue audit, then submit**

Before submission inspect `squeue`, `sacct`, partition availability, and existing output locks. Request a reasonable number of CPUs after queue inspection; no Dell compute.

**Step 5: Commit**

```text
feat(nylc): add activated NAC MM ensemble workflow
```

## Task 5: Independently audit and rank the A1 ensemble

**Files**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/audit_nylc_a1_ensemble.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_a1_final_audit.sbatch`
- Test: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_a1_final_audit_contract.py`
- Update: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/RUNBOOK.md`

**Step 1: Write failing analysis tests**

Use synthetic trajectories/tables to test:

- burn-in exclusion;
- NAC criterion distance <= 0.35 nm and angle 95–115 degrees;
- gate opening computed from residues 261–266 only;
- per-replica and pooled occupancy;
- contiguous residence with explicit frame spacing;
- recurrent cluster membership across at least two replicas;
- substrate-bound condition;
- LINCS/SETTLE/NaN/FATAL classification;
- separation of technical PASS from scientific PASS.

**Step 2: Implement independent audit**

A parent passes only if all required numerical/stability checks pass, all three replicas remain bound, at least two replicas contain post-burn-in NAC frames, and a NAC cluster recurs across at least two replicas.

Rank passing recurrent-cluster medoids by:

1. replica recurrence;
2. NAC occupancy/residence;
3. NAC geometry;
4. same-Hamiltonian MM potential-energy distribution as a descriptor only.

Do not call MM energy a barrier or free energy.

**Step 3: Emit compact handoff**

Write `FINAL_A1_AUDIT.json`, `failure_ledger.tsv`, `qmmm_eligible_a1_medoids.json`, and text checksums. Preserve every failure reason.

**Step 4: Commit**

```text
feat(nylc): audit activated NAC ensemble
```

## Task 6: Prepare the DFTB3/3OB-3-1 QM/MM preflight

**Files**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/prepare_nylc_a1_qmmm_preflight.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_a1_qmmm_preflight.sbatch`
- Test: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_a1_qmmm_preflight_contract.py`
- Update: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/RUNBOOK.md`

**Step 1: Write failing QM-region tests**

Assert that the QM region includes:

- Thr267 Nα, H1, H2, transferred HG1, Cα/Cβ/Oγ;
- reactive L2 carbonyl C/O and attacking geometry atoms;
- required adjacent substrate atoms;
- Asp306 and/or explicit waters only when selected by the frozen region definition;
- link atoms and total QM charge consistent with A1.

Assert that the proton coordinate is not restrained during the decisive optimization.

**Step 2: Build two region sizes**

Prepare a core region and a water-network sensitivity region using DFTB3/3OB-3-1. Freeze exact atom lists, charge/multiplicity, boundary bonds, input hashes, and executable/runtime versions.

**Step 3: Run smoke checks then unconstrained optimization**

For each eligible A1 medoid:

1. input/charge/link-atom validation;
2. short finite-energy smoke calculation;
3. unconstrained optimization with the activation proton free;
4. post-optimization NAC and microstate audit.

Reject collapsed structures with explicit labels such as `FAIL_QMMM_A1_PROTON_RETURN`, `FAIL_QMMM_NAC_COLLAPSE`, or `NOT_EVALUATED_TECHNICAL`.

**Step 4: Label surviving structures correctly**

A survivor is `GS_act`, not the resting-state ground state. Any later Step1 TS/PMF barrier is conditional on prior activation and excludes the activation free energy.

**Step 5: Commit and verify**

Run all workflow tests and perform an independent remote artifact audit before reporting completion.

```bash
python3 -m pytest -q workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests
```

Expected: PASS. Then commit:

```text
feat(nylc): add activated-state QM/MM preflight
```

---

## Execution checkpoints

1. **Parameter checkpoint:** no classical A1 MD submission until `PASS_A1_SCREENING_PATCH`.
2. **Build checkpoint:** no EM until atom count, total charge, valence, minimum-distance, and `grompp` audits pass.
3. **Scientific MM checkpoint:** restrained stages never count; only the fully unrestrained 1 ns windows determine A1 NAC stability.
4. **QM/MM checkpoint:** only recurrent A1 NAC medoids enter QM/MM; an unconstrained optimization must retain both A1 and NAC before the structure is called `GS_act`.
5. **Reporting checkpoint:** technical completion, scientific PASS/FAIL, and `NOT_EVALUATED_*` are always reported separately.
