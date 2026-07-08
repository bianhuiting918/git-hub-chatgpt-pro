# Folddisco Structure Recovery Runbook

This runbook records the exact recovery path used when Folddisco online hits have structures in the Folddisco/Foldseek result backend, but direct public downloads by accession fail. Keep this file lightweight: do not commit recovered PDB files, raw databases, or large hit tables.

## When to use this

Use this recovery path when a candidate appears in direct Folddisco results but the normal structure fetchers fail, for example:

- AFDB URL lookup such as `https://alphafold.ebi.ac.uk/files/AF-<accession>-F1-model_v*.pdb` returns nothing.
- RCSB/PDB or ESMAtlas direct lookup does not resolve the hit name.
- The target was emitted by Folddisco as an internal database hit, BFVD hit, AFDB subset hit, or metagenome structure hit.

Do not report these cases as biological failures. They are `NOT_EVALUATED` until the Folddisco backend structure route has been tried.

## Core lesson

Folddisco online result files contain enough information to request the matched structure from the backend. The key is that the API usually does not use the visible accession as the `id`. For non-PDB Folddisco databases, it uses the internal `dbkey` field from the `alis_*` result row.

Result API pattern:

```text
https://search.foldseek.com/api/result/folddisco/{ticket}?database={database}&id={id}
```

Where:

- `{ticket}` comes from the Folddisco online `tickets.tsv` for the query run.
- `{database}` comes from the result filename, e.g. `alis_BFVD_folddisco` -> `BFVD_folddisco`.
- `{id}` depends on the database:
  - for `pdb_folddisco`, use the target path/name from the first field;
  - for `BFVD_folddisco`, `afdb-proteome_folddisco`, `afdb50_folddisco`, `esm30_folddisco`, use the Folddisco internal `dbkey` field.

Observed `alis_*` field order:

```text
target_raw
matched_count
score
motif_rmsd
target_motif_residues
rotation
translation
coordinates
dbkey
query_motif
```

Example BFVD row, shortened:

```text
A0A345MUI4  3  32.1006  0.7977  A246,A243,A150  ...  79432  A189,A219,A267
```

The correct API id for this BFVD hit is `79432`, not `A0A345MUI4`.

## Target ID normalization

Normalize candidate IDs before matching Folddisco rows to the missing-structure manifest.

Recommended rules:

```text
AF-ACC-F1-model_vN(.pdb)       -> AFDB:ACC
AFDB:ACC or AFDB_ACC           -> AFDB:ACC
/opt/mmseqs-web/pdb100/xxxx.ent -> PDB:XXXX
MGYP...                        -> ESMATLAS:MGYP...
UPI...                         -> RAW:UPI...
plain accession-like token      -> AFDB:ACC and raw token aliases
```

Keep all aliases in the matching index. A target may appear as an accession in one table and as an AFDB-style model name in another.

## Server-side procedure used in enzyme scaffold v2

CPU project path:

```bash
cd /Dell/Dell14/bianht/enzyme_scaffold_search_v2
```

Relevant scripts created during the monomer-corrected nylonase rescue run:

```bash
python3 scripts/33_rescue_monomer_corrected_structures_from_folddisco_api.py --limit 0
python3 scripts/34_standardize_monomer_corrected_folddisco_api_pdbs.py
```

Main outputs:

```text
results/nylonase_monomer_corrected_20260708/folddisco_api_structure_rescue/folddisco_api_structure_rescue_manifest.tsv
results/nylonase_monomer_corrected_20260708/folddisco_api_structure_rescue/folddisco_api_structure_rescue_summary.tsv
results/nylonase_monomer_corrected_20260708/folddisco_api_structure_rescue/pdb_standardization_rescue_only_manifest.tsv
results/nylonase_monomer_corrected_20260708/specified_pocket_embedding/direct_folddisco_structure_manifest.tsv
results/nylonase_monomer_corrected_20260708/specified_pocket_embedding/direct_folddisco_specified_pocket_embedding_input.tsv
```

If these scripts are not present in a fresh checkout, reimplement them from the API rules above before falling back to ESMFold. The important logic is:

```python
api_id = target_raw if database.startswith("pdb") else dbkey
url = f"https://search.foldseek.com/api/result/folddisco/{ticket}?database={database}&id={api_id}"
```

Validate that the response is a structure file before marking it rescued. Some valid responses are served as `application/octet-stream`; inspect the text for PDB records such as `TITLE`, `ATOM`, or `HETATM`.

## PDB standardization requirement

Many Folddisco API structures are PDB-like but not fixed-column PDB. Example:

```text
ATOM      1  N    MET A   1       34.312 -14.523 -11.297  1.00 44.00           N
```

Bio.PDB and Pythia/Pythia-pocket can fail on these files because fixed-column residue parsing sees `line[22:26] == "A   "` instead of a residue number. After API download, standardize ATOM/HETATM rows into fixed-width PDB before building the Pythia manifest.

Audit labels:

```text
RESCUED                         API returned a valid structure and it is parseable or standardizable
SKIP_EXISTS                     an exact structure file already existed
NOT_RESCUED                     Folddisco occurrence was tried but no valid structure was returned
NO_FOLDDISCO_OCCURRENCE          missing target was not found in parsed alis rows; inspect ID normalization
NOT_EVALUATED_DOWNLOAD_FAILED    technical fetch failure, not biology
NOT_EVALUATED_PARSE_FAILED       structure exists but parsing/standardization failed
```

Never call download failure, parse failure, proxy-only status, or missing direct URL a biological negative.

## 2026-07-09 monomer-corrected nylonase checkpoint

This is the concrete case that motivated this runbook.

```text
direct Folddisco unique targets: 15,795
direct URL structures before API rescue: 11,634
missing before API rescue: 4,161
```

Folddisco API rescue over the 4,161 missing targets:

```text
RESCUED: 3,607
SKIP_EXISTS: 19
NO_FOLDDISCO_OCCURRENCE: 229
NOT_RESCUED: 306
```

After rebuilding the direct Folddisco structure manifest:

```text
structures OK: 15,262 / 15,795
still not OK: 533
```

PDB standardization over API-rescued structure paths:

```text
total Folddisco API structure paths: 3,626
already parseable: 637
standardized: 2,891
standardization failed: 98
```

## Next-step rule

Only after this Folddisco API recovery route, direct accession remapping, ESMAtlas/AFDB/RCSB lookup, and homology/proxy search are exhausted should the remaining exact-sequence candidates enter ESMFold. Mark proxy structures as `HOMOLOGY_PROXY`; do not mix them with exact target structures in RMSD, pocket embedding, or final evidence tables.