#!/usr/bin/env python3
"""Build one PA66-L2 system from the rare unbiased NylC true-Thr267 NAC."""
import hashlib
import json
import pathlib
import shutil
import sys

SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from prepare_candidates import prepare_one

TASK = pathlib.Path("/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723")
CANDIDATE_ID = "nylc_C18_trueT267_freeGS"
SOURCE_ROOT = TASK / "candidates" / CANDIDATE_ID
SOURCE_GRO_SHA = "0f0c4932e591084feba457819daf1c0bf5cd6ec249778f0673005917b82c4c22"
L2_ITP = pathlib.Path("/work/home/acshdt1dks/nylon_pa66_scnet_20260708/nylc_gyaq_pa66_l2_nac_qmmm_20260723/inputs/parameterized/PA66_L2_GMX.itp")
L2_ITP_SHA = "b0e753c60fd4b71c282d21cc6106a15e73d91d12a20d80e92dd01516162eb301"
BRANCH = pathlib.Path("/work/home/acshdt1dks/nylon_pa66_scnet_20260708/apo_gate_l4_three_carbonyl_20260715/cyclic_gate_nac_20260719/branches/nylc_gyaq/C18")


def sha256(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    selection_path = SOURCE_ROOT / "source_manifest.json"
    source_path = SOURCE_ROOT / "source.gro"
    selection = json.loads(selection_path.read_text())
    if selection["status"] != "VERIFIED_RARE_UNRESTRAINED_NAC_GS_CANDIDATE":
        raise RuntimeError("source is not the verified rare unbiased NAC candidate")
    if sha256(source_path) != SOURCE_GRO_SHA:
        raise RuntimeError("source GRO SHA mismatch")
    if sha256(L2_ITP) != L2_ITP_SHA:
        raise RuntimeError("audited L2 ITP SHA mismatch")

    input_dir = TASK / "inputs" / CANDIDATE_ID
    input_path = input_dir / "source.gro"
    build_dir = SOURCE_ROOT / "build"
    audit_dir = SOURCE_ROOT / "audit"
    if input_path.exists() or build_dir.exists() or audit_dir.exists():
        raise FileExistsError("refusing to overwrite free-GS L2 inputs/build/audit")
    input_dir.mkdir(parents=True)
    shutil.copy2(source_path, input_path)

    candidate = {
        "id": CANDIDATE_ID,
        "status": "ACTIVE",
        "enzyme": "NylC-GYAQ",
        "site": "C18_trueThr267",
        "branch_dir": str(BRANCH),
        "source_itp": "l4_nfree_GMX.itp",
        "source_molecule_type": "l4_nfree",
        "source_atom_count": 160,
        "ligand_first_global_atom": 10273,
        "reactive_local_atoms": {
            "carbonyl_c": 25,
            "carbonyl_o": 26,
            "amide_n": 24,
        },
        "reactive_global_atoms": {
            "carbonyl_c": 10297,
            "carbonyl_o": 10298,
            "amide_n": 10296,
            "thr_og1": 8961,
        },
        "source_geometry": {
            "distance_nm": selection["geometry"]["true_thr267_og1_to_carbonyl_c_nm"],
            "angle_deg": selection["geometry"]["true_thr267_attack_angle_deg"],
        },
        "gate_definition": "NylC residues 261-266; Thr267 excluded",
    }
    audit = prepare_one(candidate, TASK, L2_ITP)
    audit["source_selection_manifest"] = {
        "path": str(selection_path),
        "sha256": sha256(selection_path),
        "status": selection["status"],
        "scientific_caveat": selection["scientific_caveat"],
        "one_ns_followup": selection["one_ns_followup"],
    }
    audit_path = audit_dir / "build_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    if audit["state"] != "BUILD_PASS":
        raise SystemExit("L4-to-L2 build audit did not pass")
    print(json.dumps(audit, sort_keys=True))


if __name__ == "__main__":
    main()
