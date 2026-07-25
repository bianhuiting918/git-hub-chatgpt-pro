import hashlib
import json
import math
import sys
from pathlib import Path

import pytest

FLOW = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(FLOW / "scripts"))

import prepare_nylc_m1_nac_ensemble as ensemble_build


def _gro_atom(resid, resname, atom_name, atom_index, xyz):
    x, y, z = xyz
    return (
        f"{resid:5d}{resname:<5s}{atom_name:>5s}{atom_index:5d}"
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
    )


def _source_gro(path, thr_resname="THR"):
    carbon = (0.300, 0.000, 0.000)
    oxygen = (
        carbon[0] + 0.120 * math.cos(math.radians(75.0)),
        carbon[1] + 0.120 * math.sin(math.radians(75.0)),
        0.000,
    )
    path.write_text(
        "synthetic NylC source\n"
        "4\n"
        + _gro_atom(267, thr_resname, "OG1", 1, (0.000, 0.000, 0.000))
        + "\n"
        + _gro_atom(356, "L2", "C12", 2, carbon)
        + "\n"
        + _gro_atom(356, "L2", "O2", 3, oxygen)
        + "\n"
        + _gro_atom(356, "L2", "N3", 4, (0.300, -0.130, 0.000))
        + "\n"
        + "   4.00000   4.00000   4.00000\n"
    )


def _authority(tmp_path, source):
    selection = tmp_path / "ensemble_selection.json"
    source_sha = hashlib.sha256(source.read_bytes()).hexdigest()
    selection.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "selected": [
                    {
                        "candidate_id": "nac_evt01_time1000ps",
                        "event_id": "event_0001",
                        "time_ps": 1000.0,
                        "source_gro": str(source),
                        "source_gro_sha256": source_sha,
                    }
                ],
            }
        )
        + "\n"
    )
    return {
        "ensemble_selection": str(selection),
        "reactive_global_index1_m0": {
            "thr267_og1": 1,
            "l2_c": 2,
            "l2_o": 3,
            "l2_n": 4,
        },
        "reactive_atom_identity_m0": {
            "thr267_og1": [267, "THR", "OG1"],
            "l2_c": [356, "L2", "C12"],
            "l2_o": [356, "L2", "O2"],
            "l2_n": [356, "L2", "N3"],
        },
        "nac_gate": {
            "distance_max_nm": 0.35,
            "angle_min_deg": 95.0,
            "angle_max_deg": 115.0,
        },
        "gate_residues": [261, 262, 263, 264, 265, 266],
        "excluded_gate_residues": [267],
    }


def test_build_candidate_freezes_source_and_microstate_contract(tmp_path, monkeypatch):
    source = tmp_path / "source.gro"
    _source_gro(source)
    authority = _authority(tmp_path, source)
    candidate_dir = tmp_path / "candidate"

    def fake_materialize(source_gro, output_dir, frozen_authority, selected):
        assert source_gro == source
        assert selected["candidate_id"] == "nac_evt01_time1000ps"
        return {
            "charge_delta_e": 0.0,
            "active_chain_atom_count": 1324,
            "full_atom_count_unchanged": True,
            "reaction_geometry": {
                "attack_distance_nm": 0.300,
                "attack_angle_deg": 105.0,
                "joint_nac": True,
            },
            "minimum_chain_rest_distance_nm": 0.10,
            "minimum_heavy_atom_contact_nm": 0.20,
            "grompp_maxwarn": 0,
        }

    monkeypatch.setattr(ensemble_build, "_materialize_candidate", fake_materialize)
    audit = ensemble_build.build_candidate(source, candidate_dir, authority)

    assert audit["microstate"] == {
        "Thr267_Nalpha": "NH2_neutral_proxy",
        "Thr267_Ogamma": "OH",
        "Asp306": "ASH_HD2_neutral",
        "Asp308": "ASP_minus",
    }
    assert audit["charge_delta_e"] == 0.0
    assert audit["reaction_geometry"]["joint_nac"] is True
    assert audit["gate_definition"] == "NylC residues 261-266; Thr267 excluded"
    assert audit["source_gro_sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert json.loads((candidate_dir / "build_audit.json").read_text()) == audit
    assert (candidate_dir / "PASS.json").is_file()


def test_build_candidate_refuses_existing_target(tmp_path):
    source = tmp_path / "source.gro"
    _source_gro(source)
    authority = _authority(tmp_path, source)
    candidate_dir = tmp_path / "candidate"
    candidate_dir.mkdir()

    with pytest.raises(FileExistsError):
        ensemble_build.build_candidate(source, candidate_dir, authority)


def test_build_candidate_rejects_source_outside_selection_manifest(tmp_path):
    source = tmp_path / "source.gro"
    other = tmp_path / "other.gro"
    _source_gro(source)
    _source_gro(other)
    authority = _authority(tmp_path, source)

    with pytest.raises(ValueError, match="not present in the frozen selection"):
        ensemble_build.build_candidate(other, tmp_path / "candidate", authority)


def test_build_candidate_rejects_wrong_reactive_atom_identity(tmp_path):
    source = tmp_path / "source.gro"
    _source_gro(source, thr_resname="SER")
    authority = _authority(tmp_path, source)

    with pytest.raises(ValueError, match="reactive identity mismatch"):
        ensemble_build.build_candidate(source, tmp_path / "candidate", authority)


def test_build_array_maps_all_twelve_frozen_candidates_without_nested_srun():
    batch = FLOW / "slurm" / "run_nylc_m1_ensemble_build_array.sbatch"
    text = batch.read_text()

    assert "#SBATCH --array=0-11%4" in text
    assert "selection_job_61813011/ensemble_selection.json" in text
    assert 'srun -n 1' not in text
    assert 'build_job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}' in text
