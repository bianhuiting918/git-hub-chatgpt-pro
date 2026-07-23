import importlib.util
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
FLOW = HERE.parent
PREP = FLOW / "scripts" / "prepare_true_thr267_recapture_pilot.py"
SBATCH = FLOW / "slurm" / "run_true_thr267_recapture_pilot_array.sbatch"
AUDIT = FLOW / "scripts" / "audit_true_thr267_recapture_pilot.py"
SBATCH2 = FLOW / "slurm" / "run_true_thr267_recapture_response2_array.sbatch"
SBATCH3 = FLOW / "slurm" / "run_true_thr267_recapture_response3_array.sbatch"
SBATCH4 = FLOW / "slurm" / "run_true_thr267_recapture_response4_array.sbatch"
SBATCH5 = FLOW / "slurm" / "run_true_thr267_recapture_response5_c18.sbatch"
SBATCH6 = FLOW / "slurm" / "run_true_thr267_recapture_response6_c18.sbatch"
SBATCH7 = FLOW / "slurm" / "run_true_thr267_recapture_response7_c18.sbatch"
SBATCH8 = FLOW / "slurm" / "run_true_thr267_recapture_response8_c18.sbatch"
SBATCH9 = FLOW / "slurm" / "run_true_thr267_recapture_response9_c18.sbatch"
SBATCH10 = FLOW / "slurm" / "run_true_thr267_recapture_response10_c18.sbatch"
SBATCH11 = FLOW / "slurm" / "run_true_thr267_recapture_response11_c18.sbatch"


def load_audit():
    spec = importlib.util.spec_from_file_location("audit_recapture", AUDIT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_prep():
    spec = importlib.util.spec_from_file_location("prepare_recapture", PREP)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_recapture_mdp_is_bounded_and_uses_fresh_velocities():
    m = load_prep()
    mdp = m.make_mdp("C18", seed=181267)
    assert "nsteps                   = 50000" in mdp
    assert "continuation             = no" in mdp
    assert "gen-vel                  = yes" in mdp
    assert "gen-seed                 = 181267" in mdp
    assert "pull-coord1-start        = yes" in mdp
    assert "pull-coord1-rate         = -0.004000" in mdp
    assert "pull-coord1-k            = 100" in mdp
    assert "pull-coord2-k            = 20" in mdp
    assert "pull-group1-name         = Thr_OG1" in mdp
    assert "pull-group2-name         = L4_C18" in mdp
    assert "Gate" not in mdp


def test_corrected_index_replaces_only_wrong_thr_group():
    m = load_prep()
    text = "[ Gate ]\n8870 8896\n[ Thr_OG1 ]\n8896\n[ L4_C18 ]\n10297\n"
    corrected = m.correct_index(text)
    assert "[ Thr_OG1 ]\n8961\n" in corrected
    assert "[ Gate ]\n8870 8896\n" in corrected
    assert corrected.count("8961") == 1


def test_slurm_requires_preflight_and_true_thr267_contract():
    text = SBATCH.read_text(encoding="utf-8")
    assert "#SBATCH --gres=dcu:1" in text
    assert "nylc_C18_trueT267_recapture" in text
    assert "nylc_C23_trueT267_recapture" in text
    assert "PASS_GROMPP_TRUE_THR267" in text
    assert "mpirun -np 1" in text
    assert "SUPERSEDED_WRONG_NUCLEOPHILE_IDENTITY" in text
    assert "8896" not in text
    assert r"\\${" not in text
    assert 'candidate="${candidates[${SLURM_ARRAY_TASK_ID:?}]}\"' in text


def test_audit_distinguishes_numerical_failure_from_response_failure():
    text = AUDIT.read_text(encoding="utf-8")
    assert "FAIL_NUMERICAL" in text
    assert "FAIL_APPROACH_RESPONSE_NUMERICALLY_STABLE" in text
    assert "PASS_TECHNICAL_BOUNDED_APPROACH" in text


def test_response2_is_stronger_but_moves_reference_more_slowly():
    m = load_prep()
    mdp = m.make_mdp("C18", seed=181268, protocol="response2")
    assert "pull-coord1-rate         = -0.002000" in mdp
    assert "pull-coord1-k            = 500" in mdp
    assert "pull-coord2-k            = 50" in mdp
    assert "nsteps                   = 50000" in mdp


def test_response2_slurm_uses_separate_outputs_and_preflight():
    text = SBATCH2.read_text(encoding="utf-8")
    assert "recapture_response2" in text
    assert "#SBATCH --gres=dcu:1" in text
    assert "PASS_GROMPP_TRUE_THR267" in text
    assert "nylc_C18_trueT267_recapture" in text
    assert "nylc_C23_trueT267_recapture" in text
    assert "pilot.gro" not in text
    assert "response2.gro" in text
    assert '-deffnm "$work/response2"' in text
    assert "audit_true_thr267_recapture_pilot.py" in text
    assert "audit_true_thr267_recapture_response2.py" not in text


def test_response3_is_short_high_stiffness_probe():
    m = load_prep()
    mdp = m.make_mdp("C18", seed=181269, protocol="response3")
    assert "nsteps                   = 12500" in mdp
    assert "pull-coord1-rate         = -0.004000" in mdp
    assert "pull-coord1-k            = 2000" in mdp
    assert "pull-coord2-k            = 100" in mdp


def test_response3_slurm_uses_005_response_gate():
    text = SBATCH3.read_text(encoding="utf-8")
    assert "recapture_response3" in text
    assert '-deffnm "$work/response3"' in text
    assert "--min-response 0.05" in text
    assert "audit_true_thr267_recapture_pilot.py" in text


def test_response4_continues_from_response3_checkpoint():
    m = load_prep()
    mdp = m.make_mdp("C18", seed=181270, protocol="response4")
    assert "nsteps                   = 50000" in mdp
    assert "continuation             = yes" in mdp
    assert "gen-vel                  = no" in mdp
    assert "gen-seed" not in mdp
    assert "pull-coord1-rate         = -0.004000" in mdp
    assert "pull-coord1-k            = 2000" in mdp


def test_response4_slurm_audits_against_response3_start():
    text = SBATCH4.read_text(encoding="utf-8")
    assert "recapture_response4" in text
    assert '-deffnm "$work/response4"' in text
    assert "--reference-gro" in text
    assert "recapture_response3/response3.gro" in text
    assert "--min-response 0.15" in text


def test_response5_continues_c18_and_strengthens_angle_only():
    m = load_prep()
    mdp = m.make_mdp("C18", seed=181271, protocol="response5")
    assert "continuation             = yes" in mdp
    assert "nsteps                   = 50000" in mdp
    assert "pull-coord1-k            = 2000" in mdp
    assert "pull-coord2-k            = 500" in mdp
    assert m.PROTOCOLS["response5"]["parent_stem"] == "response4"


def test_response5_slurm_is_c18_only_and_checkpointed():
    text = SBATCH5.read_text(encoding="utf-8")
    assert "#SBATCH --array" not in text
    assert "nylc_C18_trueT267_recapture" in text
    assert "nylc_C23_trueT267_recapture" not in text
    assert '-deffnm "$work/response5"' in text
    assert "recapture_response4/response4.gro" in text
    assert "--min-response 0.20" in text


def test_response6_continues_c18_from_response5():
    m = load_prep()
    mdp = m.make_mdp("C18", seed=181272, protocol="response6")
    assert "continuation             = yes" in mdp
    assert "pull-coord1-k            = 2000" in mdp
    assert "pull-coord2-k            = 500" in mdp
    assert m.PROTOCOLS["response6"]["parent_stem"] == "response5"


def test_response6_slurm_uses_response5_reference():
    text = SBATCH6.read_text(encoding="utf-8")
    assert "nylc_C18_trueT267_recapture" in text
    assert '-deffnm "$work/response6"' in text
    assert "recapture_response5/response5.gro" in text
    assert "--min-response 0.15" in text


def test_contact_audit_reports_ligand_protein_minimum():
    m = load_audit()
    atoms = {
        1: {"resid": 10, "resname": "ALA", "atomname": "CA", "xyz_nm": (0.0, 0.0, 0.0)},
        2: {"resid": 356, "resname": "UNL", "atomname": "C18", "xyz_nm": (0.2, 0.0, 0.0)},
    }
    result = m.minimum_ligand_protein_distance(atoms, (1.0, 1.0, 1.0))
    assert result["distance_nm"] == 0.2
    assert result["protein_global_index"] == 1
    assert result["ligand_global_index"] == 2


def test_response7_continues_c18_from_response6():
    m = load_prep()
    mdp = m.make_mdp("C18", seed=181273, protocol="response7")
    assert "continuation             = yes" in mdp
    assert "pull-coord1-rate         = -0.004000" in mdp
    assert m.PROTOCOLS["response7"]["parent_stem"] == "response6"


def test_response7_slurm_uses_response6_reference():
    text = SBATCH7.read_text(encoding="utf-8")
    assert '-deffnm "$work/response7"' in text
    assert "recapture_response6/response6.gro" in text
    assert "--min-response 0.15" in text


def test_response8_slows_final_approach():
    m = load_prep()
    mdp = m.make_mdp("C18", seed=181274, protocol="response8")
    assert "continuation             = yes" in mdp
    assert "pull-coord1-rate         = -0.001500" in mdp
    assert "pull-coord1-k            = 2000" in mdp
    assert m.PROTOCOLS["response8"]["parent_stem"] == "response7"


def test_response8_slurm_uses_response7_reference():
    text = SBATCH8.read_text(encoding="utf-8")
    assert '-deffnm "$work/response8"' in text
    assert "recapture_response7/response7.gro" in text
    assert "--min-response 0.05" in text


def test_response9_targets_final_040_nm_gate():
    m = load_prep()
    mdp = m.make_mdp("C18", seed=181275, protocol="response9")
    assert "pull-coord1-rate         = -0.002000" in mdp
    assert m.PROTOCOLS["response9"]["parent_stem"] == "response8"
    text = SBATCH9.read_text(encoding="utf-8")
    assert '-deffnm "$work/response9"' in text
    assert "recapture_response8/response8.gro" in text
    assert "--min-response 0.03" in text
    assert "--max-end-distance 0.40" in text


def test_audit_supports_maximum_end_distance_gate():
    text = AUDIT.read_text(encoding="utf-8")
    assert "--max-end-distance" in text
    assert "maximum_allowed_end_distance_nm" in text


def test_response10_uses_fixed_035_nm_target():
    m = load_prep()
    mdp = m.make_mdp("C18", seed=181276, protocol="response10")
    assert "pull-coord1-start        = no" in mdp
    assert "pull-coord1-init         = 0.350000" in mdp
    assert "pull-coord1-rate         = 0.000000" in mdp
    assert "pull-coord1-k            = 5000" in mdp
    assert m.PROTOCOLS["response10"]["parent_stem"] == "response9"
    text = SBATCH10.read_text(encoding="utf-8")
    assert '-deffnm "$work/response10"' in text
    assert "recapture_response9/response9.gro" in text
    assert "--max-end-distance 0.40" in text


def test_response11_tightens_fixed_target_hold():
    m = load_prep()
    mdp = m.make_mdp("C18", seed=181277, protocol="response11")
    assert "pull-coord1-start        = no" in mdp
    assert "pull-coord1-init         = 0.350000" in mdp
    assert "pull-coord1-k            = 7500" in mdp
    assert m.PROTOCOLS["response11"]["parent_stem"] == "response10"
    text = SBATCH11.read_text(encoding="utf-8")
    assert '-deffnm "$work/response11"' in text
    assert "recapture_response10/response10.gro" in text
    assert "--max-end-distance 0.36" in text
