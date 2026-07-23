import importlib.util
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
FLOW = HERE.parent
PREP = FLOW / "scripts" / "prepare_true_thr267_recapture_pilot.py"
SBATCH = FLOW / "slurm" / "run_true_thr267_recapture_pilot_array.sbatch"


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
