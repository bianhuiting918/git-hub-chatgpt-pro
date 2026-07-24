from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREP = ROOT / "scripts" / "prepare_nylc_true_thr267_dftb3_preflight.py"
AUDIT = ROOT / "scripts" / "audit_nylc_dftb3_smoke.py"
SBATCH = ROOT / "slurm" / "run_nylc_true_thr267_dftb3_preflight_after_rebalance.sbatch"


def test_preparer_requires_free_l2_nac_and_selects_lowest_potential_frame():
    text = PREP.read_text()
    assert "PASS_UNRESTRAINED_L2_NAC_PRESENT" in text
    assert "lowest_potential_nac_frame" in text
    assert "source.tmp.gro" in text
    assert "selected_free_l2_nac.gro" in text
    assert "NOT_EVALUATED_NO_FREE_L2_NAC" in text


def test_preparer_audits_true_thr267_and_complete_ligand_qm_region():
    text = PREP.read_text()
    assert "THR267_OG1 = 8961" in text
    assert "PROTEIN_ATOMS = 10272" in text
    assert "L2_FIRST = 10273" in text
    assert "L2_LAST = 10351" in text
    assert "L2_REACTIVE_C = 10287" in text
    assert "L2_REACTIVE_O = 10288" in text
    assert "expected one QM/MM boundary bond" in text
    assert "electron_parity" in text


def test_preparer_uses_declared_dftb3_charge_spin_and_short_smoke():
    text = PREP.read_text()
    assert "qm_theory='DFTB3'" in text
    assert "QMCHARGE = 1" in text
    assert "LINK_ATOMIC_NUMBER = 1" in text
    assert "qmcharge={QMCHARGE}" in text
    assert " + LINK_ATOMIC_NUMBER - QMCHARGE" in text
    assert "spin=1" in text
    assert "maxcyc=1" in text
    assert "maxcyc=20" in text


def test_smoke_audit_is_numerical_only():
    text = AUDIT.read_text()
    assert "PASS_DFTB3_PREFLIGHT" in text
    assert "numerical preflight only" in text
    for token in ("SANDER BOMB", "Convergence could not be achieved", "vlimit", "NaN", "FATAL"):
        assert token in text


def test_runner_waits_for_postprocess_gate_and_uses_scnet_cpu_amber18():
    text = SBATCH.read_text()
    assert "#SBATCH -p xahcnormal" in text
    assert "--gres=dcu" not in text
    assert 'POSTJOB=${POSTJOB:?set POSTJOB to the successful postprocess job id}' in text
    assert "module load amber/2018" in text
    assert "export GMXDATA=/public/software/apps/Gromacs-DCU2/2022.1/mpi/share/gromacs" in text
    assert "3ob-3-1/C-C.skf" in text
    assert "01_qmmm_one_step.in" in text
    assert "02_qmmm_20_step.in" in text
    assert "NOT_EVALUATED.json" in text
    assert 'step1_dftb3_preflight_post${POSTJOB}_job${job}' in text
    assert "TS search" in text


def test_runner_converts_failed_postprocess_dependency_to_not_evaluated():
    text = SBATCH.read_text()
    assert "PASS_POSTPROCESS_TECHNICAL" in text
    assert "NOT_EVALUATED_POSTPROCESS_GATE" in text


def test_minimal_smoke_qm_region_is_not_reused_as_production_ts_region():
    text = PREP.read_text()
    assert "minimal_numerical_smoke" in text
    for residue in ("Tyr146", "Lys189", "Asn219", "Asp306", "Asp308"):
        assert residue in text
