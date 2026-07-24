from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
POST = ROOT / "scripts" / "postprocess_nylc_true_thr267_rebalance.py"
SBATCH = ROOT / "slurm" / "run_postprocess_nylc_true_thr267_rebalance.sbatch"
FREE_AUDIT = ROOT / "scripts" / "audit_nylc_true_thr267_l2_free.py"
RECAPTURE_AUDIT = ROOT / "scripts" / "audit_true_thr267_recapture_pilot.py"


def test_postprocessor_uses_existing_final_md_without_rerunning_mdrun():
    text = POST.read_text()
    assert "npt300free" in text
    assert "run.xtc" in text and "run.tpr" in text and "run.edr" in text
    assert "mdrun" not in text
    assert "postprocess_job_" in text


def test_postprocessor_uses_correct_true_thr267_geometry_and_gate():
    text = POST.read_text()
    assert "atomnr 8961 plus atomnr 10287" in text
    assert "8961 10287 10288" in text
    assert 'group "Core" plus com of group "Gate"' in text
    assert "residues 261-266" in text
    assert "audit_nylc_true_thr267_l2_free.py" in text


def test_postprocessor_repairs_only_text_history_with_backup():
    text = POST.read_text()
    assert "pre_repair_61708900" in text
    assert "run_history.tsv" in text
    assert "61708900" in text
    assert "run.gro" not in text.split("repair_history", 1)[1].split("def ", 1)[0]


def test_postprocessor_and_wrapper_have_no_double_escaped_sequences():
    for path in (POST, SBATCH):
        text = path.read_text()
        assert r"\\s" not in text
        assert r"\\n" not in text
        assert "\\\\t" not in text


def test_wrapper_is_cpu_only_and_records_dependency_target():
    text = SBATCH.read_text()
    assert "--gres=dcu" not in text
    assert "61708900" in text
    assert "postprocess_nylc_true_thr267_rebalance.py" in text


def test_wrapper_uses_cpu_partition_available_without_dcu_gres():
    text = SBATCH.read_text()
    assert "#SBATCH -p xahcnormal" in text


def test_free_l2_audit_identifies_actual_l2_residue_without_breaking_l4_default():
    free_text = FREE_AUDIT.read_text()
    recapture_text = RECAPTURE_AUDIT.read_text()
    assert 'ligand_resname="L2"' in free_text
    assert 'def minimum_ligand_protein_distance(atoms, box, ligand_resname="UNL")' in recapture_text
    assert 'atom["resname"] == ligand_resname' in recapture_text
