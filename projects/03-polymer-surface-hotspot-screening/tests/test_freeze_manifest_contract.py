from pathlib import Path


SBATCH = (
    Path(__file__).resolve().parents[1]
    / "slurm"
    / "freeze_structure_manifests.sbatch"
)


def test_freeze_job_is_commit_pinned_and_authority_aware():
    text = SBATCH.read_text(encoding="utf-8")
    assert "DEPLOY_COMMIT=${DEPLOY_COMMIT:?" in text
    assert "OUTPUT_DIR=${OUTPUT_DIR:?" in text
    assert "NYLON_EXACT_MANIFEST=${NYLON_EXACT_MANIFEST:?" in text
    assert "--nylon-exact-manifest" in text
    assert "--nylon-authority-manifest" in text
    assert "preflight_status.tsv" in text
    assert 'if [ -e "$OUTPUT_DIR" ]' in text
    assert 'sha256sum "$OUTPUT_DIR"/*.tsv "$OUTPUT_DIR"/*.json' in text
