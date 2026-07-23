import json
from pathlib import Path

FLOW = Path(__file__).resolve().parents[1]


def test_nylc_wrong_atom_candidates_are_superseded_in_manifest():
    manifest = json.loads((FLOW / "manifests/candidates.json").read_text())
    nylc = [x for x in manifest["candidates"] if x["enzyme"] == "NylC-GYAQ"]
    assert len(nylc) == 2
    for candidate in nylc:
        assert candidate["status"] == "SUPERSEDED_WRONG_NUCLEOPHILE_IDENTITY"
        assert candidate["identity_correction"]["wrong_thr262_og1"] == 8896
        assert candidate["identity_correction"]["true_thr267_og1"] == 8961


def test_prepare_entrypoints_refuse_superseded_candidates_before_writes():
    text = (FLOW / "scripts/prepare_candidates.py").read_text()
    assert 'candidate.get("status", "ACTIVE")' in text
    assert "Refusing superseded candidate" in text
    assert "NOT_EVALUATED_SUPERSEDED" in text
