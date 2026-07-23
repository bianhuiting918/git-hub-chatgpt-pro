from pathlib import Path

FLOW = Path(__file__).resolve().parents[1]


def _read(relative: str) -> str:
    return (FLOW / relative).read_text()


def test_free_1ns_audit_uses_exact_nac_atoms_and_gate_group():
    text = _read("scripts/audit_nylc_c18_free_1ns.sh")
    assert "atomnr 8896 plus atomnr 10287" in text
    assert "8896 10287 10288" in text
    assert 'group "Core"' in text
    assert 'group "Gate"' in text
    assert "source_cycle.ndx" in text
    assert "261-266" in text
    assert "Thr267" in text


def test_free_1ns_audit_uses_common_gates_and_full_thermo():
    text = _read("scripts/audit_nylc_c18_free_1ns.sh")
    assert "analyze_nac_series.py" in text
    assert "--distance-max-nm 0.35" in text
    assert "--angle-min-deg 95" in text
    assert "--angle-max-deg 115" in text
    for term in ("Potential", "Temperature", "Pressure", "Volume"):
        assert term in text


def test_free_1ns_audit_does_not_auto_approve_qmmm():
    text = _read("scripts/audit_nylc_c18_free_1ns.sh")
    assert "NOT_EVALUATED_PENDING_SCIENTIFIC_INTERPRETATION" in text
    assert "FINAL_SCIENTIFIC_PASS" not in text
    assert "PASS_QMMM" not in text
