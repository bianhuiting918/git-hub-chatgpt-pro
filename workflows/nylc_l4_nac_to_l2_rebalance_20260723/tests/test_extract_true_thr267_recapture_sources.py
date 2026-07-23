import importlib.util
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
SCRIPT = HERE.parent / "scripts" / "extract_true_thr267_recapture_sources.py"


def load_module():
    spec = importlib.util.spec_from_file_location("extract_true_sources", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_corrected_source_contract():
    m = load_module()
    assert m.TRUE_THR267_OG1 == 8961
    assert m.SUPERSEDED_THR262_OG1 == 8896
    assert set(m.CANDIDATES) == {"nylc_C18_trueT267_recapture", "nylc_C23_trueT267_recapture"}
    assert m.CANDIDATES["nylc_C18_trueT267_recapture"]["local_time_ps"] == 16.0
    assert m.CANDIDATES["nylc_C23_trueT267_recapture"]["local_time_ps"] == 98.0
    assert m.CANDIDATES["nylc_C18_trueT267_recapture"]["carbonyl_c"] == 10297
    assert m.CANDIDATES["nylc_C18_trueT267_recapture"]["carbonyl_c_atomname"] == "C18"
    assert m.CANDIDATES["nylc_C23_trueT267_recapture"]["carbonyl_c"] == 10303
    assert m.CANDIDATES["nylc_C23_trueT267_recapture"]["carbonyl_c_atomname"] == "C23"


def test_temporary_name_is_mandatory():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "source.tmp.gro" in text
    assert "os.replace" in text
    assert "SUPERSEDED_WRONG_NUCLEOPHILE_IDENTITY" in text
    assert "source.gro" in text


def test_identity_parser_uses_sequential_global_index(tmp_path):
    m = load_module()
    gro = tmp_path / "tiny.gro"
    gro.write_text(
        "tiny t= 16.000\n"
        "2\n"
        "  262THR    OG1 8896   0.000   0.000   0.000\n"
        "  267THR    OG1 8961   0.100   0.200   0.300\n"
        "   1.0 1.0 1.0\n",
        encoding="utf-8",
    )
    atoms, box, title = m.read_gro(gro)
    assert atoms[1]["resid"] == 262
    assert atoms[2]["resid"] == 267
    assert atoms[2]["atomname"] == "OG1"
    assert box == (1.0, 1.0, 1.0)
    assert "16.000" in title
