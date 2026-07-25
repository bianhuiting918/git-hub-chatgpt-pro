import pathlib
import sys
from types import SimpleNamespace

import pytest

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from generate_nylc_m1_ensemble_primitives import (
    _local_resname_matches_m1,
    gate_opening_nm,
    parse_ndx,
    pocket_retained,
)


def test_parse_ndx_keeps_one_based_membership_order(tmp_path):
    path = tmp_path / "groups.ndx"
    path.write_text("[ Core ]\n1 2 3\n[ Gate ]\n10 11\n12\n")

    groups = parse_ndx(path)

    assert groups["Core"] == [1, 2, 3]
    assert groups["Gate"] == [10, 11, 12]


def test_gate_opening_uses_frozen_axis_and_baseline():
    assert gate_opening_nm((0.0, 0.0, 0.0)) == pytest.approx(1.51109)
    assert gate_opening_nm((0.4904295935403325, 0.813080787402625, -0.3136533866174433)) == pytest.approx(2.51109)


@pytest.mark.parametrize(
    "contact_count,com_nm,expected",
    [(3, 1.2, True), (2, 1.2, False), (3, 1.200001, False)],
)
def test_pocket_retained_contract(contact_count, com_nm, expected):
    assert pocket_retained(contact_count, com_nm) is expected


@pytest.mark.parametrize(
    "record,residue,expected",
    [
        ({"resid": 306, "residue_instance": 661, "resname": "ASP"}, SimpleNamespace(resid=999, resname="ASP"), True),
        ({"resid": 306, "residue_instance": 661, "resname": "ASP"}, SimpleNamespace(resid=999, resname="ASH"), True),
        ({"resid": 308, "residue_instance": 663, "resname": "ASP"}, SimpleNamespace(resid=999, resname="ASH"), False),
        ({"resid": 306, "residue_instance": 662, "resname": "ASP"}, SimpleNamespace(resid=999, resname="ASH"), False),
        ({"resid": 306, "residue_instance": 661, "resname": "GLU"}, SimpleNamespace(resid=999, resname="GLH"), False),
    ],
)
def test_local_resname_equivalence_is_limited_to_m1_ash306(record, residue, expected):
    assert _local_resname_matches_m1(record, residue) is expected
