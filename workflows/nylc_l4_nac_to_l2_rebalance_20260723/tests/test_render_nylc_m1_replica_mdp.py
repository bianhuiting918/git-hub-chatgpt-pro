import sys
from pathlib import Path

import pytest

FLOW = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(FLOW / "scripts"))

from render_nylc_m1_replica_mdp import parse_mdp, render_mdp


RESTRAINED_TEMPLATE = """define = -DPOSRES -DPOSRES_L2_1000
integrator = md
gen-vel = yes
gen-seed = 20260723
"""


def test_renders_only_existing_key_with_approved_velocity_seed():
    rendered = render_mdp(RESTRAINED_TEMPLATE, {"gen-seed": "26711"})
    values = parse_mdp(rendered)

    assert values["gen-seed"] == "26711"
    assert values["define"] == "-DPOSRES -DPOSRES_L2_1000"


@pytest.mark.parametrize("seed", ["0", "-1", "26712", "random"])
def test_rejects_nonapproved_velocity_seed(seed):
    with pytest.raises(ValueError, match="approved"):
        render_mdp(RESTRAINED_TEMPLATE, {"gen-seed": seed})


def test_rejects_missing_or_duplicate_keys():
    with pytest.raises(ValueError, match="not present"):
        render_mdp(RESTRAINED_TEMPLATE, {"ref-t": "50"})
    with pytest.raises(ValueError, match="duplicate"):
        render_mdp("gen-seed = 1\ngen-seed = 2\n", {"gen-seed": "26711"})


def test_rejects_restraint_injection_into_free_mdp():
    free_template = "integrator = md\ngen-vel = no\ncontinuation = yes\n"
    with pytest.raises(ValueError, match="define"):
        render_mdp(free_template, {"define": "-DPOSRES"})
