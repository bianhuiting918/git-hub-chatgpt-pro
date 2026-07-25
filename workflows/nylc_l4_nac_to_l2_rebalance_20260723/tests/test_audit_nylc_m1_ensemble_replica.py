import json
import pathlib
import sys

import pytest

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from audit_nylc_m1_ensemble_replica import audit_replica


def write_xvg(path, rows):
    path.write_text(
        "\n".join(" ".join(str(value) for value in row) for row in rows) + "\n"
    )


def make_replica(root, *, lincs=False, unbound=False):
    root.mkdir()
    times = [0.0, 20.0, 40.0, 60.0, 80.0, 100.0]
    write_xvg(root / "nac_distance.xvg", zip(times, [0.34, 0.36, 0.33, 0.32, 0.37, 0.34]))
    write_xvg(root / "nac_angle.xvg", zip(times, [105.0, 104.0, 108.0, 110.0, 100.0, 112.0]))
    write_xvg(root / "gate_opening.xvg", zip(times, [1.0, 1.1, 1.2, 1.1, 1.0, 1.2]))
    write_xvg(
        root / "thermo.xvg",
        zip(
            times,
            [299.0, 300.0, 301.0, 300.0, 299.0, 301.0],
            [1.0, -5.0, 8.0, -3.0, 2.0, 4.0],
            [1000.0, 1002.0, 1001.0, 1003.0, 1002.0, 1001.0],
        ),
    )
    contacts = [4, 4, 4, 4, 4, 0 if unbound else 4]
    com = [0.8, 0.8, 0.9, 0.8, 0.9, 1.5 if unbound else 0.8]
    write_xvg(root / "pocket_state.xvg", zip(times, contacts, com))
    (root / "minimum_contact.json").write_text(
        json.dumps(
            {
                "minimum_ligand_protein_heavy_nm": 0.22,
                "minimum_ligand_water_heavy_nm": 0.20,
            }
        )
        + "\n"
    )
    with (root / "proton_network.jsonl").open("w") as handle:
        for time in times:
            handle.write(
                json.dumps(
                    {
                        "time_ps": time,
                        "channels": {
                            "OgH_to_Nalpha_preorg": {"favorable": time in (40.0, 60.0)},
                            "OgH_to_Asp306": {"favorable": False},
                            "OgH_to_Asp308": {"favorable": False},
                            "OgH_to_water": {"favorable": time >= 20.0},
                            "OgH_via_water_to_Asp306": {"favorable": time == 60.0},
                            "OgH_via_water_to_Asp308": {"favorable": False},
                        },
                    }
                )
                + "\n"
            )
    log = "Finished mdrun\n"
    if lincs:
        log += "LINCS WARNING\n"
    (root / "run.log").write_text(log)


def manifest():
    return {
        "candidate_id": "candidate_01",
        "velocity_seed": 26711,
        "microstate": "M1_NalphaH2_Asp306H",
        "gate_definition": "NylC residues 261-266; Thr267 excluded",
        "fully_unrestrained": True,
    }


def test_pass_replica_reports_denominator_nac_gate_thermo_and_proton_channels(tmp_path):
    root = tmp_path / "pass"
    make_replica(root)

    audit = audit_replica(root, manifest(), (0.0, 100.0))

    assert audit["technical_status"] == "PASS"
    assert audit["scientific_status"] == "PASS_REPLICA_NAC_PRESENT"
    assert audit["analysis"]["frame_count"] == 6
    assert audit["analysis"]["sampling_interval_ps"] == pytest.approx(20.0)
    assert audit["nac"]["nac_frame_count"] == 4
    assert audit["nac"]["nac_frame_count_after_20ps"] == 3
    assert audit["nac"]["longest_event"]["duration_ps"] == pytest.approx(20.0)
    assert audit["gate_opening_nm"]["frame_count"] == 6
    assert audit["thermodynamics"]["temperature_K"]["mean"] == pytest.approx(300.0)
    assert audit["bound_state"]["final_frame_retained"] is True
    assert audit["proton_preorganization"]["OgH_to_Nalpha_preorg"]["nac_favorable_frames"] == 2


def test_lincs_makes_science_not_evaluated(tmp_path):
    root = tmp_path / "lincs"
    make_replica(root, lincs=True)

    audit = audit_replica(root, manifest(), (0.0, 100.0))

    assert audit["technical_status"] == "FAIL"
    assert audit["scientific_status"] == "NOT_EVALUATED_NUMERICAL_FAIL"
    assert audit["numerical_issue_counts"]["lincs_warning"] == 1


def test_final_pocket_departure_is_fail_unbound(tmp_path):
    root = tmp_path / "unbound"
    make_replica(root, unbound=True)

    audit = audit_replica(root, manifest(), (0.0, 100.0))

    assert audit["technical_status"] == "PASS"
    assert audit["scientific_status"] == "FAIL_UNBOUND"
    assert audit["bound_state"]["final_frame_retained"] is False
