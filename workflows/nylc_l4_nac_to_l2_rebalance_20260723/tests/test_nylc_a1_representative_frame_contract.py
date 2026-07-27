#!/usr/bin/env python3
import importlib.util
import json
import pathlib
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock

import numpy as np

HERE = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = HERE / "scripts" / "audit_nylc_a1_representative_frame.py"


def load_module():
    spec = importlib.util.spec_from_file_location(
        "audit_nylc_a1_representative_frame", SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class FakeAtom:
    def __init__(
        self,
        index,
        name,
        resname,
        resid,
        position=(0.0, 0.0, 0.0),
        segid="H",
        bonded_names=(),
    ):
        self.index = index
        self.name = name
        self.resname = resname
        self.resid = resid
        self.position = np.asarray(position, dtype=float)
        self.segid = segid
        self.bonded_atoms = SimpleNamespace(names=np.asarray(bonded_names, dtype=object))


class FakeGroup(list):
    @property
    def positions(self):
        return np.asarray([atom.position for atom in self], dtype=float)

    @property
    def names(self):
        return np.asarray([atom.name for atom in self], dtype=object)

    @property
    def resids(self):
        return np.asarray([atom.resid for atom in self], dtype=int)

    def center_of_mass(self):
        return np.mean(self.positions, axis=0)


class A1RepresentativeFrameContractTests(unittest.TestCase):
    def test_frozen_source_and_frame_constants_are_exact(self):
        module = load_module()
        self.assertEqual(module.SELECTED_TIME_PS, 354.000)
        self.assertEqual(module.TIME_TOLERANCE_PS, 0.001)
        self.assertEqual(module.COORDINATE_TOLERANCE_NM, 0.0015)
        self.assertEqual(module.NAC_DISTANCE_MAX_NM, 0.35)
        self.assertEqual(module.NAC_ANGLE_MIN_DEG, 95.0)
        self.assertEqual(module.NAC_ANGLE_MAX_DEG, 115.0)
        self.assertEqual(module.SEVERE_CLASH_CUTOFF_NM, 0.18)
        self.assertEqual(module.EXPECTED_ATOM_COUNT, 133589)
        self.assertEqual(module.EXPECTED_L2_ATOM_COUNT, 79)
        self.assertEqual(module.EXPECTED_L2_HEAVY_COUNT, 33)
        self.assertEqual(
            module.EXPECTED_SOURCE_HASHES,
            {
                "tpr": "c60078a92c2ace51facde4ef64e453f690177fc88b4d6363427f935944fa2e43",
                "xtc": "1a54f1b5b9f139b746c22d9e0f7e9a4a94eb8154bf2b881888986eedca933d89",
            },
        )
        self.assertEqual(
            module.REACTIVE_INDEX1,
            {"thr267_og1": 8960, "l2_c": 10287, "l2_o": 10288, "l2_n": 10289},
        )

    def test_source_hashes_are_pinned_and_each_file_is_hashed_once(self):
        module = load_module()
        fake_hash = mock.Mock(
            side_effect=[
                module.EXPECTED_SOURCE_HASHES["tpr"],
                module.EXPECTED_SOURCE_HASHES["xtc"],
            ]
        )
        observed = module._validate_source_hashes(
            pathlib.Path("run.tpr"),
            pathlib.Path("run.xtc"),
            dict(module.EXPECTED_SOURCE_HASHES),
            hash_file=fake_hash,
        )
        self.assertEqual(observed, module.EXPECTED_SOURCE_HASHES)
        self.assertEqual(fake_hash.call_count, 2)
        self.assertEqual(
            [call.args[0] for call in fake_hash.call_args_list],
            [pathlib.Path("run.tpr"), pathlib.Path("run.xtc")],
        )

    def test_wrong_declared_hash_is_rejected_before_frame_loading(self):
        module = load_module()
        declared = dict(module.EXPECTED_SOURCE_HASHES)
        declared["xtc"] = "0" * 64
        with self.assertRaisesRegex(module.AuditError, "declared source hashes"):
            module._validate_source_hashes(
                pathlib.Path("run.tpr"),
                pathlib.Path("run.xtc"),
                declared,
                hash_file=mock.Mock(),
            )

    def test_time_selection_requires_exactly_one_matching_frame(self):
        module = load_module()
        frames = [
            SimpleNamespace(time=352.0),
            SimpleNamespace(time=354.0004),
            SimpleNamespace(time=356.0),
        ]
        selected = module._select_unique_time(frames, 354.000)
        self.assertIs(selected, frames[1])
        with self.assertRaisesRegex(module.AuditError, "exactly one"):
            module._select_unique_time(
                [SimpleNamespace(time=353.9995), SimpleNamespace(time=354.0005)],
                354.000,
            )
        with self.assertRaisesRegex(module.AuditError, "exactly one"):
            module._select_unique_time([SimpleNamespace(time=353.0)], 354.000)

    def test_requested_time_is_frozen_before_source_access(self):
        module = load_module()
        self.assertEqual(module._require_frozen_time(354.0004), 354.000)
        with mock.patch.object(module, "_validate_source_hashes") as validate_hashes:
            with self.assertRaisesRegex(module.AuditError, "frozen selected time"):
                module.audit_frame(
                    "run.tpr",
                    "run.xtc",
                    "frame.gro",
                    "groups.ndx",
                    355.0,
                    dict(module.EXPECTED_SOURCE_HASHES),
                )
        validate_hashes.assert_not_called()

    def test_coordinate_identity_uses_minimum_image_and_nm_tolerance(self):
        module = load_module()
        box = np.asarray((10.0, 10.0, 10.0, 90.0, 90.0, 90.0))
        source = np.asarray(((0.001, 5.0, 5.0),))
        wrapped = np.asarray(((9.999, 5.0, 5.0),))
        self.assertAlmostEqual(
            module._maximum_pbc_displacement_nm(source, wrapped, box), 0.0002
        )
        with self.assertRaisesRegex(module.AuditError, "coordinate identity"):
            module._validate_coordinate_identity(
                source, np.asarray(((9.90, 5.0, 5.0),)), box
            )

    def test_reactive_atom_identities_are_fail_closed(self):
        module = load_module()
        identities = (
            (FakeAtom(8959, "OG1", "THR", 267), 8960, "THR", "OG1", 267),
            (FakeAtom(10286, "C12", "L2", 1), 10287, "L2", "C12", None),
            (FakeAtom(10287, "O2", "L2", 1), 10288, "L2", "O2", None),
            (FakeAtom(10288, "N3", "L2", 1), 10289, "L2", "N3", None),
        )
        for atom, index1, resname, name, resid in identities:
            module._require_identity(atom, index1, resname, name, resid)
        with self.assertRaisesRegex(module.AuditError, "expected THR267:OG1"):
            module._require_identity(
                FakeAtom(8959, "HG1", "THR", 267), 8960, "THR", "OG1", 267
            )

    def test_a1_bonds_require_n_h1_h2_hg1_and_og1_only_cb(self):
        module = load_module()
        nalpha = FakeAtom(8948, "N", "THR", 267, bonded_names=("CA", "H1", "H2", "HG1"))
        og1 = FakeAtom(8959, "OG1", "THR", 267, bonded_names=("CB",))
        result = module._validate_a1_bonds(nalpha, og1)
        self.assertEqual(result["n_bonded_hydrogens"], ["H1", "H2", "HG1"])
        self.assertEqual(result["og1_bonded_atoms"], ["CB"])
        with self.assertRaisesRegex(module.AuditError, "N bonded hydrogens"):
            module._validate_a1_bonds(
                FakeAtom(8948, "N", "THR", 267, bonded_names=("CA", "H1", "H2")),
                og1,
            )
        with self.assertRaisesRegex(module.AuditError, "OG1 bonded atoms"):
            module._validate_a1_bonds(
                nalpha,
                FakeAtom(8959, "OG1", "THR", 267, bonded_names=("CB", "HG1")),
            )

    def test_joint_nac_uses_pbc_distance_and_angle(self):
        module = load_module()
        box = np.asarray((30.0, 30.0, 30.0, 90.0, 90.0, 90.0))
        result = module._joint_nac(
            oxygen_A=np.asarray((-1.0, 0.0, 0.0)),
            carbon_A=np.asarray((0.0, 0.0, 0.0)),
            og1_A=np.asarray((0.0, 1.0, 0.0)),
            box=box,
        )
        self.assertAlmostEqual(result["c_og1_distance_nm"], 0.1)
        self.assertAlmostEqual(result["o_c_og1_angle_deg"], 90.0)
        self.assertFalse(result["joint_pass"])
        passing = module._joint_nac(
            oxygen_A=np.asarray((-1.0, 0.0, 0.0)),
            carbon_A=np.asarray((0.0, 0.0, 0.0)),
            og1_A=np.asarray((0.7764571353, 2.8977774789, 0.0)),
            box=box,
        )
        self.assertAlmostEqual(passing["c_og1_distance_nm"], 0.3)
        self.assertAlmostEqual(passing["o_c_og1_angle_deg"], 105.0)
        self.assertTrue(passing["joint_pass"])

    def test_gate_is_exactly_261_to_266_and_excludes_thr267(self):
        module = load_module()
        gate = FakeGroup(
            [FakeAtom(i, "CA", "THR", resid) for i, resid in enumerate(range(261, 267))]
        )
        self.assertEqual(module._validate_gate_membership(gate), [261, 262, 263, 264, 265, 266])
        gate.append(FakeAtom(99, "CA", "THR", 267))
        with self.assertRaisesRegex(module.AuditError, "Gate residues"):
            module._validate_gate_membership(gate)

    def test_gate_record_uses_core_and_gate_and_excludes_thr267(self):
        module = load_module()
        core = FakeGroup([FakeAtom(0, "CA", "ALA", 100, (0.0, 0.0, 0.0))])
        gate = FakeGroup(
            [
                FakeAtom(i + 1, "CA", "THR", resid, (1.0, 0.0, 0.0))
                for i, resid in enumerate(range(261, 267))
            ]
        )
        primitive = mock.Mock(return_value=0.25)
        record = module._gate_opening_record(
            core,
            gate,
            np.asarray((10.0, 10.0, 10.0, 90.0, 90.0, 90.0)),
            gate_opening=primitive,
        )
        self.assertEqual(record["core_atom_count"], 1)
        self.assertEqual(record["gate_atom_count"], 6)
        self.assertEqual(record["gate_resids"], [261, 262, 263, 264, 265, 266])
        self.assertTrue(record["thr267_excluded"])
        self.assertAlmostEqual(record["opening_nm"], 0.25)
        np.testing.assert_allclose(primitive.call_args.args[0], (0.1, 0.0, 0.0))

    def test_l2_and_full_system_counts_are_exact(self):
        module = load_module()
        ligand = FakeGroup(
            [FakeAtom(i, "C", "L2", 1) for i in range(33)]
            + [FakeAtom(i + 33, "H1", "L2", 1) for i in range(46)]
        )
        self.assertEqual(
            module._validate_counts(module.EXPECTED_ATOM_COUNT, ligand),
            {"system_atoms": 133589, "l2_atoms": 79, "l2_heavy_atoms": 33},
        )
        with self.assertRaisesRegex(module.AuditError, "system atom count"):
            module._validate_counts(module.EXPECTED_ATOM_COUNT - 1, ligand)

    def test_minimum_contact_reports_exact_partner_identities_under_pbc(self):
        module = load_module()
        ligand = FakeGroup(
            [
                FakeAtom(10286, "C12", "L2", 1, (0.1, 0.0, 0.0), "L"),
                FakeAtom(10287, "O2", "L2", 1, (5.0, 5.0, 5.0), "L"),
            ]
        )
        partners = FakeGroup(
            [
                FakeAtom(20, "CB", "ALA", 12, (4.0, 4.0, 4.0), "H"),
                FakeAtom(21, "NZ", "LYS", 15, (9.9, 0.0, 0.0), "H"),
            ]
        )
        record = module._minimum_contact(
            ligand,
            partners,
            np.asarray((10.0, 10.0, 10.0, 90.0, 90.0, 90.0)),
        )
        self.assertAlmostEqual(record["distance_nm"], 0.02)
        self.assertEqual(record["severe_clash_cutoff_nm"], 0.18)
        self.assertFalse(record["contact_pass"])
        self.assertEqual(record["ligand"]["index1"], 10287)
        self.assertEqual(record["ligand"]["name"], "C12")
        self.assertEqual(record["partner"]["index1"], 22)
        self.assertEqual(record["partner"]["resname"], "LYS")
        self.assertEqual(record["partner"]["resid"], 15)
        self.assertEqual(record["partner"]["name"], "NZ")

    def test_both_contact_classes_must_clear_the_same_severe_clash_cutoff(self):
        module = load_module()
        safe = {
            "distance_nm": 0.18,
            "severe_clash_cutoff_nm": module.SEVERE_CLASH_CUTOFF_NM,
            "contact_pass": True,
        }
        self.assertTrue(
            module._validate_minimum_contacts(
                {
                    "ligand_protein_heavy": dict(safe),
                    "ligand_water_heavy": dict(safe),
                }
            )
        )
        clashing = dict(safe)
        clashing["distance_nm"] = 0.179
        clashing["contact_pass"] = False
        with self.assertRaisesRegex(module.AuditError, "severe-clash cutoff"):
            module._validate_minimum_contacts(
                {
                    "ligand_protein_heavy": dict(safe),
                    "ligand_water_heavy": clashing,
                }
            )

    def test_cli_writes_failure_json_for_audit_error(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary:
            output = pathlib.Path(temporary) / "audit.json"
            failure = module.AuditError(
                "water contact below severe-clash cutoff",
                stage="minimum_contacts",
                gate="minimum_contacts",
            )
            with mock.patch.object(module, "audit_frame", side_effect=failure):
                return_code = module.main(
                    [
                        "--tpr",
                        "run.tpr",
                        "--xtc",
                        "run.xtc",
                        "--gro",
                        "frame.gro",
                        "--ndx",
                        "groups.ndx",
                        "--output",
                        str(output),
                    ]
                )
            self.assertNotEqual(return_code, 0)
            payload = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(
                payload["scientific_status"],
                "FAIL_A1_REPRESENTATIVE_NAC_FRAME",
            )
            self.assertEqual(payload["failed_stage"], "minimum_contacts")
            self.assertEqual(payload["failed_gate"], "minimum_contacts")
            self.assertIn("water contact", payload["error"])

    def test_pass_status_is_emitted_only_when_every_gate_is_true(self):
        module = load_module()
        gates = {
            "source_hashes": True,
            "source_time_unique": True,
            "coordinate_identity": True,
            "atom_mapping": True,
            "a1_bonds": True,
            "joint_nac": True,
            "gate_definition": True,
            "counts_and_box": True,
            "minimum_contacts": True,
        }
        self.assertEqual(
            module._scientific_status(gates), "PASS_A1_REPRESENTATIVE_NAC_FRAME"
        )
        self.assertEqual(
            module.SCIENTIFIC_SCOPE,
            "classical_fixed_topology_preorganization_not_proton_transfer",
        )
        for key in gates:
            failed = dict(gates)
            failed[key] = False
            self.assertEqual(
                module._scientific_status(failed),
                "FAIL_A1_REPRESENTATIVE_NAC_FRAME",
                key,
            )


    def test_runner_contract_is_immutable_and_fail_closed(self):
        runner_path = HERE / "scripts" / "run_nylc_a1_representative_frame.sh"
        slurm_path = HERE / "slurm" / "run_nylc_a1_representative_frame.sbatch"
        self.assertTrue(runner_path.is_file(), runner_path)
        self.assertTrue(slurm_path.is_file(), slurm_path)
        runner = runner_path.read_text(encoding="utf-8")
        self.assertIn(
            "attempt_61970146_4_61970151/npt300free", runner
        )
        self.assertIn(
            "em/attempt_61968026_1_61968026/nac_evt25_time1462ps/input",
            runner,
        )
        self.assertIn(
            "ensemble/candidates/nac_evt25_time1462ps/"
            "build_job_61813799_11/build/source_cycle.ndx",
            runner,
        )
        self.assertIn(
            "c60078a92c2ace51facde4ef64e453f690177fc88b4d6363427f935944fa2e43",
            runner,
        )
        self.assertIn(
            "1a54f1b5b9f139b746c22d9e0f7e9a4a94eb8154bf2b881888986eedca933d89",
            runner,
        )
        self.assertIn('[[ ! -e "$OUT" ]]', runner)
        self.assertIn("source.tmp.gro", runner)
        self.assertIn("printf 'System\\n'", runner)
        self.assertIn("trjconv", runner)
        self.assertIn("-dump 354", runner)
        self.assertIn("audit_nylc_a1_representative_frame.py", runner)
        self.assertIn("flock -x 9", runner)
        self.assertIn("run_history.tsv", runner)
        self.assertIn("run_history.jsonl", runner)
        self.assertIn("STATE=STARTED", runner)
        self.assertIn("NOT_EVALUATED.json", runner)
        self.assertIn("PASS.json", runner)
        self.assertIn('mv "$OUT/source.tmp.gro"', runner)
        self.assertLess(
            runner.index("PASS_A1_REPRESENTATIVE_NAC_FRAME"),
            runner.index('mv "$OUT/source.tmp.gro"'),
        )
        self.assertNotIn("analyze_nylc_m1_proton_geometry.py", runner)

    def test_runner_removes_restrained_define_and_audits_preflight(self):
        runner = (
            HERE / "scripts" / "run_nylc_a1_representative_frame.sh"
        ).read_text(encoding="utf-8")
        self.assertIn("em_cg_flexible_m1.mdp", runner)
        self.assertIn("preflight.unrestrained.mdp", runner)
        self.assertIn('if key == "define":', runner)
        self.assertIn("continue", runner)
        self.assertNotIn(
            'cp "$CODE_ROOT/mdp/em_cg_flexible_m1.mdp"', runner
        )
        self.assertIn("grompp", runner)
        self.assertIn("-maxwarn 0", runner)
        self.assertIn("-po", runner)
        self.assertIn("gmx_mpi", runner)
        self.assertIn("dump -s", runner)
        self.assertIn("#posres_xA", runner)
        self.assertIn("DISRES", runner)
        self.assertIn("nonempty_defines", runner)
        self.assertNotIn("mdrun", runner)

    def test_runner_hashes_inputs_and_promoted_outputs(self):
        runner = (
            HERE / "scripts" / "run_nylc_a1_representative_frame.sh"
        ).read_text(encoding="utf-8")
        for token in (
            '"$SOURCE_TPR"',
            '"$SOURCE_XTC"',
            '"$TOPOLOGY_ROOT/topol.top"',
            '"$OUT/A1_REPRESENTATIVE_FRAME_AUDIT.json"',
            '"$OUT/representative_354ps.gro"',
            '"$OUT/representative_354ps.pdb"',
        ):
            self.assertIn(token, runner)
        self.assertIn('find "$TOPOLOGY_ROOT"', runner)
        self.assertIn("-name '*.itp'", runner)
        self.assertIn("sha256sum", runner)

    def test_slurm_wrapper_is_cpu_only_modest_and_snapshotted(self):
        slurm = (
            HERE / "slurm" / "run_nylc_a1_representative_frame.sbatch"
        ).read_text(encoding="utf-8")
        self.assertIn("#SBATCH -p xahcnormal", slurm)
        self.assertIn("#SBATCH -c 4", slurm)
        self.assertIn("#SBATCH --mem=16G", slurm)
        self.assertIn("#SBATCH -t 00:30:00", slurm)
        self.assertNotIn("#SBATCH --gres", slurm)
        self.assertNotIn("mdrun", slurm)
        self.assertIn("code_snapshots", slurm)
        self.assertIn("SNAPSHOT_SHA256.tsv", slurm)
        self.assertIn("run_nylc_a1_representative_frame.sh", slurm)
        self.assertIn("audit_nylc_a1_representative_frame.py", slurm)
        self.assertIn("em_cg_flexible_m1.mdp", slurm)


if __name__ == "__main__":
    unittest.main()
