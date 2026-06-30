#!/usr/bin/env python3
"""Generate blind PETase Stage 3/4 mechanism and QM/MM exploration templates."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


PATHS = [
    {
        "path_id": "AC_DIRECT_HIS_GENERAL_BASE",
        "reaction_stage": "acylation",
        "description": "Ser O-gamma attacks the ester carbonyl while catalytic His acts as the general base.",
        "proton_relay": "Ser_to_His_direct",
        "intermediate_model": "tetrahedral_or_concerted_to_be_resolved",
        "initial_gate": "accepted_acylation_michaelis_pose",
        "reject_if": "Ser attack geometry, His accessibility, or leaving group path fails Stage 1 filters",
    },
    {
        "path_id": "AC_WATER_ASSISTED_LEAVING_GROUP",
        "reaction_stage": "acylation",
        "description": "Ser activation is His-mediated; leaving-group proton delivery may be water-assisted.",
        "proton_relay": "Ser_to_His_and_His_or_water_to_leaving_O",
        "intermediate_model": "tetrahedral_or_stepwise",
        "initial_gate": "accepted_acylation_pose_with_relay_water",
        "reject_if": "No persistent relay water or leaving oxygen access after restrained relaxation",
    },
    {
        "path_id": "AC_ASP_COUPLED_HIS_BRANCH",
        "reaction_stage": "acylation",
        "description": "Catalytic Asp/His protonation sensitivity branch tested only if pKa/protonation gate supports it.",
        "proton_relay": "His_Asp_network_sensitivity",
        "intermediate_model": "same_as_parent_path",
        "initial_gate": "protonation_gate_disagrees_or_flags_sensitivity",
        "reject_if": "pKa/protonation evidence does not support an alternate Asp/His state",
    },
    {
        "path_id": "DE_WATER_ATTACK_HIS_GENERAL_BASE",
        "reaction_stage": "deacylation",
        "description": "A bound water attacks the acyl-enzyme carbonyl while His acts as the general base.",
        "proton_relay": "water_to_His_then_SerO_protonation",
        "intermediate_model": "tetrahedral_or_concerted_to_be_resolved",
        "initial_gate": "accepted_acyl_enzyme_pose_with_nucleophilic_water",
        "reject_if": "No water has productive attack geometry and His access",
    },
    {
        "path_id": "DE_WATER_NETWORK_ASSISTED",
        "reaction_stage": "deacylation",
        "description": "Water attack is coupled to a local water network before product release.",
        "proton_relay": "water_network_to_His_or_leaving_SerO",
        "intermediate_model": "stepwise_water_network_branch",
        "initial_gate": "stable_water_network_in_classical_ensemble",
        "reject_if": "Water network is not persistent across equilibrated conformers",
    },
    {
        "path_id": "DE_CARBOXYLATE_PROTONATION_BRANCH",
        "reaction_stage": "deacylation",
        "description": "Distal fragment carboxylate protonation sensitivity branch if protonation gate supports it.",
        "proton_relay": "fragment_carboxylate_state_sensitivity",
        "intermediate_model": "same_as_parent_path",
        "initial_gate": "pH_or_local_environment_supports_neutral_fragment",
        "reject_if": "No protonation evidence for neutral distal carboxylate",
    },
]


CVS = [
    ("acylation", "AC_DIRECT_HIS_GENERAL_BASE", "acyl_ser_o_to_carbonyl_c", "bond_formation", "distance", "SerOG Csc_acyl", "generic_serine_hydrolase_chemistry"),
    ("acylation", "AC_DIRECT_HIS_GENERAL_BASE", "ser_h_to_his_n", "proton_transfer", "distance_difference", "SerHG-HisN and SerOG-SerHG", "generic_serine_hydrolase_chemistry"),
    ("acylation", "AC_DIRECT_HIS_GENERAL_BASE", "carbonyl_c_to_leaving_o", "bond_breaking", "distance", "Csc_acyl Olg_ethylene_glycol", "generic_serine_hydrolase_chemistry"),
    ("acylation", "AC_WATER_ASSISTED_LEAVING_GROUP", "his_or_water_to_leaving_o", "proton_transfer", "distance_difference", "HisH/WatH to leaving O", "generic_serine_hydrolase_chemistry"),
    ("acylation", "AC_ASP_COUPLED_HIS_BRANCH", "his_to_asp_proton_network", "protonation_sensitivity", "distance_or_protonation_state", "HisN AspOD", "generic_serine_hydrolase_chemistry"),
    ("deacylation", "DE_WATER_ATTACK_HIS_GENERAL_BASE", "water_o_to_acyl_c", "bond_formation", "distance", "Wat_nuc_O Cacyl_ser_ester", "generic_serine_hydrolase_chemistry"),
    ("deacylation", "DE_WATER_ATTACK_HIS_GENERAL_BASE", "water_h_to_his_n", "proton_transfer", "distance_difference", "WatH-HisN and WatO-WatH", "generic_serine_hydrolase_chemistry"),
    ("deacylation", "DE_WATER_ATTACK_HIS_GENERAL_BASE", "acyl_c_to_ser_o", "bond_breaking", "distance", "Cacyl_ser_ester Ser160_OG", "generic_serine_hydrolase_chemistry"),
    ("deacylation", "DE_WATER_NETWORK_ASSISTED", "water_network_relay", "proton_transfer", "distance_network", "Wat_nuc Wat_relay HisN", "generic_serine_hydrolase_chemistry"),
    ("deacylation", "DE_CARBOXYLATE_PROTONATION_BRANCH", "fragment_carboxylate_state", "protonation_sensitivity", "state_label", "distal carboxylate", "generic_serine_hydrolase_chemistry"),
]


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_mechanism_yaml(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "blind_stage3_boundary:",
        "  allowed_inputs:",
        "    - accepted Stage 1 Michaelis/acyl-enzyme poses",
        "    - substrate chemistry",
        "    - generic serine-hydrolase chemistry",
        "    - protonation-gate evidence",
        "  forbidden_inputs:",
        "    - paper TS coordinates",
        "    - paper reaction-coordinate formulas",
        "    - paper selected CVs",
        "    - paper trajectories or umbrella windows",
        "    - paper barrier or rate-limiting conclusions",
        "mechanism_hypotheses:",
    ]
    for path_row in PATHS:
        lines.extend(
            [
                f"  - path_id: {path_row['path_id']}",
                f"    reaction_stage: {path_row['reaction_stage']}",
                f"    description: \"{path_row['description']}\"",
                f"    proton_relay: {path_row['proton_relay']}",
                f"    intermediate_model: {path_row['intermediate_model']}",
                f"    initial_gate: {path_row['initial_gate']}",
                f"    reject_if: \"{path_row['reject_if']}\"",
            ]
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate(out_root: Path) -> list[Path]:
    stage3 = out_root / "03_mechanism_tree"
    stage4 = out_root / "04_qmmm_exploration"
    generated = [
        stage3 / "mechanism_hypotheses.yaml",
        stage3 / "candidate_cv_sets.tsv",
        stage4 / "path_screening_table.tsv",
        stage4 / "ts_like_guess_manifest.tsv",
    ]

    write_mechanism_yaml(generated[0])
    write_tsv(
        generated[1],
        ["reaction_stage", "path_id", "cv_id", "cv_role", "cv_type", "atom_or_state_definition", "source"],
        [
            {
                "reaction_stage": stage,
                "path_id": path_id,
                "cv_id": cv_id,
                "cv_role": role,
                "cv_type": cv_type,
                "atom_or_state_definition": definition,
                "source": source,
            }
            for stage, path_id, cv_id, role, cv_type, definition, source in CVS
        ],
    )
    write_tsv(
        generated[2],
        [
            "path_id",
            "reaction_stage",
            "starting_pose_manifest",
            "qm_region_rule",
            "low_cost_method",
            "scan_or_string_status",
            "status",
            "notes",
        ],
        [
            {
                "path_id": row["path_id"],
                "reaction_stage": row["reaction_stage"],
                "starting_pose_manifest": "blind_work/01_system_setup/gs_pose_manifest.tsv",
                "qm_region_rule": "include_reacting_ester_or_acyl_group_Ser_His_Asp_if_supported_and_relay_water",
                "low_cost_method": "DFTB3_MM_or_xTB_MM_or_low_cost_DFT_MM",
                "scan_or_string_status": "not_started",
                "status": "not_started",
                "notes": row["initial_gate"],
            }
            for row in PATHS
        ],
    )
    write_tsv(
        generated[3],
        [
            "ts_guess_id",
            "path_id",
            "source_pose_id",
            "source_method",
            "structure_path",
            "imaginary_mode_status",
            "committor_status",
            "validation_status",
            "notes",
        ],
        [
            {
                "ts_guess_id": "pending",
                "path_id": "pending",
                "source_pose_id": "pending",
                "source_method": "scan_maximum_or_string_image_after_stage4_runs",
                "structure_path": "pending",
                "imaginary_mode_status": "not_generated",
                "committor_status": "not_generated",
                "validation_status": "not_generated",
                "notes": "Do not populate from paper TS coordinates or paper trajectories.",
            }
        ],
    )
    return generated


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-root", default="work/blind_work")
    args = parser.parse_args()
    generated = generate(Path(args.out_root))
    for path in generated:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
