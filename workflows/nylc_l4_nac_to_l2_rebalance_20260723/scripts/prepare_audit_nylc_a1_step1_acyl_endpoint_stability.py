#!/usr/bin/env python3
"""Construct, release, and audit dual-seed NylC A1 Step1 acyl endpoints."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import pathlib
import re
import warnings
from typing import Any, Iterable, Mapping, Sequence

BASE_PATH = pathlib.Path(__file__).with_name("prepare_audit_nylc_a1_step1_pt2_cn_scout.py")
_SPEC = importlib.util.spec_from_file_location("_nylc_a1_pt2_cn_authority", BASE_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"cannot load frozen authority module {BASE_PATH}")
BASE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(BASE)

SOURCES = BASE.SOURCES
PRMTOP = BASE.PRMTOP
EXPECTED_SYSTEM_ATOMS = BASE.EXPECTED_SYSTEM_ATOMS
NO_QM_WATER = True
REACTIVE_ATOMS = {
    "nalpha": 8949,
    "og1": 8960,
    "hg1": 8961,
    "c11": 10286,
    "c12": 10287,
    "o2": 10288,
    "n3": 10289,
}
INTERMEDIATE_TARGETS_A = {"attack": 1.60, "hg1_n3": 1.40, "cn": 1.75}
PRODUCT_TARGETS_A = {"attack": 1.50, "hg1_n3": 1.05, "cn": 2.20}
BUILD_FORCE_KCAL_MOL_A2 = (25.0, 50.0)
LOCAL_RELEASE_STEPS = 800
FULL_RELEASE_STEPS = 800
RELEASE_MD_STEPS = 500
RELEASE_DT_PS = 0.0005
RELEASE_NTWX = 10
PRODUCT_GATE = {
    "attack_A": (1.40, 1.65),
    "hg1_n3_A": (0.95, 1.20),
    "nalpha_hg1_A_min": 1.55,
    "qPT_A_min": 0.45,
    "c12_n3_A_min": 2.05,
    "c12_o2_A": (1.18, 1.30),
    "product_out_of_plane_A_max": 0.12,
    "product_angle_sum_deg_min": 350.0,
    "last_fraction": 0.40,
    "occupancy_min": 0.80,
}
SEED_CLASSIFICATIONS = (
    "PERSISTS_ACYL_PRODUCT",
    "RETURNS_TETRAHEDRAL",
    "RETURNS_REACTANT",
    "ZWITTERIONIC_CLEAVAGE",
    "MISROUTED_PROTON",
    "RESTRAINT_DEPENDENT_PRODUCT",
    "UNCLASSIFIED_RELEASE_ENDPOINT",
    "NOT_EVALUATED_TECHNICAL_FAILURE",
)
AUTO_START_PATH_SAMPLING = False
TERMINAL_BOUNDARIES = (
    "NOT_EVALUATED_TS_PMF_BARRIER_MECHANISM",
    "DO_NOT_START_PMF",
)
TECHNICAL_PASS = "PASS_TECHNICAL_A1_ACYL_ENDPOINT_STABILITY"
TECHNICAL_FAIL = "NOT_EVALUATED_A1_ACYL_ENDPOINT_STABILITY"
PERSISTABLE = {
    "ENDPOINT_MANIFEST.json",
    "RESULT.json",
    "PASS.json",
    "NOT_EVALUATED.json",
    "SHA256.tsv",
    "run_history.tsv",
    "run_history.jsonl",
    "constructed_product.rst7",
    "released_endpoint.rst7",
}
STAGE_ORDER = ("intermediate", "product", "local_release", "full_release", "release_md")


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: pathlib.Path, payload: Mapping[str, Any]) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def validate_source_authority(source: Mapping[str, Any]) -> Any:
    """Delegate frozen-source and QM authority validation without reconstruction."""
    return BASE.validate_authority(source)


def source_from_index(index: int) -> Mapping[str, Any]:
    if int(index) not in (0, 1):
        raise ValueError("seed index must be 0 or 1")
    return SOURCES[int(index)]


def validate_n3_reactant_graph(neighbors: Mapping[int, str]) -> bool:
    normalized = {int(index): str(element).upper() for index, element in neighbors.items()}
    values = list(normalized.values())
    return (
        len(normalized) == 3
        and normalized.get(REACTIVE_ATOMS["c12"]) == "C"
        and values.count("C") == 2
        and values.count("H") == 1
    )


def require_fresh_candidate_dir(path: pathlib.Path) -> None:
    if path.exists():
        raise FileExistsError(f"candidate attempt already exists: {path}")


def is_persistable_artifact(name: str) -> bool:
    return pathlib.Path(name).name in PERSISTABLE


def stage_spec(stage: str, seed: str) -> dict[str, Any]:
    if seed not in {"seed26723", "seed26737"}:
        raise ValueError(f"unknown seed {seed}")
    specs = {
        "intermediate": {
            "maxcyc": 400, "ncyc": 100, "targets": INTERMEDIATE_TARGETS_A,
            "force": BUILD_FORCE_KCAL_MOL_A2[0], "reactive_restraints": ("attack", "hg1_n3", "cn"),
            "environment_restraint": True, "ntr": 1, "nmropt": 1, "disang": "restraints.RST",
        },
        "product": {
            "maxcyc": 800, "ncyc": 200, "targets": PRODUCT_TARGETS_A,
            "force": BUILD_FORCE_KCAL_MOL_A2[1], "reactive_restraints": ("attack", "hg1_n3", "cn"),
            "environment_restraint": True, "ntr": 1, "nmropt": 1, "disang": "restraints.RST",
        },
        "local_release": {
            "maxcyc": LOCAL_RELEASE_STEPS, "ncyc": 200, "targets": None, "force": None,
            "reactive_restraints": (), "environment_restraint": True,
            "ntr": 1, "nmropt": 0, "disang": None,
        },
        "full_release": {
            "maxcyc": FULL_RELEASE_STEPS, "ncyc": 200, "targets": None, "force": None,
            "reactive_restraints": (), "environment_restraint": False,
            "ntr": 0, "nmropt": 0, "disang": None,
        },
    }
    if stage not in specs:
        raise ValueError(f"not a minimization stage: {stage}")
    return dict(specs[stage])


def release_md_spec(seed: str, full_release: Mapping[str, Any]) -> dict[str, Any]:
    if seed not in {"seed26723", "seed26737"}:
        raise ValueError(f"unknown seed {seed}")
    restart = str(full_release.get("restart_path", ""))
    restart_sha = str(full_release.get("restart_sha256", ""))
    if (
        full_release.get("stage") != "full_release"
        or full_release.get("technical_pass") is not True
        or not restart
        or not re.fullmatch(r"[0-9a-f]{64}", restart_sha)
    ):
        raise ValueError("release_md requires a technical-PASS SHA-fixed full_release restart")
    return {
        "stage": "release_md", "input_restart_sha256": restart_sha,
        "ig": 26723 if seed == "seed26723" else 26737,
        "nstlim": RELEASE_MD_STEPS, "dt": RELEASE_DT_PS, "temp0": 300.0,
        "ntwx": RELEASE_NTWX, "ntt": 3, "gamma_ln": 2.0, "ntc": 1, "ntf": 1, "ntpr": 10,
        "ntr": 0, "nmropt": 0, "disang": None, "reactive_restraints": (),
    }


def _xyz(value: Any) -> tuple[float, float, float]:
    if hasattr(value, "xx"):
        return float(value.xx), float(value.xy), float(value.xz)
    return tuple(float(component) for component in value)  # type: ignore[return-value]


def _distance(coordinates: Mapping[int, Any], left: int, right: int) -> float:
    a, b = _xyz(coordinates[left]), _xyz(coordinates[right])
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def _angle(coordinates: Mapping[int, Any], first: int, center: int, third: int) -> float:
    a, b, c = _xyz(coordinates[first]), _xyz(coordinates[center]), _xyz(coordinates[third])
    u = tuple(x - y for x, y in zip(a, b))
    v = tuple(x - y for x, y in zip(c, b))
    denominator = math.sqrt(sum(x * x for x in u)) * math.sqrt(sum(x * x for x in v))
    if denominator == 0.0:
        raise ValueError("undefined angle from coincident atoms")
    cosine = sum(x * y for x, y in zip(u, v)) / denominator
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def _point_plane_distance(
    coordinates: Mapping[int, Any], point: int, first: int, second: int, third: int
) -> float:
    origin = _xyz(coordinates[first])
    left = _xyz(coordinates[second])
    right = _xyz(coordinates[third])
    target = _xyz(coordinates[point])
    u = tuple(x - y for x, y in zip(left, origin))
    v = tuple(x - y for x, y in zip(right, origin))
    normal = (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )
    norm = math.sqrt(sum(value * value for value in normal))
    if norm == 0.0:
        raise ValueError("reference plane is degenerate")
    displacement = tuple(x - y for x, y in zip(target, origin))
    return abs(sum(x * y for x, y in zip(displacement, normal))) / norm


def measure_product_geometry(coordinates: Mapping[int, Any]) -> dict[str, float]:
    r = REACTIVE_ATOMS
    angle_sum = (
        _angle(coordinates, r["o2"], r["c12"], r["og1"])
        + _angle(coordinates, r["og1"], r["c12"], r["c11"])
        + _angle(coordinates, r["c11"], r["c12"], r["o2"])
    )
    return {
        "product_out_of_plane_A": _point_plane_distance(
            coordinates, r["c12"], r["o2"], r["og1"], r["c11"]
        ),
        "product_angle_sum_deg": angle_sum,
    }


def geometry_from_coordinates(
    coordinates: Mapping[int, Any], qm_heavy_atoms: Iterable[int]
) -> dict[str, Any]:
    r = REACTIVE_ATOMS
    nalpha_hg1 = _distance(coordinates, r["nalpha"], r["hg1"])
    hg1_n3 = _distance(coordinates, r["hg1"], r["n3"])
    nearest = min(qm_heavy_atoms, key=lambda index: _distance(coordinates, r["hg1"], int(index)))
    result: dict[str, Any] = {
        "attack_A": _distance(coordinates, r["og1"], r["c12"]),
        "nalpha_hg1_A": nalpha_hg1,
        "hg1_n3_A": hg1_n3,
        "qPT_A": nalpha_hg1 - hg1_n3,
        "c12_n3_A": _distance(coordinates, r["c12"], r["n3"]),
        "c12_o2_A": _distance(coordinates, r["c12"], r["o2"]),
        "attack_angle_deg": _angle(coordinates, r["o2"], r["c12"], r["og1"]),
        "c12_reactant_plane_out_of_plane_A": _point_plane_distance(
            coordinates, r["c12"], r["o2"], r["n3"], r["c11"]
        ),
        "hg1_nearest_qm_heavy_atom": int(nearest),
    }
    result.update(measure_product_geometry(coordinates))
    return result

def is_product_like(g: Mapping[str, Any]) -> bool:
    return bool(
        PRODUCT_GATE["attack_A"][0] <= g["attack_A"] <= PRODUCT_GATE["attack_A"][1]
        and PRODUCT_GATE["hg1_n3_A"][0] <= g["hg1_n3_A"] <= PRODUCT_GATE["hg1_n3_A"][1]
        and g["nalpha_hg1_A"] >= PRODUCT_GATE["nalpha_hg1_A_min"]
        and g["qPT_A"] >= PRODUCT_GATE["qPT_A_min"]
        and g["c12_n3_A"] >= PRODUCT_GATE["c12_n3_A_min"]
        and PRODUCT_GATE["c12_o2_A"][0] <= g["c12_o2_A"] <= PRODUCT_GATE["c12_o2_A"][1]
        and g["product_out_of_plane_A"] <= PRODUCT_GATE["product_out_of_plane_A_max"]
        and g["product_angle_sum_deg"] >= PRODUCT_GATE["product_angle_sum_deg_min"]
        and g["hg1_nearest_qm_heavy_atom"] == REACTIVE_ATOMS["n3"]
    )


def release_md_persists(frames: Sequence[Mapping[str, Any]]) -> bool:
    if not frames:
        return False
    count = max(1, int(math.ceil(len(frames) * PRODUCT_GATE["last_fraction"])))
    tail = frames[-count:]
    occupancy = sum(is_product_like(frame) for frame in tail) / len(tail)
    return occupancy >= PRODUCT_GATE["occupancy_min"]


def cross_seed_status(outcomes: Sequence[str]) -> str:
    if len(outcomes) != 2:
        raise ValueError("cross-seed denominator must be exactly two")
    if any(outcome == "NOT_EVALUATED_TECHNICAL_FAILURE" for outcome in outcomes):
        return "NOT_EVALUATED_TECHNICAL_A1_ACYL_PRODUCT_ENDPOINT"
    count = sum(outcome == "PERSISTS_ACYL_PRODUCT" for outcome in outcomes)
    if count == 2:
        return "PASS_ACYL_PRODUCT_ENDPOINT_REPRODUCED"
    if count == 1:
        return "NOT_REPRODUCED_A1_ACYL_PRODUCT_ENDPOINT"
    return "FAIL_NO_RELEASE_STABLE_A1_ACYL_PRODUCT_ENDPOINT"


def distance_restraint(first: int, second: int, target: float, force: float) -> str:
    return (
        f"&rst iat={first},{second}, r1={max(0.1, target - 0.35):.3f}, "
        f"r2={target - 0.05:.3f}, r3={target + 0.05:.3f}, r4=4.500, "
        f"rk2={force:.1f}, rk3={force:.1f}, /\n"
    )


def reactive_restraints(targets: Mapping[str, float], force: float) -> str:
    r = REACTIVE_ATOMS
    return "".join((
        distance_restraint(r["og1"], r["c12"], targets["attack"], force),
        distance_restraint(r["hg1"], r["n3"], targets["hg1_n3"], force),
        distance_restraint(r["c12"], r["n3"], targets["cn"], force),
    ))


def minimization_input(stage: str, seed: str, qmmask: str) -> str:
    spec = stage_spec(stage, seed)
    disang = "\nDISANG=restraints.RST\nDUMPAVE=restraint.dat" if spec["disang"] else ""
    return f"""NylC A1 acyl endpoint {stage} {seed}
&cntrl
  imin=1, ntmin=2, maxcyc={spec['maxcyc']}, ncyc={spec['ncyc']}, dx0=0.005,
  ntb=1, cut=10.0, ntpr=5, ntxo=1, ifqnt=1,
  ntr={spec['ntr']}, nmropt={spec['nmropt']},
  restraint_wt=1.0,
  restraintmask='{BASE.NON_QM_SOLUTE_HEAVY_MASK}',
/
{BASE.qmmm_block(qmmask)}&wt type='END' /{disang}
"""


def md_input(seed: str, qmmask: str, full_release: Mapping[str, Any]) -> str:
    spec = release_md_spec(seed, full_release)
    return f"""NylC A1 acyl endpoint release MD {seed}
&cntrl
  imin=0, irest=0, ntx=1, nstlim={spec['nstlim']}, dt={spec['dt']},
  ntb=1, cut=10.0, ntt={spec['ntt']}, gamma_ln={spec['gamma_ln']},
  temp0={spec['temp0']}, tempi={spec['temp0']}, ig={spec['ig']},
  ntc={spec['ntc']}, ntf={spec['ntf']}, ntwx={spec['ntwx']}, ntpr={spec['ntpr']},
  ntxo=1, ifqnt=1, ntr=0, nmropt=0,
/
{BASE.qmmm_block(qmmask)}
"""


def _coordinate_map(structure: Any) -> dict[int, Any]:
    return {atom.idx + 1: atom for atom in structure.atoms}


def _qm_indices(qmmask: str) -> list[int]:
    tokens = [token.strip() for token in qmmask.split(",") if token.strip()]
    values = [int(token[1:] if token.startswith("@") else token) for token in tokens]
    if (
        len(values) != BASE.EXPECTED_QM_ATOMS
        or len(set(values)) != len(values)
        or any(index < 1 or index > EXPECTED_SYSTEM_ATOMS for index in values)
    ):
        raise ValueError("frozen qmmask atom identity/count/range changed")
    return values


def _qm_heavy_indices(structure: Any, qmmask: str) -> list[int]:
    return [
        index for index in _qm_indices(qmmask)
        if int(getattr(structure.atoms[index - 1], "atomic_number", 0)) > 1
    ]


def _validate_structure(structure: Any) -> dict[int, str]:
    if len(structure.atoms) != EXPECTED_SYSTEM_ATOMS or structure.box is None:
        raise ValueError("source atom count or periodic box changed")
    expected = (
        (8949, "THR", "N"), (8960, "THR", "OG1"), (8961, "THR", "HG1"),
        (10286, "L2", "C11"), (10287, "L2", "C12"),
        (10288, "L2", "O2"), (10289, "L2", "N3"),
    )
    for index, residue, name in expected:
        atom = structure.atoms[index - 1]
        if atom.residue.name != residue or atom.name != name:
            raise ValueError(f"reactive atom identity changed at {index}")
    n3 = structure.atoms[REACTIVE_ATOMS["n3"] - 1]
    neighbors: dict[int, str] = {}
    for bond in n3.bonds:
        other = bond.atom2 if bond.atom1 is n3 else bond.atom1
        atomic_number = int(getattr(other, "atomic_number", 0))
        element = "H" if atomic_number == 1 else "C" if atomic_number == 6 else str(atomic_number)
        neighbors[other.idx + 1] = element
    if not validate_n3_reactant_graph(neighbors):
        raise ValueError(f"N3 reactant bond graph changed: {neighbors}")
    return neighbors


def initialize(seed_index: int, output: pathlib.Path, start_rst7: pathlib.Path, commit: str) -> None:
    source = source_from_index(seed_index)
    authority = validate_source_authority(source)
    if not isinstance(authority, tuple) or len(authority) != 2:
        raise ValueError("base authority did not return qmmask and source coordinate")
    qmmask, source_gro = authority
    require_fresh_candidate_dir(output)
    if start_rst7.exists() or not start_rst7.parent.is_dir():
        raise ValueError("scratch start restart target is not fresh")
    import parmed as pmd
    structure = pmd.load_file(str(PRMTOP), xyz=str(source_gro))
    neighbors = _validate_structure(structure)
    qm_heavy = _qm_heavy_indices(structure, qmmask)
    structure.save(str(start_rst7), overwrite=False)
    output.mkdir(parents=True, exist_ok=False)
    manifest = {
        "schema_version": 1, "status": "READY_A1_ACYL_ENDPOINT_STABILITY",
        "scientific_status": TERMINAL_BOUNDARIES[0], "next_action": TERMINAL_BOUNDARIES[1],
        "github_commit": commit, "seed_index": int(seed_index), "seed": source["seed"],
        "source": dict(source), "source_coordinate": str(source_gro),
        "source_coordinate_sha256": sha256(pathlib.Path(source_gro)),
        "prmtop": str(PRMTOP), "prmtop_sha256": sha256(PRMTOP),
        "qm_contract": {
            "qm_atom_count": BASE.EXPECTED_QM_ATOMS, "qmcharge": BASE.QMCHARGE,
            "electron_count_including_link_h": BASE.EXPECTED_ELECTRONS,
            "link_atom_count": BASE.EXPECTED_LINK_ATOMS, "step1_qm_water_count": 0,
            "qmmask": qmmask, "qm_heavy_atom_indices": qm_heavy,
        },
        "n3_reactant_neighbors": {str(k): v for k, v in sorted(neighbors.items())},
        "start_restart_sha256": sha256(start_rst7), "stages": [],
    }
    write_json(output / "ENDPOINT_MANIFEST.json", manifest)


def prepare_stage(
    stage: str, output: pathlib.Path, scratch: pathlib.Path,
    input_restart: pathlib.Path, previous_result: pathlib.Path | None,
) -> None:
    manifest_path = output / "ENDPOINT_MANIFEST.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected_index = len(manifest["stages"])
    if STAGE_ORDER[expected_index] != stage:
        raise ValueError(f"stage order violation: expected {STAGE_ORDER[expected_index]}, got {stage}")
    input_sha = sha256(input_restart)
    expected_sha = (
        manifest["start_restart_sha256"] if not manifest["stages"]
        else manifest["stages"][-1]["restart_sha256"]
    )
    if input_sha != expected_sha:
        raise ValueError("stage input restart SHA does not inherit preceding output")
    scratch.mkdir(parents=True, exist_ok=False)
    seed = manifest["seed"]
    qmmask = manifest["qm_contract"]["qmmask"]
    if stage == "release_md":
        if previous_result is None:
            raise ValueError("release_md requires full_release stage result")
        prior = json.loads(previous_result.read_text(encoding="utf-8"))
        spec = release_md_spec(seed, prior)
        if spec["input_restart_sha256"] != input_sha:
            raise ValueError("release_md full_release restart SHA mismatch")
        (scratch / "stage.in").write_text(md_input(seed, qmmask, prior), encoding="utf-8")
    else:
        spec = stage_spec(stage, seed)
        (scratch / "stage.in").write_text(minimization_input(stage, seed, qmmask), encoding="utf-8")
        if spec["reactive_restraints"]:
            (scratch / "restraints.RST").write_text(
                reactive_restraints(spec["targets"], spec["force"]), encoding="utf-8"
            )
    write_json(scratch / "PREPARED.json", {
        "stage": stage, "seed": seed, "input_restart": str(input_restart),
        "input_restart_sha256": input_sha, "spec": spec,
    })


def _hard_hits(text: str) -> dict[str, int]:
    return {
        name: len(re.findall(pattern, text, re.I))
        for name, pattern in BASE.HARD_PATTERNS.items()
    }


def _read_structure_geometry(restart: pathlib.Path, manifest: Mapping[str, Any]) -> dict[str, Any]:
    import parmed as pmd
    structure = pmd.load_file(str(PRMTOP), xyz=str(restart))
    _validate_structure(structure)
    return geometry_from_coordinates(
        _coordinate_map(structure), manifest["qm_contract"]["qm_heavy_atom_indices"]
    )


def _read_md_frames(
    trajectory: pathlib.Path, manifest: Mapping[str, Any]
) -> tuple[list[dict[str, Any]], tuple[str, ...]]:
    import parmed as pmd
    captured: list[str] = []
    reader = None
    with warnings.catch_warnings(record=True) as observed_warnings:
        warnings.simplefilter("always")
        topology = pmd.load_file(str(PRMTOP))
        try:
            if trajectory_format(trajectory) == "NETCDF":
                reader = pmd.amber.NetCDFTraj.open_old(str(trajectory))
            else:
                reader = pmd.amber.AmberMdcrd(
                    str(trajectory), len(topology.atoms), hasbox=True, mode="r"
                )
            coordinate_frames = reader.coordinates
        finally:
            if reader is not None and hasattr(reader, "close"):
                reader.close()
        captured.extend(str(item.message) for item in observed_warnings)
    heavy = manifest["qm_contract"]["qm_heavy_atom_indices"]
    results = []
    for frame in coordinate_frames:
        coordinates = {index + 1: frame[index] for index in range(len(topology.atoms))}
        results.append(geometry_from_coordinates(coordinates, heavy))
    return results, tuple(captured)


def _last_md_nstep(stage_output: str) -> Any:
    matches = re.findall(r"\bNSTEP\s*=\s*(\d+)", stage_output, re.I)
    return int(matches[-1]) if matches else None


def trajectory_format(path: pathlib.Path) -> str:
    with path.open("rb") as handle:
        return "NETCDF" if handle.read(3) == b"CDF" else "AMBER_MDCRD"


def stage_output_complete(stage: str, text: str) -> bool:
    if stage == "release_md":
        return _last_md_nstep(text) == RELEASE_MD_STEPS
    return "FINAL RESULTS" in text and bool(re.search(r"Run\s+done", text))


def release_md_evidence_complete(
    frames: Sequence[Mapping[str, Any]],
    stage_output: str,
    parse_warnings: Sequence[str],
) -> bool:
    expected_frames = RELEASE_MD_STEPS // RELEASE_NTWX
    if len(frames) != expected_frames or parse_warnings:
        return False
    if _last_md_nstep(stage_output) != RELEASE_MD_STEPS:
        return False
    for frame in frames:
        for key, value in frame.items():
            if key == "hg1_nearest_qm_heavy_atom":
                continue
            if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
                return False
    return True


def audit_stage(stage: str, output: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    manifest_path = output / "ENDPOINT_MANIFEST.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    prepared = json.loads((scratch / "PREPARED.json").read_text(encoding="utf-8"))
    if prepared["stage"] != stage:
        raise ValueError("prepared stage identity mismatch")
    stage_out = scratch / "stage.out"
    restart = scratch / "stage.rst7"
    text = stage_out.read_text(encoding="utf-8", errors="replace") if stage_out.is_file() else ""
    hard = _hard_hits(text)
    scc = len(re.findall(r"Convergence could not be achieved", text, re.I))
    vlimit = len(re.findall(r"vlimit\s+exceeded", text, re.I))
    overflow = len(re.findall(r"BOND\s*=\s*\*+", text, re.I))
    complete = stage_output_complete(stage, text)
    geometry: dict[str, Any] = {}
    frames: list[dict[str, Any]] = []
    md_warnings: list[str] = []
    if restart.is_file() and restart.stat().st_size:
        geometry = _read_structure_geometry(restart, manifest)
    if stage == "release_md" and (scratch / "release.mdcrd").is_file():
        try:
            observed = _read_md_frames(scratch / "release.mdcrd", manifest)
            if isinstance(observed, tuple) and len(observed) == 2:
                frames = list(observed[0])
                md_warnings = [str(message) for message in observed[1]]
            else:
                frames = list(observed)
        except Exception as error:
            md_warnings = [f"{type(error).__name__}: {error}"]
            frames = []
    finite = bool(geometry) and all(
        isinstance(value, (int, float)) and math.isfinite(float(value))
        for key, value in geometry.items() if key != "hg1_nearest_qm_heavy_atom"
    )
    md_complete = (
        release_md_evidence_complete(frames, text, md_warnings)
        if stage == "release_md" else True
    )
    technical = bool(
        complete and restart.is_file() and restart.stat().st_size
        and finite and scc == 0 and vlimit == 0 and overflow == 0
        and sum(hard.values()) == 0 and md_complete
    )
    expected_frames = RELEASE_MD_STEPS // RELEASE_NTWX if stage == "release_md" else None
    result = {
        "stage": stage, "technical_pass": technical,
        "input_restart_path": prepared["input_restart"],
        "input_restart_sha256": prepared["input_restart_sha256"],
        "restart_path": str(restart), "restart_sha256": sha256(restart) if restart.is_file() else "",
        "geometry": geometry, "frame_count": len(frames),
        "frames": frames if stage == "release_md" else [],
        "diagnostics": {
            "complete": complete, "scc_warnings": scc, "vlimit_warnings": vlimit,
            "bond_overflow": overflow, "hard_error_hits": hard,
            "release_md_expected_frames": expected_frames,
            "release_md_last_nstep": _last_md_nstep(text) if stage == "release_md" else None,
            "release_md_parse_warnings": md_warnings,
            "release_md_evidence_complete": md_complete if stage == "release_md" else None,
        },
    }
    if len(manifest["stages"]) != STAGE_ORDER.index(stage):
        raise ValueError("manifest stage append order changed")
    manifest["stages"].append(result)
    manifest["status"] = TECHNICAL_PASS if technical else TECHNICAL_FAIL
    write_json(manifest_path, manifest)
    write_json(scratch / "STAGE_RESULT.json", result)
    if not technical:
        raise RuntimeError(f"{stage} did not technically PASS")
    return result

def _classify(manifest: Mapping[str, Any]) -> tuple[str, float]:
    stages = manifest.get("stages", [])
    if len(stages) != len(STAGE_ORDER) or not all(stage.get("technical_pass") for stage in stages):
        return "NOT_EVALUATED_TECHNICAL_FAILURE", 0.0
    md = stages[-1]
    frames = md.get("frames", [])
    count = max(1, int(math.ceil(len(frames) * PRODUCT_GATE["last_fraction"]))) if frames else 1
    tail = frames[-count:]
    occupancy = sum(is_product_like(frame) for frame in tail) / len(tail) if tail else 0.0
    if release_md_persists(frames):
        return "PERSISTS_ACYL_PRODUCT", occupancy
    final = frames[-1] if frames else stages[-1]["geometry"]
    if final.get("hg1_nearest_qm_heavy_atom") != REACTIVE_ATOMS["n3"] and final.get("qPT_A", -99) >= 0.45:
        return "MISROUTED_PROTON", occupancy
    if final.get("c12_n3_A", 0) >= 2.05 and final.get("attack_A", 0) > 1.65:
        return "ZWITTERIONIC_CLEAVAGE", occupancy
    tetrahedral = bool(
        final.get("attack_A", math.inf) <= 1.70
        and 1.28 <= final.get("c12_o2_A", math.inf) <= 1.45
        and final.get("c12_reactant_plane_out_of_plane_A", -math.inf) >= 0.20
        and 1.30 <= final.get("c12_n3_A", math.inf) <= 1.60
        and 90.0 <= final.get("attack_angle_deg", math.inf) <= 130.0
        and final.get("nalpha_hg1_A", math.inf) <= 1.20
        and final.get("qPT_A", math.inf) <= -0.40
    )
    if tetrahedral:
        return "RETURNS_TETRAHEDRAL", occupancy
    if is_product_like(stages[1]["geometry"]):
        return "RESTRAINT_DEPENDENT_PRODUCT", occupancy
    reactant = bool(
        final.get("attack_A", -math.inf) >= 2.50
        and 1.18 <= final.get("c12_o2_A", math.inf) <= 1.30
        and final.get("c12_reactant_plane_out_of_plane_A", math.inf) <= 0.12
        and final.get("c12_n3_A", math.inf) <= 1.50
        and final.get("nalpha_hg1_A", math.inf) <= 1.20
        and final.get("qPT_A", math.inf) <= -0.40
        and final.get("hg1_nearest_qm_heavy_atom") == REACTIVE_ATOMS["nalpha"]
    )
    if reactant:
        return "RETURNS_REACTANT", occupancy
    return "UNCLASSIFIED_RELEASE_ENDPOINT", occupancy

def finalize_seed(output: pathlib.Path, technical_failure: bool = False) -> dict[str, Any]:
    manifest_path = output / "ENDPOINT_MANIFEST.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    classification, occupancy = (
        ("NOT_EVALUATED_TECHNICAL_FAILURE", 0.0)
        if technical_failure else _classify(manifest)
    )
    technical = classification != "NOT_EVALUATED_TECHNICAL_FAILURE"
    result = {
        "schema_version": 1,
        "status": TECHNICAL_PASS if technical else TECHNICAL_FAIL,
        "scientific_status": TERMINAL_BOUNDARIES[0], "next_action": TERMINAL_BOUNDARIES[1],
        "seed": manifest["seed"], "seed_index": manifest["seed_index"],
        "classification": classification, "technical_complete": technical,
        "release_md_final_40pct_product_occupancy": occupancy,
        "product_gate": PRODUCT_GATE,
        "stage_sha_chain": [
            {
                "stage": stage["stage"], "input_restart_sha256": stage["input_restart_sha256"],
                "output_restart_sha256": stage["restart_sha256"],
            }
            for stage in manifest.get("stages", [])
        ],
        "final_geometry": (
            manifest["stages"][-1]["frames"][-1]
            if manifest.get("stages") and manifest["stages"][-1].get("frames")
            else manifest["stages"][-1].get("geometry", {}) if manifest.get("stages") else {}
        ),
    }
    selected = output / ("PASS.json" if technical else "NOT_EVALUATED.json")
    opposite = output / ("NOT_EVALUATED.json" if technical else "PASS.json")
    write_json(output / "RESULT.json", result)
    write_json(selected, result)
    if opposite.exists():
        opposite.unlink()
    return result

def merge_if_ready(output_root: pathlib.Path, array_job: str) -> bool:
    attempts = [output_root / f"attempt_{array_job}_{index}" for index in (0, 1)]
    paths = [attempt / "RESULT.json" for attempt in attempts]
    if not all(path.is_file() for path in paths):
        return False
    results = [json.loads(path.read_text(encoding="utf-8")) for path in paths]
    if {result["seed_index"] for result in results} != {0, 1}:
        raise ValueError("cross-seed audit does not contain the exact two-seed denominator")
    outcomes = [result["classification"] for result in sorted(results, key=lambda item: item["seed_index"])]
    payload = {
        "schema_version": 1, "denominator": 2, "classifications": outcomes,
        "status": cross_seed_status(outcomes),
        "scientific_status": TERMINAL_BOUNDARIES[0], "next_action": TERMINAL_BOUNDARIES[1],
        "per_seed": sorted(results, key=lambda item: item["seed_index"]),
    }
    audit = output_root / "audit" / f"nylc_a1_acyl_endpoint_{array_job}.json"
    audit.parent.mkdir(parents=True, exist_ok=True)
    if audit.exists():
        existing = json.loads(audit.read_text(encoding="utf-8"))
        if existing != payload:
            raise FileExistsError(f"cross-seed audit already exists with different content: {audit}")
        return True
    write_json(audit, payload)
    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=(
        "initialize", "prepare", "audit-stage", "finalize-seed", "merge-if-ready"
    ))
    parser.add_argument("--seed-index", type=int)
    parser.add_argument("--stage", choices=STAGE_ORDER)
    parser.add_argument("--output", type=pathlib.Path)
    parser.add_argument("--output-root", type=pathlib.Path)
    parser.add_argument("--scratch", type=pathlib.Path)
    parser.add_argument("--start-rst7", type=pathlib.Path)
    parser.add_argument("--input-rst7", type=pathlib.Path)
    parser.add_argument("--previous-result", type=pathlib.Path)
    parser.add_argument("--github-commit", default="unknown")
    parser.add_argument("--array-job")
    parser.add_argument("--technical-failure", action="store_true")
    args = parser.parse_args()
    if args.mode == "initialize":
        initialize(args.seed_index, args.output, args.start_rst7, args.github_commit)
    elif args.mode == "prepare":
        prepare_stage(args.stage, args.output, args.scratch, args.input_rst7, args.previous_result)
    elif args.mode == "audit-stage":
        audit_stage(args.stage, args.output, args.scratch)
    elif args.mode == "finalize-seed":
        finalize_seed(args.output, args.technical_failure)
    else:
        merge_if_ready(args.output_root, args.array_job)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
