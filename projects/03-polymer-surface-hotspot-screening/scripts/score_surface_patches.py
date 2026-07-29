#!/usr/bin/env python3
"""Score equal-area catalytic and strongest off-target surface patches."""

from __future__ import annotations

import argparse
import csv
import hashlib
import heapq
import json
import math
import sys
from collections import deque
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


ALL_CHANNELS = ("A", "C", "OA", "HD")
MATERIAL_CHANNELS = {
    "PET": ("A", "C", "OA"),
    "NYLON": ("C", "OA", "HD"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_residue_label(label: str) -> str:
    fields = str(label).split(":")
    if len(fields) >= 3:
        return f"{fields[0]}:{fields[2]}"
    if len(fields) == 2:
        return f"{fields[0]}:{fields[1]}"
    return str(label)


def summarize_values(values: np.ndarray, top_fraction: float = 0.20) -> dict:
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("patch values must be a nonempty vector")
    if not np.isfinite(values).all():
        raise ValueError("patch values contain nonfinite entries")
    if not 0 < top_fraction <= 1:
        raise ValueError("top fraction must be in (0, 1]")
    top_count = max(1, int(math.ceil(values.size * top_fraction)))
    top = np.sort(values)[-top_count:]
    return {
        "n": int(values.size),
        "mean": float(values.mean()),
        "median": float(np.median(values)),
        "maximum": float(values.max()),
        "top_fraction": float(top_fraction),
        "top_count": int(top_count),
        "top_fraction_mean": float(top.mean()),
        "top_fraction_cutoff": float(top.min()),
    }


def robust_z(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    median = float(np.median(values))
    mad = float(np.median(np.abs(values - median)))
    scale = 1.4826 * mad
    if scale <= 1e-12:
        scale = float(values.std())
    if scale <= 1e-12:
        return np.zeros(values.shape, dtype=float)
    return (values - median) / scale


def deterministic_farthest_point_seeds(
    coordinates: np.ndarray, eligible_indices: np.ndarray, n_seeds: int
) -> np.ndarray:
    coordinates = np.asarray(coordinates, dtype=float)
    eligible = np.unique(np.asarray(eligible_indices, dtype=np.int64))
    if eligible.size == 0 or n_seeds <= 0:
        return np.array([], dtype=np.int64)
    n_seeds = min(int(n_seeds), eligible.size)
    order = np.lexsort(
        (
            eligible,
            coordinates[eligible, 2],
            coordinates[eligible, 1],
            coordinates[eligible, 0],
        )
    )
    first = int(eligible[order[0]])
    selected = [first]
    minimum_squared_distance = np.sum(
        (coordinates[eligible] - coordinates[first]) ** 2, axis=1
    )
    while len(selected) < n_seeds:
        unselected_mask = ~np.isin(eligible, np.asarray(selected))
        candidates = eligible[unselected_mask]
        candidate_distances = minimum_squared_distance[unselected_mask]
        maximum = float(candidate_distances.max())
        tied = candidates[np.isclose(candidate_distances, maximum)]
        next_seed = int(tied.min())
        selected.append(next_seed)
        squared_distance = np.sum(
            (coordinates[eligible] - coordinates[next_seed]) ** 2, axis=1
        )
        minimum_squared_distance = np.minimum(
            minimum_squared_distance, squared_distance
        )
    return np.asarray(selected, dtype=np.int64)


def grow_connected_patch(
    seed: int,
    neighbor_indptr: np.ndarray,
    neighbor_indices: np.ndarray,
    eligible_mask: np.ndarray,
    target_size: int,
) -> np.ndarray:
    eligible_mask = np.asarray(eligible_mask, dtype=bool)
    if target_size <= 0 or seed < 0 or seed >= eligible_mask.size:
        raise ValueError("invalid connected patch request")
    if not eligible_mask[seed]:
        return np.array([], dtype=np.int64)
    queue = deque([int(seed)])
    visited = {int(seed)}
    selected = []
    while queue and len(selected) < target_size:
        node = queue.popleft()
        selected.append(node)
        start = int(neighbor_indptr[node])
        end = int(neighbor_indptr[node + 1])
        for neighbor in neighbor_indices[start:end]:
            neighbor = int(neighbor)
            if eligible_mask[neighbor] and neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)
    if len(selected) != target_size:
        return np.array([], dtype=np.int64)
    return np.asarray(selected, dtype=np.int64)


def multi_source_graph_distances(
    coordinates: np.ndarray,
    neighbor_indptr: np.ndarray,
    neighbor_indices: np.ndarray,
    source_indices: np.ndarray,
) -> np.ndarray:
    coordinates = np.asarray(coordinates, dtype=float)
    n_nodes = coordinates.shape[0]
    distances = np.full(n_nodes, np.inf, dtype=float)
    heap = []
    for source in np.unique(np.asarray(source_indices, dtype=np.int64)):
        distances[source] = 0.0
        heapq.heappush(heap, (0.0, int(source)))
    while heap:
        distance, node = heapq.heappop(heap)
        if distance != distances[node]:
            continue
        start = int(neighbor_indptr[node])
        end = int(neighbor_indptr[node + 1])
        for neighbor in neighbor_indices[start:end]:
            neighbor = int(neighbor)
            weight = float(np.linalg.norm(coordinates[node] - coordinates[neighbor]))
            candidate = distance + weight
            if candidate + 1e-12 < distances[neighbor]:
                distances[neighbor] = candidate
                heapq.heappush(heap, (candidate, neighbor))
    return distances


def candidate_seeds(
    coordinates: np.ndarray,
    eligible_mask: np.ndarray,
    scores: np.ndarray,
    max_seeds: int,
) -> np.ndarray:
    eligible = np.flatnonzero(eligible_mask)
    if eligible.size == 0:
        return np.array([], dtype=np.int64)
    max_seeds = min(max(1, int(max_seeds)), eligible.size)
    n_hot = max_seeds // 2
    hot_order = np.lexsort((eligible, -np.asarray(scores)[eligible]))
    hot = eligible[hot_order[:n_hot]]
    farthest = deterministic_farthest_point_seeds(
        coordinates, eligible, max_seeds - len(hot)
    )
    ordered = []
    seen = set()
    for node in np.concatenate((hot, farthest)):
        node = int(node)
        if node not in seen:
            seen.add(node)
            ordered.append(node)
    if len(ordered) < max_seeds:
        for node in eligible:
            node = int(node)
            if node not in seen:
                seen.add(node)
                ordered.append(node)
            if len(ordered) == max_seeds:
                break
    return np.asarray(ordered, dtype=np.int64)


def strongest_equal_size_patch(
    coordinates: np.ndarray,
    neighbor_indptr: np.ndarray,
    neighbor_indices: np.ndarray,
    eligible_mask: np.ndarray,
    target_size: int,
    scores: np.ndarray,
    top_fraction: float,
    max_seeds: int,
) -> dict | None:
    best = None
    seen_patches = set()
    candidate_summaries = []
    for seed in candidate_seeds(coordinates, eligible_mask, scores, max_seeds):
        patch = grow_connected_patch(
            int(seed),
            neighbor_indptr,
            neighbor_indices,
            eligible_mask,
            target_size,
        )
        if patch.size != target_size:
            continue
        key = tuple(sorted(patch.tolist()))
        if key in seen_patches:
            continue
        seen_patches.add(key)
        summary = summarize_values(np.asarray(scores)[patch], top_fraction)
        candidate_summaries.append(summary["top_fraction_mean"])
        rank = (
            summary["top_fraction_mean"],
            summary["mean"],
            -int(seed),
        )
        if best is None or rank > best["rank"]:
            best = {
                "seed": int(seed),
                "indices": patch,
                "summary": summary,
                "rank": rank,
            }
    if best is None:
        return None
    best["candidate_top_fraction_means"] = candidate_summaries
    return best


def _validate_shell_arrays(
    coordinates: np.ndarray,
    neighbor_indptr: np.ndarray,
    neighbor_indices: np.ndarray,
    channels: dict[str, np.ndarray],
    nearest_residue: np.ndarray,
) -> None:
    n_nodes = len(coordinates)
    if np.asarray(coordinates).shape != (n_nodes, 3):
        raise ValueError("coordinates must be N by 3")
    if len(neighbor_indptr) != n_nodes + 1:
        raise ValueError("invalid CSR indptr")
    if int(neighbor_indptr[-1]) != len(neighbor_indices):
        raise ValueError("invalid CSR indices")
    if len(nearest_residue) != n_nodes:
        raise ValueError("nearest residue count differs")
    for channel in ALL_CHANNELS:
        if channel not in channels or len(channels[channel]) != n_nodes:
            raise ValueError(f"missing or mismatched shell channel {channel}")
        if not np.isfinite(channels[channel]).all():
            raise ValueError(f"nonfinite shell channel {channel}")


def score_shell_data(
    coordinates: np.ndarray,
    neighbor_indptr: np.ndarray,
    neighbor_indices: np.ndarray,
    nearest_residue: np.ndarray,
    channels: dict[str, np.ndarray],
    catalytic_residues: set[str],
    material_family: str,
    region_radii: dict[str, float],
    off_target_buffer: float,
    top_fraction: float,
    max_seeds: int,
    catalytic_anchor_indices: np.ndarray | None = None,
    missing_catalytic_residues: list[str] | None = None,
) -> dict:
    material_family = material_family.upper()
    if material_family not in MATERIAL_CHANNELS:
        raise ValueError("material family must be PET or NYLON")
    _validate_shell_arrays(
        coordinates, neighbor_indptr, neighbor_indices, channels, nearest_residue
    )
    catalytic_residues = {
        normalize_residue_label(label) for label in catalytic_residues
    }
    shell_labels = np.asarray(
        [normalize_residue_label(label) for label in nearest_residue]
    )
    if catalytic_anchor_indices is None:
        present = set(shell_labels) & catalytic_residues
        missing = sorted(catalytic_residues - present)
        anchor_indices = np.flatnonzero(np.isin(shell_labels, sorted(present)))
    else:
        anchor_indices = np.unique(
            np.asarray(catalytic_anchor_indices, dtype=np.int64)
        )
        missing = sorted(missing_catalytic_residues or [])
    if missing or anchor_indices.size == 0:
        return {
            "status": "NOT_EVALUATED_CATALYTIC_MAPPING",
            "material_family": material_family,
            "material_channels": list(MATERIAL_CHANNELS[material_family]),
            "missing_catalytic_residues": missing or sorted(catalytic_residues),
        }

    coordinates = np.asarray(coordinates, dtype=float)
    raw_channels = {
        channel: np.asarray(channels[channel], dtype=float)
        for channel in ALL_CHANNELS
    }
    affinities = {channel: -raw_channels[channel] for channel in ALL_CHANNELS}
    material_channels = MATERIAL_CHANNELS[material_family]
    composite = np.mean(
        np.vstack([robust_z(affinities[channel]) for channel in material_channels]),
        axis=0,
    )
    catalytic_distances = multi_source_graph_distances(
        coordinates, neighbor_indptr, neighbor_indices, anchor_indices
    )
    result = {
        "status": "PATCH_SCORE_PASS",
        "scientific_scope": "equal_area_surface_affinity_proxy_only",
        "material_family": material_family,
        "material_channels": list(material_channels),
        "catalytic_residues": sorted(catalytic_residues),
        "n_shell_points": int(len(coordinates)),
        "region_radii": {key: float(value) for key, value in region_radii.items()},
        "off_target_buffer": float(off_target_buffer),
        "top_fraction": float(top_fraction),
        "max_seeds": int(max_seeds),
        "regions": {},
        "_memberships": {},
    }
    for region_name, radius in region_radii.items():
        catalytic_indices = np.flatnonzero(catalytic_distances <= float(radius))
        if catalytic_indices.size == 0:
            return {
                "status": "NOT_EVALUATED_CATALYTIC_PATCH_EMPTY",
                "material_family": material_family,
                "region": region_name,
            }
        eligible_mask = catalytic_distances > float(radius) + float(off_target_buffer)
        target_size = int(catalytic_indices.size)
        region = {
            "radius": float(radius),
            "n_catalytic_points": target_size,
            "n_off_target_points": target_size,
            "channels": {},
        }
        result["_memberships"][f"catalytic__{region_name}"] = catalytic_indices

        composite_catalytic = summarize_values(
            composite[catalytic_indices], top_fraction
        )
        composite_off = strongest_equal_size_patch(
            coordinates,
            neighbor_indptr,
            neighbor_indices,
            eligible_mask,
            target_size,
            composite,
            top_fraction,
            max_seeds,
        )
        if composite_off is None:
            return {
                "status": "NOT_EVALUATED_EQUAL_AREA_OFF_TARGET",
                "material_family": material_family,
                "region": region_name,
                "target_size": target_size,
                "n_eligible_points": int(eligible_mask.sum()),
            }
        composite_candidates = np.asarray(
            composite_off["candidate_top_fraction_means"], dtype=float
        )
        region["composite"] = {
            "catalytic": composite_catalytic,
            "off_target": composite_off["summary"],
            "delta_top_fraction_mean": float(
                composite_catalytic["top_fraction_mean"]
                - composite_off["summary"]["top_fraction_mean"]
            ),
            "catalytic_candidate_percentile": float(
                np.mean(
                    composite_candidates
                    <= composite_catalytic["top_fraction_mean"]
                )
            ),
        }
        result["_memberships"][
            f"off_target__{region_name}__composite"
        ] = composite_off["indices"]
        threshold = composite_catalytic["top_fraction_cutoff"]
        region["global_sticky_fraction"] = float(np.mean(composite >= threshold))

        for channel in ALL_CHANNELS:
            catalytic_summary = summarize_values(
                affinities[channel][catalytic_indices], top_fraction
            )
            off = strongest_equal_size_patch(
                coordinates,
                neighbor_indptr,
                neighbor_indices,
                eligible_mask,
                target_size,
                affinities[channel],
                top_fraction,
                max_seeds,
            )
            if off is None:
                return {
                    "status": "NOT_EVALUATED_EQUAL_AREA_OFF_TARGET",
                    "material_family": material_family,
                    "region": region_name,
                    "channel": channel,
                    "target_size": target_size,
                }
            candidate_values = np.asarray(
                off["candidate_top_fraction_means"], dtype=float
            )
            region["channels"][channel] = {
                "catalytic": catalytic_summary,
                "off_target": off["summary"],
                "delta_mean": float(
                    catalytic_summary["mean"] - off["summary"]["mean"]
                ),
                "delta_top_fraction_mean": float(
                    catalytic_summary["top_fraction_mean"]
                    - off["summary"]["top_fraction_mean"]
                ),
                "catalytic_candidate_percentile": float(
                    np.mean(
                        candidate_values
                        <= catalytic_summary["top_fraction_mean"]
                    )
                ),
            }
            result["_memberships"][
                f"off_target__{region_name}__{channel}"
            ] = off["indices"]
        result["regions"][region_name] = region
    return result


def catalytic_anchors_from_pdbqt(
    receptor_pdbqt: Path,
    catalytic_residues: set[str],
    shell_coordinates: np.ndarray,
) -> tuple[np.ndarray, list[str]]:
    atoms_by_residue = {}
    with receptor_pdbqt.open(encoding="ascii", errors="strict") as handle:
        for line in handle:
            if line[:6] not in {"ATOM  ", "HETATM"}:
                continue
            chain = line[21:22].strip() or "_"
            residue_number = line[22:27].strip()
            key = f"{chain}:{residue_number}"
            if key not in catalytic_residues:
                continue
            coordinate = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                dtype=float,
            )
            atoms_by_residue.setdefault(key, []).append(coordinate)
    missing = sorted(catalytic_residues - set(atoms_by_residue))
    if missing:
        return np.array([], dtype=np.int64), missing
    tree = cKDTree(np.asarray(shell_coordinates, dtype=float))
    anchors = []
    for residue in sorted(catalytic_residues):
        centroid = np.mean(np.vstack(atoms_by_residue[residue]), axis=0)
        anchors.append(int(tree.query(centroid, k=1)[1]))
    return np.unique(np.asarray(anchors, dtype=np.int64)), []


def parse_region_radii(text: str) -> dict[str, float]:
    output = {}
    for item in text.split(","):
        name, value = item.split("=", maxsplit=1)
        output[name.strip()] = float(value)
    if not output or any(value <= 0 for value in output.values()):
        raise ValueError("invalid region radii")
    return output


def write_json_exclusive(path: Path, payload: dict) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shell-dir", type=Path, required=True)
    parser.add_argument("--receptor-pdbqt", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--record-id", required=True)
    parser.add_argument("--material-family", choices=("PET", "NYLON"), required=True)
    parser.add_argument("--catalytic-residues", required=True)
    parser.add_argument(
        "--region-radii",
        default="catalytic_core=6,catalytic_neighborhood=10,catalytic_extended=14",
    )
    parser.add_argument("--off-target-buffer", type=float, default=4.0)
    parser.add_argument("--top-fraction", type=float, default=0.20)
    parser.add_argument("--max-seeds", type=int, default=128)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path("/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    shell_gate_path = args.shell_dir / "SHELL_PASS.json"
    shell_npz_path = args.shell_dir / "surface_shell.npz"
    if any(
        not path.is_file()
        for path in (shell_gate_path, shell_npz_path, args.receptor_pdbqt)
    ):
        print("one or more required patch inputs are missing", file=sys.stderr)
        return 2
    root = args.project_root.resolve()
    output = args.output_dir.resolve()
    if root != output and root not in output.parents:
        print("refusing output outside project root", file=sys.stderr)
        return 2
    if args.output_dir.exists():
        pass_gate = args.output_dir / "PATCH_PASS.json"
        if pass_gate.is_file():
            print(f"existing PATCH_PASS retained: {pass_gate}")
            return 0
        print("refusing partial or unverified patch output", file=sys.stderr)
        return 3

    args.output_dir.mkdir(parents=True)
    try:
        shell_gate = json.loads(shell_gate_path.read_text(encoding="utf-8"))
        if shell_gate.get("status") != "SHELL_PASS":
            raise ValueError("incompatible shell pass gate")
        catalytic_residues = {
            normalize_residue_label(value.strip())
            for value in args.catalytic_residues.split(",")
            if value.strip()
        }
        if not catalytic_residues:
            raise ValueError("no catalytic residues")
        with np.load(shell_npz_path, allow_pickle=False) as archive:
            coordinates = archive["coordinates"]
            indptr = archive["neighbor_indptr"]
            indices = archive["neighbor_indices"]
            nearest_residue = archive["nearest_residue"]
            channels = {channel: archive[channel] for channel in ALL_CHANNELS}
        anchors, missing = catalytic_anchors_from_pdbqt(
            args.receptor_pdbqt, catalytic_residues, coordinates
        )
        result = score_shell_data(
            coordinates,
            indptr,
            indices,
            nearest_residue,
            channels,
            catalytic_residues,
            args.material_family,
            parse_region_radii(args.region_radii),
            args.off_target_buffer,
            args.top_fraction,
            args.max_seeds,
            catalytic_anchor_indices=anchors,
            missing_catalytic_residues=missing,
        )
        memberships = result.pop("_memberships", {})
        result.update(
            {
                "record_id": args.record_id,
                "shell_gate_sha256": sha256(shell_gate_path),
                "shell_npz_sha256": sha256(shell_npz_path),
                "receptor_pdbqt_sha256": sha256(args.receptor_pdbqt),
            }
        )
        if result["status"] != "PATCH_SCORE_PASS":
            write_json_exclusive(
                args.output_dir / "PATCH_NOT_EVALUATED.json", result
            )
            print(json.dumps(result, sort_keys=True), file=sys.stderr)
            return 6

        np.savez_compressed(
            args.output_dir / "patch_membership.npz",
            **{key: np.asarray(value, dtype=np.int32) for key, value in memberships.items()},
        )
        rows = []
        for region_name, region in result["regions"].items():
            for channel, metrics in {
                **region["channels"],
                "COMPOSITE": region["composite"],
            }.items():
                rows.append(
                    {
                        "record_id": args.record_id,
                        "material_family": args.material_family,
                        "region": region_name,
                        "channel": channel,
                        "n_catalytic_points": region["n_catalytic_points"],
                        "n_off_target_points": region["n_off_target_points"],
                        "catalytic_mean": metrics["catalytic"]["mean"],
                        "catalytic_top_fraction_mean": metrics["catalytic"]["top_fraction_mean"],
                        "off_target_mean": metrics["off_target"]["mean"],
                        "off_target_top_fraction_mean": metrics["off_target"]["top_fraction_mean"],
                        "delta_top_fraction_mean": metrics["delta_top_fraction_mean"],
                        "catalytic_candidate_percentile": metrics["catalytic_candidate_percentile"],
                        "global_sticky_fraction": region["global_sticky_fraction"],
                    }
                )
        with (args.output_dir / "patches.tsv").open(
            "x", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(rows[0]), delimiter="\t"
            )
            writer.writeheader()
            writer.writerows(rows)
        write_json_exclusive(args.output_dir / "PATCH_PASS.json", result)
        checksums = []
        for path in sorted(args.output_dir.iterdir()):
            if path.is_file() and path.name != "SHA256SUMS":
                checksums.append(f"{sha256(path)}  {path.name}")
        (args.output_dir / "SHA256SUMS").write_text(
            "\n".join(checksums) + "\n", encoding="ascii"
        )
        print(json.dumps(result, sort_keys=True))
        return 0
    except Exception as error:
        write_json_exclusive(
            args.output_dir / "PATCH_FAIL.json",
            {
                "status": "NOT_EVALUATED_PATCH_GENERATION",
                "record_id": args.record_id,
                "error": str(error),
            },
        )
        print(str(error), file=sys.stderr)
        return 5


if __name__ == "__main__":
    raise SystemExit(main())
