"""Export evenly sampled, solvent-free GROMACS trajectories for PyMOL."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


DEFAULT_SELECTION = "not resname SOL WAT HOH NA CL K MG CA"


def sample_indices(total_frames: int, requested_frames: int) -> list[int]:
    if total_frames < 1:
        raise ValueError("trajectory contains no frames")
    if requested_frames < 1:
        raise ValueError("requested frame count must be positive")
    if requested_frames > total_frames:
        raise ValueError("requested frame count exceeds available frames")
    return np.linspace(0, total_frames - 1, requested_frames, dtype=int).tolist()


def export_trajectory(
    topology: Path,
    trajectory: Path,
    output_prefix: Path,
    requested_frames: int,
    selection: str,
) -> dict[str, int | float]:
    import MDAnalysis as mda
    from MDAnalysis.coordinates.DCD import DCDWriter

    universe = mda.Universe(str(topology), str(trajectory), refresh_offsets=True)
    atoms = universe.select_atoms(selection)
    if not atoms.n_atoms:
        raise ValueError("solvent-free selection has no atoms")

    indices = sample_indices(len(universe.trajectory), requested_frames)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    universe.trajectory[indices[0]]
    atoms.write(f"{output_prefix}.pdb")

    with DCDWriter(f"{output_prefix}.dcd", atoms.n_atoms) as writer:
        for index in indices:
            universe.trajectory[index]
            writer.write(atoms)

    return {
        "total_frames": len(universe.trajectory),
        "exported_frames": len(indices),
        "selected_atoms": atoms.n_atoms,
        "start_ps": float(universe.trajectory[indices[0]].time),
        "end_ps": float(universe.trajectory[indices[-1]].time),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--frames", type=int, default=50)
    parser.add_argument("--selection", default=DEFAULT_SELECTION)
    args = parser.parse_args()
    report = export_trajectory(
        args.topology,
        args.trajectory,
        args.output_prefix,
        args.frames,
        args.selection,
    )
    for key, value in report.items():
        print(f"{key}={value}")


if __name__ == "__main__":
    main()

