from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def available(module: str) -> str:
    return "yes" if importlib.util.find_spec(module) is not None else "no"


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: write_ligand_build_environment_report.py <out-md>", file=sys.stderr)
        return 2
    out = Path(sys.argv[1])
    out.parent.mkdir(parents=True, exist_ok=True)
    rdkit = available("rdkit")
    openbabel = available("openbabel")
    bio = available("Bio")
    status = "ready" if rdkit == "yes" or openbabel == "yes" else "blocked_missing_chemistry_toolkit"
    out.write_text(
        "\n".join(
            [
                "# Ligand Build Environment Report",
                "",
                "Date: 2026-06-30",
                "",
                "## Boundary",
                "",
                "This report concerns only Stage 1 ligand construction environment readiness. It uses no paper-derived TS, RC, CV, trajectory, barrier, or mechanism result.",
                "",
                "## Python Runtime",
                "",
                f"- Python: `{sys.version.split()[0]}`",
                f"- RDKit import available: `{rdkit}`",
                f"- Open Babel Python import available: `{openbabel}`",
                f"- Biopython import available: `{bio}`",
                "",
                "## Status",
                "",
                f"`{status}`",
                "",
                "Ligand 3D conformer generation and atom-label validation require RDKit or Open Babel. The current bundled Python environment has neither, so this stage is intentionally stopped before producing unreliable ligand coordinates.",
                "",
                "## Next Accepted Routes",
                "",
                "1. Install or expose RDKit/Open Babel in the compute environment, then generate SDF/MOL2 from `ligand_model_definitions.md`.",
                "2. Use an existing cluster chemistry stack, but record the exact executable path, version, command, and generated atom-label table.",
                "3. If ligand coordinates are built manually in a GUI, export SDF/MOL2 plus an atom-label TSV and mark the build as manual in `ligand_model_manifest.tsv`.",
                "",
                "## Grill Gate",
                "",
                "Do not start docking or QM/MM with ligand coordinates unless the scissile ester labels in `ligand_model_manifest.tsv` are mapped to concrete atom names in the generated topology.",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
