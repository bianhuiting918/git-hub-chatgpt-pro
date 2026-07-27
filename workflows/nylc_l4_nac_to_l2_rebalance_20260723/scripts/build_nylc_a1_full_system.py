#!/usr/bin/env python3
"""Patch full NylC Thr267 from NalphaH2/OgammaH to NalphaH3+/Ogammaminus."""

import argparse
import hashlib
import json
import math
import pathlib
import re

SECTION = re.compile(r"(?m)^\s*\[\s*([^\]]+?)\s*\]\s*(?:;.*)?$")
ARITY = {"bonds": 2, "pairs": 2, "angles": 3, "dihedrals": 4}
REACTIVE = {"N", "H1", "H2", "HG1", "OG1"}
THR_NAMES = {"N","H1","H2","CA","HA","CB","HB","CG2","HG21","HG22","HG23","OG1","HG1","C","O"}
CAP_MAP = {"NCAP": (268,"N"), "HCAP": (268,"H"), "CCAP": (268,"CA"), "HC1": (268,"HA")}


class BuildError(RuntimeError):
    pass


def _sections(text):
    matches = list(SECTION.finditer(text))
    if not matches:
        raise BuildError("no topology sections")
    blocks = []
    for i, match in enumerate(matches):
        end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        blocks.append({
            "name": match.group(1).strip().lower(),
            "header": text[match.start():match.end()],
            "body": text[match.end():end],
        })
    return text[:matches[0].start()], blocks


def _join(preamble, blocks):
    return preamble + "".join(x["header"] + x["body"] for x in blocks)


def _fields(line):
    return line.split(";", 1)[0].split()


def _atoms(blocks):
    targets = [x for x in blocks if x["name"] == "atoms"]
    if len(targets) != 1:
        raise BuildError("expected exactly one atoms section")
    result = {}
    for line in targets[0]["body"].splitlines():
        f = _fields(line)
        if len(f) >= 8 and f[0].isdigit():
            result[int(f[0])] = {
                "id": int(f[0]), "type": f[1], "resnr": int(f[2]),
                "residue": f[3], "name": f[4], "cgnr": int(f[5]),
                "charge": float(f[6]), "mass": float(f[7]), "extra": f[8:],
            }
    if not result:
        raise BuildError("empty atoms section")
    return result


def _atom_line(a):
    out = (
        f"{a['id']:6d} {a['type']:>10s} {a['resnr']:6d} {a['residue']:>6s} "
        f"{a['name']:>6s} {a['cgnr']:6d} {a['charge']: .8f} {a['mass']: .6f}"
    )
    if a["extra"]:
        out += " " + " ".join(a["extra"])
    return out + "\n"


def _replace_atoms(body, replacement):
    output, seen = [], set()
    for line in body.splitlines(keepends=True):
        f = _fields(line)
        if len(f) >= 8 and f[0].isdigit() and int(f[0]) in replacement:
            atom_id = int(f[0])
            output.append(_atom_line(replacement[atom_id]))
            seen.add(atom_id)
        else:
            output.append(line)
    if seen != set(replacement):
        raise BuildError("not all Thr267 atom rows were replaced")
    return "".join(output)


def _interaction(line, arity):
    f = _fields(line)
    if len(f) < arity + 1:
        return None
    try:
        ids = tuple(int(x) for x in f[:arity])
        int(f[arity])
    except ValueError:
        return None
    return ids, f


def _model_mapping(model_atoms, parent_atoms):
    lookup = {(x["resnr"], x["name"]): i for i, x in parent_atoms.items()}
    mapping = {}
    for i, atom in model_atoms.items():
        if atom["name"] in CAP_MAP:
            target = CAP_MAP[atom["name"]]
        elif atom["name"] in THR_NAMES:
            target = (267, atom["name"])
        else:
            continue
        if target in lookup:
            mapping[i] = lookup[target]
    return mapping


def _model_terms(model_blocks, mapping, reactive_model):
    terms = {}
    for section, arity in ARITY.items():
        rows = []
        for block in model_blocks:
            if block["name"] != section:
                continue
            for line in block["body"].splitlines():
                parsed = _interaction(line, arity)
                if parsed is None:
                    continue
                ids, f = parsed
                if any(i in reactive_model for i in ids) and all(i in mapping for i in ids):
                    rows.append(" ".join([str(mapping[i]) for i in ids] + f[arity:]) + "\n")
        if not rows:
            raise BuildError(f"no mapped reactive terms in [{section}]")
        terms[section] = rows
    return terms


def _replace_terms(blocks, section, reactive_parent, additions):
    indexes = [i for i, x in enumerate(blocks) if x["name"] == section]
    if not indexes:
        raise BuildError(f"missing parent [{section}]")
    arity = ARITY[section]
    for index in indexes:
        kept = []
        for line in blocks[index]["body"].splitlines(keepends=True):
            parsed = _interaction(line, arity)
            if parsed is not None and any(i in reactive_parent for i in parsed[0]):
                continue
            kept.append(line)
        blocks[index]["body"] = "".join(kept)
    target = indexes[0]
    if not blocks[target]["body"].endswith("\n"):
        blocks[target]["body"] += "\n"
    blocks[target]["body"] += "; A1 terms from audited AM1-BCC/GAFF2 capped model\n"
    blocks[target]["body"] += "".join(additions)


def _bonds(blocks):
    result = set()
    for block in blocks:
        if block["name"] == "bonds":
            for line in block["body"].splitlines():
                parsed = _interaction(line, 2)
                if parsed:
                    result.add(tuple(sorted(parsed[0])))
    return result


def patch_chain_itp(parent_text, model_top_text):
    preamble, parent_blocks = _sections(parent_text)
    _, model_blocks = _sections(model_top_text)
    parent = _atoms(parent_blocks)
    model = _atoms(model_blocks)
    thr = {i: x for i, x in parent.items() if x["resnr"] == 267}
    by_name = {x["name"]: i for i, x in thr.items()}
    model_by_name = {x["name"]: x for x in model.values() if x["name"] in THR_NAMES}
    if len(thr) != 15 or set(by_name) != THR_NAMES or set(model_by_name) != THR_NAMES:
        raise BuildError("Thr267/model atom-name contract failed")

    count_before = len(parent)
    charge_before = sum(x["charge"] for x in parent.values())
    raw_charge = sum(model_by_name[name]["charge"] for name in THR_NAMES)
    correction = -raw_charge
    replacements = {}
    for name, atom_id in by_name.items():
        row = dict(parent[atom_id])
        row["charge"] = model_by_name[name]["charge"]
        if name in REACTIVE:
            row["type"] = model_by_name[name]["type"]
        replacements[atom_id] = row
    replacements[by_name["C"]]["charge"] += correction / 2
    replacements[by_name["O"]]["charge"] += correction / 2
    [x for x in parent_blocks if x["name"] == "atoms"][0]["body"] = _replace_atoms(
        [x for x in parent_blocks if x["name"] == "atoms"][0]["body"], replacements
    )

    mapping = _model_mapping(model, parent)
    reactive_model = {i for i, x in model.items() if x["name"] in REACTIVE}
    reactive_parent = {by_name[name] for name in REACTIVE}
    terms = _model_terms(model_blocks, mapping, reactive_model)
    for section in ARITY:
        _replace_terms(parent_blocks, section, reactive_parent, terms[section])

    patched = _join(preamble, parent_blocks)
    _, checked_blocks = _sections(patched)
    checked = _atoms(checked_blocks)
    bonds = _bonds(checked_blocks)
    n, og, hg = by_name["N"], by_name["OG1"], by_name["HG1"]
    degree = {i: 0 for i in checked}
    for first, second in bonds:
        degree[first] += 1
        degree[second] += 1
    charge_after = sum(x["charge"] for x in checked.values())
    thr_after = sum(x["charge"] for x in checked.values() if x["resnr"] == 267)
    if tuple(sorted((og, hg))) in bonds or tuple(sorted((n, hg))) not in bonds:
        raise BuildError("A1 proton-transfer bond graph is wrong")
    if (degree[n], degree[og], degree[hg]) != (4, 1, 1):
        raise BuildError(f"A1 valence is wrong: N={degree[n]} OG1={degree[og]} HG1={degree[hg]}")
    if len(checked) != count_before or abs(charge_after - charge_before) > 1e-6 or abs(thr_after) > 1e-6:
        raise BuildError("atom-count or integer-charge invariant failed")
    return patched, {
        "schema_version": 1, "status": "PASS_A1_CHAIN_PATCH",
        "atom_count_before": count_before, "atom_count_after": len(checked),
        "chain_charge_before_e": charge_before, "chain_charge_after_e": charge_after,
        "raw_a1_thr267_charge_e": raw_charge,
        "thr267_charge_correction_e": correction, "thr267_charge_after_e": thr_after,
        "removed_bond": [og, hg], "added_bond": [n, hg],
        "reactive_valence": {"N": degree[n], "OG1": degree[og], "HG1": degree[hg]},
        "reactive_atom_types": {name: checked[by_name[name]]["type"] for name in sorted(REACTIVE)},
        "model_to_parent_atom_map": {str(k): v for k, v in sorted(mapping.items())},
    }


def parse_gro(text):
    lines = text.splitlines()
    if len(lines) < 3:
        raise BuildError("truncated GRO")
    count = int(lines[1].strip())
    if len(lines) != count + 3:
        raise BuildError("GRO atom-count mismatch")
    atoms = {}
    for sequential_id, line in enumerate(lines[2:2 + count], 1):
        if len(line) < 44:
            raise BuildError(f"short GRO atom line {sequential_id}")
        atoms[sequential_id] = {
            "line": line, "display_atom_id": int(line[15:20]),
            "xyz_nm": [float(line[20:28]), float(line[28:36]), float(line[36:44])],
        }
    return {"title": lines[0], "count": count, "atoms": atoms, "box": lines[-1]}


def _mol2_xyz(text, atom_name):
    in_atoms, hits = False, []
    for line in text.splitlines():
        if line == "@<TRIPOS>ATOM":
            in_atoms = True
            continue
        if in_atoms and line.startswith("@<TRIPOS>"):
            break
        if in_atoms:
            f = line.split()
            if len(f) >= 5 and f[1] == atom_name:
                hits.append([float(f[2]), float(f[3]), float(f[4])])
    if len(hits) != 1:
        raise BuildError(f"expected one MOL2 atom {atom_name}")
    return hits[0]


def patch_gro_coordinates(source_text, model_mol2_text, chain_first_global_atom):
    before = parse_gro(source_text)
    target = chain_first_global_atom + 12
    if target not in before["atoms"]:
        raise BuildError("HG1 global atom is outside GRO")
    xyz = [x / 10 for x in _mol2_xyz(model_mol2_text, "HG1")]
    if not all(math.isfinite(x) for x in xyz):
        raise BuildError("non-finite HG1 coordinate")
    lines = source_text.splitlines()
    old = lines[target + 1]
    lines[target + 1] = old[:20] + "".join(f"{x:8.3f}" for x in xyz) + old[44:]
    patched = "\n".join(lines) + ("\n" if source_text.endswith("\n") else "")
    after = parse_gro(patched)
    changed = [i for i in before["atoms"] if before["atoms"][i]["xyz_nm"] != after["atoms"][i]["xyz_nm"]]
    if changed != [target] or after["count"] != before["count"]:
        raise BuildError(f"unexpected coordinate changes: {changed}")
    return patched, {
        "schema_version": 1, "status": "PASS_A1_COORDINATE_PATCH",
        "atom_count_before": before["count"], "atom_count_after": after["count"],
        "transferred_global_atom_id": target, "old_xyz_nm": before["atoms"][target]["xyz_nm"],
        "new_xyz_nm": after["atoms"][target]["xyz_nm"], "changed_sequential_atom_ids": changed,
    }


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-gro", required=True)
    ap.add_argument("--parent-chain-itp", required=True)
    ap.add_argument("--model-top", required=True)
    ap.add_argument("--model-mol2", required=True)
    ap.add_argument("--chain-first-global-atom", required=True, type=int)
    ap.add_argument("--candidate-id", required=True)
    ap.add_argument("--output-dir", required=True)
    a = ap.parse_args()
    out = pathlib.Path(a.output_dir)
    if out.exists():
        raise SystemExit(f"refusing to overwrite {out}")
    out.mkdir(parents=True)
    gro, itp, mt, mm = map(pathlib.Path, [a.source_gro, a.parent_chain_itp, a.model_top, a.model_mol2])
    patched_itp, chain_audit = patch_chain_itp(itp.read_text(), mt.read_text())
    patched_gro, coordinate_audit = patch_gro_coordinates(gro.read_text(), mm.read_text(), a.chain_first_global_atom)
    itp_out, gro_out = out / "topol_Protein_chain_H.itp", out / "system_A1.gro"
    itp_out.write_text(patched_itp)
    gro_out.write_text(patched_gro)
    audit = {
        "schema_version": 1, "status": "PASS_A1_FULL_SYSTEM_BUILD",
        "candidate_id": a.candidate_id,
        "scientific_scope": "fixed-topology classical A1 NAC stability screening only",
        "chain_patch": chain_audit, "coordinate_patch": coordinate_audit,
        "inputs": {str(p): _sha(p) for p in [gro, itp, mt, mm]},
        "outputs": {itp_out.name: _sha(itp_out), gro_out.name: _sha(gro_out)},
    }
    (out / "A1_FULL_SYSTEM_BUILD.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
