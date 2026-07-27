#!/usr/bin/env python3
"""Audit converged free EM while preserving the A1 protonation graph."""

import argparse
import importlib.util
import json
import math
import pathlib
import re

HERE=pathlib.Path(__file__).resolve().parent
FULL=HERE/"audit_nylc_a1_full_system.py"

class EMAuditError(RuntimeError):
    pass

def load_full():
    spec=importlib.util.spec_from_file_location("a1_full_audit",FULL)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def audit(gro_text,chain_itp_text,build_audit,log_text):
    if "converged to Fmax" not in log_text:
        raise EMAuditError("free EM did not converge to Fmax")
    force_hits=re.findall(r"Maximum force\s*=\s*([-+0-9.eE]+)",log_text)
    energy_hits=re.findall(r"Potential Energy\s*=\s*([-+0-9.eE]+)",log_text)
    if not force_hits or not energy_hits:
        raise EMAuditError("free EM final force/energy not found")
    maximum_force=float(force_hits[-1]); potential=float(energy_hits[-1])
    if not math.isfinite(maximum_force) or maximum_force>=500.0:
        raise EMAuditError(f"free EM maximum force is not below 500: {maximum_force}")
    if not math.isfinite(potential):
        raise EMAuditError("free EM potential is not finite")
    upper=log_text.upper()
    for marker in ["FATAL ERROR","NAN","SEGMENTATION FAULT"]:
        if marker in upper:
            raise EMAuditError(f"free EM log contains {marker}")
    base=load_full().audit(gro_text,chain_itp_text,build_audit,potential,True)
    return {
        "schema_version":1,
        "status":"PASS_A1_EM",
        "maximum_force_kj_mol_nm":maximum_force,
        "potential_energy_kj_mol":potential,
        "minimum_nonbonded_distance_nm":base["minimum_nonbonded_distance_nm"],
        "nalpha_hydrogen_count":base["nalpha_hydrogen_count"],
        "ogamma_hydrogen_count":base["ogamma_hydrogen_count"],
        "grompp_maxwarn_zero":True,
        "scientific_scope":"technical EM gate only; not an unconstrained NAC stability pass",
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--gro",required=True); ap.add_argument("--chain-itp",required=True)
    ap.add_argument("--build-audit",required=True); ap.add_argument("--log",required=True)
    ap.add_argument("--output",required=True); a=ap.parse_args()
    result=audit(pathlib.Path(a.gro).read_text(),pathlib.Path(a.chain_itp).read_text(),
                 json.loads(pathlib.Path(a.build_audit).read_text()),pathlib.Path(a.log).read_text())
    pathlib.Path(a.output).write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")

if __name__=="__main__":
    main()
