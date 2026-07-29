#!/usr/bin/env python3
"""Run one auditable whole-receptor A/C/OA/HD surface-field task."""
from __future__ import annotations
import argparse,csv,hashlib,json,subprocess,sys
from datetime import datetime,timezone
from pathlib import Path

def sha256_file(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1048576),b""): h.update(chunk)
    return h.hexdigest()

def load_manifest_row(path,zero_based_index):
    if zero_based_index<0: raise IndexError("manifest row index must be non-negative")
    with Path(path).open(newline="",encoding="utf-8") as f:
        for i,row in enumerate(csv.DictReader(f,delimiter="\t")):
            if i==zero_based_index: return {k:(v or "").strip() for k,v in row.items()}
    raise IndexError(f"manifest row index out of range: {zero_based_index}")

def validate_manifest_row(row,family):
    if row.get("material_family","").upper()!=family.upper(): return "NOT_EVALUATED_FAMILY_MISMATCH"
    if row.get("input_status")!="READY_FOR_STRUCTURE_EVALUATION": return "NOT_EVALUATED_INPUT_STATUS"
    if not row.get("candidate_id") or len(row.get("sequence_md5",""))!=32: return "NOT_EVALUATED_MISSING_IDENTITY"
    if not row.get("selected_chain"): return "NOT_EVALUATED_SELECTED_CHAIN_MISSING"
    receptor=Path(row.get("receptor_path",""))
    if not receptor.is_file(): return "NOT_EVALUATED_RECEPTOR_MISSING"
    expected=row.get("receptor_sha256","").lower()
    if len(expected)!=64: return "NOT_EVALUATED_RECEPTOR_CHECKSUM_MISSING"
    try: int(expected,16)
    except ValueError: return "NOT_EVALUATED_RECEPTOR_CHECKSUM_INVALID"
    if sha256_file(receptor)!=expected: return "NOT_EVALUATED_RECEPTOR_CHECKSUM_MISMATCH"
    return ""

def require_under_root(path,root):
    p,r=Path(path).resolve(),Path(root).resolve()
    if p!=r and r not in p.parents: raise ValueError(f"output outside project root: {p}")
    return p

def existing_status(target):
    target=Path(target); gate=target/"FIELD_PASS.json"
    if gate.is_file():
        try:
            if json.loads(gate.read_text(encoding="utf-8")).get("status")=="FIELD_PASS": return "PASS"
        except (OSError,json.JSONDecodeError): pass
    return "PARTIAL" if target.exists() and any(target.iterdir()) else "ABSENT"

def write_json_exclusive(path,payload):
    with Path(path).open("x",encoding="utf-8") as f:
        json.dump(payload,f,indent=2,sort_keys=True); f.write("\n")

def run_logged(command,out,err):
    with Path(out).open("x",encoding="utf-8") as o,Path(err).open("x",encoding="utf-8") as e:
        subprocess.run(command,check=True,stdout=o,stderr=e,text=True)

def parse_args():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",type=Path,required=True); p.add_argument("--row-index",type=int,required=True)
    p.add_argument("--family",choices=("PET","NYLON"),required=True); p.add_argument("--output-root",type=Path,required=True)
    p.add_argument("--project-root",type=Path,required=True); p.add_argument("--tool-root",type=Path,required=True)
    p.add_argument("--adfr-root",type=Path,required=True)
    return p.parse_args()

def main():
    a=parse_args(); output_root=require_under_root(a.output_root,a.project_root)
    row=load_manifest_row(a.manifest,a.row_index); md5=row.get("sequence_md5") or f"ROW_{a.row_index:06d}"
    target=require_under_root(output_root/md5,a.project_root); state=existing_status(target)
    if state=="PASS": print(f"existing FIELD_PASS retained: {target}"); return 0
    if state=="PARTIAL": print(f"refusing partial output: {target}",file=sys.stderr); return 3
    target.mkdir(parents=True); started=datetime.now(timezone.utc).isoformat()
    reason=validate_manifest_row(row,a.family)
    if reason:
        write_json_exclusive(target/"NOT_EVALUATED.json",{"status":reason,"scientific_scope":"technical_input_validation_only","manifest":str(a.manifest.resolve()),"row_index":a.row_index,"candidate_id":row.get("candidate_id",""),"sequence_md5":row.get("sequence_md5",""),"started_at_utc":started})
        return 4
    rec=target/"receptor"; rec.mkdir(); selected=rec/"selected.pdb"; selmeta=rec/"selection_metadata.json"
    pdbqt=rec/"receptor.pdbqt"; recmeta=rec/"metadata.json"; plan=target/"grid_plan.json"; maps=target/"maps"
    py=sys.executable; scripts=a.tool_root/"scripts"; binroot=a.adfr_root/"bin"; rid=f"{a.family}_{row['sequence_md5']}"
    commands=[
      ([py,str(scripts/"check_receptor.py"),"select","--input",row["receptor_path"],"--output",str(selected),"--metadata",str(selmeta),"--chains",row["selected_chain"],"--model","1","--catalytic-residues","","--retain-hetero",""],"01_select"),
      ([str(binroot/"prepare_receptor"),"-r",str(selected),"-o",str(pdbqt),"-A","checkhydrogens","-U","nphs_lps_waters"],"02_prepare"),
      ([py,str(scripts/"check_receptor.py"),"finalize","--selection-metadata",str(selmeta),"--prepared-pdbqt",str(pdbqt),"--output-metadata",str(recmeta)],"03_finalize"),
      ([py,str(scripts/"build_grid_plan.py"),"--receptor-metadata",str(recmeta),"--output",str(plan)],"04_grid_plan"),
      ([py,str(scripts/"run_grid_maps.py"),"--receptor-pdbqt",str(pdbqt),"--grid-plan",str(plan),"--output-dir",str(maps),"--record-id",rid,"--autogrid",str(binroot/"autogrid4"),"--autosite",str(binroot/"autosite"),"--parameter-file",str(a.adfr_root/"CCSBpckgs"/"AutoDockTools"/"AD4.1_bound.dat"),"--project-root",str(a.project_root)],"05_maps")]
    try:
        for command,label in commands: run_logged(command,target/f"{label}.stdout.log",target/f"{label}.stderr.log")
        gate=maps/"MAPS_PASS.json"
        payload=json.loads(gate.read_text(encoding="utf-8"))
        if payload.get("status")!="MAPS_PASS": raise RuntimeError("map gate is not MAPS_PASS")
        write_json_exclusive(target/"FIELD_PASS.json",{"status":"FIELD_PASS","scientific_scope":"whole_receptor_ACOAHD_field_generation_only","material_family":a.family,"manifest":str(a.manifest.resolve()),"manifest_sha256":sha256_file(a.manifest),"row_index":a.row_index,"candidate_id":row["candidate_id"],"sequence_md5":row["sequence_md5"],"selected_chain":row["selected_chain"],"receptor_path":row["receptor_path"],"receptor_sha256":row["receptor_sha256"],"maps_pass_sha256":sha256_file(gate),"channels":payload.get("channels",[]),"started_at_utc":started,"completed_at_utc":datetime.now(timezone.utc).isoformat()})
    except Exception as error:
        write_json_exclusive(target/"NOT_EVALUATED.json",{"status":"NOT_EVALUATED_SURFACE_FIELD_TECHNICAL","scientific_scope":"technical_execution_failure_only","manifest":str(a.manifest.resolve()),"row_index":a.row_index,"candidate_id":row.get("candidate_id",""),"sequence_md5":row.get("sequence_md5",""),"error_type":type(error).__name__,"error":str(error),"started_at_utc":started,"failed_at_utc":datetime.now(timezone.utc).isoformat()})
        raise
    return 0
if __name__=="__main__": raise SystemExit(main())
