"""Extract six existing PET motifs and native LCC fragments; no modeling or sampling."""
from pathlib import Path
import numpy as np,json,hashlib,datetime
from scipy.spatial.distance import cdist
P=Path("/data/bht2/polymer_material_reference_20260827/simulation_slabs/pet_dp10_400chain_direct_v1")
out=P/"motif_top6_vs_native_lcc_v1";out.mkdir(exist_ok=False)
base=P/"donor_geometry_batch_singleNH_v1/joint_transfer_map_v1"
files=[base/"joint_transfer_results.npz",base/"motifs.json",P/"catalytic_motif_compare_v1/4EB0.pdb",P/"catalytic_motif_compare_v1/5XJH.pdb",P/"native_motif_surface_compare_v1/native_3d_data.json"]
hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
assert hashes[str(files[2])]=="952288d4a49893ee259638bc72e17868b6903f0027c4bd4d207f76b9ddf85718"
assert hashes[str(files[3])]=="a0410e1407aa9ccf5fa16396a3fe329c46504fc4d6182fb438607ec23c697680"
z=np.load(files[0]);meta=json.loads(files[1].read_text());scene=json.loads(files[4].read_text())
ids=["M449","M346","M218","M318","M329","M450"]
def atoms(path):
 d={}
 for l in path.read_text().splitlines():
  if l.startswith("ATOM  ") and l[21]=="A" and l[16] in (" ","A") and l[76:78].strip()!="H":
   key=(int(l[22:26]),l[12:16].strip())
   d[key]=dict(name=key[1],resid=key[0],resn=l[17:20],chain="A",element=l[76:78].strip(),xyz=[float(l[30:38]),float(l[38:46]),float(l[46:54])])
 return d
lcc=atoms(files[2]);isp=atoms(files[3])
nkeys=[k for k in lcc if k[0] in [165,242,210]]
for d in [95,166]:
 for k in [(d-1,"C"),(d-1,"O"),(d,"N"),(d,"CA"),(d,"C"),(d,"O")]:
  if k not in nkeys:nkeys.append(k)
native=[lcc[k].copy() for k in nkeys];nq=np.array([a["xyz"] for a in native])
def nh(d,aa):
 n=np.array(aa[d,"N"]["xyz"]);a=np.array(aa[d-1,"C"]["xyz"])-n;b=np.array(aa[d,"CA"]["xyz"])-n
 v=-(a/np.linalg.norm(a)+b/np.linalg.norm(b));return n+1.01*v/np.linalg.norm(v)
nhs=np.array([nh(d,lcc) for d in [95,166]])
assert len(native)==34
def fitkeys(ser,his,asp):
 return [(ser,"CB"),(ser,"OG"),(his,"CG"),(his,"ND1"),(his,"CD2"),(his,"CE1"),(his,"NE2"),(asp,"CG"),(asp,"OD1"),(asp,"OD2")]
native_fit=[nkeys.index(k) for k in fitkeys(165,242,210)]
origin=np.array(lcc[165,"OG"]["xyz"]);ref=nq-origin
def bonds(records,q):
 b=[]
 for i in range(len(q)):
  for j in range(i+1,len(q)):
   a,c=records[i],records[j]
   if a["chain"]!=c["chain"]:continue
   same=a["resid"]==c["resid"]
   peptide=(a["name"]=="C" and c["name"]=="N" and c["resid"]==a["resid"]+1) or (c["name"]=="C" and a["name"]=="N" and a["resid"]==c["resid"]+1)
   if same or peptide:
    dd=np.linalg.norm(q[i]-q[j]);hydro=a["element"]=="H" or c["element"]=="H"
    if .75<dd<(1.15 if hydro else 1.9):b.append((i+1,j+1))
 return b
def writepdb(name,records,q,h,donor_res,remarks):
 rec=[a.copy() for a in records]
 for chain,resid in donor_res:
  orig=next(a for a in rec if a["chain"]==chain and a["resid"]==resid and a["name"]=="N")
  rec.append(dict(name="HN",resid=resid,resn=orig["resn"],chain=chain,element="H"))
 qq=np.concatenate([q,h])
 lines=["REMARK "+r+"\n" for r in remarks+["VIRTUAL NH HYDROGENS: geometrical 1.01 A bisector; not crystallographic H.","FRAGMENT ONLY; not full enzyme, not MD ready."]]
 for i,(a,v) in enumerate(zip(rec,qq),1):
  lines.append(f"ATOM  {i:5d} {a['name']:^4s} {a['resn']:>3s} {a['chain']}{a['resid']:4d}    {v[0]:8.3f}{v[1]:8.3f}{v[2]:8.3f}{1.:6.2f}{0.:6.2f}          {a['element']:>2s}\n")
 for i,j in bonds(rec,qq):lines.append(f"CONECT{i:5d}{j:5d}\n")
 lines.append("END\n");(out/name).write_text("".join(lines))
 parsed=np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in lines if l.startswith("ATOM")])
 assert len(parsed)==len(qq) and np.max(abs(parsed-qq))<=.000501
 return rec
def esterfile(name,e):
 lines=["REMARK Three target ester atoms only, a coordinate reference not complete PET.\n"]
 for i,(an,el,v) in enumerate(zip(["C","O","OE"],["C","O","O"],e),1):
  lines.append(f"HETATM{i:5d} {an:^4s} EST Z   1    {v[0]:8.3f}{v[1]:8.3f}{v[2]:8.3f}{1.:6.2f}{0.:6.2f}          {el:>2s}\n")
 lines+=["CONECT    1    2    3\n","END\n"];(out/name).write_text("".join(lines))
manifest={"status":"EXTRACTED_NO_MODELING","source_sha256":hashes,"native":"4EB0 chain A: Ser165 His242 Asp210; Tyr95 and Met166 backbone NH; 34 unique heavy atoms","candidate":"Each original 36-heavy-atom motif: own triad plus two independent copies of 5XJH Tyr87 peptide; not native donor placement","alignment":"10 catalytic sidechain atoms: Ser CB/OG, His CG/ND1/CD2/CE1/NE2, Asp CG/OD1/OD2. Proper rigid Kabsch fit, no minimization. Native reference Ser OG at origin.","items":[]}
writepdb("native_LCC_crystal_fragment.pdb",native,nq,nhs,[("A",95),("A",166)],["Native LCC 4EB0 A; original crystal coordinates."])
writepdb("native_LCC_triad_aligned.pdb",native,ref,nhs-origin,[("A",95),("A",166)],["Native LCC 4EB0 A; Ser165 OG translated to origin."])
natdata=next(d for d in scene["native"] if d["name"]=="LCC")
assert natdata["labels"]==[f"{a}:{b}" for a,b in nkeys]
example=next(d for d in scene["surface_examples"] if d["name"]=="LCC")
eq=np.array(example["q"]);eh=np.array(example["H"])
assert np.max(abs(cdist(eq,eq)-cdist(nq,nq)))<1e-7
writepdb("native_LCC_ester_frame.pdb",native,eq,eh,[("A",95),("A",166)],["Native LCC rigid motif, previously validated PET pose at "+example["site_id"],"This is an example placement, not a new optimization."])
esterfile("native_LCC_ester_reference.pdb",example["ester"])
for name in ids:
 i=next(j for j,m in enumerate(meta) if m["id"]==name);m=meta[i];q=z["motif_xyz_A"][i];h=z["motif_H_A"][i]
 nr={"SER":160,"HIS":237,"ASP":206} if m["template"]=="IsPETase" else {"SER":165,"HIS":242,"ASP":210}
 rec=[dict(name=an,resid=nr[r],resn=r,chain="A",element=el) for an,r,el in zip(m["atom_names"][:24],m["residues"][:24],m["elements"][:24])]
 for chain in ["B","C"]:
  for k in [(86,"C"),(86,"O"),(87,"N"),(87,"CA"),(87,"C"),(87,"O")]:
   a=isp[k].copy();a["chain"]=chain;rec.append(a)
 keys=[(a["resid"],a["name"]) for a in rec[:24]]
 fi=[keys.index(k) for k in fitkeys(nr["SER"],nr["HIS"],nr["ASP"])]
 A=q[fi];T=ref[native_fit];ac=A.mean(0);tc=T.mean(0)
 u,sv,vt=np.linalg.svd((A-ac).T@(T-tc));rot=u@vt
 if np.linalg.det(rot)<0:u[:,-1]*=-1;rot=u@vt
 aq=(q-ac)@rot+tc;ah=(h-ac)@rot+tc
 rmsd=float(np.sqrt(np.mean(np.sum((aq[fi]-T)**2,axis=1))))
 assert np.max(abs(cdist(aq,aq)-cdist(q,q)))<1e-8 and abs(np.linalg.det(rot)-1)<1e-8
 rem=[name+"; triad source "+m["template"]+"; sampled donor fragments, NOT native full motif.","NH1 chain B residue87; NH2 chain C residue87. Triad chain A."]
 writepdb(name+"_ester_frame.pdb",rec,q,h,[("B",87),("C",87)],rem+["Original ester registration, unchanged."])
 writepdb(name+"_triad_aligned.pdb",rec,aq,ah,[("B",87),("C",87)],rem+[f"Rigid 10-atom fit to native LCC; RMSD {rmsd:.6f} A."])
 esterfile(name+"_ester_reference.pdb",m["source_ester_local_A"])
 dn=aq[[26,32]];nn=ref[[nkeys.index((95,"N")),nkeys.index((166,"N"))]]
 dd=cdist(dn,nn);match=[0,1] if dd[0,0]**2+dd[1,1]**2<=dd[0,1]**2+dd[1,0]**2 else [1,0]
 manifest["items"].append(dict(id=name,triad_source=m["template"],center_coverage=int(z["coverage"][i].sum()),heavy_atoms=36,virtual_H=2,triad_fit_RMSD_A=rmsd,donor_N_distance_matrix_A=dd.tolist(),best_unordered_N_matching=match,donor_N_offsets_after_triad_fit_A=[float(dd[j,match[j]]) for j in range(2)],ester_N_A=q[[26,32]].tolist(),aligned_N_A=dn.tolist()))
(out/"manifest.json").write_text(json.dumps(manifest,indent=2))
(out/"RUNBOOK.md").write_text("# Top six candidate motifs versus native LCC\n\nExtraction only; no new sampling, fitting of internal bonds, or energetic minimization. Run saved export script on the CPU environment; output folder must not exist. Read manifest.json for exact atom identities, source hashes and alignment audit. Two coordinates: ester_frame preserves original substrate-registered placement; triad_aligned rigidly fits catalytic sidechain atoms to native LCC. Native ester-frame structure is one existing PET-compatible pose, not the native crystal substrate complex. Each ester reference is only three atoms, not full PET. Artificial NH fragments are chain B/C; native own NH donors are chain A Tyr95 and Met166. Two virtual H atoms are included. No duplicate atoms in native motif.\n")
for p,sha in hashes.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
(out/"RUN_LOG.jsonl").write_text(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=__file__,status=manifest["status"],source_inputs_unchanged=True,candidates=ids))+"\n")
(out/"SHA256.json").write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()},indent=2))
print(json.dumps({"output":str(out),"manifest":manifest,"file_count":len(list(out.iterdir()))}),flush=True)
