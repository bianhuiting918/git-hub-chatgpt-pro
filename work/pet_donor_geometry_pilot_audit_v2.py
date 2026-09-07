from pathlib import Path
import json,hashlib,datetime
import numpy as np
from scipy.spatial.distance import cdist
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
P=Path("/data/bht2/polymer_material_reference_20260827/simulation_slabs/pet_dp10_400chain_direct_v1")
out=P/"donor_geometry_pilot_v1"
folder=out/"P244_C42_O18_IsPETase"
r=json.loads((folder/"summary.json").read_text());a=np.load(folder/"conditional_space.npz")
k=r["example_indices"]["pose"];ia=r["example_indices"]["donor_A"];ib=r["example_indices"]["donor_B"]
tri=a["triads_xyz"][k];dA=a["A_q"][ia];dB=a["B_q"][ib]
rad={"C":1.70,"O":1.52,"N":1.55,"S":1.80}
pdb=P/"pet_dp10_400chain/water_slab_v1/dry_export/PET_DP10_400chain_water_equilibrated_dry.pdb"
rec=[l for l in pdb.read_text().splitlines() if l.startswith(("ATOM  ","HETATM")) and l[76:78].strip() in rad]
q=np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in rec]);pe=np.array([rad[l[76:78].strip()] for l in rec])
center=a["ester_xyz"][0];box=np.array([125.697,125.697]);q[:,:2]-=np.round((q[:,:2]-center[:2])/box)*box
def radii(path):
 return np.array([rad[l[76:78].strip()] for l in path.read_text().splitlines() if l.startswith(("ATOM  ","HETATM"))])
tr=radii(folder/"triad.pdb");ar=radii(folder/"donor_A.pdb");br=radii(folder/"donor_B.pdb")
def margin(x,xr,y,yr):return float(np.min(cdist(x,y)-xr[:,None]-yr[None,:]+.4))
checks={"triad_PET_clearance_A":margin(tri,tr,q,pe),"A_PET_clearance_A":margin(dA,ar,q,pe),"B_PET_clearance_A":margin(dB,br,q,pe),"A_triad_clearance_A":margin(dA,ar,tri,tr),"B_triad_clearance_A":margin(dB,br,tri,tr),"A_B_clearance_A":margin(dA,ar,dB,br)}
assert min(checks.values())>=-1e-8
for tag,idx in [("A",ia),("B",ib)]:
 N=a[tag+"_q"][idx,2];H=a[tag+"_H"][idx];O=a["ester_xyz"][1]
 u=N-H;v=O-H;angle=float(np.degrees(np.arccos(np.clip(np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v),-1,1))))
 dist=float(np.linalg.norm(N-O));assert 2.7<=dist<=3.2 and angle>=140
 checks[tag+"_N_O_A"]=dist;checks[tag+"_N_H_O_deg"]=angle
# Rotating coordinates must preserve directly computed clash margins.
rot=np.array([[0.,-1,0],[1,0,0],[0,0,1]])
assert abs(margin(dA@rot,ar,dB@rot,br)-checks["A_B_clearance_A"])<1e-10
info={"status":"INDEPENDENT_DIRECT_DISTANCE_REPLAY_PASS","checks":checks,"scope":"one representative successful pair; complete original PET heavy atoms with XY minimum images; declared 0.4 A overlap allowance"}
(folder/"independent_audit.json").write_text(json.dumps(info,indent=2))
cloud={}
for tag in ("A","B"):
 ids=np.flatnonzero(a["mask_"+tag][k])
 cloud[tag]={"N":a[tag+"_q"][ids,2].tolist(),"H":a[tag+"_H"][ids].tolist(),"count":len(ids)}
cloud["example_H"]={"A":a["A_H"][ia].tolist(),"B":a["B_H"][ib].tolist()}
cloud["carbonyl_O_name"]=[l[12:16].strip() for l in (folder/"target_ester.pdb").read_text().splitlines() if l.startswith(("ATOM  ","HETATM"))][1]
cloud["ester_C"]=a["ester_xyz"][0].tolist();cloud["ester_O"]=a["ester_xyz"][1].tolist()
(folder/"donor_cloud.json").write_text(json.dumps(cloud))
fig,axes=plt.subplots(1,2,figsize=(11,4.8),constrained_layout=True)
for ax,tag,color in zip(axes,["A","B"],["teal","darkorange"]):
 d=a[tag+"_direction"];az=np.degrees(np.arctan2(d[:,2],d[:,1]));mask=a["mask_"+tag][k]
 ax.scatter(az, d[:,0],s=14,c="lightgray",label="PET-clear, before triad")
 ax.scatter(az[mask],d[mask,0],s=22,c=color,label="PET + selected triad clear")
 for pair in a["pair_examples"]:
  if pair[0]==k:
   idx=pair[1 if tag=="A" else 2]
   ax.scatter(az[idx],d[idx,0],s=95,facecolors="none",edgecolors="red")
 ax.set(xlim=(-180,180),ylim=(-1,1),xlabel="Azimuth around C=O (deg)",ylabel="cos(polar angle to C->O)",title="Donor "+tag+" | "+str(mask.sum())+" feasible fragments")
 ax.legend(fontsize=8)
fig.suptitle("P244_C42_O18 | IsPETase triad pose 5\nRed rings: found mutually compatible pair samples; not probability or activity")
fig.savefig(folder/"conditional_donor_directions.png",dpi=150)
local="D:/Codex/polymer_material_reference_20260827/simulation_slabs/pet_dp10_400chain_direct_v1/pymol/donor_geometry_pilot_v1"
viewer="""from pathlib import Path
import json
from pymol import cmd
from pymol.cgo import BEGIN,END,LINES,COLOR,VERTEX,SPHERE
ROOT=Path(LOCAL)
cloud=json.loads((ROOT/'donor_cloud.json').read_text())
cmd.disable('all')
for filename,obj,color in [('PET_context_20A.pdb','pilot_PET','gray70'),('triad.pdb','pilot_triad','forest'),('donor_A.pdb','pilot_donor_A','cyan'),('donor_B.pdb','pilot_donor_B','orange'),('target_ester.pdb','pilot_ester','yellow')]:
 cmd.delete(obj);cmd.load(str(ROOT/filename),obj);cmd.hide('everything',obj)
 cmd.color(color,obj)
 if obj=='pilot_PET':
  cmd.show('spheres',obj);cmd.set('sphere_scale',1,obj);cmd.set('sphere_transparency',.7,obj)
 else:
  cmd.show('sticks',obj);cmd.show('spheres',obj);cmd.set('sphere_scale',.2,obj)
  cmd.color('red',obj+' and elem O');cmd.color('blue',obj+' and elem N')
for tag,color in [('A',(0.,.7,.8)),('B',(1.,.5,0.))]:
 cgo=[]
 for n,h in zip(cloud[tag]['N'],cloud[tag]['H']):
  cgo.extend([COLOR,*color,SPHERE,*n,.12,BEGIN,LINES,VERTEX,*n,VERTEX,*h,END])
 name='pilot_donor_cloud_'+tag;cmd.delete(name)
 if cgo:cmd.load_cgo(cgo,name)
 cmd.delete('pilot_H_'+tag);cmd.pseudoatom('pilot_H_'+tag,pos=cloud['example_H'][tag],elem='H',color='white')
 cmd.show('spheres','pilot_H_'+tag);cmd.set('sphere_scale',.2,'pilot_H_'+tag)
 cmd.distance('pilot_Hbond_'+tag,'pilot_donor_'+tag+' and name N','pilot_ester and name '+cloud['carbonyl_O_name'])
cmd.bg_color('white');cmd.set('orthoscopic',1);cmd.orient('pilot_triad or pilot_donor_A or pilot_donor_B or pilot_ester');cmd.zoom('pilot_triad or pilot_donor_A or pilot_donor_B or pilot_ester',5)
print('Green=triad; cyan/orange=independent donor fragments. Small dots and lines are sampled N positions and N-H directions for this fixed triad pose. Gray spheres=PET steric context. These are geometry candidates, NOT validated catalysts.')
""".replace("LOCAL",repr(local))
(folder/"show_donor_pilot.py").write_text(viewer)
with (out/"RUN_LOG.jsonl").open("a") as f:f.write(json.dumps({"time":datetime.datetime.now(datetime.timezone.utc).isoformat(),"postprocessing":info})+"\n")
print(json.dumps(info),flush=True)
print("CLOUD_COUNTS",cloud["A"]["count"],cloud["B"]["count"],flush=True)
