#!/usr/bin/env python3
"""PET surface chemical-proxy and 15 A shape statistics using PyMOL SES.
CPU-only. All chains retained. Not hydration free energy or enzyme activity.
--project ORIGINAL_PROJECT --output NEW_RESULTS --pymol-python ISOLATED_PYTHON
"""
import argparse, csv, hashlib, json, math, os, sys, time, subprocess, traceback
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
CATEGORIES=["aromatic_CH","aliphatic_CH","carbonyl_C","carbonyl_O","ester_O"]
COLORS=np.array([[.72,.52,.22],[.92,.80,.47],[.60,.45,.72],[.82,.12,.15],[1.,.48,.32]])
RADII={"C":1.70,"O":1.52,"H":1.20}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def now():return datetime.now(timezone.utc).isoformat()
def dump(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")
def table(p,rows):
    if rows:
        with Path(p).open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def weighted_quantile(x,w,qs):
    order=np.argsort(x);xx=x[order];ww=w[order]
    return np.interp(qs,(np.cumsum(ww)-.5*ww)/ww.sum(),xx)
def parse_obj(s):
    verts=[];faces=[]
    for line in s.splitlines():
        if line.startswith("v "):verts.append(list(map(float,line.split()[1:4])))
        elif line.startswith("f "):
            ids=[int(t.split("/")[0])-1 for t in line.split()[1:]]
            assert len(ids)==3,"Unexpected nontriangle OBJ"
            faces.append(ids)
    v=np.array(verts,dtype=np.float64);f=np.array(faces,dtype=np.int32)
    assert v.ndim==2 and f.ndim==2 and f.min()>=0 and f.max()<len(v)
    return v,f
def mesh_geometry(v,f):
    tri=v[f];cross=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    area=np.linalg.norm(cross,axis=1)/2
    return tri.mean(1),area,cross/np.maximum(2*area[:,None],1e-30)

def export_mesh(inp,out,quality,selftest=False):
    import pymol2
    from chempy import Atom
    from chempy.models import Indexed
    with pymol2.PyMOL() as pm:
        c=pm.cmd
        c.set("max_threads",8)
        c.set("surface_solvent",0);c.set("solvent_radius",1.4)
        c.set("surface_quality",quality);c.set("surface_mode",1)
        c.set("surface_cavity_mode",0)
        c.set("auto_zoom",0)
        if selftest:
            c.pseudoatom("pet",pos=[11.,22.,33.],vdw=2.)
        else:
            data=np.load(inp);model=Indexed()
            for pos,r,k in zip(data["halo_xyz"],data["halo_vdw"],data["halo_category"]):
                a=Atom();a.coord=pos.tolist();a.vdw=float(r);a.symbol="C";a.name="C";a.resn="PET";a.b=float(k);model.add_atom(a)
            c.load_model(model,"pet")
        c.hide("everything");c.show("surface","pet")
        for k,rgb in enumerate(COLORS):
            c.set_color("chem"+str(k),rgb.tolist());c.color("chem"+str(k),f"pet and b > {k-.1} and b < {k+.1}")
        print("PYMOL_SURFACE_START",c.count_atoms(),quality,flush=True)
        mtl,obj=c.get_mtl_obj()
        v,f=parse_obj(obj);cent,areas,norm=mesh_geometry(v,f)
        if selftest:
            centre=(v.min(0)+v.max(0))/2
            assert np.linalg.norm(centre-[11,22,33])<.1,("Export coordinate transform",centre)
            err=abs(areas.sum()/(16*np.pi)-1)
            assert err<.1,("Single sphere SES area",areas.sum())
            print("PYMOL_SELFTEST_PASS",json.dumps({"version":c.get_version()[0],"vertices":len(v),"triangles":len(f),"sphere_area_A2":float(areas.sum()),"relative_error":float(err)}),flush=True)
            return
        np.savez_compressed(out,vertices_A=v.astype(np.float32),faces=f)
        metadata={"pymol_version":c.get_version()[0],"atoms_with_xy_halo":c.count_atoms(),"quality":quality,
             "settings":{s:c.get(s) for s in ["surface_solvent","solvent_radius","surface_quality","surface_mode","surface_cavity_mode","max_threads"]},
             "vertices":len(v),"faces":len(f),"coordinate_units":"Angstrom"}
        dump(str(out)+".json",metadata)
        print("PYMOL_SURFACE_DONE",json.dumps(metadata),flush=True)
def parse_inputs(project,probe):
    pdb=project/"pet_dp10_400chain/water_slab_v1/dry_export/PET_DP10_400chain_water_equilibrated_dry.pdb"
    sdf=project/"00_chain/PET_DP10.sdf"
    itp=project/"pet_dp10_400chain/parameterization/am1bcc_singlepoint_v1/PET_DP10.acpype/PET_DP10_GMX.itp"
    lines=pdb.read_text().splitlines()
    c=next(l for l in lines if l.startswith("CRYST1"))
    box=np.array([float(c[6:15]),float(c[15:24]),float(c[24:33])])
    angles=[float(c[33:40]),float(c[40:47]),float(c[47:54])]
    assert np.allclose(angles,90),"Only orthorhombic input supported"
    atomlines=[l for l in lines if l.startswith(("ATOM  ","HETATM"))]
    xyz=np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in atomlines])
    names=np.array([l[12:16].strip() for l in atomlines])
    segs=np.array([l[72:76].strip() for l in atomlines])
    elems=np.array([l[76:78].strip() for l in atomlines])
    mol=Chem.SDMolSupplier(str(sdf),removeHs=False)[0]
    assert mol is not None
    n=mol.GetNumAtoms();base_e=np.array([a.GetSymbol() for a in mol.GetAtoms()])
    section="";topatoms=[];bonds=set()
    for line in itp.read_text().splitlines():
        line=line.split(";")[0].strip()
        if line.startswith("["):section=line.strip("[] ");continue
        if not line:continue
        t=line.split()
        if section=="atoms":topatoms.append(t)
        elif section=="bonds":bonds.add(tuple(sorted([int(t[0])-1,int(t[1])-1])))
    sbonds={tuple(sorted([b.GetBeginAtomIdx(),b.GetEndAtomIdx()])) for b in mol.GetBonds()}
    assert sbonds==bonds,"SDF vs ITP connectivity mismatch"
    expected=np.array([a[4] for a in topatoms])
    assert len(expected)==n
    segments=list(dict.fromkeys(segs));groups=[np.flatnonzero(segs==s) for s in segments]
    for g in groups:
        assert len(g)==n and np.array_equal(names[g],expected)
        assert np.array_equal(elems[g],base_e)
    graph=[[] for _ in range(n)]
    for i,j in bonds:graph[i].append(j);graph[j].append(i)
    maxbond=0
    # Bond-graph unwrapping then place each chain z-centre near box centre.
    for g in groups:
        raw=xyz[g].copy();un=np.zeros_like(raw);un[0]=raw[0];seen={0};stack=[0]
        while stack:
            i=stack.pop()
            for j in graph[i]:
                d=raw[j]-raw[i];d-=box*np.round(d/box)
                maxbond=max(maxbond,float(np.linalg.norm(d)))
                if j not in seen:un[j]=un[i]+d;seen.add(j);stack.append(j)
        assert len(seen)==n
        un[:,2]-=box[2]*np.round((un[:,2].mean()-box[2]/2)/box[2])
        xyz[g]=un
    assert maxbond<2.0,"Unexpected bond lengths or bad atom mapping"
    xyz[:,:2]%=box[:2]
    assert xyz[:,2].min()>probe+RADII["C"] and xyz[:,2].max()<box[2]-probe-RADII["C"],"Slab intersects z cell edge"
    bcat=[]
    for a in mol.GetAtoms():
        t=topatoms[a.GetIdx()][1]
        if a.GetSymbol()=="H":
            parent=next(iter(a.GetNeighbors()));t=topatoms[parent.GetIdx()][1]
        mapping={"ca":0,"c3":1,"c":2,"o":3,"os":4}
        assert t in mapping,("Unexpected chemical type",t)
        bcat.append(mapping[t])
    cat=np.zeros(len(xyz),dtype=np.int8)
    for g in groups:cat[g]=bcat
    esters=[];caps=[]
    for ci,co,eo in mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")):
        aromatic=any(a.GetIsAromatic() for a in mol.GetAtomWithIdx(ci).GetNeighbors())
        (esters if aromatic else caps).append((ci,co,eo))
    assert len(esters)==20 and len(caps)==1 and len(groups)==400 and n==231
    sites=[]
    for seg,g in zip(segments,groups):
        for k,(ci,co,eo) in enumerate(esters,1):
            sites.append({"site_id":f"{seg}_{names[g[ci]]}_{names[g[eo]]}",
                          "segment":seg,"ester_index_in_chain":k,
                          "carbonyl_C_name":str(names[g[ci]]),"carbonyl_O_name":str(names[g[co]]),
                          "ester_O_name":str(names[g[eo]]),
                          "atom_indices":[int(g[ci]),int(g[co]),int(g[eo])]})
    radii=np.array([RADII[e]+probe for e in elems])
    audit={"pdb":str(pdb),"pdb_sha256":sha(pdb),"sdf_sha256":sha(sdf),"itp_sha256":sha(itp),
           "atoms":len(xyz),"chains":len(groups),"atoms_per_chain":n,"backbone_esters":len(sites),
           "excluded_artificial_cap_esters":len(groups)*len(caps),"box_A":box.tolist(),
           "max_minimum_image_bond_A":maxbond,"z_extent_after_unwrap_A":[float(xyz[:,2].min()),float(xyz[:,2].max())],
           "atom_name_element_and_bond_mapping":"PASS"}
    return xyz,radii,box,cat,segs,sites,audit

def halo_input(xyz,vdw,cat,box,margin=8.):
    positions=[];radii=[];categories=[];ids=[]
    for i in [-1,0,1]:
        for j in [-1,0,1]:
            pos=xyz+np.array([i*box[0],j*box[1],0])
            take=np.all((pos[:,:2]>=-margin)&(pos[:,:2]<=box[:2]+margin),axis=1)
            positions.append(pos[take]);radii.append(vdw[take]);categories.append(cat[take]);ids.append(np.flatnonzero(take))
    return {"halo_xyz":np.concatenate(positions),"halo_vdw":np.concatenate(radii),
        "halo_category":np.concatenate(categories),"halo_original_atom_index":np.concatenate(ids).astype(np.int32)}



def exterior_mask(v,faces):
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree
    unique,inv=np.unique(np.round(v,4),axis=0,return_inverse=True)
    ff=inv[faces];nv=len(unique)
    row=np.concatenate([ff[:,0],ff[:,1]]);col=np.concatenate([ff[:,1],ff[:,2]])
    graph=coo_matrix((np.ones(len(row),np.int8),(row,col)),shape=(nv,nv)).tocsr()
    nc,labels=connected_components(graph,directed=False);lab=labels[ff[:,0]]
    tri=v[faces];cent=tri.mean(1);area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    areas=np.bincount(lab,weights=area,minlength=nc)
    components=[np.flatnonzero(lab==i) for i in range(nc)]
    shells=[];retained=[];excluded=[];unresolved=[];checks=[]
    def make_shell(i):
        tt=tri[components[i]];cc=cent[components[i]]
        radius=float(np.linalg.norm(tt-cc[:,None,:],axis=2).max())+1e-6
        return {"id":i,"tri":tt,"trees":[cKDTree(cc[:,[j for j in range(3) if j!=ax]]) for ax in range(3)],
            "min":tt.min((0,1)),"max":tt.max((0,1)),"radius":radius}
    def inside(shell,p):
        if np.any(p<shell["min"]) or np.any(p>shell["max"]):return False,[0,0,0]
        parities=[]
        for ax in range(3):
            xy=[j for j in range(3) if j!=ax]
            hit=shell["trees"][ax].query_ball_point(p[xy],shell["radius"])
            t=shell["tri"][hit]
            a=t[:,0][:,xy];b=t[:,1][:,xy]-a;c=t[:,2][:,xy]-a;d=p[xy]-a
            det=b[:,0]*c[:,1]-b[:,1]*c[:,0];ok=np.abs(det)>1e-12
            t=t[ok];b=b[ok];c=c[ok];d=d[ok];det=det[ok]
            u=(d[:,0]*c[:,1]-d[:,1]*c[:,0])/det
            w=(b[:,0]*d[:,1]-b[:,1]*d[:,0])/det
            yes=(u>=-1e-10)&(w>=-1e-10)&(u+w<=1+1e-10)
            z=t[:,0,ax]+u*(t[:,1,ax]-t[:,0,ax])+w*(t[:,2,ax]-t[:,0,ax])-p[ax]
            hits=np.unique(np.round(z[yes&(z>1e-6)],4))
            parities.append(int(len(hits)%2))
        if len(set(parities))!=1:return None,parities
        return bool(parities[0]),parities
    for i in np.argsort(areas)[::-1]:
        i=int(i);cids=components[i]
        # Component test point is on this shell, not on the enclosing shell.
        p=cent[cids[np.argmax(area[cids])]]
        parent=None;ambiguous=False
        for sh in shells:
            within,parities=inside(sh,p);checks.append({"component":i,"against":sh["id"],"parity_xyz":parities})
            if within is None:ambiguous=True;break
            if within:parent=sh["id"];break
        if ambiguous:unresolved.append(i)
        elif parent is None:
            retained.append(i);shells.append(make_shell(i))
        else:excluded.append(i)
    keep=np.isin(lab,retained)
    audit={"components":nc,"external_components":retained,"excluded_nested_components":excluded,
        "unresolved_components":unresolved,"unresolved_face_indices":np.flatnonzero(np.isin(lab,unresolved)).tolist(),
        "merged_vertex_precision_A":1e-4,"component_area_A2":areas.tolist(),"containment_checks":checks,
        "kept_mesh_area_A2":float(area[keep].sum()),"excluded_cavity_mesh_area_A2":float(area[np.isin(lab,excluded)].sum()),
        "unresolved_mesh_area_A2":float(area[np.isin(lab,unresolved)].sum()),
        "rule":"Weld mesh vertices, find connected shells, retain non-nested external shells. Test nesting by agreement of +x,+y,+z ray parity; no reliance on inconsistent OBJ triangle winding.",
        "containment_gate":"PASS_ALL_THREE_AXES_AGREE" if not unresolved else "NOT_EVALUATED_AMBIGUOUS_COMPONENTS"}
    return keep,audit

def exterior_selftest():
    from scipy.spatial import ConvexHull
    base=np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],float)
    f=ConvexHull(base).simplices
    for i,t in enumerate(f):
        if np.dot(np.cross(base[t[1]]-base[t[0]],base[t[2]]-base[t[0]]),base[t].mean(0))<0:f[i]=t[::-1]
    v=np.vstack([base*4+[10,10,10],base+[10,10,10],base+[30,10,10]])
    faces=np.vstack([f,f[:,::-1]+6,f+12])
    faces[0]=faces[0,::-1] # Export winding inconsistency must not alter containment.
    keep,a=exterior_mask(v,faces)
    assert keep[:8].all() and not keep[8:16].any() and keep[16:].all()
    translated,_=exterior_mask(v+[1000,-2000,3000],faces)
    assert np.array_equal(keep,translated)
    print("EXTERIOR_SELFTEST_PASS: cavity removed; external disconnected particle retained; winding and translation invariant",flush=True)

def process_mesh(mesh,halo,xyz,cat,segs,sites,box,out,save):
    from scipy.spatial import cKDTree
    dat=np.load(mesh);v=dat["vertices_A"].astype(float);faces=dat["faces"]
    cent,area,norm=mesh_geometry(v,faces)
    external,exterior_audit=exterior_mask(v,faces)
    ufull=np.array(exterior_audit.pop("unresolved_face_indices"),dtype=np.int64)
    uvalid=(area[ufull]>1e-10)&np.all((cent[ufull,:2]>=0)&(cent[ufull,:2]<box[:2]),axis=1)
    ufull=ufull[uvalid];up=cent[ufull];uw=area[ufull]
    np.savez_compressed(out/(Path(mesh).stem+"_unresolved_faces.npz"),full_mesh_face_indices=ufull,
        centroids_A=up.astype(np.float32),areas_A2=uw.astype(np.float32))
    dump(out/(Path(mesh).stem+"_exterior_audit.json"),exterior_audit)
    keep=external&(area>1e-10)&np.all((cent[:,:2]>=0)&(cent[:,:2]<box[:2]),axis=1)
    ids=np.flatnonzero(keep);cent=cent[keep];area=area[keep];norm=norm[keep]
    tree=cKDTree(halo["halo_xyz"])
    owners=np.empty(len(cent),np.int32)
    # Face-centroid chemistry is assigned to nearest VDW sphere surface.
    # This is an explicit chemical mapping rule, not a water-affinity potential.
    for start in range(0,len(cent),100000):
        block=cent[start:start+100000];dist,ix=tree.query(block,k=16,workers=8)
        owner=ix[np.arange(len(block)),np.argmin(np.abs(dist-halo["halo_vdw"][ix]),axis=1)]
        owners[start:start+len(block)]=halo["halo_original_atom_index"][owner]
    uowners=np.zeros(len(up),np.int32)
    if len(up):
        ud,ui=tree.query(up,k=16,workers=8)
        uh=ui[np.arange(len(up)),np.argmin(np.abs(ud-halo["halo_vdw"][ui]),axis=1)]
        uowners=halo["halo_original_atom_index"][uh]
    chemical=cat[owners]
    mid=float(np.median(xyz[:,2]))
    sidecode=(cent[:,2]<mid).astype(np.int8)
    summaries={};rows_by_side={};atom_areas={};uncertain_atom_areas={}
    for code,side in enumerate(["top","bottom"]):
        mask=sidecode==code;pts=cent[mask];w=area[mask];nn=norm[mask];ow=owners[mask];cc=chemical[mask];meshids=ids[mask]
        atomarea=np.bincount(ow,weights=w,minlength=len(xyz));atom_areas[side]=atomarea
        um=(up[:,2]>=mid) if side=="top" else (up[:,2]<mid)
        upts=up[um];uww=uw[um];uow=uowners[um]
        ua=np.bincount(uow,weights=uww,minlength=len(xyz));uncertain_atom_areas[side]=ua
        utree=cKDTree(np.vstack([upts+np.array([i*box[0],j*box[1],0]) for i in [-1,0,1] for j in [-1,0,1]])) if len(upts) else None
        # Periodic XY point tree; z is nonperiodic and never queried across water.
        treexyz=np.vstack([pts+np.array([i*box[0],j*box[1],0]) for i in [-1,0,1] for j in [-1,0,1]])
        t=cKDTree(treexyz)
        rows=[];patchids=[];offsets=[0]
        for site in sites:
            atomids=site["atom_indices"];exposed=float(atomarea[atomids].sum())
            if exposed<=0:continue
            origin=xyz[atomids[0]]
            hits=np.array(t.query_ball_point(origin,15.),dtype=np.int64)
            local=np.unique(hits%len(pts))
            if len(local)<6:raise ValueError("Insufficient surface faces for patch "+site["site_id"])
            delta=pts[local]-origin;delta[:,:2]-=box[:2]*np.round(delta[:,:2]/box[:2])
            assert np.max(np.linalg.norm(delta,axis=1))<=15.00001
            ww=w[local];pa=float(ww.sum());mean=np.average(delta,weights=ww,axis=0);dr=delta-mean
            cov=(dr.T@(ww[:,None]*dr))/pa
            eig,axes=np.linalg.eigh(cov);normal=axes[:,0]
            if normal[2]*(1 if side=="top" else -1)<0:normal=-normal
            frac=np.bincount(cc[local],weights=ww,minlength=5)/pa
            carbonyl=xyz[atomids[1]]-origin;carbonyl[:2]-=box[:2]*np.round(carbonyl[:2]/box[:2]);carbonyl/=np.linalg.norm(carbonyl)
            row={k:v for k,v in site.items() if k!="atom_indices"}
            unresolved_patch_area=0.
            if utree is not None:
                uhit=np.array(utree.query_ball_point(origin,15.),dtype=np.int64)
                uloc=np.unique(uhit%len(upts))
                unresolved_patch_area=float(uww[uloc].sum())
            row.update({"shape_scope_status":"COMPLETE_FOR_CLASSIFIED_MESH" if unresolved_patch_area==0 else "NOT_EVALUATED_COMPLETE_SHAPE",
                "unresolved_patch_area_A2":unresolved_patch_area,
                "unresolved_ester_group_area_A2":float(ua[atomids].sum())})
            row.update({"side":side,"center_x_A":float(origin[0]),"center_y_A":float(origin[1]),"center_z_A":float(origin[2]),
                "ester_group_surface_area_A2":exposed,"carbonyl_C_surface_area_A2":float(atomarea[atomids[0]]),
                "carbonyl_O_surface_area_A2":float(atomarea[atomids[1]]),"ester_O_surface_area_A2":float(atomarea[atomids[2]]),
                "patch_surface_area_A2":pa,"patch_triangles":len(local),"patch_contributing_chains":len(set(segs[ow[local]])),
                "roughness_plane_normal_rms_A":float(np.sqrt(max(0,eig[0]))),
                "height_p95_minus_p05_A":float(np.diff(weighted_quantile(delta[:,2],ww,[.05,.95]))[0]),
                "plane_tilt_deg":float(np.degrees(np.arccos(np.clip(abs(normal[2]),0,1)))),
                "plane_variance_fraction":float(eig[0]/eig.sum()),
                "carbonyl_to_outward_patch_normal_deg":float(np.degrees(np.arccos(np.clip(np.dot(carbonyl,normal),-1,1)))),
                "hydrophilic_proxy_O_area_fraction":float(frac[3]+frac[4]),
                "hydrophobic_proxy_CH_area_fraction":float(frac[0]+frac[1]),
                "carbonyl_C_separate_area_fraction":float(frac[2])})
            row.update({k+"_area_fraction":float(val) for k,val in zip(CATEGORIES,frac)})
            rows.append(row)
            if save:patchids.append(local.astype(np.int32));offsets.append(offsets[-1]+len(local))
        rows_by_side[side]=rows
        chem=np.bincount(cc,weights=w,minlength=5)/w.sum()
        summary={"surface_esters":len(rows),"surface_carbonyl_C":sum(r["carbonyl_C_surface_area_A2"]>0 for r in rows),
            "surface_esters_area_ge_1_A2":sum(r["ester_group_surface_area_A2"]>=1 for r in rows),
            "surface_area_A2":float(w.sum()),"surface_triangles":len(w),
            "hydrophilic_proxy_O_area_fraction":float(chem[3]+chem[4]),"hydrophobic_proxy_CH_area_fraction":float(chem[0]+chem[1]),
            "carbonyl_C_separate_area_fraction":float(chem[2]),"chemical_area_fractions":dict(zip(CATEGORIES,map(float,chem)))}
        for key in ["patch_surface_area_A2","roughness_plane_normal_rms_A","height_p95_minus_p05_A","plane_tilt_deg",
                    "hydrophilic_proxy_O_area_fraction","hydrophobic_proxy_CH_area_fraction","patch_contributing_chains"]:
            vals=np.array([r[key] for r in rows],float)
            summary["per_ester_"+key]={"mean":float(vals.mean()),"median":float(np.median(vals)),"p10":float(np.quantile(vals,.1)),"p90":float(np.quantile(vals,.9))}
        summary["unresolved_surface_area_A2"]=float(uww.sum())
        summary["unresolved_fraction_of_known_plus_unresolved_area"]=float(uww.sum()/(w.sum()+uww.sum()))
        summary["uncertain_only_ester_membership_count"]=sum(atomarea[z["atom_indices"]].sum()==0 and ua[z["atom_indices"]].sum()>0 for z in sites)
        summary["complete_patch_count"]=sum(r["unresolved_patch_area_A2"]==0 for r in rows)
        summary["partial_patch_count"]=len(rows)-summary["complete_patch_count"]
        complete=[r for r in rows if r["unresolved_patch_area_A2"]==0]
        for key in ["roughness_plane_normal_rms_A","height_p95_minus_p05_A","patch_surface_area_A2"]:
            val=[r[key] for r in complete]
            summary["complete_patches_"+key]={"n":len(val),"median":float(np.median(val)) if val else None,
                "p10":float(np.quantile(val,.1)) if val else None,"p90":float(np.quantile(val,.9)) if val else None}
        summaries[side]=summary
        if save:
            table(out/f"esters_{side}_15A.csv",rows)
            np.savez_compressed(out/f"surface_{side}_faces.npz",centroids_A=pts.astype(np.float32),areas_A2=w.astype(np.float32),
                normals=nn.astype(np.float32),owner_atom_index=ow,chemical_category=cc,full_mesh_face_indices=meshids,
                category_names=np.array(CATEGORIES),box_A=box)
            np.savez_compressed(out/f"patches_{side}_15A.npz",site_ids=np.array([r["site_id"] for r in rows]),
                offsets=np.array(offsets,np.int64),surface_face_indices=np.concatenate(patchids) if patchids else np.array([],np.int32))
            # Genuine mesh subset for 3D inspection, using original vertex coordinates.
            used=np.unique(faces[meshids]);newix=np.full(len(v),-1,np.int32);newix[used]=np.arange(len(used))
            with (out/f"surface_{side}.obj").open("w") as f:
                f.write("# PyMOL molecular surface, Angstrom; central periodic XY cell\n")
                for p in v[used]:f.write("v %.6f %.6f %.6f\n"%tuple(p))
                for tri in newix[faces[meshids]]+1:f.write("f %d %d %d\n"%tuple(tri))
    if save:
        allrows=[]
        for site in sites:
            row={k:v for k,v in site.items() if k!="atom_indices"}
            for side in ["top","bottom"]:
                ar=float(atom_areas[side][site["atom_indices"]].sum())
                row[side+"_ester_group_surface_area_A2"]=ar
                uar=float(uncertain_atom_areas[side][site["atom_indices"]].sum())
                row[side+"_unresolved_surface_area_A2"]=uar
                row[side+"_status"]="SURFACE_ASSIGNED" if ar>0 else ("NOT_EVALUATED_AMBIGUOUS_COMPONENTS" if uar>0 else "NOT_ON_RENDERED_SURFACE")
            allrows.append(row)
        table(out/"all_8000_esters.csv",allrows)
        plot_results(out,rows_by_side,box)
    summaries["z_midplane_A"]=mid
    summaries["surface_definition"]="PyMOL SES triangles with centroid in central periodic xy cell; top/bottom split at polymer median z. Closed internal cavities excluded. All chain ends retained."
    return summaries,rows_by_side

def plot_results(out,rows,box):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap,BoundaryNorm
    fig,axs=plt.subplots(2,2,figsize=(12,10),constrained_layout=True)
    for i,side in enumerate(["top","bottom"]):
        dat=np.load(out/f"surface_{side}_faces.npz");p=dat["centroids_A"];c=dat["chemical_category"]
        # Back-to-front layering supplies a top/bottom visible face projection.
        order=np.argsort(p[:,2] if side=="top" else -p[:,2])
        stride=max(1,len(order)//180000);order=order[::stride];pp=p[order]
        im=axs[i,0].scatter(pp[:,0],pp[:,1],c=pp[:,2],s=1.1,linewidths=0,cmap="viridis",rasterized=True)
        fig.colorbar(im,ax=axs[i,0],label="Surface z (A)")
        im=axs[i,1].scatter(pp[:,0],pp[:,1],c=c[order],s=1.1,linewidths=0,cmap=ListedColormap(COLORS),norm=BoundaryNorm(np.arange(-.5,5.5),5),rasterized=True)
        cb=fig.colorbar(im,ax=axs[i,1],ticks=range(5));cb.ax.set_yticklabels(CATEGORIES)
        axs[i,0].set_title(side+": PyMOL molecular surface height")
        axs[i,1].set_title(side+": chemical hydrophilic/hydrophobic proxies")
        for ax in axs[i]:
            ax.set(xlabel="x (A), periodic",ylabel="y (A), periodic",xlim=(0,box[0]),ylim=(0,box[1]),aspect="equal")
    fig.savefig(out/"top_bottom_surface_maps.png",dpi=160);plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
    for side,color in [("top","#c54b42"),("bottom","#377eb8")]:
        for ax,key,label in zip(axs,["roughness_plane_normal_rms_A","hydrophilic_proxy_O_area_fraction","patch_surface_area_A2"],
                ["15 A patch roughness RMS (A)","Hydrophilic proxy: oxygen area fraction","15 A patch surface area (A2)"]):
            ax.hist([r[key] for r in rows[side]],bins=25,histtype="step",color=color,linewidth=1.7,label=f"{side}, n={len(rows[side])}")
            ax.set(xlabel=label,ylabel="Ester sites");ax.legend()
    fig.savefig(out/"patch_distributions.png",dpi=170);plt.close(fig)

def run_analysis(args):
    global Chem
    from rdkit import Chem
    project=args.project.resolve();out=args.output.resolve()
    assert out.is_relative_to(project) and out!=project
    out.mkdir(parents=True,exist_ok=False)
    log=out/"RUN_LOG.tsv"
    def record(status,detail):
        with log.open("a") as f:f.write(now()+"\t"+status+"\t"+detail.replace("\n"," | ")+"\n")
    record("STARTED",json.dumps({"command":sys.argv,"script_sha256":sha(__file__)}))
    try:
        xyz,radii,box,cat,segs,sites,audit=parse_inputs(project,1.4)
        vdw=radii-1.4
        dump(out/"input_audit.json",audit)
        halo=halo_input(xyz,vdw,cat,box)
        np.savez_compressed(out/"mesh_input.npz",**halo)
        exterior_selftest()
        base_summary=[];base_rows=[]
        for q in [1,2]:
            mesh=args.reuse_mesh.resolve()/f"pymol_surface_q{q}.npz"
            assert mesh.is_file(),"Required precomputed mesh missing"
            record("REUSE_MESH",json.dumps({"path":str(mesh),"sha256":sha(mesh)}))
            print("MESH_REUSED",q,flush=True)
            summary,rows=process_mesh(mesh,halo,xyz,cat,segs,sites,box,out,q==2)
            dump(out/f"summary_q{q}.json",summary);base_summary.append(summary);base_rows.append(rows)
            print("SUMMARY",q,json.dumps(summary),flush=True)
        convergence={}
        for side in ["top","bottom"]:
            a,b=[s[side] for s in base_summary]
            aa,bb=[{r["site_id"]:r for r in rr[side]} for rr in base_rows]
            common=set(aa)&set(bb)
            cross=[]
            for k in sorted(common):
                cross.append({"site_id":k,"patch_area_relative_difference":abs(aa[k]["patch_surface_area_A2"]/bb[k]["patch_surface_area_A2"]-1),
                    "O_fraction_absolute_difference":abs(aa[k]["hydrophilic_proxy_O_area_fraction"]-bb[k]["hydrophilic_proxy_O_area_fraction"]),
                    "roughness_absolute_difference_A":abs(aa[k]["roughness_plane_normal_rms_A"]-bb[k]["roughness_plane_normal_rms_A"])})
            table(out/f"mesh_quality_comparison_{side}.csv",cross)
            ar=abs(a["surface_area_A2"]/b["surface_area_A2"]-1)
            ce=max(abs(a["chemical_area_fractions"][k]-b["chemical_area_fractions"][k]) for k in CATEGORIES)
            reliable_a={k for k,r in aa.items() if r["ester_group_surface_area_A2"]>=1}
            reliable_b={k for k,r in bb.items() if r["ester_group_surface_area_A2"]>=1}
            rj=len(reliable_a&reliable_b)/len(reliable_a|reliable_b)
            convergence[side]={"surface_area_relative_difference":ar,"max_chemical_fraction_absolute_difference":ce,
                "surface_site_jaccard":len(common)/len(set(aa)|set(bb)),"area_ge_1_A2_site_jaccard":rj,
                "q1_only_sites":sorted(set(aa)-set(bb)),"q2_only_sites":sorted(set(bb)-set(aa)),
                "numerical_gate":"PASS" if ar<.02 and ce<.01 and rj>=.98 else "NOT_CONVERGED",
                "patch_quality_differences":{k:{"median":float(np.median([r[k] for r in cross])),"p95":float(np.quantile([r[k] for r in cross],.95)),
                    "max":float(max(r[k] for r in cross))} for k in ["patch_area_relative_difference","O_fraction_absolute_difference","roughness_absolute_difference_A"]}}
        dump(out/"convergence.json",convergence)
        membership={str(q):json.loads((out/f"pymol_surface_q{q}_exterior_audit.json").read_text())["containment_gate"] for q in [1,2]}
        audit2={"technical_status":"COMPLETE","surface_membership_status":membership,"numerical_status":{s:convergence[s]["numerical_gate"] for s in convergence},
            "input":audit,"script":str(Path(__file__).resolve()),"script_sha256":sha(__file__),
            "method":"PyMOL 3.2.0a0 solvent-excluded molecular surface, solvent_radius 1.4 A, explicit-atom radii C=1.70 O=1.52 H=1.20 A, surface_mode=1. Quality 1 and 2 compared. Periodic xy halo 8 A; z is open water side.",
            "top_bottom":"Nested cavity shells excluded by three-axis ray containment; non-nested exterior shell(s) retained. Central-cell exterior face centroids split at median polymer z. No chain-tail removal. Not vertical-only envelope.",
            "uncertainty":"Ambiguous shell containment is NOT_EVALUATED, never forced inside/outside. Global fractions and shape metrics describe confirmed exterior faces only. Per-site unresolved area flags incomplete patches; complete-patch summaries are separate. Membership counts of uncertain-only esters are reported separately.",
            "chemistry":"Area-weighted nearest VDW-sphere mapping at each face centroid; H inherits bonded carbon class. Oxygen = hydrophilic proxy, aromatic/aliphatic CH = hydrophobic proxy; carbonyl C separate. Not a calibrated hydrophobic potential or electrostatic potential.",
            "patch":"Same-side face centroids within 15 A of carbonyl carbon under xy minimum-image distance; full face areas summed. NPZ stores face index lists; boundary triangle approximation assessed by mesh-quality comparison.",
            "shape":"Area-weighted PCA plane-normal RMS; weighted z height p95-p05; plane tilt; C=O angle to fitted outward normal. No enzyme accessibility or catalytic activity conclusion.",
            "statistics":"Overlapping patches, one input snapshot. Descriptive statistics only; no independent replicate, dynamics or significance claims.",
            "software":{"python":sys.version,"numpy":np.__version__,"rdkit":Chem.rdBase.rdkitVersion},"command":sys.argv}
        dump(out/"analysis_audit.json",audit2)
        runbook="Reproduce with the original scientific Python environment:\n"+" ".join(sys.argv)+"\n\n"+json.dumps(audit2,indent=2)+"\n\nPreserve failed/accepted outputs. Use a new output directory for any rerun. mesh_input.npz is generated from validated input; no original input is changed. pymol_surface_q2.npz is the complete halo mesh. surface_top/bottom_faces.npz contains the central-cell triangles. patches files index these per-side face tables. All lengths A and areas A2.\n"
        (out/"RUNBOOK.md").write_text(runbook)
        # Viewer uses stored triangles directly; no recomputation or hidden atom deletion.
        viewer='''from pathlib import Path\nimport numpy as np\nfrom pymol import cmd\nfrom pymol.cgo import BEGIN,END,TRIANGLES,COLOR,NORMAL,VERTEX\nroot=Path(__file__).resolve().parent\nmesh=np.load(root.parent/"results"/"pymol_surface_q2.npz");v=mesh["vertices_A"];f=mesh["faces"]\ncolors='''+repr(COLORS.tolist())+'''\nfor side in ["top","bottom"]:\n d=np.load(root/("surface_"+side+"_faces.npz"));cgo=[BEGIN,TRIANGLES]\n for fi,cat,n in zip(d["full_mesh_face_indices"],d["chemical_category"],d["normals"]):\n  cgo.extend([COLOR,*colors[int(cat)],NORMAL,*map(float,n)])\n  for p in v[f[fi]]:cgo.extend([VERTEX,*map(float,p)])\n cgo.append(END);cmd.load_cgo(cgo,"PET_"+side+"_chemical_surface")\ncmd.bg_color("white");cmd.orient("PET_*_chemical_surface")\n'''
        (out/"view_surfaces.py").write_text(viewer)
        assert sha(Path(audit["pdb"]))==audit["pdb_sha256"]
        dump(out/"SHA256SUMS.json",{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!="RUN_LOG.tsv"})
        record("COMPLETE",json.dumps(audit2["numerical_status"]))
        print("FINAL",json.dumps({"output":str(out),"numerical":audit2["numerical_status"],"membership":membership,"convergence":convergence}),flush=True)
    except Exception:
        record("FAILED",traceback.format_exc());raise
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--reuse-mesh",type=Path);ap.add_argument("--selftest-exterior",action="store_true")
    ap.add_argument("--export",nargs=3);ap.add_argument("--selftest-pymol",action="store_true")
    ap.add_argument("--project",type=Path);ap.add_argument("--output",type=Path);ap.add_argument("--pymol-python",type=Path)
    a=ap.parse_args()
    if a.selftest_exterior:exterior_selftest()
    elif a.selftest_pymol:export_mesh(None,None,2,True)
    elif a.export:export_mesh(a.export[0],a.export[1],int(a.export[2]))
    else:run_analysis(a)
if __name__=="__main__":main()
