from pathlib import Path
import csv, json, math, datetime
ROOT=Path(__file__).parent; IN=ROOT/'inputs'
AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
VDW={'H':1.20,'C':1.70,'N':1.55,'O':1.52,'S':1.80,'P':1.80}; COV={'H':0.31,'C':0.76,'N':0.71,'O':0.66,'S':1.05,'P':1.07}
def elem(l):
    e=l[76:78].strip(); return e.upper() if e else ''.join(c for c in l[12:16].strip() if c.isalpha())[0].upper()
def xyz(l): return [float(l[30:38]),float(l[38:46]),float(l[46:54])]
def vadd(a,b): return [a[i]+b[i] for i in range(3)]
def vsub(a,b): return [a[i]-b[i] for i in range(3)]
def vdot(a,b): return sum(a[i]*b[i] for i in range(3))
def norm(a): return math.sqrt(vdot(a,a))
def mean(P): return [sum(p[i] for p in P)/len(P) for i in range(3)]
def matvec(M,v): return [sum(M[i][j]*v[j] for j in range(len(v))) for i in range(len(M))]
def rotmat_from_q(q):
    w,x,y,z=q
    return [[1-2*y*y-2*z*z,2*x*y-2*z*w,2*x*z+2*y*w],[2*x*y+2*z*w,1-2*x*x-2*z*z,2*y*z-2*x*w],[2*x*z-2*y*w,2*y*z+2*x*w,1-2*x*x-2*y*y]]
def det3(R): return R[0][0]*(R[1][1]*R[2][2]-R[1][2]*R[2][1])-R[0][1]*(R[1][0]*R[2][2]-R[1][2]*R[2][0])+R[0][2]*(R[1][0]*R[2][1]-R[1][1]*R[2][0])
def apply_rt(v,R,t): return vadd(matvec(R,v),t)
def rmsd_after_transform(P,Q,R,t): return math.sqrt(sum(vdot(vsub(apply_rt(p,R,t),q),vsub(apply_rt(p,R,t),q)) for p,q in zip(P,Q))/len(P))
def jacobi_symmetric(A):
    n=len(A); a=[r[:] for r in A]; V=[[1.0 if i==j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(100):
        p,q=max(((i,j) for i in range(n) for j in range(i+1,n)), key=lambda ij: abs(a[ij[0]][ij[1]]))
        if abs(a[p][q])<1e-12: break
        tau=(a[q][q]-a[p][p])/(2*a[p][q]); tt=(1 if tau>=0 else -1)/(abs(tau)+math.sqrt(1+tau*tau)); c=1/math.sqrt(1+tt*tt); s=tt*c
        app,aqq,apq=a[p][p],a[q][q],a[p][q]; a[p][p]=app-tt*apq; a[q][q]=aqq+tt*apq; a[p][q]=a[q][p]=0.0
        for k in range(n):
            if k not in (p,q):
                akp,akq=a[k][p],a[k][q]; a[k][p]=a[p][k]=c*akp-s*akq; a[k][q]=a[q][k]=s*akp+c*akq
        for k in range(n):
            vkp,vkq=V[k][p],V[k][q]; V[k][p]=c*vkp-s*vkq; V[k][q]=s*vkp+c*vkq
    return [a[i][i] for i in range(n)],V
def kabsch(P,Q):
    Pc,Qc=mean(P),mean(Q); A=[vsub(p,Pc) for p in P]; B=[vsub(q,Qc) for q in Q]
    S=[[sum(a[i]*b[j] for a,b in zip(A,B)) for j in range(3)] for i in range(3)]
    Sxx,Sxy,Sxz=S[0]; Syx,Syy,Syz=S[1]; Szx,Szy,Szz=S[2]
    K=[[Sxx+Syy+Szz,Syz-Szy,Szx-Sxz,Sxy-Syx],[Syz-Szy,Sxx-Syy-Szz,Sxy+Syx,Szx+Sxz],[Szx-Sxz,Sxy+Syx,-Sxx+Syy-Szz,Syz+Szy],[Sxy-Syx,Szx+Sxz,Syz+Szy,-Sxx-Syy+Szz]]
    vals,V=jacobi_symmetric(K); imax=max(range(4), key=lambda i: vals[i]); q=[V[i][imax] for i in range(4)]; n=math.sqrt(sum(x*x for x in q)); q=[x/n for x in q]
    R=rotmat_from_q(q); t=vsub(Qc,matvec(R,Pc)); return R,t,rmsd_after_transform(P,Q,R,t)
def parse_atoms(path):
    return [{'line':l,'rec':l[:6].strip(),'serial':int(l[6:11]),'name':l[12:16].strip(),'resn':l[17:20].strip(),'chain':l[21].strip() or ' ','resi':int(l[22:26]),'coord':xyz(l),'element':elem(l)} for l in path.read_text().splitlines() if l.startswith(('ATOM','HETATM'))]
def residues(atoms):
    d={}
    for a in atoms:
        if a['rec']=='ATOM': d.setdefault((a['chain'],a['resi'],a['resn']),[]).append(a)
    return sorted(d.items(), key=lambda kv:(kv[0][0],kv[0][1]))
def nw(a,b):
    n,m=len(a),len(b); S=[[0]*(m+1) for _ in range(n+1)]; T={}
    for i in range(1,n+1): S[i][0]=-2*i; T[i,0]=1
    for j in range(1,m+1): S[0][j]=-2*j; T[0,j]=2
    for i in range(1,n+1):
        for j in range(1,m+1):
            opts=[(S[i-1][j-1]+(2 if a[i-1]==b[j-1] else -1),0),(S[i-1][j]-2,1),(S[i][j-1]-2,2)]
            S[i][j],T[i,j]=max(opts, key=lambda x:(x[0],-x[1]))
    i,j=n,m; pairs=[]
    while i or j:
        k=T[i,j]
        if k==0: pairs.append((i-1,j-1)); i-=1; j-=1
        elif k==1: pairs.append((i-1,None)); i-=1
        else: pairs.append((None,j-1)); j-=1
    return list(reversed(pairs))
def infer_bonds(atoms):
    bonds=set()
    for i,a in enumerate(atoms):
        for j,b in enumerate(atoms[i+1:],i+1):
            d=norm(vsub(a['coord'],b['coord']))
            if 0.4<d<=COV.get(a['element'],0.77)+COV.get(b['element'],0.77)+0.45: bonds.add((i,j))
    return bonds
def rings6_carbon(atoms,bonds):
    carb=[i for i,a in enumerate(atoms) if a['element']=='C']; adj={i:set() for i in carb}
    for i,j in bonds:
        if i in adj and j in adj: adj[i].add(j); adj[j].add(i)
    rings=set()
    def dfs(start,path):
        for nb in adj[path[-1]]:
            if nb==start and len(path)==6: rings.add(tuple(sorted(path)))
            elif nb not in path and len(path)<6 and nb>=start: dfs(start,path+[nb])
    for c in carb: dfs(c,[c])
    return rings
def carbon_rings_from_pdb(path):
    atoms=[a for a in parse_atoms(Path(path)) if a['rec']=='HETATM' and a['resn']=='LG1']; return rings6_carbon(atoms,infer_bonds(atoms))
def fmt_atom(rec,serial,name,resn,chain,resi,coord,occ,b,element): return f"{rec:<6}{serial:5d} {name:<4} {resn:>3} {chain:1}{resi:4d}    {coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}{occ:6.2f}{b:6.2f}          {element:>2}"
def compute_model():
    iccg=parse_atoms(IN/'ICCG_active_chainA_heavy.pdb'); ref=parse_atoms(IN/'IsPETase_WT_LG1_reference.pdb'); lig=[dict(a,coord=a['coord'][:]) for a in ref if a['rec']=='HETATM' and a['resn']=='LG1']
    rr,ir=residues(ref),residues(iccg); rseq=''.join(AA3.get(k[2],'X') for k,_ in rr); iseq=''.join(AA3.get(k[2],'X') for k,_ in ir)
    pairs=nw(rseq,iseq); amap={rr[i][0][1]:ir[j][0][1] for i,j in pairs if i is not None and j is not None}; gaps={str(rr[i][0][1]):None for i,j in pairs if i is not None and j is None}
    lig_heavy=[a for a in lig if a['element']!='H']; rdict={k[1]:v for k,v in rr}; idict={k[1]:v for k,v in ir}
    ref_pocket={key[1] for key,ats in rr if any(a['element']!='H' and any(norm(vsub(a['coord'],l['coord']))<=5.0 for l in lig_heavy) for a in ats)}; anchors=ref_pocket|{159,160,206,237,238}
    P=[]; Q=[]; fitmap={}
    for r in sorted(anchors):
        if r in amap:
            added=0; ra={a['name']:a for a in rdict.get(r,[])}; ia={a['name']:a for a in idict.get(amap[r],[])}
            for nm in ('N','CA','C'):
                if nm in ra and nm in ia: P.append(ra[nm]['coord']); Q.append(ia[nm]['coord']); added+=1
            if added: fitmap[str(r)]=amap[r]
    R,t,rmsd=kabsch(P,Q); original_lig=[dict(a,coord=a['coord'][:]) for a in lig]
    for a in lig: a['coord']=apply_rt(a['coord'],R,t)
    lig_heavy_t=[a for a in lig if a['element']!='H']; newp={key[1] for key,ats in ir if any(a['element']!='H' and any(norm(vsub(a['coord'],l['coord']))<=5.0 for l in lig_heavy_t) for a in ats)}
    high=set(); cys=set(); catalytic={165,210,242}
    for row in csv.DictReader(open(IN/'ICCG_legacy_mask_components.tsv'),delimiter='\t'):
        r=int(row['resi'])
        if row['high_conservation']=='1': high.add(r)
        if row['native_cysteine']=='1': cys.add(r)
        if row['catalytic']=='1': catalytic.add(r)
    fixed=sorted(catalytic|newp|high|cys); bonds=infer_bonds(lig); assert len(rings6_carbon(lig,bonds))==2
    mind=999.0; maxov=-999.0; worst=''
    for pa in iccg:
        if pa['element']=='H': continue
        for la in lig_heavy_t:
            d=norm(vsub(pa['coord'],la['coord'])); ov=VDW.get(pa['element'],1.7)+VDW.get(la['element'],1.7)-d
            if d<mind: mind=d
            if ov>maxov: maxov=ov; worst=f"A:{pa['resi']}:{pa['name']}({pa['element']})--L:301:{la['name']}({la['element']}) distance={d:.3f} overlap={ov:.3f}"
    status='PASS' if maxov<=0.8 else 'PROVISIONAL_DISPLAY_ONLY_LG1_CLASH'; keymap={str(k):(amap[k] if k in amap else None) for k in (159,160,206,237,238)}
    audit={'geometry_status':status,'min_nonbonded_protein_lg1_heavy_distance':round(mind,3),'max_vdw_overlap':round(maxov,3),'worst_contact':worst,'alignment_anchor_residue_count':len(fitmap),'alignment_backbone_atom_count':len(P),'backbone_fit_rmsd':round(rmsd,6),'rotation_det':round(det3(R),9),'rotation_matrix':R,'translation':t,'alignment_fit_mappings':fitmap,'alignment_key_mappings':keymap,'alignment_reference_gaps':gaps,'catalytic_residues':sorted(catalytic),'new_LG1_5A_pocket_residues':sorted(newp),'high_conservation_residues':sorted(high),'native_cysteine_residues':sorted(cys),'fixed_residues':fixed,'counts':{'protein_residues':len(ir),'lg1_atoms':len(lig),'lg1_heavy_atoms':len(lig_heavy_t),'fixed_residues':len(fixed),'tsv_rows':len(ir)},'note':'Viewer mask is for visualization only and does not replace any production sequence-generation mask.'}
    return {'audit':audit,'rotation':R,'translation':t,'fit_reference_points':P,'fit_iccg_points':Q,'iccg':iccg,'ir':ir,'ligand_original':original_lig,'ligand_transformed':lig,'bonds':bonds,'fixed':set(fixed),'catalytic':catalytic,'newp':newp,'high':high,'cys':cys}
def write_outputs(model, append_history=False, red_line=None, green_status=None):
    audit=model['audit']; iccg=model['iccg']; lig=model['ligand_transformed']; bonds=model['bonds']; fixed=model['fixed']; out=(IN/'6THT_secondary_structure_records.pdb').read_text().splitlines()
    if audit['geometry_status']!='PASS': out.append('REMARK PROVISIONAL_DISPLAY_ONLY_LG1_CLASH: display-only ligand pose; not computation-ready')
    serial=1
    for a in iccg:
        out.append(fmt_atom('ATOM',serial,a['name'],a['resn'],'A',a['resi'],a['coord'],1.0,0.0 if a['resi'] in fixed else 99.0,a['element'])); serial+=1
    for a in lig:
        a['newserial']=serial; out.append(fmt_atom('HETATM',serial,a['name'],'LG1','L',301,a['coord'],1.0,50.0,a['element'])); serial+=1
    for i,j in sorted(bonds): out.append(f"CONECT{lig[i]['newserial']:5d}{lig[j]['newserial']:5d}"); out.append(f"CONECT{lig[j]['newserial']:5d}{lig[i]['newserial']:5d}")
    out.append('END'); (ROOT/'ICCG_active_cartoon_fixed_sites_with_LG1.pdb').write_text('\n'.join(out)+'\n')
    with open(ROOT/'ICCG_active_fixed_sites.tsv','w') as f:
        f.write('resi\taa\tfixed\treasons\n')
        for k,ats in model['ir']:
            r=k[1]; reasons=[]
            if r in model['catalytic']: reasons.append('catalytic')
            if r in model['newp']: reasons.append('new_LG1_5A_pocket')
            if r in model['high']: reasons.append('high_conservation')
            if r in model['cys']: reasons.append('native_cysteine')
            f.write(f"{r}\t{k[2]}\t{1 if r in fixed else 0}\t{','.join(reasons)}\n")
    (ROOT/'ICCG_active_fixed_sites_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    (ROOT/'RUNBOOK.md').write_text('# RUNBOOK\n\nRun `python build_iccg_lg1_viewer.py` to regenerate all viewer artifacts, then run `python test_iccg_lg1_viewer.py`. The geometry gate controls whether the viewer is computation-ready or display-only; the visualization mask never replaces a production sequence-generation mask.\n')
def main():
    model=compute_model(); write_outputs(model); print(json.dumps(model['audit'],indent=2)); return model['audit']
if __name__=='__main__': main()
