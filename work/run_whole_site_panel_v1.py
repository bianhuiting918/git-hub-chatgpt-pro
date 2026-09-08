import sys,pathlib,subprocess,json,datetime
R=pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(R))
import chemical_contact_surface_v3 as m
cases=[(mat,cid) for mat,ids in m.CASES.items() for cid in ids if not(mat=='PET' and cid=='M329')]
out=R/'whole_site_panel_v1';out.mkdir(exist_ok=False)
records=[]
for mat,cid in cases:
 name='whole_site_'+mat+'_'+cid+'_v1'
 with (out/(mat+'_'+cid+'.log')).open('x') as log:
  cp=subprocess.run([sys.executable,'-u',str(R/'whole_site_contact_v1.py'),'--material',mat,'--core',cid,'--output',name],cwd=R,stdout=log,stderr=subprocess.STDOUT)
 rec=dict(material=mat,core=cid,exit_status=cp.returncode,output=str(R/name));records.append(rec)
 print(json.dumps(rec),flush=True)
 if cp.returncode:break
with (out/'runs.json').open('x') as f:json.dump(records,f,indent=2)
if any(x['exit_status'] for x in records):sys.exit(1)
