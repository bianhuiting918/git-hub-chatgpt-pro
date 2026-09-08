import pathlib,json,ast,datetime
R=pathlib.Path(__file__).resolve().parent
base=(R/'render_whole_site_M329_water_v1.py').read_text()
records=json.load(open(R/'whole_site_panel_v1/publication/all_core_results.json'))
outputs=[]
for rec in records:
 mat=rec['material'];cid=rec['core']
 if mat not in ['PA6','PA66']:continue
 code=base.replace("P=R/'whole_site_M329_pilot_v1';O=P/'water_preview'","P=R/'whole_site_"+mat+"_"+cid+"_v1';O=P/'nylon_water_preview_v1'")
 code=code.replace("figsize=(16,12)","figsize=(20,12)")
 code=code.replace("titles=['Actual water-facing PET: lipo','Selected PET consensus: lipo','Selected PET consensus: HBA','Candidate enzyme contact face']","titles=['Actual water-facing "+mat+": lipo','Selected "+mat+" consensus: lipo','Selected "+mat+" consensus: HBA','Selected "+mat+" consensus: HBD','Candidate enzyme contact face']")
 code=code.replace('range(4)','range(5)').replace('col==3','col==4')
 code=code.replace("if col==2:\n   values=z['material_HBA_weight'];colors=.93+(np.array([.15,.55,.85])-.93)*values[:,None];colors[~z['chemical_supported']]=.55","if col in [2,3]:\n   values=z['material_HBA_weight' if col==2 else 'material_HBD_weight'];target=np.array([.15,.55,.85]) if col==2 else np.array([.7,.25,.7]);colors=.93+(target-.93)*values[:,None];colors[~z['chemical_supported']]=.55")
 code=code.replace('fig.add_subplot(3,4,row*4+col+1','fig.add_subplot(3,5,row*5+col+1')
 code=code.replace("+'/29)\\n'","+'/'+str(summary['denominator'])+')\\n'")
 code=code.replace('PET M329 |',mat+' '+cid+' |').replace('fixed PET site subset','fixed '+mat+' site subset')
 code=code.replace('HBA: white 0 to blue 1.','HBA: white 0 to blue 1; HBD: white 0 to purple 1.')
 code=code.replace('M329_paired_water_material_enzyme',mat+'_'+cid+'_paired_water_material_enzyme')
 ast.parse(code);exec(compile(code,str(R/'render_nylon_whole_site_v1.py'),'exec'),{'__file__':str(R/'render_nylon_whole_site_v1.py'),'__name__':'__main__'})
 outputs.append(dict(material=mat,core=cid,denominator=rec['denominator'],directory=str(R/('whole_site_'+mat+'_'+cid+'_v1')/'nylon_water_preview_v1')))
 print('RENDERED',mat,cid,flush=True)
with (R/'whole_site_panel_v1/nylon_figures_v1.json').open('x') as f:json.dump(outputs,f,indent=2)
with (R/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now().astimezone().isoformat(),event='render_nylon_whole_site_v1',figures=len(outputs),exit_status=0))+'\n')
