"""Bin-specific CPU shape preparation; pinned existing protocol, no GPU submission."""
from pathlib import Path
import hashlib,sys
r=Path(__file__).resolve().parent
source=r/'prepare_top50_fragment_shapes_v2.py'
assert hashlib.sha256(source.read_bytes()).hexdigest()=='45b4fc8cdd953ec136efbb0af9e6c652b24a7341d5922f0be9c3b8d52e940402'
s=source.read_text()
start=s.index(" for mat in ['PA6','PA66']:",s.index("for row in good:"))
end=s.index("  del fields,local",start)+len("  del fields,local")
block=s[start:end]
first,body=block.split('\n',1)
body=body.replace("sites[int(i)]['material']==mat]","sites[int(i)]['material']==mat and sites[int(i)]['bin']==exposure_bin]")
body=body.replace("s['material']==mat and coverage[pi,i]","s['material']==mat and s['bin']==exposure_bin and coverage[pi,i]")
body=body.replace("name=f'{mat}_{cid}_coverage","name=f'{mat}_{cid}_bin{exposure_bin}_coverage")
body=body.replace("rec=dict(name=name,","rec=dict(exposure_bin=exposure_bin,selected_site_count=k,multi_site=(k>=2),all_exposure_group_denominator=sum(s['material']==mat and s['bin']==exposure_bin for s in sites),name=name,")
body=body.replace("material=mat,original_denominator=","material=mat,exposure_bin=exposure_bin,original_denominator=")
assert "s['bin']==exposure_bin and coverage" in body
assert "sites[int(i)]['bin']==exposure_bin]" in body
s=s[:start]+first+"\n  for exposure_bin in [0,1,2,3]:\n"+'\n'.join(' '+line for line in body.split('\n'))+s[end:]
s=s.replace("candidate_specific_shapes_v2","candidate_specific_binned_shapes_v1")
s=s.replace("route6_top50_local_shapes_v1","route6_top50_local_binned_shapes_v1")
s=s.replace("Pooled first batch; bin labels retained, bin-specific shapes not yet prepared","Independent original exposure bins 0,1,2,3; no pooled selections")
s=s.replace("Original core eligible sites, NOT reduced after adding local backbone","Original core eligible sites WITHIN THE SAME EXPOSURE BIN; NOT reduced after adding local backbone")
s=s.replace("Pooled exposure cohort first; bin-specific jobs remain pending.","Exposure bins are fitted independently. Fractions are conditional on original core-eligible sites in each bin. selected_site_count=1 is explicitly single-site, not multisite consensus. Original pooled outputs are preserved.")
compile(s,str(source),'exec')
if '--check' in sys.argv:
 assert s.count('for exposure_bin in [0,1,2,3]')==1
 assert "rec=dict(exposure_bin=exposure_bin" in s
 assert "bin{exposure_bin}_coverage" in s
 print('BINNED_DERIVED_SCRIPT_COMPILES; original parameters retained; CPU only; source hash verified')
 raise SystemExit(0)
exec(compile(s,str(source),'exec'),{'__file__':str(source),'__name__':'__main__'})
