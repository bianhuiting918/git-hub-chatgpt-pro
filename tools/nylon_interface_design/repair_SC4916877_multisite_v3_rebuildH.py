"""Repair initialization bug: rebuild hydrogen coordinates after motif restoration."""
from pathlib import Path
import hashlib
r=Path(__file__).resolve().parent
p=r/'repair_SC4916877_multisite_v2_trans.py'
src=p.read_text()
assert hashlib.sha256(p.read_bytes()).hexdigest()=='3ebfb01477acdb3139385a7ff5a882d23b57bce44b61ed11cc426e4072dab3d8'
prefix,tail=src.rsplit("\nexec(compile(s,str(p),'exec'),",1)
ns=dict(__file__=str(Path(__file__).resolve()),__name__='adapter_only')
exec(compile(prefix,str(p),'exec'),ns)
s=ns['s'].replace('SC4916877_seed2_multisite_repair_v2_trans','SC4916877_seed2_multisite_repair_v3_rebuildH')
needle='fixedkeys={(a.rsd(),a.atomno()) for a in fixed}'
assert s.count(needle)==1
s=s.replace(needle,'''# Rebuild all hydrogens from current heavy atoms, not stale pre-restoration atoms.
for ii in range(1,pose.size()+1):
 res=pose.residue(ii)
 heavy_before=[R.numeric.xyzVector_double_t(res.xyz(an)) for an in range(1,res.nheavyatoms()+1)]
 R.core.conformation.idealize_hydrogens(res,pose.conformation())
 assert max(res.xyz(an).distance(heavy_before[an-1]) for an in range(1,res.nheavyatoms()+1))<1e-10
fixedkeys={(a.rsd(),a.atomno()) for a in fixed}''')
s=s.replace("Second CPU diagnostic adds trans omega", "Third diagnostic fixes stale hydrogen initialization after restoring heavy atoms; all heavy coordinates verified unchanged by hydrogen rebuilding. Same trans omega")
compile(s,str(p),'exec')
exec(compile(s,str(p),'exec'),dict(__file__=str(Path(__file__).resolve()),__name__='__main__'))
