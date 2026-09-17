"""Same repair protocol on a second donor/backbone, not relaxed acceptance gates."""
from pathlib import Path
import hashlib
r=Path(__file__).resolve().parent
p=r/'repair_SC4916877_multisite_v3_rebuildH.py'
src=p.read_text()
assert hashlib.sha256(p.read_bytes()).hexdigest()=='0fea027ea780586e0cd9cd8858b8ba36b628148ade66a2f3ceff26548aea58a9'
prefix,tail=src.rsplit("\nexec(compile(s,str(p),'exec'),",1)
ns=dict(__file__=str(Path(__file__).resolve()),__name__='adapter_only')
exec(compile(prefix,str(p),'exec'),ns)
s=ns['s'].replace('SC4916877','SC4242958')
assert s.count('[204,358,460,543,549,648]')==1
s=s.replace('[204,358,460,543,549,648]','[204,358,460,549,648]')
s=s.replace('Six original PA6 site constraints','Five original PA6 site constraints')
compile(s,str(p),'exec')
exec(compile(s,str(p),'exec'),dict(__file__=str(Path(__file__).resolve()),__name__='__main__'))
