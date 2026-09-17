"""Second bounded diagnostic: add trans omega restraints, same starting backbone."""
from pathlib import Path
import hashlib
r=Path(__file__).resolve().parent
p=r/'repair_SC4916877_multisite_v1.py'
src=p.read_text()
assert hashlib.sha256(p.read_bytes()).hexdigest()=='621d300676cf38ed0c7ff2fa35e8b525431c717b6e861015ac925678587cc0d1'
marker="exec(compile(s,str(p),'exec'),dict(__file__=str(Path(__file__).resolve()),__name__='__main__'))"
assert src.count(marker)==1
ns=dict(__file__=str(Path(__file__).resolve()),__name__='diagnostic_adapter')
exec(compile(src.replace(marker,''),str(p),'exec'),ns)
s=ns['s'].replace('SC4916877_seed2_multisite_repair_v1','SC4916877_seed2_multisite_repair_v2_trans')
needle="def vec(i,n):"
assert s.count(needle)==1
s=s.replace(needle,'''# Only constrain peptide torsions containing movable atoms.
omega_constraints=0
for ii in range(1,pose.size()):
 ids=[R.core.id.AtomID(pose.residue(rr).atom_index(nn),rr) for rr,nn in [(ii,'CA'),(ii,'C'),(ii+1,'N'),(ii+1,'CA')]]
 if all((aid.rsd(),aid.atomno()) in fixedkeys for aid in ids):continue
 pose.add_constraint(R.core.scoring.constraints.DihedralConstraint(*ids,R.core.scoring.func.CircularHarmonicFunc(float(np.pi),float(np.radians(10)))))
 omega_constraints+=1
sf.set_weight(R.core.scoring.dihedral_constraint,1.0)
def vec(i,n):''')
s=s.replace("report=dict(status='RUNNING_DIAGNOSTIC',","report=dict(omega_constraint_count=omega_constraints,omega_target_deg=180,omega_sd_deg=10,status='RUNNING_DIAGNOSTIC',")
s=s.replace("Single CPU diagnostic, same ref2015_cart", "Second CPU diagnostic adds trans omega restraints only (180 degrees, SD 10 degrees, weight 1) on movable peptide torsions, same starting RFD backbone; same ref2015_cart")
compile(s,str(p),'exec')
exec(compile(s,str(p),'exec'),dict(__file__=str(Path(__file__).resolve()),__name__='__main__'))
