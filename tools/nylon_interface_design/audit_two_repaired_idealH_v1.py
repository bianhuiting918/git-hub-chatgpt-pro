"""Independent audit using idealized ARG H; not minimized-H or all-gates PASS."""
from pathlib import Path
import hashlib
r=Path(__file__).resolve().parent
p=r/'audit_supported_arg_repair_v1.py'
src=p.read_text()
assert hashlib.sha256(p.read_bytes()).hexdigest()=='a7699591fd7e294fa86c24fc1316bf6b92ae4bc70d29d29040dc47f79beed017'
for cid in ['SC4916877','SC4242958']:
 s=src.replace("o=d/'route6_supported_ARG_SC4559268_v1';path=o/'fixed_core_repair_v1/after.pdb'", "o=d/'route6_top50_local_generation_v1';path=o/'"+cid+"_seed2_multisite_repair_v3_rebuildH/after.pdb'")
 start=s.index("z=np.load(o/");end=s.index("sites=json.loads",start)
 s=s[:start]+"ca=[]\nfor line in (o/'candidate_specific_shapes_v2/"+cid+"/motif.pdb').read_text().splitlines():\n if line.startswith('ATOM') and line[12:16].strip()=='CA':ca.append([float(line[i:i+8]) for i in [30,38,46]])\ncenter=np.mean(ca,axis=0);rot=np.eye(3);shift=np.zeros(3);chemical=q+center\n"+s[end:]
 s=s.replace("v['id']=='SC4559268'","v['id']=='"+cid+"'")
 needle="assert sorted(a)==list(range(1,127))"
 s=s.replace(needle,needle+"\nimport pyrosetta as py\npy.init('-mute all')\npose=py.pose_from_pdb(str(path));py.rosetta.core.conformation.idealize_hydrogens(pose.residue(63),pose.conformation())\nfor h in ['1HH1','2HH1','1HH2','2HH2']:\n v=pose.residue(63).xyz(h);a[63][h]=np.array([v.x,v.y,v.z])\n")
 s=s.replace("dest=o/'fixed_core_repair_v1'/", "dest=o/'"+cid+"_seed2_multisite_repair_v3_rebuildH'/")
 s=s.replace("'independent_audit_'","'idealH_independent_audit_'")
 s=s.replace("'INDEPENDENT_GEOMETRY_AUDIT_NOT_FINAL_SEQUENCE_PASS'","'IDEALIZED_ARG_H_GEOMETRY_ONLY_NOT_FINAL_PASS'")
 s=s.replace("'122 Gly placeholder residues; no designed sequence or independent prediction'","'Gly placeholder scaffold, no designed sequence or independent prediction','ARG hydrogens geometrically rebuilt for audit only; not claimed optimized or stable','Internal contacts retain original generic-radius screen and are not MolProbity clashes','Complete core H, attack geometry and hydrogen-clash checks remain NOT_EVALUATED'")
 compile(s,str(p),'exec')
 exec(compile(s,str(p),'exec'),dict(__file__=str(Path(__file__).resolve()),__name__='__main__'))
