#!/usr/bin/env python3
"""Chemical eligibility of the 20 internal amino acids for a two-donor design."""
from pathlib import Path
import json,datetime,hashlib
from rdkit import Chem,rdBase
ROOT=Path(__file__).resolve().parent
SIDE={'ALA':'C','ARG':'CCCNC(=[NH2+])N','ASN':'CC(=O)N','ASP':'CC(=O)[O-]','CYS':'CS','GLN':'CCC(=O)N','GLU':'CCC(=O)[O-]','HIS':'Cc1cnc[nH]1','ILE':'[C@@H](C)CC','LEU':'CC(C)C','LYS':'CCCC[NH3+]','MET':'CCSC','PHE':'Cc1ccccc1','SER':'CO','THR':'[C@H](O)C','TRP':'Cc1c[nH]c2ccccc12','TYR':'Cc1ccc(O)cc1','VAL':'C(C)C'}
def entry(aa,smiles,state):
 m=Chem.MolFromSmiles(smiles);assert m is not None,(aa,state)
 Chem.SanitizeMol(m);h=Chem.AddHs(m)
 donors=[]
 for a in h.GetAtoms():
  if a.GetAtomMapNum()==5:continue # N-methyl cap NH is NOT the target residue
  hs=[b.GetIdx() for b in a.GetNeighbors() if b.GetAtomicNum()==1]
  if a.GetSymbol() in ('N','O','S') and hs:
   donors.append(dict(atom_index=a.GetIdx(),element=a.GetSymbol(),hydrogen_indices=hs,role='backbone_NH' if a.GetAtomMapNum()==1 else 'sidechain',formal_charge=a.GetFormalCharge()))
 eligible=len(donors)>=2
 return dict(amino_acid=aa,microstate=state,smiles=Chem.MolToSmiles(m,isomericSmiles=True),donor_groups=donors,n_distinct_donor_atoms=len(donors),status='ELIGIBLE_FOR_GEOMETRY_SCAN' if eligible else 'NOT_APPLICABLE_FEWER_THAN_TWO_DONOR_ATOMS',note='Cys S-H must not be treated as equivalent to N-H/O-H' if aa=='CYS' else '',hydrogens_on_same_atom='one donor atom, not independent donors')
if __name__=='__main__':
 rows=[]
 for aa in sorted(set(SIDE)|{'GLY','PRO'}):
  if aa=='GLY':sm='CC(=O)[NH:1][CH2:2][C:3](=[O:4])[NH:5]C'
  elif aa=='PRO':sm='CC(=O)[N:1]1CCC[C@H:2]1[C:3](=[O:4])[NH:5]C'
  else:sm='CC(=O)[NH:1][C@@H:2]('+SIDE[aa]+')[C:3](=[O:4])[NH:5]C'
  rows.append(entry(aa,sm,'neutral except Lys+/Arg+/Asp-/Glu-; His neutral tautomer 1'))
 assert len(rows)==20
 assert next(r for r in rows if r['amino_acid']=='PRO')['n_distinct_donor_atoms']==0
 assert next(r for r in rows if r['amino_acid']=='GLY')['n_distinct_donor_atoms']==1
 assert sum(r['status']=='ELIGIBLE_FOR_GEOMETRY_SCAN' for r in rows)==10
 for state,side in [('neutral tautomer 2','Cc1c[nH]cn1'),('protonated +1','Cc1c[nH+]c[nH]1')]:
  rows.append(entry('HIS','CC(=O)[NH:1][C@@H:2]('+side+')[C:3](=[O:4])[NH:5]C',state))
 report=dict(rdkit=rdBase.rdkitVersion,standard_amino_acids=20,microstate_records=len(rows),scope='Ac-AA-NMe geometry template; cap donor excluded; not a protein terminal residue; chemical eligibility only, no material scan yet',results=rows)
 (ROOT/'chemical_eligibility.json').open('x').write(json.dumps(report,indent=2))
 (ROOT/'RUNBOOK.md').open('x').write('Single-residue dual-donor comparison on PET 578, PA6 126, PA66 136 existing sites.\nFirst stage: CPU Python chemical_eligibility.py; records capped internal-residue donor chemistry. No substrate pose sampling performed in this stage.\nDo not count N-methyl cap NH, and do not equate two H on one N with two distinct donor atoms. Protonation and Cys SH are explicit categories. Keep previous outputs unchanged.\n')
 with (ROOT/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__).resolve()),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),command='CPU Python chemical_eligibility.py',inputs='explicit 20 Ac-AA-NMe molecular graphs and His microstates',outputs='chemical_eligibility.json',exit_code=0))+'\n')
 print(json.dumps(dict(standard_residues=20,eligible=[r['amino_acid'] for r in rows[:20] if r['status']=='ELIGIBLE_FOR_GEOMETRY_SCAN'],not_applicable=[r['amino_acid'] for r in rows[:20] if r['status']!='ELIGIBLE_FOR_GEOMETRY_SCAN'],status='CHEMISTRY_AUDIT_ONLY')))
