import unittest,importlib.util
import numpy as np
class Tests(unittest.TestCase):
 def test_nylon_rules(self):
  for mat in ["pa6","pa66"]:
   name=mat+"_tangent_figures_all_sites_v2"
   self.assertIsNotNone(importlib.util.find_spec(name),"nylon plotting implementation missing")
   m=__import__(name)
   from rdkit import Chem
   self.assertEqual([m.bin_index(x) for x in [0,.0025,.01,.02,.020001]],[None,0,1,2,3])
   for z in [[0,0,1],[0,0,-1],[1,0,0]]:
    f=m.frame_from_opening(np.array(z,float),np.array([1.,0,0]))
    self.assertTrue(np.allclose(f.T@f,np.eye(3)))
    self.assertAlmostEqual(np.linalg.det(f),1)
    self.assertTrue(np.allclose(f[:,2],z))
   for smi,expected in [("CC(=O)NCC","ACYL_CAP"),("CCC(=O)NC","ALKYL_CAP"),("CCC(=O)NCC","BACKBONE")]:
    mol=Chem.MolFromSmiles(smi)
    c=next(a for a in mol.GetAtoms() if a.GetSymbol()=="C" and any(b.GetBondTypeAsDouble()==2 and b.GetOtherAtom(a).GetSymbol()=="O" for b in a.GetBonds()))
    li=next(a.GetIdx() for a in c.GetNeighbors() if a.GetSymbol()=="N")
    self.assertEqual(m.cap_kind(mol,c.GetIdx(),li),expected)
if __name__=="__main__":unittest.main()
