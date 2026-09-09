import unittest,importlib.util
import numpy as np
class Tests(unittest.TestCase):
 def setUp(self):
  spec=importlib.util.find_spec("pet_tangent_figures_v1")
  self.assertIsNotNone(spec,"plot analysis implementation must exist")
  import pet_tangent_figures_v1 as m
  self.m=m
 def test_bins(self):
  self.assertEqual([self.m.bin_index(v) for v in [0,.0025,.0025001,.01,.010001,.02,.020001]],[None,0,1,1,2,2,3])
 def test_frame(self):
  for z in [[0,0,1],[0,0,-1],[1,0,0]]:
   f=self.m.frame_from_opening(np.array(z,float),np.array([1.,0,0]))
   self.assertTrue(np.allclose(f.T@f,np.eye(3)))
   self.assertAlmostEqual(np.linalg.det(f),1)
   self.assertTrue(np.allclose(f[:,2],z))
 def test_caps(self):
  from rdkit import Chem
  for smi,kind in [("CC(=O)OCC","ACYL_CAP"),("COC(=O)c1ccccc1","ALKYL_CAP"),("CCOC(=O)c1ccccc1","BACKBONE")]:
   mol=Chem.MolFromSmiles(smi)
   c=next(a for a in mol.GetAtoms() if a.GetSymbol()=="C" and any(b.GetBondTypeAsDouble()==2 and b.GetOtherAtom(a).GetSymbol()=="O" for b in a.GetBonds()))
   leave=next(b.GetOtherAtomIdx(c.GetIdx()) for b in c.GetBonds() if b.GetBondTypeAsDouble()==1 and b.GetOtherAtom(c).GetSymbol()=="O")
   self.assertEqual(self.m.cap_kind(mol,c.GetIdx(),leave),kind)
if __name__=="__main__":unittest.main()
