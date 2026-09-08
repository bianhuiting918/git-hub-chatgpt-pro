import unittest,importlib.util,numpy as np
from pathlib import Path
class T(unittest.TestCase):
 def setUp(self):
  p=Path(__file__).with_name('chemical_contact_surface_v1.py')
  self.assertTrue(p.exists(),'chemical surface implementation not yet present')
  s=importlib.util.spec_from_file_location('c',p);self.c=importlib.util.module_from_spec(s);s.loader.exec_module(self.c)
 def test_atom_fragment_conservation(self):
  from rdkit import Chem
  from rdkit.Chem import Crippen
  mol=Chem.AddHs(Chem.MolFromSmiles('CC(=O)NC'))
  a=self.c.fragment_props(mol)
  self.assertAlmostEqual(float(a['logp'].sum()),Crippen.MolLogP(mol),places=6)
  n=next(i for i,x in enumerate(mol.GetAtoms()) if x.GetSymbol()=='N')
  self.assertEqual(a['hba'][n],0);self.assertEqual(a['hbd'][n],1)
 def test_mesh_open_clip(self):
  a=np.arange(-18,19);x,y,z=np.meshgrid(a,a,a,indexing='ij')
  v,f,n=self.c.mesh_patch(-z.astype(float),15,False)
  self.assertLessEqual(np.linalg.norm(v,axis=1).max(),15.00001)
  self.assertLess(np.abs(v[:,2]).max(),1e-5)
 def test_core_union(self):
  grid=np.array([[0.,0.,0.],[10,10,10]])
  phi=np.array([5.,-1.]);q=np.array([[0.,0.,0.]])
  p=self.c.protein_field(phi,grid,q,np.array([1.7]),.5)
  self.assertGreater(p[0],0);self.assertGreater(p[1],0)
 def test_chemical_colors(self):
  c=self.c.lipo_colors(np.array([-1.,0.,1.]))
  self.assertEqual(c.shape,(3,3));self.assertTrue(np.all((c>=0)&(c<=1)))
if __name__=='__main__':unittest.main()
