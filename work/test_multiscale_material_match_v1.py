import unittest, importlib.util, pathlib
import numpy as np
class Tests(unittest.TestCase):
 def setUp(self):
  p=pathlib.Path(__file__).with_name("multiscale_material_match_v1.py")
  self.assertTrue(p.exists(),"multiscale material shape metric implementation missing")
  spec=importlib.util.spec_from_file_location("ms",p);self.m=importlib.util.module_from_spec(spec);spec.loader.exec_module(self.m)
 def plane(self,shift=0):
  v=np.array([[x,y,shift] for x in range(4) for y in range(4)],float)
  f=[]
  for i in range(3):
   for j in range(3):
    a=i*4+j;f.extend([[a,a+4,a+1],[a+1,a+4,a+5]])
  return v,np.array(f)
 def test_identity(self):
  a=self.m.mesh_rep(*self.plane());self.assertLess(self.m.distance(a,a)[0],1e-10)
 def test_two_angstrom_shift(self):
  a=self.m.mesh_rep(*self.plane());b=self.m.mesh_rep(*self.plane(2))
  self.assertAlmostEqual(self.m.distance(a,b)[0],2,places=6)
 def test_symmetric_extra_patch(self):
  v,f=self.plane();a=self.m.mesh_rep(v,f)
  b=self.m.mesh_rep(np.vstack([v,v+[20,0,0]]),np.vstack([f,f+len(v)]))
  self.assertGreater(self.m.distance(a,b)[0],10)
  self.assertAlmostEqual(self.m.distance(a,b)[0],self.m.distance(b,a)[0])
 def test_ceil_includes_unknown_denominator(self):
  self.assertEqual(self.m.required_tolerance([0,1,2,float("inf")],.7),2)
  self.assertTrue(np.isinf(self.m.required_tolerance([0,1,2,float("inf")],.9)))
 def test_empty(self):
  a=self.m.mesh_rep(*self.plane())
  b=self.m.mesh_rep(np.empty((0,3)),np.empty((0,3),int))
  self.assertTrue(np.isinf(self.m.distance(a,b)[0]))
if __name__=="__main__":unittest.main()
