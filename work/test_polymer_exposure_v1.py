import unittest, importlib.util
import numpy as np
from pathlib import Path
class Tests(unittest.TestCase):
 def setUp(self):
  p=Path(__file__).with_name('polymer_exposure_v1.py')
  self.assertTrue(p.exists(),'polymer exposure implementation not yet present')
  s=importlib.util.spec_from_file_location('m',p);self.m=importlib.util.module_from_spec(s);s.loader.exec_module(self.m)
 def test_cavity(self):
  a=np.ones((7,7,7),bool);a[1:6,1:6,1:6]=False;a[3,3,3]=True
  e=self.m.external_components(a);self.assertFalse(e[3,3,3]);self.assertTrue(e[0,0,0])
 def test_periodic(self):
  a=np.zeros((5,5,5),bool);a[0,2,:]=True;a[4,2,2]=True
  e=self.m.external_components(a);self.assertTrue(e[4,2,2])
 def test_frame(self):
  f=self.m.make_frame(np.array([0.,0.,-1.]),np.array([1.,0.,0.]),np.array([0.,1.,0.]))
  np.testing.assert_allclose(f.T@f,np.eye(3),atol=1e-12);self.assertAlmostEqual(np.linalg.det(f),1)
  np.testing.assert_allclose(f[:,2],[0,0,-1])
 def test_isolated_sasa(self):
  env=self.m.Env(np.array([[10.,10.,10.]]),np.array(['C']),np.array([20.,20.,20.]))
  u=self.m.fibonacci(4096);p=np.array([10.,10.,10.])+3.1*u
  self.assertTrue(np.all(env.clearance(p)>=1.4-1e-7))
 def test_bins(self):
  self.assertEqual(self.m.exposure_bin(0),'sampled_zero')
  self.assertEqual(self.m.exposure_bin(.05),'low_0_5pct')
  self.assertEqual(self.m.exposure_bin(.15),'mid_5_15pct')
  self.assertEqual(self.m.support_bin(.8),'80_90')
 def test_periodic_distance(self):
  env=self.m.Env(np.array([[.1,10.,10.]]),np.array(['C']),np.array([20.,20.,20.]))
  self.assertAlmostEqual(float(env.clearance(np.array([[19.9,10.,10.]]))[0]),-1.5,places=7)
if __name__=='__main__': unittest.main()
