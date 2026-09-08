import unittest,importlib.util,numpy as np
from pathlib import Path
class TestCoverage(unittest.TestCase):
 def setUp(self):
  p=Path(__file__).with_name('whole_site_contact_v1.py')
  self.assertTrue(p.exists(),'implementation missing')
  spec=importlib.util.spec_from_file_location('whole',p);self.m=importlib.util.module_from_spec(spec);spec.loader.exec_module(self.m)
 def test_targets(self):
  self.assertEqual([self.m.target_count(29,x) for x in [.9,.7,.5]],[27,21,15])
 def test_fixed_subset(self):
  f=np.array([[3.,0,0],[0,3,0],[0,0,3],[3,0,0]])
  sets=self.m.subsets(f,2,np.ones(3))
  self.assertIn((0,3),[tuple(sorted(s)) for s in sets])
 def test_whole_site(self):
  f=np.array([[0,0,0],[0,.5,0],[.1,0,.1]])
  self.assertEqual(self.m.whole_pass(f,.4).tolist(),[True,False,True])
if __name__=='__main__':unittest.main()
