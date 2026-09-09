import unittest,importlib.util,pathlib,numpy as np
class Test(unittest.TestCase):
 def setUp(self):
  p=pathlib.Path(__file__).with_name("tangent_probe_v1.py")
  self.assertTrue(p.exists(),"tangent-probe implementation missing")
  s=importlib.util.spec_from_file_location("tangent",p);self.m=importlib.util.module_from_spec(s);s.loader.exec_module(self.m)
 def env(self,xyz,els,box=(100,100,100)):
  return self.m.old.Env(np.array(xyz,float),np.array(els),np.array(box,float))
 def test_isolated_tangent_and_free_path(self):
  e=self.env([[0,0,0]],["C"],(100000,100000,100));u=self.m.old.fibonacci(64)
  d=self.m.trace(e,np.zeros(3),u,2.,max_path=1000)
  self.assertTrue(d["endpoint"].all());self.assertTrue(d["reached"].all())
 def test_free_endpoint_but_blocked_path(self):
  e=self.env([[0,0,0],[0,0,10]],["C","C"])
  d=self.m.trace(e,np.zeros(3),np.array([[0,0,1.]]),1.)
  self.assertTrue(d["endpoint"][0]);self.assertFalse(d["reached"][0]);self.assertTrue(d["blocked"][0])
 def test_target_neighbors_not_deleted(self):
  e=self.env([[0,0,0],[0,0,1.2]],["C","O"])
  d=self.m.trace(e,np.zeros(3),np.array([[0,0,1.]]),1.)
  self.assertFalse(d["endpoint"][0])
 def test_periodic_segment_collision(self):
  e=self.env([[19.8,0,0]],["C"],(20,20,40))
  hit=self.m.segment_hits(e,np.array([[0.,0,0]]),np.array([[0.,1,0]]),np.array([1.]),.5)
  self.assertTrue(hit[0])
 def test_horizontal_unresolved_not_pass(self):
  e=self.env([[0,0,0]],["C"],(1000,1000,100))
  d=self.m.trace(e,np.zeros(3),np.array([[1.,0,0]]),1.,max_path=10)
  self.assertFalse(d["reached"][0]);self.assertTrue(d["unresolved"][0])
 def test_score_not_density(self):
  self.assertAlmostEqual(self.m.score([1,2,4],[1,1,1]),1)
  self.assertAlmostEqual(self.m.score([1,2,4],[0,0,0]),0)
if __name__=="__main__":unittest.main()
