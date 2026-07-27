import numpy as np
from analysis.nylon_pose_mpnn_saprot_activity.clusters import sample_one_per_cluster

def test_cluster_sampler_returns_one_per_cluster():
    assignments={"a":"c1","b":"c1","c":"c2"}
    sampled=sample_one_per_cluster(assignments,np.random.default_rng(7))
    assert len(sampled)==2
    assert set(sampled)<=set(assignments)
    assert len({assignments[x] for x in sampled})==2
