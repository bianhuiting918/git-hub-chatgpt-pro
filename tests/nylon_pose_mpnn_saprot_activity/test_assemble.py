import numpy as np
import pandas as pd
from analysis.nylon_pose_mpnn_saprot_activity.assemble import assemble_authority

def test_exact_md5_and_no_zero_fill():
    activity=pd.DataFrame({
      "enzyme":["Nyl01","Nyl02","Nyl03"],"sequence_md5":["a"*32,"b"*32,"c"*32],
      "sequence":["ACDE","ACDF","ACDG"],"pa66_l1_color_span_px":[1.,0.,np.nan],
      "pa66_l2_color_span_px":[3.,0.,2.],"pa6_max_color_span_px":[2.,1.,4.]})
    pose=pd.DataFrame({"enzyme":["Nyl01","Nyl02"],"sequence_md5":["a"*32,"b"*32],
      "strict_nac_rate":[.1,.0],"retention_rate_loose_gate":[.2,.0]})
    strata=pd.DataFrame({"sequence_md5":["a"*32,"a"*32],"stratum":["weight=5","weight=10"],"retention_rate":[.11,.12]})
    mpnn=pd.DataFrame({"family":["Nylonase","Nylonase"],"canonical_sequence_md5":["a"*32,"b"*32],
      "scoring_sequence_md5":["a"*32,"x"*32],"scoring_sequence_scope":["CANONICAL_EXACT"]*2,
      "proteinmpnn_status":["PASS"]*2,"vanilla_mean_log_likelihood":[-1.,-2.],
      "soluble_mean_log_likelihood":[-1.1,-2.1]})
    saprot=pd.DataFrame({"family":["Nylonase"],"canonical_sequence_md5":["a"*32],
      "scoring_sequence_md5":["a"*32],"scoring_sequence_scope":["CANONICAL_EXACT"],
      "status_FULL_PROTEIN":["PASS"],"mean_log_likelihood_FULL_PROTEIN":[-1.5]})
    joined, excluded, summary=assemble_authority(activity,pose,strata,mpnn,saprot)
    assert summary["activity_authority_unique_md5"]==3
    assert summary["four_way_intersection"]==1
    assert joined.loc[joined.sequence_md5=="a"*32,"l1_fraction"].iat[0]==.25
    assert np.isnan(joined.loc[joined.sequence_md5=="b"*32,"l1_fraction"].iat[0])
    assert set(excluded.reason).issubset({x for x in excluded.reason if x.startswith("NOT_EVALUATED_")})
    assert joined.loc[joined.sequence_md5=="b"*32,"mpnn_vanilla_mean"].isna().all()
