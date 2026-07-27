import numpy as np
import pandas as pd
from analysis.nylon_pose_mpnn_saprot_activity.statistics import analyze_endpoint

def test_monotone_statistics_and_holm():
    x=np.arange(1,13,dtype=float)
    df=pd.DataFrame({"pose_strict_nac":x,"mpnn_vanilla_mean":x*2,"saprot_full_mean":x*3,
                     "pa66_l1_color_span_px":x})
    out=analyze_endpoint(df,"pa66_l1_color_span_px",seed=20260727,n_boot=200)
    row=out["correlations"].query("predictor == 'pose_strict_nac'").iloc[0]
    assert np.isclose(row.spearman_rho,1)
    assert out["correlations"].p_holm.between(0,1).all()
    assert out["bootstrap"].n_success.min()>=180
    assert out["model_audit"]["n_complete"]==12
