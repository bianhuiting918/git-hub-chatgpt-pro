import json
import pandas as pd
from analysis.nylon_pose_mpnn_saprot_activity.audit_run import audit_run

def test_audit_reconstructs_small_run(tmp_path):
    joined=pd.DataFrame({"enzyme":["N1","N2"],"sequence_md5":["a"*32,"b"*32],
      "pose_strict_nac":[.1,.2],"mpnn_vanilla_mean":[-1.,-2.],"saprot_full_mean":[-1.5,-2.5],
      "pa66_l1_color_span_px":[1.,2.],"pa66_l2_color_span_px":[2.,3.],
      "pa6_max_color_span_px":[3.,4.],"l1_fraction":[1/3,2/5],
      "four_way_evaluable":[True,True]})
    joined.to_csv(tmp_path/"authority_joined.tsv",sep="\t",index=False)
    pd.DataFrame({"endpoint":["pa66_l1_color_span_px"],"n_complete":[2]}).to_csv(tmp_path/"endpoint_denominators.tsv",sep="\t",index=False)
    report=audit_run(tmp_path)
    assert report["tests"]["no_missing_filled_zero"]
    assert report["tests"]["four_way_denominator_reconciles"]
