import hashlib
import pandas as pd
from analysis.nylon_pose_mpnn_saprot_activity.plotting import render_endpoint_figures

def test_figure_audit_matches_rows(tmp_path):
    md5s=["a"*32,"b"*32,"c"*32,"d"*32]
    df=pd.DataFrame({"enzyme":["N1","N2","N3","N4"],"sequence_md5":md5s,
      "pose_strict_nac":[.1,.2,.3,.4],"mpnn_vanilla_mean":[-2,-1.5,-1,-.5],
      "saprot_full_mean":[-3,-2.5,-2,-1.5],"pa66_l2_color_span_px":[1,2,3,4]})
    audit=render_endpoint_figures(df,"pa66_l2_color_span_px",tmp_path)
    expected=hashlib.sha256("\n".join(sorted(md5s)).encode()).hexdigest()
    assert audit["endpoint"]=="pa66_l2_color_span_px"
    assert audit["n_points"]==4
    assert audit["sequence_md5_sha256"]==expected
    assert audit["y_axis_label"].startswith("PA66 L2")
    assert (tmp_path/"pa66_l2_color_span_px_3d.html").exists()
