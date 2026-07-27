import re
import numpy as np
import pandas as pd

MD5_RE=re.compile(r"^[0-9a-f]{32}$")

def _require(df, cols, name):
    missing=[c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{name} missing columns: {missing}")

def _unique(df, key, name):
    dup=df[key].astype(str).duplicated(keep=False)
    if dup.any():
        vals=sorted(df.loc[dup,key].astype(str).unique())[:5]
        raise ValueError(f"{name} duplicate {key}: {vals}")

def _numeric(df, cols):
    for c in cols:
        if c in df:
            df[c]=pd.to_numeric(df[c],errors="coerce")
    return df

def _valid_md5(series):
    return series.astype(str).str.lower().str.match(MD5_RE)

def valid_mpnn(rows):
    _require(rows,["family","canonical_sequence_md5","scoring_sequence_md5",
                   "scoring_sequence_scope","proteinmpnn_status",
                   "vanilla_mean_log_likelihood","soluble_mean_log_likelihood"],"ProteinMPNN")
    x=rows.copy()
    mask=(x.family.eq("Nylonase") & x.proteinmpnn_status.eq("PASS")
          & x.scoring_sequence_scope.eq("CANONICAL_EXACT")
          & x.canonical_sequence_md5.eq(x.scoring_sequence_md5))
    x=x.loc[mask,["canonical_sequence_md5","vanilla_mean_log_likelihood",
                   "soluble_mean_log_likelihood"]].rename(columns={
        "canonical_sequence_md5":"sequence_md5",
        "vanilla_mean_log_likelihood":"mpnn_vanilla_mean",
        "soluble_mean_log_likelihood":"mpnn_soluble_mean"})
    x=_numeric(x,["mpnn_vanilla_mean","mpnn_soluble_mean"])
    x=x[np.isfinite(x.mpnn_vanilla_mean)]
    _unique(x,"sequence_md5","valid ProteinMPNN")
    return x

def valid_saprot(rows):
    _require(rows,["family","canonical_sequence_md5","scoring_sequence_md5",
                   "scoring_sequence_scope","status_FULL_PROTEIN",
                   "mean_log_likelihood_FULL_PROTEIN"],"SaProt")
    x=rows.copy()
    mask=(x.family.eq("Nylonase") & x.status_FULL_PROTEIN.eq("PASS")
          & x.scoring_sequence_scope.eq("CANONICAL_EXACT")
          & x.canonical_sequence_md5.eq(x.scoring_sequence_md5))
    x=x.loc[mask,["canonical_sequence_md5","mean_log_likelihood_FULL_PROTEIN"]].rename(
        columns={"canonical_sequence_md5":"sequence_md5",
                 "mean_log_likelihood_FULL_PROTEIN":"saprot_full_mean"})
    x=_numeric(x,["saprot_full_mean"])
    x=x[np.isfinite(x.saprot_full_mean)]
    _unique(x,"sequence_md5","valid SaProt")
    return x

def assemble_authority(activity, pose, pose_strata, mpnn, saprot, pose_not_evaluated=None):
    _require(activity,["sequence_md5","sequence","pa66_l1_color_span_px",
             "pa66_l2_color_span_px","pa6_max_color_span_px"],"activity")
    _require(pose,["sequence_md5","strict_nac_rate","retention_rate_loose_gate"],"pose")
    _require(pose_strata,["sequence_md5","stratum","retention_rate"],"pose strata")
    a=activity.copy()
    if "enzyme" not in a and "table_s1_id" in a:
        a["enzyme"]=a["table_s1_id"]
    _require(a,["enzyme"],"activity")
    a["sequence_md5"]=a.sequence_md5.astype(str).str.lower()
    if not _valid_md5(a.sequence_md5).all():
        raise ValueError("activity contains invalid sequence_md5")
    _unique(a,"sequence_md5","activity authority")
    p=pose.copy()
    p["sequence_md5"]=p.sequence_md5.astype(str).str.lower()
    _unique(p,"sequence_md5","pose")
    p=p.rename(columns={"strict_nac_rate":"pose_strict_nac",
                        "retention_rate_loose_gate":"pose_loose_gate"})
    p=_numeric(p,["pose_strict_nac","pose_loose_gate"])
    s=pose_strata.copy()
    s["sequence_md5"]=s.sequence_md5.astype(str).str.lower()
    s=_numeric(s,["retention_rate"])
    piv=s[s.stratum.isin(["weight=5","weight=10"])].pivot(
        index="sequence_md5",columns="stratum",values="retention_rate").reset_index()
    piv=piv.rename(columns={"weight=5":"pose_weight5_loose","weight=10":"pose_weight10_loose"})
    p=p.merge(piv,on="sequence_md5",how="left",validate="one_to_one")
    m=valid_mpnn(mpnn)
    q=valid_saprot(saprot)
    cols=["enzyme","sequence_md5","sequence","pa66_l1_color_span_px",
          "pa66_l2_color_span_px","pa6_max_color_span_px"]
    for optional in ["pa66_l1_figure_state","pa66_l2_figure_state","pa6_figure_state"]:
        if optional in a: cols.append(optional)
    j=a[cols].merge(p,on="sequence_md5",how="left",validate="one_to_one")
    j=j.merge(m,on="sequence_md5",how="left",validate="one_to_one")
    j=j.merge(q,on="sequence_md5",how="left",validate="one_to_one")
    j=_numeric(j,["pa66_l1_color_span_px","pa66_l2_color_span_px",
                  "pa6_max_color_span_px"])
    total=j.pa66_l1_color_span_px+j.pa66_l2_color_span_px
    j["l1_fraction"]=np.where(total>0,j.pa66_l1_color_span_px/total,np.nan)
    required=["pose_strict_nac","mpnn_vanilla_mean","saprot_full_mean"]
    j["four_way_evaluable"]=j[required].notna().all(axis=1)
    reason_map={}
    if pose_not_evaluated is not None and len(pose_not_evaluated):
        ne=pose_not_evaluated.copy()
        if "sequence_md5" in ne and "reason" in ne:
            reason_map=dict(zip(ne.sequence_md5.astype(str).str.lower(),ne.reason.astype(str)))
    excluded=[]
    for row in j.itertuples(index=False):
        reasons=[]
        if pd.isna(row.pose_strict_nac):
            reasons.append(reason_map.get(row.sequence_md5,"NOT_EVALUATED_POSE"))
        if pd.isna(row.mpnn_vanilla_mean): reasons.append("NOT_EVALUATED_PROTEINMPNN")
        if pd.isna(row.saprot_full_mean): reasons.append("NOT_EVALUATED_SAPROT_FULL_PROTEIN")
        if reasons:
            excluded.append({"enzyme":row.enzyme,"sequence_md5":row.sequence_md5,
                             "reason":";".join(reasons)})
    excluded=pd.DataFrame(excluded,columns=["enzyme","sequence_md5","reason"])
    summary={
      "activity_authority_unique_md5":int(len(a)),
      "pose_evaluated_unique_md5":int(len(p)),
      "mpnn_exact_pass_unique_md5":int(len(m)),
      "saprot_exact_full_pass_unique_md5":int(len(q)),
      "four_way_intersection":int(j.four_way_evaluable.sum()),
      "not_four_way_evaluable":int((~j.four_way_evaluable).sum()),
    }
    return j,excluded,summary
