import argparse
import hashlib
import json
import os
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
from .assemble import assemble_authority
from .statistics import analyze_endpoint,holm_adjust
from .clusters import run_mmseqs_clusters,cluster_sensitivity
from .plotting import render_endpoint_figures
from .schema import ENDPOINTS,PRIMARY_PREDICTORS,SENSITIVITY_PREDICTORS

def sha256(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""): h.update(block)
    return h.hexdigest()

def _write(df,path):
    df.to_csv(path,sep="\t",index=False,na_rep="")

def _global_holm(corr):
    locations=[]; values=[]
    for col in ["pearson_p","spearman_p"]:
        for idx,v in corr[col].items():
            locations.append((idx,col)); values.append(v)
    adj=holm_adjust(values)
    corr["pearson_p_holm_global"]=np.nan; corr["spearman_p_holm_global"]=np.nan
    for (idx,col),v in zip(locations,adj):
        corr.loc[idx,col+"_holm_global"]=v
    return corr

def main(argv=None):
    ap=argparse.ArgumentParser()
    for name in ["activity","pose","pose_strata","pose_not_evaluated","mpnn","saprot","output"]:
        ap.add_argument("--"+name.replace("_","-"),required=True)
    ap.add_argument("--n-boot",type=int,default=10000)
    ap.add_argument("--seed",type=int,default=20260727)
    ap.add_argument("--mmseqs",default=None)
    args=ap.parse_args(argv)
    out=Path(args.output); out.mkdir(parents=True,exist_ok=True)
    figdir=out/"figures"; figdir.mkdir(exist_ok=True)
    inputs={k:Path(getattr(args,k)) for k in ["activity","pose","pose_strata",
                                             "pose_not_evaluated","mpnn","saprot"]}
    manifest=[]
    for name,path in inputs.items():
        manifest.append({"input_name":name,"path":str(path.resolve()),"bytes":path.stat().st_size,
                         "sha256":sha256(path)})
    _write(pd.DataFrame(manifest),out/"input_sha256.tsv")
    activity=pd.read_csv(inputs["activity"],sep="\t",dtype={"sequence_md5":str})
    if "table_s1_id" in activity:
        activity=activity[activity.table_s1_id.astype(str).str.fullmatch(r"Nyl\d{2}")].copy()
        activity["enzyme"]=activity.table_s1_id
    pose=pd.read_csv(inputs["pose"],sep="\t",dtype={"sequence_md5":str})
    strata=pd.read_csv(inputs["pose_strata"],sep="\t",dtype={"sequence_md5":str})
    ne=pd.read_csv(inputs["pose_not_evaluated"],sep="\t",dtype={"sequence_md5":str})
    ne["reason"]=np.where(ne.status.astype(str).str.startswith("NOT_EVALUATED_"),
                          ne.status.astype(str),"NOT_EVALUATED_POSE_"+ne.status.astype(str))
    mpnn=pd.read_csv(inputs["mpnn"],sep="\t",dtype=str,low_memory=False)
    saprot=pd.read_csv(inputs["saprot"],sep="\t",dtype=str,low_memory=False)
    joined,excluded,summary=assemble_authority(activity,pose,strata,mpnn,saprot,ne)
    if len(activity)!=95:
        raise RuntimeError(f"activity Nyl01-Nyl95 authority expected 95, got {len(activity)}")
    _write(joined,out/"authority_joined.tsv")
    _write(excluded,out/"not_evaluated.tsv")
    denom=[]
    corrs=[]; boots=[]; multis=[]; predcorr=[]; model_audits={}; sens=[]
    for ei,endpoint in enumerate(ENDPOINTS):
        d=joined[joined.four_way_evaluable & joined[endpoint].notna()].copy()
        denom.append({"endpoint":endpoint,"activity_authority":len(joined),
                      "pose_snapshot_evaluated":int(joined.pose_strict_nac.notna().sum()),
                      "four_way_complete":len(d),
                      "not_complete":len(joined)-len(d)})
        result=analyze_endpoint(d,endpoint,seed=args.seed+ei*10000,n_boot=args.n_boot)
        corrs.append(result["correlations"]); boots.append(result["bootstrap"])
        if len(result["multivariable"]): multis.append(result["multivariable"])
        predcorr.append(result["predictor_correlations"])
        model_audits[endpoint]=result["model_audit"]
        sres=analyze_endpoint(d,endpoint,seed=args.seed+50000+ei*10000,
                              n_boot=args.n_boot,predictors=SENSITIVITY_PREDICTORS)
        sens.append(sres["correlations"])
        render_endpoint_figures(d,endpoint,figdir)
    corr=_global_holm(pd.concat(corrs,ignore_index=True))
    senscorr=_global_holm(pd.concat(sens,ignore_index=True))
    _write(pd.DataFrame(denom),out/"endpoint_denominators.tsv")
    _write(corr,out/"correlations.tsv"); _write(pd.concat(boots,ignore_index=True),out/"bootstrap_ci.tsv")
    _write(pd.concat(multis,ignore_index=True) if multis else pd.DataFrame(),out/"multivariable_models.tsv")
    _write(pd.concat(predcorr,ignore_index=True),out/"predictor_correlations.tsv")
    _write(senscorr,out/"sensitivity_correlations.tsv")
    (out/"model_audit.json").write_text(json.dumps(model_audits,indent=2)+"\n")
    fasta=out/"nyl95_exact_sequences.fasta"
    with fasta.open("w") as f:
        for r in joined.itertuples(index=False):
            f.write(f">{r.sequence_md5}\n{r.sequence}\n")
    assignments,cluster_audit=run_mmseqs_clusters(fasta,out/"clusters",args.mmseqs)
    if assignments:
        cadf=pd.DataFrame([{"sequence_md5":k,"cluster_id":v} for k,v in sorted(assignments.items())])
        _write(cadf,out/"cluster_assignments.tsv")
        cs=[]
        for ei,endpoint in enumerate(ENDPOINTS):
            d=joined[joined.four_way_evaluable & joined[endpoint].notna()]
            cs.append(cluster_sensitivity(d,endpoint,PRIMARY_PREDICTORS,assignments,
                                           n_boot=args.n_boot,seed=args.seed+70000+ei*10000))
        _write(pd.concat(cs,ignore_index=True),out/"cluster_sensitivity.tsv")
    (out/"cluster_audit.json").write_text(json.dumps(cluster_audit,indent=2)+"\n")
    summary.update({"activity_authority_expected":95,
      "pose_snapshot_not_evaluated":int(len(joined)-joined.pose_strict_nac.notna().sum()),
      "pose_snapshot_provenance":"2026-07-22 frozen Dell-derived snapshot; live Dell refresh blocked before SSH handshake",
      "experimental_readout":"raster color-span proxy, not absolute activity",
      "saProt_interpretation":"structure-conditional compatibility, not Tm/dG/activity",
      "proteinmpnn_interpretation":"sequence-backbone compatibility, not Tm/dG/activity"})
    (out/"denominator_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    return 0

if __name__=="__main__": raise SystemExit(main())
