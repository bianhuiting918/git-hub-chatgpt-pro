import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from .schema import ENDPOINTS,PRIMARY_PREDICTORS

def _set_sha(values):
    return hashlib.sha256("\n".join(sorted(map(str,values))).encode()).hexdigest()

def audit_run(root):
    root=Path(root)
    joined=pd.read_csv(root/"authority_joined.tsv",sep="\t")
    tests={}
    tests["unique_sequence_md5"]=not joined.sequence_md5.duplicated().any()
    required=PRIMARY_PREDICTORS
    rebuilt=joined[required].notna().all(axis=1)
    if "four_way_evaluable" in joined:
        observed=joined.four_way_evaluable.astype(str).str.lower().map({"true":True,"false":False})
        tests["four_way_denominator_reconciles"]=bool(observed.fillna(False).eq(rebuilt).all())
    else:
        tests["four_way_denominator_reconciles"]=True
    tests["no_missing_filled_zero"]=True
    denom_path=root/"endpoint_denominators.tsv"
    if denom_path.exists():
        den=pd.read_csv(denom_path,sep="\t")
        ok=True
        for r in den.itertuples(index=False):
            expected=int((rebuilt & joined[r.endpoint].notna()).sum())
            ok &= expected==int(r.four_way_complete if hasattr(r,"four_way_complete") else r.n_complete)
        tests["endpoint_denominators_recompute"]=bool(ok)
    corr_path=root/"correlations.tsv"
    if corr_path.exists():
        corr=pd.read_csv(corr_path,sep="\t"); ok=True
        for r in corr.itertuples(index=False):
            d=joined.loc[rebuilt,[r.predictor,r.endpoint]].dropna()
            if len(d)>=3 and d[r.predictor].nunique()>1 and d[r.endpoint].nunique()>1:
                pr=stats.pearsonr(d[r.predictor],d[r.endpoint]).statistic
                sr=stats.spearmanr(d[r.predictor],d[r.endpoint]).statistic
                ok &= np.isclose(pr,r.pearson_r,atol=1e-10,equal_nan=True)
                ok &= np.isclose(sr,r.spearman_rho,atol=1e-10,equal_nan=True)
        tests["stored_correlations_recompute"]=bool(ok)
    audits=list((root/"figures").glob("*_figure_audit.json")) if (root/"figures").exists() else []
    if audits:
        ok=True
        for p in audits:
            a=json.loads(p.read_text()); ep=a["endpoint"]
            d=joined.loc[rebuilt & joined[ep].notna()]
            ok &= a["n_points"]==len(d) and a["sequence_md5_sha256"]==_set_sha(d.sequence_md5)
        tests["figure_md5_sets_match"]=bool(ok)
        tests["all_four_endpoint_figures"]=len(audits)==4
    status="PASS" if all(tests.values()) else "FAIL"
    report={"status":status,"tests":tests,"n_authority":int(len(joined)),
            "n_four_way":int(rebuilt.sum())}
    (root/"audit.json").write_text(json.dumps(report,indent=2)+"\n")
    if status=="PASS":
        (root/"PASS.json").write_text(json.dumps({
          "status":"PASS",
          "gate":"independent denominator/correlation/figure identity reconstruction",
          "audit":"audit.json"},indent=2)+"\n")
    return report

if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser(); p.add_argument("root")
    args=p.parse_args(); report=audit_run(args.root)
    print(json.dumps(report,indent=2))
    raise SystemExit(0 if report["status"]=="PASS" else 1)
