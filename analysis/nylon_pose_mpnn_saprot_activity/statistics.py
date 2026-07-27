from dataclasses import dataclass
import numpy as np
import pandas as pd
from scipy import stats
from .schema import PRIMARY_PREDICTORS

def holm_adjust(pvalues):
    p=np.asarray(pvalues,dtype=float)
    out=np.full(len(p),np.nan)
    valid=np.flatnonzero(np.isfinite(p))
    if not len(valid): return out
    order=valid[np.argsort(p[valid])]
    running=0.0
    m=len(order)
    for rank,idx in enumerate(order):
        running=max(running,min(1.0,(m-rank)*p[idx]))
        out[idx]=running
    return out

def _corr(x,y,method):
    if len(x)<3 or np.nanstd(x)==0 or np.nanstd(y)==0:
        return np.nan,np.nan
    f=stats.pearsonr if method=="pearson" else stats.spearmanr
    r,p=f(x,y)
    return float(r),float(p)

def _bootstrap_corr(x,y,method,n_boot,seed):
    rng=np.random.default_rng(seed)
    vals=[]
    n=len(x)
    for _ in range(n_boot):
        idx=rng.integers(0,n,n)
        r,_=_corr(x[idx],y[idx],method)
        if np.isfinite(r): vals.append(r)
    if vals:
        lo,hi=np.percentile(vals,[2.5,97.5])
    else:
        lo=hi=np.nan
    return float(lo),float(hi),len(vals)

def _vif(z):
    vals={}
    for i,name in enumerate(z.columns):
        y=z.iloc[:,i].to_numpy(float)
        others=np.delete(z.to_numpy(float),i,axis=1)
        if others.shape[1]==0:
            vals[name]=1.0
            continue
        X=np.column_stack([np.ones(len(y)),others])
        beta=np.linalg.lstsq(X,y,rcond=None)[0]
        resid=y-X@beta
        denom=((y-y.mean())**2).sum()
        r2=1-resid.dot(resid)/denom if denom>0 else 1
        vals[name]=float(np.inf if r2>=1-1e-12 else 1/(1-r2))
    return vals

def _standardized_ols(frame,endpoint,predictors):
    d=frame[[endpoint]+predictors].dropna()
    audit={"n_complete":int(len(d)),"status":"PASS","max_vif":np.nan,
           "condition_number":np.nan,"withhold_reason":""}
    rows=[]
    if len(d)<=len(predictors)+2:
        audit.update(status="WITHHELD",withhold_reason="N_LE_P_PLUS_2")
        return pd.DataFrame(),audit
    z=(d[predictors]-d[predictors].mean())/d[predictors].std(ddof=0)
    y=(d[endpoint]-d[endpoint].mean())/d[endpoint].std(ddof=0)
    if not np.isfinite(z.to_numpy()).all() or not np.isfinite(y).all():
        audit.update(status="WITHHELD",withhold_reason="ZERO_VARIANCE")
        return pd.DataFrame(),audit
    vifs=_vif(z)
    X=np.column_stack([np.ones(len(z)),z.to_numpy()])
    cond=float(np.linalg.cond(X))
    audit["max_vif"]=float(max(vifs.values()))
    audit["condition_number"]=cond
    if audit["max_vif"]>10 or cond>1e8:
        audit.update(status="WITHHELD",withhold_reason="COLLINEARITY_GATE")
        return pd.DataFrame(),audit
    beta=np.linalg.lstsq(X,y,rcond=None)[0]
    resid=y-X@beta
    df=len(y)-X.shape[1]
    sigma2=resid.dot(resid)/df
    cov=sigma2*np.linalg.inv(X.T@X)
    se=np.sqrt(np.diag(cov))
    t=beta/se
    p=2*stats.t.sf(np.abs(t),df)
    for i,name in enumerate(predictors,1):
        rows.append({"predictor":name,"standardized_beta":beta[i],"se":se[i],
                     "t":t[i],"p":p[i],"vif":vifs[name],"n":len(d)})
    return pd.DataFrame(rows),audit

def analyze_endpoint(frame,endpoint,seed=20260727,n_boot=10000,predictors=None):
    predictors=list(predictors or PRIMARY_PREDICTORS)
    corr=[]
    boot=[]
    for pi,predictor in enumerate(predictors):
        d=frame[[predictor,endpoint]].dropna()
        x=d[predictor].to_numpy(float); y=d[endpoint].to_numpy(float)
        row={"endpoint":endpoint,"predictor":predictor,"n":len(d)}
        for mi,method in enumerate(["pearson","spearman"]):
            r,p=_corr(x,y,method)
            row[f"{method}_r" if method=="pearson" else "spearman_rho"]=r
            row[f"{method}_p"]=p
            lo,hi,n_success=_bootstrap_corr(x,y,method,n_boot,seed+1000*pi+mi)
            boot.append({"endpoint":endpoint,"predictor":predictor,"method":method,
                         "ci_low":lo,"ci_high":hi,"n_success":n_success,
                         "n_boot":n_boot,"n":len(d)})
        corr.append(row)
    corr=pd.DataFrame(corr)
    allp=[]
    locations=[]
    for col in ["pearson_p","spearman_p"]:
        for idx,val in corr[col].items():
            allp.append(val); locations.append((idx,col))
    adj=holm_adjust(allp)
    corr["pearson_p_holm"]=np.nan; corr["spearman_p_holm"]=np.nan
    for (idx,col),val in zip(locations,adj):
        corr.loc[idx,col+"_holm"]=val
    corr["p_holm"]=corr[["pearson_p_holm","spearman_p_holm"]].min(axis=1)
    multi,audit=_standardized_ols(frame,endpoint,predictors)
    if len(multi): multi.insert(0,"endpoint",endpoint)
    predcorr=frame[predictors].corr(method="spearman").stack().reset_index()
    predcorr.columns=["predictor_a","predictor_b","spearman_rho"]
    predcorr.insert(0,"endpoint",endpoint)
    return {"correlations":corr,"bootstrap":pd.DataFrame(boot),
            "multivariable":multi,"predictor_correlations":predcorr,
            "model_audit":audit}
