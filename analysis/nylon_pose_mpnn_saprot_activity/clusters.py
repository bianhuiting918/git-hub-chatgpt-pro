import shutil
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

def sample_one_per_cluster(assignments,rng):
    groups={}
    for item,cluster in assignments.items():
        groups.setdefault(cluster,[]).append(item)
    return [rng.choice(sorted(items)) for cluster,items in sorted(groups.items())]

def run_mmseqs_clusters(fasta_path,out_dir,mmseqs=None):
    out_dir=Path(out_dir); out_dir.mkdir(parents=True,exist_ok=True)
    binary=mmseqs or shutil.which("mmseqs")
    if not binary:
        return {},{"status":"NOT_EVALUATED_CLUSTER_TOOL","reason":"MMSEQS2_NOT_FOUND"}
    prefix=out_dir/"nyl95_mmseqs30"
    tmp=out_dir/"tmp"
    cmd=[binary,"easy-cluster",str(fasta_path),str(prefix),str(tmp),
         "--min-seq-id","0.30","-c","0.80","--cov-mode","0","--cluster-mode","2"]
    cp=subprocess.run(cmd,text=True,capture_output=True)
    (out_dir/"mmseqs.stdout.log").write_text(cp.stdout)
    (out_dir/"mmseqs.stderr.log").write_text(cp.stderr)
    if cp.returncode:
        return {},{"status":"NOT_EVALUATED_CLUSTER_TOOL","reason":f"MMSEQS2_EXIT_{cp.returncode}"}
    tsv=Path(str(prefix)+"_cluster.tsv")
    if not tsv.exists():
        return {},{"status":"NOT_EVALUATED_CLUSTER_TOOL","reason":"MMSEQS2_CLUSTER_TSV_MISSING"}
    assignments={}
    for line in tsv.read_text().splitlines():
        rep,member=line.split("\t")[:2]
        assignments[member]=rep
    return assignments,{"status":"PASS","n_clusters":len(set(assignments.values())),
                        "command":" ".join(cmd)}

def cluster_sensitivity(frame,endpoint,predictors,assignments,n_boot=10000,seed=20260727):
    if not assignments:
        return pd.DataFrame()
    rng=np.random.default_rng(seed)
    rows=[]
    base=frame.set_index("sequence_md5",drop=False)
    for predictor in predictors:
        vals_p=[]; vals_s=[]
        for _ in range(n_boot):
            chosen=sample_one_per_cluster(assignments,rng)
            chosen=[x for x in chosen if x in base.index]
            d=base.loc[chosen,[predictor,endpoint]].dropna()
            if len(d)<3 or d[predictor].nunique()<2 or d[endpoint].nunique()<2:
                continue
            vals_p.append(stats.pearsonr(d[predictor],d[endpoint]).statistic)
            vals_s.append(stats.spearmanr(d[predictor],d[endpoint]).statistic)
        for method,vals in [("pearson",vals_p),("spearman",vals_s)]:
            lo,med,hi=(np.percentile(vals,[2.5,50,97.5]) if vals else [np.nan]*3)
            rows.append({"endpoint":endpoint,"predictor":predictor,"method":method,
                         "median_r":med,"ci_low":lo,"ci_high":hi,
                         "n_success":len(vals),"n_boot":n_boot})
    return pd.DataFrame(rows)
