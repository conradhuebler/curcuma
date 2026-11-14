# fit(dsname, tag, seed, params(list of str), idx=None) -> dict ; model dir UC4/fits/<ds>/<tag>_s<seed>
import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase4/scripts")
from common4 import *
import json,subprocess,time,shutil
def parse(txt):
    t=ANSI.sub('',txt); r={}
    f=lambda p,s=t:[tuple(map(float,m)) for m in re.findall(p,s)]
    m=re.search(r'held-out set \(energy model.*?\((\d+) structures\)(.*?)held-out set \(joint model.*?\((\d+) structures\)(.*)',t,re.S)
    if not m: return r
    for role,blk in (('E',m.group(2)),('J',m.group(4))):
        a=re.search(r'relative energy MAE : model\s+([\d.]+) kcal/mol\s+base\s+([\d.]+)',blk)
        b=re.search(r'gradient RMSE\s+: model\s+([\d.]+) Eh/A\s+base\s+([\d.]+)',blk)
        c=re.search(r'rank correlation\s+: model\s+(-?[\d.]+)\s+base\s+(-?[\d.]+)',blk)
        r[role]=dict(mae=float(a.group(1)),mae_b=float(a.group(2)),grms=float(b.group(1)),grms_b=float(b.group(2)),sp=float(c.group(1)),sp_b=float(c.group(2)))
    r['n_test']=int(m.group(3))
    tr=re.search(r'training set \(joint model\).*?MAE : model\s+([\d.]+).*?RMSE\s+: model\s+([\d.]+)',t,re.S)
    if tr: r['J']['tr_mae']=float(tr.group(1)); r['J']['tr_grms']=float(tr.group(2))
    s=re.search(r'\(joint model\): fitted (\d+) centres from (\d+) environments, (\d+) features, sigma = ([\d.eE+-]+)',t)
    if s: r['ncen']=int(s.group(1)); r['nenv']=int(s.group(2)); r['nfeat']=int(s.group(3)); r['sigma']=float(s.group(4))
    g=re.search(r'(\d+) energy rows, (\d+) gradient rows',t)
    if g: r['erows']=int(g.group(1)); r['grows']=int(g.group(2))
    return r
def fit(dsname,tag,seed,params,idx=None,test_fraction='0.2',force=False,wdroot=None):
    wd=f'{wdroot or UC4+"/fits"}/{os.path.basename(dsname).replace(".json","")}/{tag}_s{seed}'
    if os.path.exists(wd+'/res.json') and not force: return json.load(open(wd+'/res.json'))
    shutil.rmtree(wd,ignore_errors=True)
    d=json.load(open(ds(dsname) if dsname in ('E1','M1','M2') else dsname))
    S=d['structures']; idx=list(range(len(S))) if idx is None else idx
    sub={k:v for k,v in d.items() if k!='structures'}; sub['structures']=[S[i] for i in idx]; sub['n_structures']=len(idx)
    os.makedirs(wd); json.dump(sub,open(wd+'/train.dataset.json','w'))
    # anchor molecule: first frame of the dataset (same for every fit of this molecule)
    open(wd+'/train.xyz','w').write(xyz_text(S,[idx[0]]))
    t=time.time()
    with open(wd+'/log.txt','w') as lf:
        rc=subprocess.run([BIN4,'-activelearn','train.xyz','-method','gfnff','-ref_method','gfn2','-activelearn.dataset',wd+'/train.dataset.json','-test_fraction',test_fraction,'-activelearn.seed',str(seed),'-threads','8','-no_bmt']+params,cwd=wd,stdout=lf,stderr=subprocess.STDOUT).returncode
    r=parse(open(wd+'/log.txt').read()) if test_fraction!='0' else {}
    r.update(ds=dsname,tag=tag,seed=seed,params=params,rc=rc,wall=time.time()-t)
    os.remove(wd+'/train.dataset.json')
    json.dump(r,open(wd+'/res.json','w'))
    return r
