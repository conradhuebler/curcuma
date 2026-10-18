import sys,glob; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase2/scripts')
from common2 import *
seeds=[int(s) for s in sys.argv[1:]] or [101]
SL=os.environ.get('DSVAR')=='sl'; PRE='sl' if SL else ''
def rk(v):
    r=np.empty(len(v)); r[np.argsort(v,kind='stable')]=np.arange(len(v)); return r
def sp(a,b):
    if len(a)<3: return np.nan
    return float(np.corrcoef(rk(a),rk(b))[0,1])
def metrics(E,ref):
    E=np.asarray(E)*HA; ref=np.asarray(ref)*HA; n=len(E)
    o=np.argsort(E,kind='stable'); orf=np.argsort(ref,kind='stable'); t=orf[0]
    dE=ref[:,None]-ref[None,:]; mE=E[:,None]-E[None,:]; iu=np.triu_indices(n,1)
    m=np.abs(dE[iu])>1.0; pa=float(np.mean(np.sign(dE[iu][m])==np.sign(mE[iu][m])))
    k10=max(1,int(round(0.1*n))); top=set(orf[:5].tolist()); topk=set(orf[:k10].tolist())
    low=orf[:max(3,int(round(0.25*n)))]
    c=E-ref
    return dict(spear=sp(E,ref),pairacc=pa,top1=float(o[0]==t),r5=len(top&set(o[:5].tolist()))/5,r10p=len(topk&set(o[:k10].tolist()))/k10,
        regret=float(ref[o[0]]-ref[t]),cmae=float(np.mean(np.abs(c-c.mean()))),spear_low=sp(E[low],ref[low]))
keys=['spear','pairacc','top1','r5','r10p','regret','cmae','spear_low']
out=[];csv=['M,seed,fold,method,n,'+','.join(keys)]
for M in ('M1','M2'):
    res={m:[] for m in ('GFNFF','A','C','S')}; bad=0; devs=[]; nrow=0
    for sd in seeds:
        for f in range(5):
            rows=json.load(open(f'{UC}/{M}/models/{PRE}cv{sd}_f{f}/result.json'))['rows']
            ref=[r['E_ref'] for r in rows]
            for m in res:
                if m=='GFNFF': E=[(r['E_A']-r['corr_A']) if SL else r['E_base_stored'] for r in rows]
                else: E=[r['E_'+m] for r in rows]
                mt=metrics(E,ref); res[m].append(mt)
                csv.append(f'{M},{sd},{f},{m},{len(rows)},'+','.join('%.5f'%mt[k] for k in keys))
            for r in rows:
                nrow+=1
                for m in 'ACS':
                    if r['rc_'+m]!=0 or not np.isfinite(r['E_'+m]) or not np.isfinite(r['corr_'+m]) : bad+=1
                    devs.append(abs(r['E_'+m]-r['corr_'+m]-r['E_base_stored']))
    nf=len(res['A'])
    out.append(f'== {M} CV seeds {seeds}: {nf} folds, {nrow} predictions; failed/NaN sp: {bad}; implicit base - stored: max {max(devs):.2e} Eh, median {np.median(devs):.2e}')
    out.append('%-6s %-11s %-11s %-9s %-11s %-11s %-12s %-11s %-11s'%('method','Spearman','pair-acc','top1','recall@5','recall@10%','regret','cMAE','Spear.low25'))
    for m in res:
        g=lambda k:(np.nanmean([x[k] for x in res[m]]),np.nanstd([x[k] for x in res[m]],ddof=1))
        out.append('%-6s '%m+' '.join('%.2f+-%.2f '%g(k) if k!='top1' else '%.2f      '%g(k)[0] for k in keys))
    for m in ('A','C','S'):
        s=[]
        for k in ('spear','pairacc','r5','r10p','regret','cmae','spear_low'):
            dl=np.array([a[k]-b[k] for a,b in zip(res[m],res['GFNFF'])])
            better=(dl<0).sum() if k in ('regret','cmae') else (dl>0).sum()
            s.append('%s %+.3f+-%.3f (%d/%d better)'%(k,np.nanmean(dl),np.nanstd(dl,ddof=1),better,nf))
        out.append(f' paired {m}-GFNFF n={nf} folds: '+'; '.join(s))
open(f'{UC}/test1_summary{PRE}.txt','w').write('\n'.join(out)+'\n');open(f'{UC}/test1{PRE}.csv','w').write('\n'.join(csv)+'\n');print('\n'.join(out))
