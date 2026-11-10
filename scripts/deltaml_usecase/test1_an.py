import sys,glob,json; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase/scripts')
from common import *
rows=[]
for f in sorted(glob.glob(UC+'/cv/s*_f*/result.json')):
    rows+=json.load(open(f))['rows']
M=['base','A','B','C']; NAME={'base':'GFN-FF','A':'A joint-def','B':'B joint-best','C':'C energy'}
def rk(x):
    r=np.empty(len(x));r[np.argsort(x)]=np.arange(len(x));return r
def metrics(er,e):
    N=len(er); t5=set(np.argsort(er)[:5]); o=np.argsort(e)
    rank=int(np.where(o==er.argmin())[0][0])+1
    return dict(rank=rank,top1=float(e.argmin()==er.argmin()),r5=len(t5&set(o[:5]))/5,r10=len(t5&set(o[:10]))/5,
      regret=(er[e.argmin()]-er.min())*HA,spear=float(np.corrcoef(rk(er),rk(e))[0,1]))
bad=[r for r in rows if any(not np.isfinite(r['E_'+m]) or r['rc_'+m]!=0 for m in 'ABC')]
out=["n rows %d, failed/NaN sp: %d"%(len(rows),len(bad))]
pools={'fold':{}, 'pooled':{}}
for r in rows:
    pools['fold'].setdefault((r['seed'],r['fold']),[]).append(r)
    pools['pooled'].setdefault(r['seed'],[]).append(r)
keys=['rank','top1','r5','r10','regret','spear']
csv=['pool,key,method,'+','.join(keys)]
for pn,P in pools.items():
    res={m:{k:[] for k in keys} for m in M}
    Ns=[]
    for key,rr in sorted(P.items()):
        er=np.array([r['E_ref'] for r in rr]); Ns.append(len(rr))
        for m in M:
            e=np.array([r['E_'+m] for r in rr]); mt=metrics(er,e)
            for k in keys: res[m][k].append(mt[k])
            csv.append('%s,%s,%s,'%(pn,key,m)+','.join('%.4f'%mt[k] for k in keys))
    n=len(Ns); N=np.mean(Ns)
    out.append("== pool %s: n=%d pools of size ~%d; random Top-1 %.4f, random Recall@5 %.3f, @10 %.3f"%(pn,n,N,1/N,5/N,10/N))
    out.append("%-13s %-12s %-12s %-12s %-12s %-14s %-12s"%('method','rank(mean/med)','top1','R@5','5in10','regret','spearman'))
    for m in M:
        r=res[m]
        out.append("%-13s %5.1f/%-6.0f %5.2f      %.2f+-%.2f  %.2f+-%.2f  %.2f+-%.2f   %.2f+-%.2f"%(NAME[m],np.mean(r['rank']),np.median(r['rank']),np.mean(r['top1']),np.mean(r['r5']),np.std(r['r5'],ddof=1),np.mean(r['r10']),np.std(r['r10'],ddof=1),np.mean(r['regret']),np.std(r['regret'],ddof=1),np.mean(r['spear']),np.std(r['spear'],ddof=1)))
    for m in M[1:]:
        s=[]
        for k in ['spear','regret','r5','r10']:
            dlt=np.array(res[m][k])-np.array(res['base'][k]); s.append("%s %+.3f+-%.3f(%d/%d better)"%(k,dlt.mean(),dlt.std(ddof=1),(dlt>0).sum() if k!='regret' else (dlt<0).sum(),n))
        out.append(" paired %s - GFN-FF, n=%d: "%(NAME[m],n)+'; '.join(s))
open(UC+'/test1_summary.txt','w').write('\n'.join(out)+'\n'); open(UC+'/test1.csv','w').write('\n'.join(csv)+'\n')
print('\n'.join(out))
