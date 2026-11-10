# usage: cv_an.py tag DS mode seed [seed...]
import sys; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase3/scripts')
from common3 import *
tag,DS,mode=sys.argv[1:4]; seeds=[int(s) for s in sys.argv[4:]]
res={m:[] for m in ('GFNFF','A','C','S')}; bad=0; devs=[]; nrow=0; fr=[]
for sd in seeds:
    for f in range(5):
        J=json.load(open(f'{UC}/{DS}/models/{mode}{sd}_f{f}/result.json')); rows=J['rows']
        ref=[r['E_ref'] for r in rows]
        for m in res:
            E=[r['E_base_stored'] for r in rows] if m=='GFNFF' else [r['E_'+m] for r in rows]
            res[m].append(metrics(E,ref))
        for r in rows:
            nrow+=1
            for m in 'ACS':
                if r['rc_'+m]!=0 or not np.isfinite(r['E_'+m]) or not np.isfinite(r['corr_'+m]): bad+=1
                devs.append(abs(r['E_'+m]-r['corr_'+m]-r['E_base_stored']))
        S=load(DS)['structures']
        for (i,e,rc) in J['fresh']: fr.append(abs(e-S[i]['energy_base']))
nf=len(res['A'])
out=[f'== {DS} {mode} CV seeds {seeds}: {nf} folds, {nrow} predictions; failed/NaN sp {bad}; implicit base - fresh label: max {max(devs):.1e} Eh; fresh gfnff sp vs label (n={len(fr)}): max {max(fr):.1e} Eh']
out.append('%-6s %-11s %-11s %-9s %-11s %-11s %-12s %-11s %-11s'%('method','Spearman','pair-acc','top1','recall@5','recall@10%','regret','cMAE','Spear.low25'))
for m in res:
    g=lambda k:(np.nanmean([x[k] for x in res[m]]),np.nanstd([x[k] for x in res[m]],ddof=1))
    out.append('%-6s '%m+' '.join('%.2f+-%.2f '%g(k) if k!='top1' else '%.2f      '%g(k)[0] for k in KEYS))
for m in ('A','C','S'):
    s=[]
    for k in ('spear','pairacc','r5','r10p','regret','cmae','spear_low'):
        dl=np.array([a[k]-b[k] for a,b in zip(res[m],res['GFNFF'])])
        better=(dl<0).sum() if k in ('regret','cmae') else (dl>0).sum()
        s.append('%s %+.3f+-%.3f (%d/%d)'%(k,np.nanmean(dl),np.nanstd(dl,ddof=1),better,nf))
    out.append(f' paired {m}-GFNFF n={nf}: '+'; '.join(s))
open(f'{UC}/{tag}_summary.txt','w').write('\n'.join(out)+'\n')
csv=['fold_idx,method,'+','.join(KEYS)]
for m in res:
    for j,x in enumerate(res[m]): csv.append(f'{j},{m},'+','.join('%.5f'%x[k] for k in KEYS))
open(f'{UC}/{tag}.csv','w').write('\n'.join(csv)+'\n'); print('\n'.join(out))
