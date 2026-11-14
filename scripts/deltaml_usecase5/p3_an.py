import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase5/scripts")
from common5 import *
import json,csv
SPL={1:([1,2,3,4],[5,6]),2:([3,4,5,6],[1,2]),3:([1,2,5,6],[3,4])}
MODS=['FF','A','C','S','B2']
rng=np.random.RandomState(0)
out=open(UC5+'/a2_summary.txt','w'); csvf=open(UC5+'/a2_metrics.csv','w'); cw=csv.writer(csvf)
cw.writerow(['mol','split','model','n','spear','pairacc','top1','recall3','regret','cmae'])
def P(*a):
    s=' '.join(str(x) for x in a); print(s); out.write(s+'\n')
for M in ('peo',):
    D=json.load(open(f'{UC5}/{M}/eval.json')); mins={(r['run'],r['frame']):r for r in D['mins'] if r['converged']}
    ev={}
    for (k,j,s,m,e,rc,c) in D['ev']: ev[(k,j,s,m)]=(e,rc,c)
    geo={}
    def G(r):
        key=(r['run'],r['frame'])
        if key not in geo: geo[key]=readxyz(r['dir']+'/start.opt.xyz')[1]
        return geo[key]
    allm={m:[] for m in MODS}; per={m:[] for m in MODS}
    bad=0; chk=[]
    for s,(tr,te) in SPL.items():
        cand=sorted([r for r in mins.values() if r['run'] in te],key=lambda r:r['E_ff'])
        kept=[]
        for r in cand:
            if any(abs(r['E_ff']-q['E_ff'])*HA<0.2 and rmsd(G(r),G(q))<0.3 for q in kept): continue
            kept.append(r)
        # drop failed evaluations
        ok=[r for r in kept if all(ev.get((r['run'],r['frame'],s,m),(np.nan,1,0))[1]==0 and np.isfinite(ev[(r['run'],r['frame'],s,m)][0]) for m in ('A','C','S','B2'))]
        bad+=len(kept)-len(ok); kept=ok
        ref=np.array([r['E_g2'] for r in kept])
        E={'FF':np.array([r['E_ff'] for r in kept])}
        for m in ('A','C','S','B2'): E[m]=np.array([ev[(r['run'],r['frame'],s,m)][0] for r in kept])
        for r in kept:
            for m in ('A','S'):
                e,_,c=ev[(r['run'],r['frame'],s,m)]; chk.append(abs((e-c)-r['E_ff']))
        n=len(kept)
        P(f'{M} split{s} test runs {te}: converged cand {len(cand)}, after dedupe {len(kept)}, eval failures dropped; GFN2 span {(ref.max()-ref.min())*HA:.2f} kcal, n={n}')
        for m in MODS:
            x=metrics(E[m],ref,nlow=3); per[m].append(x)
            cw.writerow([M,s,m,n]+[x[k] for k in ('spear','pairacc','top1','r5','regret','cmae')])
        # bootstrap
        for m in ('A','C','S','B2'):
            bs=[];bp=[]
            for _ in range(500):
                ii=rng.randint(0,n,n)
                if len(set(ii))<5: continue
                a=metrics(E[m][ii],ref[ii],3); b=metrics(E['FF'][ii],ref[ii],3)
                bs.append(a['spear']-b['spear']); bp.append(a['pairacc']-b['pairacc'])
            lo,hi=np.nanpercentile(bs,[2.5,97.5]); lo2,hi2=np.nanpercentile(bp,[2.5,97.5])
            P(f'   boot {m}-FF: dSpearman {per[m][-1]["spear"]-per["FF"][-1]["spear"]:+.3f} [{lo:+.3f},{hi:+.3f}]  dPairAcc {per[m][-1]["pairacc"]-per["FF"][-1]["pairacc"]:+.3f} [{lo2:+.3f},{hi2:+.3f}]')
    P(f'{M}: max |(E_model - corr) - E_FF_fresh| over minima for A,S = {max(chk):.2e} Eh (n={len(chk)}); dropped {bad}')
    P(f'{M}  model | Spearman | pairacc | top1 | recall@3 | regret | cMAE (mean over 3 splits; per split in csv)')
    for m in MODS:
        a=lambda k:np.mean([x[k] for x in per[m]])
        P(f'   {m}: {a("spear"):.3f} | {a("pairacc"):.3f} | {a("top1")*3:.0f}/3 | {a("r5"):.2f} | {a("regret"):.2f} | {a("cmae"):.2f}')
    for m in ('A','C','S','B2'):
        for k in ('spear','pairacc','regret','cmae','r5'):
            d=np.array([per[m][i][k]-per['FF'][i][k] for i in range(3)])
            better=int((d>0).sum()) if k in('spear','pairacc','r5') else int((d<0).sum())
            P(f'   paired {m}-FF {k}: {d.mean():+.3f} +- {d.std(ddof=1):.3f} per split {np.round(d,3).tolist()} better {better}/3')
out.close()
