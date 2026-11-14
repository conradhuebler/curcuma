import sys,pickle,itertools
sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase7/scripts')
from common7 import *
MOL=sys.argv[1];PART=sys.argv[2]
D,tors,Eref=pickle.load(open(f'{U}/analysis_{PART}_{MOL}/D.pkl','rb'))
tags=sorted({k[0] for k in D});seeds=sorted({k[1] for k in D})
def pool(tag,ss): return {q:np.concatenate([D[(tag,s)][q] for s in ss]) for q in ('Rg','ee','E','ang')}
def met(a,b):
    j=np.mean([jsd(hist(a['ang'][:,i]),hist(b['ang'][:,i])) for i in range(a['ang'].shape[1])])
    return np.array([j,w1(a['E']-Eref,b['E']-Eref)*HA,w1(a['Rg'],b['Rg']),w1(a['ee'],b['ee'])])
h=len(seeds)//2
L=[f'{MOL}: half-vs-complement-half (seeds {seeds}); method pooled over {h} seeds vs GFN2 pooled over the complementary {len(seeds)-h} seeds; all {len(list(itertools.combinations(seeds,h)))} splits (not independent)','tag    JSD    W1_E[kcal]  W1_Rg[A]  W1_ee[A]   (mean+-std over splits)']
res={}
for t in tags:
    r=[]
    for c in itertools.combinations(seeds,h):
        o=[s for s in seeds if s not in c]
        r.append(met(pool(t,c),pool('gfn2',o)))
    r=np.array(r);res[t]=r
    L.append(f"{t:6s} "+'  '.join(f"{r[:,k].mean():.4f}+-{r[:,k].std(ddof=1):.4f}" for k in range(4)))
for t in tags:
    if t in ('gfn2','gfnff'): continue
    d=res[t]-res['gfnff'];L.append(f"{t} minus gfnff over splits: mean "+' '.join(f"{d[:,k].mean():+.4f}" for k in range(4))+"; splits better "+' '.join(f"{int((d[:,k]<0).sum())}/{len(d)}" for k in range(4)))
# per torsion JSD (pooled all seeds)
REF=pool('gfn2',seeds)
for t in tags:
    if t=='gfn2': continue
    P=pool(t,seeds);pj=[jsd(hist(P['ang'][:,i]),hist(REF['ang'][:,i])) for i in range(P['ang'].shape[1])]
    L.append(f"per-torsion JSD pooled {t}: median {np.median(pj):.3f} max {max(pj):.3f} (idx {int(np.argmax(pj))}: {tors[int(np.argmax(pj))]})")
# trans fraction of OCCO (gauche) 
for t in tags:
    P=pool(t,seeds);a=P['ang']
    g=[i for i,q in enumerate(tors) if False]
open(f'{U}/analysis_{PART}_{MOL}/pooled.txt','w').write('\n'.join(L)+'\n');print('\n'.join(L))
