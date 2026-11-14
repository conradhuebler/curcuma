import sys; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase2/scripts')
from common2 import *
def rk(v):
    r=np.empty(len(v)); r[np.argsort(v,kind='stable')]=np.arange(len(v)); return r
def dist(name,x):
    q=np.percentile(x,[10,25,50,75,90])
    return "%s min %.2f q10 %.2f q25 %.2f q50 %.2f q75 %.2f q90 %.2f max %.2f std %.2f"%(name,x.min(),*q,x.max(),x.std(ddof=1))
out=[]
for M in ('M1','M2'):
    d=load(M);S=d['structures'];n=len(S)
    er=np.array([s['energy_ref'] for s in S])*HA; eb=np.array([s['energy_base'] for s in S])*HA
    er-=er.min(); eb-=eb.min()
    out.append(f'== {M} n={n} natoms={len(S[0]["elements"])}')
    out.append(dist('GFN2 rel [kcal/mol]',er)); out.append(dist('GFN-FF fresh rel',eb))
    pr=np.corrcoef(er,eb)[0,1]; sp=np.corrcoef(rk(er),rk(eb))[0,1]
    a,b=np.polyfit(eb,er,1); rm=np.sqrt(np.mean((a*eb+b-er)**2))
    out.append("Pearson %.3f Spearman %.3f RMSE after linear fit %.3f kcal/mol (slope %.3f)"%(pr,sp,rm,a))
    out.append("within GFN2 min: "+', '.join("%g kcal: %d"%(t,(er<=t).sum()) for t in (1,2,5,10)))
    X=[np.array(s['geometry']).reshape(-1,3) for s in S]
    D=np.zeros((n,n))
    for i in range(n):
        for j in range(i+1,n): D[i,j]=D[j,i]=rmsd(X[i],X[j])
    np.save(f'{UC}/{M}/rmsd_matrix.npy',D)
    iu=np.triu_indices(n,1); p=D[iu]
    out.append("pairwise RMSD [A]: "+dist('',p).strip()+" mean %.3f"%p.mean())
    Dn=D+np.eye(n)*1e9; nn=Dn.min(1)
    out.append("nearest-neighbour RMSD per frame (distance to nearest other frame): mean %.3f max %.3f median %.3f"%(nn.mean(),nn.max(),np.median(nn)))
    # nearest-neighbour distance when 20% held out (CV-like): mean over folds
    perm=np.random.RandomState(101).permutation(n); folds=np.array_split(perm,5); nd=[]
    for f in folds:
        tr=np.setdiff1d(np.arange(n),f); nd+=list(D[np.ix_(f,tr)].min(1))
    out.append("distance test frame -> nearest training frame (5-fold, seed 101): mean %.3f max %.3f"%(np.mean(nd),np.max(nd)))
open(UC+'/test0_summary.txt','w').write('\n'.join(out)+'\n'); print('\n'.join(out))
