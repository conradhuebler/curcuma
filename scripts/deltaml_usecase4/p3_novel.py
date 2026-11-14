import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase4/scripts")
from common4 import *
import json
SPL={1:([1,2,3,4],[5,6]),2:([3,4,5,6],[1,2]),3:([1,2,5,6],[3,4])}
o=open(UC4+'/p3_novelty.txt','w')
for M in ('M1','M2'):
    D=json.load(open(f'{UC4}/p3/{M}/eval.json')); mins=[r for r in D['mins'] if r['converged']]
    ev={(k,j,s,m):e for (k,j,s,m,e,rc,c) in D['ev']}
    for r in mins: r['g']=readxyz(r['dir']+'/start.opt.xyz')[1]
    for s,(tr,te) in SPL.items():
        T=[r for r in mins if r['run'] in tr]; X=[r for r in mins if r['run'] in te]
        nov=[]; mr=[]
        for r in X:
            d=min(rmsd(r['g'],q['g']) for q in T); mr.append(d)
            same=any(abs(r['E_ff']-q['E_ff'])*HA<0.2 and rmsd(r['g'],q['g'])<0.3 for q in T); nov.append(not same)
        nov=np.array(nov); mr=np.array(mr)
        # pair accuracy FF vs A on pairs where both minima are novel
        ref=np.array([r['E_g2'] for r in X])*HA; ff=np.array([r['E_ff'] for r in X])*HA; A=np.array([ev[(r['run'],r['frame'],s,'A')] for r in X])*HA
        idx=np.where(nov)[0]; iu=np.triu_indices(len(idx),1)
        def pa(E):
            E=E[idx]; R=ref[idx]; dE=(R[:,None]-R[None,:])[iu]; mE=(E[:,None]-E[None,:])[iu]; m=abs(dE)>1
            return np.mean(np.sign(dE[m])==np.sign(mE[m]))
        s_=f'{M} split{s}: test minima {len(X)}, novel (no GFN-FF training minimum within 0.3 A & 0.2 kcal) {nov.sum()}, median RMSD to nearest training minimum {np.median(mr):.2f} A; pair acc on novel-novel pairs FF {pa(ff):.3f} A {pa(A):.3f}'
        print(s_); o.write(s_+'\n')
