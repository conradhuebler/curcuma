import sys; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase/scripts')
from common import *
exec(open(UC+'/scripts/t2_an.py').read().split("Ms=[")[0])
sel_set=set(sel); tr=[np.array(S[j]['geometry']).reshape(-1,3) for j in range(len(S)) if j not in sel_set]
out=["nearest-training-structure RMSD (A), mean+-std over 20 [min over 333 training conformers]:"]
for m in ['start','gfn2','gfnff','A','B','C']:
    v=[]
    for i in sel:
        X=np.array(S[i]['geometry']).reshape(-1,3) if m=='start' else rd('%s/t2/s%03d/%s/start.opt.xyz'%(UC,i,m))[1]
        v.append(min(rmsd(X,T) for T in tr))
    out.append("  %-6s %.3f+-%.3f (n=20)"%(m,np.mean(v),np.std(v,ddof=1)))
open(UC+'/test2_domain.txt','w').write('\n'.join(out)+'\n');print('\n'.join(out))
