import sys,pickle
sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase7/scripts')
from common7 import *
OUT='/home/conrad/src/curcuma/build_deltaml/usecase8/analysis'
D=pickle.load(open(OUT+'/D8.pkl','rb'))
def pool(tag,ss): return {q:np.concatenate([D[(tag,s)][q] for s in ss]) for q in ('Rg','ee','E','ang')}
REF=pool('gfn2',[2,3,4]);Eref=REF['E'].mean()
def met(a,b,e):
    return np.array([np.mean([jsd(hist(a['ang'][:,i]),hist(b['ang'][:,i])) for i in range(21)]),w1(a['E']-e,b['E']-e)*HA,w1(a['Rg'],b['Rg']),w1(a['ee'],b['ee'])])
n3=[]
for s in (1,2,3,4):
    rest=pool('gfn2',[x for x in (1,2,3,4) if x!=s]);n3.append(met(pool('gfn2',[s]),rest,rest['E'].mean()))
n3=np.array(n3);print('noise3 (GFN2 seed vs pool of the other 3 of 1-4; same sizes as models) JSD per seed',np.round(n3[:,0],3),'mean %.3f+-%.3f'%(n3[:,0].mean(),n3[:,0].std(ddof=1)),'W1E %.2f'%n3[:,1].mean(),'Rg %.2f ee %.2f'%(n3[:,2].mean(),n3[:,3].mean()))
tags=['gfnff','V0_A','V0_B2','V1_A','V1_B2','V2_A','V2_B2','V3_v3half','V3_v3a','V3_v3b','V3_v3c','V3_v3d']
M={t:np.array([met(pool(t,[s]),REF,Eref) for s in (1,2,3,4)]) for t in tags}
print('per-seed JSD:')
for t in tags: print(f'{t:10s}',np.round(M[t][:,0],3),f'mean {M[t][:,0].mean():.3f}+-{M[t][:,0].std(ddof=1):.3f}  seeds JSD<=max noise3 ({n3[:,0].max():.3f}): {(M[t][:,0]<=n3[:,0].max()).sum()}/4')
def pr(a,b):
    o=[]
    for k,n in enumerate(('JSD','W1E','Rg','ee')):
        d=M[a][:,k]-M[b][:,k];o.append(f'{n} {d.mean():+.3f}+-{d.std(ddof=1)/2:.3f} ({(d<0).sum()}/4)')
    return f'{a} - {b}: '+'; '.join(o)
for a,b in [('V1_A','V0_A'),('V1_B2','V0_B2'),('V2_A','V0_A'),('V2_B2','V0_B2'),('V3_v3a','V0_A'),('V3_v3a','V0_B2'),('V3_v3a','gfnff'),('V3_v3b','V0_B2'),('V3_v3b','gfnff'),('V3_v3b','V1_B2'),('V3_v3b','V2_B2'),('V3_v3half','gfnff'),('V3_v3a','V3_v3half'),('V3_v3c','V0_B2'),('V3_v3d','V0_B2'),('V3_v3c','gfnff'),('V3_v3d','gfnff'),('V1_B2','gfnff'),('V2_B2','gfnff')]: print(pr(a,b))
# per-torsion of V2_B2 and V1_B2 worst
tors=None
import json
for t in ('V1_B2','V2_B2','V3_v3b','V3_v3c','V3_v3d','V0_B2','gfnff'):
    P=pool(t,[1,2,3,4]);pj=np.array([jsd(hist(P['ang'][:,i]),hist(REF['ang'][:,i])) for i in range(21)])
    print(t,'worst torsion idx',np.argsort(pj)[-3:][::-1],np.round(np.sort(pj)[-3:][::-1],3),'n torsions JSD>0.2:',int((pj>0.2).sum()))
