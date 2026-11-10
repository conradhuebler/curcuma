# usage: an10.py tag ...  torsion JSD of final models vs pooled GFN2 seeds 2-4, as usecase9 final9.py (Claude Generated Oct 2026)
import sys,pickle,os,json
sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase7/scripts')
import common7
from common7 import *
U10='/home/conrad/src/curcuma/build_deltaml/usecase10';OUT=U10+'/analysis';os.makedirs(OUT,exist_ok=True)
common7.BIN=U10+'/curcuma_uc10'
D9=pickle.load(open('/home/conrad/src/curcuma/build_deltaml/usecase9/analysis/D9.pkl','rb'))
st=read_trj(U+'/peo8_start.xyz')[0];tors=backbone_torsions(st[2],st[3])
REF={q:np.concatenate([D9[('gfn2',s)][q] for s in (2,3,4)]) for q in ('ang','E')};Eref=REF['E'].mean()
L=[]
for tag in sys.argv[1:]:
    J=[];hm=[];hn=[]
    for s in (1,2,3,4):
        wd=f'{U10}/runs/{tag}_s{s}'
        seen=set();F=[]
        for f in read_trj(wd+'/md.snapshots/md.trj.xyz'):
            if f[0]>=5000 and f[0] not in seen: seen.add(f[0]);F.append(f)
        ang=np.array([[dihedral(f[3],t) for t in tors] for f in F])
        J.append(np.mean([jsd(hist(ang[:,i]),hist(REF['ang'][:,i])) for i in range(ang.shape[1])]))
        cf=f'{OUT}/E_{tag}_s{s}.npy';idx=np.unique(np.round(np.linspace(0,len(F)-1,150)).astype(int))
        if os.path.exists(cf): E=np.load(cf)
        else: E=gfn2_sp([(F[k][2],F[k][3]) for k in idx],OUT,f'sp_{tag}_s{s}');np.save(cf,E)
        d=(np.array([F[k][1] for k in idx])-E)*HA;hm.append(d.mean());hn.append(d.min())
    J=np.array(J);L.append(f"{tag}: JSD per MD seed {np.round(J,3)} mean {J.mean():.3f}+-{J.std(ddof=1):.3f}; model-GFN2 own frames mean {np.mean(hm):+.1f} deepest {min(hn):.1f}")
open(OUT+'/an10.txt','a').write('\n'.join(L)+'\n');print('\n'.join(L))
