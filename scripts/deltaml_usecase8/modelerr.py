# model energy (MD comment) minus GFN2 SP on the SP frames, per tag, seed
import sys,os
sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase7/scripts')
from common7 import *
U8='/home/conrad/src/curcuma/build_deltaml/usecase8'
for tag in sys.argv[1:]:
    out=[]
    for s in (1,2,3,4):
        wd=f'{U}/runs/A/peo8_{tag[3:]}_s{s}' if tag.startswith('V0_') else f'{U8}/runs/{tag}_s{s}'
        cf=f'{U}/analysis_A_peo8/E_{tag[3:]}_s{s}.npy' if tag.startswith('V0_') else f'{U8}/analysis/E_{tag}_s{s}.npy'
        if not os.path.exists(cf): continue
        seen=set();F=[]
        for f in read_trj(wd+'/md.snapshots/md.trj.xyz'):
            if f[0]>=5000 and f[0] not in seen: seen.add(f[0]);F.append(f)
        idx=np.unique(np.round(np.linspace(0,len(F)-1,150)).astype(int))
        em=np.array([F[k][1] for k in idx]);eg=np.load(cf)
        d=(em-eg)*HA;out.append((d.mean(),d.std(),np.abs(d).max()))
    print(tag,' | '.join(f'm{a:+.2f} sd{b:.2f} max{c:.1f}' for a,b,c in out))
