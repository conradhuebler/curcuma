import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase5/scripts")
from fit import *
import subprocess
CUT=1
os.environ['REFM']='r2scan-3c'
NC=['-activelearn.n_centers','4000']
# cross-run: train on one MD run (all frames), evaluate on the other
for tr,te in (('CAF','CAF2'),('CAF2','CAF')):
    for tag,par in (('full_std',[]),('full_nc4000',NC),('full_se1e-5',['-activelearn.sigma_energy','1e-5'])):
        r=fit(tr,tag,1,par,test_fraction='0'); print(tr,tag,r['rc'],'%.0f'%r['wall'],flush=True)
        for mod,mf in (('J','train.deltaml.json'),('E','train.deltaml.energy.json')):
            out=f'{UC5}/xeval_{tr}_{tag}_{mod}.json'
            if os.path.exists(out): continue
            subprocess.run(['python3','-I',UC5+'/scripts/eval_on.py',f'{UC5}/fits/{tr}/{tag}_s1/{mf}',ds(te),out,f'{UC5}/xeval_wd'],capture_output=True)
            print(tr,'->',te,tag,mod,open(out).read(),flush=True)
# combined 501 frames, held-out seeds 1-5
for tag,par in (('std',[]),('nc4000',NC),('se1e-5',['-activelearn.sigma_energy','1e-5'])):
    for s in range(1,6):
        r=fit('CAFALL',tag,s,par); print('CAFALL',tag,s,r['rc'],'%.0f'%r['wall'],flush=True)
