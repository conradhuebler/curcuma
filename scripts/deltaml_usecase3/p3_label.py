# dedup converged minima, then label (fresh gfnff + gfn2 with gradients) -> M/p3/min_dataset.json
import sys,json,subprocess; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase3/scripts')
from common3 import *
from concurrent.futures import ThreadPoolExecutor
M=sys.argv[1]; B=f'{UC}/{M}/p3'
res=[json.load(open(f'{B}/opt/{d}/res.json')) for d in sorted(os.listdir(B+'/opt'))]
conv=[r for r in res if r['converged']]
print(M,'opt runs',len(res),'converged',len(conv),'nonconv idx',[r['idx'] for r in res if not r['converged']])
def lab(r):
    i=r['idx']; wd=f'{B}/lab/{i:03d}'; os.makedirs(wd+'/base',exist_ok=True); os.makedirs(wd+'/ref',exist_ok=True)
    for s in ('base','ref'): shutil.copy(f'{B}/opt/{i:03d}/start.opt.xyz',f'{wd}/{s}/x.xyz')
    rb=subprocess.run([BIN,'-sp','x.xyz','-method','gfnff','-gfnff.cache_topology','false','-dump_gradient','g.txt','-threads','1','-no_bmt'],cwd=wd+'/base',capture_output=True)
    rr=subprocess.run([BIN,'-sp','x.xyz','-method','gfn2','-dump_gradient','g.txt','-threads','1','-no_bmt'],cwd=wd+'/ref',capture_output=True)
    el,X=readxyz(f'{wd}/base/x.xyz'); eb,gb=dumpg(wd+'/base/g.txt'); er,gr=dumpg(wd+'/ref/g.txt')
    return dict(idx=i,rc=(rb.returncode,rr.returncode),el=el,X=X,eb=eb,er=er,gb=gb,gr=gr)
import shutil
with ThreadPoolExecutor(24) as ex: L=list(ex.map(lab,conv))
assert all(l['rc']==(0,0) for l in L)
L.sort(key=lambda l:l['eb'])
keep=[]; dup=0
for l in L:
    if any(rmsd(l['X'],k['X'])<0.3 and abs(l['eb']-k['eb'])*HA<0.2 for k in keep): dup+=1
    else: keep.append(l)
print(M,'unique minima',len(keep),'duplicates removed',dup)
S=[dict(charge=0,cycle=-1,elements=[Z[e] for e in l['el']],energy_base=l['eb'],energy_ref=l['er'],geometry=l['X'].flatten().tolist(),
  gradient_base=l['gb'].flatten().tolist(),gradient_ref=l['gr'].flatten().tolist(),has_gradients=True,spin=0.0,start_idx=l['idx']) for l in keep]
json.dump({'format':'curcuma-deltaml-dataset','n_structures':len(S),'structures':S,'version':1,'n_converged':len(conv),'n_dup':dup},open(B+'/min_dataset.json','w'))
er=np.array([s['energy_ref'] for s in S])*HA; eb=np.array([s['energy_base'] for s in S])*HA
print(M,'GFN2 rel span %.2f std %.2f | GFN-FF rel span %.2f std %.2f | Spearman %.3f'%(np.ptp(er),er.std(ddof=1),np.ptp(eb),eb.std(ddof=1),sp(er,eb)))
