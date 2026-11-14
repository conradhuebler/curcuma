# usage: cv3.py DS mode(rand|block) seed fold   -> DS/models/{mode}{seed}_f{fold}/result.json
import sys,json,shutil,subprocess; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase3/scripts')
from common3 import *
from train3 import train
from concurrent.futures import ThreadPoolExecutor
DS,mode,seed,fold=sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4])
d=load(DS);S=d['structures'];n=len(S)
if mode=='rand':
    perm=np.random.RandomState(seed).permutation(n); folds=np.array_split(perm,5)
    test=sorted(folds[fold].tolist()); train_i=sorted(np.concatenate([folds[k] for k in range(5) if k!=fold]).tolist())
else:
    blocks=np.array_split(np.arange(n),5); test=blocks[fold].tolist(); a,b=test[0],test[-1]
    excl=set(range(a-5,b+6)); train_i=[i for i in range(n) if i not in excl]
tag=f'{mode}{seed}_f{fold}'; base=f'{UC}/{DS}/models/{tag}'
os.makedirs(base,exist_ok=True); json.dump(train_i,open(base+'/train_idx.json','w'))
for v in 'AS':
    rc,t=train(DS,f'{base}/{v}',v,train_i)
    print(DS,tag,v,'rc',rc,'%.0fs'%t,len(train_i),flush=True)
models={'A':f'{base}/A/train.deltaml.json','C':f'{base}/A/train.deltaml.energy.json','S':f'{base}/S/train.deltaml.json'}
def one(i):
    wd=f'{base}/sp/{i}'; os.makedirs(wd,exist_ok=True); open(wd+'/x.xyz','w').write(xyz_text(S,[i]))
    r={'idx':i,'E_ref':S[i]['energy_ref'],'E_base_stored':S[i]['energy_base']}
    for m,p in models.items():
        e,rc,t=sp_energy('x.xyz','deltaml:'+p,wd,1,['-verbosity','3'])
        mm=re.search(r'deltaml_correction:\s*(-?\d+\.\d+(?:[eE][-+]?\d+)?)',t)
        r['E_'+m]=e; r['rc_'+m]=rc; r['corr_'+m]=float(mm.group(1)) if mm else float('nan')
    return r
with ThreadPoolExecutor(8) as ex: rows=list(ex.map(one,test))
# fresh GFN-FF check on first 3 test structures
fresh=[]
for i in test[:3]:
    wd=f'{base}/fr/{i}'; os.makedirs(wd,exist_ok=True); open(wd+'/x.xyz','w').write(xyz_text(S,[i]))
    e,rc,t=sp_energy('x.xyz','gfnff',wd,1,['-gfnff.cache_topology','false']); fresh.append((i,e,rc))
json.dump({'rows':rows,'train_n':len(train_i),'fresh':fresh},open(base+'/result.json','w'))
shutil.rmtree(base+'/sp',ignore_errors=True); shutil.rmtree(base+'/fr',ignore_errors=True)
print(DS,tag,'done',len(rows),flush=True)
