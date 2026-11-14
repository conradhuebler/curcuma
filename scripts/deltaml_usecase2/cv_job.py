# usage: cv_job.py M fold  (seed 101): trains A,S on other folds, predicts test fold with A,C,S
import sys,json,shutil,subprocess; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase2/scripts')
from common2 import *
from concurrent.futures import ThreadPoolExecutor
M=sys.argv[1]; fold=int(sys.argv[2]); seed=int(sys.argv[3]) if len(sys.argv)>3 else 101
d=load(M);S=d['structures'];n=len(S)
perm=np.random.RandomState(seed).permutation(n); folds=np.array_split(perm,5)
test=sorted(folds[fold].tolist()); train=sorted(np.concatenate([folds[k] for k in range(5) if k!=fold]).tolist())
tag=('sl' if os.environ.get('DSVAR')=='sl' else '')+f'cv{seed}_f{fold}'; base=f'{UC}/{M}/models/{tag}'
os.makedirs(base,exist_ok=True); json.dump(train,open(base+'/train_idx.json','w'))
for v in 'AS':
    subprocess.run(['python3','-I',UC+'/scripts/train.py',M,tag,v,base+'/train_idx.json'],check=True)
models={'A':f'{base}/A/train.deltaml.json','C':f'{base}/A/train.deltaml.energy.json','S':f'{base}/S/train.deltaml.json'}
def one(i):
    wd=f'{base}/sp/{i}'; os.makedirs(wd,exist_ok=True); open(wd+'/x.xyz','w').write(xyz_text(S,[i]))
    r={'idx':i,'E_ref':S[i]['energy_ref'],'E_base_stored':S[i]['energy_base']}
    for m,p in models.items():
        e,rc,t=sp_energy('x.xyz','deltaml:'+p,wd,1,['-verbosity','3'])
        mm=re.search(r'deltaml_correction:\s*(-?\d+\.\d+(?:[eE][-+]?\d+)?)',t)
        r['E_'+m]=e; r['rc_'+m]=rc; r['corr_'+m]=float(mm.group(1)) if mm else float('nan')
        r['nan_'+m]=bool(re.search(r'\bnan\b',t.lower()))
    return r
with ThreadPoolExecutor(8) as ex: rows=list(ex.map(one,test))
json.dump({'rows':rows,'train_n':len(train)},open(base+'/result.json','w'))
shutil.rmtree(base+'/sp',ignore_errors=True)
print(M,fold,'done',len(rows))
