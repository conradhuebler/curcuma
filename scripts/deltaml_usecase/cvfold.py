# usage: cvfold.py rep seed fold   -> trains A,B on other folds, predicts test fold with A,B,C
import sys,shutil,time,os,json,subprocess; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase/scripts')
from common import *
seed=int(sys.argv[1]); fold=int(sys.argv[2])
d=load();S=d['structures'];n=len(S)
perm=np.random.RandomState(seed).permutation(n); folds=np.array_split(perm,5)
test=sorted(folds[fold].tolist()); train=sorted(np.concatenate([folds[k] for k in range(5) if k!=fold]).tolist())
base=UC+'/cv/s%d_f%d'%(seed,fold)
B_ARGS=['-activelearn.n_centers','4000','-activelearn.lambda','1e-7','-activelearn.acsf_rcut_radial','7','-activelearn.acsf_zeta','1,2,4,8,16']
fit={}
for var,extra in [('A',[]),('B',B_ARGS)]:
    dd=base+'/'+var; write_train(d,train,dd)
    t=time.time()
    with open(dd+'/log.txt','w') as lf:
        rc=subprocess.run([BIN,'-activelearn','train.xyz','-method','gfnff','-ref_method','gfn2','-activelearn.dataset',dd+'/train.dataset.json','-test_fraction','0','-activelearn.seed','1','-threads','8','-no_bmt']+extra,cwd=dd,stdout=lf,stderr=subprocess.STDOUT).returncode
    fit[var]=(rc,time.time()-t)
models={'A':base+'/A/train.deltaml.json','B':base+'/B/train.deltaml.json','C':base+'/A/train.deltaml.energy.json'}
rows=[]
for i in test:
    wd=base+'/sp'; shutil.rmtree(wd,ignore_errors=True); os.makedirs(wd)
    open(wd+'/x.xyz','w').write(xyz_text(S,[i]))
    r={'idx':i,'fold':fold,'seed':seed,'E_ref':S[i]['energy_ref'],'E_base':S[i]['energy_base']}
    for m,p in models.items():
        e,rc=sp_energy('x.xyz','deltaml:'+p,wd,1); r['E_'+m]=e; r['rc_'+m]=rc
    rows.append(r)
shutil.rmtree(base+'/sp',ignore_errors=True)
json.dump({'rows':rows,'fit':fit,'train_n':len(train)},open(base+'/result.json','w'))
