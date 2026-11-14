# usage: train.py M tag variant(A|S) idxfile|all
import sys,subprocess,time,json; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase2/scripts')
from common2 import *
M,tag,var,sel=sys.argv[1:5]
d=load(M); n=len(d['structures'])
idx=list(range(n)) if sel=='all' else json.load(open(sel))
dd=f'{UC}/{M}/models/{tag}/{var}'; write_train(d,idx,dd)
extra=['-activelearn.lambda','1e-3'] if var=='S' else []
t=time.time()
with open(dd+'/log.txt','w') as lf:
    rc=subprocess.run([BIN,'-activelearn','train.xyz','-method','gfnff','-ref_method','gfn2','-activelearn.dataset',dd+'/train.dataset.json','-test_fraction','0','-activelearn.seed','1','-threads','8','-no_bmt']+extra,cwd=dd,stdout=lf,stderr=subprocess.STDOUT).returncode
print(M,tag,var,'rc',rc,'%.0fs'%(time.time()-t),len(idx))
