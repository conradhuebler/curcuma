# usage: train3.py DS modeldir variant(A|S) idxfile   (DS = dataset name dir under usecase3)
import sys,subprocess,time,json; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase3/scripts')
from common3 import *
def train(DS,dd,var,idx,dsfile=None):
    d=json.load(open(dsfile)) if dsfile else load(DS)
    write_train(d,idx,dd)
    extra=['-activelearn.lambda','1e-3'] if var=='S' else []
    t=time.time()
    with open(dd+'/log.txt','w') as lf:
        rc=subprocess.run([BIN,'-activelearn','train.xyz','-method','gfnff','-ref_method','gfn2','-activelearn.dataset',dd+'/train.dataset.json','-test_fraction','0','-activelearn.seed','1','-threads','8','-no_bmt']+extra,cwd=dd,stdout=lf,stderr=subprocess.STDOUT).returncode
    return rc,time.time()-t
