import sys; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase6/scripts')
from common6 import *
from concurrent.futures import ThreadPoolExecutor
jobs=[]
for M in ('M1','M2','PEO'):
    N=len(load(M)['structures']); te=set(testidx(M)); tr=[i for i in range(N) if i not in te]
    for tag,par in (('A',[]),('S',['-activelearn.lambda','1e-3']),('B2',['-activelearn.n_centers','4000'])): jobs.append((M,'p1'+tag,tr,par))
for M in ('M1','PEO'):
    N=len(load(M)['structures']); tr=[i for i in range(N) if i!=200]
    for tag,par in (('A',[]),('B2',['-activelearn.n_centers','4000'])): jobs.append((M,'p2'+tag,tr,par))
def run(j):
    r=fit(*j); print(j[0],j[1],r['rc'],'%.0f'%r['wall'],r['model'],flush=True)
with ThreadPoolExecutor(3) as ex: list(ex.map(run,jobs))
print('DONE')
