import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase4/scripts")
from fit import *
from concurrent.futures import ThreadPoolExecutor
jobs=[(d,t,s,p) for d in ('E1','M1','M2') for s in (11,12,13,14,15) for t,p in (('std',[]),('nc4000',['-activelearn.n_centers','4000']))]
def run(j):
    r=fit(*j); print(j[0],j[1],j[2],r.get('rc'),'%.0f'%r.get('wall',0),flush=True)
with ThreadPoolExecutor(3) as ex: list(ex.map(run,jobs))
print('DONE')
