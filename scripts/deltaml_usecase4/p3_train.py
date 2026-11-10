# p3_train.py M : train A (+C), S, B2 for the 3 splits
import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase4/scripts")
from fit import *
from concurrent.futures import ThreadPoolExecutor
M=sys.argv[1]
SPL={1:([1,2,3,4],[5,6]),2:([3,4,5,6],[1,2]),3:([1,2,5,6],[3,4])}
jobs=[]
for s,(tr,te) in SPL.items():
    S=[]
    for k in tr: S+=json.load(open(f'{UC4}/p3/{M}/md{k}/frames.deltaml.dataset.json'))['structures']
    d=json.load(open(f'{UC4}/p3/{M}/md{tr[0]}/frames.deltaml.dataset.json')); d['structures']=S; d['n_structures']=len(S)
    p=f'{UC4}/p3/{M}/{M}_split{s}_train.json'; json.dump(d,open(p,'w'))
    for tag,par in (('A',[]),('S',['-activelearn.lambda','1e-3']),('B2',['-activelearn.n_centers','4000'])): jobs.append((p,tag,1,par))
def run(j):
    r=fit(j[0],j[1],j[2],j[3],test_fraction='0'); print(j[0].split('/')[-1],j[1],r['rc'],'%.0f'%r['wall'],flush=True)
with ThreadPoolExecutor(3) as ex: list(ex.map(run,jobs))
print('DONE')
