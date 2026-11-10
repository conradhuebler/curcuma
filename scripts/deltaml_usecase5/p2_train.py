import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase5/scripts")
from fit import *
from concurrent.futures import ThreadPoolExecutor
jobs=[]
for M in ('PEO',):
    te=sorted(np.random.RandomState(7).choice(351,20,replace=False).tolist()); tr=[i for i in range(351) if i not in te]
    json.dump({'test':te,'train':tr},open(f'{UC5}/p2_split_{M}.json','w'))
    for tag,par in (('p2A',[]),('p2S',['-activelearn.lambda','1e-3']),('p2B2',['-activelearn.n_centers','4000'])): jobs.append((M,tag,1,par,tr))
def run(j):
    r=fit(j[0],j[1],j[2],j[3],j[4],test_fraction='0'); print(j[0],j[1],r['rc'],'%.0f'%r['wall'],flush=True)
with ThreadPoolExecutor(int(os.environ.get('WORKERS','2'))) as ex: list(ex.map(run,jobs))
