import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase4/scripts")
from p2_opt import *
from concurrent.futures import ThreadPoolExecutor
jobs=[]
for M in ('M1','M2'):
    te=json.load(open(UC4+f'/p2_split_{M}.json'))['test']
    for i in te:
        jobs.append((M,i,'gfn2','gfn2',()))
        jobs.append((M,i,'gfnff','gfnff',()))
        for tag in ('A','S','B2'):
            mp='deltaml:'+UC4+f'/fits/{M}/p2{tag}_s1/train.deltaml.json'
            jobs.append((M,i,tag,mp,()))
            for th in ('0.5','0.9'): jobs.append((M,i,f'{tag}_ab{th}',mp,('-deltaml.uncertainty_abort',th)))
print(len(jobs),flush=True)
def run(j):
    r=job(*j); return r
with ThreadPoolExecutor(3) as ex: list(ex.map(run,jobs))
print('DONE',flush=True)
