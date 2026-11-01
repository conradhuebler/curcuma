# p3_eval.py M : energies of all models at all converged minima of the 6 runs, ranking per split
import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase4/scripts")
from common4 import *
import json,glob,shutil
from concurrent.futures import ThreadPoolExecutor
M=sys.argv[1]
SPL={1:([1,2,3,4],[5,6]),2:([3,4,5,6],[1,2]),3:([1,2,5,6],[3,4])}
mins=[]
for p in sorted(glob.glob(f'{UC4}/p3/{M}/opt/r*/res.json')):
    r=json.load(open(p)); r['dir']=os.path.dirname(p)
    mins.append(r)
print(M,'opt results',len(mins),'converged',sum(r['converged'] for r in mins),flush=True)
def mpath(s,tag,en=False):
    return f'{UC4}/fits/{M}_split{s}_train/{tag}_s1/train.deltaml{".energy" if en else ""}.json'
models={'A':lambda s:mpath(s,'A'),'C':lambda s:mpath(s,'A',True),'S':lambda s:mpath(s,'S'),'B2':lambda s:mpath(s,'B2')}
out={}
def ev(a):
    r,s,m=a; sd=r['dir']+f'/ev_{m}_{s}'; os.makedirs(sd,exist_ok=True); shutil.copy(r['dir']+'/start.opt.xyz',sd+'/x.xyz')
    e,rc,t=sp_energy('x.xyz','deltaml:'+models[m](s),sd,1,['-verbosity','3'])
    mm=re.search(r'deltaml_correction:\s*(-?\d+\.\d+(?:[eE][-+]?\d+)?)',t)
    return (r['run'],r['frame'],s,m,e,rc,float(mm.group(1)) if mm else float('nan'))
jobs=[(r,s,m) for r in mins if r['converged'] for s,(tr,te) in SPL.items() if r['run'] in te for m in models]
print(len(jobs),'model evaluations',flush=True)
with ThreadPoolExecutor(3) as ex: res=list(ex.map(ev,jobs))
json.dump({'mins':mins,'ev':res},open(f'{UC4}/p3/{M}/eval.json','w'))
print('DONE')
