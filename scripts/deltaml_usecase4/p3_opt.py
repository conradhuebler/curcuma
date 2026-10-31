# p3_opt.py M k j : GFN-FF opt of frame j of run k, then fresh GFN-FF + GFN2 sp at result
import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase4/scripts")
from common4 import *
import json,shutil,subprocess,time
def job(M,k,j):
    S=json.load(open(f'{UC4}/p3/{M}/md{k}/frames.deltaml.dataset.json'))['structures']
    wd=f'{UC4}/p3/{M}/opt/r{k}_{j:03d}'
    if os.path.exists(wd+'/res.json'): return
    shutil.rmtree(wd,ignore_errors=True); os.makedirs(wd)
    open(wd+'/start.xyz','w').write(xyz_text(S,[j]))
    with open(wd+'/log.txt','w') as lf:
        rc=subprocess.run([BIN4,'-opt','start.xyz','-method','gfnff','-threads','8','-no_bmt'],cwd=wd,stdout=lf,stderr=subprocess.STDOUT).returncode
    txt=ANSI.sub('',open(wd+'/log.txt').read())
    rows=re.findall(r'^\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)',txt,re.M)
    res={'run':k,'frame':j,'rc':rc,'steps':int(rows[-1][0]) if rows else -1,'grad_last':float(rows[-1][4]) if rows else float('nan'),'has_opt':os.path.exists(wd+'/start.opt.xyz')}
    res['converged']=bool(rc==0 and res['has_opt'] and res['grad_last']<=5e-4)
    if res['has_opt']:
        sd=wd+'/sp'; os.makedirs(sd); shutil.copy(wd+'/start.opt.xyz',sd+'/x.xyz')
        res['E_ff'],r1,_=sp_energy('x.xyz','gfnff',sd,8,['-gfnff.cache_topology','false'])
        res['E_g2'],r2,_=sp_energy('x.xyz','gfn2',sd,8)
        res['rc_sp']=[r1,r2]
    json.dump(res,open(wd+'/res.json','w'))
if __name__=='__main__':
    from concurrent.futures import ThreadPoolExecutor
    M=sys.argv[1]; J=np.linspace(0,350,40).round().astype(int).tolist()
    jobs=[(M,k,j) for k in range(1,7) for j in J]
    with ThreadPoolExecutor(3) as ex: list(ex.map(lambda a:job(*a),jobs))
    print(M,'opt done',flush=True)
