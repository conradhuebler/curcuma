# usage: p3_opt.py M idx : GFN-FF opt (default optimizer) from MD frame idx
import sys,time,shutil,json,subprocess; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase3/scripts')
from common3 import *
M,i=sys.argv[1],int(sys.argv[2])
S=load(M)['structures']
wd=f'{UC}/{M}/p3/opt/{i:03d}'; shutil.rmtree(wd,ignore_errors=True); os.makedirs(wd)
open(wd+'/start.xyz','w').write(xyz_text(S,[i]))
t=time.time()
with open(wd+'/log.txt','w') as lf:
    rc=subprocess.run([BIN,'-opt','start.xyz','-method','gfnff','-threads','8','-no_bmt'],cwd=wd,stdout=lf,stderr=subprocess.STDOUT).returncode
txt=ANSI.sub('',open(wd+'/log.txt').read())
rows=re.findall(r'^\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)',txt,re.M)
res={'idx':i,'rc':rc,'wall':time.time()-t,'steps':int(rows[-1][0]) if rows else -1,'grad_last':float(rows[-1][4]) if rows else float('nan'),'has_opt':os.path.exists(wd+'/start.opt.xyz'),'nan_in_log':bool(re.search(r'\bnan\b',txt.lower()))}
res['converged']=bool(rc==0 and res['has_opt'] and res['grad_last']<=5e-4)
json.dump(res,open(wd+'/res.json','w'))
