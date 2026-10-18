# usage: t2_job.py M idx method(gfn2|gfnff|A|S)
import sys,time,shutil,json,subprocess; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase2/scripts')
from common2 import *
M,i,m=sys.argv[1],int(sys.argv[2]),sys.argv[3]
OPT=os.environ.get('OPT','lbfgs'); RUN='run' if OPT=='lbfgs' else 'run_'+OPT
d=load(M); S=d['structures']
T2=f'{UC}/{M}/models/t2'
meth={'gfn2':'gfn2','gfnff':'gfnff','A':'deltaml:'+T2+'/A/train.deltaml.json','S':'deltaml:'+T2+'/S/train.deltaml.json'}[m]
wd=f'{UC}/{M}/t2/{RUN}/{i:03d}/{m}'; shutil.rmtree(wd,ignore_errors=True); os.makedirs(wd)
open(wd+'/start.xyz','w').write(xyz_text(S,[i]))
t=time.time()
with open(wd+'/log.txt','w') as lf:
    rc=subprocess.run([BIN,'-opt','start.xyz','-method',meth,*(['-opt.optimizer',OPT] if OPT!='auto' else []),'-threads','8','-no_bmt'],cwd=wd,stdout=lf,stderr=subprocess.STDOUT).returncode
wall=time.time()-t
txt=ANSI.sub('',open(wd+'/log.txt').read())
rows=re.findall(r'^\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)',txt,re.M)
res={'idx':i,'M':m,'rc':rc,'wall':wall,'steps':int(rows[-1][0]) if rows else -1,'grad_last':float(rows[-1][4]) if rows else float('nan'),'E_last':float(rows[-1][1]) if rows else float('nan'),'has_opt':os.path.exists(wd+'/start.opt.xyz'),'nan_in_log':bool(re.search(r'\bnan\b',txt.lower()))}
if res['has_opt']:
    sd=wd+'/sp'; os.makedirs(sd); shutil.copy(wd+'/start.opt.xyz',sd+'/x.xyz')
    res['E_gfn2'],res['rc_sp'],_=sp_energy('x.xyz','gfn2',sd,1)
    if m!='gfn2':
        res['E_model'],res['rc_spm'],_=sp_energy('x.xyz',meth,sd,1)
json.dump(res,open(wd+'/res.json','w'))
