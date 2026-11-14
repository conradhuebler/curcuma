import sys,time,shutil,os,json,subprocess; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase/scripts')
from common import *
sid=sys.argv[1]; M=sys.argv[2]
T2=UC+'/t2'
meth={'gfn2':'gfn2','gfnff':'gfnff','A':'deltaml:'+T2+'/A/train.deltaml.json','B':'deltaml:'+T2+'/B/train.deltaml.json','C':'deltaml:'+T2+'/A/train.deltaml.energy.json'}[M]
wd='%s/%s/%s'%(T2,sid,M); shutil.rmtree(wd,ignore_errors=True); os.makedirs(wd)
shutil.copy('%s/%s/start.xyz'%(T2,sid),wd+'/start.xyz')
t=time.time()
with open(wd+'/log.txt','w') as lf:
    rc=subprocess.run([BIN,'-opt','start.xyz','-method',meth,'-opt.optimizer','lbfgs','-threads','8','-no_bmt'],cwd=wd,stdout=lf,stderr=subprocess.STDOUT).returncode
wall=time.time()-t
txt=ANSI.sub('',open(wd+'/log.txt').read())
rows=re.findall(r'^\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)',txt,re.M)
res={'sid':sid,'M':M,'rc':rc,'wall':wall,'steps':int(rows[-1][0]) if rows else -1,'grad_last':float(rows[-1][4]) if rows else float('nan'),'E_last':float(rows[-1][1]) if rows else float('nan'),'has_opt':os.path.exists(wd+'/start.opt.xyz'),'nan_in_log':'nan' in txt.lower()}
if res['has_opt']:
    sd=wd+'/sp'; os.makedirs(sd); shutil.copy(wd+'/start.opt.xyz',sd+'/x.xyz')
    res['E_gfn2'],res['rc_sp']=sp_energy('x.xyz','gfn2',sd,1)
json.dump(res,open(wd+'/res.json','w'))
