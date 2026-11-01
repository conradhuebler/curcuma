# p2_opt.py M idx label method [extra args...] : -opt from MD frame idx (default optimizer), then GFN2 sp (+model sp) at result
import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase4/scripts")
from common4 import *
import json,shutil,subprocess,time
def job(M,i,label,meth,extra=(),root='p2',xyz_override=None):
    S=json.load(open(ds(M)))['structures']
    wd=f'{UC4}/{root}/{M}/{label}/{i:03d}'
    if os.path.exists(wd+'/res.json'): return json.load(open(wd+'/res.json'))
    shutil.rmtree(wd,ignore_errors=True); os.makedirs(wd)
    open(wd+'/start.xyz','w').write(xyz_text(S,[i]))
    t=time.time()
    with open(wd+'/log.txt','w') as lf:
        rc=subprocess.run([BIN4,'-opt','start.xyz','-method',meth,'-threads','8','-no_bmt',*extra],cwd=wd,stdout=lf,stderr=subprocess.STDOUT).returncode
    wall=time.time()-t
    txt=ANSI.sub('',open(wd+'/log.txt').read())
    rows=re.findall(r'^\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)',txt,re.M)
    res={'idx':i,'label':label,'rc':rc,'wall':wall,'steps':int(rows[-1][0]) if rows else -1,'grad_last':float(rows[-1][4]) if rows else float('nan'),
         'has_opt':os.path.exists(wd+'/start.opt.xyz'),'nan':bool(re.search(r'\bnan\b',txt.lower())),
         'abort_msg':bool(re.search(r'refusing to return an unsupported energy',txt))}
    if res['has_opt']:
        sd=wd+'/sp'; os.makedirs(sd); shutil.copy(wd+'/start.opt.xyz',sd+'/x.xyz')
        res['E_gfn2'],res['rc_sp'],_=sp_energy('x.xyz','gfn2',sd,8)
        if meth.startswith('deltaml'):
            res['E_model'],res['rc_spm'],_=sp_energy('x.xyz',meth,sd,8)
        if meth=='gfnff': res['E_gfnff'],_,_=sp_energy('x.xyz','gfnff',sd,8,['-gfnff.cache_topology','false'])
    json.dump(res,open(wd+'/res.json','w'))
    return res
if __name__=='__main__':
    M,i,label,meth=sys.argv[1],int(sys.argv[2]),sys.argv[3],sys.argv[4]
    print(job(M,i,label,meth,sys.argv[5:]))
