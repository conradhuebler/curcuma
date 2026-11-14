import sys; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase6/scripts')
from common6 import *
from concurrent.futures import ThreadPoolExecutor
OVH=0.015
def opt(wd,xyz,meth):
    """run -opt in wd on xyz (file in wd); returns dict"""
    t=time.time()
    with open(wd+'/log.txt','w') as lf:
        rc=subprocess.run([BIN,'-opt',xyz,'-method',meth,'-threads','8','-no_bmt'],cwd=wd,stdout=lf,stderr=subprocess.STDOUT).returncode
    wall=time.time()-t
    txt=ANSI.sub('',open(wd+'/log.txt').read())
    rows=re.findall(r'^\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$',txt,re.M)
    base=xyz[:-4]
    return dict(rc=rc,wall=wall,wall_net=wall-OVH,steps=int(rows[-1][0]) if rows else -1,grad_last=float(rows[-1][4]) if rows else float('nan'),
        tsum=sum(float(r[5]) for r in rows),out=f'{wd}/{base}.opt.xyz' if os.path.exists(f'{wd}/{base}.opt.xyz') else None,nan=bool(re.search(r'\bnan\b',txt.lower())))
def sp(wd,xyz):
    r=subprocess.run([BIN,'-sp',xyz,'-method','gfn2','-threads','8','-no_bmt'],cwd=wd,capture_output=True,text=True)
    m=re.search(r'Single Point Energy\s*=\s*(-?\d+\.\d+)',ANSI.sub('',r.stdout+r.stderr)); return float(m.group(1)) if m else float('nan')
def route(M,i,rt,S):
    wd=f'{UC6}/p1/{M}/{i:03d}/{rt}'
    if os.path.exists(wd+'/res.json'): return json.load(open(wd+'/res.json'))
    shutil.rmtree(wd,ignore_errors=True); os.makedirs(wd+'/g2'); os.makedirs(wd+'/pre')
    res=dict(M=M,idx=i,route=rt)
    open(wd+'/x0.xyz','w').write(xyz_text(S,[i]))
    pre={'R0':None,'R0b':None,'R1':'gfnff','R2':'deltaml:'+modelpath(M,'p1A'),'R3':'deltaml:'+modelpath(M,'p1S'),'R4':'deltaml:'+modelpath(M,'p1B2')}[rt]
    src=wd+'/x0.xyz'
    if pre:
        shutil.copy(src,wd+'/pre/start.xyz'); p=opt(wd+'/pre','start.xyz',pre); res['pre']=p
        if not p['out']: res['fail']='pre'; json.dump(res,open(wd+'/res.json','w')); return res
        src=p['out']
    shutil.copy(src,wd+'/g2/start.xyz'); g=opt(wd+'/g2','start.xyz','gfn2'); res['g2']=g
    if g['out']:
        res['E']=sp(wd+'/g2','start.opt.xyz')
    res['t_total']=(res['pre']['wall_net'] if pre else 0)+g['wall_net']
    json.dump(res,open(wd+'/res.json','w')); return res
def frame(a):
    M,i=a; S=load(M)['structures']
    for rt in ('R0','R1','R2','R3','R4','R0b'):
        r=route(M,i,rt,S); print(M,i,rt,r.get('g2',{}).get('steps'),'%.1f'%r.get('t_total',-1),flush=True)
if __name__=='__main__':
    jobs=[(M,i) for M in ('M1','M2','PEO') for i in testidx(M)]
    with ThreadPoolExecutor(int(os.environ.get('WORKERS','2'))) as ex: list(ex.map(frame,jobs))
    print('DONE')
