import sys,os,json,csv,shutil,time
sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase7/scripts')
from common7 import *
from concurrent.futures import ThreadPoolExecutor
OUT=U+'/analysis_B_peo8';os.makedirs(OUT,exist_ok=True)
RUNS=[('gfn2',1,10),('gfnff',1,40),('gfnff',2,40),('gfnff',3,40),('gfnff',4,40),('B2',1,40),('B2',2,40),('B2',3,40),('B2',4,40)]
def kabsch(P,Q):
    P=P-P.mean(0);Q=Q-Q.mean(0);H=P.T@Q;Uu,s,Vt=np.linalg.svd(H);d=np.sign(np.linalg.det(Vt.T@Uu.T))
    D=np.diag([1,1,d]);R=Vt.T@D@Uu.T;return float(np.sqrt(((P@R.T-Q)**2).sum()/len(P)))
def quench(job):
    tag,s,k,el,X=job
    wd=f'{OUT}/opt/{tag}_s{s}/f{k:02d}';rf=wd+'/res.json'
    if os.path.exists(rf): return json.load(open(rf))
    shutil.rmtree(wd,ignore_errors=True);os.makedirs(wd)
    open(wd+'/f.xyz','w').write(f'{len(el)}\nf\n'+''.join('%s %.8f %.8f %.8f\n'%(e,*x) for e,x in zip(el,X)))
    t=time.time()
    with open(wd+'/log.txt','w') as lf: rc=subprocess.run([BIN,'-opt','f.xyz','-method','gfnff','-threads','8','-no_bmt'],cwd=wd,stdout=lf,stderr=subprocess.STDOUT).returncode
    r=dict(tag=tag,seed=s,k=k,rc=rc,wall=time.time()-t,ok=False)
    if rc==0 and os.path.exists(wd+'/f.opt.xyz'):
        f=read_trj(wd+'/f.opt.xyz')[-1] if False else None
        L=open(wd+'/f.opt.xyz').read().split('\n');n=int(L[0]);c=L[1].split()
        r['E']=float(c[c.index('=')+1]);r['el']=[l.split()[0] for l in L[2:2+n]];r['X']=[[float(v) for v in l.split()[1:4]] for l in L[2:2+n]];r['ok']=True
    json.dump(r,open(rf,'w'));return r
jobs=[]
for tag,s,ps in RUNS:
    wd=f'{U}/runs/B/peo8_{tag}_s{s}'
    if not os.path.exists(wd+'/done.json'): print('missing',wd);continue
    seen=set()
    for f in read_trj(wd+'/md.snapshots/md.trj.xyz'):
        if f[0]>=1000 and f[0] not in seen and abs(f[0]/1000-round(f[0]/1000))<1e-6:
            seen.add(f[0]);jobs.append((tag,s,int(round(f[0]/1000)),f[2],f[3]))
print(len(jobs),'quench jobs')
with ThreadPoolExecutor(3) as ex: res=list(ex.map(quench,jobs))
bad=[r for r in res if not r['ok']];print('failed opts',len(bad))
# dedupe
def same(a,b): return abs(a['E']-b['E'])*HA<0.2 and kabsch(np.array(a['X']),np.array(b['X']))<0.3
def dedupe(items):
    kept=[]
    for r in sorted(items,key=lambda r:r['E']):
        for q in kept:
            if same(r,q): q['count']+=1;break
        else: r=dict(r);r['count']=1;kept.append(r)
    return kept
M={}
for tag in ('gfn2','gfnff','B2'):
    M[tag]=dedupe([r for r in res if r['ok'] and r['tag']==tag])
    per={s:len(dedupe([r for r in res if r['ok'] and r['tag']==tag and r['seed']==s])) for s in (1,2,3,4)}
    print(tag,'frames',sum(1 for r in res if r['tag']==tag),'distinct',len(M[tag]),'per seed',per)
# union with membership
U_=[]
for tag in M:
    for r in M[tag]:
        for u in U_:
            if same(r,u['rep']): u['found'].add(tag);u['E2']=None;break
        else: U_.append(dict(rep=r,found={tag}))
print('union',len(U_))
frames=[(u['rep']['el'],np.array(u['rep']['X'])) for u in U_]
E2=gfn2_sp(frames,OUT,'union_sp')
for u,e in zip(U_,E2): u['E2']=float(e)
emin=min(E2)
walls={}
for tag,s,ps in RUNS:
    p=f'{U}/runs/B/peo8_{tag}_s{s}/done.json'
    if os.path.exists(p): walls[tag]=walls.get(tag,0)+json.load(open(p))['wall']
rows=[]
for tag in M:
    mine=[u for u in U_ if tag in u['found']]
    # all minima of this tag are in union; GFN2 energies
    es=sorted(u['E2'] for u in mine)
    new=[u for u in mine if u['found']=={tag}]
    rows.append(dict(tag=tag,distinct=len(mine),emin_kcal=(min(es)-emin)*HA,low2=sum(1 for e in es if (e-emin)*HA<=2.0),unique_to_this=len(new),unique_frac=len(new)/len(mine),
        wall_s=walls[tag],wall_per_min=walls[tag]/len(mine),frames=sum(1 for r in res if r['tag']==tag)))
with open(OUT+'/minima.csv','w') as f:
    f.write('tag,distinct,emin_rel_all_kcal,n_within_2kcal_of_global_lowest,unique_to_tag,unique_frac,sampling_wall_s,wall_s_per_minimum,frames\n')
    for r in rows: f.write(f"{r['tag']},{r['distinct']},{r['emin_kcal']:.3f},{r['low2']},{r['unique_to_this']},{r['unique_frac']:.2f},{r['wall_s']:.0f},{r['wall_per_min']:.1f},{r['frames']}\n")
with open(OUT+'/union.csv','w') as f:
    f.write('E_gfn2_Eh,E_rel_kcal,E_gfnff,found,count_per_tag\n')
    for u in sorted(U_,key=lambda u:u['E2']): f.write(f"{u['E2']:.8f},{(u['E2']-emin)*HA:.3f},{u['rep']['E']:.8f},{'+'.join(sorted(u['found']))},\n")
print(open(OUT+'/minima.csv').read())
# also: GFN2 energy of every quenched structure distribution (diagnostic)
E=[(u['E2']-emin)*HA for u in U_];print('union GFN2 rel energies: min 0, median %.2f, max %.2f'%(np.median(E),max(E)))
