import os,re,json,subprocess,time,shutil,sys
import numpy as np
UC6='/home/conrad/src/curcuma/build_deltaml/usecase6'
BIN=UC6+'/curcuma_uc6'
UC2='/home/conrad/src/curcuma/build_deltaml/usecase2'
UC5='/home/conrad/src/curcuma/build_deltaml/usecase5'
DS={'M1':UC2+'/M1/dataset.json','M2':UC2+'/M2/dataset.json','PEO':UC5+'/peo/md0/frames.deltaml.dataset.json'}
SYM={1:'H',6:'C',7:'N',8:'O'}
ANSI=re.compile(r'\x1b\[[0-9;]*m')
HA=627.5094740631
def load(M): return json.load(open(DS[M]))
def xyz_text(S,idx,comment=None):
    out=[]
    for i in idx:
        s=S[i]; g=s['geometry']; el=s['elements']
        out.append(f"{len(el)}\n{comment or 'struct %d'%i}\n")
        for k,e in enumerate(el): out.append("%s %.10f %.10f %.10f\n"%(SYM[e],g[3*k],g[3*k+1],g[3*k+2]))
    return ''.join(out)
def testidx(M):
    N=len(load(M)['structures']); return sorted(np.random.RandomState(7).choice(N,20,replace=False).tolist())
def readxyz(p):
    L=open(p).read().split('\n'); n=int(L[0]); el=[];X=[]
    for l in L[2:2+n]:
        a=l.split(); el.append(a[0]); X.append([float(v) for v in a[1:4]])
    return el,np.array(X)
def rmsd(P,Q):
    P=P-P.mean(0);Q=Q-Q.mean(0);H=P.T@Q;U,s,Vt=np.linalg.svd(H);dd=np.sign(np.linalg.det(Vt.T@U.T))
    D=np.diag([1,1,dd]);R=Vt.T@D@U.T;return float(np.sqrt(((P@R.T-Q)**2).sum()/len(P)))
def fit(M,tag,idx,params,wdroot=UC6+'/fits'):
    wd=f'{wdroot}/{M}_{tag}'
    if os.path.exists(wd+'/res.json'): return json.load(open(wd+'/res.json'))
    shutil.rmtree(wd,ignore_errors=True); os.makedirs(wd)
    d=load(M); S=d['structures']
    sub={k:v for k,v in d.items() if k!='structures'}; sub['structures']=[S[i] for i in idx]; sub['n_structures']=len(idx)
    json.dump(sub,open(wd+'/train.dataset.json','w'))
    open(wd+'/train.xyz','w').write(xyz_text(S,[idx[0]]))
    t=time.time()
    with open(wd+'/log.txt','w') as lf:
        rc=subprocess.run([BIN,'-activelearn','train.xyz','-method','gfnff','-ref_method','gfn2','-activelearn.dataset',wd+'/train.dataset.json','-test_fraction','0','-activelearn.seed','1','-threads','8','-no_bmt']+params,cwd=wd,stdout=lf,stderr=subprocess.STDOUT).returncode
    r=dict(M=M,tag=tag,rc=rc,wall=time.time()-t,n_train=len(idx),params=params,model=os.path.exists(wd+'/train.deltaml.json'))
    os.remove(wd+'/train.dataset.json'); json.dump(r,open(wd+'/res.json','w')); return r
def modelpath(M,tag): return f'{UC6}/fits/{M}_{tag}/train.deltaml.json'
