import json, numpy as np, os, re, subprocess
UC='/home/conrad/src/curcuma/build_deltaml/usecase3'
BIN=UC+'/curcuma_uc3'
HA=627.5094740631
Z={'H':1,'C':6,'N':7,'O':8}
SYM={v:k for k,v in Z.items()}
ANSI=re.compile(r'\x1b\[[0-9;]*m')
def readxyz(p):
    L=open(p).read().split('\n'); n=int(L[0]); el=[];X=[]
    for l in L[2:2+n]:
        a=l.split(); el.append(a[0]); X.append([float(v) for v in a[1:4]])
    return el,np.array(X)
def dumpg(p):
    t=open(p).read().strip().split('\n'); e=float(t[1].split()[2])
    g=np.array([[float(v) for v in l.split()] for l in t[2:]]); return e,g
def load(M):
    return json.load(open(f'{UC}/{M}/dataset.json'))
def xyz_text(S,idx_list):
    out=[]
    for i in idx_list:
        s=S[i]; g=s['geometry']; el=s['elements']
        out.append(str(len(el))+"\nstruct %d\n"%i)
        for k,e in enumerate(el): out.append("%s %.10f %.10f %.10f\n"%(SYM[e],g[3*k],g[3*k+1],g[3*k+2]))
    return ''.join(out)
def write_train(d,idx,dirn):
    os.makedirs(dirn,exist_ok=True); S=d['structures']
    sub={k:v for k,v in d.items() if k!='structures'}
    sub['structures']=[S[i] for i in idx]; sub['n_structures']=len(idx)
    json.dump(sub,open(dirn+'/train.dataset.json','w'))
    open(dirn+'/train.xyz','w').write(xyz_text(S,idx))
def sp_energy(xyz,method,cwd,threads=1,extra=()):
    r=subprocess.run([BIN,'-sp',xyz,'-method',method,'-threads',str(threads),'-no_bmt',*extra],cwd=cwd,capture_output=True,text=True)
    t=ANSI.sub('',r.stdout+r.stderr)
    m=re.search(r'Single Point Energy\s*=\s*(-?\d+\.\d+)',t)
    return (float(m.group(1)) if m else float('nan')), r.returncode, t
def rmsd(P,Q):
    P=P-P.mean(0);Q=Q-Q.mean(0);H=P.T@Q;U,s,Vt=np.linalg.svd(H);dd=np.sign(np.linalg.det(Vt.T@U.T))
    D=np.diag([1,1,dd]);R=Vt.T@D@U.T;return float(np.sqrt(((P@R.T-Q)**2).sum()/len(P)))
def rk(v):
    r=np.empty(len(v)); r[np.argsort(v,kind='stable')]=np.arange(len(v)); return r
def sp(a,b):
    if len(a)<3: return np.nan
    return float(np.corrcoef(rk(a),rk(b))[0,1])
def metrics(E,ref,nlow=5):
    E=np.asarray(E)*HA; ref=np.asarray(ref)*HA; n=len(E)
    o=np.argsort(E,kind='stable'); orf=np.argsort(ref,kind='stable'); t=orf[0]
    dE=ref[:,None]-ref[None,:]; mE=E[:,None]-E[None,:]; iu=np.triu_indices(n,1)
    m=np.abs(dE[iu])>1.0; pa=float(np.mean(np.sign(dE[iu][m])==np.sign(mE[iu][m]))) if m.any() else np.nan
    k10=max(1,int(round(0.1*n))); top=set(orf[:nlow].tolist()); topk=set(orf[:k10].tolist())
    low=orf[:max(3,int(round(0.25*n)))]
    c=E-ref
    return dict(spear=sp(E,ref),pairacc=pa,top1=float(o[0]==t),r5=len(top&set(o[:nlow].tolist()))/nlow,r10p=len(topk&set(o[:k10].tolist()))/k10,
        regret=float(ref[o[0]]-ref[t]),cmae=float(np.mean(np.abs(c-c.mean()))),spear_low=sp(E[low],ref[low]))
KEYS=['spear','pairacc','top1','r5','r10p','regret','cmae','spear_low']
