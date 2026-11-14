import os,re,json,subprocess,glob
import numpy as np
U='/home/conrad/src/curcuma/build_deltaml/usecase7'
BIN=U+'/curcuma_uc7'
COV={'H':0.31,'C':0.76,'O':0.66,'N':0.71}
HA=627.5094740631
def read_trj(p):
    """returns list of (t_fs, E_comment, el, X)"""
    out=[];L=open(p).read().split('\n');i=0
    while i<len(L) and L[i].strip():
        n=int(L[i]);c=L[i+1].split();t=float(c[0]) if c[0]!='**' else 0.0;E=float(c[c.index('=')+1]) if '=' in c else float('nan')
        el=[];X=[]
        for l in L[i+2:i+2+n]:
            a=l.split();el.append(a[0]);X.append([float(v) for v in a[1:4]])
        out.append((t,E,el,np.array(X)));i+=n+2
    return out
def graph(el,X,f=1.25):
    n=len(el);adj={i:[] for i in range(n)}
    for i in range(n):
        for j in range(i+1,n):
            if np.linalg.norm(X[i]-X[j])<f*(COV[el[i]]+COV[el[j]]): adj[i].append(j);adj[j].append(i)
    return adj
def backbone_torsions(el,X):
    """all heavy-atom a-b-c-d paths (a<d), from the bond graph"""
    adj=graph(el,X);hv=lambda i:el[i]!='H'
    T=set()
    for b in range(len(el)):
        if not hv(b): continue
        for c in adj[b]:
            if c<=b or not hv(c): continue
            for a in adj[b]:
                if a==c or not hv(a): continue
                for d in adj[c]:
                    if d==b or d==a or not hv(d): continue
                    T.add((a,b,c,d))
    return sorted(T)
def rotatable_torsions(el,X):
    """M1: heavy-atom torsions about single bonds that are acyclic (not in a ring) and both ends have heavy neighbours; one per bond"""
    adj=graph(el,X);n=len(el)
    def in_ring(b,c):
        # bond b-c in ring if c reachable from b without the bond
        seen={b};st=[x for x in adj[b] if x!=c]
        seen.update(st)
        while st:
            u=st.pop()
            if u==c: return True
            for w in adj[u]:
                if w not in seen: seen.add(w);st.append(w)
        return False
    T=[]
    for b in range(n):
        for c in adj[b]:
            if c<=b or el[b]=='H' or el[c]=='H': continue
            if in_ring(b,c): continue
            A=[a for a in adj[b] if a!=c and el[a]!='H'];D=[d for d in adj[c] if d!=b and el[d]!='H']
            if A and D: T.append((A[0],b,c,D[0]))
    return T
def dihedral(X,t):
    a,b,c,d=[X[i] for i in t]
    b0=a-b;b1=c-b;b2=d-c
    b1n=b1/np.linalg.norm(b1);v=b0-np.dot(b0,b1n)*b1n;w=b2-np.dot(b2,b1n)*b1n
    x=np.dot(v,w);y=np.dot(np.cross(b1n,v),w)
    return np.degrees(np.arctan2(y,x))
def rg(el,X):
    m=np.array([{'H':1.008,'C':12.011,'O':15.999,'N':14.007}[e] for e in el]);c=(m[:,None]*X).sum(0)/m.sum()
    return float(np.sqrt((m*((X-c)**2).sum(1)).sum()/m.sum()))
def w1(a,b):
    a=np.sort(np.asarray(a,float));b=np.sort(np.asarray(b,float))
    z=np.sort(np.concatenate([a,b]));dz=np.diff(z)
    Fa=np.searchsorted(a,z[:-1],side='right')/len(a);Fb=np.searchsorted(b,z[:-1],side='right')/len(b)
    return float(np.sum(np.abs(Fa-Fb)*dz))
def jsd(p,q):
    p=p/p.sum();q=q/q.sum();m=0.5*(p+q)
    def kl(x,y):
        k=x>0;return float(np.sum(x[k]*np.log2(x[k]/y[k])))
    return 0.5*kl(p,m)+0.5*kl(q,m)
BINS=np.arange(-180,181,30)
def hist(angles): return np.histogram(angles,bins=BINS)[0].astype(float)
def gfn2_sp(frames,wd,tag,threads=8):
    """frames: list of (el,X); returns energies (Eh) via curcuma multi-structure -sp"""
    os.makedirs(wd,exist_ok=True)
    f=f'{wd}/{tag}.xyz'
    with open(f,'w') as o:
        for k,(el,X) in enumerate(frames):
            o.write(f'{len(el)}\nf{k}\n'+''.join('%s %.8f %.8f %.8f\n'%(e,*x) for e,x in zip(el,X)))
    log=f'{wd}/{tag}.sp.log'
    with open(log,'w') as lf: rc=subprocess.run([BIN,'-sp',f,'-method','gfn2','-threads',str(threads),'-no_bmt'],cwd=wd,stdout=lf,stderr=subprocess.STDOUT).returncode
    txt=re.sub(r'\x1b\[[0-9;]*m','',open(log).read())
    E=[float(m.group(2)) for m in re.finditer(r'Structure\s+(\d+)\s+Single Point Energy =\s*(-?[\d.]+)',txt)]
    assert rc==0 and len(E)==len(frames),(rc,len(E),len(frames))
    return np.array(E)
