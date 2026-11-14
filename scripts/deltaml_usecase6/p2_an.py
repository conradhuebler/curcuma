import sys; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase6/scripts')
from common6 import *
KC=HA
RC={'H':0.31,'C':0.76,'N':0.71,'O':0.66}
def series(wd):
    txt=ANSI.sub('',open(wd+'/md.log').read()); R=[]
    for l in txt.split('\n'):
        a=l.split()
        if len(a)>=9 and re.fullmatch(r'\d+\.\d+',a[0]):
            try: R.append([float(x) for x in a[:9]])
            except: pass
    return np.array(R),txt   # cols: t(ps) Epot <Epot> Ekin <Ekin> Etot <Etot> T <T>
def traj(wd):
    L=open(wd+'/md.snapshots/md.trj.xyz').read().split('\n'); n=int(L[0]); fr=[]
    for k in range(len(L)//(n+2)):
        b=L[k*(n+2):(k+1)*(n+2)]
        if len(b)<n+2: break
        t=float(b[1].split()[0]); el=[x.split()[0] for x in b[2:]]; X=np.array([[float(v) for v in x.split()[1:4]] for x in b[2:]])
        fr.append((t,el,X,b))
    return fr
def analyse(M,lab,fs):
    wd=f'{UC6}/p2/{M}/{lab}_{fs}'; d=json.load(open(wd+'/done.json')); S,txt=series(wd); fr=traj(wd)
    r=dict(M=M,label=lab,fs=fs,rc=d['rc'],wall=d['wall'],nrows=len(S),nan=bool(re.search(r'\bnan\b',txt.lower().replace('-nan','')) ))
    t=S[:,0]; E=S[:,5]*KC; Tt=S[:,7]; n=len(fr[0][1])
    ok=np.isfinite(E).all() and np.isfinite(S[:,1]).all()
    r['finite']=bool(ok); r['t_end_ps']=float(t[-1])
    p=np.polyfit(t,E,1); r['drift_kcal_ps']=float(p[0]); r['drift_kcal_ps_atom']=float(p[0]/n)
    nb=int(round(t[-1]-t[0])) if t[-1]-t[0]>=2 else 1
    tb=[];eb=[]
    edges=np.linspace(t[0],t[-1]+1e-9,nb+1)
    for a,b in zip(edges[:-1],edges[1:]):
        m=(t>=a)&(t<b); tb.append(t[m].mean()); eb.append(E[m].mean())
    tb=np.array(tb); eb=np.array(eb)
    if len(tb)>=4:
        pb=np.polyfit(tb,eb,1); res=eb-np.polyval(pb,tb); r['drift_se']=float(np.sqrt((res**2).sum()/(len(tb)-2)/((tb-tb.mean())**2).sum()))
    else: r['drift_se']=float('nan')
    r['etot_std_fit']=float(np.std(E-np.polyval(p,t))); r['maxdev']=float(np.max(np.abs(E-E[0])))
    r['T_mean']=float(Tt.mean()); r['T_std']=float(Tt.std()); r['T0']=float(Tt[0]); r['Etot0']=float(E[0])
    r['Epot_mean']=float(S[:,1].mean()*KC); r['Epot_std']=float(S[:,1].std()*KC); r['Ekin_std']=float(S[:,3].std()*KC)
    # bonds from first frame
    t0,el,X0,_=fr[0]; bonds=[(i,j) for i in range(n) for j in range(i+1,n) if np.linalg.norm(X0[i]-X0[j])<1.3*(RC[el[i]]+RC[el[j]])]
    d0=np.array([np.linalg.norm(X0[i]-X0[j]) for i,j in bonds]); r['nbonds']=len(bonds)
    bs=set(bonds); nonb=[(i,j) for i in range(n) for j in range(i+1,n) if (i,j) not in bs and np.linalg.norm(X0[i]-X0[j])>=1.3*(RC[el[i]]+RC[el[j]])]
    newc=0
    for tt,e,X,_ in fr:
        newc=max(newc,sum(1 for i,j in nonb if np.linalg.norm(X[i]-X[j])<1.2*(RC[el[i]]+RC[el[j]])))
    r['max_new_contacts']=newc
    allr=[]; 
    for tt,e,X,_ in fr:
        dd=np.array([np.linalg.norm(X[i]-X[j]) for i,j in bonds]); allr.append(dd)
    allr=np.array(allr); fin=allr[-1]
    r['final_dmin']=float(fin.min()); r['final_dmax']=float(fin.max()); r['traj_dmin']=float(allr.min()); r['traj_dmax']=float(allr.max())
    r['traj_max_ratio']=float((allr/d0).max()); r['traj_min_ratio']=float((allr/d0).min()); r['start_dmin']=float(d0.min()); r['start_dmax']=float(d0.max())
    r['nframes']=len(fr); r['frame_t']=[f[0] for f in fr[:3]]+[fr[-1][0]]
    return r
def pick(wd,times_fs):
    fr=traj(wd); ts=np.array([f[0] for f in fr]); sel=[]
    # trj time unit: detect (ps or fs) by max
    scale=1000.0 if ts.max()<=12 else 1.0
    for tf in times_fs:
        k=int(np.argmin(np.abs(ts*scale-tf))); sel.append(fr[k])
    return sel
def gfn2_sets(M,lab,fs,times,tag):
    wd=f'{UC6}/p2/{M}/{lab}_{fs}'; sel=pick(wd,times); es=[]
    sd=f'{UC6}/p2/{M}/sp_{lab}_{tag}'; os.makedirs(sd,exist_ok=True)
    cache=sd+'/E.json'
    if os.path.exists(cache): return json.load(open(cache))
    for k,(t,el,X,b) in enumerate(sel):
        open(f'{sd}/s{k}.xyz','w').write('\n'.join([b[0],f'{t}']+b[2:])+'\n')
        rr=subprocess.run([BIN,'-sp',f's{k}.xyz','-method','gfn2','-threads','8','-no_bmt'],cwd=sd,capture_output=True,text=True)
        m=re.search(r'Single Point Energy\s*=\s*(-?\d+\.\d+)',ANSI.sub('',rr.stdout+rr.stderr)); es.append(float(m.group(1)) if m else float('nan'))
    json.dump(es,open(cache,'w')); return es
if __name__=='__main__':
    out=[];res=[]
    for M in ('M1','PEO'):
        for lab,fs in (('ff',10000),('dmlA',10000),('dmlB2',10000),('gfn2',2000),('ff',2000),('dmlA',2000),('dmlB2',2000),('gfn2',1000)):
            if os.path.exists(f'{UC6}/p2/{M}/{lab}_{fs}/done.json'): res.append(analyse(M,lab,fs))
    json.dump(res,open(UC6+'/p2_nve.json','w'),indent=1)
    # ms per step
    ms={}
    for M in ('M1','PEO'):
        for lab,(a,b) in (('ff',(10000,2000)),('dmlA',(10000,2000)),('dmlB2',(10000,2000)),('gfn2',(2000,1000))):
            try:
                wa=json.load(open(f'{UC6}/p2/{M}/{lab}_{a}/done.json'))['wall']; wb=json.load(open(f'{UC6}/p2/{M}/{lab}_{b}/done.json'))['wall']
                ms[(M,lab)]=(wa-wb)/(a-b)*1000
            except Exception as e: pass
    for M in ('M1','PEO'):
        out.append(f'== {M}: NVE 1 fs, identical start state (equil 5 ps GFN-FF 300 K csvr), n atoms from trj')
        out.append('run          fs   n_rows finite  drift kcal/mol/ps  per atom      Etot std(fit)  maxdev  T mean+-std   (T0)   ms/step   bonds d_min/d_max end | traj | max d/d0')
        for r in res:
            if r['M']!=M: continue
            k=(M,r['label']); 
            out.append(f"{r['label']:6s} {r['fs']:6d} {r['nrows']:6d} {str(r['finite'] and not r['nan']):5s} {r['drift_kcal_ps']:+7.3f}+-{r['drift_se']:5.3f} {r['drift_kcal_ps_atom']:+9.4f} {r['etot_std_fit']:8.3f} {r['maxdev']:7.3f} {r['T_mean']:6.1f}+-{r['T_std']:5.1f} ({r['T0']:5.1f}) {ms.get(k,float('nan')) if r['fs'] in (10000,) or r['label']=='gfn2' and r['fs']==2000 else float('nan'):8.2f}  {r['final_dmin']:.3f}/{r['final_dmax']:.3f} | {r['traj_dmin']:.3f}/{r['traj_dmax']:.3f} | {r['traj_max_ratio']:.2f}  newc={r['max_new_contacts']} rc={r['rc']} start {r['start_dmin']:.2f}/{r['start_dmax']:.2f} nb={r['nbonds']}")
    # GFN2 energies of 20 structures
    out.append('')
    out.append('GFN2 single-point energy of 20 uniformly spaced MD structures [kcal/mol, relative to the GFN2-MD mean of the same window]')
    for M in ('M1','PEO'):
        full=list(range(500,10001,500)); w2=list(range(100,2001,100))
        g2=np.array(gfn2_sets(M,'gfn2',2000,w2,'w2'))*KC
        out.append(f'{M}: GFN2-MD (0-2 ps, 20 str.) mean {g2.mean():.2f} std {g2.std(ddof=1):.2f} [abs kcal/mol]')
        for lab in ('ff','dmlA','dmlB2'):
            for tag,tm in (('w2',w2),('full',full)):
                e=np.array(gfn2_sets(M,lab,10000,tm,tag))*KC
                se=np.sqrt(e.var(ddof=1)/20+g2.var(ddof=1)/20)
                out.append(f'  {lab:6s} {tag:4s}: mean {e.mean()-g2.mean():+7.2f} +- {e.std(ddof=1):5.2f} (std of the 20; SE of mean diff {se:.2f}; n=20 vs 20, autocorrelated)')
    open(UC6+'/p2_nve_summary.txt','w').write('\n'.join(out)+'\n'); print('\n'.join(out))
