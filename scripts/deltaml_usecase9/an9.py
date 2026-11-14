# usage: an8.py [tag ...]   tags: V0_* baselines from usecase7, new from usecase8/runs
import sys,os,json,csv,pickle
sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase7/scripts')
import common7
from common7 import *
U8='/home/conrad/src/curcuma/build_deltaml/usecase9'
common7.BIN=U8+'/curcuma_uc9'
OUT=U8+'/analysis';os.makedirs(OUT,exist_ok=True)
DISC=5000.0;NSP=150
st=read_trj(U+'/peo8_start.xyz')[0];el0,X0=st[2],st[3]
tors=backbone_torsions(el0,X0);adj=graph(el0,X0);bonds=[(i,j) for i in adj for j in adj[i] if j>i]
SEEDS=[1,2,3,4];REFS=[2,3,4]
def rundir(tag,s):
    if tag in('gfn2','gfnff'): return f'{U}/runs/A/peo8_{tag}_s{s}'
    if tag in('V0_A','V0_B2'): return f'{U}/runs/A/peo8_{tag[3:]}_s{s}'
    return f'{U8}/runs/{tag}_s{s}'
def cachepath(tag,s):
    if tag=='gfnff': return f'{U}/analysis_A_peo8/E_gfnff_s{s}.npy'
    if tag in('V0_A','V0_B2'): return f'{U}/analysis_A_peo8/E_{tag[3:]}_s{s}.npy'
    return f'{OUT}/E_{tag}_s{s}.npy'
D={}
tags=sys.argv[1:] or ['gfn2','gfnff','V0_A','V0_B2']
for tag in tags:
  for s in SEEDS:
    wd=rundir(tag,s)
    if not os.path.exists(wd+'/done.json'): continue
    meta=json.load(open(wd+'/done.json'))
    seen=set();F=[]
    for f in read_trj(wd+'/md.snapshots/md.trj.xyz'):
        if f[0]>=DISC and f[0] not in seen: seen.add(f[0]);F.append(f)
    ang=np.array([[dihedral(f[3],t) for t in tors] for f in F]);R=np.array([rg(f[2],f[3]) for f in F])
    ee=np.array([np.linalg.norm(f[3][0]-f[3][23]) for f in F])
    stretch=max(max(np.linalg.norm(f[3][i]-f[3][j])/np.linalg.norm(X0[i]-X0[j]) for i,j in bonds) for f in F)
    if tag=='gfn2': E=np.array([f[1] for f in F])
    else:
        cf=cachepath(tag,s)
        if os.path.exists(cf): E=np.load(cf)
        else:
            idx=np.unique(np.round(np.linspace(0,len(F)-1,NSP)).astype(int))
            E=gfn2_sp([(F[k][2],F[k][3]) for k in idx],OUT,f'sp_{tag}_s{s}');np.save(cf,E)
    D[(tag,s)]=dict(ang=ang,Rg=R,ee=ee,E=E,wall=meta['wall'],ps=meta['ps'],rc=meta['rc'],n=len(F),stretch=stretch)
pk=OUT+'/D9.pkl'
old=pickle.load(open(pk,'rb')) if os.path.exists(pk) else {}
old.update(D);D=old;pickle.dump(D,open(pk,'wb'))
def pool(tag,ss):
    ks=[(tag,s) for s in ss]
    return {q:np.concatenate([D[k][q] for k in ks]) for q in ('Rg','ee','E','ang')}
Eref=pool('gfn2',REFS)['E'].mean()
def metrics(a,b,eref):
    j=np.mean([jsd(hist(a['ang'][:,i]),hist(b['ang'][:,i])) for i in range(a['ang'].shape[1])])
    return dict(JSD=j,W1_E=w1(a['E']-eref,b['E']-eref)*HA,W1_Rg=w1(a['Rg'],b['Rg']),W1_ee=w1(a['ee'],b['ee']))
REF=pool('gfn2',REFS);rows=[]
for (tag,s),d in sorted(D.items()):
    if tag=='gfn2':
        if s in REFS:
            rest=pool('gfn2',[x for x in REFS if x!=s]);m=metrics(d,rest,rest['E'].mean());m['kind']='noise2v2'
        else: m=metrics(d,REF,Eref);m['kind']='noise1v3'
    else: m=metrics(d,REF,Eref);m['kind']='vs234'
    m.update(tag=tag,seed=s,Erel=float((d['E']-Eref).mean()*HA),Rg=float(d['Rg'].mean()),ee=float(d['ee'].mean()),wall_ps=d['wall']/d['ps'],stretch=float(d['stretch']),rc=d['rc'])
    rows.append(m)
keys=['tag','seed','kind','JSD','W1_E','W1_Rg','W1_ee','Erel','Rg','ee','wall_ps','stretch','rc']
with open(OUT+'/metrics9.csv','w',newline='') as f:
    w=csv.DictWriter(f,keys);w.writeheader();[w.writerow({k:r[k] for k in keys}) for r in rows]
def ms(rs,k): 
    v=[r[k] for r in rs]; return f"{np.mean(v):.3f}+-{np.std(v,ddof=1) if len(v)>1 else 0:.3f}"
L=['tag            n  JSD            W1_E[kcal]    W1_Rg[A]      W1_ee[A]     Erel[kcal]   Rg[A]       ee[A]      s/ps  maxstretch']
allt=sorted({r['tag'] for r in rows},key=lambda t:(t!='gfn2',t!='gfnff',t))
for t in allt:
    for kind in sorted({r['kind'] for r in rows if r['tag']==t}):
        rs=[r for r in rows if r['tag']==t and r['kind']==kind]
        L.append(f"{t+('('+kind+')' if t=='gfn2' else ''):15s}{len(rs)} {ms(rs,'JSD'):14s} {ms(rs,'W1_E'):13s} {ms(rs,'W1_Rg'):13s} {ms(rs,'W1_ee'):12s} {ms(rs,'Erel'):12s} {ms(rs,'Rg'):11s} {ms(rs,'ee'):10s} {np.mean([r['wall_ps'] for r in rs]):5.1f} {max(r['stretch'] for r in rs):.2f}")
def pair(a,b):
    A={r['seed']:r for r in rows if r['tag']==a and r['kind']=='vs234'};B={r['seed']:r for r in rows if r['tag']==b and r['kind']=='vs234'};ss=sorted(set(A)&set(B))
    if not ss: return None
    o=[]
    for k in('JSD','W1_E','W1_Rg','W1_ee'):
        d=np.array([A[s][k]-B[s][k] for s in ss]);o.append(f"{k} {d.mean():+.3f}+-{(d.std(ddof=1)/np.sqrt(len(d)) if len(d)>1 else 0):.3f} ({int((d<0).sum())}/{len(ss)})")
    return f"{a} minus {b}: "+'; '.join(o)
L.append('paired by seed (negative = first closer to GFN2 pool 2-4; mean+-SE; closer in k/n seeds)')
for t in allt:
    if t in('gfn2','gfnff'): continue
    for b in ('gfnff',)+(('V0_'+t.split('_',1)[1],) if t.startswith('V') and not t.startswith('V0') else ()):
        if b!=t:
            p=pair(t,b)
            if p:L.append(p)
# per-torsion pooled JSD
for t in allt:
    if t=='gfn2' or not all((t,s) in D for s in SEEDS): continue
    P=pool(t,SEEDS);pj=[jsd(hist(P['ang'][:,i]),hist(REF['ang'][:,i])) for i in range(21)]
    L.append(f"pooled-4-seeds per-torsion JSD {t}: median {np.median(pj):.3f} max {max(pj):.3f} mean {np.mean(pj):.3f}")
open(OUT+'/summary9.txt','w').write('\n'.join(L)+'\n');print('\n'.join(L))
