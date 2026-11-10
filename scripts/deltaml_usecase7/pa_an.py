import sys,os,json,csv
sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase7/scripts')
from common7 import *
MOL=sys.argv[1] if len(sys.argv)>1 else 'peo8'
PART=sys.argv[2] if len(sys.argv)>2 else 'A'
DISC=float(sys.argv[3]) if len(sys.argv)>3 else 5000.0
NSP=int(sys.argv[4]) if len(sys.argv)>4 else 150
START=U+'/'+{'peo8':'peo8_start.xyz','M1':'M1_start.xyz'}[MOL]
tags=['gfn2','gfnff','A','B2'] if MOL=='peo8' else ['gfn2','gfnff','M1A','M1B2']
seeds=[1,2,3,4] if MOL=='peo8' else [1,2]
st=read_trj(START)[0];el0,X0=st[2],st[3]
tors=backbone_torsions(el0,X0) if MOL=='peo8' else rotatable_torsions(el0,X0)
adj=graph(el0,X0);bonds=[(i,j) for i in adj for j in adj[i] if j>i]
heavy=[i for i in range(len(el0)) if el0[i]!='H']
ends=(0,23) if MOL=='peo8' else None
OUT=f'{U}/analysis_{PART}_{MOL}';os.makedirs(OUT,exist_ok=True)
D={}
for tag in tags:
  for s in seeds:
    wd=f'{U}/runs/{PART}/{MOL}_{tag}_s{s}'
    if not os.path.exists(wd+'/done.json'): continue
    meta=json.load(open(wd+'/done.json'))
    fr=[f for f in read_trj(wd+'/md.snapshots/md.trj.xyz') if f[0]>=DISC]
    # drop duplicated t entries
    seen=set();F=[]
    for f in fr:
        if f[0] in seen: continue
        seen.add(f[0]);F.append(f)
    ang=np.array([[dihedral(f[3],t) for t in tors] for f in F])
    R=np.array([rg(f[2],f[3]) for f in F])
    ee=np.array([np.linalg.norm(f[3][ends[0]]-f[3][ends[1]]) for f in F]) if ends else np.zeros(len(F))
    stretch=max(max(np.linalg.norm(f[3][i]-f[3][j])/np.linalg.norm(X0[i]-X0[j]) for i,j in bonds) for f in F)
    # new close contacts
    if tag=='gfn2': E=np.array([f[1] for f in F])
    else:
        cf=f'{OUT}/E_{tag}_s{s}.npy'
        idx=np.unique(np.round(np.linspace(0,len(F)-1,NSP)).astype(int))
        if os.path.exists(cf): E=np.load(cf)
        else:
            E=gfn2_sp([(F[k][2],F[k][3]) for k in idx],OUT,f'sp_{tag}_s{s}');np.save(cf,E)
    D[(tag,s)]=dict(ang=ang,Rg=R,ee=ee,E=E,wall=meta['wall'],ps=meta['ps'],rc=meta['rc'],n=len(F),stretch=stretch,Emd=np.array([f[1] for f in F]))
if not D: sys.exit('no runs')
# check: GFN2 SP vs MD comment energy on 5 frames of one gfn2 run
Eall=np.concatenate([D[k]['E'] for k in D if k[0]=='gfn2']);Eref=Eall.mean()
ref={}
def pool(tag,ss): 
    ks=[(tag,s) for s in ss if (tag,s) in D]
    return {q:np.concatenate([D[k][q] for k in ks]) for q in ('Rg','ee','E')}|{'ang':np.concatenate([D[k]['ang'] for k in ks])}
def metrics(a,b,eref):
    j=np.mean([jsd(hist(a['ang'][:,i]),hist(b['ang'][:,i])) for i in range(a['ang'].shape[1])])
    return dict(JSD=j,W1_E=w1(a['E']-eref,b['E']-eref),W1_Rg=w1(a['Rg'],b['Rg']),W1_ee=w1(a['ee'],b['ee']))
rows=[]
gs=[s for s in seeds if ('gfn2',s) in D]
REF=pool('gfn2',gs)
for (tag,s),d in sorted(D.items()):
    if tag=='gfn2':
        rest=pool('gfn2',[x for x in gs if x!=s]);m=metrics(d,rest,rest['E'].mean())  # noise: seed vs rest
        m['kind']='noise'
    else:
        m=metrics(d,REF,Eref);m['kind']='vsGFN2'
    m.update(tag=tag,seed=s,Erel_mean=float((d['E']-Eref).mean()),Erel_std=float(d['E'].std()),Rg_mean=float(d['Rg'].mean()),Rg_std=float(d['Rg'].std()),ee_mean=float(d['ee'].mean()),ee_std=float(d['ee'].std()),
      wall_per_ps=d['wall']/d['ps'],n_frames=d['n'],n_E=len(d['E']),max_stretch=float(d['stretch']),rc=d['rc'])
    rows.append(m)
keys=['tag','seed','kind','JSD','W1_E','W1_Rg','W1_ee','Erel_mean','Erel_std','Rg_mean','Rg_std','ee_mean','ee_std','wall_per_ps','n_frames','n_E','max_stretch','rc']
with open(f'{OUT}/metrics.csv','w',newline='') as f:
    w=csv.DictWriter(f,keys);w.writeheader();[w.writerow({k:r[k] for k in keys}) for r in rows]
# summary
L=[]
Et=lambda tag:[r for r in rows if r['tag']==tag]
L.append(f'{MOL} part {PART}: discard {DISC} fs; torsions {len(tors)}; seeds present per tag: '+', '.join(f"{t}:{len(Et(t))}" for t in tags))
L.append('mean+-std over seeds (n per tag above) of per-seed metric vs GFN2 pooled (gfn2 row = leave-one-seed-out noise)')
L.append('tag     JSD            W1_E[kcal]     W1_Rg[A]       W1_ee[A]      Erel[kcal]   Rg[A]        ee[A]       ms/ps(s)  maxstretch')
for t in tags:
    R=Et(t)
    if not R: continue
    f=lambda k,sc=1:f"{np.mean([r[k] for r in R])*sc:.4f}+-{np.std([r[k] for r in R],ddof=1)*sc:.4f}" if len(R)>1 else f"{R[0][k]*sc:.4f}"
    L.append(f"{t:6s} {f('JSD'):15s}{f('W1_E',HA):15s}{f('W1_Rg'):15s}{f('W1_ee'):15s}{f('Erel_mean',HA):13s}{f('Rg_mean'):13s}{f('ee_mean'):13s}{np.mean([r['wall_per_ps'] for r in R]):7.1f} {max(r['max_stretch'] for r in R):.2f}")
# pooled
for t in tags[1:]:
    if any(r['tag']==t for r in rows):
        P=pool(t,seeds);m=metrics(P,REF,Eref);L.append(f"pooled {t}: JSD {m['JSD']:.4f} W1_E {m['W1_E']*HA:.3f} W1_Rg {m['W1_Rg']:.4f} W1_ee {m['W1_ee']:.4f}")
# paired
def pair(a,b):
    A={r['seed']:r for r in rows if r['tag']==a};B={r['seed']:r for r in rows if r['tag']==b};ss=sorted(set(A)&set(B))
    out=[]
    for k in ('JSD','W1_E','W1_Rg','W1_ee'):
        d=np.array([A[s][k]-B[s][k] for s in ss]);sc=HA if k=='W1_E' else 1
        out.append(f"{k} {d.mean()*sc:+.4f}+-{(d.std(ddof=1)/np.sqrt(len(d)) if len(d)>1 else 0)*sc:.4f} (better {int((d<0).sum())}/{len(ss)})")
    return f"paired {a} minus {b} (negative = {a} closer to GFN2): "+'; '.join(out)
for t in tags[2:]: 
    if any(r['tag']==t for r in rows) and any(r['tag']=='gfnff' for r in rows): L.append(pair(t,'gfnff'))
open(f'{OUT}/summary.txt','w').write('\n'.join(L)+'\n')
print('\n'.join(L))
import pickle;pickle.dump((D,tors,Eref),open(f'{OUT}/D.pkl','wb'))
