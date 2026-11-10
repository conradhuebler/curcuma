import sys,glob,json;sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase7/scripts')
from common7 import *
OUT=U+'/analysis_B_peo8'
R=[json.load(open(f)) for f in sorted(glob.glob(OUT+'/opt/*/*/res.json'))]
cf=OUT+'/E2_all.npy'
if os.path.exists(cf): E2=np.load(cf)
else: E2=gfn2_sp([(r['el'],np.array(r['X'])) for r in R],OUT,'all_sp');np.save(cf,E2)
emin=E2.min();rel=(E2-emin)*HA
tors=backbone_torsions(R[0]['el'],np.array(R[0]['X']))
def state(X):
    a=[dihedral(np.array(X),t) for t in tors]
    return ''.join('t' if abs(x)>120 else ('p' if x>0 else 'm') for x in a)
S=[state(r['X']) for r in R]
L=['Part B supplement (not in spec): GFN2 energy of each quench structure and distinct backbone-torsion states (21 torsions, t/p/m bins at 120 deg)','tag   n  GFN2rel[kcal]: min  p10  median | distinct states | states shared with other tags | lowest-5 mean']
tags=['gfn2','gfnff','B2']
st={t:set(s for s,r in zip(S,R) if r['tag']==t) for t in tags}
for t in tags:
    i=[k for k,r in enumerate(R) if r['tag']==t];e=np.sort(rel[i])
    oth=set().union(*[st[u] for u in tags if u!=t])
    L.append(f"{t:6s} {len(i):4d} {e[0]:7.2f} {np.percentile(e,10):7.2f} {np.median(e):7.2f} | {len(st[t]):4d} | {len(st[t]&oth):3d} | {e[:5].mean():.2f}")
for t in ('gfnff','B2'):
    per=[]
    for s in (1,2,3,4):
        i=[k for k,r in enumerate(R) if r['tag']==t and r['seed']==s];e=np.sort(rel[i]);per.append((e[0],np.median(e),len(set(S[k] for k in i)),(e<=2.0).sum()))
    L.append(f"{t} per seed (lowest, median, distinct states, n<=2kcal of global lowest): "+'; '.join(f"({a:.2f},{b:.2f},{c},{d})" for a,b,c,d in per))
    # seed-wise: lowest of each seed within the pooled
L.append('rel energies measured against global lowest GFN2 energy over all %d quench structures (GFN-FF minima, GFN2 single points); n per tag as listed'%len(R))
open(OUT+'/supplement.txt','w').write('\n'.join(L)+'\n');print('\n'.join(L))
