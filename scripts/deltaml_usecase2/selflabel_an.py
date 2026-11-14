import sys; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase2/scripts')
from common2 import *
out=[]
for M in ('M1','M2'):
    d=load(M)['structures']; s=json.load(open(f'{UC}/{M}/selflabel/all.deltaml.dataset.json'))['structures']
    n=len(d)
    er=np.array([x['energy_ref'] for x in d]); es=np.array([x['energy_ref'] for x in s]); perm=[int(np.argmin(np.abs(es-er[i]))) for i in range(n)]
    assert len(set(perm))==n
    s=[s[j] for j in perm]
    fr=np.array([d[i]['energy_base'] for i in range(n)]); st=np.array([s[i]['energy_base'] for i in range(n)])
    rf=np.array([d[i]['energy_ref'] for i in range(n)]); rs=np.array([s[i]['energy_ref'] for i in range(n)])
    g=np.array([np.abs(np.array(d[i]['gradient_base'])-np.array(s[i]['gradient_base'])).max() for i in range(n)])
    x=np.array([np.abs(np.array(d[i]['geometry'])-np.array(s[i]['geometry'])).max() for i in range(n)])
    dd=(fr-st)*HA
    out.append(f'{M}: n={n} order same (max geom diff {x.max():.1e} A); self-labelled energy_base vs fresh GFN-FF: fresh-stored mean {dd.mean():+.3f} std {dd.std(ddof=1):.3f} min {dd.min():+.3f} max {dd.max():+.3f} kcal/mol (|d|>1e-6 Eh: {(np.abs(fr-st)>1e-6).sum()}/{n}); energy_ref diff max {np.abs(rf-rs).max()*HA:.1e} kcal/mol; max |dgradient_base| {g.max():.2e} Eh/A')
open(UC+'/selflabel_summary.txt','w').write('\n'.join(out)+'\n');print('\n'.join(out))
