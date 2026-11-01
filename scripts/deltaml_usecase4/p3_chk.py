import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase4/scripts")
from common4 import *
import json
for M in ('M1','M2'):
    D=[json.load(open(f'{UC4}/p3/{M}/md{k}/frames.deltaml.dataset.json'))['structures'] for k in range(1,7)]
    print(M,[len(d) for d in D])
    G=[np.array([s['geometry'] for s in d]) for d in D]
    E=[np.array([s['energy_ref'] for s in d])*HA for d in D]
    print(' GFN2 E mean/std kcal per run',[(round(e.mean()-E[0].mean(),1),round(e.std(),1)) for e in E])
    print(' max|dx| between run1 and run k at frame 100/200',[(round(float(np.abs(G[0][100]-g[100]).max()),2),round(float(np.abs(G[0][200]-g[200]).max()),2)) for g in G[1:]])
    allE=np.concatenate(E); print(' GFN2 span',allE.max()-allE.min())
