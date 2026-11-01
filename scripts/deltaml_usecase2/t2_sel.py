import sys,json; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase2/scripts')
from common2 import *
for M in ('M1','M2'):
    n=len(load(M)['structures']); sel=sorted(np.random.RandomState(7).choice(n,20,replace=False).tolist())
    rest=[i for i in range(n) if i not in sel]; os.makedirs(f'{UC}/{M}/t2',exist_ok=True)
    json.dump(sel,open(f'{UC}/{M}/t2/sel.json','w')); json.dump(rest,open(f'{UC}/{M}/t2/rest.json','w'))
