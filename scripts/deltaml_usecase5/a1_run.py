# a1_run.py DS seeds tag... : std and nc4000 fits for seeds
import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase5/scripts")
from fit import *
dsn=sys.argv[1]; seeds=[int(x) for x in sys.argv[2].split(',')]
V={'se1e-5':['-activelearn.sigma_energy','1e-5'],'se1e-6':['-activelearn.sigma_energy','1e-6'],'se1e-3':['-activelearn.sigma_energy','1e-3'],'lam1e-4':['-activelearn.lambda','1e-4'],'lam1e-8':['-activelearn.lambda','1e-8'],'noforce':['-activelearn.use_forces','false'],'std':[],'nc4000':['-activelearn.n_centers','4000'],'lam1e-3':['-activelearn.lambda','1e-3']}
for tag in sys.argv[3:]:
    for s in seeds:
        r=fit(dsn,tag,s,V[tag]); print(dsn,tag,s,r.get('rc'),'%.0f'%r.get('wall',0),flush=True)
