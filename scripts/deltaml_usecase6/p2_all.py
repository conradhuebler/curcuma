import sys; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase6/scripts')
from p2_run import *
for M in ('M1','PEO'):
    equil(M)
    runs=[('ff','gfnff',10000),('ff','gfnff',2000),
          ('dmlA','deltaml:'+modelpath(M,'p2A'),10000),('dmlA','deltaml:'+modelpath(M,'p2A'),2000),
          ('dmlB2','deltaml:'+modelpath(M,'p2B2'),10000),('dmlB2','deltaml:'+modelpath(M,'p2B2'),2000),
          ('gfn2','gfn2',2000),('gfn2','gfn2',1000)]
    for lab,m,t in runs: nve(M,lab,m,t)
print('DONE',flush=True)
