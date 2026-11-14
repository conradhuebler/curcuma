import sys,glob,os,re,json
import numpy as np
U='/home/conrad/src/curcuma/build_deltaml/usecase9/al'
def clean(t): return re.sub(r'\x1b\[[0-9;]*m','',t)
rows=[]
for d in sorted(os.listdir(U)):
    cf=f'{U}/{d}/start.deltaml.learning_curve.dat'
    if not os.path.exists(cf): continue
    A=np.loadtxt(cf);A=np.atleast_2d(A)
    log=clean(open(f'{U}/{d}/log.txt').read())
    conv='not converged' not in log and ('converged' in log.lower())
    m=re.findall(r'converged[^\n]*',log,flags=re.I)
    dj=json.load(open(f'{U}/{d}/done.json')) if os.path.exists(f'{U}/{d}/done.json') else {}
    v=A[:,7:11]; vv=v[~np.isnan(v[:,0])]
    ncyc=len(A)
    print(f"{d:8s} cycles {ncyc:2d} labels(last n_train) {int(A[-1,1])} wall {dj.get('wall')} conv={conv} | umax first/mean/last {A[0,5]:.2f}/{A[1:,5].mean():.2f}/{A[-1,5]:.2f}")
    print(f"         val MAE first/mean/last {vv[0,0]:.2f}/{vv[:,0].mean():.2f}/{vv[-1,0]:.2f}; signed {vv[0,1]:.2f}/{vv[:,1].mean():.2f}/{vv[-1,1]:.2f}; hole {vv[0,2]:.2f}/{vv[:,2].mean():.2f}/{vv[-1,2]:.2f}; mean over last5 MAE {vv[-5:,0].mean():.2f} hole {vv[-5:,2].mean():.2f}")
    print('         ',[l[:150] for l in m[-1:]])
