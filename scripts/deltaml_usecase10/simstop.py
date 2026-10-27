# usage: simstop.py run ...  replay stop rules on recorded pooled curves (Claude Generated Oct 2026)
import sys,numpy as np
U='/home/conrad/src/curcuma/build_deltaml/usecase10/al'
runs=sys.argv[1:]
C={d:np.atleast_2d(np.loadtxt(f'{U}/{d}/start.deltaml.learning_curve.dat')) for d in runs}
def stop(A,te,th,mode,pat=2,uthr=0.3):
    st=0
    for r in A:
        if r[11]<30: st=0; continue
        v=r[12]<=te and r[14]>=-th; u=r[5]<=uthr
        if mode=='validation': met=v
        elif mode=='both': met=v and u
        st=st+1 if met else 0
        if mode=='both_now':
            vs=getattr(stop,'vs',0)
        if st>=pat: return int(r[0])+1
    return None
def stop_now(A,te,th,pat=2,uthr=0.3):
    vs=0
    for r in A:
        if r[11]<30: vs=0; continue
        vs=vs+1 if (r[12]<=te and r[14]>=-th) else 0
        if vs>=pat and r[5]<=uthr: return int(r[0])+1
    return None
for te,th in [(2.0,3.0),(2.0,3.5),(2.5,3.5),(2.5,4.0),(3.0,4.0),(2.2,3.5)]:
    print(f"te {te} th {th}: validation {[stop(C[d],te,th,'validation') for d in runs]}  both {[stop(C[d],te,th,'both') for d in runs]}  val-streak+unc-now {[stop_now(C[d],te,th) for d in runs]}")
