import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase4/scripts")
from fit import *
from concurrent.futures import ThreadPoolExecutor
seeds=[int(x) for x in sys.argv[1].split(',')]; stage=sys.argv[2]
V={'std':[]}
for n in (500,2000,4000): V[f'nc{n}']=['-activelearn.n_centers',str(n)]
for l in ('1e-8','1e-4','1e-2'): V[f'lam{l}']=['-activelearn.lambda',l]
for s in ('1e-5','1e-3'): V[f'se{s}']=['-activelearn.sigma_energy',s]
for s in ('1e-4','1e-2'): V[f'sf{s}']=['-activelearn.sigma_force',s]
V['noforce']=['-activelearn.use_forces','false']
for n in (8,12): V[f'nsh{n}']=['-activelearn.acsf_n_shell',str(n)]
for n in (8,16): V[f'eta{n}']=['-activelearn.acsf_radial_eta',str(n)]
V['rc7']=['-activelearn.acsf_rcut_radial','7']
V['zeta14']=['-activelearn.acsf_zeta','1,4']; V['zeta5']=['-activelearn.acsf_zeta','1,2,4,8,16']
jobs=[]
if stage=='std':
    jobs=[(d,'std',s,[]) for d in ('E1','M1','M2') for s in seeds]
elif stage=='ks':
    for d in ('E1','M1','M2'):
        for s in seeds:
            sg=json.load(open(f'{UC4}/fits/{d}/std_s{s}/res.json'))['sigma']
            for f in (0.5,2): jobs.append((d,f'ks{f}',s,['-activelearn.kernel_sigma','%.6f'%(sg*f)]))
else:
    jobs=[(d,k,s,v) for d in ('E1','M1','M2') for k,v in V.items() if k!='std' for s in seeds]
print(len(jobs),flush=True)
def run(j):
    r=fit(*j); print(j[0],j[1],j[2],r.get('rc'),'%.0f'%r.get('wall',0),flush=True); return r
with ThreadPoolExecutor(3) as ex: list(ex.map(run,jobs))
print('DONE',flush=True)
