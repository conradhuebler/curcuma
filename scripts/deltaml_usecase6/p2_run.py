# Part 2 driver: equilibrate (GFN-FF 300 K csvr 5 ps), then NVE runs from the same restart state for each method
import sys; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase6/scripts')
from common6 import *
def equil(M):
    wd=f'{UC6}/p2/{M}/equil'
    if os.path.exists(wd+'/md.snapshots/md.final.json'): return wd
    shutil.rmtree(wd,ignore_errors=True); os.makedirs(wd)
    S=load(M)['structures']; open(wd+'/start.xyz','w').write(xyz_text(S,[200]))
    with open(wd+'/md.log','w') as lf:
        rc=subprocess.run([BIN,'-md','start.xyz','-method','gfnff','-md.temperature','300','-md.max_time','5000','-md.time_step','1','-md.thermostat','csvr','-md.dump_frequency','1000','-md.print_frequency','500','-md.no_restart','true','-threads','8','-no_bmt'],cwd=wd,stdout=lf,stderr=subprocess.STDOUT).returncode
    print(M,'equil rc',rc,flush=True); return wd
def nve(M,label,method,maxfs):
    wd=f'{UC6}/p2/{M}/{label}_{maxfs}'
    if os.path.exists(wd+'/done.json'): return json.load(open(wd+'/done.json'))
    shutil.rmtree(wd,ignore_errors=True); os.makedirs(wd)
    eq=json.load(open(f'{UC6}/p2/{M}/equil/start.snapshots/start.final.json'))
    st=eq['MD'] if 'MD' in eq else eq; st['method']=method; st['thermostat']='none'; st['MaxTime']=float(maxfs); st['currentStep']=0.0; st['T']=300.0
    json.dump(st,open(wd+'/rs.json','w'))
    shutil.copy(f'{UC6}/p2/{M}/equil/start.xyz',wd+'/md.xyz')
    t=time.time()
    with open(wd+'/md.log','w') as lf:
        rc=subprocess.run([BIN,'-md','md.xyz','-method',method,'-md.restart_file','rs.json','-md.thermostat','none','-md.temperature','300','-md.max_time',str(maxfs),'-md.time_step','1','-md.dump_frequency','10','-md.print_frequency','10','-threads','8','-no_bmt'],cwd=wd,stdout=lf,stderr=subprocess.STDOUT).returncode
    r=dict(M=M,label=label,method=method,maxfs=maxfs,rc=rc,wall=time.time()-t)
    json.dump(r,open(wd+'/done.json','w')); print(r,flush=True); return r
