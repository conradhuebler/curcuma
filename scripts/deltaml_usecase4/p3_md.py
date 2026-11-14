# p3_md.py M k : GFN-FF MD 400 K, 40 ps, seed 20+k, dump 100 fs; frames t>=5000 fs -> p3/M/mdk/frames.xyz ; then internal labelling
import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase4/scripts")
from common4 import *
import shutil,subprocess,time,json
M,k=sys.argv[1],int(sys.argv[2])
wd=f'{UC4}/p3/{M}/md{k}'; shutil.rmtree(wd,ignore_errors=True); os.makedirs(wd)
S0=json.load(open(ds(M)))['structures']; START=[30,100,170,240,310,380]
open(wd+'/md.xyz','w').write(xyz_text(S0,[START[k-1]]))
t=time.time()
rc=subprocess.run([BIN4,'-md','md.xyz','-method','gfnff','-md.temperature','400','-md.max_time','40000','-md.time_step','1','-md.dump_frequency','100','-md.thermostat','csvr','-md.seed',str(20+k),'-md.print_frequency','5000','-threads','8','-no_bmt'],cwd=wd,stdout=open(wd+'/md.log','w'),stderr=subprocess.STDOUT).returncode
print(M,k,'md rc',rc,'%.0fs'%(time.time()-t),flush=True)
L=open(wd+'/md.snapshots/md.trj.xyz').read().split('\n'); n=int(L[0]); out=[];seen=set();nf=0
for i in range(len(L)//(n+2)):
    blk=L[i*(n+2):(i+1)*(n+2)]
    if len(blk)<n+2: break
    tt=float(blk[1].split()[0])
    if tt<5000 or tt in seen: continue
    seen.add(tt); out.append(blk[0]+f'\nt={tt}\n'+'\n'.join(blk[2:])+'\n'); nf+=1
open(wd+'/frames.xyz','w').write(''.join(out)); print(M,k,'frames',nf,flush=True)
t=time.time()
rc=subprocess.run([BIN4,'-activelearn','frames.xyz','-method','gfnff','-ref_method','gfn2','-test_fraction','0','-threads','8','-no_bmt'],cwd=wd,stdout=open(wd+'/label.log','w'),stderr=subprocess.STDOUT).returncode
print(M,k,'label rc',rc,'%.0fs'%(time.time()-t),flush=True)
