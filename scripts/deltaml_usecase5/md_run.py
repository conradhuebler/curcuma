# md_run.py wd startxyz T maxtime_fs seed nthreads [tmin_fs] [nmax] : GFN-FF MD (csvr, 1 fs, dump 100 fs); frames t>=tmin -> wd/frames.xyz (first nmax)
import sys,os,subprocess,time,shutil
wd,start,T,mt,seed,nt=sys.argv[1:7]; tmin=float(sys.argv[7]) if len(sys.argv)>7 else 5000; nmax=int(sys.argv[8]) if len(sys.argv)>8 else 10**9
BIN='/home/conrad/src/curcuma/build_deltaml/usecase5/curcuma_uc5'
os.makedirs(wd,exist_ok=True); shutil.copy(start,wd+'/md.xyz')
t=time.time()
rc=subprocess.run([BIN,'-md','md.xyz','-method','gfnff','-md.temperature',T,'-md.max_time',mt,'-md.time_step','1','-md.dump_frequency','100','-md.thermostat','csvr','-md.seed',seed,'-md.print_frequency','5000','-threads',nt,'-no_bmt'],cwd=wd,stdout=open(wd+'/md.log','w'),stderr=subprocess.STDOUT).returncode
print(wd,'md rc',rc,'%.0fs'%(time.time()-t),flush=True)
L=open(wd+'/md.snapshots/md.trj.xyz').read().split('\n'); n=int(L[0]); out=[];seen=set()
for i in range(len(L)//(n+2)):
    blk=L[i*(n+2):(i+1)*(n+2)]
    if len(blk)<n+2: break
    tt=float(blk[1].split()[0])
    if tt<tmin or tt in seen: continue
    seen.add(tt); out.append(blk[0]+f'\nt={tt}\n'+'\n'.join(blk[2:])+'\n')
    if len(out)>=nmax: break
open(wd+'/frames.xyz','w').write(''.join(out)); print(wd,'frames',len(out),flush=True)
