import sys,time,shutil,subprocess; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase/scripts')
from common import *
T2=UC+'/t2'
cfg=[('gfnff','gfnff',8,200,1000),('A','deltaml:'+T2+'/A/train.deltaml.json',8,200,1000),('B','deltaml:'+T2+'/B/train.deltaml.json',8,100,500),
     ('gfn2','gfn2',8,20,60),('A_1thr','deltaml:'+T2+'/A/train.deltaml.json',1,100,500),('gfnff_1thr','gfnff',1,100,500)]
REP=5; rows=[]; csv=['method,threads,rep,n_short,n_long,t_short,t_long,ms_per_step,final_time_ps_long']
def run(meth,thr,n,rep,tag):
    wd=UC+'/t3/w'; shutil.rmtree(wd,ignore_errors=True); os.makedirs(wd); shutil.copy(UC+'/t2/s003/start.xyz',wd)
    t=time.time()
    r=subprocess.run([BIN,'-md','start.xyz','-method',meth,'-md.max_time',str(n),'-md.dump_frequency','1000000','-md.print_frequency',str(n),'-md.write_xyz','false','-threads',str(thr),'-no_bmt'],cwd=wd,capture_output=True,text=True)
    dt=time.time()-t; txt=ANSI.sub('',r.stdout)
    tt=re.findall(r'^\s+(\d+\.\d+)\s+-?\d+\.\d+\s',txt,re.M)
    return dt,(float(tt[-1]) if tt else -1),r.returncode
res={}
for rep in range(REP):
    for name,meth,thr,ns,nl in cfg:
        ts,fs,rc1=run(meth,thr,ns,rep,'s'); tl,fl,rc2=run(meth,thr,nl,rep,'l')
        ms=(tl-ts)/(nl-ns)*1000
        res.setdefault(name,[]).append(ms)
        csv.append('%s,%d,%d,%d,%d,%.3f,%.3f,%.3f,%.4f'%(name,thr,rep,ns,nl,ts,tl,ms,fl))
        if rc1 or rc2 or abs(fl-nl/1000)>1e-6: csv.append('# WARNING rc %d %d final_time %.4f expected %.4f'%(rc1,rc2,fl,nl/1000))
open(UC+'/test3.csv','w').write('\n'.join(csv)+'\n')
out=[]
for k,v in res.items(): out.append("%-11s ms/step(E+grad+MD) %.2f +- %.2f (n=%d)"%(k,np.mean(v),np.std(v,ddof=1),len(v)))
open(UC+'/test3_summary.txt','w').write('\n'.join(out)+'\n');print('\n'.join(out))
