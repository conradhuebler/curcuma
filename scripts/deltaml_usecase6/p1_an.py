import sys; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase6/scripts')
from common6 import *
import csv
RT=['R0','R0b','R1','R2','R3','R4']
rows=[]; out=[]
def get(M,i,rt):
    p=f'{UC6}/p1/{M}/{i:03d}/{rt}/res.json'
    return json.load(open(p)) if os.path.exists(p) else None
for M in ('M1','M2','PEO'):
    ids=testidx(M); R={}
    for i in ids:
        r0=get(M,i,'R0'); ref=readxyz(r0['g2']['out']) [1] if r0 and r0['g2']['out'] else None
        for rt in RT:
            r=get(M,i,rt)
            if r is None or 'g2' not in r: continue
            g=r['g2']; conv=(g['rc']==0 and g['grad_last']<=5e-4 and g['out'] is not None)
            d=dict(M=M,idx=i,route=rt,steps_g2=g['steps'],t_g2=g['wall_net'],t_pre=r['pre']['wall_net'] if 'pre' in r else 0.0,steps_pre=r['pre']['steps'] if 'pre' in r else 0,
                   t_total=r.get('t_total',float('nan')),conv=conv,grad=g['grad_last'],E=r.get('E',float('nan')))
            if conv and ref is not None and r0.get('E') is not None:
                d['dE']=(r['E']-r0['E'])*HA; d['rmsd']=rmsd(readxyz(g['out'])[1],ref)
                d['basin']=(d['rmsd']<0.3 and abs(d['dE'])<0.2)
            R[(i,rt)]=d; rows.append(d)
    out.append(f'== {M} (n frames={len(ids)}; test idx seed 7)')
    out.append('route  n  conv  g2steps med/mean   t_pre mean  t_g2 med/mean   t_tot med/mean   basin  dE med[kcal]  rmsd med')
    for rt in RT:
        L=[R[(i,rt)] for i in ids if (i,rt) in R]
        if not L: continue
        f=lambda k:np.array([x[k] for x in L],float)
        bas=[x['basin'] for x in L if 'basin' in x]
        out.append(f"{rt:4s} {len(L):3d} {sum(x['conv'] for x in L):3d}  {np.median(f('steps_g2')):6.0f}/{f('steps_g2').mean():6.1f}  {f('t_pre').mean():7.2f}  {np.median(f('t_g2')):6.1f}/{f('t_g2').mean():6.1f}  {np.median(f('t_total')):6.1f}/{f('t_total').mean():6.1f}  {sum(bas)}/{len(bas)}  {np.nanmedian([x.get('dE',np.nan) for x in L]):7.3f}  {np.nanmedian([x.get('rmsd',np.nan) for x in L]):6.3f}")
    out.append('endpoint vs R0: n  mean dE  | dE<-0.2 | |dE|<0.2 | dE>+0.2 | rmsd<0.3 | mean rmsd | pre-opt steps mean')
    for rt in RT[1:]:
        L=[R[(i,rt)] for i in ids if (i,rt) in R and 'dE' in R[(i,rt)]]
        dE=np.array([x['dE'] for x in L]); rm=np.array([x['rmsd'] for x in L])
        out.append(f"  {rt:4s} {len(L):3d} {dE.mean():+7.3f} | {(dE<-0.2).sum():3d} | {(abs(dE)<0.2).sum():3d} | {(dE>0.2).sum():3d} | {(rm<0.3).sum():3d} | {rm.mean():5.2f} | {np.mean([x['steps_pre'] for x in L]):6.1f}")
    out.append('paired (route - ref): d_ttot mean+-se [s] (n faster/n), d_g2steps mean+-se (n fewer/n)')
    for rt,ref in (('R0b','R0'),('R1','R0'),('R2','R0'),('R3','R0'),('R4','R0'),('R2','R1'),('R3','R1'),('R4','R1')):
        ii=[i for i in ids if (i,rt) in R and (i,ref) in R]
        if not ii: continue
        dt=np.array([R[(i,rt)]['t_total']-R[(i,ref)]['t_total'] for i in ii]); ds=np.array([R[(i,rt)]['steps_g2']-R[(i,ref)]['steps_g2'] for i in ii])
        se=lambda x:x.std(ddof=1)/np.sqrt(len(x)) if len(x)>1 else np.nan
        out.append(f"  {rt} - {ref}: dt {dt.mean():+7.2f}+-{se(dt):5.2f} ({int((dt<0).sum())}/{len(ii)})  dsteps {ds.mean():+7.1f}+-{se(ds):5.1f} ({int((ds<0).sum())}/{len(ii)})  ratio mean t {np.mean([R[(i,rt)]['t_total'] for i in ii])/np.mean([R[(i,ref)]['t_total'] for i in ii]):.3f}")
open(UC6+'/p1_summary.txt','w').write('\n'.join(out)+'\n')
keys=['M','idx','route','steps_pre','steps_g2','t_pre','t_g2','t_total','conv','grad','E','dE','rmsd','basin']
with open(UC6+'/p1_results.csv','w',newline='') as f:
    w=csv.DictWriter(f,keys,extrasaction='ignore'); w.writeheader(); w.writerows(rows)
print('\n'.join(out))
