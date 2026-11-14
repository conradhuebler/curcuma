import sys,glob; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase2/scripts')
from common2 import *
OPT=sys.argv[1] if len(sys.argv)>1 else 'lbfgs'; RUN='run' if OPT=='lbfgs' else 'run_'+OPT
out=[f'#### Test 2, optimizer={OPT}  (geometry metrics only over runs with exit 0 AND final gradient <= 5e-4 Eh/Bohr; R2 must converge too)'];csv=['M,idx,method,rc,steps,converged,grad_last,rmsd_to_R2,dE_gfn2_kcal,model_minus_gfn2_kcal,wall']
for M in ('M1','M2'):
    sel=json.load(open(f'{UC}/{M}/t2/sel.json')); S=load(M)['structures']
    res={};R={};cv={}
    for i in sel:
        for m in ('gfn2','gfnff','A','S'):
            p=f'{UC}/{M}/t2/{RUN}/{i:03d}/{m}/res.json'
            if os.path.exists(p):
                r=json.load(open(p)); res[(i,m)]=r
                cv[(i,m)]=bool(r['rc']==0 and r['has_opt'] and r['grad_last']<=5e-4)
                if r['has_opt']: R[(i,m)]=readxyz(f'{UC}/{M}/t2/{RUN}/{i:03d}/{m}/start.opt.xyz')[1]
    ok2=[i for i in sel if cv.get((i,'gfn2'))]
    out.append(f'== {M}: runs finished {len(res)}/80; R2 (GFN2) converged for {len(ok2)}/20 starts')
    tab={}
    for m in ('gfn2','gfnff','A','S'):
        rr=[];de=[];sm=[];stp=[];ex=0;c=0;N=0;nn=0
        for i in sel:
            r=res.get((i,m))
            if r is None: continue
            N+=1; ex+=r['rc']==0; c+=cv[(i,m)]; stp.append(r['steps']); nn+=r['nan_in_log']
            a=b=s=np.nan
            if cv[(i,m)] and i in ok2:
                a=rmsd(R[(i,m)],R[(i,'gfn2')]); b=(r['E_gfn2']-res[(i,'gfn2')]['E_gfn2'])*HA
                if m in 'AS': s=(r['E_model']-r['E_gfn2'])*HA
                rr.append(a);de.append(b);sm.append(s)
            tab[(i,m)]=(a,b,s)
            csv.append(f'{M},{i},{m},{r["rc"]},{r["steps"]},{int(cv[(i,m)])},{r["grad_last"]},{a:.4f},{b:.4f},{s:.3f},{r["wall"]:.0f}')
        f=lambda x:('%.2f+-%.2f (n=%d)'%(np.mean(x),np.std(x,ddof=1),len(x))) if len(x)>1 else 'n<2'
        ps=[x for x in stp if x>0]
        out.append(f'{m:6s} exit0 {ex}/{N}, converged {c}/{N}, NaN-log {nn}, steps med {np.median(ps) if ps else -1:.0f} (max {max(stp)}); RMSD to R2 {f(rr)} A; dE_GFN2 {f(de)} kcal/mol, median {np.median(de) if de else np.nan:.2f}; model-GFN2 at R_m {f([x for x in sm if np.isfinite(x)]) if m in "AS" else "-"} kcal/mol')
    for m in ('A','S'):
        for k,nm in ((0,'RMSD [A]'),(1,'dE_GFN2 [kcal/mol]')):
            dl=np.array([tab[(i,m)][k]-tab[(i,'gfnff')][k] for i in sel if (i,m) in tab and (i,'gfnff') in tab])
            dl=dl[np.isfinite(dl)]
            out.append(f' paired {m}-GFNFF {nm}: '+(f'{dl.mean():+.2f}+-{dl.std(ddof=1):.2f}, n={len(dl)} (both+R2 converged), {m} better in {(dl<0).sum()}' if len(dl)>1 else f'n={len(dl)}, not computable'))
open(f'{UC}/test2_summary_{OPT}.txt','w').write('\n'.join(out)+'\n');open(f'{UC}/test2_{OPT}.csv','w').write('\n'.join(csv)+'\n');print('\n'.join(out))
