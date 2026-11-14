import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase4/scripts")
from common4 import *
import json,csv
out=open(UC4+'/p2_summary.txt','w'); cf=open(UC4+'/p2_results.csv','w'); cw=csv.writer(cf)
cw.writerow(['mol','idx','label','rc','steps','has_opt','abort_msg','E_gfn2_kcal_rel_R2','rmsd_R2','model_minus_gfn2_kcal','grad_last'])
def P(s): print(s); out.write(s+'\n')
LAB=['gfnff','A','A_ab0.5','A_ab0.9','S','S_ab0.5','S_ab0.9','B2','B2_ab0.5','B2_ab0.9']
for M in ('M1','M2'):
    te=json.load(open(UC4+f'/p2_split_{M}.json'))['test']
    R={}
    for l in ['gfn2']+LAB:
        for i in te:
            p=f'{UC4}/p2/{M}/{l}/{i:03d}/res.json'
            if os.path.exists(p): R[(l,i)]=json.load(open(p))
    ref={}
    for i in te:
        r=R.get(('gfn2',i))
        if r and r['has_opt'] and r['rc']==0: ref[i]=(r['E_gfn2'],readxyz(f'{UC4}/p2/{M}/gfn2/{i:03d}/start.opt.xyz')[1])
    P(f'== {M}: 20 frames, GFN2 reference opt (R2) ok {len(ref)}/20, rc!=0 {sum(1 for i in te if R.get(("gfn2",i),{}).get("rc",9)!=0)}, median steps {np.median([R[("gfn2",i)]["steps"] for i in te]):.0f}')
    D={}
    for l in LAB:
        rows=[]
        for i in te:
            r=R.get((l,i))
            if r is None: continue
            ok=r['has_opt'] and r['rc']==0 and i in ref
            dE=rm=mm=np.nan
            if ok:
                dE=(r['E_gfn2']-ref[i][0])*HA; rm=rmsd(readxyz(f'{UC4}/p2/{M}/{l}/{i:03d}/start.opt.xyz')[1],ref[i][1])
                mm=(r['E_model']-r['E_gfn2'])*HA if 'E_model' in r else np.nan
            cw.writerow([M,i,l,r['rc'],r['steps'],r['has_opt'],r['abort_msg'],dE,rm,mm,r['grad_last']])
            rows.append(dict(i=i,ok=ok,rc=r['rc'],steps=r['steps'],dE=dE,rm=rm,mm=mm,ab=r['abort_msg']))
        D[l]=rows
    P('label | ok(rc0,result)/20 | abort msg | median steps | GFN2 E(R_m)-E(R2) kcal mean/median | RMSD to R2 mean | model-GFN2 at result mean(min) | paired vs GFN-FF (common ok n): dE mean+-sd, #better | dRMSD mean #better')
    ff={r['i']:r for r in D['gfnff']}
    for l in LAB:
        rows=D[l]; ok=[r for r in rows if r['ok']]
        c=[r for r in ok if ff[r['i']]['ok']]
        dd=np.array([r['dE']-ff[r['i']]['dE'] for r in c]); dr=np.array([r['rm']-ff[r['i']]['rm'] for r in c])
        mm=np.array([r['mm'] for r in ok])
        P(f"{l} | {len(ok)}/20 | {sum(r['ab'] for r in rows)} | {np.median([r['steps'] for r in rows]):.0f} | {np.mean([r['dE'] for r in ok]) if ok else np.nan:.2f}/{np.median([r['dE'] for r in ok]) if ok else np.nan:.2f} | {np.mean([r['rm'] for r in ok]) if ok else np.nan:.2f} | {mm.mean() if len(mm) else np.nan:+.2f}({mm.min() if len(mm) else np.nan:+.1f}) | n={len(c)}: {dd.mean() if len(dd) else np.nan:+.2f}+-{dd.std(ddof=1) if len(dd)>1 else np.nan:.2f} {int((dd<0).sum())}/{len(c)} | {dr.mean() if len(dr) else np.nan:+.2f} {int((dr<0).sum())}/{len(c)}")
    # fake minima count: model energy more than 2 kcal below GFN2 at result
    for l in LAB:
        mm=np.array([r['mm'] for r in D[l] if r['ok']])
        if l!='gfnff': P(f'   {l}: spurious (model < GFN2 by > 2 kcal) {int((mm<-2).sum())}/{len(mm)}; |model-GFN2|>2: {int((abs(mm)>2).sum())}')
out.close()
