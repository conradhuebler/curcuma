import sys,glob; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase/scripts')
from common import *
d=load();S=d['structures'];sel=json.load(open(UC+'/t2_sel.json'))['sel']
def rd(p):
    L=open(p).read().split('\n');n=int(L[0]);el=[];X=[]
    for l in L[2:2+n]: a=l.split();el.append(a[0]);X.append([float(v) for v in a[1:4]])
    return el,np.array(X)
def rmsd(P,Q):
    P=P-P.mean(0);Q=Q-Q.mean(0);H=P.T@Q;U,s,Vt=np.linalg.svd(H);dd=np.sign(np.linalg.det(Vt.T@U.T))
    D=np.diag([1,1,dd]);R=Vt.T@D@U.T;return float(np.sqrt(((P@R.T-Q)**2).sum()/len(P)))
Ms=['gfn2','gfnff','A','B','C']; NAME={'gfn2':'GFN2(R2)','gfnff':'GFN-FF','A':'A','B':'B','C':'C'}
R={};res={}
for i in sel:
    sid='s%03d'%i; el0=S[i]['elements']; X0=np.array(S[i]['geometry']).reshape(-1,3)
    for m in Ms:
        p='%s/t2/%s/%s/res.json'%(UC,sid,m)
        if not os.path.exists(p): continue
        r=json.load(open(p)); res[(i,m)]=r
        if r['has_opt']:
            el,X=rd('%s/t2/%s/%s/start.opt.xyz'%(UC,sid,m)); assert el==[{1:"H",6:"C",7:"N"}[z] for z in el0]; R[(i,m)]=X
out=[];csv=['idx,method,rc,steps,converged,grad_last,rmsd_to_R2,dE_gfn2_kcal,wall']
g=[np.linalg.norm(S[i]['gradient_ref']) for i in sel]  # unit as stored
rs=[rmsd(X0,R[(i,'gfn2')]) if (i,'gfn2') in R else np.nan for i,X0 in [(i,np.array(S[i]['geometry']).reshape(-1,3)) for i in sel]]
out.append("Start: |g_GFN2| from dataset (stored unit, norm over 270 comps) mean %.4f max %.4f min %.4f; RMSD(start,R2) mean %.3f +- %.3f A (min %.3f max %.3f), n=%d"%(np.mean(g),np.max(g),np.min(g),np.nanmean(rs),np.nanstd(rs,ddof=1),np.nanmin(rs),np.nanmax(rs),sum(np.isfinite(rs))))
tab={}
for m in Ms:
    rm=[];de=[];ex=[];cv=[];st=[]
    for i in sel:
        r=res.get((i,m))
        if r is None: continue
        ex.append(r['rc']==0); st.append(r['steps']); cv.append(r['rc']==0 and r['has_opt'] and r['steps']<5000 and r['grad_last']<=5e-4)
        rr=rmsd(R[(i,m)],R[(i,'gfn2')]) if (i,m) in R and (i,'gfn2') in R else np.nan
        e2=(r.get('E_gfn2',np.nan)-res[(i,'gfn2')].get('E_gfn2',np.nan))*HA if (i,'gfn2') in res else np.nan
        tab[(i,m)]=(rr,e2); rm.append(rr);de.append(e2)
        csv.append('%d,%s,%d,%d,%d,%.6f,%.4f,%.4f,%.1f'%(i,m,r['rc'],r['steps'],cv[-1],r['grad_last'],rr,e2,r['wall']))
    rm=np.array(rm);de=np.array(de)
    out.append("%-9s n=%d exit0 %d/%d, converged(grad<=5e-4,steps<5000) %d, steps med %d (min %d max %d); RMSD to R2 %.3f+-%.3f A; dE(GFN2 at R_m - at R2) %.2f+-%.2f kcal/mol (median %.2f)"%(NAME[m],len(rm),sum(ex),len(ex),sum(cv),np.median(st),min(st),max(st),np.nanmean(rm),np.nanstd(rm,ddof=1),np.nanmean(de),np.nanstd(de,ddof=1),np.nanmedian(de)))
for m in ['A','B','C']:
    for k,nm in [(0,'RMSD'),(1,'dE')]:
        dl=np.array([tab[(i,m)][k]-tab[(i,'gfnff')][k] for i in sel if (i,m) in tab and (i,'gfnff') in tab]);
        out.append("paired %s - GFN-FF, %s: %+.3f+-%.3f, n=%d, %s better in %d"%(NAME[m],nm,np.nanmean(dl),np.nanstd(dl,ddof=1),np.isfinite(dl).sum(),NAME[m],(dl<0).sum()))
open(UC+'/test2_summary.txt','w').write('\n'.join(out)+'\n');open(UC+'/test2.csv','w').write('\n'.join(csv)+'\n');print('\n'.join(out))
