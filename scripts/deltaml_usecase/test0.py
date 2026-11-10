import sys; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase/scripts')
from common import *
d=load();S=d['structures']
er=np.array([s['energy_ref'] for s in S]);eb=np.array([s['energy_base'] for s in S])
def rk(x): 
    r=np.empty(len(x));r[np.argsort(x)]=np.arange(len(x));return r
out=[]
for n,e in [('GFN2',er),('GFN-FF',eb)]:
    r=(e-e.min())*HA
    q=np.percentile(r,[10,25,50,75,90])
    out.append("%s rel E kcal/mol: min %.2f q10 %.2f q25 %.2f q50 %.2f q75 %.2f q90 %.2f max %.2f std %.2f; within 0.5/1/2 kcal: %d/%d/%d"%(n,r.min(),*q,r.max(),r.std(),(r<=.5).sum(),(r<=1).sum(),(r<=2).sum()))
pe=np.corrcoef(er,eb)[0,1];sp=np.corrcoef(rk(er),rk(eb))[0,1]
a,b=np.polyfit(eb*HA,er*HA,1);res=er*HA-(a*eb*HA+b)
out.append("Pearson %.3f Spearman %.3f; linear fit slope %.3f RMSE after fit %.3f kcal/mol; raw RMSE of relative energies (no fit) %.3f"%(pe,sp,a,np.sqrt((res**2).mean()),np.sqrt((((er-er.min())-(eb-eb.min()))*HA)**2).mean()) if False else "Pearson %.3f Spearman %.3f; linear fit slope %.3f RMSE after fit %.3f kcal/mol"%(pe,sp,a,np.sqrt((res**2).mean())))
out.append("argmin GFN2 idx %d, its GFN-FF rank %d (of 353)"%(er.argmin(),rk(eb)[er.argmin()]+1))
open(UC+'/test0.txt','w').write('\n'.join(out)+'\n');print('\n'.join(out))
