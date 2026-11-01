# eval_on.py model.json dataset.json outjson workdir : model (total, base+corr) vs stored reference on every structure of dataset
import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase5/scripts")
from common5 import *
import json,shutil
mp,dsf,outp,wd=sys.argv[1:5]
S=json.load(open(dsf))['structures']
os.makedirs(wd,exist_ok=True)
E=[];G=[]
for i,s in enumerate(S):
    open(wd+'/x.xyz','w').write(xyz_text(S,[i]))
    r=subprocess.run([BIN5,'-sp','x.xyz','-method','deltaml:'+mp,'-threads','1','-dump_gradient','g.txt','-no_bmt'],cwd=wd,capture_output=True)
    e,g=dumpg(wd+'/g.txt'); E.append(e); G.append(g.ravel().tolist())
Er=np.array([s['energy_ref'] for s in S]); Eb=np.array([s['energy_base'] for s in S]); E=np.array(E)
Gr=np.array([s['gradient_ref'] for s in S]); Gb=np.array([s['gradient_base'] for s in S]); G=np.array(G)
def relmae(x): d=(x-Er)*HA; return float(np.mean(np.abs(d-d.mean())))
res=dict(n=len(S),mae_model=relmae(E),mae_base=relmae(Eb),grms_model=float(np.sqrt(np.mean((G-Gr)**2))),grms_base=float(np.sqrt(np.mean((Gb-Gr)**2))),
  sp_model=sp(E,Er),sp_base=sp(Eb,Er),offset_model=float(np.mean((E-Er))*HA),offset_base=float(np.mean((Eb-Er))*HA))
json.dump(res,open(outp,'w')); print(res)
