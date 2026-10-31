import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase4/scripts")
from common4 import *
import json
a=json.load(open(ds('M2')))['structures']; b=json.load(open(UC4+'/chk/x.deltaml.dataset.json'))['structures']
idx=list(range(0,400,20)); print(len(b))
dE=[];dEb=[];dg=[];dgb=[];dgeo=[]
for k,i in enumerate(idx):
    s,t=a[i],b[k]
    dgeo.append(np.abs(np.array(s['geometry'])-np.array(t['geometry'])).max())
    dE.append(abs(s['energy_ref']-t['energy_ref']));dEb.append(abs(s['energy_base']-t['energy_base']))
    dg.append(np.abs(np.array(s['gradient_ref'])-np.array(t['gradient_ref'])).max());dgb.append(np.abs(np.array(s['gradient_base'])-np.array(t['gradient_base'])).max())
print('geo',max(dgeo),'Eref',max(dE),'Ebase',max(dEb),'gref',max(dg),'gbase',max(dgb))
