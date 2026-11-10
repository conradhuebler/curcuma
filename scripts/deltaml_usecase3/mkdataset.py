import sys,json; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase3/scripts')
from common3 import *
M=sys.argv[1]
S=[];dirs=sorted(os.listdir(f'{UC}/{M}/fr'))
for k in dirs:
    p=f'{UC}/{M}/fr/{k}'; el,X=readxyz(p+'/s.xyz')
    assert open(p+'/base/rc.txt').read().strip()=='0' and open(p+'/ref/rc.txt').read().strip()=='0', p
    eb,gb=dumpg(p+'/base/g_base.txt'); er,gr=dumpg(p+'/ref/g_ref.txt')
    S.append(dict(charge=0,cycle=-1,elements=[Z[e] for e in el],energy_base=eb,energy_ref=er,geometry=X.flatten().tolist(),
      gradient_base=gb.flatten().tolist(),gradient_ref=gr.flatten().tolist(),has_gradients=True,spin=0.0))
d={'format':'curcuma-deltaml-dataset','n_structures':len(S),'structures':S,'version':1}
json.dump(d,open(f'{UC}/{M}/dataset.json','w'))
print(M,len(S),np.isfinite([s['energy_ref'] for s in S]).all(),np.isfinite([s['energy_base'] for s in S]).all())
