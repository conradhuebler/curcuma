import sys,json; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase3/scripts')
from common3 import *
from concurrent.futures import ThreadPoolExecutor
M=sys.argv[1]; B=f'{UC}/{M}/p3'; S=json.load(open(B+'/min_dataset.json'))['structures']
models={'A':f'{B}/models/A/train.deltaml.json','C':f'{B}/models/A/train.deltaml.energy.json','S':f'{B}/models/S/train.deltaml.json'}
def one(i):
    wd=f'{B}/pred/{i}'; os.makedirs(wd,exist_ok=True); open(wd+'/x.xyz','w').write(xyz_text(S,[i])); r={'i':i}
    for m,p in models.items():
        e,rc,t=sp_energy('x.xyz','deltaml:'+p,wd,1,['-verbosity','3'])
        mm=re.search(r'deltaml_correction:\s*(-?\d+\.\d+(?:[eE][-+]?\d+)?)',t)
        r['E_'+m]=e; r['rc_'+m]=rc; r['corr_'+m]=float(mm.group(1)) if mm else float('nan')
    return r
with ThreadPoolExecutor(24) as ex: rows=list(ex.map(one,range(len(S))))
json.dump(rows,open(B+'/pred.json','w')); print(M,'pred',len(rows))
