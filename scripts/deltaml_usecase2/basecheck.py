import sys,shutil; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase2/scripts')
from common2 import *
out=[]
for M in ('M1','M2'):
    d=load(M);S=d['structures']
    for i in (0,57,133,250,399):
        wd=f'{UC}/{M}/bc/{i}'; shutil.rmtree(wd,ignore_errors=True); os.makedirs(wd)
        open(wd+'/x.xyz','w').write(xyz_text(S,[i]))
        for name in ('joint','energy'):
            mp=f'{UC}/{M}/models/full/A/train.deltaml.'+('json' if name=='joint' else 'energy.json')
            e,rc,t=sp_energy('x.xyz','deltaml:'+mp,wd,1,['-verbosity','3'])
            open(f'{wd}/log_{name}.txt','w').write(t)
            m=re.search(r'deltaml_correction:\s*(-?\d+\.\d+(?:[eE][-+]?\d+)?)',t)
            c=float(m.group(1)) if m else float('nan')
            out.append("%s idx %d %s: total %.8f corr %.8f implicit base %.8f stored %.8f diff %.3e Eh (%.4f kcal/mol)"%(M,i,name,e,c,e-c,S[i]['energy_base'],e-c-S[i]['energy_base'],(e-c-S[i]['energy_base'])*HA))
open(UC+'/basecheck.txt','w').write('\n'.join(out)+'\n');print('\n'.join(out))
