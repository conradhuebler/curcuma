# chk_labels.py dataset.json ref_method workdir : 5 structures, fresh -sp -dump_gradient vs stored labels
import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase5/scripts")
from common5 import *
import json,shutil
dsf,ref,wd=sys.argv[1:4]
S=json.load(open(dsf))['structures']; n=len(S)
idx=[int(round(i*(n-1)/4)) for i in range(5)]
mx=dict(geo=0,Eb=0,Er=0,gb=0,gr=0)
for i in idx:
    for role,meth,extra in (('base','gfnff',['-gfnff.cache_topology','false']),('ref',ref,[])):
        d=f'{wd}/{role}_{i}'; shutil.rmtree(d,ignore_errors=True); os.makedirs(d)
        open(d+'/x.xyz','w').write(xyz_text(S,[i]))
        r=subprocess.run([BIN5,'-sp','x.xyz','-method',meth,'-threads','1','-dump_gradient','g.txt','-no_bmt',*extra],cwd=d,capture_output=True,text=True)
        e,g=dumpg(d+'/g.txt')
        k='b' if role=='base' else 'r'
        key='energy_base' if role=='base' else 'energy_ref'; gk='gradient_base' if role=='base' else 'gradient_ref'
        mx['E'+k]=max(mx['E'+k],abs(e-S[i][key])); mx['g'+k]=max(mx['g'+k],np.abs(g.ravel()-np.array(S[i][gk])).max())
        print(i,role,'dE',e-S[i][key],'dg',np.abs(g.ravel()-np.array(S[i][gk])).max(),flush=True)
print('MAX',mx)
