import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase5/scripts")
from common5 import *
import json,shutil
S=json.load(open(UC5+'/caf/md/frames.deltaml.dataset.json'))['structures']
d=f'{UC5}/chk/repro2'; shutil.rmtree(d,ignore_errors=True); os.makedirs(d)
for i in (36,37):
    open(d+'/x.xyz','w').write(xyz_text(S,[i]))
    subprocess.run([BIN5,'-sp','x.xyz','-method','r2scan-3c','-threads','1','-dump_gradient','g.txt','-no_bmt'],cwd=d,capture_output=True)
    e,g=dumpg(d+'/g.txt'); print(i,'dE vs stored',e-S[i]['energy_ref'],'dg',np.abs(g.ravel()-np.array(S[i]['gradient_ref'])).max())
