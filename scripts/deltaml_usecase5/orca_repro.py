import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase5/scripts")
from common5 import *
import json,shutil
S=json.load(open(UC5+'/caf/md/frames.deltaml.dataset.json'))['structures']; i=37
G=[]
for r in range(2):
    d=f'{UC5}/chk/repro{r}'; shutil.rmtree(d,ignore_errors=True); os.makedirs(d)
    open(d+'/x.xyz','w').write(xyz_text(S,[i]))
    subprocess.run([BIN5,'-sp','x.xyz','-method','r2scan-3c','-threads','1','-dump_gradient','g.txt','-no_bmt'],cwd=d,capture_output=True)
    G.append(dumpg(d+'/g.txt'))
print('fresh vs fresh dE',G[0][0]-G[1][0],'dg',np.abs(G[0][1]-G[1][1]).max())
print('fresh vs stored dg',np.abs(G[0][1].ravel()-np.array(S[i]['gradient_ref'])).max(),'|g| max',np.abs(G[0][1]).max())
