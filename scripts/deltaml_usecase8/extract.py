# usage: extract.py trj.xyz out.xyz tmin  -> frames with t>=tmin (unique t), written without energy comment
import sys
sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase7/scripts')
from common7 import *
fr=read_trj(sys.argv[1]);tmin=float(sys.argv[3]);seen=set();F=[]
for f in fr:
    if f[0]>=tmin and f[0] not in seen: seen.add(f[0]);F.append(f)
with open(sys.argv[2],'w') as o:
    for t,E,el,X in F: o.write(f'{len(el)}\nt={t:.0f}\n'+''.join('%s %.8f %.8f %.8f\n'%(e,*x) for e,x in zip(el,X)))
print(len(F),'frames',F[0][0],F[-1][0])
