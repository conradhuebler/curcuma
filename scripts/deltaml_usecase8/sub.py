import sys
sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase7/scripts')
from common7 import *
fr=read_trj(sys.argv[1]);seen=set();F=[]
for f in fr:
    if f[0]>=5000 and f[0] not in seen: seen.add(f[0]);F.append(f)
F=F[::3]
open(sys.argv[2],'w').write(''.join(f'{len(el)}\nt={t:.0f}\n'+''.join('%s %.8f %.8f %.8f\n'%(e,*x) for e,x in zip(el,X)) for t,E,el,X in F));print(len(F))
