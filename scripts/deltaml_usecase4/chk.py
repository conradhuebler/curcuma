import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase4/scripts")
import json,subprocess,os
from common4 import *
d=json.load(open(ds('M2')));S=d['structures']
idx=list(range(0,400,20))
wd=UC4+'/chk'
open(wd+'/x.xyz','w').write(xyz_text(S,idx))
rc=subprocess.run([BIN4,'-activelearn','x.xyz','-method','gfnff','-ref_method','gfn2','-test_fraction','0','-threads','8','-no_bmt'],cwd=wd,stdout=open(wd+'/log.txt','w'),stderr=subprocess.STDOUT).returncode
print('rc',rc,os.listdir(wd))
