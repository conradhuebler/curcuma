#!/bin/bash
cd /home/conrad/src/curcuma/build_deltaml/usecase3
python3 -I - > logs/p3_jobs.txt <<'PY'
import numpy as np
for M in ('M1','M2'):
    for i in np.unique(np.round(np.linspace(0,400,60)).astype(int)): print(M,i)
PY
wc -l logs/p3_jobs.txt
xargs -P 3 -L1 bash -c 'python3 -I scripts/p3_opt.py $0 $1 >> logs/p3_opt.log 2>&1' < logs/p3_jobs.txt
echo done > logs/p3_opt_done.flag
