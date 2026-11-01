#!/bin/bash
cd /home/conrad/src/curcuma/build_deltaml/usecase2
python3 -I - <<'PY' > t2_jobs.txt
import json
for M in ('M1','M2'):
    for i in json.load(open(f'{M}/t2/sel.json')):
        for m in ('gfnff','A','S','gfn2'): print(M,i,m)
PY
export OPT=${OPT:-lbfgs}
xargs -P 3 -L1 bash -c 'python3 -I scripts/t2_job.py $0 $1 $2 >> t2_jobs_$OPT.log 2>&1' < t2_jobs.txt
echo done > t2_done_$OPT.flag
