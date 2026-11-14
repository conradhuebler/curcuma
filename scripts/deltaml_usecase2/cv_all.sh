#!/bin/bash
export SEED=${SEED:-101}
cd /home/conrad/src/curcuma/build_deltaml/usecase2
for M in M1 M2; do for f in 0 1 2 3 4; do echo "$M $f"; done; done | xargs -P 3 -L1 bash -c 'python3 -I scripts/cv_job.py $0 $1 $SEED > cvlog_$0_$1_$SEED.txt 2>&1'
echo done > cv_done.flag
