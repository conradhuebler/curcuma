#!/bin/bash
# usage: cv_all3.sh "DS1 DS2" mode seed
cd /home/conrad/src/curcuma/build_deltaml/usecase3
MODE=$2; SEED=$3
for M in $1; do for f in 0 1 2 3 4; do echo "$M $MODE $SEED $f"; done; done | xargs -P 3 -L1 bash -c 'python3 -I scripts/cv3.py $0 $1 $2 $3 > logs/cv_$0_$1$2_f$3.txt 2>&1'
echo done > logs/cvdone_$2_$3_$(echo $1 | tr ' ' _).flag
