#!/bin/bash
# usage: al.sh dir seed binary extra-args...   (setup of usecase9 al.sh; Claude Generated Oct 2026)
U=/home/conrad/src/curcuma/build_deltaml/usecase10
d=$1; seed=$2; bin=$3; shift 3
mkdir -p $U/al/$d; cd $U/al/$d; cp $U/peo8_start.xyz start.xyz
s=$(date +%s)
$U/$bin -activelearn start.xyz -al_mode learn -method gfnff -ref_method gfn2 -activelearn.seed $seed -activelearn.al_cycles 20 -activelearn.al_temperature 300 -activelearn.al_md_time 10000 -activelearn.al_md_runs 4 -activelearn.al_dump_frequency 100 -activelearn.al_candidates 15 -activelearn.al_initial 15 -activelearn.al_uncertainty 0.3 -test_fraction 0 -threads 8 -no_bmt "$@" > log.txt 2>&1
rc=$?
echo "{\"dir\":\"$d\",\"rc\":$rc,\"wall\":$(( $(date +%s)-s ))}" > done.json
echo "$(date +%H:%M:%S) AL $d rc=$rc wall=$(( $(date +%s)-s ))" >> $U/STATUS.txt
