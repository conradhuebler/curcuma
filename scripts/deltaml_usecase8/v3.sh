#!/bin/bash
# usage: v3.sh dir extra-args...
U=/home/conrad/src/curcuma/build_deltaml/usecase8
d=$1; shift
cd $U/$d
s=$(date +%s)
$U/curcuma_uc8 -activelearn start.xyz -al_mode learn -method gfnff -ref_method gfn2 -activelearn.seed 1 -activelearn.al_cycles 8 -activelearn.al_temperature 300 -activelearn.al_md_time 10000 -activelearn.al_md_runs 4 -activelearn.al_dump_frequency 100 -activelearn.al_candidates 15 -activelearn.al_initial 15 -test_fraction 0 -threads 8 -no_bmt "$@" > log.txt 2>&1
echo "$(date +%H:%M:%S) V3 $d rc=$? wall=$(( $(date +%s)-s ))" >> $U/STATUS.txt
