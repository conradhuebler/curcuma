#!/bin/bash
# usage: train_var.sh V frames.xyz   -> models/V/A, models/V/B2 (A labels internally, B2 reuses the dataset)
U=/home/conrad/src/curcuma/build_deltaml/usecase8
V=$1; X=$2
mkdir -p $U/models/$V/A $U/models/$V/B2
cd $U/models/$V/A; cp $X start.xyz
s=$(date +%s)
$U/curcuma_uc8 -activelearn start.xyz -method gfnff -ref_method gfn2 -test_fraction 0 -activelearn.seed 1 -threads 8 -no_bmt > log.txt 2>&1
echo "rc=$? wall=$(( $(date +%s)-s ))" > time.txt
echo "$(date +%H:%M:%S) $V A trained $(cat time.txt)" >> $U/STATUS.txt
DS=$U/models/$V/A/start.deltaml.dataset.json
cd $U/models/$V/B2; cp $X start.xyz
s=$(date +%s)
$U/curcuma_uc8 -activelearn start.xyz -method gfnff -ref_method gfn2 -activelearn.dataset $DS -test_fraction 0 -activelearn.seed 1 -threads 8 -no_bmt -activelearn.n_centers 4000 > log.txt 2>&1
echo "rc=$? wall=$(( $(date +%s)-s ))" > time.txt
echo "$(date +%H:%M:%S) $V B2 trained $(cat time.txt)" >> $U/STATUS.txt
