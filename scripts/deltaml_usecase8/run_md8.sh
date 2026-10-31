#!/bin/bash
# usage: run_md8.sh tag model.json seed   (model "gfnff" or "gfn2" => plain method)
U=/home/conrad/src/curcuma/build_deltaml/usecase8
tag=$1; mod=$2; seed=$3
M=deltaml:$mod
wd=$U/runs/${tag}_s$seed
[ -f $wd/done.json ] && exit 0
rm -rf $wd; mkdir -p $wd; cd $wd; cp $U/../usecase7/peo8_start.xyz md.xyz
s=$(date +%s.%N)
$U/curcuma_uc8 -md md.xyz -method $M -md.temperature 300 -md.max_time 40000 -md.time_step 1 -md.thermostat csvr -md.seed $seed -md.dump_frequency 100 -md.print_frequency 1000 -md.no_restart true -threads 8 -no_bmt > md.log 2>&1
rc=$?
w=$(echo "$(date +%s.%N)-$s"|bc)
echo "{\"tag\":\"$tag\",\"ps\":40,\"seed\":$seed,\"rc\":$rc,\"wall\":$w}" > done.json
echo "$(date +%H:%M:%S) MD done $tag s$seed rc=$rc wall=$w" >> $U/STATUS.txt
