#!/bin/bash
# usage: run_md.sh part mol method_tag T ps seed dump_steps start_xyz
U=/home/conrad/src/curcuma/build_deltaml/usecase7
part=$1; mol=$2; tag=$3; T=$4; ps=$5; seed=$6; dump=$7; start=$8
case $tag in
 gfn2) M=gfn2;; gfnff) M=gfnff;;
 A) M=deltaml:$U/models/A/start.deltaml.json;;
 B2) M=deltaml:$U/models/B2/start.deltaml.json;;
 M1A) M=deltaml:$U/models/M1A/start.deltaml.json;;
 M1B2) M=deltaml:$U/models/M1B2/start.deltaml.json;;
esac
wd=$U/runs/$part/${mol}_${tag}_s$seed
[ -f $wd/done.json ] && exit 0
rm -rf $wd; mkdir -p $wd; cd $wd; cp $start md.xyz
s=$(date +%s.%N)
$U/curcuma_uc7 -md md.xyz -method $M -md.temperature $T -md.max_time $((ps*1000)) -md.time_step 1 -md.thermostat csvr -md.seed $seed -md.dump_frequency $dump -md.print_frequency 1000 -md.no_restart true -threads 8 -no_bmt > md.log 2>&1
rc=$?
w=$(echo "$(date +%s.%N)-$s"|bc)
echo "{\"part\":\"$part\",\"mol\":\"$mol\",\"tag\":\"$tag\",\"T\":$T,\"ps\":$ps,\"seed\":$seed,\"rc\":$rc,\"wall\":$w}" > done.json
echo "$(date +%H:%M:%S) done $part $mol $tag s$seed rc=$rc wall=$w" >> $U/STATUS.txt
