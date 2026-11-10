#!/bin/bash
U=/home/conrad/src/curcuma/build_deltaml/usecase8
cd $U/v1md; cp $U/peo8_start.xyz md.xyz
s=$(date +%s)
$U/curcuma_uc8 -md md.xyz -method gfnff -md.temperature 300 -md.max_time 40000 -md.time_step 1 -md.thermostat csvr -md.seed 11 -md.dump_frequency 100 -md.print_frequency 1000 -md.no_restart true -threads 8 -no_bmt > md.log 2>&1
echo "$(date +%H:%M:%S) V1 gfnff MD rc=$? wall=$(( $(date +%s)-s ))" >> $U/STATUS.txt
