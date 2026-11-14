#!/bin/bash
cd /home/conrad/src/curcuma/build_deltaml/usecase2
for M in M1 M2; do ls -d $M/fr/* ; done | xargs -P 24 -I{} bash scripts/label_one.sh {}
echo done > label_done.flag
