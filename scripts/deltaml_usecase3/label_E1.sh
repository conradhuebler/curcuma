#!/bin/bash
cd /home/conrad/src/curcuma/build_deltaml/usecase3
ls -d E1/fr/* | xargs -P 24 -I{} bash scripts/label_one.sh {}
echo done > E1_label_done.flag
