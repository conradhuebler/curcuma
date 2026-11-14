#!/bin/bash
cd /home/conrad/src/curcuma/build_deltaml/usecase3
bash scripts/cv_all3.sh E1 rand 101
bash scripts/cv_all3.sh E1 rand 202
echo done > partI_done.flag
