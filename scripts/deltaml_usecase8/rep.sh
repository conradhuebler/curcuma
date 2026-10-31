cd /home/conrad/src/curcuma/build_deltaml/usecase8
(bash scripts/v3.sh v3c -activelearn.al_cycles 12 -activelearn.al_uncertainty 0.3) &
(bash scripts/v3.sh v3d -activelearn.al_cycles 12 -activelearn.al_uncertainty 0.3) &
wait
U=$PWD; rm -f jobs3.txt; for d in v3c v3d; do for s in 1 2 3 4; do echo "V3_$d $U/$d/start.deltaml.json $s" >> jobs3.txt; done; done
xargs -P 3 -L 1 bash scripts/run_md8.sh < jobs3.txt
python3 -I scripts/an8.py V3_v3c V3_v3d > an_v3cd.log 2>&1
echo "$(date +%H:%M:%S) replicates analysed" >> STATUS.txt
