cd /home/conrad/src/curcuma/build_deltaml/usecase8
for d in v3a v3half v3b; do echo "== $d"; sed 's/\x1b\[[0-9;]*m//g' $d/log.txt | grep -i "sampled\|converged" | sed 's/activelearn: //' | cut -c1-140; echo "label calls: $(sed 's/\x1b\[[0-9;]*m//g' $d/log.txt | grep -c 'labelling')"; sed 's/\x1b\[[0-9;]*m//g' $d/log.txt | grep "model written"|cut -c1-120; done
U=$PWD; rm -f jobs2.txt; for d in v3a v3half v3b; do for s in 1 2 3 4; do echo "V3_$d $U/$d/start.deltaml.json $s" >> jobs2.txt; done; done
nohup sh -c "while pgrep -f 'xargs -P 2 -L 1 bash scripts/run_md8.sh' >/dev/null; do sleep 5; done; xargs -P 3 -L 1 bash scripts/run_md8.sh < jobs2.txt" > xargs2.log 2>&1 &
