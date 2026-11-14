cd /home/conrad/src/curcuma/build_deltaml/usecase8
for d in v3b v3c v3d; do echo "== $d"; sed 's/\x1b\[[0-9;]*m//g' $d/log.txt | grep -i "sampled\|converged" | sed 's/.*largest/largest/' | cut -c1-70 | tr '\n' ';'; echo; sed 's/\x1b\[[0-9;]*m//g' $d/log.txt | grep "model written"|cut -c1-110; done
