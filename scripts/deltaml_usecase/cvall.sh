cd /home/conrad/src/curcuma/build_deltaml/usecase
for seed in 101 202; do for f in 0 1 2 3 4; do echo "$seed $f"; done; done | xargs -P 3 -L 1 sh -c 'python3 -I scripts/cvfold.py $0 $1 > cv_$0_$1.log 2>&1'
echo done > cv_done.flag
