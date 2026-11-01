cd /home/conrad/src/curcuma/build_deltaml/usecase4
for M in M1 M2; do python3 -I scripts/p3_opt.py $M; done
for M in M1 M2; do python3 -I scripts/p3_train.py $M; done
for M in M1 M2; do python3 -I scripts/p3_eval.py $M; done
touch p3_chain_done.flag
