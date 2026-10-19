cd /home/conrad/src/curcuma/build_deltaml/usecase4
for M in M1 M2; do for k in 1 2 3 4 5 6; do echo "$M $k"; done; done | xargs -P 3 -L 1 sh -c 'python3 -I scripts/p3_md.py $0 $1' > p3_md.log 2>&1
touch p3_md_done.flag
