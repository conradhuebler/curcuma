cd /home/conrad/src/curcuma/build_deltaml/usecase5
python3 -I - <<'PY'
L=open('peo/md0/frames.xyz').read().split('\n'); n=58
for k,i in enumerate([20,80,140,200,260,320],1):
    import os; os.makedirs(f'peo/md{k}',exist_ok=True)
    open(f'peo/md{k}_start.xyz','w').write('\n'.join(L[i*(n+2):(i+1)*(n+2)])+'\n')
PY
for k in 1 2 3 4 5 6; do python3 -I scripts/md_run.py peo/md$k peo/md${k}_start.xyz 400 40000 $((100+k)) 8; done
for k in 1 2 3 4 5 6; do bash scripts/label.sh /home/conrad/src/curcuma/build_deltaml/usecase5/peo/md$k gfn2; done
touch peo_mds_done.flag
