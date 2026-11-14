import sys,os
M=sys.argv[1]
L=open(f'/home/conrad/src/curcuma/build_deltaml/usecase2/{M}/md/opt.snapshots/opt.trj.xyz').read().split('\n'); n=int(L[0])
seen=set();k=0
os.makedirs(f'/home/conrad/src/curcuma/build_deltaml/usecase2/{M}/fr',exist_ok=True)
nfr=len(L)//(n+2)
for i in range(nfr):
    blk=L[i*(n+2):(i+1)*(n+2)]
    if len(blk)<n+2: break
    t=float(blk[1].split()[0])
    if t<5000 or t in seen: continue
    seen.add(t)
    d=f'/home/conrad/src/curcuma/build_deltaml/usecase2/{M}/fr/{k:03d}'; os.makedirs(d,exist_ok=True)
    open(d+'/s.xyz','w').write(blk[0]+f'\nt={t}\n'+'\n'.join(blk[2:])+'\n'); k+=1
print(M,k)
