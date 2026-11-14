#!/usr/bin/env python3
"""Claude Generated (Oct 2026): compare atomic partial charges of several models for a PDMAEMA (AB)n chain.

Usage (in a directory holding q_<model>.txt written by  curcuma -sp chain.xyz -method <model> -charge N -dump_atomic_charges q_<model>.txt):
    polymer_charge_compare.py chain.xyz gfnff gfn1 gfn2
Prints the charge of the protonated and neutral side chains, N and N-H charges, and per-atom / per-chain correlations.
"""
import sys, numpy as np, collections
r={'H':.31,'C':.76,'N':.71,'O':.66}
xyz,models=sys.argv[1],sys.argv[2:]
L=open(xyz).read().split('\n')[2:]
at=[(l.split()[0],np.array(l.split()[1:4],float)) for l in L if l.strip()]
el=[a[0] for a in at];X=np.array([a[1] for a in at]);n=len(el)
D=np.linalg.norm(X[:,None]-X[None],axis=2)
adj={i:[j for j in range(n) if j!=i and D[i,j]<1.2*(r[el[i]]+r[el[j]])] for i in range(n)}
q={m:np.loadtxt(f'q_{m}.txt') for m in models}
# side chain = component of carbonyl C (2 O neighbours) after cutting the bond to its quaternary C
carb=[i for i in range(n) if el[i]=='C' and sum(el[k]=='O' for k in adj[i])==2]
side=[]
for c in carb:
    quat=[k for k in adj[c] if el[k]=='C'][0]
    seen={c};st=[c]
    while st:
        k=st.pop()
        for m in adj[k]:
            if m not in seen and not (k==c and m==quat): seen.add(m);st.append(m)
    N=[k for k in seen if el[k]=='N'][0]
    side.append((sorted(seen),len(adj[N])==4))
prot=[s for s,p in side if p]; neut=[s for s,p in side if not p]
print('side chains: %d protonated, %d neutral, atoms per chain %s/%s'%(len(prot),len(neut),len(prot[0]),len(neut[0])))
print('%-6s %24s %24s %12s'%('model','side chain (NH+)  mean [min,max]','side chain (NMe2) mean [min,max]','backbone sum'))
allside=set(sum([s for s,_ in side],[]))
for m in models:
    a=[q[m][s].sum() for s in prot];b=[q[m][s].sum() for s in neut]
    bb=sum(q[m][i] for i in range(n) if i not in allside)
    print('%-6s %8.3f [%6.3f,%6.3f]   %8.3f [%6.3f,%6.3f]   %8.3f'%(m,np.mean(a),min(a),max(a),np.mean(b),min(b),max(b),bb))
# charge on the N-H hydrogen and N
print('\nN charge / N-H charge (protonated), N charge (neutral):')
for m in models:
    Np=[k for s in prot for k in s if el[k]=='N']; Nn=[k for s in neut for k in s if el[k]=='N']
    Hn=[h for k in Np for h in adj[k] if el[h]=='H']
    print('%-6s N+ %6.3f  H(N) %6.3f  N0 %6.3f'%(m,np.mean(q[m][Np]),np.mean(q[m][Hn]),np.mean(q[m][Nn])))
def rank(v): return np.argsort(np.argsort(v)).astype(float)
print('\nper-atom correlation over %d atoms (Pearson r / Spearman rho / RMSD [e] / slope):'%n)
for i,a in enumerate(models):
    for b in models[i+1:]:
        x,y=q[a],q[b]
        print('%-6s vs %-6s r=%.3f rho=%.3f rmsd=%.3f slope=%.2f'%(a,b,np.corrcoef(x,y)[0,1],np.corrcoef(rank(x),rank(y))[0,1],np.sqrt(np.mean((x-y)**2)),np.polyfit(x,y,1)[0]))
print('\nper-side-chain sum correlation (%d chains):'%len(side))
sc={m:np.array([q[m][s].sum() for s,_ in side]) for m in models}
for i,a in enumerate(models):
    for b in models[i+1:]:
        print('%-6s vs %-6s r=%.3f'%(a,b,np.corrcoef(sc[a],sc[b])[0,1]))
