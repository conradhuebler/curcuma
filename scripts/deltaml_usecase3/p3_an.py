import sys,json; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase3/scripts')
from common3 import *
out=[];csv=['M,idx,start_idx,E_ref,E_base,E_A,E_C,E_S']
rng=np.random.RandomState(0)
for M in ('M1','M2'):
    B=f'{UC}/{M}/p3'; D=json.load(open(B+'/min_dataset.json')); S=D['structures']; P=json.load(open(B+'/pred.json')); n=len(S)
    ref=np.array([s['energy_ref'] for s in S]); base=np.array([s['energy_base'] for s in S])
    bad=sum(1 for r in P for m in 'ACS' if r['rc_'+m]!=0 or not np.isfinite(r['E_'+m]))
    dev=max(abs(r['E_'+m]-r['corr_'+m]-S[r['i']]['energy_base']) for r in P for m in 'ACS')
    E={'GFNFF':base,**{m:np.array([r['E_'+m] for r in P]) for m in 'ACS'}}
    for i in range(n): csv.append(f"{M},{i},{S[i]['start_idx']},{ref[i]:.9f},{base[i]:.9f},{E['A'][i]:.9f},{E['C'][i]:.9f},{E['S'][i]:.9f}")
    rr=(ref-ref.min())*HA
    out.append(f"== {M}: 60 starts, converged {D['n_converged']}, duplicates {D['n_dup']}, unique minima n={n}; failed sp {bad}; implicit base - fresh label max {dev:.1e} Eh")
    out.append(f"Test0 minima: GFN2 rel span {np.ptp(rr):.2f} std {rr.std(ddof=1):.2f} kcal/mol; GFN-FF span {np.ptp(base)*HA:.2f} std {base.std(ddof=1)*HA:.2f}; Spearman GFN-FF vs GFN2 {sp(base,ref):.3f}; within 1/2/5 kcal of GFN2 min: {(rr<=1).sum()}/{(rr<=2).sum()}/{(rr<=5).sum()}")
    mt={m:metrics(E[m],ref,nlow=3) for m in E}
    out.append('%-6s %-8s %-8s %-6s %-9s %-11s %-9s'%('method','Spearman','pair-acc','top1','recall@3','regret[kc]','cMAE[kc]'))
    for m in E: x=mt[m]; out.append('%-6s %-8.3f %-8.3f %-6.0f %-9.2f %-11.2f %-9.2f'%(m,x['spear'],x['pairacc'],x['top1'],x['r5'],x['regret'],x['cmae']))
    # bootstrap over minima (1000 resamples) of paired differences
    bs={m:{k:[] for k in ('spear','pairacc','cmae')} for m in 'ACS'}
    for _ in range(1000):
        ix=rng.randint(0,n,n)
        if len(set(ix))<5: continue
        b0=metrics(E['GFNFF'][ix],ref[ix],nlow=3)
        for m in 'ACS':
            b=metrics(E[m][ix],ref[ix],nlow=3)
            for k in bs[m]: bs[m][k].append(b[k]-b0[k])
    c0=E['GFNFF']-ref; c0=np.abs(c0-c0.mean())*HA
    for m in 'ACS':
        d=mt[m]['spear']-mt['GFNFF']['spear']; dp=mt[m]['pairacc']-mt['GFNFF']['pairacc']
        c=E[m]-ref; c=np.abs(c-c.mean())*HA; better=(c<c0).sum()
        # pairwise (|dE|>1) paired: pairs where model right, GFN-FF wrong etc.
        iu=np.triu_indices(n,1); dE=(ref[:,None]-ref[None,:])[iu]*HA; mk=np.abs(dE)>1
        okm=np.sign((E[m][:,None]-E[m][None,:])[iu][mk])==np.sign(dE[mk]); okb=np.sign((E['GFNFF'][:,None]-E['GFNFF'][None,:])[iu][mk])==np.sign(dE[mk])
        lo,hi=np.percentile(bs[m]['spear'],[2.5,97.5]); frac=np.mean(np.array(bs[m]['spear'])>0)
        out.append(f" paired {m}-GFNFF n={n} minima: dSpearman {d:+.3f} (bootstrap 95% [{lo:+.3f},{hi:+.3f}], >0 in {100*frac:.0f}% resamples); dPairAcc {dp:+.3f} ({mk.sum()} pairs: only-{m}-right {(okm&~okb).sum()}, only-GFNFF-right {(~okm&okb).sum()}); dRecall@3 {mt[m]['r5']-mt['GFNFF']['r5']:+.2f}; dTop1 {mt[m]['top1']-mt['GFNFF']['top1']:+.0f}; dRegret {mt[m]['regret']-mt['GFNFF']['regret']:+.2f} kcal; dcMAE {mt[m]['cmae']-mt['GFNFF']['cmae']:+.2f} kcal, per-minimum |centered err| smaller in {better}/{n}")
open(UC+'/partIII_summary.txt','w').write('\n'.join(out)+'\n');open(UC+'/partIII.csv','w').write('\n'.join(csv)+'\n');print('\n'.join(out))
