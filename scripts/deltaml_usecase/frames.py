import json,math,glob,numpy as np
K=627.5094740631
d=json.load(open('/home/conrad/src/curcuma/build_deltaml/ablation_run/conf.deltaml.dataset.json'))['structures']
fresh={int(l.split()[0]):float(l.split()[1]) for l in open('fresh_gfnff.txt')}
folds=[]
for f in sorted(glob.glob('../cv/s*_f*/result.json')):
    folds.append(json.load(open(f))['rows'])
def spear(a,b):
    def rk(v):
        o=np.argsort(v,kind='stable'); r=np.empty(len(v)); r[o]=np.arange(len(v)); return r
    ra,rb=rk(a),rk(b); ra-=ra.mean(); rb-=rb.mean(); return float(ra@rb/np.sqrt((ra@ra)*(rb@rb)))
def metrics(E,ref):
    E=np.asarray(E); ref=np.asarray(ref); n=len(E)
    tbest=int(np.argmin(ref)); order=np.argsort(E,kind='stable')
    rank=int(np.where(order==tbest)[0][0])+1
    t5=set(np.argsort(ref,kind='stable')[:5].tolist())
    r5=len(t5&set(order[:5].tolist()))/5; r10=len(t5&set(order[:10].tolist()))/5
    reg=(ref[order[0]]-ref[tbest])*K
    return dict(rank=rank,top1=float(order[0]==tbest),r5=r5,r10=r10,regret=reg,spear=spear(E,ref))
def run(frame):
    res={m:[] for m in ('base','A','B','C')}
    for rows in folds:
        idx=[r['idx'] for r in rows]; ref=[r['E_ref'] for r in rows]
        for m in res:
            if frame=='L':   # labelling-consistent: dataset base + learned correction
                if m=='base': E=[r['E_base'] for r in rows]
                else: E=[r['E_base']+(r['E_'+m]-fresh[r['idx']]) for r in rows]
            else:            # deployment: fresh GFN-FF (+ correction) as returned by -sp
                if m=='base': E=[fresh[r['idx']] for r in rows]
                else: E=[r['E_'+m] for r in rows]
            res[m].append(metrics(E,ref))
    return res
for frame,name in (('L','labelling-consistent (dataset base + correction)'),('D','deployment (fresh GFN-FF + correction)')):
    res=run(frame); print('\n== frame',frame,name,' n folds =',len(folds))
    for m in res:
        g=lambda k:(np.mean([x[k] for x in res[m]]),np.std([x[k] for x in res[m]]))
        print('%-5s rank %5.1f  top1 %.2f  R@5 %.2f+-%.2f  R@10 %.2f  regret %.2f+-%.2f  spearman %.2f+-%.2f'%(m,g('rank')[0],g('top1')[0],g('r5')[0],g('r5')[1],g('r10')[0],g('regret')[0],g('regret')[1],g('spear')[0],g('spear')[1]))
    for m in ('A','B','C'):
        ds=np.array([a['spear']-b['spear'] for a,b in zip(res[m],res['base'])]); dr=np.array([a['r5']-b['r5'] for a,b in zip(res[m],res['base'])])
        print('   %s - base: dSpearman %+.3f+-%.3f (%d/%d better)  dR@5 %+.2f+-%.2f (%d/%d better)'%(m,ds.mean(),ds.std(),int((ds>0).sum()),len(ds),dr.mean(),dr.std(),int((dr>0).sum()),len(dr)))
# correction itself: how large is  E_X - E_fresh  vs the true remaining error  E_ref - E_fresh ?
rows=[r for f in folds for r in f]
print('\nstd of true delta in deployment frame  (E_ref - E_fresh):  %.3f kcal/mol'%(np.std([ (r['E_ref']-fresh[r['idx']]) for r in rows])*K))
print('std of true delta in labelling frame (E_ref - E_ds):       %.3f kcal/mol'%(np.std([ (r['E_ref']-r['E_base']) for r in rows])*K))
