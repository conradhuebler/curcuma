import sys,pickle,os
sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase7/scripts')
from common7 import *
U9='/home/conrad/src/curcuma/build_deltaml/usecase9';OUT=U9+'/analysis'
D=pickle.load(open(OUT+'/D9.pkl','rb'))
def pool(tag,ss): return {q:np.concatenate([D[(tag,s)][q] for s in ss]) for q in ('Rg','ee','E','ang')}
REF=pool('gfn2',[2,3,4]);Eref=REF['E'].mean()
def met(a,b,e): return np.array([np.mean([jsd(hist(a['ang'][:,i]),hist(b['ang'][:,i])) for i in range(21)]),w1(a['E']-e,b['E']-e)*HA,w1(a['Rg'],b['Rg']),w1(a['ee'],b['ee'])])
def holes(tag):
    m=[];mn=[]
    for s in (1,2,3,4):
        if tag.startswith('V0_'): wd=f'/home/conrad/src/curcuma/build_deltaml/usecase7/runs/A/peo8_{tag[3:]}_s{s}';cf=f'/home/conrad/src/curcuma/build_deltaml/usecase7/analysis_A_peo8/E_{tag[3:]}_s{s}.npy'
        elif tag.startswith('V3_'): wd=f'/home/conrad/src/curcuma/build_deltaml/usecase8/runs/{tag}_s{s}';cf=f'/home/conrad/src/curcuma/build_deltaml/usecase8/analysis/E_{tag}_s{s}.npy'
        else: wd=f'{U9}/runs/{tag}_s{s}';cf=f'{OUT}/E_{tag}_s{s}.npy'
        seen=set();F=[]
        for f in read_trj(wd+'/md.snapshots/md.trj.xyz'):
            if f[0]>=5000 and f[0] not in seen: seen.add(f[0]);F.append(f)
        idx=np.unique(np.round(np.linspace(0,len(F)-1,150)).astype(int))
        d=(np.array([F[k][1] for k in idx])-np.load(cf))*HA;m.append(d.mean());mn.append(d.min())
    return np.mean(m),np.min(mn),np.array(m),np.array(mn)
tags=['gfnff','V0_B2','V3_v3b','V3_v3c','V3_v3d']+[f'N{a}_s{b}' for a,bs in (('0',(1,2,3)),('1',(1,2,3)),('2',(1,2)),('3',(1,2,3))) for b in bs]
M={}
for t in tags: M[t]=np.array([met(pool(t,[s]),REF,Eref) for s in (1,2,3,4)])
H={t:holes(t) for t in tags if t!='gfnff'}
L=[]
L.append('tag     JSD per MD seed 1-4          mean+-std      <=0.271  W1E   W1Rg  W1ee | hole: mean(Em-E2) / deepest frame [kcal] | closer(k/4) vs V0_B2, GFN-FF')
for t in tags:
    j=M[t][:,0];h=H.get(t)
    pb=(j<M['V0_B2'][:,0]).sum();pg=(j<M['gfnff'][:,0]).sum()
    L.append(f"{t:7s} {np.round(j,3)} {j.mean():.3f}+-{j.std(ddof=1):.3f} {(j<=0.271).sum()}/4 {M[t][:,1].mean():5.1f} {M[t][:,2].mean():.2f} {M[t][:,3].mean():.2f} | "+(f"{h[0]:+.1f} / {h[1]:.1f}" if h else '-')+f" | {pb}/4 {pg}/4")
for g,ts in (('N0',['N0_s1','N0_s2','N0_s3']),('N1',['N1_s1','N1_s2','N1_s3']),('N2',['N2_s1','N2_s2']),('N3',['N3_s1','N3_s2','N3_s3']),('V3old',['V3_v3b','V3_v3c','V3_v3d'])):
    r=np.array([M[t][:,0].mean() for t in ts]);hh=[H[t][0] for t in ts];hm=[H[t][1] for t in ts]
    allj=np.concatenate([M[t][:,0] for t in ts])
    pbv=np.mean([(M[t][:,0]-M['V0_B2'][:,0]).mean() for t in ts])
    L.append(f"GROUP {g}: run-mean JSD {np.round(r,3)} mean {r.mean():.3f} sd-between-runs {r.std(ddof=1):.3f}; all MD seeds mean {allj.mean():.3f} sd {allj.std(ddof=1):.3f}; hole mean {np.mean(hh):+.1f}, deepest {min(hm):.1f}; JSD minus V0_B2 {pbv:+.3f}; W1E {np.mean([M[t][:,1].mean() for t in ts]):.1f}; runs <=0.198: {(r<=0.198).sum()}")
open(OUT+'/final9.txt','w').write('\n'.join(L)+'\n');print('\n'.join(L))
