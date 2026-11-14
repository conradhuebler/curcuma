import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase5/scripts")
from common5 import *
import glob,json
seeds=[int(x) for x in sys.argv[1].split(',')]; out=sys.argv[2]
def load(d,tag,s):
    p=f'{UC5}/fits/{d}/{tag}_s{s}/res.json'
    return json.load(open(p)) if os.path.exists(p) else None
def ms(x): x=np.array(x); return f'{x.mean():.4g}+-{x.std(ddof=1) if len(x)>1 else 0:.2g}'
rows=[]
for d in ('PEO','CAF','CAFALL','PEODFT'):
    tags=sorted({os.path.basename(p).rsplit('_s',1)[0] for p in glob.glob(f'{UC5}/fits/{d}/*_s*')})
    for role in ('J','E'):
        for tag in ['std']+[t for t in tags if t!='std']:
            R=[load(d,tag,s) for s in seeds]; S0=[load(d,'std',s) for s in seeds]
            if any(r is None or 'J' not in r for r in R): continue
            def v(r,k): return r[role][k]
            mae=[v(r,'mae') for r in R]; mb=[v(r,'mae_b') for r in R]; g=[v(r,'grms') for r in R]; gb=[v(r,'grms_b') for r in R]; sp=[v(r,'sp') for r in R]; spb=[v(r,'sp_b') for r in R]
            dmae=np.array(mae)-np.array(mb); dg=np.array(g)-np.array(gb); dsp=np.array(sp)-np.array(spb)
            # paired to std
            ps=lambda k,sign=1:[(v(r,k)-v(q,k)) for r,q in zip(R,S0)]
            d_mae_s=ps('mae'); d_g_s=ps('grms'); d_sp_s=ps('sp')
            def pm(x): x=np.array(x); return f'{x.mean():+.3g}+-{x.std(ddof=1):.2g}' if len(x)>1 else f'{x.mean():+.3g}'
            nb=lambda x,sgn:int(sum(np.array(x)*sgn<0))
            rows.append(dict(ds=d,role=role,tag=tag,n=len(R),mae=ms(mae),mae_b=ms(mb),g=ms(g),gb=ms(gb),sp=ms(sp),spb=ms(spb),
              dmae=pm(dmae),dmae_better=nb(dmae,1),dg=pm(dg),dg_better=nb(dg,1),dsp=pm(dsp),dsp_better=nb(dsp,-1),
              s_dmae=pm(d_mae_s),s_dmae_b=nb(d_mae_s,1),s_dg=pm(d_g_s),s_dg_b=nb(d_g_s,1),s_dsp=pm(d_sp_s),s_dsp_b=nb(d_sp_s,-1),
              m_dmae=float(np.mean(d_mae_s)),sd_dmae=float(np.std(d_mae_s,ddof=1)) if len(R)>1 else 0,m_dg=float(np.mean(d_g_s)),sd_dg=float(np.std(d_g_s,ddof=1)) if len(R)>1 else 0,
              rel_mae=float(np.mean(mae)/np.mean(mb)),rel_g=float(np.mean(g)/np.mean(gb)),
              ratio_g=ms(np.array(gb)/np.array(g)),sigma=ms([r['sigma'] for r in R]),ncen=R[0]['ncen'],nfeat=R[0]['nfeat'],grows=R[0]['grows'],tr_mae=ms([r['J']['tr_mae'] for r in R]),wall=ms([r['wall'] for r in R])))
json.dump(rows,open(out+'.json','w'))
with open(out+'.txt','w') as f:
    for d in ('PEO','CAF','CAFALL','PEODFT'):
        for role in ('J','E'):
            f.write(f'== {d} role {role} seeds {seeds}\n')
            f.write('tag | n | MAEmod(kcal) | MAEbase | dMAE(m-b) [better] | gRMSE mod | base | dG [better] | Spearman mod | base | dSp [better] | paired-to-std dMAE [std-better] dG [better] dSp | sigma nc nfeat grows trMAE\n')
            for r in rows:
                if r['ds']==d and r['role']==role:
                    f.write(f"{r['tag']} | {r['n']} | {r['mae']} | {r['mae_b']} | {r['dmae']} [{r['dmae_better']}/{r['n']}] | {r['g']} | {r['gb']} | {r['dg']} [{r['dg_better']}/{r['n']}] | {r['sp']} | {r['spb']} | {r['dsp']} [{r['dsp_better']}/{r['n']}] | {r['s_dmae']} [{r['s_dmae_b']}] {r['s_dg']} [{r['s_dg_b']}] {r['s_dsp']} [{r['s_dsp_b']}] | {r['sigma']} {r['ncen']} {r['nfeat']} {r['grows']} {r['tr_mae']}\n")
