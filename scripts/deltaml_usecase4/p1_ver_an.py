import sys; sys.path.insert(0,"/home/conrad/src/curcuma/build_deltaml/usecase4/scripts")
from common4 import *
import json
def L(d,t,s): return json.load(open(f'{UC4}/fits/{d}/{t}_s{s}/res.json'))
for d in ('E1','M1','M2'):
  for role in ('J','E'):
    S=range(11,16); a=[L(d,'std',s) for s in S]; b=[L(d,'nc4000',s) for s in S]
    f=lambda R,k:np.array([r[role][k] for r in R])
    def ms(x): return f'{x.mean():.4g}+-{x.std(ddof=1):.2g}'
    dm=f(b,'mae')-f(a,'mae'); dg=f(b,'grms')-f(a,'grms'); dsp=f(b,'sp')-f(a,'sp')
    print(d,role,'n=5 seeds 11-15 | std MAE',ms(f(a,'mae')),'base',ms(f(a,'mae_b')),'| nc4000 MAE',ms(f(b,'mae')),'| std g',ms(f(a,'grms')),'base',ms(f(a,'grms_b')),'nc4000 g',ms(f(b,'grms')),'| Sp std',ms(f(a,'sp')),'base',ms(f(a,'sp_b')),'nc4000',ms(f(b,'sp')),
      '| paired nc4000-std dMAE',ms(dm),f'[{int((dm<0).sum())}/5 better] dG',ms(dg),f'[{int((dg<0).sum())}/5] dSp',ms(dsp),f'[{int((dsp>0).sum())}/5]',
      '| std-base: dMAE',ms(f(a,'mae')-f(a,'mae_b')),'dG',ms(f(a,'grms')-f(a,'grms_b')),'dSp',ms(f(a,'sp')-f(a,'sp_b')), '| nc4000-base: dMAE',ms(f(b,'mae')-f(b,'mae_b')),'dG',ms(f(b,'grms')-f(b,'grms_b')),'dSp',ms(f(b,'sp')-f(b,'sp_b')),'| gradient factor base/model std',ms(f(a,'grms_b')/f(a,'grms')),'nc4000',ms(f(b,'grms_b')/f(b,'grms')))
