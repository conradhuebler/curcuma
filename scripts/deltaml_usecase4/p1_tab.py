import json,sys,numpy as np
rows=json.load(open(sys.argv[1]+'.json'))
role=sys.argv[2] if len(sys.argv)>2 else 'J'
import math
for d in ('E1','M1','M2'):
    R=[r for r in rows if r['ds']==d and r['role']==role]
    std=[r for r in R if r['tag']=='std'][0]
    print(f'== {d} {role}: std MAE {std["mae"]} (base {std["mae_b"]}) gRMSE {std["g"]} (base {std["gb"]}) Sp {std["sp"]} (base {std["spb"]})')
    print('tag | MAE | gRMSE | Sp | dMAE_vs_std(sd) | dG_vs_std(sd)e-3 | dSp_vs_std(sd) | noise flags(MAE,G,Sp) | score')
    for r in R:
        if r['tag']=='std': continue
        t=2.776/math.sqrt(r['n'])
        fl=''.join('N' if abs(m)<t*s else '-' for m,s in ((r['m_dmae'],r['sd_dmae']),(r['m_dg'],r['sd_dg'])))
        print(f"{r['tag']} | {r['mae'].split('+')[0]} | {r['g'].split('+')[0]} | {r['sp'].split('+')[0]} | {r['s_dmae']} | {r['s_dg']} | {r['s_dsp']} | {fl}")
