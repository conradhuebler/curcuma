# usage: curves10.py run ...  per-cycle pooled validation curve and stop summary (Claude Generated Oct 2026)
import sys,os,re,json
import numpy as np
U='/home/conrad/src/curcuma/build_deltaml/usecase10/al'
def clean(t): return re.sub(r'\x1b\[[0-9;]*m','',t)
for d in sys.argv[1:]:
    A=np.atleast_2d(np.loadtxt(f'{U}/{d}/start.deltaml.learning_curve.dat'))
    log=clean(open(f'{U}/{d}/log.txt').read())
    m=re.search(r'converged in cycle (\d+)',log)
    nlab=len(json.load(open(f'{U}/{d}/start.deltaml.dataset.json'))['structures']) if os.path.exists(f'{U}/{d}/start.deltaml.dataset.json') else -1
    dj=json.load(open(f'{U}/{d}/done.json'))
    print(f"{d}: stop {'cycle '+m.group(1) if m else 'none (20 cycles)'}, labels {nlab}, wall {dj['wall']} s, rc {dj['rc']}")
    rows=[r for r in A if r[11]>0]
    print('  cycle(1-based) '+' '.join(f'{int(r[0])+1:5d}' for r in rows))
    print('  pool_n         '+' '.join(f'{int(r[11]):5d}' for r in rows))
    print('  cMAE           '+' '.join(f'{r[12]:5.2f}' for r in rows))
    print('  offset         '+' '.join(f'{r[13]:5.2f}' for r in rows))
    print('  hole_q10       '+' '.join(f'{r[14]:5.2f}' for r in rows))
    print('  umax           '+' '.join(f'{r[5]:5.2f}' for r in rows))
