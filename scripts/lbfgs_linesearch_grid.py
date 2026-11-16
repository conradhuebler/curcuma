#!/usr/bin/env python3
"""Grid of -opt runs over (method, structure, lbfgs_max_line_search). Claude Generated (Nov 2026).

usage: lbfgs_linesearch_grid.py BINARY OUTDIR RESULT_JSON   (needs numpy; about 140 runs, 5 in parallel)
Measurement behind the decision on lbfgs_max_line_search (TODO.md, docs/debt_2026-11/PARAMS.md H5).
Reports per (method, structure): convergence, steps and final energy for each line-search value.
A run counts as converged when the process exits 0 and prints no 'did not converge'.
Perturbed variants: coordinates plus gaussian noise (sigma 0.10 A, fixed seed), written to OUTDIR/pert.
"""
import concurrent.futures as cf
import json
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np

BIN, OUT, RES = sys.argv[1], sys.argv[2], sys.argv[3]
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'test_cases', 'structures') + '/'
VALUES = ['2', '3', '5', '20']
os.makedirs(OUT + '/pert', exist_ok=True)


def perturbed(rel, sigma=0.10, seed=7):
    src = ROOT + rel
    L = open(src).read().split('\n')
    n = int(L[0])
    rng = np.random.RandomState(seed)
    X = np.array([[float(v) for v in L[2 + i].split()[1:4]] for i in range(n)]) + rng.normal(0, sigma, (n, 3))
    name = OUT + '/pert/' + os.path.basename(rel).replace('.xyz', '') + '.pert.xyz'
    with open(name, 'w') as f:
        f.write('%d\nperturbed %.2f\n' % (n, sigma))
        for i in range(n):
            f.write('%s %.10f %.10f %.10f\n' % (L[2 + i].split()[0], *X[i]))
    return name


small = ['small/ethane.xyz', 'small/benzene.xyz', 'medium/benzoic-acid.xyz', 'medium/caffeine.xyz',
         'medium/methyl-a-d-galactopyranoside.xyz', 'medium/fi-pyridine.xyz', 'medium/peo8-chain.h-capped.xyz',
         'clusters/water8-cluster.xyz']
big = ['clusters/aaa-bglc.conf96.xyz', 'clusters/macrocycle-bgal.xyz']
combos = []  # (method, path, extra args, timeout)
for s in small + big:
    combos.append(('gfnff', ROOT + s, [], 900))
for s in ['small/ethane.xyz', 'medium/caffeine.xyz', 'medium/methyl-a-d-galactopyranoside.xyz', 'clusters/macrocycle-bgal.xyz']:
    combos.append(('gfnff', perturbed(s), [], 900))
for s in ['small/ethane.xyz', 'small/benzene.xyz', 'medium/benzoic-acid.xyz', 'medium/caffeine.xyz', 'clusters/water8-cluster.xyz']:
    combos.append(('gfn1', ROOT + s, [], 900))
combos.append(('gfn1', perturbed('medium/caffeine.xyz'), [], 900))
for s in ['small/ethane.xyz', 'small/benzene.xyz', 'medium/benzoic-acid.xyz', 'medium/caffeine.xyz',
          'medium/methyl-a-d-galactopyranoside.xyz', 'clusters/water8-cluster.xyz']:
    combos.append(('gfn2', ROOT + s, [], 1200))
combos.append(('gfn2', perturbed('medium/caffeine.xyz'), [], 1200))
combos.append(('gfn2', ROOT + 'clusters/aaa-bglc.conf96.xyz', [], 2400))
for s in ['small/ethane.xyz', 'small/benzene.xyz', 'medium/benzoic-acid.xyz', 'medium/caffeine.xyz',
          'medium/methyl-a-d-galactopyranoside.xyz']:
    combos.append(('uff', ROOT + s, [], 900))
for s in ['small/ethane.xyz', 'small/benzene.xyz']:
    combos.append(('pm3', ROOT + s, [], 900))
# a large system: bounded to 300 steps, values 2 and 20 only
large = ('gfnff', ROOT + 'bulk/peo201-chain.xyz', ['-max_iterations', '300'], 3000)


def run(method, path, extra, timeout, v):
    tag = '%s_%s' % (method, os.path.basename(path).replace('.xyz', ''))
    d = os.path.join(OUT, tag, 'ls' + v)
    os.makedirs(d, exist_ok=True)
    shutil.copy(path, d + '/in.xyz')
    t = time.time()
    try:
        p = subprocess.run([BIN, '-opt', 'in.xyz', '-method', method, '-lbfgs_max_line_search', v, '-threads', '1',
                            '-no_bmt', '-verbosity', '1'] + extra,
                           cwd=d, capture_output=True, text=True, timeout=timeout)
        out = re.sub(r'\x1b\[[0-9;]*m', '', p.stdout + p.stderr)
        rc = p.returncode
    except subprocess.TimeoutExpired:
        out, rc = '', -9
    rows = re.findall(r'^\s+(\d+)\s+(-?\d+\.\d+)\s+(-?[\d.]+)\s+([\d.]+)\s+([\d.eE+-]+)', out, re.M)
    return dict(method=method, structure=os.path.basename(path), v=v, rc=rc,
                conv=(rc == 0 and 'did not converge' not in out),
                steps=int(rows[-1][0]) if rows else -1,
                E=float(rows[-1][1]) if rows else None, t=round(time.time() - t, 1))


jobs = [(m, p, e, to, v) for (m, p, e, to) in combos for v in VALUES]
jobs += [(large[0], large[1], large[2], large[3], v) for v in ('2', '20')]
results = []
with cf.ThreadPoolExecutor(5) as ex:
    futs = [ex.submit(run, *j) for j in jobs]
    for f in futs:
        results.append(f.result())
        json.dump(results, open(RES, 'w'))
print('done', len(results), 'runs')
