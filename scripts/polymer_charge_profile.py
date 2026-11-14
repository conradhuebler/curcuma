#!/usr/bin/env python3
"""Claude Generated (Oct 2026): charge profile around the protonated nitrogens of a PDMAEMA chain, by bond distance.

Usage:  polymer_charge_profile.py chain.xyz label=q_file [label=q_file ...]
        q files are written by  curcuma -sp chain.xyz -method M -charge N -dump_atomic_charges FILE
Prints (1) the mean charge per atom type around the N+ (N with four neighbours), (2) the cumulative charge within d bonds of
the N+ (mean over all centres). Bonds: distance < 1.2 x sum of covalent radii. Used for docs/POLYMERBUILD_CHARGE.md.
"""
import sys
import numpy as np

radius = {'H': .31, 'C': .76, 'N': .71, 'O': .66}
xyz, pairs = sys.argv[1], [a.split('=', 1) for a in sys.argv[2:]]
lines = open(xyz).read().split('\n')[2:]
atoms = [(l.split()[0], np.array(l.split()[1:4], float)) for l in lines if l.strip()]
el = [a[0] for a in atoms]
X = np.array([a[1] for a in atoms])
n = len(el)
D = np.linalg.norm(X[:, None] - X[None], axis=2)
adj = {i: [j for j in range(n) if j != i and D[i, j] < 1.2 * (radius[el[i]] + radius[el[j]])] for i in range(n)}
centres = [i for i in range(n) if el[i] == 'N' and len(adj[i]) == 4]
Q = {label: np.loadtxt(path) for label, path in pairs}

def bond_distance(start):
    d, front = {start: 0}, [start]
    while front:
        nxt = []
        for a in front:
            for b in adj[a]:
                if b not in d:
                    d[b] = d[a] + 1
                    nxt.append(b)
        front = nxt
    return d

names = ['N+', 'H on N', 'C bonded to N', 'H on those C', 'next C (beta)', 'H on beta C', 'O (next)']
classes = {k: [] for k in names}
for c in centres:
    for a, v in bond_distance(c).items():
        if v == 0: classes['N+'].append(a)
        elif v == 1 and el[a] == 'H': classes['H on N'].append(a)
        elif v == 1: classes['C bonded to N'].append(a)
        elif v == 2 and el[a] == 'H': classes['H on those C'].append(a)
        elif v == 2 and el[a] == 'C': classes['next C (beta)'].append(a)
        elif v == 3 and el[a] == 'H': classes['H on beta C'].append(a)
        elif v == 3 and el[a] == 'O': classes['O (next)'].append(a)
print('%d centres. Mean charge per position (e):' % len(centres))
print('%-16s %4s' % ('position', 'n') + ''.join('%12s' % k for k in Q))
for k in names:
    idx = classes[k]
    print('%-16s %4d' % (k, len(idx) // max(1, len(centres))) + ''.join('%12.3f' % q[idx].mean() for q in Q.values()))
print('\nCumulative charge within d bonds of N+ (mean over centres, e):')
print('%-26s' % 'model' + ''.join('%8s' % ('d<=%d' % k) for k in range(9)))
for label, q in Q.items():
    row = [np.mean([sum(q[a] for a, v in bond_distance(c).items() if v <= k) for c in centres]) for k in range(9)]
    print('%-26s' % label + ''.join('%8.2f' % x for x in row))
print('atoms within d bonds       ' + ' '.join('%7d' % int(np.mean([sum(1 for v in bond_distance(c).values() if v <= k) for c in centres])) for k in range(9)))
