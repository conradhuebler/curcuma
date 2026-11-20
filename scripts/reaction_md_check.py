#!/usr/bin/env python3
"""
Bond analysis of a reactive-MD trajectory (multi-frame XYZ) - Claude Generated (Nov 2026).

Pure Python, no numpy. A bond i-j exists when r_ij < 1.25 * (rc_i + rc_j) (covalent radii, Cordero 2008, as used
in docs/REACTION_MD_EXAMPLES.md); fragments are the connected components; a fragment is named by its gross
formula (C, H first, then alphabetical). The 1.25 factor is a geometric criterion, not a chemical bond order:
hot, stretched H2 (above 0.775 A) and loosely bound intermediates are missed or split.

Usage:
  reaction_md_check.py <scenario> <trajectory.xyz | run directory> [--dt 1.0] [--dump 10]
  scenario: h2-metal | soup | wgs | haber | fe4 | doeb
  A directory is searched for the newest *.trj.xyz (BMT output directory or CWD).
Time of frame k is k * dump * dt (fs); frame 0 is the start structure.
Exit status 0 = analysis done (it does not judge whether the expected event happened; the report does).
"""
import argparse
import collections
import glob
import itertools
import math
import os
import sys

RC = {'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'Fe': 1.32, 'Pt': 1.36}  # Angstrom, Cordero 2008 (Fe low spin)
FACTOR = 1.25


def read_xyz(path):
    lines = open(path).read().split('\n')
    frames, i = [], 0
    while i < len(lines) and lines[i].strip():
        n = int(lines[i].split()[0])
        block = [l.split() for l in lines[i + 2:i + 2 + n]]
        frames.append(([b[0] for b in block], [tuple(float(v) for v in b[1:4]) for b in block]))
        i += n + 2
    return frames


def bonds(sym, xyz):
    out = set()
    for i, j in itertools.combinations(range(len(sym)), 2):
        if math.dist(xyz[i], xyz[j]) < FACTOR * (RC[sym[i]] + RC[sym[j]]):
            out.add((i, j))
    return out


def components(n, bset):
    adj = [[] for _ in range(n)]
    for i, j in bset:
        adj[i].append(j)
        adj[j].append(i)
    seen, comps = [False] * n, []
    for s in range(n):
        if seen[s]:
            continue
        seen[s] = True
        stack, comp = [s], []
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if not seen[v]:
                    seen[v] = True
                    stack.append(v)
        comps.append(sorted(comp))
    return comps


def formula(sym, idx):
    cnt = collections.Counter(sym[k] for k in idx)
    order = [e for e in ('C', 'H') if cnt[e]] + sorted(e for e in cnt if e not in ('C', 'H'))
    return ''.join(e + (str(cnt[e]) if cnt[e] > 1 else '') for e in order)


def species(sym, bset):
    return sorted(formula(sym, c) for c in components(len(sym), bset))


def fmt_time(fs):
    return f"{fs:.0f} fs" if fs < 1000 else f"{fs / 1000:.2f} ps"


def find_traj(path):
    if os.path.isdir(path):
        c = glob.glob(os.path.join(path, '**', '*.trj.xyz'), recursive=True)
        if not c:
            sys.exit(f"no *.trj.xyz below {path}")
        return max(c, key=os.path.getmtime)
    return path


def summary(sym, B, dt, final_species_n=1):
    last = collections.Counter(species(sym, B[-1]))
    print("final species :", dict(last))


def scen_h2_metal(sym, X, B, dt):
    # last two atoms are H2 (library order: metal atoms first, then H, H)
    n = len(sym)
    h = [n - 2, n - 1]
    metal = [k for k in range(n) if sym[k] == 'Pt']
    hh = [math.dist(x[h[0]], x[h[1]]) for x in X]
    first = next((k for k, d in enumerate(hh) if d > 2.0), None)
    print(f"H-H start/min/final/max : {hh[0]:.2f} / {min(hh):.2f} / {hh[-1]:.2f} / {max(hh):.2f} A")
    print("H-H > 2.0 A first at    :", fmt_time(first * dt) if first is not None else "never")
    nearest = [min(math.dist(X[-1][a], X[-1][m]) for m in metal) for a in h]
    print(f"final nearest Pt-H      : {nearest[0]:.2f}, {nearest[1]:.2f} A")
    ptpt = max(math.dist(x[a], x[b]) for x in X for a, b in itertools.combinations(metal, 2)
               if math.dist(X[0][a], X[0][b]) < 3.2)
    print(f"max Pt-Pt (start nbrs) : {ptpt:.2f} A (cluster intact if < ~3.3)")
    print("final species           :", dict(collections.Counter(species(sym, B[-1]))))
    print("EVENT dissociated       :", first is not None and all(d < 1.8 for d in nearest))


def scen_soup(sym, X, B, dt):
    b0 = B[0]
    ever = set().union(*B)
    broken = {b for b in b0 if any(b not in bs for bs in B)}
    new = {b for b in ever if b not in b0}
    print(f"start bonds {len(b0)}, broken at least once {len(broken)}, new atom pairs bonded at least once {len(new)}")
    newtypes = collections.Counter('-'.join(sorted((sym[i], sym[j]))) for i, j in new)
    print("new pairs by element    :", dict(newtypes))
    first_cc = next((k for k, bs in enumerate(B) if any(sym[i] == 'C' and sym[j] == 'C' for i, j in bs)), None)
    first_new_cx = next((k for k, bs in enumerate(B)
                         if any(b in new and 'C' in (sym[b[0]], sym[b[1]]) and 'H' not in (sym[b[0]], sym[b[1]]) for b in bs)), None)
    print("first C-C bond          :", fmt_time(first_cc * dt) if first_cc is not None else "never")
    print("first new C-X (X!=H)    :", fmt_time(first_new_cx * dt) if first_new_cx is not None else "never")
    seen = collections.defaultdict(lambda: None)
    for k, bs in enumerate(B):
        for s in species(sym, bs):
            if seen[s] is None:
                seen[s] = k * dt
    interest = ['C2H2', 'C2H4', 'C2H6', 'CH3', 'CH2', 'H2N', 'HO', 'H', 'CH2O', 'CHO', 'CH3O', 'CH4O', 'C2H4O', 'H2', 'NH', 'CN', 'CHN', 'CH3N', 'CO', 'CO2']
    print("first seen (selected)   :", {s: fmt_time(seen[s]) for s in interest if seen[s] is not None})
    print("all species ever (n)    :", len(seen))
    print("final species           :", dict(collections.Counter(species(sym, B[-1]))))
    rmax = max(math.dist(x[a], (0, 0, 0)) for x in X for a in range(len(sym)))
    print(f"max distance from origin: {rmax:.2f} A (wall R 6.5)")
    print("EVENT bonds broken      :", len(broken) > 0, "| new C-C or C-O/C-N pair:",
          any({sym[i], sym[j]} <= {'C', 'O', 'N'} and ({sym[i], sym[j]} != {'N', 'O'}) for i, j in new))


def scen_wgs(sym, X, B, dt):
    first, inter = {}, {}
    cnt_co2 = 0
    for k, bs in enumerate(B):
        sp = species(sym, bs)
        for s in sp:
            if s in ('CO2', 'H2'):
                first.setdefault(s, k * dt)
            if s in ('CHO2', 'CH2O2', 'CH2O3', 'CHO3', 'CH2O', 'CHO'):
                inter.setdefault(s, k * dt)
        cnt_co2 += 'CO2' in sp
    print("first CO2 / H2          :", {k: fmt_time(v) for k, v in first.items()} or "none")
    print("intermediates first     :", {k: fmt_time(v) for k, v in inter.items()} or "none")
    print(f"fraction of frames with CO2 : {cnt_co2 / len(B):.2f}")
    print("final species           :", dict(collections.Counter(species(sym, B[-1]))))
    print("EVENT CO2 formed        :", 'CO2' in first, "| H2 formed:", 'H2' in first)


def nh_stats(sym, bset):
    n_h = collections.Counter()
    for i, j in bset:
        for a, b in ((i, j), (j, i)):
            if sym[a] == 'N' and sym[b] == 'H':
                n_h[a] += 1
    return n_h


def scen_haber(sym, X, B, dt, fe=False):
    nidx = [k for k, s in enumerate(sym) if s == 'N']
    nn = [math.dist(x[nidx[0]], x[nidx[1]]) for x in X]
    first_nh = next((k for k, bs in enumerate(B) if nh_stats(sym, bs)), None)
    # free NH3: component with exactly N1 H3 and no other element
    def nh3(bs):
        return sum(1 for c in components(len(sym), bs) if formula(sym, c) == 'H3N')
    n3 = [nh3(bs) for bs in B]
    first_nh3 = next((k for k, v in enumerate(n3) if v), None)
    tail = max(1, len(B) // 10)
    print(f"N-N start/min/max/final : {nn[0]:.2f} / {min(nn):.2f} / {max(nn):.2f} / {nn[-1]:.2f} A")
    print("first N-H bond          :", fmt_time(first_nh * dt) if first_nh is not None else "never")
    print("first free NH3          :", fmt_time(first_nh3 * dt) if first_nh3 is not None else "never")
    print(f"frames with >=1 free NH3 in last 10 % ({tail} frames): {sum(1 for v in n3[-tail:] if v)}/{tail}")
    print("free NH3 in last frame  :", n3[-1], "| max N-H bonds per N (any frame):",
          max((max(nh_stats(sym, bs).values(), default=0) for bs in B)))
    print("final species           :", dict(collections.Counter(species(sym, B[-1]))))
    print("final N-H bonds per N   :", dict(nh_stats(sym, B[-1])))
    if fe:
        fei = [k for k, s in enumerate(sym) if s == 'Fe']
        fe0 = [(a, b) for a, b in itertools.combinations(fei, 2)]
        print(f"max Fe-Fe over run      : {max(math.dist(x[a], x[b]) for x in X for a, b in fe0):.2f} A (start {max(math.dist(X[0][a], X[0][b]) for a, b in fe0):.2f})")
        partners = lambda n: [sym[j if i == n else i] for i, j in B[-1] if n in (i, j)]
        print("neighbours of N (last)  :", {n + 1: sorted(partners(n)) for n in nidx})
    print("EVENT NH3 (free) formed :", first_nh3 is not None, "| N-N max > 2.5 A:", max(nn) > 2.5)


def scen_doeb(sym, X, B, dt):
    """2 H2 + O2 with or without a Pt cluster (Doebereiner lighter), Claude Generated (Nov 2026).

    Events (geometric criterion only): H-H pairs of the start structure broken, O-O broken, first O-H bond,
    first O bound to two H (H2O unit, free or on Pt), first free H2O (component of formula H2O), first Pt-O/Pt-H contact.
    """
    h = [k for k, s in enumerate(sym) if s == 'H']
    o = [k for k, s in enumerate(sym) if s == 'O']
    pt = [k for k, s in enumerate(sym) if s == 'Pt']
    hh0 = [(a, b) for a, b in itertools.combinations([k for k, s in enumerate(sym) if s == 'H'], 2) if math.dist(X[0][a], X[0][b]) < 0.9]
    oo0 = [b for b in B[0] if sym[b[0]] == 'O' and sym[b[1]] == 'O']

    def oh_count(bs):
        c = collections.Counter()
        for i, j in bs:
            if {sym[i], sym[j]} == {'O', 'H'}:
                c[i if sym[i] == 'O' else j] += 1
        return c

    def first(pred):
        k = next((k for k, bs in enumerate(B) if pred(bs)), None)
        return fmt_time(k * dt) if k is not None else "never"

    def n_oh(bs):
        return sum(oh_count(bs).values())

    def h2o_unit(bs):
        return any(v >= 2 for v in oh_count(bs).values())

    def free_h2o(bs):
        return any(formula(sym, c) == 'H2O' for c in components(len(sym), bs))

    # H-H "intact" = distance < 1.0 A (a stretched, hot H2 stays an H2; the 1.25 criterion would call it broken)
    hh_pairs = [(a, b) for a, b in itertools.combinations(h, 2) if math.dist(X[0][a], X[0][b]) < 0.9]
    hh_n = [sum(1 for a, b in hh_pairs if math.dist(x[a], x[b]) < 1.0) for x in X]
    hh_id = {id(bs): k for k, bs in enumerate(B)}

    def hh_left(bs):
        return hh_n[hh_id[id(bs)]]

    def pt_contact(bs, el):
        return any({sym[i], sym[j]} == {'Pt', el} for i, j in bs)

    print(f"start H-H pairs {len(hh0)}, O-O pairs {len(oo0)}, Pt atoms {len(pt)}")
    print("first H-H pair broken   :", first(lambda bs: hh_left(bs) < len(hh0)))
    print("both H-H pairs broken   :", first(lambda bs: hh_left(bs) == 0))
    print("first O-O broken        :", first(lambda bs: oo0 and not all(b in bs for b in oo0)))
    print("first O-H bond          :", first(lambda bs: n_oh(bs) > 0))
    print("first O with 2 H (H2O)  :", first(h2o_unit))
    print("first free H2O          :", first(free_h2o))
    print("first Pt-H / Pt-O bond  :", first(lambda bs: pt_contact(bs, 'H')), "/", first(lambda bs: pt_contact(bs, 'O')))
    print("frames with an O-H bond :", f"{sum(1 for bs in B if n_oh(bs)) / len(B):.2f}",
          "| with an H2O unit:", f"{sum(1 for bs in B if h2o_unit(bs)) / len(B):.2f}")
    print("max O-H bonds (any frame)/final :", max(n_oh(bs) for bs in B), "/", n_oh(B[-1]))
    print("final H-H pairs (<1.0 A):", hh_left(B[-1]), "of", len(hh0), "| final O-O intact:",
          sum(1 for b in oo0 if b in B[-1]), "of", len(oo0))
    print("final species           :", dict(collections.Counter(species(sym, B[-1]))))
    if pt:
        pt0 = [(a, b) for a, b in itertools.combinations(pt, 2) if math.dist(X[0][a], X[0][b]) < 3.2]
        print(f"max Pt-Pt (start nbrs)  : {max(math.dist(x[a], x[b]) for x in X for a, b in pt0):.2f} A (start {max(math.dist(X[0][a], X[0][b]) for a, b in pt0):.2f})")
    rmax = max(math.dist(x[a], (0, 0, 0)) for x in X for a in range(len(sym)))
    print(f"max distance from origin: {rmax:.2f} A")
    print("EVENT H2O unit formed   :", any(h2o_unit(bs) for bs in B), "| O-H bond:", any(n_oh(bs) for bs in B))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('scenario', choices=['h2-metal', 'soup', 'wgs', 'haber', 'fe4', 'doeb'])
    ap.add_argument('path')
    ap.add_argument('--dt', type=float, default=1.0, help='time step in fs')
    ap.add_argument('--dump', type=int, default=10, help='dump_frequency in steps')
    a = ap.parse_args()
    path = find_traj(a.path)
    fr = read_xyz(path)
    sym = fr[0][0]
    X = [f[1] for f in fr]
    B = [bonds(sym, x) for x in X]
    dt = a.dt * a.dump
    print(f"trajectory              : {path}")
    print(f"frames                  : {len(fr)} ({fmt_time((len(fr) - 1) * dt)} of simulated time), atoms {len(sym)}")
    {'h2-metal': scen_h2_metal, 'soup': scen_soup, 'wgs': scen_wgs,
     'haber': scen_haber, 'fe4': lambda *p: scen_haber(*p, fe=True), 'doeb': scen_doeb}[a.scenario](sym, X, B, dt)


if __name__ == '__main__':
    main()
