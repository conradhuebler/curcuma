#!/usr/bin/env python3
"""Print the markdown tables of docs/SCHWARZ_LEARNED_OPERATOR_2026-10.md from the json results in LMP2_WORK
(free_<sys>.json, transfer_<family>.json, lrsweep_<sys>.json, rank_needed.json).

Usage: LMP2_WORK=<dir> ~/.venvs/qmref/bin/python -I scripts/schwarz_learn/tables.py
Claude Generated (Oct 2026).
"""
import json
import os
from collections import defaultdict

import numpy as np

W = os.environ.get("LMP2_WORK", "local_mp2_work")


def load(f):
    p = os.path.join(W, f)
    return json.load(open(p)) if os.path.exists(p) else None


def e(x, fmt="{:.1e}"):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "-"
    return fmt.format(x)


def it(x):
    return "n.c." if x in (400, 2000) else ("indef." if x == -1 else str(x))


def seed_group(rows, keyf):
    g = defaultdict(list)
    for r in rows:
        g[keyf(r)].append(r)
    return g


def med(rows, k):
    v = [r[k] for r in rows if k in r]
    return float(np.median(v)) if v else None


def span(rows, k):
    v = [r[k] for r in rows if k in r]
    return f"{e(min(v))}..{e(max(v))}" if v else "-"


def free_tables():
    for name in ("alkane8", "water4", "acene2"):
        o = load(f"free_{name}.json")
        if not o:
            continue
        nd = o["naux"]
        pd = nd * (nd + 1) // 2
        print(f"\n### {name}: natm {o['natm']}, naux {nd}, E_J {o['E_J']:.3f} Eh, E_corr {o['E_mp2']:.6f} Eh, "
              f"{o['steps']} Adam steps, lr {o['lr']}, Cholesky of J {o['t_cholesky_s'] * 1e3:.1f} ms; local density "
              f"fit R=8 (outside the operator class): dE_J {o['local_density_fit_R8_EJ_err']:.1e} Eh\n")
        print("| operator | params / dense | flops / dense | val. rel. J-norm err | Frob. err (norm.) | dE_J robust "
              "(Eh/atom) | dE_MP2 robust (Eh/atom) | dE_MP2 plain (Eh/atom) | PCG it. | GMRES it. | train s |")
        print("|---|---|---|---|---|---|---|---|---|---|---|")
        rows = o["rows"]
        na = o["natm"]

        def line(lbl, rs):
            r0 = rs[0]
            val = med(rs, "val_rel_jnorm")
            vtxt = e(val) if len(rs) == 1 else f"{e(val)} ({span(rs, 'val_rel_jnorm')})"
            ej = med(rs, "EJ_robust_err")
            mp = med(rs, "Emp2_robust_err")
            mpp = med(rs, "Emp2_plain_err")
            tt = med(rs, "train_time_s")
            print(f"| {lbl} | {r0['params'] / pd:.3f} | {r0['flops_apply'] / (2 * nd * nd):.3f} | {vtxt} | "
                  f"{e(med(rs, 'frob_rel_norm'))} | {e(ej / na if ej is not None else None)} | "
                  f"{e(mp / na if mp is not None else None)} | {e(mpp / na if mpp is not None else None)} | "
                  f"{it(r0.get('pcg_iter'))} | {it(r0.get('gmres_iter'))} | {e(tt, '{:.0f}') if tt else '-'} |")

        for r in rows:
            if r["kind"] in ("exact", "reference"):
                line(r["label"], [r])
        for r in rows:
            if r["kind"] == "global-closed":
                line(r["label"], [r])
        for r in rows:
            if r["kind"] == "schwarz-exact" and (r["variant"], r["coarse"]) in (("asm", None), ("asm-pu", None),
                                                                                 ("asm-pu", "bal"), ("ras-core", "bal")):
                line(r["label"], [r])
        for r in rows:
            if r["kind"] == "schwarz-rrr":
                line(r["label"], [r])
        g = seed_group([r for r in rows if r["kind"] == "schwarz-learned"], lambda r: (r["R"], r["r_l"], r["coarse"]))
        for (R, rl, co), rs in g.items():
            line(f"learned Schwarz R={R} r_l={rl} coarse={co} ({len(rs)} seeds)", rs)
        g = seed_group([r for r in rows if r["kind"] == "global-learned"], lambda r: (r["rank"], r["match"]))
        for (k, mt), rs in g.items():
            line(f"learned global rank {k} ({mt}, {len(rs)} seeds)", rs)


def transfer_tables():
    for fam in ("alkanes", "acenes"):
        o = load(f"transfer_{fam}.json")
        if not o:
            continue
        print(f"\n### Transfer {fam}: train {', '.join(o['train'])}; test {', '.join(o['test'])}; {o['steps']} steps\n")
        for k, c in o["classes"].items():
            print(f"- R, r_l = {k}: {c['n_classes']} classes for {c['n_subdomains']} training subdomains, "
                  f"max |Jn_ii - Jn_jj| within a class {c['max_dev_equal_keys']:.1e}")
        print("\n| system | role | operator | params | val. rel. J-norm err | dE_J robust (Eh/atom) | "
              "dE_MP2 robust (Eh/atom) | GMRES it. | unseen subdomains |")
        print("|---|---|---|---|---|---|---|---|---|")
        rows = o["rows"]
        names = o["train"] + o["test"]
        for n in names:
            rs = [r for r in rows if r["system"] == n]
            if not rs:
                continue
            na = rs[0]["natm"]
            for r in rs:
                if r["kind"] in ("reference", "schwarz-exact", "global-oracle"):
                    if r["kind"] == "schwarz-exact" and r["R"] != 6.0:
                        continue
                    if r["kind"] == "global-oracle" and (r["r_l"] != 32 or r["R"] != 6.0):
                        continue
                    mp = r.get("Emp2_robust_err")
                    print(f"| {n} | {r['role']} | {r['label']} | {r['params']} | {e(r.get('val_rel_jnorm'))} | "
                          f"{e(r['EJ_robust_err'] / na)} | {e(mp / na if mp is not None else None)} | "
                          f"{it(r['gmres_iter'])} | - |")
            g = seed_group([r for r in rs if r["kind"] == "schwarz-learned-shared"], lambda r: (r["R"], r["r_l"]))
            for (R, rl), gs in g.items():
                mp = med(gs, "Emp2_robust_err")
                print(f"| {n} | {gs[0]['role']} | learned shared R={R} r_l={rl} + Galerkin ({len(gs)} seeds) | "
                      f"{gs[0]['params']} | {e(med(gs, 'val_rel_jnorm'))} ({span(gs, 'val_rel_jnorm')}) | "
                      f"{e(med(gs, 'EJ_robust_err') / na)} | {e(mp / na if mp is not None else None)} | "
                      f"{it(gs[0]['gmres_iter'])} | {gs[0]['unseen']} (dev {gs[0]['max_dev_vs_train_blocks']:.0e}) |")


def lrsweep_table():
    o = load("lrsweep_alkane4.json")
    if not o:
        return
    print(f"\n### Learning-rate sweep, {o['system']} (naux {o['naux']}), {o['steps']} steps, 3 seeds; "
          f"Schwarz params {o['schwarz_params']}, ranks {o['ranks']}\n")
    for tag, v in o["optimum"].items():
        print(f"- closed-form optimum, global rank {v['rank']} ({tag}): train {v['train']:.2e}, val {v['val']:.2e}")
    print("\n| lr | model | params | final train loss median (min..max) | val loss median | Frob. err (norm.) | "
          "diverged |")
    print("|---|---|---|---|---|---|---|")
    g = seed_group(o["runs"], lambda r: (r["lr"], r["model"]))
    for (lr, m), rs in g.items():
        print(f"| {lr:.0e} | {m} | {rs[0]['params']} | {e(med(rs, 'train'))} ({span(rs, 'train')}) | "
              f"{e(med(rs, 'val'))} | {e(med(rs, 'frob_rel_norm'))} | {sum(r['diverged'] for r in rs)}/{len(rs)} |")


def rank_table():
    o = load("rank_needed.json")
    if not o:
        return
    print("\n### Rank needed by the best global low-rank operator (closed form)\n")
    print("| system | naux | k needed, dE_J < 1e-6 Eh/atom | k/naux | val. err (ia) at k=0.5 naux | at 0.9 naux | "
          "at 0.98 naux | dE_J/atom at 0.9 naux |")
    print("|---|---|---|---|---|---|---|---|")
    for r in o:
        c = {x["frac"]: x for x in r["curve"]}
        k = r["k_needed_EJ_1e-6_per_atom"]
        print(f"| {r['system']} | {r['naux']} | {k} | {k / r['naux']:.3f} | {e(c[0.5]['val_rel_jnorm'])} | "
              f"{e(c[0.9]['val_rel_jnorm'])} | {e(c[0.98]['val_rel_jnorm'])} | {e(c[0.9]['EJ_err_atom'])} |")


if __name__ == "__main__":
    rank_table()
    free_tables()
    transfer_tables()
    lrsweep_table()
