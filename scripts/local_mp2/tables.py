#!/usr/bin/env python3
"""Print the markdown tables of docs/LOCAL_MP2_FEASIBILITY_2026-10.md from the json results in LMP2_WORK.

Reads every q*_*.json (tagged and untagged) written by the q1..q6 scripts and check_solver.py.
Usage: LMP2_WORK=<dir> ~/.venvs/qmref/bin/python -I scripts/local_mp2/tables.py [q1|q2|q3|q4|q5|q6 ...]
Claude Generated (Oct 2026).
"""
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lmp2_common as lc  # first: sets the thread environment before numpy loads

W = lc.work_dir()
ORDER = ["alkane4", "alkane8", "alkane12", "alkane16", "alkane20", "water4", "water8", "water12", "water16",
         "acene1", "acene2", "acene3", "acene4", "acene5"]


def load(prefix):
    rows = []
    for f in sorted(glob.glob(os.path.join(W, prefix + "*.json"))):
        rows += json.load(open(f))
    key = lambda r: ORDER.index(r["system"]) if r["system"] in ORDER else 99
    return sorted(rows, key=key)


def sci(x, d=1):
    return f"{x:.{d}e}"


def q1():
    R = load("q1_locality")
    print("\n### Q1a pair counts and pair-cutoff error (Pipek-Mezey, exact pair energies)\n")
    print("| system | atoms | o_val | pairs | n(|e_ij|>1e-5) | n(>1e-6) | n(>1e-7) | err R>8 | err R>10 | err R>12 |"
          " err R>15 | err R>15 /atom | slope log|e| vs log R (R>8) |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in R:
        if r["loc"] != "pm":
            continue
        c = r["cut_err"]
        print(f"| {lc.label(r['system'])} | {r['natm']} | {r['nval']} | {r['npairs']} | {r['nsig']['1e-05']} |"
              f" {r['nsig']['1e-06']} | {r['nsig']['1e-07']} | {sci(c['8'])} | {sci(c['10'])} | {sci(c['12'])} |"
              f" {sci(c['15'])} | {sci(c['15'] / r['natm'])} | {r.get('slope_logR', float('nan')):.2f} |")
    print("\n### Q1b localisation methods (mean/max spread in Bohr, BP(0.02) atoms per orbital, n(|e|>1e-6), err R>10)\n")
    print("| system | PM spread | Boys spread | IBO spread | PM BP | Boys BP | IBO BP | PM n>1e-6 | Boys | IBO |"
          " PM err R>10 | Boys | IBO |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    by = {}
    for r in R:
        by.setdefault(r["system"], {})[r["loc"]] = r
    for sname, d in by.items():
        if not all(k in d for k in ("pm", "boys", "ibo")):
            continue
        g = lambda k, f: f(d[k])
        print(f"| {lc.label(sname)} | " + " | ".join(f"{d[k]['mean_spread_bohr']:.2f}/{d[k]['max_spread_bohr']:.2f}"
                                                   for k in ("pm", "boys", "ibo")) + " | "
              + " | ".join(f"{d[k]['bp002_mean']:.1f}/{d[k]['bp002_max']}" for k in ("pm", "boys", "ibo")) + " | "
              + " | ".join(str(d[k]["nsig"]["1e-06"]) for k in ("pm", "boys", "ibo")) + " | "
              + " | ".join(sci(d[k]["cut_err"]["10"]) for k in ("pm", "boys", "ibo")) + " |")
    print("\n### Q1c decay: max |e_ij| per centroid-distance bin (Bohr), PM, largest system per class\n")
    for sname in ("alkane20", "water16", "acene5"):
        rr = [r for r in R if r["system"] == sname and r["loc"] == "pm"]
        if not rr:
            continue
        print(f"{lc.label(sname)}: " + ", ".join(f"[{a},{b}) n={n} max {m:.1e}" for a, b, n, m, s in rr[0]["decay"]))
        print()


def q2():
    R = load("q2_domains")
    cfgs = []
    for r in R:
        if r["cfg"] not in cfgs:
            cfgs.append(r["cfg"])
    print("\n### Q2a LMP2 domain error per atom (Eh/atom, E(LMP2) - E(RI-MP2)); recovery in % in brackets\n")
    print("| system | atoms | " + " | ".join(cfgs) + " |")
    print("|---|---|" + "---|" * len(cfgs))
    by = {}
    for r in R:
        by.setdefault(r["system"], {})[r["cfg"]] = r
    for sname, d in by.items():
        print(f"| {lc.label(sname)} | {list(d.values())[0]['natm']} | " + " | ".join(
            (f"{sci(d[c]['err_per_atom'])} ({d[c]['recovery_pct']:.2f})" if c in d else "-") for c in cfgs) + " |")
    print("\n### Q2b domain sizes: atoms per orbital domain [i] (mean/max), pair space dimension m (mean) / nvir,"
          " iterations\n")
    print("| system | " + " | ".join(cfgs) + " |")
    print("|---|" + "---|" * len(cfgs))
    for sname, d in by.items():
        print(f"| {lc.label(sname)} | " + " | ".join(
            (f"{d[c]['dom_atoms_mean']:.1f}/{d[c]['dom_atoms_max']}, m {d[c]['m_mean']:.0f}/{d[c]['nvir']}, "
             f"it {d[c]['iterations']}" if c in d else "-") for c in cfgs) + " |")


def q3(base):
    R = [r for r in load("q3_pno") if r["base"] == base]
    if not R:
        return
    tcs = sorted({r["tcut"] for r in R}, reverse=True)
    by = {}
    for r in R:
        by.setdefault(r["system"], {})[r["tcut"]] = r
    print(f"\n### Q3a ({base}) error per atom with the semicanonical PNO correction (Eh/atom); without it in brackets\n")
    print("| system | atoms | " + " | ".join(f"T_CutPNO {t:.0e}" for t in tcs) + " |")
    print("|---|---|" + "---|" * len(tcs))
    for sname, d in by.items():
        print(f"| {lc.label(sname)} | {list(d.values())[0]['natm']} | " + " | ".join(
            (f"{sci(d[t]['err_corr_per_atom'])} ({sci(d[t]['err'] / d[t]['natm'])})" if t in d else "-") for t in tcs)
              + " |")
    print(f"\n### Q3b ({base}) PNOs per pair: mean over all pairs / mean over pairs with centroid distance < 10 Bohr"
          f" / max; iterations\n")
    print("| system | pairs | nvir | " + " | ".join(f"{t:.0e}" for t in tcs) + " |")
    print("|---|---|---|" + "---|" * len(tcs))
    for sname, d in by.items():
        r0 = list(d.values())[0]
        print(f"| {lc.label(sname)} | {r0['npairs']} | {r0['nvir']} | " + " | ".join(
            (f"{d[t]['pno_mean']:.1f} / {d[t]['pno_mean_close']:.1f} / {d[t]['pno_max']}, it {d[t]['iterations']}"
             if t in d else "-") for t in tcs) + " |")
    print(f"\n### Q3c ({base}) flop proxies relative to the canonical RI-MP2 energy step o^2 v^2 naux"
          f" (K_PNO / B->PNO transform / one iteration all k / one iteration |F_ik|>1e-5 / overlaps S^pq all k)\n")
    print("| system | canonical flops | " + " | ".join(f"{t:.0e}" for t in tcs) + " |")
    print("|---|---|" + "---|" * len(tcs))
    for sname, d in by.items():
        r0 = list(d.values())[0]
        can = r0["cost"]["canonical"]
        cell = lambda c: (f"{c['k_pno'] / can:.2g} / {c['b_transform'] / can:.2g} / {c['iter_all_k'] / can:.2g} /"
                          f" {c['iter_screened'] / can:.2g} / {c['overlaps_all_k'] / can:.2g}")
        print(f"| {lc.label(sname)} | {can:.2e} | " + " | ".join((cell(d[t]["cost"]) if t in d else "-") for t in tcs)
              + " |")


def q4():
    R = load("q4_metric")
    print("\n### Q4a inverse metric: max |element| per atom-distance bin relative to the global max"
          " (J^-1 / Coulomb-normalised Jn^-1 / J itself)\n")
    for r in R:
        print(f"{lc.label(r['system'])} (naux {r['naux']}, cond(J) {r['cond_J']:.1e}): " + ", ".join(
            f"[{a},{b}) {x:.1e}/{y:.1e}/{z:.1e}" for a, b, x, y, z in r["inv_decay"][1:]))
        print()
    radii = list(R[0]["rij"].keys())
    print("\n### Q4b RI-J Coulomb energy, local atom-centred fit domains, robust: error per atom (Eh/atom);"
          " plain (non-robust) error per atom in brackets\n")
    print("| system | atoms | naux | " + " | ".join(f"R={x}" for x in radii) + " |")
    print("|---|---|---|" + "---|" * len(radii))
    for r in R:
        print(f"| {lc.label(r['system'])} | {r['natm']} | {r['naux']} | " + " | ".join(
            f"{sci(r['rij'][x]['err_robust_per_atom'])} ({sci(r['rij'][x]['err_plain'] / r['natm'])})" for x in radii)
              + " |")
    print("\n### Q4c RI-MP2 with LMO-centred local fit domains (robust): error per atom (Eh/atom); mean domain size"
          " (aux functions) in brackets\n")
    print("| system | atoms | naux | " + " | ".join(f"R={x}" for x in radii) + " |")
    print("|---|---|---|" + "---|" * len(radii))
    for r in R:
        print(f"| {lc.label(r['system'])} | {r['natm']} | {r['naux']} | " + " | ".join(
            f"{sci(r['rimp2'][x]['err_per_atom'])} ({r['rimp2'][x]['ndom_mean']:.0f})" for x in radii) + " |")
    T = load("q4_timing")
    if T:
        print("\n### Q4d Cholesky timing (ms, min of 3, n=1): global vs all atom-centred local blocks\n")
        print("| system | atoms | naux | global | R=6 local (ndom) | R=8 local (ndom) | R=10 local (ndom) |"
              " flop ratio local/global R=8 |")
        print("|---|---|---|---|---|---|---|---|")
        for r in T:
            L = r["local"]
            print(f"| {r['system']} | {r['natm']} | {r['naux']} | {r['t_global_s'] * 1e3:.1f} | " + " | ".join(
                f"{L[x]['t_local_s'] * 1e3:.1f} ({L[x]['ndom_mean']:.0f})" for x in ("6", "8", "10"))
                  + f" | {L['8']['flop_local'] / r['flop_global']:.2f} |")


def q5():
    R = load("q5_schwarz")
    if not R:
        return
    print("\n### Q5 iterations to convergence (|dE|<1e-8, max|R|<1e-6) [iterations until within 1e-8 of the tight"
          " solution] and cost in matvec equivalents (outer + inner flops / flops of one full A T)\n")
    vars_ = []
    for r in R:
        if r["variant"] not in vars_:
            vars_.append(r["variant"])
    print("| system | pairs | PNO mean | " + " | ".join(vars_) + " |")
    print("|---|---|---|" + "---|" * len(vars_))
    by = {}
    for r in R:
        by.setdefault((r["system"], r["tcut"]), {})[r["variant"]] = r
    for (sname, tc), d in by.items():
        r0 = list(d.values())[0]
        print(f"| {lc.label(sname)} (T {tc:.0e}) | {r0['npairs']} | {r0['pno_mean']:.1f} | " + " | ".join(
            (f"{d[v]['iterations']} [{d[v]['it_to_tight_1e-8']}], {(d[v]['flop_outer'] + d[v]['flop_inner']) / d[v]['flop_per_matvec']:.0f}"
             + (f", sub {d[v]['sub_mean']:.0f}" if "sub_mean" in d[v] else "") if v in d else "-") for v in vars_) + " |")


def q6():
    R = load("q6_orca")
    if not R:
        return
    print("\n### Q6 ORCA 6.1 DLPNO-MP2 vs ORCA RI-MP2 (def2-SVP/C), frozen core\n")
    print("| system | E_corr RI-MP2 | NormalPNO error (recovery) | PNOs/pair | TightPNO error (recovery) | PNOs/pair |")
    print("|---|---|---|---|---|---|")
    for r in R:
        print(f"| {lc.label(r['system'])} | {r['rimp2']['e_corr']:.8f} | {sci(r['err_normal'], 2)}"
              f" ({r['recovery_normal_pct']:.3f}%) | {r['dlpno_normal'].get('pno_mean')} |"
              f" {sci(r['err_tight'], 2)} ({r['recovery_tight_pct']:.3f}%) | {r['dlpno_tight'].get('pno_mean')} |")


def main():
    what = sys.argv[1:] or ["q1", "q2", "q3", "q4", "q5", "q6"]
    for w in what:
        if w == "q3":
            for base in sorted({r["base"] for r in load("q3_pno")}):
                q3(base)
        else:
            globals()[w]()


if __name__ == "__main__":
    main()
