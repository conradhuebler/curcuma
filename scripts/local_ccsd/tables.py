#!/usr/bin/env python3
"""Print the markdown tables of docs/LOCAL_CCSD_PROTOTYPE_2026-10.md from the json results in LCCSD_WORK.

Usage: LCCSD_WORK=<dir> python3 scripts/local_ccsd/tables.py [stageA check accuracy schwarz orca]
Claude Generated (Oct 2026).
"""
import json
import os
import re
import sys

W = os.environ.get("LCCSD_WORK", os.path.join(os.getcwd(), "local_ccsd_work"))


def load(name):
    p = os.path.join(W, name)
    return json.load(open(p)) if os.path.exists(p) else []


def e(x, f="{:.1e}"):
    return "-" if x is None else f.format(x)


def stageA():
    print("| molecule | nbf | o / v | E_corr CCSD (Eh) | it (ORCA) | - PySCF | - ORCA | E(it 0) - MP2 PySCF | E(it 0) - MP2 curcuma "
          "| E_HF - curcuma | E_HF - ORCA | residual - PySCF update_amps | PM - canonical (it) |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in load("stageA_validation.json"):
        o, c = r["orca"] or {}, r["curcuma"] or {}
        print(f"| {r['system']} | {r['nbf']} | {r['o']} / {r['v']} | {r['e_ccsd']:.10f} | {r['iterations']} ({o.get('iterations', '-')}) "
              f"| {e(r['e_ccsd'] - r['pyscf_ccsd'])} | {e(r['e_ccsd'] - o['e_ccsd'] if o else None)} "
              f"| {e(r['e_iter0'] - r['pyscf_mp2'])} | {e(r['e_iter0'] - c['e_mp2'] if c else None)} "
              f"| {e(r['e_hf'] - c['e_hf'] if c else None)} | {e(r['e_hf'] - o['e_hf'] if o else None)} "
              f"| {e(r['d_update_amps'])} | {e(r['e_ccsd_pm'] - r['e_ccsd'])} ({r['iterations_pm']}) |")


def check():
    print("| system | full-space local CCSD - canonical | it | LMP2 - MP2 | it |")
    print("|---|---|---|---|---|")
    for r in load("check_full.json"):
        print(f"| {r['system']} | {e(r['d_ccsd'])} | {r['iterations']} | {e(r['d_mp2'])} | {r['iterations_mp2']} |")


def accuracy():
    A = [r for f in sorted(os.listdir(W)) if f.startswith("accuracy") and f.endswith(".json") for r in load(f)]
    order = ["alkane", "water", "acene"]
    A = list({(r["system"], r["space"]): r for r in A}.values())  # later files win (restarted runs)
    A.sort(key=lambda r: (order.index(r["system"].rstrip("0123456789")), int(re.sub(r"\D", "", r["system"]))))
    systems = list(dict.fromkeys(r["system"] for r in A))
    spaces = list(dict.fromkeys(r["space"] for r in A if r["space"] != "full"))
    hdr = lambda sp: sp.replace("|", " & ")  # a literal | would split the markdown cell
    get = {(r["system"], r["space"]): r for r in A}
    first = {}
    for r in A:
        first.setdefault(r["system"], r)
    print("CCSD error per atom (Eh/atom) [LMP2 error per atom, same spaces]:\n")
    print("| system | atoms | " + " | ".join(hdr(x) for x in spaces) + " |")
    print("|---|---|" + "---|" * len(spaces))
    for s in systems:
        row = []
        for sp in spaces:
            r = get.get((s, sp))
            row.append("-" if r is None else f"{r['err_ccsd'] / r['natm']:.1e} [{r['err_mp2'] / r['natm']:.1e}]")
        r0 = first[s]
        print(f"| {r0['label']} | {r0['natm']} | " + " | ".join(row) + " |")
    print("\nCCSD error per atom after correction: + (E_MP2 canonical - E_LMP2 same space) [+ dE_PNO semicanonical]:\n")
    print("| system | " + " | ".join(hdr(x) for x in spaces) + " |")
    print("|---|" + "---|" * len(spaces))
    for s in systems:
        row = []
        for sp in spaces:
            r = get.get((s, sp))
            if r is None:
                row.append("-")
            elif sp.startswith("bp") and "pno" not in sp:
                row.append(f"{r['err_ccsd_mp2corr'] / r['natm']:.1e}")
            else:
                row.append(f"{r['err_ccsd_mp2corr'] / r['natm']:.1e} [{r['err_ccsd_pnocorr'] / r['natm']:.1e}]")
        print(f"| {first[s]['label']} | " + " | ".join(row) + " |")
    print("\nPair-space size, kept t2 amplitudes and ladder flop proxy relative to canonical; iterations (local / full space"
          " in local orbitals, same criterion):\n")
    print("| system | o / v | " + " | ".join(hdr(x) for x in spaces) + " |")
    print("|---|---|" + "---|" * len(spaces))
    for s in systems:
        row = []
        for sp in spaces:
            r = get.get((s, sp))
            row.append("-" if r is None else
                       f"m {r['m_mean']:.0f}, {100 * r['frac_t2']:.1f} %, {r['flops']['ladder'] / r['flops']['canonical_ladder']:.3f}, it {r['iterations']}")
        rf = get.get((s, "full"))
        r0 = first[s]
        print(f"| {r0['label']} | {r0['o']} / {r0['v']} (it full {rf['iterations'] if rf else '-'}) | " + " | ".join(row) + " |")
    print("\nAbsolute numbers (kept t2 amplitudes, canonical o^2 v^2, ladder proxy flops):\n")
    print("| system | canonical t2 | canonical ladder | " + " | ".join(hdr(x) for x in spaces) + " |")
    print("|---|---|---|" + "---|" * len(spaces))
    for s in systems:
        r0 = first[s]
        row = [("-" if get.get((s, sp)) is None else
                f"{get[(s, sp)]['n_t2']:.2e} / {get[(s, sp)]['flops']['ladder']:.1e}") for sp in spaces]
        print(f"| {r0['label']} | {r0['n_t2_canonical']:.2e} | {r0['flops']['canonical_ladder']:.1e} | " + " | ".join(row) + " |")


def schwarz():
    R = []
    for fn in sorted(f for f in os.listdir(W) if f.startswith("schwarz_") and f.endswith(".json")):
        R += load(fn)
    get = {(r["system"], r["space"], r["variant"]): r for r in R}  # later files win
    R = list(get.values())
    variants = list(dict.fromkeys(v for v in ["jacobi+diis", "jacobi"] + [r["variant"] for r in R]))
    variants = [v for v in variants if any(r["variant"] == v for r in R)]
    order = ["alkane", "water", "acene"]
    systems = sorted(dict.fromkeys(r["system"] for r in R),
                     key=lambda x: (order.index(x.rstrip("0123456789")), int(re.sub(r"\D", "", x))))
    for space in dict.fromkeys(r["space"] for r in R):
        print(f"\nspace {space}: outer iterations until the energy stays within 1e-6 / 1e-8 Eh of the tight solution;"
              " Schwarz: + inner work in ladder units, mean subdomain size in pairs; last column fock_all / ladder:\n")
        print("| system | pairs | m | " + " | ".join(variants) + " | fock_all/ladder |")
        print("|---|---|---|" + "---|" * len(variants) + "---|")
        for sy in systems:
            rs = [r for r in R if r["system"] == sy and r["space"] == space]
            if not rs:
                continue
            row = []
            for v in variants:
                r = get.get((sy, space, v))
                if r is None:
                    row.append("-")
                    continue
                c = f"{r.get('it_to_tight_1e-6', '-')} / {r['it_to_tight_1e-8']}"
                if r["flop_inner"]:
                    c += f", {r['flop_inner'] / r['flops']['ladder']:.0f}, {r.get('sub_mean', 0):.0f}"
                row.append(c)
            r0 = rs[0]
            print(f"| {sy} | {r0['npairs']} | {r0['m_mean']:.1f} | " + " | ".join(row)
                  + f" | {r0['flops']['fock_all'] / r0['flops']['ladder']:.2f} |")


def orca():
    print("| system | atoms | E_corr CCSD (ORCA) | NormalPNO error (/atom) | PNOs/pair | TightPNO error (/atom) | PNOs/pair | RI-MP2(SVP/C) - MP2 "
          "| complete-PNO DLPNO - CCSD | Normal - complete (/atom) | Tight - complete (/atom) |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in load("orca_dlpno.json"):
        n = r["natm"]
        print(f"| {r['system']} | {n} | {r['ccsd']['e_corr']:.8f} | {r['err_normal']:.2e} ({r['err_normal'] / n:.1e}) "
              f"| {r['dlpno_normal'].get('pno_mean', '-')} | {r['err_tight']:.2e} ({r['err_tight'] / n:.1e}) "
              f"| {r['dlpno_tight'].get('pno_mean', '-')} | {r['aux_err_mp2']:.1e} | {e(r['complete_minus_ccsd'])} "
              f"| {e(r['err_normal_vs_complete'])} ({e(r['err_normal_vs_complete'] and r['err_normal_vs_complete'] / n)}) "
              f"| {e(r['err_tight_vs_complete'])} ({e(r['err_tight_vs_complete'] and r['err_tight_vs_complete'] / n)}) |")


if __name__ == "__main__":
    for k in (sys.argv[1:] or ["stageA", "check", "accuracy", "schwarz", "orca"]):
        print(f"\n## {k}\n")
        globals()[k]()
