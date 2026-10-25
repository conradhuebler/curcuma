#!/usr/bin/env python3
"""Markdown tables for docs/QM_TIMING_ORCA_PYSCF_2026-10.md from results.jsonl (Claude Generated, Oct 2026).

  report.py [results.jsonl]      (default ~/.cache/qmtime/results.jsonl)

Wall times in s (whole program run). Energy deviations in Eh relative to the reference column.
'cap' = killed at the wall-time cap, 'pred' = skipped because the previous size predicted more
than 1.5 x cap, '-' = not run. The last run of a job in the file wins.
"""
import json
import os
import sys

PATH = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/.cache/qmtime/results.jsonl")
R = {}
for ln in open(PATH):
    if ln.strip():
        r = json.loads(ln)
        R[(r["program"], r["method"], r["variant"], r["mol"], r["threads"])] = r

MOLS = {
    "hf": ["benzene", "anthracene", "alkane20", "alkane40", "alkane80"],
    "rimp2": ["benzene", "anthracene", "alkane10", "alkane12", "alkane16", "alkane20", "alkane40", "alkane80"],
    "ccsd": ["water4", "benzene", "water6", "alkane6", "water8", "alkane8", "water10", "water12"],
    "lmp2": ["alkane20", "alkane40", "alkane80"],
}
LABEL = {"alkane%d" % n: "C%dH%d" % (n, 2 * n + 2) for n in range(1, 200)}
LABEL.update({"water%d" % n: "(H2O)%d" % n for n in range(1, 50)})


def lab(m):
    return LABEL.get(m, m)


def nbf(method, mol):
    for k, r in R.items():
        if k[3] == mol and r.get("nbf"):
            return r["nbf"]
    return "-"


def cell(r):
    if r is None:
        return "-"
    if r.get("ok"):
        return f"{r['wall']:.1f}" if r["wall"] < 100 else f"{r['wall']:.0f}"
    k = r.get("killed") or ""
    if k == "timeout":
        return f"cap ({r['wall']:.0f})"
    if k.startswith("skipped"):
        return "pred"
    return k or f"fail rc={r.get('rc')}"


def get(p, m, v, mol, t):
    return R.get((p, m, v, mol, t))


def wall_table(method, cols, threads_list, mols=None):
    mols = mols or MOLS[method]
    hdr = ["system", "nbf"] + [f"{p} {v} {t}T" for (p, v) in cols for t in threads_list]
    out = ["| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)]
    for mol in mols:
        row = [lab(mol), str(nbf(method, mol))]
        any_run = False
        for (p, v) in cols:
            for t in threads_list:
                r = get(p, method if v != "dlpno" else "dlpno_ccsd", v if v != "dlpno" else "normalpno", mol, t)
                any_run |= r is not None
                row.append(cell(r))
        if any_run:
            out.append("| " + " | ".join(row) + " |")
    return "\n".join(out)


def energy_table(method, key, ref, cols, threads=24, mols=None, alt_threads=1):
    """Deviation of `key` (e_hf | e_corr) from the reference column, taken from the 24-thread run
    (or the 1-thread run when the 24-thread run is missing)."""
    mols = mols or MOLS[method]

    def val(p, v, mol):
        for t in (threads, alt_threads):
            r = get(p, method, v, mol, t)
            if r and r.get("ok") and r.get(key) is not None:
                return r[key]
        return None

    hdr = ["system", f"{ref[0]} {ref[1]} {key} (Eh)"] + [f"{p} {v} - ref" for (p, v) in cols]
    out = ["| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)]
    for mol in mols:
        e0 = val(ref[0], ref[1], mol)
        row = [lab(mol), f"{e0:.9f}" if e0 is not None else "-"]
        have = e0 is not None
        for (p, v) in cols:
            e = val(p, v, mol)
            row.append(f"{e - e0:+.1e}" if (e is not None and e0 is not None) else "-")
            have |= e is not None
        if have:
            out.append("| " + " | ".join(row) + " |")
    return "\n".join(out)


def ratio_table(method, pairs, threads_list, mols=None):
    """Wall-time ratio curcuma / other for each (curcuma col, other col) pair."""
    mols = mols or MOLS[method]
    hdr = ["system"] + [f"{a[1]}/{b[0]} {b[1]} {t}T" for (a, b) in pairs for t in threads_list]
    out = ["| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)]
    for mol in mols:
        row, have = [lab(mol)], False
        for (a, b) in pairs:
            for t in threads_list:
                ra, rb = get(a[0], method, a[1], mol, t), get(b[0], method, b[1], mol, t)
                if ra and rb and ra.get("ok") and rb.get("ok"):
                    row.append(f"{ra['wall'] / rb['wall']:.2f}")
                    have = True
                else:
                    row.append("-")
        if have:
            out.append("| " + " | ".join(row) + " |")
    return "\n".join(out)


def lmp2_error_table():
    """LMP2 / DLPNO-MP2 correlation energy minus canonical RI-MP2 (curcuma, same aux, 24T or gpu)."""
    hdr = ["system", "nbf", "RI-MP2 E_corr (Eh)", "curcuma lmp2", "ORCA NormalPNO", "ORCA TightPNO",
           "per atom: curcuma / Normal / Tight (uEh)"]
    out = ["| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)]
    for mol in MOLS["lmp2"]:
        ref = None
        for v, t in (("cpu", 24), ("gpu", 24), ("gpu_bigri", 24), ("cpu", 1)):
            r = get("curcuma", "rimp2", v, mol, t)
            if r and r.get("ok"):
                ref = r["e_corr"]
                break
        if ref is None:  # ORCA RI-MP2 as fallback reference
            r = get("orca", "rimp2", "nori", mol, 24)
            ref = r["e_corr"] if r and r.get("ok") else None
        nat = 3 * int(mol[6:]) + 2
        vals = []
        for p, v in (("curcuma", "cpu"), ("orca", "normalpno"), ("orca", "tightpno")):
            e = None
            for vv in ((v, "gpu") if p == "curcuma" else (v,)):
                for t in (24, 1):
                    r = get(p, "lmp2", vv, mol, t)
                    if r and r.get("ok"):
                        e = r["e_corr"]
                        break
                if e is not None:
                    break
            vals.append(e)
        row = [lab(mol), str(nbf("lmp2", mol)), f"{ref:.8f}" if ref is not None else "-"]
        row += [f"{e - ref:+.2e}" if (e is not None and ref is not None) else (f"{e:.8f} (abs)" if e is not None else "-")
                for e in vals]
        row.append(" / ".join(f"{(e - ref) / nat * 1e6:+.1f}" if (e is not None and ref is not None) else "-" for e in vals))
        out.append("| " + " | ".join(row) + " |")
    return "\n".join(out)


def loads():
    ls = [r.get("load_before") for r in R.values() if r.get("load_before") is not None]
    return f"load average (1 min) before each job: min {min(ls):.2f}, median {sorted(ls)[len(ls)//2]:.2f}, max {max(ls):.2f} over {len(ls)} jobs" if ls else ""


if __name__ == "__main__":
    print("## HF wall time (s)\n")
    hf_cols = [("curcuma", "cpu"), ("curcuma", "cpu_rijk"), ("orca", "nori"), ("orca", "rijk"), ("pyscf", "exact"), ("pyscf", "df")]
    print("### 24 threads\n")
    print(wall_table("hf", hf_cols[:2] + [("curcuma", "gpu_direct"), ("curcuma", "gpu_rijk")] + hf_cols[2:], [24]))
    print("\n### 24 threads, RI-JK with the stored (CPU) / streamed (GPU) tensor, -qm.ri_max_memory_mb 200000 (not default)\n")
    print(wall_table("hf", [("curcuma", "cpu_rijk_stored"), ("curcuma", "gpu_rijk_stored")], [24], mols=["alkane40", "alkane80"]))
    print("\n### 1 thread\n")
    print(wall_table("hf", hf_cols, [1]))
    print("\n### HF ratios curcuma / other\n")
    print(ratio_table("hf", [(("curcuma", "cpu"), ("orca", "nori")), (("curcuma", "cpu_rijk"), ("orca", "rijk")),
                             (("curcuma", "gpu_direct"), ("orca", "nori")), (("curcuma", "cpu"), ("pyscf", "exact")),
                             (("curcuma", "cpu_rijk"), ("pyscf", "df"))], [24]))
    print()
    print(ratio_table("hf", [(("curcuma", "cpu"), ("orca", "nori")), (("curcuma", "cpu_rijk"), ("orca", "rijk")),
                             (("curcuma", "cpu"), ("pyscf", "exact")), (("curcuma", "cpu_rijk"), ("pyscf", "df"))], [1]))
    print("\n### HF energies\n")
    print(energy_table("hf", "e_hf", ("pyscf", "exact"), [("curcuma", "cpu"), ("curcuma", "gpu_direct"), ("orca", "nori")]))
    print()
    print(energy_table("hf", "e_hf", ("pyscf", "df"), [("curcuma", "cpu_rijk"), ("curcuma", "gpu_rijk"), ("orca", "rijk")]))

    print("\n## RI-MP2 wall time (s)\n")
    mp_cols = [("curcuma", "cpu"), ("orca", "nori"), ("pyscf", "exact")]
    print(wall_table("rimp2", mp_cols[:1] + [("curcuma", "gpu")] + mp_cols[1:], [24]))
    print()
    print(wall_table("rimp2", mp_cols, [1]))
    print()
    print(ratio_table("rimp2", [(("curcuma", "cpu"), ("orca", "nori")), (("curcuma", "cpu"), ("pyscf", "exact")),
                                (("curcuma", "gpu"), ("orca", "nori"))], [24, 1]))
    print()
    print(energy_table("rimp2", "e_corr", ("pyscf", "exact"), [("curcuma", "cpu"), ("curcuma", "gpu"), ("orca", "nori")]))

    print("\n## Local MP2 wall time (s)\n")
    print(wall_table("lmp2", [("curcuma", "cpu"), ("curcuma", "gpu"), ("orca", "normalpno"), ("orca", "tightpno")], [24]))
    print()
    print(wall_table("lmp2", [("curcuma", "cpu"), ("orca", "normalpno"), ("orca", "tightpno")], [1]))
    print("\nCanonical RI-MP2 reference runs (s):\n")
    print(wall_table("rimp2", [("curcuma", "gpu"), ("curcuma", "gpu_bigri"), ("curcuma", "cpu"), ("orca", "nori")], [24], mols=MOLS["lmp2"]))
    print()
    print(lmp2_error_table())

    print("\n## CCSD wall time (s)\n")
    cc_cols = [("curcuma", "cpu"), ("orca", "nori"), ("pyscf", "exact")]
    print(wall_table("ccsd", cc_cols + [("orca", "dlpno")], [24]))
    print()
    print(wall_table("ccsd", cc_cols, [1]))
    print()
    print(ratio_table("ccsd", [(("curcuma", "cpu"), ("orca", "nori")), (("curcuma", "cpu"), ("pyscf", "exact"))], [24, 1]))
    print()
    print(energy_table("ccsd", "e_corr", ("pyscf", "exact"), [("curcuma", "cpu"), ("orca", "nori")]))
    print()
    print(loads())
