#!/usr/bin/env python3
"""Claude Generated (Oct 2026): GFN-FF -gfnff.cn_zero_guard (Coulomb TERM 1b chain rule at CN = 0).

Uses only the curcuma binary; every structure is built at run time from the library water (structlib id "water"),
no structure file is added to the repository.
  1. The switch never changes the energy: printed single-point energy and the energy line of -dump_gradient are equal
     with and without it, for an isolated Cl- ... H2O (H-Cl 3.2 A), a bound Cl- ... H2O (H-Cl 2.2 A) and a water dimer.
  2. Isolated Cl- (H-Cl 3.2 A, CN(Cl) = 0 exactly): with the switch the analytic gradient agrees with a central finite
     difference of the own energy (h = 2e-3 A, 12 components incl. the Cl ones) within 2e-5 Eh/A; so does the default
     without the switch (charged, two fragments: frag_charge_model ensemble is a fragment-charge override, and every
     override applies the guard since Oct 31, 2026); with -gfnff.frag_charge_model reference and the Cl- as atom 1 the
     reference formula (switch off) is off by more than 1e-4 Eh/A, the switch repairs it.
  3. Bound Cl- (H-Cl 2.2 A, CN > 0) and the neutral water dimer: gradients bit-identical with and without the switch,
     and within 2e-5 Eh/A of the finite difference.
Usage: check_cn_zero_guard.py CURCUMA_BINARY SOURCE_DIR
"""
import os, re, subprocess, sys, tempfile
import numpy as np

curcuma, src = sys.argv[1], sys.argv[2]
failed = 0
E_RE = re.compile(r"Single Point Energy = (-?\d+\.\d+)")
ON = ["-gfnff.cn_zero_guard", "true"]


def check(ok, what):
    global failed
    print(("  PASS  " if ok else "  FAIL  ") + what)
    failed += 0 if ok else 1


work = tempfile.mkdtemp(prefix="curcuma_cn_zero_guard_")
os.chdir(work)
count = [0]


def sp(el, X, charge, extra=(), grad=False):
    """single point in a fresh directory; returns (printed energy string, gradient Eh/A or None, dump energy line)"""
    count[0] += 1
    d = os.path.join(work, f"r{count[0]}")
    os.mkdir(d)
    with open(os.path.join(d, "in.xyz"), "w") as f:
        f.write(f"{len(el)}\nx\n" + "".join("%s %.10f %.10f %.10f\n" % (e, *x) for e, x in zip(el, X)))
    cmd = [curcuma, "-sp", "in.xyz", "-method", "gfnff", "-charge", str(charge), "-verbosity", "0", "-no_bmt",
           "-threads", "1", "-gfnff.cache_topology", "false"] + list(extra)
    if grad:
        cmd += ["-dump_gradient", "g.txt"]
    p = subprocess.run(cmd, cwd=d, capture_output=True, text=True)
    m = E_RE.findall(p.stdout + p.stderr)
    G, eline = None, None
    if grad and os.path.exists(os.path.join(d, "g.txt")):
        txt = open(os.path.join(d, "g.txt")).read()
        eline = [l for l in txt.split("\n") if l.startswith("# energy")][0].split(",")[0]
        G = np.loadtxt(txt.split("\n"), comments="#").reshape(len(el), 3)
    return (m[-1] if m else None), G, eline


def fd_error(el, X, charge, extra, h=0.002):
    _, G, _ = sp(el, X, charge, extra, grad=True)
    err = []
    for i in range(len(el)):
        for k in range(3):
            Y = X.copy(); Y[i, k] += h; ep = float(sp(el, Y, charge, extra)[0])
            Y = X.copy(); Y[i, k] -= h; em = float(sp(el, Y, charge, extra)[0])
            err.append(abs((ep - em) / (2 * h) - G[i, k]))
    return max(err)


def library_water():
    path = subprocess.run([sys.executable, src + "/scripts/structlib.py", "path", "water"], capture_output=True,
                          text=True).stdout.strip()
    L = open(path).read().split("\n")
    el = [l.split()[0] for l in L[2:5]]
    X = np.array([[float(v) for v in l.split()[1:4]] for l in L[2:5]])
    o = el.index("O")
    order = [o] + [i for i in range(3) if i != o]          # O first, then the two H
    return [el[i] for i in order], X[order] - X[o]


W_EL, W_X = library_water()


def water(c):
    return list(W_EL), W_X + c


def cl_water(d_hcl, seed):
    el, X = water(np.zeros(3))
    v = X[1] - X[0]
    X = np.vstack([X, X[1] + d_hcl * v / np.linalg.norm(v)]) + np.random.RandomState(seed).normal(0, 0.06, (4, 3))
    return el + ["Cl"], X


el_iso, X_iso = cl_water(3.2, 11)   # CN(Cl) = 0 exactly
el_cn, X_cn = cl_water(2.2, 12)     # CN(Cl) > 0
e1, x1 = water(np.zeros(3)); e2, x2 = water(np.array([2.9, 0.0, 0.3]))
el_w, X_w = e1 + e2, np.vstack([x1, x2]) + np.random.RandomState(5).normal(0, 0.06, (6, 3))
SYS = {"Cl- ... H2O 3.2 A": (el_iso, X_iso, -1), "Cl- ... H2O 2.2 A": (el_cn, X_cn, -1), "water dimer": (el_w, X_w, 0)}

print("1. energy unchanged by the switch")
G_off, G_on = {}, {}
for name, (el, X, q) in SYS.items():
    e_off = sp(el, X, q)[0]
    e_on = sp(el, X, q, ON)[0]
    _, G_off[name], l_off = sp(el, X, q, grad=True)
    _, G_on[name], l_on = sp(el, X, q, ON, grad=True)
    check(e_off is not None and e_off == e_on and l_off == l_on,
          f"{name}: printed E {e_off} / {e_on}, dump '{l_off}' / '{l_on}'")

print("2. isolated Cl- (CN = 0): gradient against central finite difference, 12 components")
err_on = fd_error(el_iso, X_iso, -1, ON)
err_off = fd_error(el_iso, X_iso, -1, [])
check(err_on <= 2e-5, f"switch on: max |FD - analytic| = {err_on:.2e} Eh/A (<= 2e-5)")
# Claude Generated (Oct 31, 2026): a charged system with several fragments runs the default frag_charge_model ensemble,
# whose variants are fragment-charge overrides, and every override applies the guard (GFNFF::cnZeroGuardEffective, the
# charge-groups Cl- defect of 1 Eh/A): the switch changes nothing here (before: 6.58e-4 Eh/A with the switch off).
check(err_off <= 2e-5, f"switch off, default charge model (ensemble override, guard always on): max |FD - analytic| = {err_off:.2e} Eh/A (<= 2e-5)")
check(np.array_equal(G_off["Cl- ... H2O 3.2 A"], G_on["Cl- ... H2O 3.2 A"]), "default charge model: the switch does not change this gradient")
# the reference charge model puts the net charge on the fragment of atom 1: Cl- first, so that it carries the -1; only
# there the reference formula (switch off) is still in force
REF = ["-gfnff.frag_charge_model", "reference"]
el_ref, X_ref = [el_iso[3]] + el_iso[:3], np.vstack([X_iso[3:4], X_iso[:3]])
err_ref_off = fd_error(el_ref, X_ref, -1, REF)
err_ref_on = fd_error(el_ref, X_ref, -1, REF + ON)
check(err_ref_off > 1e-4, f"reference charge model (Cl- first), switch off (reference formula): max |FD - analytic| = {err_ref_off:.2e} Eh/A (> 1e-4)")
check(err_ref_on <= 2e-5, f"reference charge model (Cl- first), switch on: max |FD - analytic| = {err_ref_on:.2e} Eh/A (<= 2e-5)")

print("3. CN > 0 and neutral: gradient unchanged and consistent")
for name in ("Cl- ... H2O 2.2 A", "water dimer"):
    el, X, q = SYS[name]
    check(np.array_equal(G_off[name], G_on[name]), f"{name}: gradient bit-identical with and without the switch")
    err = fd_error(el, X, q, ON)
    check(err <= 2e-5, f"{name}: switch on, max |FD - analytic| = {err:.2e} Eh/A (<= 2e-5)")

print(f"{count[0]} curcuma runs, {failed} failure(s)")
sys.exit(1 if failed else 0)
