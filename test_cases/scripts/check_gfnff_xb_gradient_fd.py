#!/usr/bin/env python3
"""Claude Generated (Nov 2026): GFN-FF halogen-bond (XB) gradient against a central finite difference (Known Issue #46).

FFWorkspace::calcHalogenBonds had the derivative of the short-range damping with the wrong sign and that of the
long-range damping with r/longcut instead of 1/r (energy unaffected). This test checks the fully assembled surface:
-dump_gradient (Eh/A) against a central difference of the -dump_gradient energy line (15 decimals), every Cartesian
component, h = 1e-4 and 2e-4 A. Gate: the Richardson-extrapolated difference (4 D(h) - D(2h)) / 3 agrees with the
analytic gradient within 1e-8 Eh/A. The plain h = 1e-4 A difference is only reported: after the fix its error is the
O(h^2) truncation of the difference itself (ratio between the two h exactly 4; 8e-9 / 1.6e-8 Eh/A on the two
structures, Nov 2026), before the fix it was 1.2e-4 / 1.3e-4 Eh/A, independent of h.

Structures are built at run time from the library (no structure file is added):
  1. fi-pyridine (F-I ... N(pyridine), structlib id "fi-pyridine"), perturbed by N(0, 0.03 A) to break the planar symmetry
  2. fi-pyridine + acetic-acid-dimer (structlib id "acetic-acid-dimer") shifted by +8 A in x, both perturbed
The XB term must be active (nonzero XB energy at verbosity 2) in both, otherwise the test fails.
Usage: check_gfnff_xb_gradient_fd.py CURCUMA_BINARY SOURCE_DIR
"""
import os, re, subprocess, sys, tempfile
import numpy as np

curcuma, src = sys.argv[1], sys.argv[2]
failed = 0
GATE = 1e-8


def check(ok, what):
    global failed
    print(("  PASS  " if ok else "  FAIL  ") + what)
    failed += 0 if ok else 1


work = tempfile.mkdtemp(prefix="curcuma_xb_fd_")
count = [0]


def run(el, X, grad=True, verbosity="0"):
    """single point in a fresh directory (no topology cache); returns (energy, gradient Eh/A, stdout)"""
    count[0] += 1
    d = os.path.join(work, f"r{count[0]}")
    os.mkdir(d)
    with open(os.path.join(d, "in.xyz"), "w") as f:
        f.write(f"{len(el)}\nx\n" + "".join("%s %.12f %.12f %.12f\n" % (e, *x) for e, x in zip(el, X)))
    cmd = [curcuma, "-sp", "in.xyz", "-method", "gfnff", "-verbosity", verbosity, "-no_bmt", "-threads", "1",
           "-gfnff.cache_topology", "false", "-dump_gradient", "g.txt"]
    p = subprocess.run(cmd, cwd=d, capture_output=True, text=True)
    txt = open(os.path.join(d, "g.txt")).read()
    E = float(re.search(r"# energy (-?\d+\.\d+)", txt).group(1))
    G = np.loadtxt(txt.split("\n"), comments="#").reshape(len(el), 3) if grad else None
    return E, G, p.stdout + p.stderr


def fd_gradient(el, X, h):
    """central difference of the dump energy, every Cartesian component (Eh/A)"""
    D = np.zeros_like(X)
    for i in range(len(el)):
        for k in range(3):
            Y = X.copy(); Y[i, k] += h; ep = run(el, Y, grad=False)[0]
            Y = X.copy(); Y[i, k] -= h; em = run(el, Y, grad=False)[0]
            D[i, k] = (ep - em) / (2 * h)
    return D


def library(sid):
    path = subprocess.run([sys.executable, src + "/scripts/structlib.py", "path", sid], capture_output=True,
                          text=True).stdout.strip()
    L = open(path).read().split("\n")
    n = int(L[0])
    el = [l.split()[0] for l in L[2:2 + n]]
    X = np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])
    return el, X


fi_el, fi_X = library("fi-pyridine")
ac_el, ac_X = library("acetic-acid-dimer")
rng = np.random.RandomState(46)
SYS = {
    "fi-pyridine (perturbed)": (fi_el, fi_X + rng.normal(0, 0.03, fi_X.shape)),
    "fi-pyridine + acetic-acid-dimer": (fi_el + ac_el,
                                        np.vstack([fi_X, ac_X - ac_X.mean(0) + np.array([8.0, 0.0, 0.0])])
                                        + rng.normal(0, 0.03, (len(fi_el) + len(ac_el), 3))),
}

for name, (el, X) in SYS.items():
    print(f"{name}: {len(el)} atoms")
    out = run(el, X, verbosity="2")[2]
    m = re.search(r"X-bonds\s+([-+]?\d+\.\d+)", re.sub(r"\x1b\[[0-9;]*m", "", out))
    exb = float(m.group(1)) if m else 0.0
    check(m is not None and exb != 0.0, f"XB term active: '{m.group(0).strip() if m else 'no XB line'}'")
    E, G, _ = run(el, X)
    D1, D2 = fd_gradient(el, X, 1e-4), fd_gradient(el, X, 2e-4)
    e1, e2 = np.abs(D1 - G).max(), np.abs(D2 - G).max()
    # Richardson: (4 D(h) - D(2h)) / 3 removes the h^2 truncation term of the central difference
    eR = np.abs((4 * D1 - D2) / 3 - G).max()
    print(f"  max |FD - analytic|: h = 1e-4 A {e1:.3e}, h = 2e-4 A {e2:.3e} (ratio {e2 / max(e1, 1e-300):.2f}),"
          f" Richardson {eR:.3e} Eh/A")
    check(eR < GATE, f"Richardson-extrapolated FD (h = 1e-4 / 2e-4 A): {eR:.3e} < {GATE:.0e} Eh/A")

print(f"{count[0]} curcuma runs, {failed} failure(s)")
sys.exit(1 if failed else 0)
