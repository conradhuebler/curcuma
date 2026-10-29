#!/usr/bin/env python3
"""Claude Generated (Oct 2026): -lbfgs_max_atom_step (cap on the displacement of one atom per L-BFGS step).

Structure: library acetic-acid dimer (structlib id "acetic-acid-dimer"), every coordinate shifted by a fixed random
vector (sigma 0.12 A, seed 5) so that the first steps are large. -opt, default optimiser (lbfgspp), GFN-FF.
  1. Cap not given and cap 0.0 give bit-identical optimised coordinates (the cap is off by default).
  2. Without the cap the largest single-atom displacement in iterations 1-4 exceeds 0.05 A, so
     the test can see the cap (otherwise checks 3-4 would pass for any code).
  3. With -lbfgs_max_atom_step 0.05 no atom moves by more than 0.05 A (+1e-6) between consecutive frames.
  4. The capped run ends at the same minimum: energy within 1e-3 Eh of the uncapped one (measured 2.1e-4 Eh: both stop at the default convergence criteria; threshold set after that run).
Usage: check_lbfgs_atom_step_cap.py CURCUMA_BINARY SOURCE_DIR
"""
import os, subprocess, sys, tempfile
import numpy as np

curcuma, src = sys.argv[1], sys.argv[2]
failed = 0


def check(ok, what):
    global failed
    print(("  PASS  " if ok else "  FAIL  ") + what)
    failed += 0 if ok else 1


path = subprocess.run([sys.executable, src + "/scripts/structlib.py", "path", "acetic-acid-dimer"],
                      capture_output=True, text=True).stdout.strip().split("\n")[-1]
lines = open(path).read().split("\n")
n = int(lines[0])
el = [l.split()[0] for l in lines[2:2 + n]]
X = np.array([[float(v) for v in l.split()[1:4]] for l in lines[2:2 + n]])
X = X + np.random.RandomState(5).normal(0.0, 0.12, X.shape)

work = tempfile.mkdtemp(prefix="curcuma_lbfgs_cap_")


def frames(fn):
    L = open(fn).read().split("\n")
    out, i = [], 0
    while i < len(L) and L[i].strip().isdigit():
        k = int(L[i])
        out.append((L[i + 1], np.array([[float(v) for v in l.split()[1:4]] for l in L[i + 2:i + 2 + k]])))
        i += k + 2
    return out


def energy(comment):
    for tok in comment.replace("=", " ").replace(",", " ").split():
        try:
            v = float(tok)
            if v < -1.0:
                return v
        except ValueError:
            pass
    return None


def opt(tag, extra, iters=None):
    """one -opt run in its own directory; returns (coordinates, energy) of in.opt.xyz or None"""
    d = os.path.join(work, tag)
    os.mkdir(d)
    with open(os.path.join(d, "in.xyz"), "w") as f:
        f.write(f"{n}\nx\n" + "".join("%s %.10f %.10f %.10f\n" % (e, *x) for e, x in zip(el, X)))
    cmd = [curcuma, "-opt", "in.xyz", "-method", "gfnff", "-verbosity", "0", "-no_bmt", "-threads", "1",
           "-gfnff.cache_topology", "false"] + extra
    if iters is not None:
        cmd += ["-opt.max_iterations", str(iters)]
    subprocess.run(cmd, cwd=d, capture_output=True, text=True)
    fn = os.path.join(d, "in.opt.xyz")
    if not os.path.exists(fn):
        return None
    last = frames(fn)[-1]
    return last[1], energy(last[0])


def run(tag, extra):
    """converged run plus the largest single-atom step of the first 4 iterations (iteration k minus k-1)"""
    full = opt(tag + "_full", extra)
    prev, maxstep = X, 0.0
    for k in (1, 2, 3, 4):
        r = opt(f"{tag}_it{k}", extra, k)
        if r is None:
            return None
        maxstep = max(maxstep, np.linalg.norm(r[0] - prev, axis=1).max())
        prev = r[0]
    return None if full is None else {"x": full[0], "e": full[1], "maxstep": maxstep}


r_none = run("none", [])
r_zero = run("zero", ["-lbfgs_max_atom_step", "0.0"])
r_cap = run("cap", ["-lbfgs_max_atom_step", "0.05"])
for tag, r in (("none", r_none), ("zero", r_zero), ("cap", r_cap)):
    check(r is not None, f"{tag}: optimisation wrote result")
if failed == 0:
    print("  info  max single-atom step in iterations 1-4: none %.4f, cap %.4f A; E none %s cap %s" % (
        r_none["maxstep"], r_cap["maxstep"], r_none["e"], r_cap["e"]))
    check(np.array_equal(r_none["x"], r_zero["x"]), "1. cap 0.0 == no flag: optimised coordinates bit-identical")
    check(r_none["maxstep"] > 0.05, "2. uncapped run has a single-atom step above 0.05 A (test is sensitive)")
    check(r_cap["maxstep"] <= 0.05 + 1e-6, "3. capped run: largest single-atom step <= 0.05 A")
    check(r_none["e"] is not None and r_cap["e"] is not None and abs(r_none["e"] - r_cap["e"]) < 1e-3,
          "4. capped run reaches the same minimum (|dE| < 1e-3 Eh)")
print("RESULT:", "FAIL (%d)" % failed if failed else "PASS")
sys.exit(1 if failed else 0)
