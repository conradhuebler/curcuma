#!/usr/bin/env python3
# Claude Generated (Oct 2026). AI-generated, machine-tested only.
"""Property tests of one structure with one method: things that must hold for any correct energy/gradient code,
without a stored reference value.

Usage:
    struct_invariants.py --curcuma BIN --structure ID --method gfn2 [--charge Q] [--spin S] [--tests a,b,...]
                         [--gpu cuda] [--tol-file test_cases/struct_test_policy.json] [--verbose]

Tests (names as in `--tests`):
  translate   energy and gradient unchanged by a rigid translation
  rotate      energy unchanged, gradient rotates with the molecule (g' = R g)
  permute     energy unchanged by a fixed shuffle of the atom order
  sumforce    sum of the forces = 0 (translation invariance of the energy)
  torque      sum of r x g = 0 (rotation invariance of the energy)
  fd          analytic gradient against a central finite difference of curcuma's own energy (Angstrom steps)
  repeat      two identical CPU runs give the same energy and gradient
  threads     -threads 1 against -threads 4
  gpu         CPU against the GPU backend named by --gpu (skipped, exit 77, when it is not available)
  compressed      the geometry scaled about its centroid (policy compressed_scale, CN far beyond every D4 reference): the
                  finite-difference gradient test and the force sum again (Oct 2026: a D4 weight underflow gave wrong gradients here)
  compressed_gpu  CPU against the GPU backend on the scaled geometry

Units: curcuma's -dump_gradient is Eh/Angstrom (the getGradient() contract, root CLAUDE.md). Exit 0 = all selected tests
passed, 1 = a test failed, 77 = nothing could be run (skip).
"""
import argparse
import json
import math
import os
import random
import re
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_TESTS = "translate,rotate,permute,sumforce,torque,fd,repeat,threads"


def read_xyz(path):
    with open(path) as f:
        lines = f.read().split("\n")
    n = int(lines[0].split()[0])
    atoms = []
    for l in lines[2:2 + n]:
        p = l.split()
        atoms.append([p[0], float(p[1]), float(p[2]), float(p[3])])
    return atoms


def write_xyz(path, atoms, comment="struct_invariants"):
    with open(path, "w") as f:
        f.write("%d\n%s\n" % (len(atoms), comment))
        for a in atoms:
            f.write("%-3s %.12f %.12f %.12f\n" % (a[0], a[1], a[2], a[3]))


def structure_path(sid):
    out = subprocess.run([sys.executable, os.path.join(REPO, "scripts", "structlib.py"), "path", sid],
                         capture_output=True, text=True)
    p = out.stdout.strip().split("\n")[0] if out.returncode == 0 else ""
    if not p or not os.path.exists(p):
        raise SystemExit("structure '%s' not found in the library: %s" % (sid, out.stderr.strip()))
    return p


class Runner:
    def __init__(self, args):
        self.curcuma = args.curcuma
        self.method = args.method
        self.charge = args.charge
        self.spin = args.spin
        self.verbose = args.verbose
        self.tmp = tempfile.mkdtemp(prefix="structinv_")
        self.n = 0
        # tight SCF so that the finite difference is not limited by the convergence threshold
        self.extra = ["-scf_threshold", "1e-10"] if args.method in ("gfn1", "gfn2") else []

    def calc(self, atoms, threads=1, gpu="cpu", tag="run"):
        """Energy (Eh) and gradient (list of [gx, gy, gz], Eh/Angstrom) at the given atoms; None when curcuma fails."""
        self.n += 1
        d = os.path.join(self.tmp, "%s_%d" % (tag, self.n))   # fresh directory per call: no topology cache is reused
        os.makedirs(d)
        write_xyz(os.path.join(d, "s.xyz"), atoms)
        cmd = [self.curcuma, "-sp", "s.xyz", "-method", self.method, "-charge", str(self.charge), "-spin", str(self.spin),
               "-gpu", gpu, "-threads", str(threads), "-dump_gradient", "g.txt", "-no_bmt", "-verbosity", "0"] + self.extra
        if gpu != "cpu":
            cmd += ["-gpu_devices", "0", "-gpu_strict", "true"]
        r = subprocess.run(cmd, cwd=d, capture_output=True, text=True, timeout=600)
        gfile = os.path.join(d, "g.txt")
        if r.returncode != 0 or not os.path.exists(gfile):
            if self.verbose:
                print("   curcuma exit %d: %s" % (r.returncode, (r.stdout + r.stderr).strip()[-300:]))
            return None
        e, g = None, []
        for l in open(gfile):
            if l.startswith("# energy"):
                e = float(l.split()[2])
            elif not l.startswith("#") and l.strip():
                g.append([float(x) for x in l.split()])
        if e is None or len(g) != len(atoms):
            return None
        return e, g


def rot_matrix(axis, angle):
    n = math.sqrt(sum(a * a for a in axis))
    x, y, z = [a / n for a in axis]
    c, s = math.cos(angle), math.sin(angle)
    C = 1 - c
    return [[c + x * x * C, x * y * C - z * s, x * z * C + y * s],
            [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
            [z * x * C - y * s, z * y * C + x * s, c + z * z * C]]


def apply(R, v):
    return [sum(R[i][j] * v[j] for j in range(3)) for i in range(3)]


def maxdiff(a, b):
    return max(abs(x - y) for ra, rb in zip(a, b) for x, y in zip(ra, rb))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--curcuma", required=True)
    ap.add_argument("--structure", required=True)
    ap.add_argument("--method", required=True)
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--spin", type=int, default=0)
    ap.add_argument("--tests", default=DEFAULT_TESTS)
    ap.add_argument("--gpu", default="")
    ap.add_argument("--policy", default=os.path.join(REPO, "test_cases", "struct_test_policy.json"))
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    policy = json.load(open(args.policy))
    tol_all = policy["tolerances"]
    tol = dict(tol_all.get("default", {}))
    tol.update(tol_all.get(args.method, {}))

    atoms = read_xyz(structure_path(args.structure))
    R = Runner(args)
    base = R.calc(atoms, tag="base")
    if base is None:
        print("SKIP %s %s: the reference run failed (see --verbose)" % (args.structure, args.method))
        return 77
    e0, g0 = base
    n = len(atoms)
    failed = 0
    results = []

    # Known deviations of the policy (structure + method + test): a looser, explicit limit and the reason, printed with the result.
    known = {k["test"]: k for k in policy.get("known_deviations", [])
             if k["structure"] == args.structure and k["method"] == args.method}

    def report(name, value, limit, note=""):
        nonlocal failed
        if name in known:
            limit = known[name]["limit"]
            note = (note + " [known deviation: %s]" % known[name]["why"]).strip()
        ok = value <= limit
        failed += 0 if ok else 1
        print("%s %-9s %s %s: %.3e (limit %.1e)%s" % ("PASS" if ok else "FAIL", name, args.structure, args.method,
                                                      value, limit, " " + note if note else ""))

    wanted = [t for t in args.tests.split(",") if t]

    if "translate" in wanted:
        sh = [0.37, -1.21, 2.03]
        r = R.calc([[a[0], a[1] + sh[0], a[2] + sh[1], a[3] + sh[2]] for a in atoms], tag="tr")
        if r is None:
            report("translate", float("inf"), tol["energy"], "run failed")
        else:
            report("translate", max(abs(r[0] - e0), 0.0), tol["energy"], "|dE| Eh")
            report("translate", maxdiff(r[1], g0), tol["gradient"], "max|dg| Eh/A")

    if "rotate" in wanted:
        M = rot_matrix([1.0, 2.0, 3.0], 0.7)
        r = R.calc([[a[0]] + apply(M, a[1:4]) for a in atoms], tag="rot")
        if r is None:
            report("rotate", float("inf"), tol["energy"], "run failed")
        else:
            report("rotate", abs(r[0] - e0), tol["energy"], "|dE| Eh")
            report("rotate", maxdiff(r[1], [apply(M, g) for g in g0]), tol["gradient"], "max|g' - R g| Eh/A")

    if "permute" in wanted:
        idx = list(range(n))
        random.Random(11).shuffle(idx)
        r = R.calc([atoms[i] for i in idx], tag="perm")
        if r is None:
            report("permute", float("inf"), tol["energy"], "run failed")
        else:
            lim = tol.get("permute_energy", tol["energy"])
            report("permute", abs(r[0] - e0), lim, "|dE| Eh" + (" (known issue #34)" if "permute_energy" in tol else ""))

    if "sumforce" in wanted:
        s = max(abs(sum(g[k] for g in g0)) for k in range(3))
        report("sumforce", s, tol["sumforce"], "max|sum g| Eh/A")

    if "torque" in wanted:
        cx = [sum(a[k + 1] for a in atoms) / n for k in range(3)]
        t = [0.0, 0.0, 0.0]
        for a, g in zip(atoms, g0):
            r_ = [a[k + 1] - cx[k] for k in range(3)]
            t[0] += r_[1] * g[2] - r_[2] * g[1]
            t[1] += r_[2] * g[0] - r_[0] * g[2]
            t[2] += r_[0] * g[1] - r_[1] * g[0]
        report("torque", max(abs(x) for x in t), tol["torque"], "max|sum r x g| Eh")

    if "fd" in wanted:
        # Central differences with Richardson extrapolation of the h^2 term, g = (4 D(h) - D(2h)) / 3: a stiff, anharmonic
        # bond (N2: third derivative) leaves a truncation error of 1e-4 Eh/A at h = 0.005 A with the plain formula.
        h = tol.get("fd_step", 0.005)
        pick = sorted({0, n // 2, n - 1})
        worst = 0.0
        bad = False
        for i in pick:
            for k in range(3):
                D = {}
                for step in (h, 2 * h):
                    e = []
                    for sgn in (+1, -1):
                        a2 = [list(a) for a in atoms]
                        a2[i][k + 1] += sgn * step
                        r = R.calc(a2, tag="fd")
                        if r is None:
                            bad = True
                            break
                        e.append(r[0])
                    if bad:
                        break
                    D[step] = (e[0] - e[1]) / (2 * step)
                if bad:
                    break
                worst = max(worst, abs((4 * D[h] - D[2 * h]) / 3 - g0[i][k]))
            if bad:
                break
        if bad:
            report("fd", float("inf"), tol["fd"], "run failed")
        else:
            report("fd", worst, tol["fd"], "max|g_fd - g| Eh/A, Richardson, %d atoms x 3, h=%g A" % (len(pick), h))

    if "repeat" in wanted:
        r = R.calc(atoms, tag="rep")
        if r is None:
            report("repeat", float("inf"), tol["repeat"], "run failed")
        else:
            report("repeat", max(abs(r[0] - e0), maxdiff(r[1], g0)), tol["repeat"], "max(|dE|, |dg|)")

    if "threads" in wanted:
        r = R.calc(atoms, threads=4, tag="thr")
        if r is None:
            report("threads", float("inf"), tol["threads"], "run failed")
        else:
            report("threads", max(abs(r[0] - e0), maxdiff(r[1], g0)), tol["threads"], "max(|dE|, |dg|) 1 vs 4 threads")

    def scaled_atoms():
        f = tol.get("compressed_scale", 0.35)
        c = [sum(a[k + 1] for a in atoms) / n for k in range(3)]
        return [[a[0]] + [c[k] + f * (a[k + 1] - c[k]) for k in range(3)] for a in atoms], f

    if "compressed" in wanted or "compressed_gpu" in wanted:
        sc, f = scaled_atoms()
        base_c = R.calc(sc, tag="comp")
        if base_c is None:
            print("SKIP compressed %s %s: the run failed" % (args.structure, args.method))
        else:
            ec, gc = base_c
            if "compressed" in wanted:
                h = tol.get("fd_step", 0.005) / 5.0   # compressed bonds are stiffer
                pick = sorted({0, n // 2, n - 1})
                worst, bad = 0.0, False
                for i in pick:
                    for k in range(3):
                        D = {}
                        for step in (h, 2 * h):
                            e = []
                            for sgn in (+1, -1):
                                a2 = [list(a) for a in sc]
                                a2[i][k + 1] += sgn * step
                                r = R.calc(a2, tag="fdc")
                                if r is None:
                                    bad = True
                                    break
                                e.append(r[0])
                            if bad:
                                break
                            D[step] = (e[0] - e[1]) / (2 * step)
                        if bad:
                            break
                        worst = max(worst, abs((4 * D[h] - D[2 * h]) / 3 - gc[i][k]))
                    if bad:
                        break
                gmax = max(abs(x) for row in gc for x in row)
                if bad:
                    report("compressed", float("inf"), tol["fd_compressed"], "run failed")
                else:
                    report("compressed", worst, tol["fd_compressed"],
                           "fd gradient at scale %.2f, max|g_fd - g| Eh/A (max|g| %.1e)" % (f, gmax))
                report("compsum", max(abs(sum(g[k] for g in gc)) for k in range(3)), tol["compressed_sumforce"],
                       "max|sum g| Eh/A at scale %.2f" % f)
            if "compressed_gpu" in wanted and args.gpu:
                r = R.calc(sc, gpu=args.gpu, tag="compgpu")
                if r is None:
                    print("SKIP compressed_gpu %s %s: -gpu %s not available or the run failed" % (args.structure, args.method, args.gpu))
                else:
                    report("compgpu", abs(r[0] - ec), tol["gpu_energy"], "|dE| Eh at scale %.2f, cpu vs %s" % (f, args.gpu))
                    report("compgpu", maxdiff(r[1], gc), tol["compressed_gpu_gradient"],
                           "max|dg| Eh/A at scale %.2f, cpu vs %s" % (f, args.gpu))

    if "gpu" in wanted and args.gpu:
        r = R.calc(atoms, gpu=args.gpu, tag="gpu")
        if r is None:
            print("SKIP gpu %s %s: -gpu %s not available or the run failed" % (args.structure, args.method, args.gpu))
        else:
            report("gpu", abs(r[0] - e0), tol["gpu_energy"], "|dE| Eh, cpu vs " + args.gpu)
            report("gpu", maxdiff(r[1], g0), tol["gpu_gradient"], "max|dg| Eh/A, cpu vs " + args.gpu)

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
