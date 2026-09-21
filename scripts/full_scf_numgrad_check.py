#!/usr/bin/env python3
"""Full-SCF numerical gradient vs curcuma analytic gradient (Claude Generated, F5).

Runs curcuma -sp -gradient at displaced geometries, differences the TOTAL
(relaxed-density) energy read from the -dump_gradient header (12 digits), and
compares to the analytic gradient. Unlike the frozen-density audit this exercises
every response term (W = C·fε·Cᵀ energy-weighted density, D4/multipole/Coulomb
self-consistent q-response). Ground-truth test of the analytic gradient.

The native xTB gradient is in Eh/Angstrom (native_xtb_method.cpp:331); the FD
is computed as dE/dx with x in Angstrom, so both sides are Eh/Angstrom directly.

Usage: full_scf_numgrad_check.py <xyz> [method] [charge] [h1,h2,...] [atoms]
"""
import subprocess, sys, re, tempfile, os, math
from pathlib import Path

CURCUMA = "./release/curcuma"
E_HDR_RE = re.compile(r"#\s*energy\s*(-?\d+\.\d+)\s*Eh")

def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    els, xyz = [], []
    for ln in lines[2:2+n]:
        p = ln.split(); els.append(p[0])
        xyz.append([float(p[1]), float(p[2]), float(p[3])])
    return els, xyz

def write_xyz(path, els, xyz):
    with open(path, "w") as f:
        f.write(f"{len(els)}\n\n")
        for e, r in zip(els, xyz):
            f.write(f"{e} {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}\n")

def run(xyzpath, method, charge, scf_thr, gradpath):
    """Return (energy_12digit, analytic_grad_rows[Eh/Ang])."""
    subprocess.run(
        [CURCUMA, "-sp", xyzpath, "-method", method, "-charge", str(charge),
         "-gradient", "-dump_gradient", gradpath, "-verbosity", "1", "-no_bmt",
         "-scf_threshold", scf_thr], capture_output=True, text=True)
    txt = Path(gradpath).read_text()
    e = None; rows = []
    for ln in txt.splitlines():
        m = E_HDR_RE.search(ln)
        if m: e = float(m.group(1))
        if ln.startswith("#"): continue
        p = ln.split()
        if len(p) == 3: rows.append([float(x) for x in p])
    return e, rows

def main():
    struct = sys.argv[1] if len(sys.argv) > 1 else "test_cases/MOR41-testset/ED39/mol.xyz"
    method = sys.argv[2] if len(sys.argv) > 2 else "gfn2"
    charge = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    hs = [float(x) for x in sys.argv[4].split(",")] if len(sys.argv) > 4 else [4e-3, 2e-3, 1e-3]
    atoms = [int(x) for x in sys.argv[5].split(",")] if len(sys.argv) > 5 else None
    scf_thr = "1e-10"

    els, xyz0 = read_xyz(struct)
    nat = len(els)
    if atoms is None: atoms = list(range(nat))
    tmp = tempfile.mkdtemp(prefix="fdgrad_")

    _, ana = run(struct, method, charge, scf_thr, os.path.join(tmp, "ana.grad"))
    print(f"# {struct}  method={method}  charge={charge}  nat={nat}  scf={scf_thr}")
    print(f"# analytic gradient in Eh/Ang (native_xtb_method.cpp:331)")
    print(f"# {'atom':>4} {'el':>3} {'c':>1} {'analytic':>16}  " +
          "  ".join(f"FD(h={h:g})".rjust(15) for h in hs) +
          f"  {'rel@min_h':>10}")

    def energy_at(disp_i, disp_j, sign, h):
        x = [r[:] for r in xyz0]; x[disp_i][disp_j] += sign*h
        pth = os.path.join(tmp, f"d_{disp_i}_{disp_j}_{sign}.xyz")
        write_xyz(pth, els, x)
        e, _ = run(pth, method, charge, scf_thr, os.path.join(tmp, "g.grad"))
        return e

    max_rel = 0.0; worst = None
    for i in atoms:
        for j in range(3):
            a = ana[i][j]
            fds = []
            for h in hs:
                ep = energy_at(i, j, +1, h)
                em = energy_at(i, j, -1, h)
                fds.append((ep - em) / (2.0*h))   # Eh/Ang
            fd_min = fds[-1]  # smallest h
            rel = abs(a - fd_min) / max(1e-6, abs(fd_min))
            if abs(fd_min) > 1e-3 and rel > max_rel:
                max_rel = rel; worst = (i, els[i], 'xyz'[j])
            print(f"  {i:>4} {els[i]:>3} {'xyz'[j]:>1} {a:>16.10f}  " +
                  "  ".join(f"{v:>15.10f}" for v in fds) +
                  f"  {rel:>10.2e}")
    print(f"# max_rel (components > 1e-3 Eh/Ang) at smallest h = {max_rel:.3e}"
          + (f" @atom {worst}" if worst else ""))

if __name__ == "__main__":
    main()
