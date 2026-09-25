#!/usr/bin/env python3
"""Run two curcuma binaries over a reference set and compare the energies per structure.

Answers "did my change alter any reference-set energy" WITHOUT needing xtb: the
old binary is the reference. Build it from the commit you started at, e.g.

    git worktree add /tmp/ref <commit> && cd /tmp/ref && mkdir release && cd release
    cmake .. -DCMAKE_BUILD_TYPE=Release && make -j16 curcuma

    python scripts/refset_regression.py --set mor41   --method gfnff \
        --old /tmp/ref/release/curcuma --new release/curcuma
    python scripts/refset_regression.py --set gmtkn55 --method gfn2  ...

Fetch the sets first with scripts/fetch_testset.py. Every single point runs
single-threaded and the structures are spread over a thread pool; each one gets
its own scratch directory, because every GMTKN55 structure file is called
struc.xyz and GFN-FF caches its topology next to that name (CLAUDE.md #28).

A change that is meant to be numerically neutral should come back with zero
differences; a rounding-level change shows up in the last printed digits. Note
that a non-converged SCF (GMTKN55 has three: the G21IP single-atom cations)
moves by tens of kcal/mol on any perturbation and says nothing about the change.

Claude Generated (Sep 2026).
"""
import argparse, concurrent.futures as cf, os, re, shutil, subprocess, sys, tempfile
from pathlib import Path

E_RE = re.compile(r"Single Point Energy\s*=\s*(-?\d+\.\d+)")

def run_one(binary, xyz, method, charge=0, spin=0):
    # Fresh scratch dir per structure: GFN-FF caches its topology next to the basename
    # and every GMTKN55 structure is called struc.xyz (see CLAUDE.md Known Issue #28).
    with tempfile.TemporaryDirectory() as td:
        local = Path(td) / "struc.xyz"
        shutil.copy(xyz, local)
        cmd = [str(binary), "-sp", str(local), "-method", method, "-threads", "1",
               "-no_bmt", "-verbosity", "1"]
        if charge: cmd += ["-charge", str(charge)]
        if spin:   cmd += ["-spin", str(spin)]
        try:
            out = subprocess.run(cmd, capture_output=True, text=True, timeout=900).stdout
        except subprocess.TimeoutExpired:
            return None
    m = E_RE.findall(re.sub(r"\x1b\[[0-9;]*m", "", out))
    return float(m[-1]) if m else None

def charge_spin(xyz):
    d = xyz.parent
    c = s = 0
    for name, conv in ((".CHRG", "c"), (".UHF", "s")):
        f = d / name
        if f.exists():
            try:
                v = int(f.read_text().split()[0])
            except Exception:
                v = 0
            if conv == "c": c = v
            else: s = v
    return c, s

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", required=True, choices=["mor41", "gmtkn55"])
    ap.add_argument("--method", default="gfnff")
    ap.add_argument("--old", required=True)
    ap.add_argument("--new", required=True)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    a = ap.parse_args()

    repo = Path(__file__).resolve().parents[1]
    if a.set == "mor41":
        root = repo / "test_cases" / "MOR41-testset"
        xyzs = sorted(p for p in root.rglob("*.xyz"))
    else:
        root = repo / "test_cases" / "GMTKN55-testset"
        xyzs = sorted(root.rglob("struc.xyz"))
    if a.limit:
        xyzs = xyzs[:a.limit]
    print(f"{a.set}/{a.method}: {len(xyzs)} structures, {a.jobs} jobs", flush=True)

    def both(xyz):
        c, s = charge_spin(xyz)
        return (xyz, run_one(a.old, xyz, a.method, c, s), run_one(a.new, xyz, a.method, c, s))

    worst, n_ok, n_fail, n_diff = [], 0, 0, 0
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for xyz, e_old, e_new in ex.map(both, xyzs):
            if e_old is None or e_new is None:
                n_fail += 1
                print(f"  FAIL {xyz.relative_to(root)}: old={e_old} new={e_new}", flush=True)
                continue
            n_ok += 1
            d = abs(e_new - e_old) * 627.5094740631  # kcal/mol
            if d > 0.0:
                n_diff += 1
            worst.append((d, str(xyz.relative_to(root)), e_old, e_new))
    worst.sort(reverse=True)
    print(f"\n{a.set}/{a.method}: {n_ok} compared, {n_fail} failed, {n_diff} with any difference")
    if worst:
        mad = sum(w[0] for w in worst) / len(worst)
        print(f"  MAD {mad:.3e} kcal/mol, max {worst[0][0]:.3e} kcal/mol")
        for d, name, eo, en in worst[:5]:
            if d == 0.0: break
            print(f"    {d:.3e}  {name}  {eo:.10f} -> {en:.10f}")
    return 1 if n_fail else 0

sys.exit(main())
