#!/usr/bin/env python3
"""Claude Generated (Oct 2026): GFN-FF fused D4 pair C6 (-gfnff.disp_c6_fused, default true) is bit-identical.

The fused path evaluates each D4 pair's C6 and dC6/dCN inside the CPU dispersion energy loop instead of storing
them per pair every step (D4ParameterGenerator::requestFusedPairC6). It performs the same floating-point
operations, so every result must equal the stored path (-gfnff.disp_c6_fused false) bit for bit:
  1. single point: -dump_gradient file (energy line + gradient) identical, 1 and 4 threads;
  2. 20 MD steps (reused calculator, per-step C6 refresh): trajectory file and MD log energies identical;
  3. the same MD with -d4param.d4_cn_cache_threshold 0.01, which exercises the weight refresh without
     derivative rebuild (the fused view then keeps the derivative state of the last rebuild).
Structures from the library (test_cases/structures), staged at run time.
Usage: check_d4_fused_c6.py CURCUMA_BINARY SOURCE_DIR
"""
import hashlib, os, re, subprocess, sys, tempfile

curcuma, src = sys.argv[1], sys.argv[2]
failed = 0
OFF = ["-gfnff.disp_c6_fused", "false"]


def check(ok, what):
    global failed
    print(("  PASS  " if ok else "  FAIL  ") + what, flush=True)
    failed += 0 if ok else 1


work = tempfile.mkdtemp(prefix="curcuma_d4_fused_c6_")
n_run = [0]


def stage(sid):
    p = subprocess.run([sys.executable, src + "/scripts/structlib.py", "path", sid], capture_output=True, text=True)
    return p.stdout.strip()


def fresh(xyz):
    n_run[0] += 1
    d = os.path.join(work, f"r{n_run[0]}")
    os.mkdir(d)
    with open(xyz) as f, open(os.path.join(d, "in.xyz"), "w") as g:
        g.write(f.read())
    return d


def sp(xyz, threads, extra):
    d = fresh(xyz)
    cmd = [curcuma, "-sp", "in.xyz", "-method", "gfnff", "-threads", str(threads), "-verbosity", "0", "-no_bmt",
           "-gfnff.cache_topology", "false", "-dump_gradient", "g.txt"] + extra
    subprocess.run(cmd, cwd=d, capture_output=True, text=True)
    p = os.path.join(d, "g.txt")
    return open(p, "rb").read() if os.path.exists(p) else None


def md(xyz, extra):
    d = fresh(xyz)
    cmd = [curcuma, "-md", "in.xyz", "-method", "gfnff", "-threads", "2", "-verbosity", "1", "-no_bmt",
           "-gfnff.cache_topology", "false", "-md.time_step", "0.5", "-md.max_time", "10", "-md.thermostat", "csvr",
           "-md.seed", "11", "-md.temperature", "400", "-md.print_frequency", "1", "-md.dump_frequency", "1",
           "-md.write_restart_frequency", "-1", "-md.no_restart", "true"] + extra
    out = subprocess.run(cmd, cwd=d, capture_output=True, text=True).stdout
    out = re.sub(r"\x1b\[[0-9;]*m", "", out)
    # MD table rows: lines whose first token is the simulation time; keep the energy columns, not the wall times
    rows = [l.split()[:6] for l in out.splitlines() if re.match(r"^\s*\d+\.\d+\s+-?\d", l)]
    trj = [os.path.join(r, f) for r, _, fs in os.walk(d) for f in fs if f.endswith(".trj.xyz")]  # <name>.snapshots/
    h = hashlib.md5(open(trj[0], "rb").read()).hexdigest() if trj else None
    return rows, h


for sid in ("caffeine", "trisaccharide-c18h32o16", "urea-water-cluster197.unrelaxed"):
    xyz = stage(sid)
    for t in (1, 4):
        a, b = sp(xyz, t, []), sp(xyz, t, OFF)
        check(a is not None and a == b, f"{sid}: -sp gradient file identical fused/stored, {t} thread(s)")

xyz = stage("urea-water-cluster197.unrelaxed")
for label, extra in (("default", []), ("d4_cn_cache_threshold 0.01", ["-d4param.d4_cn_cache_threshold", "0.01"])):
    (ra, ha), (rb, hb) = md(xyz, extra), md(xyz, extra + OFF)
    check(len(ra) >= 20 and ra == rb, f"MD 20 steps ({label}): {len(ra)} log rows, energies identical fused/stored")
    check(ha is not None and ha == hb, f"MD 20 steps ({label}): trajectory file identical fused/stored")

print(f"{failed} check(s) failed")
sys.exit(1 if failed else 0)
