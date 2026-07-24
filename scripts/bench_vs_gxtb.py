#!/usr/bin/env python3
"""Wall-clock benchmark: curcuma vs gxtb for gfnff/gfn1/gfn2 on CPU (Claude Generated).

Sweeps core counts. Times a single-point (+ optional gradient) as full-process wall
clock (what the user experiences), median of N runs after a warmup. curcuma is pinned
to K threads via -threads K + OMP/MKL_NUM_THREADS=K; gxtb via OMP_NUM_THREADS=K.

Usage: bench_vs_gxtb.py <xyz> [methods=gfnff,gfn1,gfn2] [cores=1,2,4,8,16] [reps=3] [grad=1]
"""
import subprocess, sys, time, os, statistics, tempfile, shutil

CURCUMA = "/home/conrad/src/curcuma/release/curcuma"
GXTB    = "/opt/bin/gxtb"
GX_FLAG = {"gfnff": "--gfnff", "gfn1": "--gfn1", "gfn2": "--gfn2"}

def median_wall(cmd, env, reps, cwd=None):
    r = subprocess.run(cmd, env=env, cwd=cwd, capture_output=True)  # warmup
    if r.returncode != 0:
        return None
    ts = []
    for _ in range(reps):
        t = time.perf_counter()
        subprocess.run(cmd, env=env, cwd=cwd, capture_output=True)
        ts.append(time.perf_counter() - t)
    return statistics.median(ts)

def main():
    xyz = os.path.abspath(sys.argv[1])
    methods = (sys.argv[2].split(",") if len(sys.argv) > 2 else ["gfnff","gfn1","gfn2"])
    cores = ([int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else [1,2,4,8,16])
    reps = int(sys.argv[4]) if len(sys.argv) > 4 else 3
    grad = (len(sys.argv) <= 5) or sys.argv[5] == "1"
    nat = int(open(xyz).readline().split()[0])
    print(f"# {xyz}  nat={nat}  grad={grad}  reps={reps}  cores={cores}  (ms, median)")
    hdr = "  ".join(f"K={k:<2}" for k in cores)
    print(f"{'method':>6} {'engine':>7} | {hdr}")

    for m in methods:
        cbase = [CURCUMA, "-sp", xyz, "-method", m, "-charge", "0", "-verbosity", "0", "-no_bmt"]
        if grad: cbase += ["-gradient"]
        # gxtb: fresh xtbrestart-free cwd for EVERY run (+ --norestart), so each is a cold
        # start matching curcuma's fresh -sp. Reusing xtbrestart warm-starts gxtb to ~3
        # iterations and grossly understates its time.
        gbase = [GXTB, "mol.xyz", GX_FLAG[m], "--norestart"] + (["--grad"] if grad else [])
        def gx_wall(k):
            genv = {**os.environ, "OMP_NUM_THREADS": str(k)}
            ts = []
            for i in range(reps + 1):                       # +1 warmup
                wd = tempfile.mkdtemp(); shutil.copy(xyz, wd+"/mol.xyz")
                t = time.perf_counter()
                r = subprocess.run(gbase, env=genv, cwd=wd, capture_output=True)
                dt = time.perf_counter() - t
                shutil.rmtree(wd, ignore_errors=True)
                if r.returncode != 0: return None
                if i > 0: ts.append(dt)
            return statistics.median(ts)
        cur, gx = [], []
        for k in cores:
            env = {**os.environ, "OMP_NUM_THREADS": str(k), "MKL_NUM_THREADS": str(k)}
            cur.append(median_wall(cbase + ["-threads", str(k)], env, reps))
            gx.append(gx_wall(k))
        def row(vals): return "  ".join(f"{v*1000:5.0f}" if v else " FAIL" for v in vals)
        print(f"{m:>6} {'curcuma':>7} | {row(cur)}")
        print(f"{m:>6} {'gxtb':>7} | {row(gx)}")
        rr = "  ".join((f"{c/g:5.2f}" if (c and g) else "  -  ") for c,g in zip(cur,gx))
        print(f"{m:>6} {'cur/gx':>7} | {rr}")

if __name__ == "__main__":
    main()
