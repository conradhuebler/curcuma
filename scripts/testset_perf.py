#!/usr/bin/env python3
# Claude Generated (Sep 2026)
"""CPU/GPU/threading performance benchmark on a fetched testset.

Generic wall-clock timing harness: picks the N largest structures (by atom
count) under a testset directory (e.g. test_cases/MOR41-testset, populated by
scripts/fetch_testset.py) and times `curcuma -sp`/`-opt` across a grid of
`-threads N` values and `-gpu backend` values, on this machine. Reports
threading speedup/parallel efficiency and GPU-vs-CPU speedup.

GPU backends are runtime dlopen plugins (release/libcurcuma_<backend>.so); by
default the script only benchmarks backends whose plugin is actually present
next to the curcuma binary on this machine (plus "cpu"), so results are never
silently CPU-only-pretending-to-be-GPU. Pass --gpu to override explicitly.

Usage:
    python scripts/testset_perf.py                                   # MOR41, defaults
    python scripts/testset_perf.py --root test_cases/MOR41-testset \\
        --method gfn2 --threads 1,2,4,8 --n-structures 5
    python scripts/testset_perf.py --gpu cpu,cuda --opt --repeats 5
"""
import argparse
import csv
import statistics
import subprocess
import time
from datetime import datetime
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CURCUMA = REPO / "release" / "curcuma"
RELEASE_DIR = REPO / "release"
KNOWN_GPU_BACKENDS = ("cuda", "rocm", "vulkan")


def count_atoms(xyz_path):
    try:
        with xyz_path.open() as f:
            return int(f.readline().strip())
    except (OSError, ValueError):
        return 0


def detect_gpu_backends():
    """GPU backends whose plugin .so actually sits next to the binary."""
    found = []
    for backend in KNOWN_GPU_BACKENDS:
        if (RELEASE_DIR / f"libcurcuma_{backend}.so").exists():
            found.append(backend)
    return found


def pick_structures(root, n):
    xyz_files = sorted(root.glob("**/*.xyz"))
    scored = [(count_atoms(p), p) for p in xyz_files]
    scored = [(a, p) for a, p in scored if a > 0]
    scored.sort(key=lambda t: -t[0])
    return scored[:n]


def run_once(xyz_path, method, threads, backend, use_opt, timeout):
    cmd = [str(CURCUMA), "-opt" if use_opt else "-sp", str(xyz_path),
           "-method", method, "-charge", "0", "-verbosity", "0", "-no_bmt",
           "-threads", str(threads)]
    if backend != "cpu":
        cmd += ["-gpu", backend]
    start = time.perf_counter()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        ok = proc.returncode == 0
        stderr = proc.stderr
    except subprocess.TimeoutExpired:
        return None, "TIMEOUT"
    wall = time.perf_counter() - start
    return (wall if ok else None), ("" if ok else stderr[-500:])


def median_run(xyz_path, method, threads, backend, use_opt, repeats, timeout):
    walls = []
    last_err = ""
    for _ in range(repeats):
        wall, err = run_once(xyz_path, method, threads, backend, use_opt, timeout)
        if wall is None:
            last_err = err
            continue
        walls.append(wall)
    if not walls:
        return None, last_err
    return statistics.median(walls), ""


def stats_line(rows, key):
    vals = [r[key] for r in rows if r[key] is not None]
    return vals


def write_csv(rows, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["structure", "atoms", "method", "backend", "threads",
                    "wall_s_median", "error"])
        for r in rows:
            w.writerow([r["structure"], r["atoms"], r["method"], r["backend"],
                        r["threads"], f"{r['wall_s']:.4f}" if r["wall_s"] is not None else "",
                        r["error"]])


def write_summary(rows, structures, threads_list, backends, path):
    with path.open("w") as f:
        f.write("# curcuma testset performance benchmark\n\n")
        f.write(f"Structures: {len(structures)} (largest by atom count), "
                f"threads={threads_list}, backends={backends}\n\n")

        f.write("## Threading scaling (per backend, largest structure)\n\n")
        largest = structures[0][1].name if structures else None
        for backend in backends:
            f.write(f"### {backend}\n\n")
            f.write("| threads | wall_s | speedup vs threads=1 | parallel efficiency |\n")
            f.write("|---:|---:|---:|---:|\n")
            base = None
            for t in threads_list:
                match = [r for r in rows if r["structure"] == largest
                         and r["backend"] == backend and r["threads"] == t]
                wall = match[0]["wall_s"] if match else None
                if t == threads_list[0]:
                    base = wall
                if wall is None:
                    f.write(f"| {t} | - | - | - |\n")
                    continue
                speedup = (base / wall) if base else None
                eff = (speedup / t) if speedup is not None else None
                f.write(f"| {t} | {wall:.3f} | "
                        f"{f'{speedup:.2f}x' if speedup is not None else '-'} | "
                        f"{f'{eff:.2f}' if eff is not None else '-'} |\n")
            f.write("\n")

        gpu_backends = [b for b in backends if b != "cpu"]
        if gpu_backends:
            f.write("## GPU vs CPU (best CPU threading vs each GPU backend, per structure)\n\n")
            f.write("| structure | atoms | best CPU wall_s (threads) | "
                    + " | ".join(f"{b} wall_s" for b in gpu_backends)
                    + " | " + " | ".join(f"{b} speedup" for b in gpu_backends) + " |\n")
            f.write("|---|---:|---:|" + "---:|" * len(gpu_backends) * 2 + "\n")
            for atoms, xyz in structures:
                s = xyz.name
                cpu_rows = [r for r in rows if r["structure"] == s
                            and r["backend"] == "cpu" and r["wall_s"] is not None]
                best_cpu = min(cpu_rows, key=lambda r: r["wall_s"]) if cpu_rows else None
                cells_wall, cells_speedup = [], []
                for b in gpu_backends:
                    grow = [r for r in rows if r["structure"] == s and r["backend"] == b
                            and r["wall_s"] is not None]
                    gwall = min(r["wall_s"] for r in grow) if grow else None
                    cells_wall.append(f"{gwall:.3f}" if gwall is not None else "-")
                    if gwall and best_cpu:
                        cells_speedup.append(f"{best_cpu['wall_s'] / gwall:.2f}x")
                    else:
                        cells_speedup.append("-")
                bc = (f"{best_cpu['wall_s']:.3f} ({best_cpu['threads']})"
                      if best_cpu else "-")
                f.write(f"| {s} | {atoms} | {bc} | "
                        + " | ".join(cells_wall) + " | " + " | ".join(cells_speedup) + " |\n")
            f.write("\n")

        errors = [r for r in rows if r["error"]]
        if errors:
            f.write(f"## Failures ({len(errors)})\n\n")
            for r in errors:
                f.write(f"- {r['structure']} method={r['method']} backend={r['backend']} "
                        f"threads={r['threads']}: {r['error'][:200]}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=REPO / "test_cases" / "MOR41-testset",
                     help="testset directory to search for *.xyz (default: MOR41-testset)")
    ap.add_argument("--method", default="gfn2")
    ap.add_argument("--threads", default="1,2,4,8",
                     help="comma-separated thread counts")
    ap.add_argument("--gpu", default="auto",
                     help="comma-separated backends (cpu,cuda,rocm,vulkan), or "
                          "'auto' = cpu + whatever plugin .so files are present")
    ap.add_argument("--n-structures", type=int, default=5,
                     help="benchmark the N largest structures by atom count")
    ap.add_argument("--repeats", type=int, default=3,
                     help="repeats per (structure, threads, backend); median reported")
    ap.add_argument("--opt", action="store_true",
                     help="time -opt instead of -sp (heavier, more realistic workload)")
    ap.add_argument("--timeout", type=int, default=1800)
    args = ap.parse_args()

    if not CURCUMA.exists():
        raise SystemExit(f"curcuma binary not found at {CURCUMA} - build release/ first")

    threads_list = [int(t) for t in args.threads.split(",")]
    if args.gpu == "auto":
        backends = ["cpu"] + detect_gpu_backends()
    else:
        requested = [b.strip() for b in args.gpu.split(",")]
        available = set(["cpu"] + detect_gpu_backends())
        backends = []
        for b in requested:
            if b != "cpu" and b not in available:
                print(f"! backend '{b}' requested but no libcurcuma_{b}.so found "
                      f"next to {CURCUMA} - skipping", flush=True)
                continue
            backends.append(b)
    if not backends:
        raise SystemExit("no usable backend selected")

    structures = pick_structures(args.root, args.n_structures)
    if not structures:
        raise SystemExit(f"no *.xyz found under {args.root} "
                          f"(fetch a testset first: python scripts/fetch_testset.py fetch ...)")

    print(f"structures: {[p.name for _, p in structures]}", flush=True)
    print(f"threads: {threads_list}  backends: {backends}  method: {args.method}  "
          f"mode: {'opt' if args.opt else 'sp'}", flush=True)

    rows = []
    for atoms, xyz in structures:
        for backend in backends:
            for t in threads_list:
                wall, err = median_run(xyz, args.method, t, backend, args.opt,
                                        args.repeats, args.timeout)
                print(f"  {xyz.name:12s} atoms={atoms:4d} backend={backend:6s} "
                      f"threads={t:3d} wall={wall}", flush=True)
                rows.append({"structure": xyz.name, "atoms": atoms, "method": args.method,
                             "backend": backend, "threads": t, "wall_s": wall, "error": err})

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    rundir = args.root / "_run" / f"perf_{stamp}"
    csv_path = rundir / "results.csv"
    md_path = rundir / "summary.md"
    write_csv(rows, csv_path)
    write_summary(rows, structures, threads_list, backends, md_path)
    print(f"\nWrote {csv_path}\nWrote {md_path}")


if __name__ == "__main__":
    main()
