#!/usr/bin/env python3
"""
perf_table.py - wall-time table of curcuma single points per backend, method and system.

Claude Generated (Oct 2026). One block per machine for docs/PERFORMANCE_TABLE.md; repeat it after
each hardware change or after a change that touches the GPU or SCF paths.

What is measured: the wall time of the whole process `curcuma -sp <system> -method <m> [-gpu <b>]`
(includes start, plugin load and device initialisation), n repeats per cell, backends alternating
inside every repeat so that clock and load drift hit all backends alike, one process at a time.
Every run starts in a fresh temporary directory with a copy of the input (GFN-FF writes
<name>.topo.json next to the input and reuses it). The energy of each run is read from the output;
cells whose energies differ between backends by more than --etol are marked.

Not measured: MD or optimisation steps, memory, power, the first-call cost separately.

Usage
  python3 scripts/perf_table.py [--curcuma BIN] [--repeats 5] [--threads 10]
                                [--backends cpu,rocm] [--systems NAME,...] [--quick]
                                [--timeout 900] [--json out.json] [--note TEXT]
Output: markdown block on stdout.
"""
import argparse
import json
import os
import re
import shutil
import statistics
import sys
import subprocess
import tempfile
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
S = ROOT / "test_cases/structures"

# name, file, atoms, methods. Fixed list: change it only together with the table in the document.
SYSTEMS = [
    ("caffeine", S / "medium/caffeine.xyz", 24, ["gfnff", "gfn2"]),
    ("trisaccharide", S / "large/trisaccharide-c18h32o16.xyz", 66, ["gfnff", "gfn2"]),
    ("aaa-bglc.conf76", S / "clusters/aaa-bglc.conf76.xyz", 117, ["gfnff", "gfn2"]),
    ("macrocycle-bgal", S / "clusters/macrocycle-bgal.xyz", 231, ["gfnff", "gfn2"]),
    ("peo201-chain", S / "bulk/peo201-chain.xyz", 1410, ["gfnff", "gfn2"]),
    ("peo201-water500", S / "bulk/peo201-water500.xyz", 2910, ["gfnff"]),
    ("urea400-water1000", S / "bulk/urea400-water1000.xyz", 6200, ["gfnff"]),
    ("peo201x2-water1500", S / "bulk/peo201x2-water1500.xyz", 7320, ["gfnff"]),
]
QUICK = {"caffeine", "peo201-chain", "peo201x2-water1500"}

ENERGY = re.compile(r"Single Point Energy\s*=\s*(-?\d+\.\d+)\s*Eh")
ANSI = re.compile(r"\x1b\[[0-9;]*m")


def sh(cmd):
    try:
        return subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=20).stdout.strip()
    except Exception:
        return ""


def machine_info(binary):
    cpu = sh("lscpu | awk -F: '/Model name/{gsub(/^ +/,\"\",$2); print $2; exit}'")
    mem = sh("awk '/MemTotal/{printf \"%.0f\", $2/1048576}' /proc/meminfo")
    gpu = sh("rocminfo 2>/dev/null | awk -F: '/Marketing Name/{gsub(/^ +/,\"\",$2); if ($2 !~ /Intel|AMD Ryzen|Core/) {print $2; exit}}'")
    arch = sh("rocminfo 2>/dev/null | awk '/Name: +gfx/{print $2; exit}'")
    nv = sh("nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -n 1")
    return {"cpu": cpu, "ram_gb": mem, "rocm_gpu": f"{gpu} ({arch})" if gpu else "", "cuda_gpu": nv,
            "load": sh("cut -d' ' -f1-3 /proc/loadavg")}


def provenance(binary):
    head = sh(f"git -C {ROOT} rev-parse --short HEAD")
    dirty = sh(f"git -C {ROOT} status --porcelain | wc -l")
    newer = sh(f"find {ROOT}/src -newer {binary} -type f \\( -name '*.cpp' -o -name '*.h' -o -name '*.cu' -o -name '*.hip' -o -name '*.hiph' \\) | wc -l")
    st = os.stat(binary)
    return {"binary": str(binary), "binary_mtime": datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M"),
            "repo_head": head, "repo_dirty_files": dirty, "sources_newer_than_binary": newer}


def run_one(binary, src, method, backend, threads, timeout):
    tmp = Path(tempfile.mkdtemp(prefix="perf_"))
    try:
        inp = tmp / src.name
        shutil.copy(src, inp)
        cmd = [str(binary), "-sp", inp.name, "-method", method, "-threads", str(threads), "-no_bmt"]
        if backend != "cpu":
            cmd += ["-gpu", backend, "-gpu_strict", "true"]   # a silent CPU fallback would be timed as GPU
        if shutil.which("systemd-run"):
            cmd = ["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=24G", "-p", "MemorySwapMax=0"] + cmd
        t0 = time.perf_counter()
        try:
            r = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            return None, None, "timeout"
        dt = time.perf_counter() - t0
        out = ANSI.sub("", r.stdout + r.stderr)
        m = ENERGY.search(out)
        if r.returncode != 0 or not m:
            return None, None, f"exit {r.returncode}" + ("" if m else ", no energy")
        head = re.search(r"Git Commit Hash:\s*(\w+)", out)
        return dt, float(m.group(1)), head.group(1) if head else ""
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--curcuma", default=str(ROOT / "release/curcuma"))
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--threads", type=int, default=min(os.cpu_count() or 1, 16))
    ap.add_argument("--backends", default="cpu,rocm")
    ap.add_argument("--systems", default="")
    ap.add_argument("--quick", action="store_true", help="three systems only (smoke test of the setup)")
    ap.add_argument("--timeout", type=int, default=900)
    ap.add_argument("--etol", type=float, default=1e-6, help="energy tolerance between backends in Eh")
    ap.add_argument("--json", default="")
    ap.add_argument("--note", default="")
    a = ap.parse_args()

    binary = Path(a.curcuma).resolve()
    backends = a.backends.split(",")
    names = set(a.systems.split(",")) if a.systems else (QUICK if a.quick else None)
    systems = [s for s in SYSTEMS if names is None or s[0] in names]
    for s in systems:
        if not s[1].exists():
            raise SystemExit(f"missing input {s[1]}")

    info, prov = machine_info(binary), provenance(binary)
    load0 = info["load"]
    results = []   # (system, atoms, method, backend, [times], [energies], errors)
    binary_hash = ""
    for name, path, atoms, methods in systems:
        for method in methods:
            times = {b: [] for b in backends}
            energies = {b: [] for b in backends}
            errors = {b: [] for b in backends}
            for rep in range(a.repeats):
                order = backends[rep % len(backends):] + backends[:rep % len(backends)]   # alternate
                for b in order:
                    dt, e, extra = run_one(binary, path, method, b, a.threads, a.timeout)
                    if dt is None:
                        errors[b].append(extra)
                    else:
                        times[b].append(dt)
                        energies[b].append(e)
                        binary_hash = binary_hash or extra
            for b in backends:
                results.append((name, atoms, method, b, times[b], energies[b], errors[b]))
            print(f"[done] {name} {method}", flush=True, file=sys.stderr)

    cpu_med = {}
    for name, atoms, method, b, t, e, err in results:
        if b == "cpu" and t:
            cpu_med[(name, method)] = statistics.median(t)

    load1 = sh("cut -d' ' -f1-3 /proc/loadavg")
    print(f"### {datetime.now():%Y-%m-%d}: {info['cpu']}" + (f", {info['rocm_gpu']}" if info['rocm_gpu'] else "") +
          (f", {info['cuda_gpu']}" if info['cuda_gpu'] else ""))
    print()
    print(f"- Binary `{prov['binary']}` (built {prov['binary_mtime']}, reports commit `{binary_hash or 'unknown'}`); "
          f"repository HEAD `{prov['repo_head']}`, {prov['repo_dirty_files']} uncommitted files, "
          f"{prov['sources_newer_than_binary']} source files newer than the binary.")
    print(f"- {a.repeats} repeats per cell, backends alternating, `-threads {a.threads}`, RAM {info['ram_gb']} GB, "
          f"load average at start {load0}, at end {load1}.")
    print("- Wall time of the whole process (start, plugin load, device initialisation included), `-sp`, fresh directory per run, "
          "GPU runs with `-gpu_strict true`.")
    if a.note:
        print(f"- {a.note}")
    print()
    print("| system | atoms | method | backend | median s | min s | max s | n | speed-up vs CPU | note |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for name, atoms, method, b, t, e, err in results:
        if not t:
            print(f"| {name} | {atoms} | {method} | {b} | - | - | - | 0 | - | {', '.join(sorted(set(err))) or 'no run'} |")
            continue
        med = statistics.median(t)
        sp = ""
        if b != "cpu" and (name, method) in cpu_med:
            sp = f"{cpu_med[(name, method)] / med:.2f}"
        notes = []
        if err:
            notes.append(f"{len(err)} failed: {', '.join(sorted(set(err)))}")
        ref = next((x for x in results if x[0] == name and x[2] == method and x[3] == "cpu" and x[5]), None)
        if b != "cpu" and ref is not None:
            d = abs(statistics.median(e) - statistics.median(ref[5]))
            if d > a.etol:
                notes.append(f"energy differs from CPU by {d:.1e} Eh")
        print(f"| {name} | {atoms} | {method} | {b} | {med:.2f} | {min(t):.2f} | {max(t):.2f} | {len(t)} | {sp} | {'; '.join(notes)} |")

    if a.json:
        Path(a.json).write_text(json.dumps({"info": info, "provenance": prov, "binary_commit": binary_hash,
                                            "repeats": a.repeats, "threads": a.threads,
                                            "cells": [{"system": r[0], "atoms": r[1], "method": r[2], "backend": r[3],
                                                       "times": r[4], "energies": r[5], "errors": r[6]} for r in results]},
                                           indent=1))


if __name__ == "__main__":
    main()
