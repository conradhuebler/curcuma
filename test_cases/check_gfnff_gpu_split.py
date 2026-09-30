#!/usr/bin/env python3
# Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Sep 2026)
"""GFN-FF multi-GPU split of one molecule: results must equal the single-device run.

The split (gpu_split_devices, docs/MULTI_GPU_GAPS.md F-3) spreads the implicit Coulomb
tiles and the projected-PCG EEQ (matrix column blocks + matvec) of one molecule over
several GPUs. Both only change the summation order, so a short NVE MD of polymer_2x.xyz
(7320 atoms, 1502 EEQ fragments -> projected PCG; polymer.xyz is ONE fragment and
takes the dense Cholesky, so it would not exercise the EEQ split) must give the same EEQ charges,
energy terms and gradients as the run with -gpu_split_devices none, to well below the
PCG tolerance (absolute 1e-10 on the residual).

Exit 0 = pass, 1 = fail, 77 = skipped (fewer than two CUDA devices or no CUDA plugin).

Usage: check_gfnff_gpu_split.py <curcuma-binary> <polymer_2x.xyz>
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

TOL_Q = 1e-9      # e
TOL_E = 1e-9      # Eh, per energy term
TOL_G = 1e-8      # gradient components


def visible_gpus():
    try:
        out = subprocess.run(["nvidia-smi", "-L"], capture_output=True, text=True, timeout=30)
        return sum(1 for l in out.stdout.splitlines() if l.startswith("GPU "))
    except Exception:
        return 0


def flat(x):
    if isinstance(x, list):
        o = []
        for y in x:
            o.extend(flat(y))
        return o
    return [x]


def run(binary, xyz, workdir, extra):
    os.makedirs(workdir, exist_ok=True)
    shutil.copy(xyz, workdir)
    base = os.path.splitext(os.path.basename(xyz))[0]
    cmd = [binary, "-md", os.path.basename(xyz), "-method", "gfnff", "-gpu", "cuda",
           "-gpu_strict", "true", "-no_bmt", "-verbosity", "1", "-maxtime", "4", "-dt", "1",
           "-thermostat", "none", "-seed", "7", "-md_diagnostics", "true", "-dump", "1",
           "-gpu_split_min_atoms", "0"] + extra
    r = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        print(f"run failed ({' '.join(extra)}): rc={r.returncode}\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
        return None
    with open(os.path.join(workdir, base + ".diag.jsonl")) as f:
        return [json.loads(l) for l in f]


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    binary, xyz = sys.argv[1], sys.argv[2]
    if visible_gpus() < 2:
        print("SKIP: fewer than two CUDA devices")
        return 77
    plugin = os.path.join(os.path.dirname(os.path.abspath(binary)), "libcurcuma_cuda.so")
    if not os.path.exists(plugin):
        print("SKIP: no CUDA plugin next to the binary")
        return 77
    tmp = tempfile.mkdtemp(prefix="gfnff_split_")
    try:
        split = run(binary, xyz, os.path.join(tmp, "split"), ["-gpu_split_devices", "all"])
        single = run(binary, xyz, os.path.join(tmp, "single"), ["-gpu_split_devices", "none"])
        if split is None or single is None or len(split) != len(single) or not split:
            print("FAIL: missing diagnostics")
            return 1
        worst = {"q": 0.0, "E": 0.0, "g": 0.0}
        for a, b in zip(split, single):
            worst["q"] = max(worst["q"], max(abs(x - y) for x, y in zip(flat(a["charges"]), flat(b["charges"]))))
            worst["E"] = max(worst["E"], max(abs(a["energy"][k] - b["energy"][k]) for k in a["energy"]))
            worst["g"] = max(worst["g"], max(abs(x - y) for x, y in zip(flat(a["gradient_norm"]), flat(b["gradient_norm"]))))
        print(f"{len(split)} MD records: max|dq| {worst['q']:.2e} e, max|dE term| {worst['E']:.2e} Eh, "
              f"max|dgrad| {worst['g']:.2e}")
        ok = worst["q"] <= TOL_Q and worst["E"] <= TOL_E and worst["g"] <= TOL_G
        print("PASS" if ok else "FAIL")
        return 0 if ok else 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
