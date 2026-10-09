#!/usr/bin/env python3
"""EEQ-solver benchmark for GFN-FF: single points vs xtb, then MD timings per solver.

Claude Generated (Oct 2026). Works for any structure(s); nothing here is specific to one molecule.

Structures are named by id (or alias, e.g. `polymer_2x` = `peo201x2-water1500`) from the structure
library test_cases/structures (see its README, `scripts/structlib.py`); a path to an .xyz file works
too. For library structures charge and spin of the manifest are passed on (`-charge/-spin` for
curcuma, `--chrg/--uhf` for xtb) when the manifest records them, and id, class, atom count and
sha256 are stored in results.json.

For every structure, solver variant and device (cpu/cuda) it runs, each in a FRESH directory
(no .topo.json cache, so the topology is rebuilt every time):
  * sp : curcuma single point, 1 thread by default, energy + wall time. With --xtb the reference
         xtb --gfnff --sp is run once per structure (OMP_NUM_THREADS=1, ulimit -s unlimited).
  * md : curcuma MD of --md-fs femtoseconds with per-step phase timings from the
         md_diagnostics_timing block of <basename>.diag.jsonl (median/min/max over steps 1..n,
         step 0 is dropped as warm-up) and the final Epot from the MD table.

Solver variants (flag spelling matters, see below):
  default        no solver flag (= auto: ppcg for N >= 500 atoms, else cholesky)
  auto, cholesky, ppcg, pcg, lu, ldlt, batched   -eeq_solver.solve_method <name>; an explicit name
                 is never overridden by the size threshold (cholesky is dense at any size).
                 lu and ldlt have no GPU implementation yet: with --devices cuda they are recorded as rc != 0.
  `batched` DROPS cross-fragment Coulomb: energies are not comparable, it is flagged in the table.
(-gfnff.solve X is an equivalent alias since Oct 2026; the benchmark uses the eeq_solver spelling.)

Usage:
  scripts/eeq_solver_benchmark.py peo201x2-water1500 [caffeine ...] --curcuma build/curcuma \\
      --solvers default,cholesky_dense,ppcg,lu --devices cpu,cuda --md-fs 50 --threads 1 \\
      --xtb /opt/bin/xtb --out bench_out
Slow solvers (pcg) can be limited with --slow-md-fs; a run that exceeds --timeout is recorded as
timeout with the steps that were written.
"""
import argparse
import json
import math
import os
import re
import shutil
import signal
import statistics as st
import subprocess
import sys
import time
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

VARIANTS = {
    "default": [],
    "cholesky": ["-eeq_solver.solve_method", "cholesky"],
    "ppcg": ["-eeq_solver.solve_method", "ppcg"],
    "auto": ["-eeq_solver.solve_method", "auto"],
    "ldlt": ["-eeq_solver.solve_method", "ldlt"],
    "lu": ["-eeq_solver.solve_method", "lu"],
    "pcg": ["-eeq_solver.solve_method", "pcg"],
    "batched": ["-eeq_solver.solve_method", "batched"],
}
INEXACT = {"batched"}
ANSI = re.compile(r"\x1b\[[0-9;]*m")


def resolve_structure(arg):
    """Library id/alias or file path -> dict(name, path, id, natoms, class, charge, spin, sha256)."""
    p = Path(arg)
    if p.suffix == ".xyz" and p.exists():
        return dict(name=p.stem, path=p.resolve(), id=None, natoms=None, cls=None,
                    charge=None, spin=None, sha256=None)
    sys.path.insert(0, str(SCRIPT_DIR))
    import structlib                                    # scripts/structlib.py, standard library only
    for e in structlib.load_manifest()["structures"]:
        if arg == e["id"] or arg in e.get("aliases", []):
            return dict(name=e["id"], path=Path(structlib.LIB) / e["file"], id=e["id"],
                        natoms=e["natoms"], cls=e["class"], charge=e.get("charge"),
                        spin=e.get("spin"), sha256=e["sha256"])
    sys.exit(f"'{arg}' is neither an existing .xyz file nor a structure id/alias of the library")


def charge_spin_flags(st_):
    """Flags for a structure whose manifest charge/spin is known and non-zero (unknown = None: no flag)."""
    cur, xtb = [], []
    if st_["charge"]:
        cur += ["-charge", str(st_["charge"])]
        xtb += ["--chrg", str(st_["charge"])]
    if st_["spin"]:
        cur += ["-spin", str(st_["spin"])]
        xtb += ["--uhf", str(st_["spin"])]
    return cur, xtb


def run(cmd, cwd, env, timeout, log):
    """Run cmd in its own process group; kill the group on timeout. Returns (rc, wall_s)."""
    t0 = time.time()
    with open(log, "w") as f:
        p = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=f, stderr=subprocess.STDOUT,
                             start_new_session=True)
        try:
            rc = p.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(p.pid, signal.SIGKILL)
            p.wait()
            rc = "timeout"
    return rc, time.time() - t0


def fresh_dir(root, name, xyz, tag):
    """Fresh run directory holding only a copy of the structure as <tag>.xyz (no .topo.json cache)."""
    d = Path(root) / name
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    shutil.copy(xyz, d / (tag + ".xyz"))
    return d


def parse_sp(log):
    m = re.search(r"Single Point Energy\s*=\s*(-?\d+\.\d+)", ANSI.sub("", Path(log).read_text(errors="ignore")))
    return float(m.group(1)) if m else math.nan


def parse_md_epot(log):
    rows = [l.split() for l in ANSI.sub("", Path(log).read_text(errors="ignore")).splitlines()
            if re.match(r"\s+\d+\.\d+\s+-?\d+\.\d+", l)]
    return float(rows[-1][1]) if rows else math.nan


def md_stats(diag):
    try:
        recs = [json.loads(l) for l in open(diag)]
    except OSError:
        return None
    recs = [r for r in recs if "timing_ms" in r][1:]          # drop step 0
    if not recs:
        return None
    out = {"steps": len(recs)}
    for key in ("step_total", "eeq_solve", "ff_total"):
        v = [r["timing_ms"].get(key, 0.0) for r in recs]
        out[key] = {"median": st.median(v), "min": min(v), "max": max(v)}
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("structures", nargs="+", metavar="STRUCTURE",
                    help="library id or alias (test_cases/structures) or path to an .xyz file")
    ap.add_argument("--curcuma", required=True, type=Path)
    ap.add_argument("--xtb", type=Path, help="xtb binary: adds a reference SP per structure")
    ap.add_argument("--solvers", default="default,cholesky,ppcg,pcg,lu,ldlt")
    ap.add_argument("--devices", default="cpu,cuda")
    ap.add_argument("--method", default="gfnff")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--modes", default="sp,md")
    ap.add_argument("--md-fs", type=float, default=50.0)
    ap.add_argument("--slow-md-fs", type=float, default=3.0, help="MD length for pcg")
    ap.add_argument("--dt", type=float, default=1.0)
    ap.add_argument("--temperature", type=float, default=300.0)
    ap.add_argument("--timeout", type=int, default=3600)
    ap.add_argument("--out", type=Path, default=Path("eeq_solver_bench"))
    ap.add_argument("--extra", default="", help="extra curcuma flags, e.g. '-gpu_strict true'")
    a = ap.parse_args()
    a.curcuma = a.curcuma.resolve()
    a.xtb = a.xtb.resolve() if a.xtb else None
    a.out = a.out.resolve()
    a.structures = [resolve_structure(x) for x in a.structures]

    env = dict(os.environ, OMP_NUM_THREADS=str(a.threads), MKL_NUM_THREADS=str(a.threads))
    a.out.mkdir(parents=True, exist_ok=True)
    modes = a.modes.split(",")
    results = []
    repo = Path(__file__).resolve().parent.parent
    commit = subprocess.run(["git", "-C", str(repo), "rev-parse", "--short", "HEAD"],
                            capture_output=True, text=True).stdout.strip()

    for stc in a.structures:
        xyz, tag = stc["path"], stc["name"]
        cflags, xflags = charge_spin_flags(stc)
        meta = {k: stc[k] for k in ("id", "natoms", "cls", "charge", "spin", "sha256")}
        if a.xtb and "sp" in modes:
            d = fresh_dir(a.out, f"{tag}__xtb", xyz, tag)
            sh = f"ulimit -s unlimited; exec {a.xtb} {tag}.xyz --{a.method} --sp {' '.join(xflags)}"
            e2 = dict(env, OMP_STACKSIZE="4G")
            rc, t = run(["bash", "-c", sh], d, e2, a.timeout, d / "out.log")
            m = re.search(r"TOTAL ENERGY\s+(-?\d+\.\d+)", (d / "out.log").read_text(errors="ignore"))
            results.append(dict(mol=tag, meta=meta, mode="sp", prog="xtb", solver="-", device="cpu", rc=rc,
                                wall_s=t, epot=float(m.group(1)) if m else math.nan))
            print(f"[{tag}] xtb sp rc={rc} {t:.1f}s", flush=True)
        for mode in modes:
            for solver in a.solvers.split(","):
                for dev in a.devices.split(","):
                    fl = list(VARIANTS[solver])
                    if dev == "cuda":
                        fl += ["-gpu", "cuda", "-gpu_strict", "true"]
                    base = [str(a.curcuma), "-sp" if mode == "sp" else "-md", tag + ".xyz", "-method", a.method,
                            "-no_bmt", "-threads", str(a.threads)] + cflags
                    md_fs = a.md_fs
                    if mode == "md":
                        if solver == "pcg":
                            md_fs = min(md_fs, a.slow_md_fs)
                        base += ["-max_time", str(md_fs), "-time_step", str(a.dt),
                                 "-temperature", str(a.temperature), "-md_diagnostics", "true",
                                 "-md_diagnostics_timing", "true", "-dump_frequency", "1"]
                    cmd = base + fl + a.extra.split()
                    d = fresh_dir(a.out, f"{tag}__{mode}_{solver}_{dev}", xyz, tag)
                    rc, t = run(cmd, d, env, a.timeout, d / "out.log")
                    r = dict(mol=tag, meta=meta, mode=mode, prog="curcuma", solver=solver, device=dev, rc=rc,
                             wall_s=t, inexact=solver in INEXACT, cmd=" ".join(cmd))
                    if mode == "sp":
                        r["epot"] = parse_sp(d / "out.log")
                    else:
                        r["md_fs"] = md_fs
                        r["epot"] = parse_md_epot(d / "out.log")
                        r["stats"] = md_stats(d / f"{tag}.diag.jsonl")
                    results.append(r)
                    print(f"[{tag}] {mode} {solver:15} {dev:5} rc={rc} {t:8.1f}s", flush=True)
                    (a.out / "results.json").write_text(json.dumps(
                        dict(commit=commit, threads=a.threads, results=results), indent=1))

    # ---- markdown tables -------------------------------------------------------------
    print(f"\ncommit {commit}, threads {a.threads}")
    for stc in a.structures:
        tag = stc["name"]
        sp = [r for r in results if r["mol"] == tag and r["mode"] == "sp"]
        if sp:
            ref = next((r["epot"] for r in sp if r["prog"] == "xtb"), None)
            print(f"\n### {tag}: single point (fresh topology)\n")
            print("| program | solver | device | wall / s | E / Eh | dE vs xtb / Eh |\n|---|---|---|---:|---:|---:|")
            for r in sp:
                de = f"{r['epot'] - ref:+.2e}" if ref is not None and not math.isnan(r["epot"]) else "-"
                note = " (inexact)" if r.get("inexact") else ""
                print(f"| {r['prog']} | {r['solver']}{note} | {r['device']} | {r['wall_s']:.1f} | "
                      f"{r['epot']:.8f} | {de} |")
        md = [r for r in results if r["mol"] == tag and r["mode"] == "md"]
        if md:
            print(f"\n### {tag}: MD per step (median [min-max] ms, steps 1..n)\n")
            print("| solver | device | fs | n | step | EEQ | FF | Epot_end / Eh | rc |\n|---|---|---:|---:|---:|---:|---:|---:|---|")
            for r in md:
                s = r.get("stats")
                note = " (inexact)" if r.get("inexact") else ""
                if s:
                    f = lambda k: f"{s[k]['median']:.1f} [{s[k]['min']:.0f}-{s[k]['max']:.0f}]"
                    print(f"| {r['solver']}{note} | {r['device']} | {r['md_fs']:g} | {s['steps']} | "
                          f"{f('step_total')} | {f('eeq_solve')} | {f('ff_total')} | {r['epot']:.6f} | {r['rc']} |")
                else:
                    print(f"| {r['solver']}{note} | {r['device']} | {r['md_fs']:g} | 0 | - | - | - | - | {r['rc']} |")


if __name__ == "__main__":
    sys.exit(main())
