#!/usr/bin/env python3
# Claude Generated (Jul 2026)
"""MOR41: curcuma native (CPU + CUDA GPU) vs gxtb (xtb 6.7.1) for gfnff/gfn1/gfn2.

For every MOR41 structure and each requested method this runs three engines on
the identical mol.xyz (Angstrom, charge 0, singlet) and captures both the total
energy and the full analytic gradient dE/dx [Eh/Bohr]:

    cur_cpu : curcuma -sp mol.xyz -method M              -dump_gradient <f>
    cur_gpu : curcuma -sp mol.xyz -method M -gpu cuda    -dump_gradient <f>
    gxtb    : gxtb    mol.xyz --sp {--gfnff|--gfn1|--gfn2} --grad   (TM gradient file)

Two comparison axes are reported per structure and per method:
  - internal   : cur_cpu vs cur_gpu   (does the CUDA path reproduce the CPU path?)
  - fidelity   : cur_cpu vs gxtb       (implementation fidelity / method divergence)

Both axes report dE [kcal/mol] and the gradient max|component| and RMS deviation
[Eh/Bohr]. Results are cached in _run_cpu_gpu_gxtb/eval.json keyed by
(structure, method), so re-runs and --only subsets never recompute.

Read-only w.r.t. the test set; writes only under _run_cpu_gpu_gxtb/.

gxtb note: this xtb 6.7.1 build honours the *attached* GFN selectors
"--gfn1"/"--gfn2"; the spaced form "--gfn 1" is silently ignored (defaults to
GFN2). Verified against the printed "Hamiltonian GFNx-xTB" line.

Usage:
    python scripts/mor41_cpu_gpu_gxtb_eval.py                  # all methods, all structures
    python scripts/mor41_cpu_gpu_gxtb_eval.py --method gfn2
    python scripts/mor41_cpu_gpu_gxtb_eval.py --only ED03 PR27 --recompute
"""
import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
TESTSET = REPO / "test_cases" / "MOR41-testset"
RUNDIR = TESTSET / "_run_cpu_gpu_gxtb"
CURCUMA = REPO / "release" / "curcuma"
GXTB = Path("/opt/bin/gxtb")

AU2KCAL = 627.509474
BOHR = 1.8897259886  # Bohr per Angstrom

# curcuma -method name -> (curcuma method arg, gxtb attached selector flag, cur_grad_unit)
# cur_grad_unit: units of curcuma's dumped Gradient(). Verified by finite difference
# against curcuma's OWN energy (Jul 2026):
#   gfnff  -> Eh/Bohr   (matches gxtb + own FD to <0.1%)
#   gfn1   -> Eh/Ang    (own FD ratio 1.000; = gxtb after /BOHR)
#   gfn2   -> Eh/Ang    (own FD ratio ~1.0 for H2O; but analytic is ~1% off its own
#                        FD / gxtb for heavier/TM systems -> real gradient residual)
# gxtb TM gradient files are always Eh/Bohr, so curcuma Eh/Ang gradients are divided
# by BOHR before the fidelity comparison.
METHODS = {
    "gfnff": ("gfnff", "--gfnff", "bohr"),
    "gfn1":  ("gfn1",  "--gfn1",  "ang"),
    "gfn2":  ("gfn2",  "--gfn2",  "ang"),
}

# ------------------------------------------------------------------ parsing


def parse_cur_gradient(path):
    """curcuma -dump_gradient file -> (energy[Eh], [[gx,gy,gz],...])."""
    if not path.exists():
        return None, None
    energy = None
    grad = []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            m = re.search(r"energy\s+(-?\d+\.\d+)", line)
            if m:
                energy = float(m.group(1))
            continue
        toks = line.split()
        if len(toks) == 3:
            try:
                grad.append([float(t) for t in toks])
            except ValueError:
                pass
    return energy, (grad or None)


def parse_tm_gradient(path):
    """xtb/gxtb TM `gradient` file -> (energy[Eh], [[gx,gy,gz],...]).

    Layout: $grad / cycle line (SCF energy) / N coordinate rows (x y z Elem) /
    N gradient rows (gx gy gz) / $end. Coordinate rows have 4 tokens (last is a
    letter) and are skipped; gradient rows have 3 numeric tokens."""
    if not path.exists():
        return None, None
    lines = path.read_text().splitlines()
    energy = None
    start = len(lines)
    for i, line in enumerate(lines):
        if "SCF energy" in line:
            m = re.search(r"SCF energy\s*=\s*(-?\d+\.\d+)", line)
            if m:
                energy = float(m.group(1))
            start = i + 1
            break
    grad = []
    for line in lines[start:]:
        s = line.strip()
        if s.startswith("$"):
            break
        toks = s.split()
        if len(toks) == 3:
            try:
                grad.append([float(t.replace("D", "E").replace("d", "e")) for t in toks])
            except ValueError:
                pass
    return energy, (grad or None)


# ------------------------------------------------------------------ engines


def run_curcuma(xyz, method_arg, gpu, grad_path, log_path):
    cmd = [str(CURCUMA), "-sp", str(xyz), "-method", method_arg,
           "-charge", "0", "-verbosity", "0", "-no_bmt",
           "-threads", "4", "-dump_gradient", str(grad_path)]
    if gpu:
        cmd += ["-gpu", "cuda"]
    grad_path.parent.mkdir(parents=True, exist_ok=True)
    if grad_path.exists():
        grad_path.unlink()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    except subprocess.TimeoutExpired:
        log_path.write_text("TIMEOUT (curcuma)\n")
        return None, None
    log_path.write_text(proc.stdout + proc.stderr)
    return parse_cur_gradient(grad_path)


def run_gxtb(xyz, gxtb_flag, workdir, log_path):
    workdir.mkdir(parents=True, exist_ok=True)
    # clean any prior scratch so we never parse a stale gradient file
    for f in workdir.iterdir():
        try:
            f.unlink()
        except OSError:
            pass
    local = workdir / "mol.xyz"
    shutil.copy(xyz, local)
    cmd = [str(GXTB), "mol.xyz", "--sp", gxtb_flag, "--grad", "-c", "0"]
    env = dict(os.environ)
    env["XTBPATH"] = ""
    env["OMP_NUM_THREADS"] = "4"
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600,
                              cwd=str(workdir), env=env)
    except subprocess.TimeoutExpired:
        log_path.write_text("TIMEOUT (gxtb)\n")
        return None, None
    log_path.write_text(proc.stdout + proc.stderr)
    return parse_tm_gradient(workdir / "gradient")


# ------------------------------------------------------------------ metrics


def grad_dev(a, b, rel_floor=1e-3):
    """Deviation metrics between two N x 3 gradient lists (same units).

    Returns (max|abs diff|, RMS abs diff, max relative diff, mean ratio) where the
    relative metrics only use components with |b| > rel_floor (avoids blow-up on
    near-zero components). Returns Nones if shapes mismatch."""
    if a is None or b is None or len(a) != len(b):
        return None, None, None, None
    mx = ss = 0.0
    n = 0
    maxrel = 0.0
    ratios = []
    for ra, rb in zip(a, b):
        for ca, cb in zip(ra, rb):
            d = abs(ca - cb)
            mx = max(mx, d)
            ss += d * d
            n += 1
            if abs(cb) > rel_floor:
                maxrel = max(maxrel, d / abs(cb))
                ratios.append(ca / cb)
    rms = math.sqrt(ss / n) if n else None
    mean_ratio = sum(ratios) / len(ratios) if ratios else None
    return mx, rms, (maxrel if ratios else None), mean_ratio


def gnorm(g):
    if g is None:
        return None
    return math.sqrt(sum(c * c for r in g for c in r))


# ------------------------------------------------------------------ cache


def load_cache():
    p = RUNDIR / "eval.json"
    return json.loads(p.read_text()) if p.exists() else {}


def save_cache(cache):
    RUNDIR.mkdir(parents=True, exist_ok=True)
    (RUNDIR / "eval.json").write_text(json.dumps(cache, indent=1, sort_keys=True))


def evaluate(structure, method, cache, recompute):
    key = f"{structure}|{method}"
    if not recompute and key in cache:
        return cache[key]
    xyz = TESTSET / structure / "mol.xyz"
    if not xyz.exists():
        cache[key] = {"status": "no_xyz"}
        return cache[key]
    logdir = RUNDIR / "logs"
    logdir.mkdir(parents=True, exist_ok=True)
    method_arg, gxtb_flag, cur_unit = METHODS[method]

    e_cpu, g_cpu = run_curcuma(xyz, method_arg, False,
                               RUNDIR / "grad" / f"{structure}.{method}.cpu.grad",
                               logdir / f"{structure}.{method}.cpu.log")
    e_gpu, g_gpu = run_curcuma(xyz, method_arg, True,
                               RUNDIR / "grad" / f"{structure}.{method}.gpu.grad",
                               logdir / f"{structure}.{method}.gpu.log")
    e_ref, g_ref = run_gxtb(xyz, gxtb_flag,
                            RUNDIR / "gxtb" / f"{structure}_{method}",
                            logdir / f"{structure}.{method}.gxtb.log")

    # curcuma-vs-GPU: same units (both curcuma dumps) -> compare directly.
    gmax_cg, grms_cg, _, _ = grad_dev(g_cpu, g_gpu)
    # curcuma-vs-gxtb: convert curcuma gradient to Eh/Bohr (gxtb's unit) first.
    scale = 1.0 if cur_unit == "bohr" else 1.0 / BOHR
    g_cpu_bohr = [[c * scale for c in r] for r in g_cpu] if g_cpu else None
    gmax_cx, grms_cx, gmaxrel_cx, gratio_cx = grad_dev(g_cpu_bohr, g_ref)
    rec = {
        "status": "ok" if (e_cpu is not None and e_gpu is not None and e_ref is not None) else "partial",
        "e_cpu": e_cpu, "e_gpu": e_gpu, "e_gxtb": e_ref,
        "gnorm_cpu": gnorm(g_cpu), "gnorm_gpu": gnorm(g_gpu), "gnorm_gxtb": gnorm(g_ref),
        "natoms": len(g_cpu) if g_cpu else None,
        # internal: CPU vs GPU (curcuma native units)
        "dE_cpu_gpu_kcal": (e_cpu - e_gpu) * AU2KCAL if (e_cpu is not None and e_gpu is not None) else None,
        "gmax_cpu_gpu": gmax_cg, "grms_cpu_gpu": grms_cg,
        # fidelity: CPU vs gxtb (Eh/Bohr, unit-corrected)
        "dE_cpu_gxtb_kcal": (e_cpu - e_ref) * AU2KCAL if (e_cpu is not None and e_ref is not None) else None,
        "gmax_cpu_gxtb": gmax_cx, "grms_cpu_gxtb": grms_cx,
        "gmaxrel_cpu_gxtb": gmaxrel_cx, "gratio_cpu_gxtb": gratio_cx,
    }
    cache[key] = rec
    return rec


# ------------------------------------------------------------------ report


def stats(vals):
    vals = [abs(v) for v in vals if v is not None]
    if not vals:
        return None
    n = len(vals)
    return {"n": n, "mad": sum(vals) / n, "max": max(vals),
            "rms": math.sqrt(sum(v * v for v in vals) / n)}


def signed_stats(vals):
    vals = [v for v in vals if v is not None]
    if not vals:
        return None
    n = len(vals)
    return {"n": n, "md": sum(vals) / n, "mad": sum(abs(v) for v in vals) / n,
            "max": max(abs(v) for v in vals),
            "rms": math.sqrt(sum(v * v for v in vals) / n)}


def write_report(methods, structures, cache):
    RUNDIR.mkdir(parents=True, exist_ok=True)
    path = RUNDIR / "mor41_cpu_gpu_gxtb.md"
    with path.open("w") as f:
        f.write("# MOR41: curcuma CPU vs GPU (CUDA) vs gxtb 6.7.1\n\n")
        f.write("Per-structure single points on identical mol.xyz (charge 0, singlet).\n")
        f.write("Energies [Eh]; dE [kcal/mol]; gradient deviations [Eh/Bohr].\n")
        f.write("- internal = cur_cpu vs cur_gpu; fidelity = cur_cpu vs gxtb.\n\n")
        f.write(f"Reference: gxtb {GXTB} (xtb 6.7.1). Device: CUDA GPU.\n\n")
        for method in methods:
            f.write(f"## {method}\n\n")
            rows = [(s, cache.get(f"{s}|{method}", {})) for s in structures]
            de_cg = [r.get("dE_cpu_gpu_kcal") for _, r in rows]
            gm_cg = [r.get("gmax_cpu_gpu") for _, r in rows]
            de_cx = [r.get("dE_cpu_gxtb_kcal") for _, r in rows]
            gm_cx = [r.get("gmax_cpu_gxtb") for _, r in rows]

            f.write("### Internal: curcuma CPU vs GPU\n\n")
            se = signed_stats(de_cg)
            sg = stats(gm_cg)
            if se:
                f.write(f"- energy dE: n={se['n']}  MD={se['md']:+.2e}  MAD={se['mad']:.2e}  "
                        f"max={se['max']:.2e} kcal/mol\n")
            if sg:
                f.write(f"- gradient max|comp|: n={sg['n']}  MAD={sg['mad']:.2e}  "
                        f"max={sg['max']:.2e} Eh/Bohr\n")
            worst = sorted([r for r in rows if r[1].get("dE_cpu_gpu_kcal") is not None],
                           key=lambda r: -abs(r[1]["dE_cpu_gpu_kcal"]))[:5]
            f.write("- largest |dE(cpu-gpu)|: "
                    + ", ".join(f"{s} {r['dE_cpu_gpu_kcal']:+.2e}" for s, r in worst) + "\n\n")

            f.write("### Fidelity: curcuma CPU vs gxtb 6.7.1\n\n")
            se = signed_stats(de_cx)
            sg = stats(gm_cx)
            gr = signed_stats([r.get("gratio_cpu_gxtb") for _, r in rows])
            gmr = stats([r.get("gmaxrel_cpu_gxtb") for _, r in rows])
            if se:
                f.write(f"- energy dE: n={se['n']}  MD={se['md']:+.4f}  MAD={se['mad']:.4f}  "
                        f"max={se['max']:.4f}  RMS={se['rms']:.4f} kcal/mol\n")
            if sg:
                f.write(f"- gradient max|comp| (unit-corrected Eh/Bohr): n={sg['n']}  "
                        f"MAD={sg['mad']:.2e}  max={sg['max']:.2e} Eh/Bohr\n")
            if gr:
                f.write(f"- gradient component ratio cur/gxtb: mean={gr['md']:.4f} "
                        f"(1.0 = identical); per-structure max relative dev "
                        f"mean={gmr['mad']*100:.2f}%  worst={gmr['max']*100:.2f}%\n")
            within = [abs(v) for v in de_cx if v is not None]
            if within:
                f.write(f"- structures within 0.1 kcal/mol of gxtb: "
                        f"{sum(1 for v in within if v < 0.1)}/{len(within)}; "
                        f"within 1.0: {sum(1 for v in within if v < 1.0)}/{len(within)}\n")
            worst = sorted([r for r in rows if r[1].get("dE_cpu_gxtb_kcal") is not None],
                           key=lambda r: -abs(r[1]["dE_cpu_gxtb_kcal"]))[:10]
            f.write("- largest |dE(cpu-gxtb)|:\n")
            for s, r in worst:
                f.write(f"    {s:8s} dE={r['dE_cpu_gxtb_kcal']:+8.4f} kcal  "
                        f"gmax={r.get('gmax_cpu_gxtb') or float('nan'):.2e}\n")

            partial = [s for s, r in rows if r.get("status") not in ("ok", None)]
            if partial:
                f.write(f"\n- incomplete/failed: {partial}\n")

            # full per-structure table
            f.write("\n<details><summary>per-structure table</summary>\n\n")
            f.write("| structure | nat | E_cpu | E_gpu | E_gxtb | dE(cpu-gpu) | gmax(cpu-gpu) | dE(cpu-gxtb) | gmax(cpu-gxtb) |\n")
            f.write("|---|--:|--:|--:|--:|--:|--:|--:|--:|\n")
            for s, r in rows:
                def fe(x):
                    return f"{x:.6f}" if isinstance(x, (int, float)) else "-"
                def fk(x):
                    return f"{x:+.2e}" if isinstance(x, (int, float)) else "-"
                f.write("| {s} | {nat} | {ec} | {eg} | {ex} | {dcg} | {gcg} | {dcx} | {gcx} |\n".format(
                    s=s, nat=r.get("natoms") or "-",
                    ec=fe(r.get("e_cpu")), eg=fe(r.get("e_gpu")), ex=fe(r.get("e_gxtb")),
                    dcg=fk(r.get("dE_cpu_gpu_kcal")), gcg=fk(r.get("gmax_cpu_gpu")),
                    dcx=fk(r.get("dE_cpu_gxtb_kcal")), gcx=fk(r.get("gmax_cpu_gxtb"))))
            f.write("\n</details>\n\n")
    return path


# ------------------------------------------------------------------ main


def all_structures():
    return sorted(d.name for d in TESTSET.iterdir()
                  if d.is_dir() and (d / "mol.xyz").exists())


def main():
    ap = argparse.ArgumentParser(description="MOR41 curcuma CPU/GPU vs gxtb 6.7.1")
    ap.add_argument("--method", choices=list(METHODS) + ["all"], default="all")
    ap.add_argument("--only", nargs="*", default=None, help="structure names")
    ap.add_argument("--recompute", action="store_true")
    args = ap.parse_args()

    if not CURCUMA.exists():
        sys.exit(f"missing curcuma binary: {CURCUMA}")
    if not GXTB.exists():
        sys.exit(f"missing gxtb binary: {GXTB}")

    methods = list(METHODS) if args.method == "all" else [args.method]
    structures = args.only if args.only else all_structures()
    cache = load_cache()

    total = len(structures) * len(methods)
    done = 0
    for method in methods:
        for s in structures:
            r = evaluate(s, method, cache, args.recompute)
            save_cache(cache)
            done += 1
            print(f"[{done:3d}/{total}] {s:8s} {method:5s} "
                  f"E_cpu={r.get('e_cpu')} dE(cpu-gpu)={r.get('dE_cpu_gpu_kcal')} "
                  f"dE(cpu-gxtb)={r.get('dE_cpu_gxtb_kcal')}", flush=True)

    path = write_report(methods, structures, cache)
    print(f"\nWrote {path}")


if __name__ == "__main__":
    main()
