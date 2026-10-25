#!/usr/bin/env python3
# Claude Generated (Sep 2026)
"""Scan curcuma's performance knobs on THIS machine and report the fastest setting.

Every knob it varies is a normal CLI flag (see docs/GPU_TUNING.md and
docs/SQM_PERFORMANCE.md), so whatever comes out of a run can be pasted onto a
production command line - the script never needs a rebuild and never sets an
environment variable.

Three stages:
  1. baseline        - plain defaults, repeated, to get a reference time and
                       the reference energy
  2. one-at-a-time   - each knob's values, everything else at its default
  3. combination     - the winning value of every knob that beat the baseline,
                       applied together and measured (a combination can be
                       slower than its parts, so this is measured, not assumed)

Every run is checked against the baseline energy. A setting that changes the
energy by more than --energy-tol kcal/mol is reported as SUSPECT and excluded
from the recommendation, however fast it was.

Usage:
    python scripts/tuning_sweep.py STRUCTURE.xyz                      # CPU knobs, gfn2
    python scripts/tuning_sweep.py mol.xyz --method gfnff --repeats 3
    python scripts/tuning_sweep.py big.xyz --gpu cuda                 # adds the GPU knobs
    python scripts/tuning_sweep.py big.xyz --gpu cuda --knobs gpu_eigensolver_devices,threads
    python scripts/tuning_sweep.py big.xyz --gradient --json out.json
    python scripts/tuning_sweep.py big.xyz --method gfnff --gpu cuda --md --json md.json   # MD step time

--md (Oct 2026): time = milliseconds per MD STEP of a thermostatted run (-T 300, fixed seed), taken
from the difference of two runs of different length (--md-steps 100,600) so that start-up, topology
and the CUDA context cancel; the "energy" checked against the baseline is the potential energy of
the last step. The MD knobs add the multi-GPU split ladder (docs/MULTI_GPU.md) and the EEQ and overlap
switches of docs/GPU_TUNING.md; knobs that change the trajectory on purpose (hb_update_force_every)
are not scanned.

On a cluster node run it exactly as the production job runs (same module load,
same CUDA_VISIBLE_DEVICES / SLURM allocation) - several knobs, above all the
multi-GPU ones, depend on what the process actually sees.
"""
import argparse
import json
import os
import platform
import re
import statistics
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_CURCUMA = REPO / "release" / "curcuma"

EH_TO_KCAL = 627.5094740631

ENERGY_RE = re.compile(r"Single Point Energy\s*=\s*(-?\d+\.\d+)")
ITER_RE = re.compile(r"SCF converged in (\d+) iterations")
ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")
# MD table row: time, Epot, ...; the last one is the potential energy of the last printed step
MD_ROW_RE = re.compile(r"^\s*\d+\.\d+\s+(-?\d+\.\d+)\s", re.M)


def cpu_count():
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


def thread_ladder(maxthreads):
    """1,2,4,8,... up to the cores this process may use, plus the core count."""
    ladder, t = [], 1
    while t <= maxthreads:
        ladder.append(t)
        t *= 2
    if maxthreads not in ladder:
        ladder.append(maxthreads)
    return ladder


def cpu_knobs(maxthreads):
    """Knobs that apply to any run. Value lists are ordered; index 0 is the default."""
    return {
        "threads": ("-threads", [str(t) for t in thread_ladder(maxthreads)]),
        "eigensolver_max_threads": ("-eigensolver_max_threads", ["0", "4", "8", "16"]),
        "scf_reduce": ("-scf_reduce", ["auto", "sygst", "trsm"]),
        "scf_mixed_precision": ("-scf_mixed_precision", ["true", "false"]),
        "scf_fp32_threshold": ("-scf_fp32_threshold", ["1e-3", "1e-4", "1e-5"]),
        "eigensolver": ("-eigensolver", ["mkl", "native"]),
        "scf_guess": ("-scf_guess", ["eeq", "h0"]),
    }


def gpu_knobs():
    """Knobs that only mean something with -gpu <backend>. See docs/MULTI_GPU.md."""
    return {
        "gpu_eigensolver_devices": ("-gpu_eigensolver_devices", ["all", "none", "0,1"]),
        "gpu_eigensolver_backend": ("-gpu_eigensolver_backend", ["auto", "mp", "mg"]),
        "gpu_eigensolver_block": ("-gpu_eigensolver_block", ["128", "64", "256", "512"]),
        "gpu_eigensolver_fp32": ("-gpu_eigensolver_fp32", ["true", "false"]),
        "gpu_eigensolver_verify": ("-gpu_eigensolver_verify", ["true", "false"]),
        "gpu_density_devices": ("-gpu_density_devices", ["all", "none"]),
        "gpu_multipole_otf": ("-gpu_multipole_otf", ["auto", "on", "off"]),
        "gpu_sparse_integrals": ("-gpu_sparse_integrals", ["auto", "on", "off"]),
    }


def count_gpus():
    """GPUs this process sees (CUDA_VISIBLE_DEVICES honoured), 0 if nvidia-smi is missing."""
    vis = os.environ.get("CUDA_VISIBLE_DEVICES")
    if vis is not None and vis.strip():
        return len([v for v in vis.split(",") if v.strip()])
    try:
        p = subprocess.run(["nvidia-smi", "-L"], capture_output=True, text=True, timeout=30)
        return len([l for l in p.stdout.splitlines() if l.startswith("GPU ")]) if p.returncode == 0 else 0
    except (OSError, subprocess.TimeoutExpired):
        return 0


def md_gpu_knobs(ngpu):
    """GFN-FF MD on a GPU (Oct 2026). Index 0 of every list is the default.
    gpu_split_devices: auto/none/all and the ladder 0,1 / 0,1,2 / ... of the GPUs seen."""
    ladder = [",".join(str(d) for d in range(k)) for k in range(2, max(ngpu, 1))]
    return {
        "gpu_split_devices": ("-gpu_split_devices", ["auto", "none", "all"] + ladder),
        "gpu_split_min_atoms": ("-gpu_split_min_atoms", ["1000", "500", "2000"]),
        "gpu_split_eeq_min_atoms": ("-gpu_split_eeq_min_atoms", ["4000", "1000", "8000"]),
        "gpu_eeq_mixed_precision": ("-gpu_eeq_mixed_precision", ["true", "false"]),
        "gpu_eeq_mixed_rho": ("-gpu_eeq_mixed_rho", ["1e-5", "1e-4", "1e-6"]),
        "gpu_eeq_overlap": ("-gpu_eeq_overlap", ["true", "false"]),
        "gpu_dispersion_overlap": ("-gpu_dispersion_overlap", ["true", "false"]),
        "gpu_staged_list_upload": ("-gpu_staged_list_upload", ["true", "false"]),
        "gpu_cn_pair_skin_bohr": ("-gpu_cn_pair_skin_bohr", ["1.0", "0.5", "2.0"]),
    }


def gfnff_knobs(gpu):
    """GFN-FF's own knobs (docs/GPU_TUNING.md section 3). The dotted scope is required."""
    k = {"gfnff_coulomb_implicit": ("-gfnff.coulomb_implicit", ["true", "false"])}
    if gpu:
        k.update({
            "gfnff_gpu_coulomb_implicit": ("-gfnff.gpu_coulomb_implicit", ["true", "false"]),
            "gfnff_gpu_disp_pairs_on_device": ("-gfnff.gpu_disp_pairs_on_device", ["false", "true"]),
            "gfnff_eeq_mixed_precision": ("-gfnff.eeq_mixed_precision", ["false", "true"]),
            "gfnff_gpu_block_size": ("-gfnff.gpu_block_size", ["0", "128", "256", "512"]),
        })
    return k


# Knobs that only exist for the native GFN1/GFN2 SCF; dropped for force fields.
SCF_ONLY = {
    "eigensolver_max_threads", "scf_reduce", "scf_mixed_precision",
    "scf_fp32_threshold", "eigensolver", "scf_guess",
    "gpu_eigensolver_devices", "gpu_eigensolver_backend", "gpu_eigensolver_block",
    "gpu_eigensolver_fp32", "gpu_eigensolver_verify", "gpu_density_devices",
    "gpu_multipole_otf", "gpu_sparse_integrals",
}


class Runner:
    def __init__(self, args):
        self.args = args
        self.md = args.md
        if self.md:
            self.md_steps = [int(x) for x in args.md_steps.split(",")]
            if len(self.md_steps) != 2 or self.md_steps[0] >= self.md_steps[1]:
                sys.exit("--md-steps needs two increasing step counts, e.g. 100,600")
        self.base_cmd = [str(args.curcuma),
                         "-opt" if args.opt else "-sp", str(args.structure),
                         "-method", args.method, "-no_bmt"]
        if self.md:
            self.base_cmd[1] = "-md"
        if args.gradient and not args.opt:
            self.base_cmd += ["-dump_gradient", "true"]
        if args.gpu:
            self.base_cmd += ["-gpu", args.gpu]
        if args.charge:
            self.base_cmd += ["-charge", str(args.charge)]
        if args.spin:
            self.base_cmd += ["-spin", str(args.spin)]
        self.base_cmd += args.extra

    def run_md_once(self, flags):
        """ms per MD step from two run lengths; energy = Epot of the last row of the longer run."""
        walls, energy = [], None
        for steps in self.md_steps:
            cmd = self.base_cmd + ["-T", "300", "-dt", "1", "-seed", "7", "-max_time", str(steps),
                                   "-print", str(steps), "-dump", "100000000", "-verbosity", "1"] + flags
            t0 = time.perf_counter()
            try:
                p = subprocess.run(cmd, capture_output=True, text=True,
                                   timeout=self.args.timeout, cwd=self.args.workdir)
            except subprocess.TimeoutExpired:
                return {"ok": False, "why": "timeout", "wall": float(self.args.timeout)}
            walls.append(time.perf_counter() - t0)
            out = ANSI_RE.sub("", p.stdout + p.stderr)
            rows = MD_ROW_RE.findall(out)
            if p.returncode != 0 or not rows:
                why = "exit %d" % p.returncode if p.returncode != 0 else "no MD row in output"
                return {"ok": False, "why": why, "wall": walls[-1], "tail": out.strip()[-400:]}
            energy = float(rows[-1])
        ms = (walls[1] - walls[0]) / (self.md_steps[1] - self.md_steps[0]) * 1000.0
        return {"ok": True, "wall": ms, "energy": energy, "iterations": None}

    def run_once(self, flags):
        if self.md:
            return self.run_md_once(flags)
        cmd = self.base_cmd + flags
        t0 = time.perf_counter()
        try:
            p = subprocess.run(cmd, capture_output=True, text=True,
                               timeout=self.args.timeout, cwd=self.args.workdir)
        except subprocess.TimeoutExpired:
            return {"ok": False, "why": "timeout", "wall": float(self.args.timeout)}
        wall = time.perf_counter() - t0
        out = p.stdout + p.stderr
        m = ENERGY_RE.search(out)
        if p.returncode != 0 or not m:
            why = "exit %d" % p.returncode if p.returncode != 0 else "no energy in output"
            return {"ok": False, "why": why, "wall": wall, "tail": out.strip()[-400:]}
        it = ITER_RE.search(out)
        return {"ok": True, "wall": wall, "energy": float(m.group(1)),
                "iterations": int(it.group(1)) if it else None}

    def measure(self, flags, repeats):
        """Repeat a setting; report the BEST wall time (least disturbed by noise)."""
        runs = [self.run_once(flags) for _ in range(repeats)]
        bad = [r for r in runs if not r["ok"]]
        if bad:
            return {"ok": False, "why": bad[0]["why"], "tail": bad[0].get("tail", ""),
                    "flags": flags}
        walls = [r["wall"] for r in runs]
        return {"ok": True, "flags": flags,
                "wall": min(walls),
                "wall_median": statistics.median(walls),
                "wall_spread": max(walls) - min(walls),
                "energy": runs[0]["energy"],
                "iterations": runs[0]["iterations"]}


def hardware_note(curcuma):
    note = {"host": platform.node(), "cpu_count": cpu_count(),
            "when": datetime.now().isoformat(timespec="seconds")}
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("model name"):
                    note["cpu"] = line.split(":", 1)[1].strip()
                    break
    except OSError:
        pass
    try:
        p = subprocess.run(["nvidia-smi", "--query-gpu=name,memory.total",
                            "--format=csv,noheader"], capture_output=True, text=True,
                           timeout=30)
        if p.returncode == 0 and p.stdout.strip():
            note["gpus"] = [g.strip() for g in p.stdout.strip().split("\n")]
    except (OSError, subprocess.TimeoutExpired):
        pass
    for env in ("CUDA_VISIBLE_DEVICES", "OMP_NUM_THREADS", "SLURM_JOB_ID"):
        if os.environ.get(env):
            note[env] = os.environ[env]
    note["curcuma"] = str(curcuma)
    return note


UNIT = ["s"]   # "ms/step" in --md mode


def fmt_row(label, res, base_wall, base_energy, tol):
    if not res["ok"]:
        return "  %-46s FAILED (%s)" % (label, res["why"])
    speed = base_wall / res["wall"] if res["wall"] > 0 else 0.0
    de = (res["energy"] - base_energy) * EH_TO_KCAL
    mark = "  SUSPECT dE=%+.4f kcal/mol" % de if abs(de) > tol else ""
    iters = "%4d" % res["iterations"] if res["iterations"] else "   -"
    return "  %-46s %8.2f %s  %5.2fx  %s it%s" % (label, res["wall"], UNIT[0], speed, iters, mark)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("structure", type=Path, help="XYZ file to benchmark on")
    ap.add_argument("--curcuma", type=Path, default=DEFAULT_CURCUMA)
    ap.add_argument("--method", default="gfn2", help="gfn2 (default), gfn1, gfnff, ...")
    ap.add_argument("--gpu", default="", help="GPU backend (cuda|rocm|vulkan); adds the GPU knobs")
    ap.add_argument("--opt", action="store_true", help="benchmark -opt instead of -sp")
    ap.add_argument("--gradient", action="store_true", help="single point WITH gradient")
    ap.add_argument("--md", action="store_true",
                    help="benchmark MD (-md, 300 K, seed 7): ms per step from two run lengths")
    ap.add_argument("--md-steps", default="100,600", help="the two MD run lengths in steps (default 100,600)")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--spin", type=int, default=0)
    ap.add_argument("--repeats", type=int, default=2, help="runs per setting (default 2)")
    ap.add_argument("--timeout", type=int, default=3600, help="seconds per run (default 3600)")
    ap.add_argument("--knobs", default="", help="comma list: only scan these knobs")
    ap.add_argument("--skip", default="", help="comma list: skip these knobs")
    ap.add_argument("--max-threads", type=int, default=0, help="upper end of the -threads ladder")
    ap.add_argument("--energy-tol", type=float, default=0.1,
                    help="kcal/mol; a setting deviating more is reported SUSPECT (default 0.1)")
    ap.add_argument("--no-combination", action="store_true", help="skip stage 3")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the settings that would be measured, run nothing")
    ap.add_argument("--json", type=Path, help="write the full result table here")
    ap.add_argument("--workdir", type=Path, default=None,
                    help="run in this directory (default: a scratch dir next to the structure)")
    ap.add_argument("extra", nargs="*", help="extra flags passed to every run, after --")
    args = ap.parse_args()

    if not args.curcuma.exists():
        sys.exit("curcuma binary not found: %s" % args.curcuma)
    if not args.structure.exists():
        sys.exit("structure not found: %s" % args.structure)
    args.structure = args.structure.resolve()
    args.curcuma = args.curcuma.resolve()
    if args.workdir is None:
        # GFN-FF caches its perceived topology next to the structure basename, and
        # -sp writes files into the CWD; keep that out of the repo.
        args.workdir = Path("tuning_sweep_run").resolve()
    args.workdir.mkdir(parents=True, exist_ok=True)

    maxthreads = args.max_threads or cpu_count()
    knobs = dict(cpu_knobs(maxthreads))
    if args.gpu:
        knobs.update(gpu_knobs())
    if args.method.lower().startswith(("gfnff", "uff", "qmdff", "cg")):
        knobs = {k: v for k, v in knobs.items() if k not in SCF_ONLY}
    if args.method.lower().startswith("gfnff"):
        knobs.update(gfnff_knobs(bool(args.gpu)))
        if args.md and args.gpu:
            knobs.update(md_gpu_knobs(count_gpus()))
    if args.md:
        UNIT[0] = "ms/step"
        # host stages (HB/XB detection, integrator) matter in MD, but 1 and 2 threads only waste time on a GPU run
        if "threads" in knobs and args.gpu:
            flag, vals = knobs["threads"]
            keep = [v for v in vals if int(v) >= 4] or vals
            knobs["threads"] = (flag, keep)
    if args.knobs:
        wanted = [k.strip() for k in args.knobs.split(",") if k.strip()]
        unknown = [k for k in wanted if k not in knobs]
        if unknown:
            sys.exit("unknown knob(s): %s\nknown: %s"
                     % (", ".join(unknown), ", ".join(sorted(knobs))))
        knobs = {k: knobs[k] for k in wanted}
    for k in (s.strip() for s in args.skip.split(",")):
        knobs.pop(k, None)

    runner = Runner(args)
    hw = hardware_note(args.curcuma)

    if args.dry_run:
        total = sum(len(v[1]) for v in knobs.values()) * args.repeats + args.repeats
        print("dry run: %d knobs, %d curcuma runs" % (len(knobs), total))
        for name, (flag, values) in knobs.items():
            print("  %-34s %s" % (name, "  ".join("%s %s" % (flag, v) for v in values)))
        return

    print("curcuma tuning sweep")
    print("  structure : %s" % args.structure)
    print("  method    : %s%s" % (args.method, "  (-gpu %s)" % args.gpu if args.gpu else ""))
    print("  mode      : %s%s" % ("md %s steps" % args.md_steps if args.md else "opt" if args.opt else "sp",
                                  " + gradient" if args.gradient and not args.opt and not args.md else ""))
    print("  repeats   : %d   timeout %d s   workdir %s" % (args.repeats, args.timeout, args.workdir))
    print("  host      : %s, %s cores%s" % (hw["host"], hw["cpu_count"],
                                            ", %d GPU(s)" % len(hw["gpus"]) if hw.get("gpus") else ""))
    print()

    if args.md:
        # The first run writes the GFN-FF topology cache next to the structure and warms the file cache and the
        # CUDA kernel cache; its short run is slower than the longer one would suggest, which would understate
        # the baseline ms/step (seen: 14.0 against 19.2 ms/step for the same setting). Discarded.
        print("[0/3] warm-up run (discarded)")
        warm = runner.run_once([])
        if not warm["ok"]:
            print(warm.get("tail", ""))
            sys.exit("warm-up run failed: %s" % warm["why"])
        print()
    print("[1/3] baseline (all defaults)")
    base = runner.measure([], args.repeats)
    if not base["ok"]:
        print(base.get("tail", ""))
        sys.exit("baseline run failed: %s" % base["why"])
    print("  %.2f %s (median %.2f, spread %.2f)   E = %.8f Eh   %s iterations"
          % (base["wall"], UNIT[0], base["wall_median"], base["wall_spread"], base["energy"],
             base["iterations"] if base["iterations"] else "-"))
    print()

    results = {"baseline": base, "knobs": {}}
    n_runs = sum(len(v[1]) for v in knobs.values()) * args.repeats
    run_s = (sum(runner.md_steps) * base["wall_median"] / 1000.0) if args.md else base["wall_median"]
    print("[2/3] one knob at a time (%d knobs, %d runs, roughly %.0f min at the baseline time)"
          % (len(knobs), n_runs, n_runs * run_s / 60.0))
    for name, (flag, values) in knobs.items():
        print(" %s" % name)
        per_knob = []
        for v in values:
            res = runner.measure([flag, v], args.repeats)
            res["value"] = v
            per_knob.append(res)
            print(fmt_row("%s %s" % (flag, v), res, base["wall"], base["energy"], args.energy_tol))
        results["knobs"][name] = per_knob
    print()

    # Winners: faster than the baseline by more than the measured noise, and
    # energy-clean. The spread of the baseline repeats is the noise floor.
    noise = max(base["wall_spread"], 0.02 * base["wall"])
    winners = {}
    for name, (flag, values) in knobs.items():
        ok = [r for r in results["knobs"][name]
              if r["ok"] and abs((r["energy"] - base["energy"]) * EH_TO_KCAL) <= args.energy_tol]
        if not ok:
            continue
        best = min(ok, key=lambda r: r["wall"])
        if best["value"] != values[0] and best["wall"] < base["wall"] - noise:
            winners[name] = (flag, best["value"], best["wall"])

    combo = None
    if winners and not args.no_combination:
        flags = []
        for name, (flag, value, _) in winners.items():
            flags += [flag, value]
        print("[3/3] combination of the %d winning knob(s)" % len(winners))
        print("  %s" % " ".join(flags))
        combo = runner.measure(flags, args.repeats)
        combo["flags"] = flags
        print(fmt_row("combined", combo, base["wall"], base["energy"], args.energy_tol))
        results["combination"] = combo
    elif not args.no_combination:
        print("[3/3] combination: nothing beat the baseline by more than the noise "
              "(%.2f %s); defaults are already the best setting here." % (noise, UNIT[0]))
    print()

    print("RECOMMENDATION")
    if combo and combo["ok"] and combo["wall"] < base["wall"] - noise \
            and abs((combo["energy"] - base["energy"]) * EH_TO_KCAL) <= args.energy_tol:
        best_single = min(winners.values(), key=lambda w: w[2])
        use = combo["flags"] if combo["wall"] <= best_single[2] else [best_single[0], best_single[1]]
        print("  %s" % " ".join(use))
        print("  %.2f %s vs %.2f %s baseline (%.2fx), energy identical to %.1e kcal/mol"
              % (min(combo["wall"], best_single[2]), UNIT[0], base["wall"], UNIT[0],
                 base["wall"] / min(combo["wall"], best_single[2]),
                 abs((combo["energy"] - base["energy"]) * EH_TO_KCAL)))
        results["recommendation"] = use
    elif winners:
        name, (flag, value, wall) = min(winners.items(), key=lambda kv: kv[1][2])
        print("  %s %s" % (flag, value))
        print("  %.2f %s vs %.2f %s baseline (%.2fx). The combination was not faster."
              % (wall, UNIT[0], base["wall"], UNIT[0], base["wall"] / wall))
        results["recommendation"] = [flag, value]
    else:
        print("  keep the defaults - no knob beat them outside the noise on this machine.")
        results["recommendation"] = []
    print()
    print("Caveat: measured for THIS structure, method and machine. A knob that wins")
    print("on a 7000-atom single point can lose on a 200-atom optimisation; re-run for")
    print("the workload you actually submit.")

    if args.json:
        results["hardware"] = hw
        results["setup"] = {"structure": str(args.structure), "method": args.method,
                            "gpu": args.gpu, "opt": args.opt, "gradient": args.gradient,
                            "repeats": args.repeats, "energy_tol": args.energy_tol}
        args.json.write_text(json.dumps(results, indent=2))
        print("\nwrote %s" % args.json)


if __name__ == "__main__":
    main()
