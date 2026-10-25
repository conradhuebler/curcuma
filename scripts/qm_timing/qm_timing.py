#!/usr/bin/env python3
"""Wall-time comparison of curcuma's native QM engine with ORCA 6.1 and PySCF (Claude Generated, Oct 2026).

Runs one job at a time, measures the wall time of the whole program run (start-up, integrals,
SCF, correlation; no internal timer is used for the headline number), parses energies and
appends one JSON line per job to WORK/results.jsonl. Large files go to WORK (default
~/.cache/qmtime, never /tmp: RAM tmpfs); a watchdog kills a job when the home disk has less than
40 GB free, its scratch exceeds 150 GB or MemAvailable drops below 8 GB. Scratch is deleted
after each job.

  qm_timing.py run PROGRAM METHOD VARIANT MOL THREADS     one job
  qm_timing.py queue FILE                                  lines "PROGRAM METHOD VARIANT MOL THREADS";
        jobs already in results.jsonl are skipped; after a timeout/OOM/failure of a series
        (program, method, variant, threads) its larger molecules are skipped (file order = size order)
  qm_timing.py table                                       markdown tables from results.jsonl

PROGRAM/METHOD/VARIANT combinations (see CMD_* below):
  curcuma hf    cpu | cpu_rijk | gpu_direct | gpu_rijk | cpu_rijk_stored | gpu_rijk_stored
  curcuma rimp2 cpu | gpu | gpu_bigri | cpu_bigri ;  lmp2 cpu | gpu ;  ccsd cpu
  orca    hf    nori | rijk ;  rimp2 nori ; lmp2 normalpno | tightpno ; ccsd nori ; dlpno_ccsd normalpno
  pyscf   hf    exact | df ;  rimp2 exact ;  ccsd exact
MOL: a name understood by geometries.py (alkaneN, benzene, anthracene, waterN).
Environment: CURCUMA (binary), ORCA (binary), PYSCF_PY (venv python), QMTIME_WORK, QMTIME_CAP (s, default 1800), QMTIME_PIN_CPU (core for 1-thread runs, default 2),
  QMTIME_CURCUMA_SCF = legacy (default; -qm.scf_threshold 1e-7, the setting of the first study, which also keeps the
  tight density test under the newer `-qm.scf_convergence auto`) | default (no SCF flag: curcuma's own convergence
  default, used for the re-measurement), QMTIME_PYSCF_BLAS = omp_nospin (default; OMP_NUM_THREADS and
  OPENBLAS/MKL_NUM_THREADS = threads, OPENBLAS_THREAD_TIMEOUT=4 so idle BLAS workers sleep instead of spinning) |
  single (OPENBLAS/MKL_NUM_THREADS = 1: one pool, but numpy-BLAS-bound steps such as DF-K run serial) |
  omp (first study: both = threads, the two pools spin at once, up to ~70 runnable threads on 24 cores).
"""
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.expanduser(os.environ.get("QMTIME_WORK", "~/.cache/qmtime"))
CURCUMA = os.environ.get("CURCUMA", os.path.join(REPO, "release/curcuma"))
ORCA = os.environ.get("ORCA", "/home/conrad/Downloads/orca_6_1_0_linux_x86-64_shared_openmpi418_avx2/orca")
PYSCF_PY = os.environ.get("PYSCF_PY", os.path.expanduser("~/.venvs/qmref/bin/python"))
CAP = float(os.environ.get("QMTIME_CAP", "1800"))
RESULTS = os.path.join(WORK, "results.jsonl")
PIN_CPU = os.environ.get("QMTIME_PIN_CPU", "2")
# cost exponent in the number of basis functions, used only to skip a size predicted to exceed the cap
GROWTH = {"hf": 2.5, "rimp2": 4.0, "lmp2": 1.5, "ccsd": 5.5, "dlpno_ccsd": 1.5}

sys.path.insert(0, HERE)
import geometries  # noqa: E402

# ---------------------------------------------------------------- command builders
CUR_SCF_MODE = os.environ.get("QMTIME_CURCUMA_SCF", "legacy")
if CUR_SCF_MODE not in ("legacy", "default"):
    sys.exit("QMTIME_CURCUMA_SCF must be legacy or default")
PYSCF_BLAS = os.environ.get("QMTIME_PYSCF_BLAS", "omp_nospin")
if PYSCF_BLAS not in ("single", "omp", "omp_nospin"):
    sys.exit("QMTIME_PYSCF_BLAS must be single, omp or omp_nospin")
CUR_SCF = ["-qm.scf_threshold", "1e-7"] if CUR_SCF_MODE == "legacy" else []
CUR_COMMON = CUR_SCF + ["-qm.frozen_core", "true", "-verbosity", "2", "-no_bmt",
              "-qm.mp2_max_memory_mb", "150000", "-qm.cc_max_memory_mb", "150000",
              "-qm.lmp2_max_memory_mb", "150000"]
CUR_METHOD = {"hf": "hf", "rimp2": "ri-mp2", "lmp2": "lmp2", "ccsd": "ccsd"}
CUR_VARIANT = {
    "cpu": [],
    "cpu_rijk": ["-qm.ri", "jk"],
    "gpu_direct": ["-gpu", "cuda", "-qm.ri", "off", "-qm.gpu_direct", "on", "-qm.gpu_scf", "on"],
    "gpu_rijk": ["-gpu", "cuda", "-qm.ri", "jk"],
    "gpu": ["-gpu", "cuda"],
    # extra, not default: raise the RI host-memory limit so the tensor is stored (CPU) / streamed (GPU)
    # instead of the default batched recomputation above 16 GB
    "cpu_rijk_stored": ["-qm.ri", "jk", "-qm.ri_max_memory_mb", "200000"],
    "gpu_bigri": ["-gpu", "cuda", "-qm.ri_max_memory_mb", "200000"],  # RI-MP2 needs the stored tensor
    # RI-MP2 C80H162: the auto limit is re-read after the 133 GB tensor exists and then refuses the run (Oct 25 2026)
    "cpu_bigri": ["-qm.ri_max_memory_mb", "200000"],
    "gpu_rijk_stored": ["-gpu", "cuda", "-qm.ri", "jk", "-qm.ri_max_memory_mb", "200000"],  # correlation methods: exact SCF (and RI-MP2) on the device where implemented
}
CUR_CC = ["-qm.cc_threshold", "1e-8", "-qm.cc_residual_threshold", "1e-7"]


def cmd_curcuma(method, variant, xyz, threads, jobdir):
    cmd = [CURCUMA, "-sp", xyz, "-method", CUR_METHOD[method], "-threads", str(threads)]
    cmd += CUR_COMMON + CUR_VARIANT[variant] + (CUR_CC if method == "ccsd" else [])
    return cmd, {}


ORCA_KW = {
    ("hf", "nori"): "! RHF def2-SVP NoRI TightSCF",
    ("hf", "rijk"): "! RHF def2-SVP RIJK def2/JK TightSCF",
    ("rimp2", "nori"): "! RHF RI-MP2 def2-SVP NoRI TightSCF",
    ("lmp2", "normalpno"): "! RHF DLPNO-MP2 def2-SVP NoRI TightSCF NormalPNO",
    ("lmp2", "tightpno"): "! RHF DLPNO-MP2 def2-SVP NoRI TightSCF TightPNO",
    ("ccsd", "nori"): "! RHF CCSD def2-SVP NoRI TightSCF",
    ("dlpno_ccsd", "normalpno"): "! RHF DLPNO-CCSD def2-SVP NoRI TightSCF NormalPNO",
}


def cmd_orca(method, variant, xyz, threads, jobdir):
    kw = ORCA_KW[(method, variant)]
    maxcore = 4000 if threads > 1 else 64000
    lines = [kw]
    if threads > 1:
        lines.append(f"%pal nprocs {threads} end")
    lines.append(f"%maxcore {maxcore}")
    if method in ("rimp2", "lmp2", "dlpno_ccsd"):
        lines.append('%basis AuxC "def2/JK" end')  # same aux basis as curcuma (def2-universal-JKFIT)
    if method == "ccsd":
        lines.append("%mdci STol 1e-7 end")
    shutil.copy(xyz, os.path.join(jobdir, "mol.xyz"))
    lines.append("* xyzfile 0 1 mol.xyz")
    with open(os.path.join(jobdir, "job.inp"), "w") as f:
        f.write("\n".join(lines) + "\n")
    env = {"LD_LIBRARY_PATH": os.path.dirname(ORCA) + ":" + os.environ.get("LD_LIBRARY_PATH", "")}
    return [ORCA, "job.inp"], env


def cmd_pyscf(method, variant, xyz, threads, jobdir):
    # PySCF's C kernels run under OpenMP (OMP_NUM_THREADS); numpy's OpenBLAS is a separate pthread pool. With both
    # set to the thread count (first study, QMTIME_PYSCF_BLAS=omp) the two pools spin at the same time (up to ~70
    # runnable threads on 24 cores); 'single' keeps one pool of `threads` workers.
    # 'omp_nospin' keeps numpy's BLAS parallel (PySCF's DF-K and MP2 contractions run through numpy.dot) but lets
    # the idle OpenBLAS workers sleep after 2^4 cycles instead of spinning (OPENBLAS_THREAD_TIMEOUT).
    blas = "1" if PYSCF_BLAS == "single" else str(threads)
    env = {"OMP_NUM_THREADS": str(threads), "OPENBLAS_NUM_THREADS": blas, "MKL_NUM_THREADS": blas,
           "PYSCF_TMPDIR": jobdir, "TMPDIR": jobdir}
    if PYSCF_BLAS == "omp_nospin":
        env["OPENBLAS_THREAD_TIMEOUT"] = "4"
    return [PYSCF_PY, os.path.join(HERE, "pyscf_job.py"), xyz, method, variant, str(threads)], env


BUILDERS = {"curcuma": cmd_curcuma, "orca": cmd_orca, "pyscf": cmd_pyscf}

# ---------------------------------------------------------------- parsers
ANSI = re.compile(r"\x1b\[[0-9;]*m")


def _f(pat, txt, cast=float, last=True):
    m = re.findall(pat, txt, re.M)
    return cast(m[-1] if last else m[0]) if m else None


def parse_curcuma(txt):
    txt = ANSI.sub("", txt)
    r = {"nbf": _f(r"-> (\d+) basis functions", txt, int),
         "e_hf": _f(r"E\(HF\)\s+=\s+(-?\d+\.\d+)", txt) or _f(r"total HF energy:\s+(-?\d+\.\d+)", txt),
         "scf_iter": _f(r"SCF converged in (\d+) iterations", txt, int),
         "e_corr": _f(r"E\(corr, ri-mp2 \)\s+=\s+(-?\d+\.\d+)", txt)
         or _f(r"E\(LMP2 corr\)\s+=\s+(-?\d+\.\d+)", txt) or _f(r"E\(CCSD corr\)\s+=\s+(-?\d+\.\d+)", txt),
         "e_final": _f(r"Single Point Energy = (-?\d+\.\d+)", txt)}
    if re.search(r"(?i)fall(ing)? ?back|fallback", txt):
        r["fallback_lines"] = [ln.strip() for ln in txt.splitlines() if re.search(r"(?i)fallback|falls back", ln)][:5]
    r["eri_stored"] = "built 4-centre ERI" in txt
    m = re.search(r"SCF loop ([\d.]+) ms", txt)
    r["t_scf_loop"] = float(m.group(1)) / 1000 if m else None
    return r


def parse_orca(txt):
    r = {"nbf": _f(r"Number of basis functions\s+\.\.\.\s+(\d+)", txt, int, last=False),
         "e_hf": _f(r"Total Energy\s+:\s+(-?\d+\.\d+) Eh", txt, last=False),
         "scf_iter": _f(r"SCF CONVERGED AFTER\s+(\d+) CYCLES", txt, int),
         "e_final": _f(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", txt),
         "e_corr": _f(r"RI-MP2 CORRELATION ENERGY:\s+(-?\d+\.\d+)", txt)
         or _f(r"DLPNO-MP2 CORRELATION ENERGY:\s+(-?\d+\.\d+)", txt)
         or _f(r"E\(CORR\)\(corrected\)\s+\.\.\.\s+(-?\d+\.\d+)", txt)
         or _f(r"E\(CORR\)\s+\.\.\.\s+(-?\d+\.\d+)", txt)
         or _f(r"MP2 CORRELATION ENERGY\s+:?\s+(-?\d+\.\d+)", txt),
         "normal_end": "ORCA TERMINATED NORMALLY" in txt}
    if r["e_corr"] is None and r["e_final"] is not None and r["e_hf"] is not None and "MP2" in txt:
        r["e_corr_from_diff"] = r["e_final"] - r["e_hf"]
    return r


def parse_pyscf(txt):
    m = re.search(r"^RESULT_JSON (.*)$", txt, re.M)
    r = json.loads(m.group(1)) if m else {}
    r["scf_iter"] = len(re.findall(r"^cycle= *\d+", txt, re.M)) or None
    return r


PARSERS = {"curcuma": parse_curcuma, "orca": parse_orca, "pyscf": parse_pyscf}

# ---------------------------------------------------------------- execution


def du_bytes(path):
    tot = 0
    for root, _, files in os.walk(path):
        for fn in files:
            try:
                tot += os.lstat(os.path.join(root, fn)).st_size
            except OSError:
                pass
    return tot


def mem_available_gb():
    with open("/proc/meminfo") as f:
        for ln in f:
            if ln.startswith("MemAvailable:"):
                return int(ln.split()[1]) / 1048576
    return 999.0


def procs_running():
    """Instantaneous number of runnable threads on the machine (/proc/loadavg field 4)."""
    with open("/proc/loadavg") as f:
        return int(f.read().split()[3].split("/")[0])


def run_job(program, method, variant, mol, threads):
    geo = os.path.join(WORK, "geo")
    os.makedirs(geo, exist_ok=True)
    xyz = os.path.join(geo, mol + ".xyz")
    if not os.path.exists(xyz):
        geometries.write(geo, mol)
    tag = f"{program}_{method}_{variant}_{mol}_t{threads}"
    jobdir = os.path.join(WORK, "scratch", tag)
    shutil.rmtree(jobdir, ignore_errors=True)
    os.makedirs(jobdir)
    os.makedirs(os.path.join(WORK, "logs"), exist_ok=True)
    log = os.path.join(WORK, "logs", tag + ".log")
    cmd, extra_env = BUILDERS[program](method, variant, xyz, threads, jobdir)
    if threads == 1:  # pin serial runs to one P-core (CPU 2 of the 8 P-cores 0-7) for all programs
        cmd = ["taskset", "-c", PIN_CPU] + cmd
    env = dict(os.environ)
    env.update(extra_env)
    load0 = os.getloadavg()
    reason = None
    peak_scratch = 0
    min_mem = 999.0
    running = []
    t0 = time.perf_counter()
    with open(log, "w") as lf:
        lf.write("# " + " ".join(cmd) + "\n")
        lf.flush()
        p = subprocess.Popen(cmd, cwd=jobdir, env=env, stdout=lf, stderr=subprocess.STDOUT,
                             start_new_session=True)
        while True:
            try:
                p.wait(timeout=5)
                break
            except subprocess.TimeoutExpired:
                pass
            el = time.perf_counter() - t0
            sc = du_bytes(jobdir) if int(el) % 30 < 5 else peak_scratch
            peak_scratch = max(peak_scratch, sc)
            ma = mem_available_gb()
            min_mem = min(min_mem, ma)
            running.append(procs_running())
            free = shutil.disk_usage(os.path.expanduser("~")).free / 1e9
            if el > CAP:
                reason = "timeout"
            elif free < 40 or sc > 150e9:
                reason = "disk"
            elif ma < 8:
                reason = "memory"
            if reason:
                os.killpg(p.pid, signal.SIGKILL)
                p.wait()
                break
    wall = time.perf_counter() - t0
    load1 = os.getloadavg()
    txt = open(log, errors="replace").read()
    if program == "orca":  # keep the ORCA output next to the log, drop the scratch
        pass
    res = {"program": program, "method": method, "variant": variant, "mol": mol, "threads": threads,
           "wall": round(wall, 2), "rc": p.returncode, "killed": reason, "load_before": load0[0],
           "load_after": load1[0], "peak_scratch_gb": round(peak_scratch / 1e9, 2),
           "min_mem_available_gb": round(min_mem, 1), "log": log, "date": time.strftime("%Y-%m-%d %H:%M:%S"),
           "running_max": max(running) if running else None,
           "running_median": sorted(running)[len(running) // 2] if running else None}
    if program == "curcuma":
        res["curcuma_scf"] = CUR_SCF_MODE
    if program == "pyscf":
        res["pyscf_blas"] = PYSCF_BLAS
    try:
        res.update(PARSERS[program](txt))
    except Exception as e:  # parsing must never lose the timing
        res["parse_error"] = str(e)
    ok = p.returncode == 0 and reason is None and res.get("e_hf") is not None
    if method != "hf" and res.get("e_corr") is None:
        ok = False
    res["ok"] = ok
    shutil.rmtree(jobdir, ignore_errors=True)
    with open(RESULTS, "a") as f:
        f.write(json.dumps(res) + "\n")
    print(f"{tag}: wall {wall:.1f} s ok={ok} rc={p.returncode} killed={reason} "
          f"E_HF={res.get('e_hf')} E_corr={res.get('e_corr')} nbf={res.get('nbf')} load={load0[0]:.2f}", flush=True)
    return res


def load_results():
    if not os.path.exists(RESULTS):
        return []
    return [json.loads(ln) for ln in open(RESULTS) if ln.strip()]


def nbf_estimate(mol):
    """def2-SVP spherical: H 5, C/N/O 14 functions."""
    return sum(5 if el == "H" else 14 for el, _ in geometries.build(mol))


def queue(path):
    done = {(r["program"], r["method"], r["variant"], r["mol"], r["threads"]) for r in load_results()}
    failed = {(r["program"], r["method"], r["variant"], r["threads"]) for r in load_results() if not r["ok"]}
    last = {}
    for r in load_results():
        if r["ok"]:
            last[(r["program"], r["method"], r["variant"], r["threads"])] = (r["wall"], nbf_estimate(r["mol"]))
    for ln in open(path):
        ln = ln.split("#")[0].strip()
        if not ln:
            continue
        pr, me, va, mo, th = ln.split()
        th = int(th)
        if (pr, me, va, mo, th) in done:
            continue
        if (pr, me, va, th) in failed:
            print(f"skip {ln} (series failed or exceeded the cap at a smaller size)", flush=True)
            continue
        if (pr, me, va, th) in last:
            w, n0 = last[(pr, me, va, th)]
            pred = w * (nbf_estimate(mo) / n0) ** GROWTH[me]
            if pred > 1.5 * CAP:
                print(f"skip {ln} (predicted {pred:.0f} s from the previous size > 1.5 x cap)", flush=True)
                with open(RESULTS, "a") as f:
                    f.write(json.dumps({"program": pr, "method": me, "variant": va, "mol": mo, "threads": th,
                                        "wall": None, "ok": False, "killed": f"skipped, predicted {pred:.0f} s",
                                        "rc": None, "nbf": nbf_estimate(mo)}) + "\n")
                failed.add((pr, me, va, th))
                continue
        r = run_job(pr, me, va, mo, th)
        if not r["ok"]:
            failed.add((pr, me, va, th))
        else:
            last[(pr, me, va, th)] = (r["wall"], nbf_estimate(mo))


# ---------------------------------------------------------------- tables


def table():
    rs = load_results()
    best = {}
    for r in rs:  # the last run of a job wins (warm repeat)
        best[(r["program"], r["method"], r["variant"], r["mol"], r["threads"])] = r
    for method in ("hf", "rimp2", "lmp2", "ccsd", "dlpno_ccsd"):
        sub = [r for r in best.values() if r["method"] == method]
        if not sub:
            continue
        print(f"\n### {method}\n")
        print("| program | variant | mol | nbf | threads | wall s | E_HF | E_corr | iter | note |")
        print("|---|---|---|---|---|---|---|---|---|---|")
        for r in sorted(sub, key=lambda r: (r.get("nbf") or 0, r["mol"], r["program"], r["variant"], r["threads"])):
            note = "" if r["ok"] else (r.get("killed") or f"rc={r['rc']}")
            ehf = f"{r['e_hf']:.8f}" if r.get("e_hf") is not None else "-"
            ec = f"{r['e_corr']:.8f}" if r.get("e_corr") is not None else "-"
            wall = f"{r['wall']:.1f}" if r.get("wall") is not None else "-"
            print(f"| {r['program']} | {r['variant']} | {r['mol']} | {r.get('nbf') or '-'} | {r['threads']} | "
                  f"{wall} | {ehf} | {ec} | {r.get('scf_iter') or '-'} | {note} |")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    if sys.argv[1] == "run":
        run_job(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], int(sys.argv[6]))
    elif sys.argv[1] == "queue":
        queue(sys.argv[2])
    elif sys.argv[1] == "table":
        table()
    else:
        sys.exit(__doc__)
