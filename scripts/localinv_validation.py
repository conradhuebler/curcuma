#!/usr/bin/env python3
"""Validation of the ppcg local-inverse preconditioner (-eeq_solver.eeq_ppcg_precond localinv) against
the exact dense solver (-eeq_solver.solve_method cholesky) and against ppcg+diag.

Claude Generated (Oct 2026).  Measurement tool, no network access, nothing is written into the repo.
Every curcuma call runs under a memory cap (systemd-run scope + timeout, or the wrapper named in $SAFE_RUN).

Sub-commands
  run      single points for a structure set and the four solver variants (chol, diag, li, li46), appended to a JSONL file
             sets: gmtkn (every STRIDE-th structure plus all charged), mor41, alpb_small, alpb_big
  grad     analytic gradients (-gradient -dump_gradient) for a selection of structures, all variants vs chol
  opt      -opt runs (variants diag/li alternating) on one structure, step count / final energy / wall time
  md       NVE MD runs (variants diag/li alternating), Etot drift, iteration statistics, fallback counts
  analyze  summary (median / P95 / max / worst 5 / flags) of a run JSONL, written as JSON and printed

Variants (reference = chol for the small sets):
  chol : -eeq_solver.solve_method cholesky
  diag : ppcg, eeq_ppcg_precond diag
  li   : ppcg, localinv, default radii (8 / 10 / rebuild 1.0 Bohr)
  li46 : ppcg, localinv, eeq_localinv_radius_bohr 4, eeq_localinv_ext_bohr 6

Fallbacks are counted from the warnings 'preconditioner rejected' (rejected build) and 'indefinite preconditioner
detected' / 'solve repeated with diagonal Jacobi' (runtime fallback in the CG) in the log (verbosity >= 1).
"""
import argparse
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

CURCUMA = os.environ.get("CURCUMA", "/home/conrad/src/curcuma_branches/master/curcuma/build-native-cuda/curcuma")
SAFE_RUN = os.environ.get("SAFE_RUN", "")
REPO = Path(__file__).resolve().parent.parent
GMTKN = Path(os.environ.get("GMTKN55", "/home/conrad/src/curcuma_branches/curcuma/test_cases/GMTKN55-testset"))
MOR41 = Path(os.environ.get("MOR41", "/home/conrad/src/curcuma_branches/curcuma/test_cases/MOR41-testset"))
STRUCTS = REPO / "test_cases" / "structures"
HARTREE_TO_KCAL = 627.509474

VARIANTS = {
    "chol": ["-eeq_solver.solve_method", "cholesky"],
    "diag": ["-eeq_solver.solve_method", "ppcg", "-eeq_solver.eeq_ppcg_precond", "diag"],
    "li": ["-eeq_solver.solve_method", "ppcg", "-eeq_solver.eeq_ppcg_precond", "localinv"],
    "li46": ["-eeq_solver.solve_method", "ppcg", "-eeq_solver.eeq_ppcg_precond", "localinv",
             "-eeq_solver.eeq_localinv_radius_bohr", "4", "-eeq_solver.eeq_localinv_ext_bohr", "6"],
}
ANSI = re.compile(r"\x1b\[[0-9;]*m")
RE_PPCG = re.compile(r"PPCG N=(\d+) nfrag=(\d+) iters=(\d+) converged=(\d)(?:.*?precond=(\w+).*?discarded=(\d+))?")
RE_SP = re.compile(r"Single Point Energy\s*=\s*(-?\d+\.\d+)")
RE_TOT = re.compile(r"^\[RESULT\]\s+Total\s+(-?\d+\.\d{6,})\s*$", re.M)
RE_REJ = re.compile(r"preconditioner rejected")
RE_RT = re.compile(r"indefinite preconditioner detected|solve repeated with diagonal")


def busy_cores(interval=3.0):
    """Instantaneous number of busy cores from /proc/stat (the 1-min load average lags behind our own runs)."""
    def snap():
        v = [int(x) for x in open("/proc/stat").readline().split()[1:]]
        return sum(v) - v[3] - v[4], sum(v)
    b0, t0 = snap()
    time.sleep(interval)
    b1, t1 = snap()
    return (b1 - b0) / max(t1 - t0, 1) * (os.cpu_count() or 1)


def wait_for_load(limit, maxwait=600):
    """Wait until fewer than `limit` cores are busy (foreign load), at most maxwait seconds."""
    if os.environ.get("LOCALINV_NOWAIT"):
        return
    t0 = time.time()
    while busy_cores() > limit and time.time() - t0 < maxwait:
        time.sleep(7)


def run_cmd(args, workdir, tmax=300, memmax="12G", extra_env=None):
    """Run curcuma through the safe wrapper in workdir; returns dict(rc, wall, out, err) with logs as text."""
    env = dict(os.environ)
    env.update({"TMAX": str(tmax), "MEMMAX": memmax, "CURCUMA_PPCG_ITERS": "1"})
    if extra_env:
        env.update(extra_env)
    # Memory cap is mandatory: an uncapped run of a large charged system once reached 155 GB and took tmux down.
    if SAFE_RUN:
        prefix = [SAFE_RUN]
    elif shutil.which("systemd-run"):
        prefix = ["systemd-run", "--user", "--scope", "-q", "-p", f"MemoryMax={memmax}", "-p", "MemorySwapMax=0",
                  "timeout", str(tmax)]
    else:
        prefix = ["timeout", str(tmax)]
    cmd = prefix + [CURCUMA] + [str(a) for a in args]
    Path(workdir).mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with open(Path(workdir) / "out.txt", "w") as fo, open(Path(workdir) / "err.txt", "w") as fe:
        p = subprocess.run(cmd, cwd=workdir, stdout=fo, stderr=fe, env=env)
    wall = time.time() - t0
    out = ANSI.sub("", (Path(workdir) / "out.txt").read_text(errors="replace"))
    err = ANSI.sub("", (Path(workdir) / "err.txt").read_text(errors="replace"))
    return {"rc": p.returncode, "wall": wall, "out": out, "err": err}


def parse_sp(res):
    out, err = res["out"], res["err"]
    r = {"rc": res["rc"], "wall": round(res["wall"], 2)}
    m = RE_SP.search(out)
    r["e8"] = float(m.group(1)) if m else None
    m = RE_TOT.search(out)
    r["e10"] = float(m.group(1)) if m else None
    pl = [(int(a), int(b), int(c), int(d), e or "diag", int(f or 0)) for a, b, c, d, e, f in RE_PPCG.findall(err + out)]
    r["ppcg_lines"] = len(pl)
    r["iters_last"] = pl[-1][2] if pl else None
    r["iters_sum"] = sum(x[2] for x in pl) if pl else None
    r["iters_max"] = max(x[2] for x in pl) if pl else None
    r["precond_last"] = pl[-1][4] if pl else None
    r["n_localinv_lines"] = sum(1 for x in pl if x[4] == "localinv")
    r["nonconv"] = sum(1 for x in pl if x[3] == 0)
    r["rejected_builds"] = len(RE_REJ.findall(out + err))
    r["runtime_fallbacks"] = len(RE_RT.findall(out + err))
    return r


# ----------------------------------------------------------------- case lists
def gmtkn_cases(stride):
    cases = []
    subsets = sorted(p.name for p in GMTKN.iterdir() if p.is_dir() and not p.name.startswith(("_", ".")))
    allc = []
    for s in subsets:
        for d in sorted((GMTKN / s).iterdir()):
            xyz = d / "struc.xyz"
            if not xyz.exists():
                continue
            chg = int((d / ".CHRG").read_text().strip()) if (d / ".CHRG").exists() else 0
            uhf_t = (d / ".UHF").read_text().strip() if (d / ".UHF").exists() else ""
            uhf = int(uhf_t) if uhf_t else 0
            allc.append((s, d.name, xyz, chg, uhf))
    for i, (s, n, xyz, chg, uhf) in enumerate(allc):
        if i % stride == 0 or chg != 0:
            cases.append({"id": f"{s}/{n}", "xyz": str(xyz), "charge": chg, "spin": uhf,
                          "flags": ["-gfnff.frag_charge_model", "reference"] if chg != 0 else []})
    return cases, len(allc)


def mor41_cases():
    cases = []
    for d in sorted(MOR41.iterdir()):
        if (d / "mol.xyz").exists():
            cases.append({"id": d.name, "xyz": str(d / "mol.xyz"), "charge": 0, "spin": 0, "flags": []})
    return cases


ALPB_FLAGS = ["-gfnff.solvent", "water", "-gfnff.solvent_model", "alpb"]
ALPB_SMALL = ["acetic-acid-dimer", "benzoic-acid", "caffeine.opt", "water8-cluster", "methyl-a-d-galactopyranoside",
              "methyl-a-d-glucopyranoside", "trisaccharide-c18h32o16.gfnff-opt", "aaa-host.rmsd-a"]
ALPB_BIG = ["peo201-water500", "peo201-chain"]


def struct_path(sid):
    m = json.load(open(STRUCTS / "manifest.json"))
    for v in m["structures"]:
        if v["id"] == sid:
            return str(STRUCTS / v["file"]), v
    raise SystemExit(f"structure {sid} not in manifest")


def alpb_cases(ids):
    cases = []
    for sid in ids:
        p, v = struct_path(sid)
        cases.append({"id": sid, "xyz": p, "charge": v.get("charge") or 0, "spin": v.get("spin") or 0,
                      "flags": list(ALPB_FLAGS), "natoms": v.get("natoms")})
    return cases


# ----------------------------------------------------------------- run
def one_sp(case, variant, work, threads, tmax, gradient=False):
    wd = Path(work) / f"{re.sub(r'[^A-Za-z0-9_.-]', '_', case['id'])}__{variant}"
    if wd.exists():
        shutil.rmtree(wd)
    wd.mkdir(parents=True)
    shutil.copy(case["xyz"], wd / "mol.xyz")
    args = ["-sp", "mol.xyz", "-method", "gfnff", "-no_bmt", "-charge", case["charge"], "-spin", case["spin"],
            "-verbosity", 2, "-threads", threads] + VARIANTS[variant] + case["flags"]
    gpath = wd / "cur.grad"
    if gradient:
        args += ["-gradient", "-dump_gradient", str(gpath)]
    res = run_cmd(args, wd, tmax=tmax, memmax=case.get("memmax", "12G"))
    r = parse_sp(res)
    r.update({"id": case["id"], "variant": variant, "charge": case["charge"]})
    r["natoms"] = case.get("natoms") or atom_count(wd / "mol.xyz")
    if gradient and gpath.exists():
        g = []
        for line in gpath.read_text().splitlines():
            if not line.startswith("#"):
                g.append([float(x) for x in line.split()])
        r["grad"] = g
    bad = (r["rc"] != 0 or r["e8"] is None)
    if bad:
        r["tail"] = (res["out"][-800:] + "\n--ERR--\n" + res["err"][-400:])
        # keep the dir for inspection
    else:
        shutil.rmtree(wd, ignore_errors=True)
    return r


def atom_count(xyz):
    try:
        return int(Path(xyz).read_text().split()[0])
    except Exception:
        return None


def cmd_run(a):
    if a.set == "gmtkn":
        cases, total = gmtkn_cases(a.stride)
        print(f"gmtkn: {total} structures found, {len(cases)} selected (stride {a.stride} + charged)")
    elif a.set == "mor41":
        cases = mor41_cases()
    elif a.set == "alpb_small":
        cases = alpb_cases(ALPB_SMALL)
    elif a.set == "alpb_big":
        cases = alpb_cases(ALPB_BIG)
    else:
        raise SystemExit("unknown set")
    if a.limit:
        cases = cases[:a.limit]
    variants = a.variants.split(",")
    done = set()
    outp = Path(a.out)
    if outp.exists():
        for line in outp.read_text().splitlines():
            d = json.loads(line)
            done.add((d["set"], d["id"], d["variant"]))
    jobs = [(c, v) for c in cases for v in variants if (a.set, c["id"], v) not in done]
    print(f"{len(jobs)} runs to do ({len(done)} already in {outp})", flush=True)
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=a.parallel) as ex, open(outp, "a") as fo:
        futs = []
        for c, v in jobs:
            if a.parallel > 1:
                wait_for_load(8)
            futs.append(ex.submit(one_sp, c, v, a.work, a.threads, a.tmax))
            # keep at most 2*parallel jobs queued so load checks stay meaningful
            while sum(1 for f in futs if not f.done()) >= a.parallel:
                time.sleep(0.2)
        n = 0
        for f in as_completed(futs):
            r = f.result()
            r["set"] = a.set
            r.pop("grad", None)
            fo.write(json.dumps(r) + "\n")
            fo.flush()
            n += 1
            if n % 50 == 0:
                print(f"  {n}/{len(jobs)} done, {time.time() - t0:.0f} s", flush=True)
    print(f"finished {len(jobs)} runs in {time.time() - t0:.0f} s")


# ----------------------------------------------------------------- analyze
def pct(vals, p):
    if not vals:
        return None
    s = sorted(vals)
    k = (len(s) - 1) * p / 100.0
    lo, hi = int(k), min(int(k) + 1, len(s) - 1)
    return s[lo] + (s[hi] - s[lo]) * (k - lo)


def energy(r):
    return r["e10"] if r.get("e10") is not None else r.get("e8")


def analyze_file(path, setname, ref_variant=None):
    rows = [json.loads(l) for l in Path(path).read_text().splitlines() if l.strip()]
    rows = [r for r in rows if r["set"] == setname]
    by = {}
    for r in rows:
        by.setdefault(r["id"], {})[r["variant"]] = r
    ref = ref_variant or ("chol" if any("chol" in v for v in by.values()) else "diag")
    out = {"set": setname, "n_structures": len(by), "reference": ref, "variants": {}, "flags": []}
    for var in ("diag", "li", "li46"):
        d = []
        for sid, v in by.items():
            if ref in v and var in v and energy(v[ref]) is not None and energy(v[var]) is not None and v[ref]["rc"] == 0 and v[var]["rc"] == 0:
                d.append((energy(v[var]) - energy(v[ref]), sid, v[var]["natoms"], v[var]["charge"]))
        ad = [abs(x[0]) for x in d]
        out["variants"][var] = {
            "n": len(d), "median_abs": statistics.median(ad) if ad else None, "p95_abs": pct(ad, 95), "max_abs": max(ad) if ad else None,
            "worst5": [{"id": x[1], "natoms": x[2], "charge": x[3], "dE": x[0]} for x in sorted(d, key=lambda t: -abs(t[0]))[:5]],
            "n_fail": sum(1 for v in by.values() if var in v and (v[var]["rc"] != 0 or v[var]["e8"] is None)),
            "n_nonconv": sum(1 for v in by.values() if var in v and v[var]["nonconv"]),
            "rejected_builds": sum(v[var]["rejected_builds"] for v in by.values() if var in v),
            "runtime_fallbacks": sum(v[var]["runtime_fallbacks"] for v in by.values() if var in v),
            "structures_with_fallback": sum(1 for v in by.values() if var in v and (v[var]["rejected_builds"] or v[var]["runtime_fallbacks"])),
            "no_localinv_line": sum(1 for v in by.values() if var in v and var != "diag" and v[var]["n_localinv_lines"] == 0),
        }
        its = [v[var]["iters_last"] for v in by.values() if var in v and v[var]["iters_last"] is not None]
        out["variants"][var]["iters_last_median"] = statistics.median(its) if its else None
        out["variants"][var]["iters_last_max"] = max(its) if its else None
        wl = [v[var]["wall"] for v in by.values() if var in v and v[var]["rc"] == 0]
        out["variants"][var]["wall_median_s"] = statistics.median(wl) if wl else None
    for sid, v in by.items():
        for var, r in v.items():
            if r["rc"] != 0 or r["e8"] is None:
                out["flags"].append({"id": sid, "variant": var, "kind": f"exit {r['rc']} / no energy"})
            if r["nonconv"]:
                out["flags"].append({"id": sid, "variant": var, "kind": "PPCG not converged"})
            if r["rejected_builds"] or r["runtime_fallbacks"]:
                out["flags"].append({"id": sid, "variant": var, "kind": f"jacobi fallback rej={r['rejected_builds']} rt={r['runtime_fallbacks']}"})
        for var in ("li", "li46"):
            if var in v and "diag" in v and ref in v and all(energy(x) is not None for x in (v[var], v["diag"], v[ref])):
                dl = abs(energy(v[var]) - energy(v[ref]))
                dd = abs(energy(v["diag"]) - energy(v[ref]))
                if dl > max(10 * dd, 1e-8):
                    out["flags"].append({"id": sid, "variant": var, "kind": "dE rule", "dE_li": dl, "dE_diag": dd})
            if var in v and "diag" in v and v[var]["iters_last"] is not None and v["diag"]["iters_last"] is not None \
                    and v[var]["iters_last"] > v["diag"]["iters_last"]:
                out["flags"].append({"id": sid, "variant": var, "kind": "more iterations than diag",
                                     "it_li": v[var]["iters_last"], "it_diag": v["diag"]["iters_last"]})
    out["flag_counts"] = {}
    for f in out["flags"]:
        k = (f["variant"], f["kind"].split(" ")[0] + " " + (f["kind"].split(" ")[1] if " " in f["kind"] else ""))
        out["flag_counts"][f"{k[0]}: {k[1]}"] = out["flag_counts"].get(f"{k[0]}: {k[1]}", 0) + 1
    return out


def cmd_analyze(a):
    sets = a.sets.split(",")
    res = {}
    for s in sets:
        r = analyze_file(a.inp, s)
        res[s] = r
        print(f"== {s}: n={r['n_structures']} reference={r['reference']}")
        for var, v in r["variants"].items():
            print(f"  {var:5s} n={v['n']} |dE| med={v['median_abs']} p95={v['p95_abs']} max={v['max_abs']} "
                  f"fail={v['n_fail']} nonconv={v['n_nonconv']} rej={v['rejected_builds']} rt={v['runtime_fallbacks']} "
                  f"it_last med={v['iters_last_median']} max={v['iters_last_max']} wall_med={v['wall_median_s']}")
            for w in v["worst5"]:
                print(f"       {w['id']} N={w['natoms']} q={w['charge']} dE={w['dE']:.3e}")
        print("  flag counts:", r["flag_counts"])
    Path(a.out).write_text(json.dumps(res, indent=1))


# ----------------------------------------------------------------- gradients
def cmd_grad(a):
    sel = json.load(open(a.selection))     # list of {id, xyz, charge, spin, flags}
    outp = {}
    for c in sel:
        row = {}
        for v in ("chol", "diag", "li", "li46"):
            r = one_sp(c, v, a.work, a.threads, a.tmax, gradient=True)
            row[v] = r
        g0 = row["chol"].get("grad")
        rec = {"id": c["id"], "natoms": row["chol"]["natoms"], "rc": {v: row[v]["rc"] for v in row}}
        if g0:
            for v in ("diag", "li", "li46"):
                g = row[v].get("grad")
                if g and len(g) == len(g0):
                    rec["maxdg_" + v] = max(abs(x - y) for p, q in zip(g, g0) for x, y in zip(p, q))
            rec["maxg_chol"] = max(abs(x) for p in g0 for x in p)
        outp[c["id"]] = rec
        print(json.dumps(rec), flush=True)
    Path(a.out).write_text(json.dumps(outp, indent=1))


# ----------------------------------------------------------------- opt / md
def cmd_opt(a):
    p = a.xyz
    rows = []
    for i, v in enumerate(a.variants.split(",")):
        wd = Path(a.work) / f"opt_{Path(p).stem}_{v}_{i}"
        if wd.exists():
            shutil.rmtree(wd)
        wd.mkdir(parents=True)
        shutil.copy(p, wd / "mol.xyz")
        args = ["-opt", "mol.xyz", "-method", "gfnff", "-no_bmt", "-charge", a.charge, "-verbosity", 1, "-threads", a.threads] \
            + VARIANTS[v] + (ALPB_FLAGS if a.alpb else []) + a.extra.split()
        wait_for_load(6)
        res = run_cmd(args, wd, tmax=a.tmax, memmax="16G")
        txt = res["out"] + res["err"]
        pl = [(int(x[2]), int(x[5] or 0)) for x in RE_PPCG.findall(txt)]
        rows_tab = re.findall(r"^\s+(\d+)\s+(-?\d+\.\d+)\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+", res["out"], re.M)
        steps = int(rows_tab[-1][0]) if rows_tab else None
        rec = {"variant": v, "run": i, "rc": res["rc"], "wall": round(res["wall"], 1), "ppcg_solves": len(pl),
               "steps": steps, "final_e": float(rows_tab[-1][1]) if rows_tab else None,
               "written": "Optimized structure written" in res["out"], "converged_heuristic": bool(steps and steps < 5000 and "Optimized structure written" in res["out"]),
               "iters_median": statistics.median([x[0] for x in pl]) if pl else None,
               "iters_max": max([x[0] for x in pl]) if pl else None,
               "rejected_builds": len(RE_REJ.findall(txt)), "runtime_fallbacks": len(RE_RT.findall(txt)),
               "log": str(wd / "out.txt")}
        rows.append(rec)
        print(json.dumps(rec), flush=True)
    Path(a.out).write_text(json.dumps(rows, indent=1))


def cmd_md(a):
    rows = []
    for i, v in enumerate(a.variants.split(",")):
        wd = Path(a.work) / f"md_{Path(a.xyz).stem}_{v}_{i}"
        if wd.exists():
            shutil.rmtree(wd)
        wd.mkdir(parents=True)
        shutil.copy(a.xyz, wd / "mol.xyz")
        args = ["-md", "mol.xyz", "-method", "gfnff", "-no_bmt", "-charge", 0, "-verbosity", 1, "-threads", a.threads,
                "-thermostat", "none", "-time_step", 1, "-temperature", 300, "-max_time", a.maxtime, "-seed", a.seed, "-print", a.print_every] \
            + VARIANTS[v] + a.extra.split()
        wait_for_load(6)
        res = run_cmd(args, wd, tmax=a.tmax, memmax="16G")
        txt = res["out"] + res["err"]
        pl = [int(x[2]) for x in RE_PPCG.findall(txt)]
        tab = re.findall(r"^\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+-?\d+\.\d+\s+(-?\d+\.\d+)\s+-?\d+\.\d+\s+(-?\d+\.\d+)\s+-?\d+\.\d+\s+(-?\d+\.\d+)", res["out"], re.M)
        t = [float(x[0]) for x in tab]
        et = [float(x[3]) for x in tab]
        drift = sd = None
        if len(t) > 3:
            n = len(t)
            mt, me = sum(t) / n, sum(et) / n
            sxx = sum((x - mt) ** 2 for x in t)
            slope = sum((x - mt) * (y - me) for x, y in zip(t, et)) / sxx          # Eh/ps
            resid = [y - (me + slope * (x - mt)) for x, y in zip(t, et)]
            drift = slope
            sd = statistics.pstdev(et)
            sd_resid = statistics.pstdev(resid)
        q = max(len(pl) // 4, 1)
        rec = {"variant": v, "run": i, "rc": res["rc"], "wall": round(res["wall"], 1), "ppcg_solves": len(pl),
               "etot_rows": len(t), "t_last_ps": t[-1] if t else None, "etot_first": et[0] if et else None,
               "T_first": float(tab[0][4]) if tab else None, "etot_last": et[-1] if et else None,
               "etot_drift_Eh_per_ps": drift, "etot_std_Eh": sd, "etot_resid_std_Eh": sd_resid if drift is not None else None,
               "iters_median": statistics.median(pl) if pl else None, "iters_max": max(pl) if pl else None,
               "iters_median_first_quarter": statistics.median(pl[:q]) if pl else None,
               "iters_median_last_quarter": statistics.median(pl[-q:]) if pl else None,
               "rejected_builds": len(RE_REJ.findall(txt)), "runtime_fallbacks": len(RE_RT.findall(txt)),
               "log": str(wd / "out.txt")}
        rows.append(rec)
        print(json.dumps(rec), flush=True)
    Path(a.out).write_text(json.dumps(rows, indent=1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    r = sp.add_parser("run", help="single points of a set")
    r.add_argument("--set", required=True, choices=["gmtkn", "mor41", "alpb_small", "alpb_big"])
    r.add_argument("--work", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--variants", default="chol,diag,li,li46")
    r.add_argument("--stride", type=int, default=8)
    r.add_argument("--limit", type=int, default=0)
    r.add_argument("--parallel", type=int, default=4)
    r.add_argument("--threads", type=int, default=2)
    r.add_argument("--tmax", type=int, default=300)
    r.set_defaults(f=cmd_run)
    an = sp.add_parser("analyze")
    an.add_argument("--inp", required=True)
    an.add_argument("--out", required=True)
    an.add_argument("--sets", default="gmtkn,mor41,alpb_small,alpb_big")
    an.set_defaults(f=cmd_analyze)
    g = sp.add_parser("grad")
    g.add_argument("--selection", required=True, help="JSON list of cases {id,xyz,charge,spin,flags}")
    g.add_argument("--work", required=True)
    g.add_argument("--out", required=True)
    g.add_argument("--threads", type=int, default=2)
    g.add_argument("--tmax", type=int, default=300)
    g.set_defaults(f=cmd_grad)
    o = sp.add_parser("opt")
    o.add_argument("--xyz", required=True)
    o.add_argument("--work", required=True)
    o.add_argument("--out", required=True)
    o.add_argument("--variants", default="diag,li")
    o.add_argument("--charge", type=int, default=0)
    o.add_argument("--alpb", action="store_true")
    o.add_argument("--threads", type=int, default=16)
    o.add_argument("--tmax", type=int, default=300)
    o.add_argument("--extra", default="", help="extra curcuma flags as ONE string (give it with =, e.g. --extra='-flag value')")
    o.set_defaults(f=cmd_opt)
    m = sp.add_parser("md")
    m.add_argument("--xyz", required=True)
    m.add_argument("--work", required=True)
    m.add_argument("--out", required=True)
    m.add_argument("--variants", default="diag,li")
    m.add_argument("--threads", type=int, default=16)
    m.add_argument("--maxtime", type=float, default=2000, help="fs")
    m.add_argument("--seed", type=int, default=1)
    m.add_argument("--print_every", type=int, default=10, help="status table every N steps")
    m.add_argument("--tmax", type=int, default=300)
    m.add_argument("--extra", default="", help="extra curcuma flags as ONE string (give it with =)")
    m.set_defaults(f=cmd_md)
    a = ap.parse_args()
    a.f(a)


if __name__ == "__main__":
    main()
