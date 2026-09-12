#!/usr/bin/env python3
# Claude Generated (Sep 2026) - rev-gfnff stage 1b: dE_jump histogram of revgfnff react MD runs (rev-gfnff stage 1b acceptance metric).
# Usage: python3 scripts/revgfnff_jump_stats.py --tag NAME [--dt 0.25] [--tscale 1.05] [--only h4_2000_rev ...] [--extra "-gfnff.rev_blend false"]
# Output: test_cases/revgfnff/jump_stats/<tag>/summary.md (+ per-run stdout.log, summary.json).
import json, re, statistics, subprocess, sys, time, argparse, shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
REPO = Path(__file__).resolve().parents[1]
CUR = REPO / "release" / "curcuma"
SYS = REPO / "test_cases" / "revgfnff" / "systems"
ANSI = re.compile(r"\x1b\[[0-9;]*m")
BASE = ["-gfnff.topology_mode", "react", "-md.no_restart",
        "-md.rattle_12", "false", "-md.thermostat", "csvr", "-md.coupling", "10", "-md.wall_type", "spheric",
        "-threads", "1", "-verbosity", "3", "-no_bmt"]
TERMS = ["bond", "angle", "tors", "inv", "brep", "nbrep", "coul", "disp", "hb", "over", "batm"]
RUNS = [  # name, xyz, method, T, maxtime fs, wall radius
    ("h4_2000_rev", "h4_square.xyz", "revgfnff", 2000, 5000, 2.5),
    ("h4_2000_gfnff", "h4_square.xyz", "gfnff", 2000, 5000, 2.5),
    ("2h2_3000_rev", "2h2.xyz", "revgfnff", 3000, 3000, 3.0),
    ("n2_3h2_3500_rev", "n2_3h2.xyz", "revgfnff", 3500, 3000, 3.5),
    ("n2_3h2_3500_gfnff", "n2_3h2.xyz", "gfnff", 3500, 3000, 3.5),
]
def run(name, xyz, method, T, tmax, wall, out, extra):
    d = out / name; shutil.rmtree(d, ignore_errors=True); d.mkdir(parents=True)
    (d / "input.xyz").write_text((SYS / xyz).read_text())
    cmd = [str(CUR), "-md", "input.xyz", "-method", method, "-temperature", str(int(T * TSCALE)), "-maxtime", str(tmax),
           "-md.wall_radius", str(wall)] + BASE + extra
    (d / "cmd.txt").write_text(" ".join(cmd) + "\n")
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=d)
    txt = ANSI.sub("", p.stdout + p.stderr); (d / "stdout.log").write_text(txt)
    jumps, terms = [], {t: [] for t in TERMS}
    for m in re.finditer(r"dE_jump = ([-+\d.a-z]+) Eh \(([-+\d.]+) kJ/mol\)", txt):
        if "nan" in m.group(1): jumps.append(float("nan")); continue
        jumps.append(float(m.group(2)))
    for m in re.finditer(r"REACT jump terms \[kJ/mol\]: (.*)", txt):
        for t, v in re.findall(r"(\w+) ([-+\d.]+)", m.group(1)):
            if t in terms: terms[t].append(float(v))
    fin = [j for j in jumps if j == j]; ab = [abs(j) for j in fin]
    rows = [l.split() for l in txt.splitlines() if re.match(r"\s+\d+\.\d+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+", l)]
    tmax = max((float(r[7]) for r in rows), default=float("nan"))
    s = {"name": name, "T_max": tmax, "method": method, "T": T, "maxtime_fs": tmax, "wall_s": round(time.time() - t0, 1),
         "rc": p.returncode, "n_rebuild": len(jumps), "n_nan": len(jumps) - len(fin),
         "formed": len(re.findall(r"REACT bond formed", txt)), "broken": len(re.findall(r"REACT bond broken", txt)),
         "median_abs": statistics.median(ab) if ab else None, "max_abs": max(ab) if ab else None,
         "frac_lt1": (sum(a < 1 for a in ab) / len(ab)) if ab else None,
         "frac_lt5": (sum(a < 5 for a in ab) / len(ab)) if ab else None,
         "term_median_abs": {t: (statistics.median([abs(x) for x in v]) if v else None) for t, v in terms.items()},
         "epot_nan": bool(re.search(r"dE_jump = nan", txt[txt.find("REACT rebuild #2"):])) if "REACT rebuild #2" in txt else False}
    (d / "summary.json").write_text(json.dumps(s, indent=1))
    return s
TSCALE = 1.0
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", required=True); ap.add_argument("--only", nargs="*")
    ap.add_argument("-j", type=int, default=4); ap.add_argument("--extra", default="", help="extra curcuma flags as ONE string"); ap.add_argument("--dt", default="0.5"); ap.add_argument("--seed", default="42"); ap.add_argument("--tscale", type=float, default=1.0)
    a = ap.parse_args(); global TSCALE; TSCALE = a.tscale
    out = REPO / "test_cases" / "revgfnff" / "jump_stats" / a.tag; out.mkdir(parents=True, exist_ok=True)
    runs = [r for r in RUNS if not a.only or any(r[0].startswith(o) for o in a.only)]
    with ThreadPoolExecutor(a.j) as ex:
        res = list(ex.map(lambda r: run(*r, out, ["-md.time_step", a.dt, "-md.seed", a.seed] + a.extra.split()), runs))
    lines = ["| run | rebuilds (form/break) | median|dE| | max|dE| | <1 kJ | <5 kJ | dominant term (median) | NaN | wall s |", "|---|---|---|---|---|---|---|---|---|"]
    for s in res:
        tm = {t: v for t, v in s["term_median_abs"].items() if v is not None}
        dom = max(tm, key=tm.get) if tm else "-"
        f = lambda x, fmt="{:.1f}": (fmt.format(x) if x is not None else "-")
        lines.append(f"| {s['name']} | {s['n_rebuild']} ({s['formed']}/{s['broken']}) | {f(s['median_abs'])} | {f(s['max_abs'])} | "
                     f"{f(s['frac_lt1'], '{:.2f}')} | {f(s['frac_lt5'], '{:.2f}')} | {dom} {f(tm.get(dom))} | {s['n_nan']} | {s['T_max']:.0f} | {s['wall_s']} |")
    (out / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
if __name__ == "__main__": main()
