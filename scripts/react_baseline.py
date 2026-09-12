#!/usr/bin/env python3
# Claude Generated (Sep 2026) - rev-gfnff WP0c
"""Re-measure the react-topology numbers of docs/GFNFF_REACT_TOPOLOGY.md.

Every number in that document was recorded before the GFN-FF gradient-unit fix
(CLAUDE.md Known Issue #28: forces were 1/au = 1.8897x too weak). This script
re-runs the documented conditions with the current binary, with fixed seeds, and
records per run: bond formations / breaks / rebuilds, dE_jump statistics, exchange
resolutions, instability, final fragments and the largest distance from the origin.

Runs (doc line references in the docstrings of RUNS):
  R1  H + H formation temperature threshold      (4 H, 2.5 A wall, 5 ps)
  R2  H2 break at break factor 1.45 vs 2.6       (2 H2, 3000 K, 3 ps)
  R3  N2 + 3 H2, 3500 K, 3.5 A, 20 ps, cap+refractory on / off
  R4  N4H4 (2 N2 + 2 H2), 3500 K, 3.2 A, 15 ps, slack radius on / off
  R5  2 N2 + 6 H2, 3000 K, 4.5 A: wall table (5 ps) and container table (20 ps)
  R6  2 H2, 2500 K, 3 ps: no spurious events

Outputs: test_cases/revgfnff/react_baseline/<run>/{input.xyz,cmd.txt,stdout.log,
summary.json} and test_cases/revgfnff/react_baseline/summary.md. Inputs are
generated deterministically (seeded random packing) and written to
test_cases/revgfnff/systems/.

Usage:
    python scripts/react_baseline.py                 # all runs, 8 in parallel
    python scripts/react_baseline.py --only R1 R6 -j 4
"""
import argparse
import json
import math
import random
import re
import shutil
import statistics
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CURCUMA = REPO / "release" / "curcuma"
SYSTEMS = REPO / "test_cases" / "revgfnff" / "systems"
OUT = REPO / "test_cases" / "revgfnff" / "react_baseline"
ANSI = re.compile(r"\x1b\[[0-9;]*m")
RCOV = {"H": 0.32, "N": 0.71}  # Angstrom, Pyykko; used only for the final fragment count

# ------------------------------------------------------------------ inputs


def pack(molecules, radius, seed, min_dist):
    """Random rigid placement of small molecules in a sphere (seeded, rejection sampling)."""
    rng = random.Random(seed)
    atoms = []
    for sym_list, coords in molecules:
        for _ in range(20000):
            # random rotation (uniform via random quaternion) + random centre inside the sphere
            u1, u2, u3 = rng.random(), rng.random(), rng.random()
            q = (math.sqrt(1 - u1) * math.sin(2 * math.pi * u2), math.sqrt(1 - u1) * math.cos(2 * math.pi * u2),
                 math.sqrt(u1) * math.sin(2 * math.pi * u3), math.sqrt(u1) * math.cos(2 * math.pi * u3))
            w, x, y, z = q[3], q[0], q[1], q[2]
            R = [[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                 [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                 [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]]
            while True:
                c = [rng.uniform(-radius, radius) for _ in range(3)]
                if math.sqrt(sum(v * v for v in c)) < radius * 0.7:
                    break
            placed = []
            for sym, (px, py, pz) in zip(sym_list, coords):
                rx = R[0][0] * px + R[0][1] * py + R[0][2] * pz + c[0]
                ry = R[1][0] * px + R[1][1] * py + R[1][2] * pz + c[1]
                rz = R[2][0] * px + R[2][1] * py + R[2][2] * pz + c[2]
                placed.append((sym, (rx, ry, rz)))
            ok = all(math.dist(p[1], a[1]) >= min_dist for p in placed for a in atoms)
            if ok:
                atoms.extend(placed)
                break
        else:
            raise RuntimeError("packing failed")
    return atoms


H2 = (["H", "H"], [(0.0, 0.0, 0.37), (0.0, 0.0, -0.37)])
N2 = (["N", "N"], [(0.0, 0.0, 0.55), (0.0, 0.0, -0.55)])


def write_xyz(path, atoms, comment):
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [str(len(atoms)), comment] + [f"{s} {x:.6f} {y:.6f} {z:.6f}" for s, (x, y, z) in atoms]
    path.write_text("\n".join(lines) + "\n")


def make_inputs():
    files = {}
    files["h4"] = SYSTEMS / "h4_square.xyz"
    write_xyz(files["h4"], [("H", (0.0, 0.0, 0.0)), ("H", (1.3, 0.0, 0.0)), ("H", (0.0, 1.3, 0.0)), ("H", (1.3, 1.3, 0.0))],
              "4 free H atoms on a 1.3 A square (cli_simplemd_14 geometry)")
    files["2h2"] = SYSTEMS / "2h2.xyz"
    write_xyz(files["2h2"], pack([H2, H2], 2.5, 3, 2.6), "2 H2 in a 2.5 A sphere, packed seed 3, min dist 2.6 A")
    files["n2_3h2"] = SYSTEMS / "n2_3h2.xyz"
    write_xyz(files["n2_3h2"], pack([N2, H2, H2, H2], 3.5, 7, 2.0), "N2 + 3 H2 in a 3.5 A sphere, packed seed 7 (rev-gfnff demo system)")
    files["n4h4"] = SYSTEMS / "2n2_2h2.xyz"
    write_xyz(files["n4h4"], pack([N2, N2, H2, H2], 3.2, 7, 2.0), "2 N2 + 2 H2 in a 3.2 A sphere, packed seed 7")
    files["2n2_6h2"] = SYSTEMS / "2n2_6h2.xyz"
    write_xyz(files["2n2_6h2"], pack([N2, N2, H2, H2, H2, H2, H2, H2], 4.5, 7, 2.2), "2 N2 + 6 H2 in a 4.5 A sphere, packed seed 7, min dist 2.2 A")
    return files


# ------------------------------------------------------------------ runs

BASE = ["-method", "gfnff", "-gfnff.topology_mode", "react", "-md.time_step", "0.5", "-md.seed", "42",
        "-md.no_restart", "-md.rattle_12", "false", "-md.thermostat", "csvr", "-md.coupling", "10",
        "-md.wall_type", "spheric", "-threads", "1", "-verbosity", "1", "-no_bmt"]


def run_defs(files):
    runs = []
    for T in (1000, 2000, 3000, 4000, 5000, 6000):
        runs.append((f"R1_h4_T{T}", files["h4"], ["-temperature", str(T), "-maxtime", "5000", "-md.wall_radius", "2.5"],
                     "R1 H+H formation threshold (doc lines 46-50: 'free H atoms recombine at 3000 K'; test 14 needs 6000 K)"))
    # break factor 1.45 alone is reset to the defaults by the form<break guard, so the early-break
    # variant lowers the formation factor with it (1.2/1.45); the default pair is 1.6/2.6
    for tag, ff, bf in (("early_1.45", "1.2", "1.45"), ("default_2.6", "1.6", "2.6")):
        runs.append((f"R2_2h2_break_{tag}", files["2h2"], ["-temperature", "3000", "-maxtime", "3000", "-md.wall_radius", "3.0",
                     "-gfnff.react_bond_form_factor", ff, "-gfnff.react_bond_break_factor", bf],
                     "R2 H2 break dE_jump (doc lines 54-59: +482 kJ/mol at 1.45 vs +21..34 at 2.6)"))
    for tag, extra in (("filters_on", []), ("filters_off", ["-gfnff.react_valence_cap", "false", "-gfnff.react_refractory_scans", "0"])):
        runs.append((f"R3_n2_3h2_{tag}", files["n2_3h2"], ["-temperature", "3500", "-maxtime", "20000", "-md.wall_radius", "3.5"] + extra,
                     "R3 N2+3H2 3500 K / 3.5 A / 20 ps (doc lines 103-108: 301 events without cap+refractory, 35 with, final N2H2 + 4 H)"))
    for tag, extra in (("slack_on", []), ("slack_off", ["-gfnff.react_slack_form_factor", "1.6"])):
        runs.append((f"R4_n4h4_{tag}", files["n4h4"], ["-temperature", "3500", "-maxtime", "15000", "-md.wall_radius", "3.2"] + extra,
                     "R4 N4H4 3500 K / 3.2 A / 15 ps (doc lines 129-135: 278 events / 115 resolutions with slack radius, 736/248 without)"))
    for pot, wt in (("harmonic", "298.15"), ("harmonic", "10000"), ("logfermi", "298.15"), ("logfermi", "10000")):
        runs.append((f"R5_wall5ps_{pot}_{wt}", files["2n2_6h2"], ["-temperature", "3000", "-maxtime", "5000", "-md.wall_radius", "4.5",
                     "-md.wall_potential", pot, "-md.wall_temp", wt], "R5 wall table 5 ps (doc lines 155-163: max|r| 9.70 / 5.20 / 4.54 / 3.73 A)"))
    for pot in ("harmonic", "logfermi", "pbc"):
        runs.append((f"R5_wall20ps_{pot}", files["2n2_6h2"], ["-temperature", "3000", "-maxtime", "20000", "-md.wall_radius", "4.5",
                     "-md.wall_potential", pot], "R5 container table 20 ps (doc lines 176-179: 11.00/177, 8.76/294, 4.91/601)"))
    runs.append(("R6_2h2_2500K", files["2h2"], ["-temperature", "2500", "-maxtime", "3000", "-md.wall_radius", "3.0"],
                 "R6 2 H2 at 2500 K, 3 ps (doc line 278: no spurious events)"))
    return runs


def fragments(atoms):
    n = len(atoms)
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if math.dist(atoms[i][1], atoms[j][1]) < 1.3 * (RCOV[atoms[i][0]] + RCOV[atoms[j][0]]):
                adj[i].append(j)
                adj[j].append(i)
    seen, frags = set(), []
    for i in range(n):
        if i in seen:
            continue
        stack, comp = [i], []
        while stack:
            a = stack.pop()
            if a in seen:
                continue
            seen.add(a)
            comp.append(a)
            stack.extend(adj[a])
        formula = "".join(f"{s}{c if c > 1 else ''}" for s, c in sorted(
            ((s, sum(1 for a in comp if atoms[a][0] == s)) for s in {atoms[a][0] for a in comp}), key=lambda t: (t[0] != "N", t[0])))
        frags.append(formula)
    return sorted(frags)


def last_frame(trj):
    if not trj.exists():
        return None
    lines = trj.read_text().splitlines()
    if not lines:
        return None
    n = int(lines[0].split()[0])
    block = len(lines) // (n + 2)
    start = (block - 1) * (n + 2) + 2
    atoms = []
    for ln in lines[start:start + n]:
        t = ln.split()
        atoms.append((t[0], (float(t[1]), float(t[2]), float(t[3]))))
    return atoms


def run_one(name, xyz, extra, doc):
    d = OUT / name
    d.mkdir(parents=True, exist_ok=True)
    for f in d.glob("*"):
        if f.is_dir():
            shutil.rmtree(f, ignore_errors=True)
        else:
            f.unlink()
    inp = d / "input.xyz"
    inp.write_text(xyz.read_text())
    cmd = [str(CURCUMA), "-md", "input.xyz"] + BASE + extra
    (d / "cmd.txt").write_text(" ".join(cmd) + "\n")
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=d)
    wall = time.time() - t0
    out = ANSI.sub("", proc.stdout + proc.stderr)
    (d / "stdout.log").write_text(out)
    # "dE_jump = n/a (...)" (Sep 12, 2026) and the older "nan" both mean "not measured": skipped here
    jumps = [float(m.group(1)) for m in re.finditer(r"dE_jump = [-+\d.a-z]+ Eh \(([-+\d.]+) kJ/mol\)", out) if "nan" not in m.group(0)]
    summary = {
        "run": name, "doc": doc, "cmd": " ".join(cmd), "exit": proc.returncode, "wall_s": round(wall, 1),
        "formed": len(re.findall(r"REACT bond formed", out)), "broken": len(re.findall(r"REACT bond broken", out)),
        "rebuilds": len(re.findall(r"REACT rebuild #", out)), "exchange_resolved": len(re.findall(r"REACT exchange", out)),
        "dE_jump_kJ": {"n": len(jumps), "median": statistics.median(jumps) if jumps else None,
                       "min": min(jumps) if jumps else None, "max": max(jumps) if jumps else None,
                       "sum": sum(jumps) if jumps else None},
        # only the MD's own instability messages count; the status table prints "-nan" for
        # an undefined average at step 0, which is not an instability
        "nan_or_unstable": bool(re.search(r"Simulation got unstable|dynamics became unstable|NaN/Inf|REACT scan skipped", out)),
    }
    trj = d / "input.trj.xyz"
    if not trj.exists():  # with -no_bmt SimpleMD writes its files into <basename>.snapshots/
        trj = d / "input.snapshots" / "input.trj.xyz"
    frame = last_frame(trj)
    if frame:
        summary["final_fragments"] = fragments(frame)
        summary["max_r_A"] = round(max(math.sqrt(sum(v * v for v in a[1])) for a in frame), 2)
    (d / "summary.json").write_text(json.dumps(summary, indent=1))
    return summary


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*", help="run name prefixes to run (default all)")
    ap.add_argument("-j", type=int, default=8)
    args = ap.parse_args()
    files = make_inputs()
    runs = run_defs(files)
    if args.only:
        runs = [r for r in runs if any(r[0].startswith(p) for p in args.only)]
    print(f"{len(runs)} runs, {args.j} parallel, binary {CURCUMA}")
    with ThreadPoolExecutor(max_workers=args.j) as ex:
        results = list(ex.map(lambda r: run_one(*r), runs))
    lines = ["# React baseline with corrected forces (post Known Issue #28)", "",
             "AI-generated (scripts/react_baseline.py), machine-measured. Same seed (42) for every run; MD with CSVR (10 fs), dt 0.5 fs.", "",
             "| run | formed | broken | rebuilds | exchange | dE_jump n / median / min / max [kJ/mol] | NaN | max r [A] | final fragments | wall [s] |",
             "|---|---:|---:|---:|---:|---|---|---:|---|---:|"]
    for s in results:
        j = s["dE_jump_kJ"]
        jt = f"{j['n']} / {j['median']:+.1f} / {j['min']:+.1f} / {j['max']:+.1f}" if j["n"] else "0"
        lines.append(f"| {s['run']} | {s['formed']} | {s['broken']} | {s['rebuilds']} | {s['exchange_resolved']} | {jt} | "
                     f"{'yes' if s['nan_or_unstable'] else 'no'} | {s.get('max_r_A', '-')} | {' + '.join(s.get('final_fragments', []))} | {s['wall_s']} |")
        print(lines[-1])
    (OUT / "summary.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {OUT / 'summary.md'}")


if __name__ == "__main__":
    main()
