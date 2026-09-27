#!/usr/bin/env python3
"""rev-gfnff X2- survey: every 2c-3e diatomic anion of stage 2 against its DLPNO-CCSD(T) curve.

Fresh single points (no topology reuse, no topology cache) on each pair's own reference grid,
energies relative to the model's own separated limit A- + B (for ClF-: Cl- + F, the physical
asymptote the reference uses), each fragment computed with the same binary and flags. The
reference is relative to its own fragment sum (fragment_energies_eh). Same construction as the
I2-/ClF- campaign harness (I2_CLF_STATUS.md section 6), ported here so the survey is reproducible.

Regions of a grid point (reported for the worst point and as separate rms columns):
  compressed  r <  r_min of the reference curve (the repulsive / sigma* wall)
  bonded      r >= r_min and the model still perceives the A-B bond (Bond term != 0)
  asymptotic  the model perceives no A-B bond (Bond term == 0), wherever that happens

Usage:
  revgfnff_x2_survey.py survey BIN OUT.json [pair ...] [--cfg name ...]
  revgfnff_x2_survey.py decomp BIN PAIR CFG r1,r2,...   (term table, kcal/mol, vs asymptote)
Environment: X2S_SCRATCH (default /var/tmp/x2s), X2S_JOBS (default 12).

Claude Generated (Sep 27, 2026). AI-generated, machine-tested only.
"""
import json, math, os, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor

import numpy as np

KCAL = 627.509474
HERE = os.path.dirname(os.path.abspath(__file__))
REFROOT = os.path.join(HERE, "..", "test_cases", "revgfnff", "ref", "E")
TMP = os.environ.get("X2S_SCRATCH", "/var/tmp/x2s")
os.makedirs(TMP, exist_ok=True)

# ---- flag sets (docs/REV_GFNFF_STAGE2.md, "Consolidated state ... (2026-09-26)") -------------
P2P3 = ["-gfnff.rev_charge_model", "sqe", "-gfnff.rev_sqe_phase1", "true",
        "-gfnff.rev_excess_electron", "true"]
HARRIS = ["-gfnff.rev_excess_mode", "harris"]
ENS = ["-gfnff.frag_charge_model", "ensemble", "-gfnff.frag_charge_s_max", "1.2"]
VP = ["-gfnff.rev_sqe_virtual_pairs", "true"]
GP = ["-gfnff.rev_sqe_group_pairs_only", "true"]
EA = ["-gfnff.frag_charge_atomic_ea", "true"]
PI = ["-gfnff.rev_pi_excess_electron", "true"]
REC = P2P3 + HARRIS + ENS + VP + GP  # the consolidated recommended halogen setting

# pair -> (A, B, reference dir, pair-specific additions to REC, reference r to exclude + why)
PAIRS = {
    "Cl2-": ("Cl", "Cl", "cl2m_Cl-Cl-_dlpno_ccsdt", [], {}),
    "F2-": ("F", "F", "f2m_F-F-_dlpno_ccsdt", [], {}),
    "Br2-": ("Br", "Br", "br2m_Br-Br-_dlpno_ccsdt", [], {}),
    "I2-": ("I", "I", "i2m_I-I-_dlpno_ccsdt", [], {}),
    "ClF-": ("Cl", "F", "clfm_Cl-F-_dlpno_ccsdt", EA, {}),
    # PI_STAR_STATUS 11.3: the 7.5 A point is in the wrong electronic state (+61.9 above O + O-)
    "O2-": ("O", "O", "o2m_O-O-_dlpno_ccsdt", PI, {7.5: "wrong electronic state (PI_STAR 11.3)"}),
    "S2-": ("S", "S", "s2m_S-S-_dlpno_ccsdt", PI, {}),
}


def configs(pair):
    """name -> (method, flags). 'rec' = the current recommended setting for this pair;
    'rec_noGP' = the same without group_pairs_only (the setting several earlier numbers used)."""
    add = PAIRS[pair][3]
    return {
        "gfnff": ("gfnff", []),
        "rev_default": ("revgfnff", []),
        "rec_noGP": ("revgfnff", P2P3 + HARRIS + ENS + VP + add),
        "rec": ("revgfnff", REC + add),
        # P-A (section 8 of the status file): rec + the 2c-3e candidate bond extension
        "recX": ("revgfnff", REC + add + ["-gfnff.rev_excess_bond_extend", os.environ.get("X2S_EXT", "1.8")]),
        "flat100": ("revgfnff", P2P3 + add),
        "flat100_raw": ("revgfnff", P2P3),               # the state the half rows were fitted on
        "harris_raw": ("revgfnff", P2P3 + HARRIS),        # harris, no window / carrier flags
    }


def ref_points(pair):
    A, B, d, _, excl = PAIRS[pair]
    R = json.load(open(os.path.join(REFROOT, d, "energies.json")))
    fe = R["fragment_energies_eh"]
    e0 = fe["asymptote_eh"] if "asymptote_eh" in fe else fe["atom"] + fe["anion"]
    out = []
    for p in R["points"]:
        r = float(p["label"].split("=")[1])
        if p["energy_eh"] is None or any(abs(r - x) < 1e-6 for x in excl):
            continue
        out.append((r, (p["energy_eh"] - e0) * KCAL))
    return out


def ref_interp(pair, r):
    """natural cubic spline through the reference points (x2lib.ref_interp); nan outside."""
    pts = ref_points(pair)
    x = np.array([p[0] for p in pts]); y = np.array([p[1] for p in pts])
    if r < x[0] - 1e-9 or r > x[-1] + 1e-9:
        return float("nan")
    n = len(x); h = np.diff(x)
    A = np.zeros((n, n)); b = np.zeros(n); A[0, 0] = A[-1, -1] = 1
    for i in range(1, n - 1):
        A[i, i - 1] = h[i - 1]; A[i, i] = 2 * (h[i - 1] + h[i]); A[i, i + 1] = h[i]
        b[i] = 3 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1])
    c = np.linalg.solve(A, b)
    i = min(max(np.searchsorted(x, r) - 1, 0), n - 2)
    t = r - x[i]
    bb = (y[i + 1] - y[i]) / h[i] - h[i] * (2 * c[i] + c[i + 1]) / 3
    d = (c[i + 1] - c[i]) / (3 * h[i])
    return float(y[i] + bb * t + c[i] * t * t + d * t ** 3)


def batch(binp, frames, charge, method, flags, env=None, reuse=False):
    d = tempfile.mkdtemp(dir=TMP)
    with open(os.path.join(d, "s.xyz"), "w") as f:
        for fr in frames:
            f.write(f"{len(fr)}\n\n")
            for a in fr:
                f.write(f"{a[0]} {a[1]:.10f} {a[2]:.10f} {a[3]:.10f}\n")
    cmd = [os.path.abspath(binp), "-sp", "s.xyz", "-method", method, "-batch", "true",
           "-batch_out", "o.jsonl", "-gfnff.cache_topology", "false", "-charge", str(charge),
           "-threads", "1", "-verbosity", "0", "-no_bmt"] + (
               ["-batch_reuse_topology", "true", "-gfnff.topology_mode", "react"] if reuse == "react" or reuse is True
               else ["-batch_reuse_topology", "true"] if reuse == "default"
               else ["-batch_reuse_topology", "false"]) + flags
    e = dict(os.environ)
    if env:
        e.update(env)
    p = subprocess.run(cmd, cwd=d, capture_output=True, text=True, timeout=3600, env=e)
    o = os.path.join(d, "o.jsonl")
    recs = [json.loads(l) for l in open(o)] if os.path.exists(o) else []
    if len(recs) != len(frames):
        raise RuntimeError(f"{len(recs)}/{len(frames)} records in {d}: {p.stderr[-400:]}")
    return recs


def curve(binp, pair, cfg, grid=None, react=False, reuse=None):
    A, B = PAIRS[pair][:2]
    method, flags = configs(pair)[cfg]
    refs = ref_points(pair)
    if grid is None:
        grid = [r for r, _ in refs]
    recs = batch(binp, [[(A, 0, 0, 0), (B, 0, 0, r)] for r in grid], -1, method, flags, reuse=(reuse or react))
    fa = batch(binp, [[(A, 0, 0, 0)]], -1, method, flags)[0]   # A- (Cl- for ClF-)
    fb = batch(binp, [[(B, 0, 0, 0)]], 0, method, flags)[0]    # B
    e0 = fa["energy_eh"] + fb["energy_eh"]
    t0 = {k: fa["terms"].get(k, 0) + fb["terms"].get(k, 0) for k in set(fa["terms"]) | set(fb["terms"])}
    rmap = dict(refs)
    rows = []
    for r, rec in zip(grid, recs):
        E = (rec["energy_eh"] - e0) * KCAL
        terms = {k: (v - t0.get(k, 0)) * KCAL for k, v in rec["terms"].items()}
        rows.append(dict(r=r, E=E, ref=(ref_interp(pair, r) if react else rmap.get(r, float("nan"))), terms=terms, q=rec.get("charges")))
    return rows


def rms(v):
    v = [x for x in v if not math.isnan(x)]
    return math.sqrt(sum(x * x for x in v) / len(v)) if v else float("nan")


def stats(pair, rows):
    rmin = min(ref_points(pair), key=lambda t: t[1])[0]
    d = []
    for x in rows:
        if math.isnan(x["ref"]):
            continue
        bonded = abs(x["terms"].get("Bond", 0.0)) > 1e-9
        reg = "compressed" if x["r"] < rmin - 1e-6 else ("bonded" if bonded else "asymptotic")
        d.append((x["r"], x["E"] - x["ref"], reg, x["E"], x["ref"]))
    w = max(d, key=lambda t: abs(t[1]))
    m = min(d, key=lambda t: t[3])
    return dict(
        n=len(d), full=rms([t[1] for t in d]),
        bonded_all=rms([t[1] for t in d if t[2] != "asymptotic"]),  # every point with a perceived bond
        compressed=rms([t[1] for t in d if t[2] == "compressed"]),
        near=rms([t[1] for t in d if t[2] == "bonded"]),
        asym=rms([t[1] for t in d if t[2] == "asymptotic"]),
        n_comp=sum(t[2] == "compressed" for t in d), n_bond=sum(t[2] != "asymptotic" for t in d),
        worst=dict(r=w[0], dev=w[1], region=w[2], E=w[3], ref=w[4]),
        model_min=dict(r=m[0], E=m[3]), ref_min_r=rmin,
        r_max=max(t[0] for t in d), r_min_grid=min(t[0] for t in d),
    )


def cmd_survey(argv):
    binp, out = argv[0], argv[1]
    rest = argv[2:]
    cfgs = None
    if "--cfg" in rest:
        i = rest.index("--cfg")
        cfgs = rest[i + 1:]
        rest = rest[:i]
    pairs = rest or list(PAIRS)
    cfgs = cfgs or ["gfnff", "rev_default", "rec_noGP", "rec"]
    jobs = {}
    with ThreadPoolExecutor(int(os.environ.get("X2S_JOBS", "12"))) as ex:
        for p in pairs:
            for c in cfgs:
                jobs[(p, c)] = ex.submit(curve, binp, p, c)
        res = {f"{p}|{c}": f.result() for (p, c), f in jobs.items()}
    summ = {k: stats(k.split("|")[0], v) for k, v in res.items()}
    json.dump({"rows": res, "stats": summ}, open(out, "w"), indent=1)
    print(f"{'pair|config':20s} {'n':>3s} {'full':>7s} {'bond+c':>7s} {'compr':>7s} {'n_c':>3s} "
          f"{'near':>7s} {'asym':>7s}  worst dev @ r (region)          model min / ref r_min")
    for k, s in summ.items():
        w = s["worst"]
        print(f"{k:20s} {s['n']:3d} {s['full']:7.2f} {s['bonded_all']:7.2f} {s['compressed']:7.2f} "
              f"{s['n_comp']:3d} {s['near']:7.2f} {s['asym']:7.2f}  {w['dev']:+8.2f} @ {w['r']:.3f} "
              f"({w['region']:10s})  {s['model_min']['E']:8.2f} @ {s['model_min']['r']:.3f} / {s['ref_min_r']:.3f}")


TERMS = ["Bond", "RepulsionBonded", "RepulsionNonbonded", "Dispersion", "Coulomb", "OverCoord",
         "SqeHardness", "XBond", "HBond"]


def cmd_decomp(argv):
    binp, pair, cfg, rs = argv[0], argv[1], argv[2], [float(x) for x in argv[3].split(",")]
    rows = curve(binp, pair, cfg, grid=rs)
    refs = ref_points(pair)
    x = np.array([r for r, _ in refs]); y = np.array([e for _, e in refs])
    print(f"{pair} {cfg}: terms in kcal/mol relative to the model's own A- + B fragments")
    print(f"{'term':20s} " + " ".join(f"{r:9.4f}" for r in rs))
    for t in TERMS:
        print(f"{t:20s} " + " ".join(f"{row['terms'].get(t, 0.0):9.2f}" for row in rows))
    print(f"{'Total':20s} " + " ".join(f"{row['E']:9.2f}" for row in rows))
    refv = [row["ref"] if not math.isnan(row["ref"]) else float(np.interp(row["r"], x, y)) for row in rows]
    print(f"{'reference':20s} " + " ".join(f"{v:9.2f}" for v in refv))
    print(f"{'Total - ref':20s} " + " ".join(f"{row['E'] - v:9.2f}" for row, v in zip(rows, refv)))
    print(f"{'q(A) / q(B)':20s} " + " ".join(
        f"{row['q'][0]:+.2f}/{row['q'][1]:+.2f}" if row["q"] else "      n/a" for row in rows))


def cmd_react(argv):
    """react-mode breaking (ascending) and forming (descending) chains, 0.05 A steps over the
    reference range (the x2lib protocol of X2_SCOPE 12 / PI_STAR 11.5), rms vs the spline."""
    binp, pairs, cfgs = argv[0], argv[1].split(","), argv[2].split(",")
    jobs = {}
    with ThreadPoolExecutor(int(os.environ.get("X2S_JOBS", "12"))) as ex:
        for p in pairs:
            pts = ref_points(p)
            grid = [float(x) for x in np.round(np.arange(pts[0][0], pts[-1][0] + 1e-9, 0.05), 4)]
            for c in cfgs:
                jobs[(p, c, "break")] = ex.submit(curve, binp, p, c, grid, True)
                jobs[(p, c, "form")] = ex.submit(curve, binp, p, c, grid[::-1], True)
    print(f"{'pair|config|dir':26s} {'n':>3s} {'rms':>7s} {'max dev @ r':>18s}")
    for (p, c, dn), f in jobs.items():
        d = [(x["r"], x["E"] - x["ref"]) for x in f.result() if not math.isnan(x["ref"])]
        w = max(d, key=lambda t: abs(t[1]))
        print(f"{p + '|' + c + '|' + dn:26s} {len(d):3d} {rms([t[1] for t in d]):7.2f} {w[1]:+9.2f} @ {w[0]:.3f}")


def cmd_updown(argv):
    """up-vs-down topology history (FRAG_CHARGE_STATUS section 9 / I2_CLF 8): 0.05 A chains with
    batch topology reuse at the DEFAULT refresh check, |E_up - E_down| at equal r, split at the
    first r where a fresh single point has no A-B bond."""
    binp, pairs, cfgs = argv[0], argv[1].split(","), argv[2].split(",")
    jobs = {}
    with ThreadPoolExecutor(int(os.environ.get("X2S_JOBS", "12"))) as ex:
        for p in pairs:
            pts = ref_points(p)
            g = [float(x) for x in np.round(np.arange(pts[0][0], pts[-1][0] + 1e-9, 0.05), 4)]
            for c in cfgs:
                jobs[(p, c)] = [ex.submit(curve, binp, p, c, g, False, "default"),
                                ex.submit(curve, binp, p, c, g[::-1], False, "default"),
                                ex.submit(curve, binp, p, c, g)]
    for (p, c), (u, d, f) in jobs.items():
        up = {x["r"]: x["E"] for x in u.result()}; dn = {x["r"]: x["E"] for x in d.result()}
        split = next((x["r"] for x in f.result() if abs(x["terms"].get("Bond", 0)) < 1e-9), 1e9)
        dd = [(r, abs(up[r] - dn[r])) for r in up if r in dn]
        lo = [v for r, v in dd if r < split]; hi = [v for r, v in dd if r >= split]
        print(f"{p + '|' + c:18s} split {split:.2f}: r<split max {max(lo, default=0):6.2f} mean "
              f"{(sum(lo) / len(lo)) if lo else 0:5.2f} | r>=split max {max(hi, default=0):6.2f} | "
              f"n>1 kcal {sum(v > 1 for _, v in dd)}/{len(dd)}")


def cmd_steps(argv):
    """fresh single points on a fine grid (default 0.01 A) over the reference range: the largest
    adjacent-point energy step and its r, against the largest analytic |dE/dr| * h (a step well
    above that is a discontinuity, not a steep wall)."""
    binp, pairs, cfgs = argv[0], argv[1].split(","), argv[2].split(",")
    h = float(os.environ.get("X2S_STEP", "0.01"))
    jobs = {}
    with ThreadPoolExecutor(int(os.environ.get("X2S_JOBS", "12"))) as ex:
        for p in pairs:
            pts = ref_points(p)
            g = [float(x) for x in np.round(np.arange(pts[0][0], min(pts[-1][0], 9.0) + 1e-9, h), 4)]
            for c in cfgs:
                jobs[(p, c)] = (g, ex.submit(curve, binp, p, c, g))
    for (p, c), (g, f) in jobs.items():
        rows = f.result()
        st = []
        for a, b in zip(rows[:-1], rows[1:]):
            # only the region past the wall (E < +20): inside the wall a 0.01 A step is legitimately large
            if a["E"] < 20 and b["E"] < 20:
                st.append((abs(b["E"] - a["E"]), a["r"], b["E"] - a["E"]))
        w = max(st)
        print(f"{p + '|' + c:18s} n {len(rows)}: largest step {w[2]:+8.3f} kcal/mol at {w[1]:.3f}->{w[1] + h:.3f} A")


if __name__ == "__main__":
    {"survey": cmd_survey, "decomp": cmd_decomp, "react": cmd_react, "updown": cmd_updown, "steps": cmd_steps}[sys.argv[1]](sys.argv[2:])
