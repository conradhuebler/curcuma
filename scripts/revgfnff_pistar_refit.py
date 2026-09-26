#!/usr/bin/env python3
"""rev-gfnff P3 pi* (O2-/S2-): refit of the pi-excess well rows and harris rows against the
RUNTIME well inputs. Claude Generated (Sep 26, 2026), test_cases/revgfnff/_log/PI_STAR_STATUS.md
section 12. The first fit reconstructed the well with CURCUMA_BONDDUMP's r0_dyn, which lacks the
stage-3a(i) pair-CN correction the bond kernel applies; this one reads CURCUMA_WELLDUMP (the
kernel's own r0/fc/alpha). Order: well (flat, kappa_x 100) -> insert + rebuild -> harris target
(kappa 0, final well). Needs the build_pi binary (or PISTAR_BIN) and the ref/E DLPNO curves.
Usage: python3 scripts/revgfnff_pistar_refit.py   (prints rows; writes refit.json in the CWD)
"""
import json, os, re, subprocess, tempfile, shutil
from pathlib import Path
os.environ.setdefault("TMPDIR", "/var/tmp")
tempfile.tempdir = os.environ["TMPDIR"]
REPO = Path("/home/conrad/src/curcuma_branches/curcuma/.claude/worktrees/agent-a997450b4245937f9")
BIN = os.environ.get("PISTAR_BIN", str(REPO / "build_pi" / "curcuma"))
K2 = 627.5094740631
ANSI = re.compile(r"\x1b\[[0-9;]*m")
BASE = ["-gfnff.rev_enabled", "true", "-gfnff.rev_charge_model", "sqe", "-gfnff.rev_sqe_phase1", "true",
        "-gfnff.rev_excess_electron", "true"]
REC = BASE + ["-gfnff.rev_sqe_virtual_pairs", "true", "-gfnff.rev_excess_mode", "harris",
              "-gfnff.frag_charge_model", "ensemble", "-gfnff.frag_charge_s_max", "1.2"]
PI = ["-gfnff.rev_pi_excess_electron", "true"]

def xyz(atoms):
    return f"{len(atoms)}\nx\n" + "".join(f"{s} {x:.10f} {y:.10f} {z:.10f}\n" for s, x, y, z in atoms)

def dimer(el, r):
    return [(el, 0, 0, 0), (el, 0, 0, r)]

def sp(atoms, flags, charge, spin, method="revgfnff", env_extra=None, grad=False, keep_stdout=False):
    d = Path(tempfile.mkdtemp(prefix="sp_"))
    try:
        (d / "m.xyz").write_text(xyz(atoms))
        cmd = [BIN, "-sp", "m.xyz", "-method", method, "-batch", "true", "-batch_out", "o.jsonl",
               "-gfnff.cache_topology", "false", "-charge", str(charge), "-spin", str(spin),
               "-no_bmt", "-threads", "1", "-verbosity", "2" if keep_stdout else "0"] + list(flags)
        if grad:
            cmd += ["-gradient", "true"]
        env = dict(os.environ); env.update(env_extra or {})
        p = subprocess.run(cmd, capture_output=True, text=True, cwd=d, env=env, timeout=300)
        out = {}
        if (d / "o.jsonl").exists():
            J = [json.loads(l) for l in (d / "o.jsonl").read_text().splitlines() if l.strip()]
            out["E"] = J[0].get("energy_eh") if J else None
            if J and "gradient_eh_ang" in J[0]:
                out["g"] = J[0]["gradient_eh_ang"]
        else:
            out["E"] = None
        if grad and (d / "g.dat").exists():
            out["g"] = [[float(v) for v in ln.split()] for ln in (d / "g.dat").read_text().splitlines()
                        if ln.strip() and not ln.startswith("#")]
        if keep_stdout:
            out["stdout"] = ANSI.sub("", p.stdout)
        out["rc"] = p.returncode
        return out
    finally:
        shutil.rmtree(d, ignore_errors=True)

def chain(frames, flags, charge, spin, method="revgfnff", reuse=True, check=None, env_extra=None, keep_stdout=False):
    d = Path(tempfile.mkdtemp(prefix="ch_"))
    try:
        (d / "s.xyz").write_text("".join(xyz(f) for f in frames))
        cmd = [BIN, "-sp", "s.xyz", "-method", method, "-batch", "true", "-batch_out", "o.jsonl",
               "-gfnff.cache_topology", "false", "-charge", str(charge), "-spin", str(spin),
               "-no_bmt", "-threads", "1", "-verbosity", "2" if keep_stdout else "0"] + list(flags)
        if reuse:
            cmd += ["-batch_reuse_topology", "true"]
        if check is not None:
            cmd += ["-gfnff.reuse_topology_check", "true" if check else "false"]
        env = dict(os.environ); env.update(env_extra or {})
        p = subprocess.run(cmd, capture_output=True, text=True, cwd=d, env=env, timeout=600)
        J = [json.loads(l) for l in (d / "o.jsonl").read_text().splitlines() if l.strip()] if (d / "o.jsonl").exists() else []
        res = [j.get("energy_eh") for j in J]
        return (res, ANSI.sub("", p.stdout)) if keep_stdout else res
    finally:
        shutil.rmtree(d, ignore_errors=True)

def ref(sysname):
    e = json.load(open(REPO / f"test_cases/revgfnff/ref/E/{sysname}/energies.json"))
    pts = [(float(p["label"].split("=")[1]), p["energy_eh"]) for p in e["points"]]
    fr = e["fragment_energies_eh"]
    return pts, fr["atom"] + fr["anion"]


# ---- refit (PI_STAR_STATUS.md section 12) ----
"""Correct-r0 refit of the pi-excess well row and the harris row (PI_STAR_STATUS.md section 12).
Uses the RUNTIME well inputs (CURCUMA_WELLDUMP: r0 after the rev pair-CN correction, fc, alpha)
instead of CURCUMA_BONDDUMP's r0_dyn. Claude Generated (Sep 26, 2026)."""
import re, math, json, sys, importlib.util
h = sys.modules[__name__]   # the helper block above lives in this same file
import numpy as np
from concurrent.futures import ThreadPoolExecutor
spec = importlib.util.spec_from_file_location("wf", str(h.REPO / "scripts/revgfnff_wellfit.py"))
WF = importlib.util.module_from_spec(spec); spec.loader.exec_module(WF)
BA = 1.8897261246257702
RW = re.compile(r"WELLDUMP \d+-\d+ r=(\S+) r0=(\S+) fc=(\S+) alpha=(\S+) form=(\d) well=(\S+) w=(\S+) E=(\S+)")
RB = re.compile(r"\[RESULT\]\s+Bond\s+(-?\d+\.\d+)")
# the pi-excess rows CURRENTLY in rev_well_table_v2.h (the binary must carry exactly these; the
# Bond-term reproduction check below verifies it). First-campaign rows, for reference:
# O (0.789528, 0.926116, 0.443642, 0.392741), S (0.581348, 0.751854, 0.114145, 0.286588).
ROWS = {"O": (0.798346, 1.188780, 1.507767, 0.375625), "S": (0.582195, 0.765328, 0.175854, 0.286630)}

def well_b(par, r, r0, kb, alpha):
    s, ca, beta, dr0 = par
    D = s * kb; a = ca * math.sqrt(alpha / s); bb = beta / BA**2
    x = r - r0 - dr0 * BA
    phi = a * x + bb * x * x
    y = math.exp(-min(max(phi, -50.0), 200.0))
    yc = float(WF._y_cap(np.array([y]))[0])
    return -D * (2 * yc - yc * yc)

def point(el, r, flags):
    o = h.sp(h.dimer(el, r), flags, -1, 1, env_extra={"CURCUMA_WELLDUMP": "1"}, keep_stdout=True)
    corners = []
    bond = None
    for l in o["stdout"].splitlines():
        m = RW.search(l)
        if m:
            rr, r0, fc, al, form, well, w, E = (float(v) for v in m.groups())
            key = (round(r0, 9), round(fc, 9))
            if key not in [c[0] for c in corners]:
                corners.append((key, dict(r=rr, r0=r0, kb=abs(fc), alpha=al, form=int(form), well=well)))
        m = RB.search(l)
        if m:
            bond = float(m.group(1))
    return dict(r=r, E=o["E"], bond=bond, corners=[c[1] for c in corners])

def lam(p):
    c = p["corners"]
    if len(c) == 1:
        return [1.0]
    w1, w2 = c[0]["well"], c[1]["well"]
    l1 = (p["bond"] - w2) / (w1 - w2) if abs(w1 - w2) > 1e-14 else 0.5
    return [l1, 1 - l1]

def bond_new(p, par):
    return sum(l * well_b(par, c["r"], c["r0"], c["kb"], c["alpha"]) for l, c in zip(lam(p), p["corners"]))

def collect(el, flags, pts):
    with ThreadPoolExecutor(8) as ex:
        P = list(ex.map(lambda rr: point(el, rr, flags), [r for r, _ in pts]))
    z1 = h.sp(h.dimer(el, 60.0), flags, -1, 1)["E"]
    return P, z1

def fit_well(el, sysn):
    pts, fsum = h.ref(sysn)
    P, z1 = collect(el, h.BASE + h.PI, pts)
    rows = []
    for (r, er), p in zip(pts, P):
        if er is None or not p["corners"]:
            continue
        rest = p["E"] - p["bond"]
        rows.append((p, (er - fsum) - (rest - z1)))     # the well must supply this [Eh]
    # sanity: the current row reproduces the runtime Bond term
    chk = max(abs(bond_new(p, ROWS[el]) - p["bond"]) for p, _ in rows)
    def loss(par):
        if par[0] <= 0 or par[1] <= 0:
            return 1e6
        return sum((bond_new(p, par) - t) ** 2 for p, t in rows) / len(rows)
    best = None
    for s0 in (0.6, 1.0, 1.4):
        for ca0 in (0.8, 1.2, 1.6):
            for b0 in (0.1, 0.5, 1.5):
                xf, fv = WF.nelder(loss, [s0, ca0, b0, 0.0], [0.2, 0.2, 0.2, 0.1], it=3000)
                xf, fv = WF.nelder(loss, xf, [0.02] * 4, it=3000)
                if best is None or fv < best[1]:
                    best = (xf, fv)
    par = [float(v) for v in best[0]]
    rms_new = math.sqrt(best[1]) * h.K2
    rms_old = math.sqrt(loss(list(ROWS[el]))) * h.K2
    res = [(p["r"], (bond_new(p, par) - t) * h.K2) for p, t in rows]
    return dict(par=par, rms=rms_new, rms_current_row_runtime=rms_old, n=len(rows), check_bond_reprod_eh=chk, resid=res)

def harris_target(el, sysn, par):
    pts, fsum = h.ref(sysn)
    P, z1 = collect(el, h.BASE + ["-gfnff.rev_excess_kappa", "0.0"] + h.PI, pts)
    out = []
    for (r, er), p in zip(pts, P):
        if er is None or not p["corners"]:
            continue
        e_new = p["E"] - p["bond"] + bond_new(p, par)
        out.append((r, ((er - fsum) - (e_new - z1)) * h.K2))
    return out

def fit_harris(tg, cgrid=np.linspace(0.02, 3.0, 1491)):
    r = np.array([a for a, _ in tg]); t = np.array([b for _, b in tg])
    best = None
    scan = []
    for c in cgrid:
        ex = np.exp(-c * r)
        A_ = np.vstack([np.ones(len(r)), -ex]).T
        sol, *_ = np.linalg.lstsq(A_, t, rcond=None)
        rms = math.sqrt(float(np.mean((A_ @ sol - t) ** 2)))
        scan.append((float(c), rms))
        if best is None or rms < best[3]:
            best = (float(sol[0]), float(sol[1]), float(c), rms)
    return best, scan

if __name__ == "__main__":
    out = {}
    for el, sysn in (("O", "o2m_O-O-_dlpno_ccsdt"), ("S", "s2m_S-S-_dlpno_ccsdt")):
        W = fit_well(el, sysn)
        print(el, "WELL", {k: v for k, v in W.items() if k != "resid"})
        print("   resid", [(round(a, 4), round(b, 2)) for a, b in W["resid"]])
        tg_new = harris_target(el, sysn, W["par"])
        tg_old = harris_target(el, sysn, list(ROWS[el]))   # current row
        for nm, tg in (("new-well", tg_new), ("old-well", tg_old)):
            b, scan = fit_harris(tg)
            print(el, "HARRIS", nm, "A=%.4f B=%.4f c=%.4f rms=%.3f" % b, "target", [round(v, 2) for _, v in tg])
            print("   rms(c):", [(round(c, 2), round(rr, 2)) for c, rr in scan[::100]])
        out[el] = dict(well=W, tg_new=tg_new, tg_old=tg_old)
    json.dump(out, open("refit.json", "w"), indent=1)
