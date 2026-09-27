#!/usr/bin/env python3
"""Joint refit of an X-Y half-order well row + harris g(r) row IN the runtime harris setting.

Why (X2_COMPRESSED_SURVEY_STATUS.md, Part 2): the shipped recipe fits the half row against the
flat100 (kappa_x = 100, P3 flat) rest and then g against the harris rest. For ClF- the flat100
charge state localises the electron on F and jumps at the pass-1 split inside the well, so the
difference of the two rests has a 12 kcal/mol step no smooth g can follow (I2_CLF_STATUS 12.2).
This script removes flat100 from the chain: both rows are fitted together against the reference,
on the energy the recommended harris setting actually produces:

    E(p) = E_cur - Bond_cur - H_cur + Bond_cur * f_p(r) / f_cur(r) + x_eff(r) * g_p(r)

with H = the SqeHardness term (= x * g at kappa 0), x_eff = H_cur / g_cur(r), and f(r) the well
per unit |fc| evaluated with the KERNEL's own r0 / alpha (CURCUMA_WELLDUMP, which includes the
stage-3a(i) pair-CN correction that CURCUMA_BONDDUMP's r0_dyn lacks - PI_STAR_STATUS 11.1(a)).
The substitution is exact because D = s|fc| scales linearly with fc, a = ca sqrt(alpha/s) does not
depend on fc, the harris term is not fed back into the charges, and the carrier weights of the
ensemble window do not depend on the well. Checked after the rebuild: fit rms == runtime rms.

Usage: revgfnff_x2_jointrefit.py BIN PAIR CFG OUT.json   (PAIR/CFG as in revgfnff_x2_survey.py)
Current rows are read from the command line: --half s,ca,beta,dr0 --harris A,B,c
Claude Generated (Sep 27, 2026). AI-generated, machine-tested only.
"""
import json, math, os, re, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import revgfnff_x2_survey as S

BOHR = 1.8897261246257702
K = 627.509474
# c range. With c free up to 10 (or 4, the shipped recipe's grid bound) the joint fit turns g into
# a short-range wall (c on the bound, B ~ -1e4 .. -1e7) and shrinks the well - measured, rejected:
# 11 bonded points cannot separate a Morse wall from an exponential g wall. Bounded instead by the
# family's own fitted range (I 0.48, O 0.63, Br 0.68, Cl 0.74, S 0.80, F 1.38); env X2J_CMAX.
C_MIN, C_MAX = 0.05, float(os.environ.get("X2J_CMAX", "1.5"))
G_CONST = os.environ.get("X2J_GCONST", "0") == "1"
B_NONNEG = os.environ.get("X2J_BNEG", "0") != "1"
C_FIX = os.environ.get("X2J_CFIX")  # fix c (a scan over c shows how well it is determined)
RE = re.compile(r"WELLDUMP 1-2 r=(\S+) r0=(\S+) fc=(\S+) alpha=(\S+)")


def welldump(binp, pair, cfg, r):
    A, B = S.PAIRS[pair][:2]
    method, flags = S.configs(pair)[cfg]
    d = tempfile.mkdtemp(dir=S.TMP)
    open(os.path.join(d, "m.xyz"), "w").write(f"2\n\n{A} 0 0 0\n{B} 0 0 {r:.10f}\n")
    env = dict(os.environ, CURCUMA_WELLDUMP="1")
    p = subprocess.run([os.path.abspath(binp), "-sp", "m.xyz", "-method", method, "-charge", "-1",
                        "-threads", "1", "-verbosity", "0", "-no_bmt", "-gfnff.cache_topology", "false"] + flags,
                       cwd=d, capture_output=True, text=True, env=env)
    m = RE.findall(p.stdout)
    if not m:
        return None
    R, r0, fc, al = (float(x) for x in m[0])
    return {"r0": r0, "alpha": al, "fc": fc, "uncap": UNCAP[pair]}


# half-order rows (sigma excess) are uncapped (ff_workspace_gfnff.cpp, P3 uncap = halfOrderWeight
# = 1 at order 0.5); the pi-excess rows of O-O / S-S have no half-order row and keep the y = 2 cap.
UNCAP = {"Cl2-": 1.0, "F2-": 1.0, "Br2-": 1.0, "I2-": 1.0, "ClF-": 1.0, "O2-": 0.0, "S2-": 0.0}


def share_min_one(x, a):
    u = x - 1.0
    if u >= 0.0:
        return 1.0
    if u <= -a:
        return x
    return 1.0 - u * u * u / (a * a) - 2.0 * u * u / a


def f_well(r, v, s, ca, beta, dr0):
    """well per unit |fc| (Eh): the kernel's form 3 (calcBonds), cap at y = 2 blended by uncap."""
    a = ca * math.sqrt(v["alpha"] / s)
    x = r * BOHR - v["r0"] - dr0 * BOHR
    phi = a * x + beta / BOHR ** 2 * x * x
    y = math.exp(-min(max(phi, -50.0), 200.0))
    yc = 2.0 * share_min_one(y / 2.0, 0.2)
    capped = -s * (2 * yc - yc * yc)
    u = v.get("uncap", 1.0)
    return (1.0 - u) * capped + u * (-s * (2 * y - y * y))


def g_of(r, A, B, c):
    return A - B * math.exp(-c * r)


def main():
    binp, pair, cfg, out = sys.argv[1:5]
    half = [float(x) for x in sys.argv[sys.argv.index("--half") + 1].split(",")]
    harris = [float(x) for x in sys.argv[sys.argv.index("--harris") + 1].split(",")]
    rows = S.curve(binp, pair, cfg)
    bonded = [x for x in rows if abs(x["terms"].get("Bond", 0)) > 1e-9 and not math.isnan(x["ref"])]
    with ThreadPoolExecutor(12) as ex:
        wd = list(ex.map(lambda x: welldump(binp, pair, cfg, x["r"]), bonded))
    pts = []
    for x, v in zip(bonded, wd):
        g0 = g_of(x["r"], *harris)
        pts.append(dict(r=x["r"], ref=x["ref"], E=x["E"], bond=x["terms"]["Bond"], H=x["terms"].get("SqeHardness", 0.0),
                        v=v, fcur=f_well(x["r"], v, *half), xeff=x["terms"].get("SqeHardness", 0.0) / g0))

    def model(p, q):
        s, ca, beta, dr0, A, B, c = p
        return q["E"] - q["bond"] - q["H"] + q["bond"] * f_well(q["r"], q["v"], s, ca, beta, dr0) / q["fcur"] \
            + q["xeff"] * g_of(q["r"], A, B, c)

    # sanity: the well formula x |fc| must reproduce the kernel's Bond term wherever the window does
    # not blend (x_eff = 1); a wrong cap / r0 / unit shows up here, not in a self-substitution.
    chk = max((abs(q["bond"] - f_well(q["r"], q["v"], *half) * abs(q["v"]["fc"]) * K)
               for q in pts if abs(q["xeff"] - 1.0) < 1e-6), default=float("nan"))
    print(f"{pair} {cfg}: {len(pts)} bonded points {pts[0]['r']:.3f}-{pts[-1]['r']:.3f} A, "
          f"x_eff {min(q['xeff'] for q in pts):.3f}-{max(q['xeff'] for q in pts):.3f}, kernel-replica check (Bond term) {chk:.2e} kcal/mol")

    def solve_lin(nl, P):
        """given (s, ca, beta, dr0, c): A, B by linear least squares."""
        s, ca, beta, dr0, c = nl
        M, y = [], []
        for q in P:
            base = q["E"] - q["bond"] - q["H"] + q["bond"] * f_well(q["r"], q["v"], s, ca, beta, dr0) / q["fcur"]
            M.append([q["xeff"], -q["xeff"] * math.exp(-c * q["r"])])
            y.append(q["ref"] - base)
        M, y = np.array(M), np.array(y)
        ab, *_ = np.linalg.lstsq(M, y, rcond=None)
        if G_CONST:
            # g = A (B = 0): a pure offset, the harris correction carries no r shape of its own
            a0 = float(np.dot(M[:, 0], y) / np.dot(M[:, 0], M[:, 0]))
            ab = np.array([a0, 0.0])
        elif B_NONNEG and ab[1] < 0.0:
            # g must rise towards its asymptote A (B >= 0), as for every other pair: the sigma*
            # compression wall belongs to the uncapped half-order well (ff_workspace_gfnff.cpp,
            # P3 uncap), not to g. With B < 0 allowed the fit makes g the wall - measured, rejected.
            a0 = float(np.dot(M[:, 0], y) / np.dot(M[:, 0], M[:, 0]))
            ab = np.array([a0, 0.0])
        res = M @ ab - y
        return ab, res

    def resid(z, P):
        if C_FIX is not None:
            z = np.array(z, float); z[4] = math.log(float(C_FIX))
        if not (-4 < z[0] < 3 and -4 < z[1] < 3 and -0.5 < z[3] < 1.5 and -5 < z[4] < 3):
            return np.full(len(P), 1e3)
        s, ca = math.exp(z[0]), math.exp(z[1])
        beta, dr0, c = max(z[2], 0.0), z[3], math.exp(z[4])
        if not (-4 < z[0] < 3 and -4 < z[1] < 3 and -0.5 < dr0 < 1.5 and C_MIN <= c <= C_MAX):
            return np.full(len(P), 1e3)
        _, res = solve_lin((s, ca, beta, dr0, c), P)
        return res

    def lm(z, P, it=200):
        z = np.array(z, float); lam = 1e-3; r = resid(z, P); cst = r @ r
        for _ in range(it):
            J = np.zeros((len(r), 5))
            for k in range(5):
                dz = np.zeros(5); dz[k] = 1e-6 * max(1, abs(z[k])); J[:, k] = (resid(z + dz, P) - r) / dz[k]
            Ah = J.T @ J; gr = J.T @ r
            while True:
                zn = z + np.linalg.solve(Ah + lam * np.diag(np.diag(Ah) + 1e-12), -gr)
                rn = resid(zn, P); cn = rn @ rn
                if cn < cst: z, r, cst = zn, rn, cn; lam = max(lam / 3, 1e-12); break
                lam *= 4
                if lam > 1e12: return z, cst
        return z, cst

    def fit(P, it=200):
        best = None
        for dr0 in (0.0, 0.2, 0.4, 0.6):
            for s0 in (0.7, 1.5, 2.5):
                for ca0 in (0.8, 1.4, 2.0):
                    for c0 in (0.1, 0.7, 2.0, 3.9):
                        z, cst = lm([math.log(s0), math.log(ca0), 0.3, dr0, math.log(c0)], P, it)
                        if best is None or cst < best[1]:
                            best = (z, cst)
        z = np.array(best[0], float)
        if C_FIX is not None:
            z[4] = math.log(float(C_FIX))
        nl = (math.exp(z[0]), math.exp(z[1]), max(z[2], 0.0), z[3], math.exp(z[4]))
        ab, res = solve_lin(nl, P)
        return nl, ab, res

    nl, ab, res = fit(pts)
    s, ca, beta, dr0, c = nl
    A, B = ab
    cur = [q["E"] - q["ref"] for q in pts]
    print(f"current rows: bonded rms {math.sqrt(np.mean(np.square(cur))):.2f}  max {max(abs(x) for x in cur):.2f}")
    print(f"joint fit   : s {s:.6f} ca {ca:.6f} beta {beta:.6f} dr0 {dr0:.6f} | A {A:.6f} B {B:.6f} c {c:.6f}")
    print(f"              bonded rms {math.sqrt(np.mean(np.square(res))):.2f}  max {max(abs(x) for x in res):.2f}")
    print("  residuals: " + " ".join(f"{q['r']:.3f}:{x:+.2f}" for q, x in zip(pts, res)))
    loo = []
    for k in range(len(pts)):
        nlk, abk, _ = fit(pts[:k] + pts[k + 1:], it=100)
        q = pts[k]
        loo.append(model(list(nlk[:4]) + [abk[0], abk[1], nlk[4]], q) - q["ref"])
    print(f"  LOO rms {math.sqrt(np.mean(np.square(loo))):.2f}: " + " ".join(f"{x:+.2f}" for x in loo))
    D0 = s * abs(-0.072340) * K
    print(f"  g(1.5) {g_of(1.5, A, B, c):.1f}  g(2.2) {g_of(2.2, A, B, c):.1f}  g(3.0) {g_of(3.0, A, B, c):.1f} kcal/mol")
    json.dump(dict(half=dict(s=s, ca=ca, beta=beta, dr0=dr0), harris=dict(A=A, B=B, c=c),
                   rms=math.sqrt(np.mean(np.square(res))), loo=math.sqrt(np.mean(np.square(loo))),
                   residuals={f"{q['r']:.4f}": x for q, x in zip(pts, res)}), open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
