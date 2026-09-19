#!/usr/bin/env python3
"""rev-gfnff stage 3a(iii): fit the two curvature-pinned bond-well forms on the class-A set.

Claude Generated (Sep 18, 2026).  AI-fitted parameters - see the header of the generated table.

For every class-A bond type under test_cases/revgfnff/ref/A/ the model is run along the reference
grid with the topology KEPT from the r_eq frame (so the stretched pair stays in the bond list and
there is exactly one topology corner), CURCUMA_SHAREDUMP=1 gives the pair's own well depth
D_p = -k_b e^{-alpha dr^2} w, weight w and share factor c per frame, and the per-frame term
decomposition comes from the -batch JSONL.  From that:

    E_rest  = E_total - E_pair              frozen (the rest of the molecule)
    E_cand  = E_rest + E_form(r - r0)       the candidate curve

with r0, alpha and k_b recovered from ln(D/w) by a quadratic least-squares fit (the same
recovery the review's offline test used; the residual is reported and is <= 1e-2).

Two forms, both with the CURVATURE AT THE MINIMUM PINNED to the delivered Gaussian's
K = 2 alpha |k_b|, so r_min and the force constant are reproduced by construction and only the
DEPTH and the TAIL are free (two parameters per bond type):

    MG        E = -D (2 y - y^2),  y = exp(-(a x + beta x^2)),  a = sqrt(K / (2 D))
    erf-Morse E = -D (2 y - y^2),  y = erfc((x - u)/sigma) / erfc(-u/sigma),
              u solved by bisection from  2 D h(u/sigma)^2 / sigma^2 = K

The objective is the RMS of (candidate - reference) over the BREAK side (from one grid point
inside the reference minimum outwards), each curve referred to its own value at the reference
minimum - the harness definition.

The rest is taken q-FROZEN (the Coulomb term held at its r_eq value) by default: 9 of 32 bonds
drift by >= 10 kcal/mol in this static protocol and the SIGN of that drift flips between the
static and the react protocol (FABLE_REVIEW_2 B.2, EEQ_DRIFT_STATUS), so a depth fitted on the
raw rest absorbs a charge-model error that is not the well's.  --rest raw fits on the raw rest.

Output: a JSON table keyed by class-A bond type and, aggregated, by ELEMENT PAIR (median over
the systems that carry that pair), with s = D / |k_b| and the tail parameter.  --header writes
the compiled C++ table the force field reads.
"""
import argparse, hashlib, importlib.util, json, math, os, re, shutil, subprocess, sys, tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("classa", REPO / "scripts" / "revgfnff_classa.py")
H = importlib.util.module_from_spec(spec)
spec.loader.exec_module(H)

K2 = 627.5094740631          # Eh -> kcal/mol
B2A = 0.529177210903         # Bohr -> Angstrom
ANSI = re.compile(r"\x1b\[[0-9;]*m")
RD = re.compile(r"shareD\s+\d+\s+(\d+)-\s*(\d+) r\s+(\S+) D\s+(\S+) w\s+(\S+) c\s+(\S+) E\s+(\S+)")
RB = re.compile(r"BONDPARAM\s+(\d+)\(\d+\)-(\d+)\(\d+\).*\border=([\d.eE+-]+)")


# --------------------------------------------------------------------------- scan
def scan(binary, base, workdir, extra=()):
    """One kept-topology scan of one bond type -> per-frame energies, terms and shareD rows."""
    ref = H.load_reference(base)
    rows = ref["rows"]
    imin = min(range(len(rows)), key=lambda j: rows[j]["energy"])
    frames = [rows[imin]["frame"]] + [e["frame"] for e in rows]     # frame 0 fixes the topology
    d = Path(tempfile.mkdtemp(dir=workdir))
    try:
        H.write_frames(d / "scan.xyz", frames, ["f%d" % k for k in range(len(frames))])
        cmd = [str(binary), "-sp", "scan.xyz", "-method", "revgfnff", "-batch", "true",
               "-batch_out", "out.jsonl", "-batch_reuse_topology", "true",
               "-gfnff.cache_topology", "false", "-charge", str(ref["charge"]),
               "-spin", str(ref["mult"] - 1), "-no_bmt", "-threads", "1", "-verbosity", "1",
               # The scan MUST run the DELIVERED Gaussian well, whatever the current default is.
               # prepare() recovers (r0, alpha, k_b) by fitting a quadratic to ln(D/w), which is
               # exact only if D IS k_b exp(-alpha dr^2) w.  When this script was first written
               # `gauss` was the default and the flag was unnecessary; the Sep 19 flip to `mg`
               # made the omission silently wrong - the recovery then fits a Gaussian to an MG
               # curve, and the whole table inherits the error.  Measured: without this flag the
               # reported "delivered" break rms is 9.4 (that is MG's) instead of 19.2, and the
               # free-curvature table that came out of it made the class-A harness WORSE, not
               # better.  Claude Generated (Sep 20, 2026).
               "-gfnff.rev_well_form", "gauss",
               "-gfnff.reuse_topology_check", "false"] + list(extra)
        env = dict(os.environ)
        env["CURCUMA_SHAREDUMP"] = "1"
        # rev-gfnff stage 3b (Claude Generated, Sep 20, 2026): the same run also dumps the
        # per-bond CONTINUOUS bond order (Bond::rev_order, printed by CURCUMA_BONDDUMP), so the
        # order table is keyed on exactly the number the force field will compute at runtime and
        # not on a chemical label parsed from the directory name. The two disagree in one place
        # and it matters: GFN-FF gives O2 hyb = 1 on both ends, so its ORDER comes out 3, not the
        # nominal 2 - keying on the runtime value is what keeps the table self-consistent.
        env["CURCUMA_BONDDUMP"] = "1"
        p = subprocess.run(cmd, capture_output=True, text=True, cwd=d, env=env)
        orders = {}
        for ln in p.stdout.splitlines():
            m = RB.match(ANSI.sub("", ln).strip())
            if m:
                key = (min(int(m.group(1)), int(m.group(2))), max(int(m.group(1)), int(m.group(2))))
                orders.setdefault(key, float(m.group(3)))
        J = [json.loads(l) for l in (d / "out.jsonl").read_text().splitlines()]
        blocks, cur = [], None
        for ln in p.stdout.splitlines():
            ln = ANSI.sub("", ln).replace("[RESULT]", "").strip()
            if ln.startswith("share dump"):
                cur = []
                blocks.append(cur)
                continue
            m = RD.match(ln)
            if m and cur is not None:
                cur.append((int(m.group(1)), int(m.group(2)), float(m.group(3)),
                            float(m.group(4)), float(m.group(5)), float(m.group(6))))
        return dict(bond=base, r=[e["r"] for e in rows], Eref=[e["energy"] for e in rows],
                    E=[j.get("energy_eh") for j in J][1:], terms=[j.get("terms") for j in J][1:],
                    blocks=blocks[1:], charge=ref["charge"], orders=orders,
                    symbols=[t[0] for t in frames[0]])
    finally:
        shutil.rmtree(d, ignore_errors=True)


def prepare(rec):
    """Pick the stretched pair, recover (r0, alpha, k_b) and build the frozen rest."""
    ok = [k for k, e in enumerate(rec["E"]) if e is not None]
    r = np.array([rec["r"][k] for k in ok])
    Eref = np.array([rec["Eref"][k] for k in ok]) * K2
    E = np.array([rec["E"][k] for k in ok]) * K2
    blocks = [rec["blocks"][k] for k in ok]
    terms = [rec["terms"][k] for k in ok]
    keys = set((p[0], p[1]) for p in blocks[0])
    best = None
    for key in keys:
        sel = [[p for p in b if (p[0], p[1]) == key] for b in blocks]
        if any(len(x) != 1 for x in sel):
            continue
        v = np.array([x[0][2] for x in sel])
        if best is None or v.std() > best[1]:
            best = (key, v.std(), sel)
    if best is None:
        return None
    key, _, sel = best
    rp = np.array([x[0][2] for x in sel]) * B2A
    D = np.array([x[0][3] for x in sel]) * K2
    w = np.array([x[0][4] for x in sel])
    c = np.array([x[0][5] for x in sel])
    rest = E - (-D * c)
    coul = np.array([t["Coulomb"] for t in terms]) * K2
    i0 = int(np.argmin(Eref))
    m = (w > 0.5) & (D / np.maximum(w, 1e-12) > 1e-3)
    if m.sum() < 4:
        return None
    A = np.vstack([np.ones(m.sum()), rp[m], rp[m] ** 2]).T
    coef = np.linalg.lstsq(A, np.log(D[m] / w[m]), rcond=None)[0]
    alpha = -coef[2]
    r0 = coef[1] / (2 * alpha)
    kb = math.exp(coef[0] + alpha * r0 * r0)          # positive depth
    res = float(np.max(np.abs(np.log(D[m] / w[m]) - A @ coef)))
    sym = rec.get("symbols") or []
    # The element pair comes from the ATOM INDICES of the stretched pair, not from the directory
    # name: the class-A tags are chemical shorthand ("ch3oh_HO-H" is an O-H bond, "hcn_HC-H" a
    # C-H bond, "c2h2_CTC" a C-C triple bond), which no name parser should have to know.
    e1 = sym[key[0] - 1] if 0 < key[0] <= len(sym) else "?"
    e2 = sym[key[1] - 1] if 0 < key[1] <= len(sym) else "?"
    order = (rec.get("orders") or {}).get((min(key), max(key)))
    return dict(r=r, Eref=Eref, E=E, rest=rest, drift=coul - coul[i0], i0=i0,
                alpha=float(alpha), r0=float(r0), kb=float(kb), gres=res, pair=key,
                order=order, elems=tuple(sorted((e1, e2))))


# --------------------------------------------------------------------------- forms
def mg(D, beta, K, x):
    a = math.sqrt(K / (2.0 * D))
    phi = a * x + beta * x * x
    y = np.exp(-np.clip(phi, -50.0, 200.0))
    return -D * (2.0 * y - y * y)



def d2_at(r, e, i):
    """3-point Lagrange second derivative at index i (the same formula as the class-A harness's
    second_derivative, on a curve that is already in kcal/mol)."""
    if i <= 0 or i >= len(r) - 1:
        return float("nan")
    h1, h2 = r[i] - r[i - 1], r[i + 1] - r[i]
    return 2.0 * (e[i - 1] * h2 - e[i] * (h1 + h2) + e[i + 1] * h1) / (h1 * h2 * (h1 + h2))



def vertex(r, e, i=None):
    """Sub-grid minimum of a sampled curve: the vertex of the parabola through the three points
    around its sampled minimum.  The class-A grid is r_eq * 1.0668^k, i.e. ~0.065 Angstrom around
    an O-H bond, so the SAMPLED minimum cannot resolve the 0.02-0.03 Angstrom shifts the r0
    re-solve is about; the vertex can."""
    if i is None:
        i = int(np.argmin(e))
    if i <= 0 or i >= len(r) - 1:
        return float(r[i])
    h1, h2 = r[i] - r[i - 1], r[i + 1] - r[i]
    d1 = (e[i + 1] - e[i]) / h2
    d0 = (e[i] - e[i - 1]) / h1
    d2 = 2.0 * (d1 - d0) / (h1 + h2)
    if d2 <= 0:
        return float(r[i])
    slope = (d1 * h1 + d0 * h2) / (h1 + h2)          # centred first derivative at i
    return float(r[i] - slope / d2)


DR0_FIXED = {}
DR0_MAX = 0.40          # Angstrom: the largest r0 RE-SOLVE accepted (see fit())


def _y_cap(y):
    """The runtime's inner-side cap on y, replicated EXACTLY (FFWorkspace::calcBonds): y is run
    through 2 * shareMinOne(y/2, 0.2), so it is the literal y below y = 1.6 and saturates at 2.
    The curvature-pinned forms never came near it on the class-A grids, but a fitted curvature of
    up to ~2.5 x the Gaussian's does at the innermost compressed points, and a fit that ignored it
    would not be the curve the force field evaluates."""
    t = y / 2.0
    u = t - 1.0
    out = np.where(u >= 0.0, 1.0,
                   np.where(u <= -0.2, t, 1.0 - u ** 3 / 0.04 - 2.0 * u * u / 0.2))
    return 2.0 * out


def mg2(D, ca, beta, dr0, K, x):
    """rev-gfnff 3a(iii) step 2 (Claude Generated, Sep 20, 2026): MG with FREE curvature and a
    re-solved minimum.  Same E = -D (2y - y^2), y = exp(-(a x' + beta x'^2)), but

        x'   = x - dr0            the minimum is moved by the fitted offset dr0 [Angstrom]
        a    = ca sqrt(K / (2 D)) the curvature at the minimum is K_well = 2 D a^2 = ca^2 K,
                                  i.e. ca = 1 is exactly the curvature-pinned form above

    so the four free parameters are (D, ca, beta, dr0) instead of (D, beta).  ca and dr0 are
    dimensionless resp. a length, so they transfer to the force field's Bohr units unchanged /
    by one conversion - see rev_well_table_v2.h."""
    a = ca * math.sqrt(K / (2.0 * D))
    xx = x - dr0
    phi = a * xx + beta * xx * xx
    y = _y_cap(np.exp(-np.clip(phi, -50.0, 200.0)))
    return -D * (2.0 * y - y * y)


def em_u(D, sigma, K):
    """u with 2 D h(u/sigma)^2/sigma^2 = K; h(z) = (2/sqrt(pi)) e^{-z^2}/erfc(z) at z = -u/sigma."""
    tgt = math.sqrt(K / (2.0 * D))

    def g(u):
        z = -u / sigma
        if z > 25.0:
            return 2.0 * z / sigma
        return 2.0 / (sigma * math.sqrt(math.pi)) * math.exp(-z * z) / math.erfc(z)

    lo, hi = -40.0 * sigma, 8.0 * sigma
    if g(lo) < tgt or g(hi) > tgt:
        return None
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if g(mid) > tgt:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def em(D, sigma, K, x):
    u = em_u(D, sigma, K)
    if u is None:
        return None
    n0 = math.erfc(-u / sigma)
    if n0 < 1e-12:
        return None
    y = np.array([math.erfc((xx - u) / sigma) for xx in x]) / n0
    return -D * (2.0 * y - y * y)


def nelder(f, x0, step, it=3000, tol=1e-10):
    n = len(x0)
    S = [np.array(x0, float)]
    for i in range(n):
        y = np.array(x0, float)
        y[i] += step[i]
        S.append(y)
    fv = [f(s) for s in S]
    for _ in range(it):
        o = np.argsort(fv)
        S = [S[i] for i in o]
        fv = [fv[i] for i in o]
        if abs(fv[-1] - fv[0]) < tol:
            break
        c = np.mean(S[:-1], axis=0)
        xr = c + (c - S[-1])
        fr = f(xr)
        if fr < fv[0]:
            xe = c + 2 * (c - S[-1])
            fe = f(xe)
            S[-1], fv[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < fv[-2]:
            S[-1], fv[-1] = xr, fr
        else:
            xc = c + 0.5 * (S[-1] - c)
            fc = f(xc)
            if fc < fv[-1]:
                S[-1], fv[-1] = xc, fc
            else:
                S = [S[0]] + [S[0] + 0.5 * (s - S[0]) for s in S[1:]]
                fv = [fv[0]] + [f(s) for s in S[1:]]
    i = int(np.argmin(fv))
    return S[i], fv[i]


def metrics(r, e):
    i = int(np.argmin(e))
    req = r[i]
    de = e[-1] - e[i]
    rel = e - e[i]

    def rise(fr):
        t = fr * de
        for k in range(len(r) - 1):
            if r[k] < req - 1e-9:
                continue
            if rel[k] <= t <= rel[k + 1]:
                return (r[k] + (r[k + 1] - r[k]) * (t - rel[k]) / (rel[k + 1] - rel[k])) / req
        return float("nan")

    return de, rise(0.5), rise(0.9)


def fit(d, form, rest, seed=None, dr0_fixed=None):
    x = d["r"] - d["r0"]
    ref = d["Eref"] - d["Eref"][d["i0"]]
    msk = np.arange(len(x)) >= d["i0"] - 1
    Kc = 2.0 * d["alpha"] * d["kb"]

    # ----------------------------------------------------------------- mg2: the curvature
    # rev-gfnff 3a(iii) step 2 (Claude Generated, Sep 20, 2026).  The curvature at the minimum is
    # RE-SOLVED FROM THE REFERENCE, not fitted and not inherited from the Gaussian.  Measured
    # reason for not fitting it: the break-side objective barely constrains it (the well's second
    # derivative at its own minimum is one point of a 15-point branch), so a free ca is used by
    # the simplex as a spare tail knob and runs to ZERO - 8 of 32 bonds came out with a flat,
    # QUARTIC bottom (ca -> 0 means E''(r_min) = 0), which looks excellent on a break-side rms
    # and would destroy every vibrational force constant.  Fitting on the whole grid instead does
    # constrain it, but then the five compressed points - where the REPULSION term carries the
    # error, not the well - drag the well off: f2 went from a 1.4 break-side rms to 13.1.
    #
    # The well's curvature is exactly 2 D a^2 = ca^2 K_gauss (beta enters only at x^3), and it is
    # the only term whose curvature the well form changes, so
    #
    #     k_model + (ca^2 - 1) K_gauss = k_ref   ->   ca^2 = 1 + (k_ref - k_model) / K_gauss
    #
    # with k_ref and k_model the measured 3-point curvatures of the reference and of the DELIVERED
    # model at their own sampled minima - the same quantity the class-A harness prints as
    # "k ref/model".  Zero knobs, one closed form, both inputs measured.
    i0 = d["i0"]
    k_ref = d2_at(d["r"], ref, i0)
    im = int(np.argmin(d["E"]))
    k_mod = d2_at(d["r"], d["E"] - d["E"][im], im)
    ca_solved, ca_clamped = 1.0, False
    dr0_state = [0.0]
    r_eq_ref = vertex(d["r"], ref, i0)
    if form == "mg2":
        if math.isnan(k_ref) or math.isnan(k_mod) or Kc <= 0.0:
            ca_solved, ca_clamped = 1.0, True
        else:
            c2 = 1.0 + (k_ref - k_mod) / Kc
            ca_solved = math.sqrt(c2) if c2 > 0.09 else 0.3
            ca_clamped = not (c2 > 0.09)
            if ca_solved > 3.0:
                ca_solved, ca_clamped = 3.0, True

    def tot(p):
        D = p[0]
        if D <= 0:
            return None
        if form == "mg":
            e = mg(D, abs(p[1]), Kc, x)
        elif form == "mg2":
            e = mg2(D, ca_solved, abs(p[1]), dr0_state[0], Kc, x)
        else:
            e = em(D, abs(p[1]) + 1e-6, Kc, x)
        if e is None or not np.all(np.isfinite(e)):
            return None
        return rest + e

    # rev-gfnff 3a(iii) step 2 (Claude Generated, Sep 20, 2026): the FREE-curvature form is fitted
    # on the WHOLE reference grid, not on the break side alone.  Measured, and the reason: with
    # only the break side the curvature at the minimum is barely constrained, so the simplex uses
    # it as a spare tail knob and drives it to ZERO - 8 of 32 bonds came out with a flat, quartic
    # bottom (ca -> 0 means E'' (r_min) = 0), which would destroy every vibrational force constant
    # while looking excellent on a break-side rms.  The class-A curves carry FIVE points inside
    # r_eq; including them is what "re-solve r0 and the curvature from the reference data" means,
    # and the full grid is also exactly the quantity the class-A harness reports as its rms.  The
    # two curvature-PINNED forms keep the break-side objective so their package-4 numbers stay
    # reproducible - their curvature is not a free parameter, so the inner points add nothing.
    def obj(p):
        t = tot(p)
        if t is None:
            return 1e9
        return float(np.sqrt(np.mean((((t - t[d["i0"]]) - ref)[msk]) ** 2)))

    def rms_full(p):
        t = tot(p)
        if t is None:
            return float("nan")
        return float(np.sqrt(np.mean(((t - t[d["i0"]]) - ref) ** 2)))

    if form == "mg":
        starts = [[d["kb"] * 1.2, 0.5], [d["kb"] * 1.5, 0.2], [d["kb"], 1.0]]
    elif form == "mg2":
        D0 = seed["D"] if seed else d["kb"]
        b0 = seed["tail"] if seed else 0.5
        starts = [[D0, b0], [D0 * 1.3, b0 * 0.6], [D0 * 0.8, b0 * 1.5]]
    else:
        starts = [[d["kb"] * 1.2, 1.2], [d["kb"] * 1.5, 0.8], [d["kb"], 2.0]]
    def solve_dr0(p):
        """rev-gfnff 3a(iii) step 2 (Claude Generated, Sep 20, 2026): the r0 RE-SOLVE.

        dr0 is not a fitted parameter either - it is SOLVED so that the candidate curve's own
        sub-grid minimum sits on the REFERENCE minimum.  Measured reason: a dr0 fitted on the
        break-side objective does not target r_eq at all, and it moved the wrong way - the
        pair-table H-O offset came out -0.052 Angstrom and the optimised water O-H went
        0.9727 -> 0.9386 against a reference r_eq of ~0.962, i.e. FURTHER from the reference
        than the curvature-pinned form, at 5x the equilibrium shift the operator's guard
        precedent allows.  Solving it instead is what "re-solve r0 from the reference data"
        means, and it costs no free parameter.

        Monotone: moving the well outward moves the total minimum outward (d r*/d dr0 =
        K_well/k_total > 0), so a bisection is exact.  D and beta barely enter - the well's
        curvature is ca^2 K_gauss independent of D, and its own minimum is at x = 0 for any
        beta - which is why one solve per (D, beta) pair converges in a couple of outer
        iterations."""
        lo, hi = -DR0_MAX, DR0_MAX
        keep = dr0_state[0]
        def rq(u):
            dr0_state[0] = u
            t_ = tot(p)
            return None if t_ is None else vertex(d["r"], t_)
        a_, b_ = rq(lo), rq(hi)
        if a_ is None or b_ is None:
            dr0_state[0] = keep
            return keep
        if not (a_ <= r_eq_ref <= b_):
            dr0_state[0] = lo if abs(a_ - r_eq_ref) < abs(b_ - r_eq_ref) else hi
            return dr0_state[0]
        for _ in range(50):
            mid = 0.5 * (lo + hi)
            if rq(mid) < r_eq_ref:
                lo = mid
            else:
                hi = mid
        dr0_state[0] = 0.5 * (lo + hi)
        return dr0_state[0]

    best = None
    for s in starts:
        if form == "mg2" and dr0_fixed is not None:
            # rev-gfnff 3a(iii) step 2 (Claude Generated, Sep 20, 2026): dr0 supplied from the
            # RELAXED-geometry solve (see WORK_STATUS package 9a).  The offline rigid-scan solve
            # does not transfer to an optimisation - it moved the optimised water O-H from 0.9727
            # to 0.9179 against a reference 0.9644 - so the offset is measured against the
            # OPTIMISED class-A bond lengths instead, and only (D, beta) are fitted here.
            dr0_state[0] = dr0_fixed
            p = list(s)
            p, v = nelder(obj, p, [0.1 * abs(q) + 0.05 for q in p])
            p, v = nelder(obj, p, [0.02 * abs(q) + 0.01 for q in p])
            if best is None or v < best[1]:
                best = (p, v, dr0_fixed)
            continue
        if form == "mg2":
            # alternate: fit (D, beta) at the current dr0, re-solve dr0, repeat.  Two rounds are
            # enough - the third moves dr0 by < 1e-4 Angstrom on every class-A bond.
            dr0_state[0] = 0.0
            p = list(s)
            for _ in range(3):
                p, v = nelder(obj, p, [0.1 * abs(q) + 0.05 for q in p])
                solve_dr0(p)
            p, v = nelder(obj, p, [0.02 * abs(q) + 0.01 for q in p])
            solve_dr0(p)
            v = obj(p)
            if best is None or v < best[1]:
                best = (p, v, dr0_state[0])
            continue
        p, v = nelder(obj, s, [0.1 * abs(q) + 0.05 for q in s])
        p, v = nelder(obj, p, [0.02 * abs(q) + 0.01 for q in p])
        if best is None or v < best[1]:
            best = (p, v, 0.0)
    p, v = best[0], best[1]
    dr0_state[0] = best[2]
    t = tot(p)
    out = dict(D=float(p[0]), tail=float(abs(p[1])), s=float(p[0] / d["kb"]), rms=float(v),
               metrics=[float(q) for q in metrics(d["r"], t)],
               join=float(t[-1] - rest[-1]))
    if form == "mg2":
        out["ca"] = float(ca_solved)
        out["ca_clamped"] = bool(ca_clamped)
        out["dr0"] = float(dr0_state[0])
        out["dr0_at_bound"] = bool(abs(abs(dr0_state[0]) - DR0_MAX) < 1e-4)
        out["r_eq_ref"] = r_eq_ref
        out["r_eq_model"] = vertex(d["r"], t)
        out["rms_full"] = rms_full(p)
    return out


ELEM = {1: "H", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 14: "Si", 15: "P", 16: "S", 17: "Cl",
        35: "Br", 53: "I"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", default=str(REPO / "build_rev" / "curcuma"))
    ap.add_argument("--json", default="wellfit.json")
    ap.add_argument("--header", default="", help="write the compiled C++ table here")
    ap.add_argument("--header-v2", dest="header_v2", default="",
                    help="write the FREE-CURVATURE table (rev_well_table_v2.h) here")
    ap.add_argument("--rest", choices=["qfrozen", "raw"], default="qfrozen")
    ap.add_argument("--dr0-file", dest="dr0_file", default="",
                    help="JSON {bond type: dr0 [Angstrom]} - fixes the mg2 r0 offset instead of "
                         "solving it on the rigid scan (see fit(), dr0_fixed)")
    ap.add_argument("--workdir", default="")
    ap.add_argument("--systems", nargs="*", default=None)
    a = ap.parse_args()
    global DR0_FIXED
    if a.dr0_file:
        DR0_FIXED = json.load(open(a.dr0_file))
    work = a.workdir or tempfile.mkdtemp(prefix="wellfit_")
    Path(work).mkdir(parents=True, exist_ok=True)
    bonds = a.systems if a.systems else H.bond_types()
    md5 = hashlib.md5(Path(a.binary).read_bytes()).hexdigest()
    print("# binary %s md5 %s  rest %s  %d bond types" % (a.binary, md5, a.rest, len(bonds)))
    out = dict(binary=a.binary, md5=md5, rest=a.rest, bonds={})
    for b in bonds:
        rec = scan(a.binary, b, work)
        d = prepare(rec)
        if d is None:
            print("%-22s SKIP (no usable pair)" % b)
            continue
        rest = d["rest"] - (d["drift"] if a.rest == "qfrozen" else 0.0)
        ref = d["Eref"]
        rms_del = float(np.sqrt(np.mean((((d["E"] - d["E"][d["i0"]]) - (ref - ref[d["i0"]]))
                                         [np.arange(len(ref)) >= d["i0"] - 1]) ** 2)))
        f_mg = fit(d, "mg", rest)
        row = dict(r0=d["r0"], alpha=d["alpha"], kb=d["kb"], gres=d["gres"],
                   rms_delivered=rms_del, pair="-".join(d["elems"]),
                   order=d["order"],
                   order_key=(None if d["order"] is None else int(round(d["order"]))),
                   mg=f_mg, em=fit(d, "em", rest), mg2=fit(d, "mg2", rest, seed=f_mg, dr0_fixed=DR0_FIXED.get(b)))
        out["bonds"][b] = row
        print("%-22s pair %-6s ord %s kb %8.2f alpha %6.3f | delivered rms %7.2f | MG s %5.3f beta %6.3f rms %6.2f"
              " | EM rms %6.2f | MG2 s %5.3f ca %5.3f%s beta %6.3f dr0 %+7.4f%s rms %6.2f"
              % (b, row["pair"], "-" if d["order"] is None else "%.2f" % d["order"],
                 d["kb"], d["alpha"], rms_del,
                 row["mg"]["s"], row["mg"]["tail"], row["mg"]["rms"],
                 row["em"]["rms"],
                 row["mg2"]["s"], row["mg2"]["ca"], "!" if row["mg2"]["ca_clamped"] else " ",
                 row["mg2"]["tail"], row["mg2"]["dr0"],
                 "!" if row["mg2"]["dr0_at_bound"] else " ", row["mg2"]["rms"]))
    # aggregate to element pairs
    agg = {}
    for b, row in out["bonds"].items():
        agg.setdefault(row["pair"], []).append(row)
    out["pairs"] = {}
    for p, rows in sorted(agg.items()):
        out["pairs"][p] = dict(
            n=len(rows),
            mg_s=float(np.median([r["mg"]["s"] for r in rows])),
            mg_tail=float(np.median([r["mg"]["tail"] for r in rows])),
            em_s=float(np.median([r["em"]["s"] for r in rows])),
            em_tail=float(np.median([r["em"]["tail"] for r in rows])),
            mg_rms=float(np.median([r["mg"]["rms"] for r in rows])),
            em_rms=float(np.median([r["em"]["rms"] for r in rows])),
            mg2_s=float(np.median([r["mg2"]["s"] for r in rows])),
            mg2_ca=float(np.median([r["mg2"]["ca"] for r in rows])),
            mg2_beta=float(np.median([r["mg2"]["tail"] for r in rows])),
            mg2_dr0=float(np.median([r["mg2"]["dr0"] for r in rows])),
            mg2_rms=float(np.median([r["mg2"]["rms"] for r in rows])))
    # stage 3b: the same aggregation with the BOND ORDER in the key.  The order is the runtime
    # value (rounded to the nearest integer for the key, see scan()); every class-A pair that
    # carries more than one order is resolved this way, the rest keep an n = 1 or n = 2 row that
    # behaves exactly like the element-pair entry.
    agg2 = {}
    for b, row in out["bonds"].items():
        if row["order_key"] is None:
            continue
        agg2.setdefault((row["pair"], row["order_key"]), []).append((b, row))
    out["pair_orders"] = {}
    for (p, o), named in sorted(agg2.items()):
        rows = [r for _, r in named]
        out["pair_orders"]["%s|%d" % (p, o)] = dict(
            n=len(rows), pair=p, order=o,
            members=sorted(nm for nm, _ in named),
            mg2_s=float(np.median([r["mg2"]["s"] for r in rows])),
            mg2_ca=float(np.median([r["mg2"]["ca"] for r in rows])),
            mg2_beta=float(np.median([r["mg2"]["tail"] for r in rows])),
            mg2_dr0=float(np.median([r["mg2"]["dr0"] for r in rows])),
            mg2_rms=float(np.median([r["mg2"]["rms"] for r in rows])))
    json.dump(out, open(a.json, "w"), indent=1)
    print("# wrote", a.json, "-", len(out["bonds"]), "bond types,", len(out["pairs"]),
          "element pairs,", len(out["pair_orders"]), "pair/order keys")
    for tag, key in (("delivered", "rms_delivered"), ("MG", "mg"), ("EM", "em"), ("MG2", "mg2")):
        v = sorted((r[key] if key == "rms_delivered" else r[key]["rms"])
                   for r in out["bonds"].values())
        if v:
            print("# per-system break rms, %-9s median %7.3f  max %8.3f  n(<3) %d/%d"
                  % (tag, v[len(v) // 2], v[-1], sum(1 for q in v if q < 3.0), len(v)))
    if a.header:
        write_header(a.header, out)
        print("# wrote", a.header)
    if a.header_v2:
        write_header_v2(a.header_v2, out)
        print("# wrote", a.header_v2)


def write_header(path, out):
    inv = {v: k for k, v in ELEM.items()}
    rows = []
    for p, v in sorted(out["pairs"].items()):
        e1, e2 = p.split("-")
        if e1 not in inv or e2 not in inv:
            continue
        z1, z2 = sorted((inv[e1], inv[e2]))
        rows.append((z1, z2, v))
    with open(path, "w") as f:
        f.write("""/*
 * rev-gfnff stage 3a(iii): per-element-pair parameters of the two curvature-pinned bond-well
 * forms.  GENERATED by scripts/revgfnff_wellfit.py - do not edit by hand.
 *
 * AI-FITTED (Claude Generated, Sep 2026), machine-tested only.  Source: the class-A reference
 * scans under test_cases/revgfnff/ref/A/ (r2SCAN-3c, merged min(RKS, UKS) with the QUALITY.md
 * exclusions), fitted on the BREAK side against the CHARGE-FROZEN rest of the molecule.  Per
 * element pair the value is the MEDIAN over the class-A systems that carry that pair.
 *
 * s    = D / |k_b|, the depth scale on the delivered Gaussian's own force constant
 * tail = beta (MG, 1/Angstrom^2) or sigma (erf-Morse, Angstrom)
 *
 * Binary that produced the scans: md5 %s, rest %s.
 */
#pragma once

#include <cstddef>

namespace RevWellTable {

struct Entry {
    int z1, z2;
    double mg_s, mg_beta;
    double em_s, em_sigma;
};

static constexpr Entry kEntries[] = {
""" % (out["md5"], out["rest"]))
        for z1, z2, v in rows:
            f.write("    { %2d, %2d, %10.6f, %10.6f, %10.6f, %10.6f },   // %s (n = %d)\n"
                    % (z1, z2, v["mg_s"], v["mg_tail"], v["em_s"], v["em_tail"],
                       ELEM[z1] + "-" + ELEM[z2], v["n"]))
        f.write("""};

static constexpr std::size_t kCount = sizeof(kEntries) / sizeof(kEntries[0]);

/// The entry for an element pair, or nullptr when the class-A set has no data for it - the
/// caller then falls back to the delivered Gaussian and says so at verbosity 2.
inline const Entry* find(int za, int zb)
{
    const int z1 = za < zb ? za : zb;
    const int z2 = za < zb ? zb : za;
    for (std::size_t k = 0; k < kCount; ++k)
        if (kEntries[k].z1 == z1 && kEntries[k].z2 == z2)
            return &kEntries[k];
    return nullptr;
}

} // namespace RevWellTable
""")



HEADER_V2_DOC = r"""/*
 * rev-gfnff stage 3a(iii) step 2 + stage 3b: parameters of the FREE-CURVATURE MG bond well.
 * GENERATED by scripts/revgfnff_wellfit.py --form mg2 - do not edit by hand.
 *
 * AI-FITTED (Claude Generated, Sep 2026), machine-tested only.  Source: the class-A reference
 * scans under test_cases/revgfnff/ref/A/ (r2SCAN-3c, merged min(RKS, UKS) with the QUALITY.md
 * exclusions), fitted on the BREAK side against the CHARGE-FROZEN rest of the molecule.
 *
 * The form is the same MG well as rev_well_table.h,
 *
 *     E = -D (2 y - y^2),   y = exp(-(a x + beta x^2)),   x = r - (r0_model + dr0),
 *
 * but with FOUR free parameters per key instead of two - the curvature at the minimum and the
 * position of the minimum are no longer inherited from the delivered Gaussian:
 *
 *   s    = D / |k_b|                     depth scale on the delivered force constant
 *   ca   = a / sqrt(alpha / s)           CURVATURE scale; the curvature at the minimum is
 *                                        K = 2 D a^2 = ca^2 * (2 alpha |k_b|), so ca = 1 is
 *                                        exactly the curvature-pinned table (rev_well_table.h)
 *   beta = the MG tail parameter [1/Angstrom^2]
 *   dr0  = the r0 RE-SOLVE [Angstrom], an offset ON TOP of the model's own dynamic r0, so the
 *          stage-3a(i) CN dependence of r0 is preserved and only its class-A offset is fitted
 *
 * ca = 1, dr0 = 0 reproduces rev_well_table.h's mg entry bit-for-bit, which is how the code path
 * was verified before the free-curvature fit was run.
 *
 * Binary that produced the scans: md5 %s, rest %s.
 *
 * kPairEntries is keyed on the ELEMENT PAIR (well form 'mg2'); kOrderEntries adds the nominal
 * bond ORDER (1/2/3) and is interpolated at runtime on the CONTINUOUS order Bond::rev_order
 * (well form 'mg3').
 */
#pragma once

#include <cstddef>

namespace RevWellTableV2 {

struct Entry {
    int z1, z2;
    double s, ca, beta, dr0;
};

struct OrderEntry {
    int z1, z2;
    int order;                 ///< nominal bond order of the class-A system this row was fitted on
    double s, ca, beta, dr0;
};

"""

HEADER_V2_TAIL = r"""
/// The element-pair entry, or nullptr when the class-A set has no data for that pair - the
/// caller then falls back to the delivered Gaussian and says so at verbosity 2.
inline const Entry* find(int za, int zb)
{
    const int z1 = za < zb ? za : zb;
    const int z2 = za < zb ? zb : za;
    for (std::size_t k = 0; k < kPairCount; ++k)
        if (kPairEntries[k].z1 == z1 && kPairEntries[k].z2 == z2)
            return &kPairEntries[k];
    return nullptr;
}

/// The bond-order-resolved parameters of an element pair at a CONTINUOUS order, linearly
/// interpolated between the nominal orders the class-A set covers for that pair and CLAMPED at
/// the lowest / highest of them.  Returns false when the pair has no order row at all (the
/// caller then uses find() above, i.e. the element-pair fit).
///
/// The interpolation is Lipschitz in \p order and \p order is a topology constant inside one
/// energy call, so this introduces no geometry derivative and no switch: a pair whose pi order
/// drifts moves its well parameters proportionally, never in a step.  A pair the class-A set
/// covers at ONE order only (e.g. C-H, which has no double or triple member) returns that row
/// for every order, i.e. it behaves exactly like the element-pair table.
inline bool findOrder(int za, int zb, double order, double& s, double& ca, double& beta, double& dr0)
{
    const int z1 = za < zb ? za : zb;
    const int z2 = za < zb ? zb : za;
    const OrderEntry* lo = nullptr;
    const OrderEntry* hi = nullptr;
    for (std::size_t k = 0; k < kOrderCount; ++k) {
        const OrderEntry& e = kOrderEntries[k];
        if (e.z1 != z1 || e.z2 != z2)
            continue;
        const double o = static_cast<double>(e.order);
        if (o <= order && (!lo || o > static_cast<double>(lo->order)))
            lo = &e;
        if (o >= order && (!hi || o < static_cast<double>(hi->order)))
            hi = &e;
    }
    if (!lo && !hi)
        return false;
    if (!lo)
        lo = hi;                       // order below the lowest fitted one -> clamp
    if (!hi)
        hi = lo;                       // order above the highest fitted one -> clamp
    double t = 0.0;
    if (hi->order != lo->order)
        t = (order - static_cast<double>(lo->order))
            / static_cast<double>(hi->order - lo->order);
    s = lo->s + t * (hi->s - lo->s);
    ca = lo->ca + t * (hi->ca - lo->ca);
    beta = lo->beta + t * (hi->beta - lo->beta);
    dr0 = lo->dr0 + t * (hi->dr0 - lo->dr0);
    return true;
}

} // namespace RevWellTableV2
"""


def write_header_v2(path, out):
    """rev-gfnff 3a(iii) step 2 + 3b (Claude Generated, Sep 20, 2026): the free-curvature table.

    Two arrays: the element-pair fit (well form 'mg2') and the bond-order-resolved fit
    ('mg3', interpolated at runtime on the continuous order).  Everything else about the file -
    the struct layout, find(), findOrder() - is fixed C++ and is written verbatim, so the
    generated header is a drop-in replacement for the previous one."""
    inv = {v: k for k, v in ELEM.items()}
    prows = []
    for p_, v in sorted(out["pairs"].items()):
        e1, e2 = p_.split("-")
        if e1 not in inv or e2 not in inv:
            continue
        z1, z2 = sorted((inv[e1], inv[e2]))
        prows.append((z1, z2, p_, v))
    orows = []
    for key, v in sorted(out.get("pair_orders", {}).items()):
        e1, e2 = v["pair"].split("-")
        if e1 not in inv or e2 not in inv:
            continue
        z1, z2 = sorted((inv[e1], inv[e2]))
        orows.append((z1, z2, v["order"], v["pair"], v))
    orows.sort(key=lambda t: (t[0], t[1], t[2]))
    with open(path, "w") as f:
        f.write(HEADER_V2_DOC % (out["md5"], out["rest"]))
        f.write("static constexpr Entry kPairEntries[] = {\n")
        for z1, z2, name, v in prows:
            f.write("    { %2d, %2d, %10.6f, %10.6f, %10.6f, %10.6f },   // %s (n = %d, fit rms %.2f)\n"
                    % (z1, z2, v["mg2_s"], v["mg2_ca"], v["mg2_beta"], v["mg2_dr0"],
                       name, v["n"], v["mg2_rms"]))
        f.write("};\n\nstatic constexpr std::size_t kPairCount = sizeof(kPairEntries) / sizeof(kPairEntries[0]);\n\n")
        f.write("static constexpr OrderEntry kOrderEntries[] = {\n")
        for z1, z2, o, name, v in orows:
            f.write("    { %2d, %2d, %d, %10.6f, %10.6f, %10.6f, %10.6f },   // %s order %d (n = %d, %s, fit rms %.2f)\n"
                    % (z1, z2, o, v["mg2_s"], v["mg2_ca"], v["mg2_beta"], v["mg2_dr0"],
                       name, o, v["n"], "+".join(v["members"]), v["mg2_rms"]))
        f.write("};\n\nstatic constexpr std::size_t kOrderCount = sizeof(kOrderEntries) / sizeof(kOrderEntries[0]);\n")
        f.write(HEADER_V2_TAIL)


if __name__ == "__main__":
    main()
