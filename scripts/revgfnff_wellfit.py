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
               "-gfnff.reuse_topology_check", "false"] + list(extra)
        env = dict(os.environ)
        env["CURCUMA_SHAREDUMP"] = "1"
        p = subprocess.run(cmd, capture_output=True, text=True, cwd=d, env=env)
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
                    blocks=blocks[1:], charge=ref["charge"],
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
    return dict(r=r, Eref=Eref, E=E, rest=rest, drift=coul - coul[i0], i0=i0,
                alpha=float(alpha), r0=float(r0), kb=float(kb), gres=res, pair=key,
                elems=tuple(sorted((e1, e2))))


# --------------------------------------------------------------------------- forms
def mg(D, beta, K, x):
    a = math.sqrt(K / (2.0 * D))
    phi = a * x + beta * x * x
    y = np.exp(-np.clip(phi, -50.0, 200.0))
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


def fit(d, form, rest):
    x = d["r"] - d["r0"]
    ref = d["Eref"] - d["Eref"][d["i0"]]
    msk = np.arange(len(x)) >= d["i0"] - 1
    Kc = 2.0 * d["alpha"] * d["kb"]

    def tot(p):
        D = p[0]
        if D <= 0:
            return None
        e = mg(D, abs(p[1]), Kc, x) if form == "mg" else em(D, abs(p[1]) + 1e-6, Kc, x)
        if e is None or not np.all(np.isfinite(e)):
            return None
        return rest + e

    def obj(p):
        t = tot(p)
        if t is None:
            return 1e9
        return float(np.sqrt(np.mean((((t - t[d["i0"]]) - ref)[msk]) ** 2)))

    starts = ([[d["kb"] * 1.2, 0.5], [d["kb"] * 1.5, 0.2], [d["kb"], 1.0]] if form == "mg"
              else [[d["kb"] * 1.2, 1.2], [d["kb"] * 1.5, 0.8], [d["kb"], 2.0]])
    best = None
    for s in starts:
        p, v = nelder(obj, s, [0.1 * abs(q) + 0.05 for q in s])
        p, v = nelder(obj, p, [0.02 * abs(q) + 0.01 for q in p])
        if best is None or v < best[1]:
            best = (p, v)
    p, v = best
    t = tot(p)
    return dict(D=float(p[0]), tail=float(abs(p[1])), s=float(p[0] / d["kb"]), rms=float(v),
                metrics=[float(q) for q in metrics(d["r"], t)],
                join=float(t[-1] - rest[-1]))


ELEM = {1: "H", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 14: "Si", 15: "P", 16: "S", 17: "Cl",
        35: "Br", 53: "I"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", default=str(REPO / "build_rev" / "curcuma"))
    ap.add_argument("--json", default="wellfit.json")
    ap.add_argument("--header", default="", help="write the compiled C++ table here")
    ap.add_argument("--rest", choices=["qfrozen", "raw"], default="qfrozen")
    ap.add_argument("--workdir", default="")
    ap.add_argument("--systems", nargs="*", default=None)
    a = ap.parse_args()
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
        row = dict(r0=d["r0"], alpha=d["alpha"], kb=d["kb"], gres=d["gres"],
                   rms_delivered=rms_del, pair="-".join(d["elems"]),
                   mg=fit(d, "mg", rest), em=fit(d, "em", rest))
        out["bonds"][b] = row
        print("%-22s pair %-6s kb %8.2f alpha %6.3f | delivered rms %7.2f | MG s %5.3f beta %6.3f rms %6.2f"
              " | EM s %5.3f sigma %6.3f rms %6.2f"
              % (b, row["pair"], d["kb"], d["alpha"], rms_del,
                 row["mg"]["s"], row["mg"]["tail"], row["mg"]["rms"],
                 row["em"]["s"], row["em"]["tail"], row["em"]["rms"]))
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
            em_rms=float(np.median([r["em"]["rms"] for r in rows])))
    json.dump(out, open(a.json, "w"), indent=1)
    print("# wrote", a.json, "-", len(out["bonds"]), "bond types,", len(out["pairs"]), "element pairs")
    if a.header:
        write_header(a.header, out)
        print("# wrote", a.header)


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


if __name__ == "__main__":
    main()
